# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 16:18:57 2016

@author: slauniai


"""
import pathlib
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import configparser
import timeit

from canopygrid_3 import CanopyGrid # CHANGED 18.12.
from bucketgrid_3 import BucketGrid
from topmodel_3 import Topmodel_Homogenous as Topmodel

from iotools import create_catchment, read_FMI_weather, read_SVE_runoff

import spathypath
#spathy_path = os.path.join('c:', 'c:\pyspace\spathy')
#spathy_path = os.path.join('./')
eps = np.finfo(float).eps  # machine epsilon

""" ************** SpatHy ************************************

Simple spatial hydrology and catchment water balance model.

CONSISTS OF THREE CLASSES, imported from spathy.modules.: 
    CanopyGrid - vegetation and snowpack water storages and flows
    BucketGrid - topsoil bucket model (root zone / topsoil water storage)
    Topmodel - integration to catchment scale using Topmodel -concept
HELPER FUNCTIONS:
    in iotools
 
MAIN PROGRAM:   
    spathy_driver is main program, call it as
    outargs = spathy_driver(spathyparamfile, args)
    
    spathyparamfile - path to parameter file, default is 'spathy_default.ini'
    soil type dependent parameters are in 'soilparam.ini'

NEEDS 2D gis rasters in ascii-grid format

CanopyGrid & BucketGrid can be initialized from gis-data or set to spatially constant properties

ToDo:
    CanopyGrid:
        -include topographic shading to radiation received at canopy top
        -radiation-based snowmelt coefficient
        -add simple GPP-model
    BucketGrid:
        -make soil hydrologic properties (soiltypes.ini) more realistic
        -kasvupaikkatyyppi (multi-NFI) --> soil properties
        -add soil frost model, simplest would be Stefan equation with coefficients modified based on snow insulation
          --> we need snow density algorithm: SWE <-----> depth
    Topmodel:
        -think of definging 'relative m & to grids' (soil-type & elevation-dependent?) and calibrate 'catchment averages'
        -topmodel gives 'saturated zone storage deficit in [m]'. This can be converted to gwl proxy (?) if: 
        local water retention characteristics are known & hydrostatic equilibrium assumes. Look which water retention model was
        analytically integrable (Campbell, brooks-corey?)
    
    Graphics and analysis of results:
        -make ready functions

Changes June 2016:
    -separated sub-models into own modules
    -seasonal cycle of LAI of deciduous trees
    -add species composition (coniferous, deciduous) and their transpiratio parameters
    -to speed up calibration, possibility to 'flatten' 2d arrays to 1d, include only cells within catchment
    
    - re-written calibration scripts    
(C) Samuli Launiainen 10/2016-->    

VERSION 01.11.2017: CHANGES:

1) Changed the way Transpiration & Efloor are computed in CanopyGrid and now returned separately
2) added CO2-response of canopy conductance and CO2 as input variable (iotools.Read_FMI_weather)
3) modified NetCdf-output: new variables Transpi and Efloor

"""


def spathy_driver(setupfile, catchment_id=None, ncf=False, cpy_outputs=False, 
                  bu_outputs=False, top_outputs=False, flatten=False):
    """ 
    ******************** spathy_driver **********************
    
    Normal Spathy run without parameter optimization
    1) reads parameters from 'spathy.ini' (setupfile) 
    2) reads GIS-data, here predefined format for Seurantaverkko -cathcments
    3) reads forcing (defined in 'spathy.ini')
    4) optionally reads runoff measurements (file defined in 'spathy.ini')
    
    5) creates CanopyGrid (cpy), BucketGrid (bu) and Topmodel (top) -objects  within Spathy-object (spa) and temporary outputs
    6) creates netCDF -file for outputs if 'ncf' = True
    
    7) loops through forging period and sequentially solves: i) aboveground water budget (cpy), ii) topsoil water budget (bu),
        iii) catchment water budget, saturated areas and streamflow generation (top)
        
        Note: drainage from topsoil bucket (bu) is allowed only if water content of a cell is above field capacity AND 
        'lower soil layer' of the same cell is unsaturated according to Topmodel. This delays drying of cells that have high 
        topographic wetness index and is thus consistent with 'topmodel' conceptualization.
    
    8) returns following outputs:
        spa - spathy object
        outf - filepath to output netCDF file.
    
    TODO:
        read_Runoff_data(args)
        netCDF writing in 10/30day steps
        implement calibration option - current structure supports that
        netCDF file: open and data access, graphics etc. netCDF can be accessed as dict or DataFrame.
    
    IN:
        setupfile - path to ini-file
        catchment_id - id of catchment, overrides that in ini-file
        ncf - True saves outputs to netCDF-file
        cpy_outputs, bu_, top_ - True saves cpy, bu and top outputs to lists within each object. 
            Use only for testing, memory issue.
        flatten - True flattens 2d arrys to 1d array containing only cells inside catchment. not working with current
                netcdf output file so for calibration only.
    OUT:
        spa - spathy object
        outf - filepath to netCDF-file. if ncf=False, returns None
    """
    start_time = timeit.default_timer()

    """read parameter file into dicts"""
    pgen, pcpy, pbu, ptop = read_setup(setupfile)

    # full path to soil_file
    soil_file=pathlib.Path(pgen['soil_file'])
    print('soil_file: ',soil_file)
    pgen['soil_file'] =  soil_file.resolve()

    # if given, override cathcment_id of setupfile
    if catchment_id:
        pgen['catchment_id'] = catchment_id

    gisdata = create_catchment(pgen['catchment_id'], fpath=pgen['gis_folder'],
                               plotgrids=False, plotdistr=False)

    """ greate SpatHy object """
    spa = SpatHy(pgen, pcpy, pbu, ptop, gisdata, cpy_outputs=cpy_outputs,
                 bu_outputs=bu_outputs, top_outputs=top_outputs, flatten=flatten)
    Nsteps = spa.Nsteps

    """ create netCDF output file """
    if ncf:
        ncf, outf = initialize_netCDF(spa.id, spa.GisData, spa.FORC,
                                      fpath=spa.pgen['output_folder'])
    else:
        outf = None
                                
    """ initialize 3D arrays for temporary outputs """
    r, c = np.shape(spa.cmask)

    #3d array indexing: dim1=time, dim2=rows(lat), dim3=cols(lon). W[1,:,:] --> grid at 1st timestep. 
    #W[:,30,30]=timeserie at cell with index 30,30

    """ ----- MAIN CALCULATION LOOP ----- """

    print('Reading files and init objects [s]: ', timeit.default_timer() - start_time )
    del start_time
    start_time = timeit.default_timer()

    print('******* Running Spathy ********')
    spa._run(0, Nsteps, calibr=False, ncf=ncf)

    print('Loops total [s]: ', timeit.default_timer() - start_time)
    print('********* done *********')

    return spa, outf


"""
******************************************************************************
            ----------- SpatHy model class --------
******************************************************************************
"""


class SpatHy():
    """
    SpatHy model class
    """
    def __init__(self, pgen, pcpy, pbu, ptop, gisdata, cpy_outputs=False,
                 bu_outputs=False, top_outputs=False, flatten=False):

        self.dt = pgen['dt']  # s
        self.id = pgen['catchment_id']
        self.spinup_end = pgen['spinup_end']
        self.pgen = pgen

        """ read forcing data and catchment runoff file """
        FORC = read_FMI_weather(pgen['catchment_id'],
                                pgen['start_date'],
                                pgen['end_date'],
                                sourcefile=pgen['forcing_file'])
        FORC['Prec'] = FORC['Prec'] / self.dt  # mms-1

        self.FORC = FORC
        self.Nsteps = len(FORC)

        # read runoff measurements
        if pgen['runoff_file'] is not '':
            self.Qmeas = read_SVE_runoff(pgen['catchment_id'],
                                         pgen['start_date'], pgen['end_date'])
        else:
            self.Qmeas = None

        self.GisData = gisdata
        self.cmask = self.GisData['cmask']

        cmask= self.cmask.copy()
        lai_conif = gisdata['LAI_conif'].copy()
        lai_decid = gisdata['LAI_decid'].copy()
        hc = gisdata['hc'].copy()
        cf = gisdata['cf'].copy()
        sdata = gisdata['soil'].copy()
        flowacc = gisdata['flowacc'].copy()
        slope = gisdata['slope'].copy()        
        
        """ --- flatten 2d arrays and omit cells outside catchment ---"""
        if flatten:
            ix = np.where(np.isfinite(cmask))
            cmask = cmask[ix].copy()
            lai_conif = lai_conif[ix].copy()
            lai_decid = lai_decid[ix].copy()
            hc = hc[ix].copy()
            cf = cf[ix].copy()
            sdata = sdata[ix].copy()
            flowacc = flowacc[ix].copy()
            slope = slope[ix].copy()

#            cmask, _, _ = matrix_to_array(cmask)
#            lai_conif, _, _ = matrix_to_array(lai_conif)
#            lai_decid, _, _ = matrix_to_array(lai_decid)
#            hc, _, _ = matrix_to_array(hc)
#            cf, _, _ = matrix_to_array(cf)
#            sdata, _, _ = matrix_to_array(sdata)
#            flowacc, _, _ = matrix_to_array(flowacc)
#            slope, _, _ = matrix_to_array(slope)
         
        """--- initialize CanopyGrid and BucketGrid ---"""
        if pgen['spatial_cpy'] == 'no':  # spatially constant stand properties
            self.cpy = CanopyGrid(pcpy, cmask=cmask, outputs=cpy_outputs)
        else:
            self.cpy = CanopyGrid(pcpy, lai_conif=lai_conif, lai_decid = lai_decid,
                                  cf=cf, hc=hc, outputs=cpy_outputs)

        if pgen['spatial_soil'] == 'no':  # spatially constant soil properties
            self.bu = BucketGrid(pbu, cmask=cmask, outputs=bu_outputs)
        else:
            self.bu = BucketGrid(pbu, soiltypefile=pgen['soil_file'], sdata=sdata, outputs=bu_outputs)

        """ --- initialize homogenous topmodel --- """
        self.top=Topmodel(ptop, self.GisData['cellsize']**2, cmask,
                          flowacc, slope, outputs=top_outputs)
                          
    def _run(self, fstep, Nsteps, calibr=False, ncf=None):
        """ 
        Runs Spathy
        IN:
            fstep - index of starting point [int]
            Nsteps - number of timesteps [int]
            calibr - set True for parameter optimization, returns Qmod
            ncf - netCDF -file handle, for outputs
        OUT:
            res - modeled streamflow at catchment outlet [mm/d]
        """
        dt = self.dt

        # for calibration run, return res
        if calibr:
             res={'Qm':[None]*self.Nsteps}  # 'RR':[None]*self.Nsteps, 'ET':[None]*self.Nsteps, 'Inter':[None]*self.Nsteps, 'Mbet':[None]*self.Nsteps}
             print('M', self.top.M, 'to', self.top.To)
        
        RR = 0.0 # initial value for recharge [m]
        for k in range(fstep, fstep + Nsteps):
            print('k=' + str(k))
            # forcing
#            doy = self.FORC['doy'].iloc[k]; ta = self.FORC['T'].iloc[k]
#            vpd = self.FORC['VPD'].iloc[k]; rg = self.FORC['Rg'].iloc[k]
#            par = self.FORC['Par'].iloc[k]; prec = self.FORC['Prec'].iloc[k]
#            co2 = self.FORC['CO2'].iloc[k] # ADDED 1.11.


            # forcing
            doy = self.FORC['doy'].iloc[k]; ta = self.FORC['T'].iloc[k]
            vpd = self.FORC['VPD'].iloc[k]; rg = self.FORC['Rg'].iloc[k]
            par = self.FORC['Par'].iloc[k]; prec = self.FORC['Prec'].iloc[k]
            co2 = self.FORC['CO2'].iloc[k]; u = 2.0            
            #u = self.FORC['U'].iloc[k];  u = 2.0;
            #beta = self.bu.Wliq / self.bu.poros
            beta0 = self.bu.WatStoTop / self.bu.MaxStoTop
            
            # run Topmodel, take recharge RR from prev. timestep
            # baseflow [m], returnflow grid [m], sat.area [-]
            qb, qr, fsat = self.top.run_timestep(RR)

            # run CanopyGrid
            potinf, trfall, interc, evap, et, transpi, efloor, mbe = \
                self.cpy.run_timestep(doy, dt, ta, prec, rg, par, vpd, U=u, CO2=co2,
                                      beta=beta0, Rew=self.bu.Rew, P=101300.0)

            # run BucketGrid water balance
            infi, infi_ex, drain, _, mbes = self.bu.watbal(dt=dt, rr=1e-3*potinf, tr=1e-3*transpi,
                                                           evap=1e-3*efloor, retflow=qr)

            # catchment average [m per unit area] saturation excess --> goes to stream
            # as surface runoff
            qs = np.nansum(infi_ex)*self.top.CellArea / self.top.CatchmentArea

            # catchment average ground water recharge [m per unit area]
            RR = np.nansum(drain)*self.top.CellArea / self.top.CatchmentArea

            """ update outputs """
            if ncf:
                # writes to netCDF -file at every timestep; bit slow - should accumulate into temporary variables and save every 10 days?                
                # canopygrid       
                ncf['cpy']['W'][k,:,:] = self.cpy.W
                ncf['cpy']['SWE'][k,:,:] = self.cpy.SWE
                ncf['cpy']['Trfall'][k,:,:] = trfall 
                ncf['cpy']['Potinf'][k,:,:] = potinf
                ncf['cpy']['ET'][k,:,:] = et
                
                # NEW 1.11.
                ncf['cpy']['Transpi'][k,:,:] = transpi
                ncf['cpy']['Efloor'][k,:,:] = efloor
                #                 
                ncf['cpy']['Evap'][k,:,:] = evap
                ncf['cpy']['Inter'][k,:,:] = interc
                ncf['cpy']['Mbe'][k,:,:] = mbe              
                
                # bucketgrid
                ncf['bu']['Drain'][k,:,:] = drain
                ncf['bu']['Infil'][k,:,:] = infi
                ncf['bu']['Wliq'][k,:,:] = self.bu.Wliq
                ncf['bu']['Wsto'][k,:,:] = self.bu.WatSto
                ncf['bu']['SurSto'][k,:,:] = self.bu.SurfSto
                ncf['bu']['Mbe'][k,:,:] = mbes
                
                # topmodel
                ncf['top']['Qb'][k] = qb
                # ncf['top']['Qr'][k] = qr
                ncf['top']['Qs'][k] = qs
                ncf['top']['Qt'][k] = qb + qs  # total runoff
                ncf['top']['R'][k] = RR
                ncf['top']['fsat'][k] = fsat
                ncf['top']['S'][k] = self.top.S

            if calibr:  # calibration run, return streamflow
                res['Qm'][k] = 1e3*(qb + qs)

        # end of time loop
        if ncf:
            ncf.close()  # close netCDF-file

        if calibr:  # return streamflow in dataframe
            res = pd.DataFrame.from_dict(res)
            res.index = self.FORC.index

            return res


""" Functions for data input and output """

def read_setup(inifile):
    """
    reads Spathy.ini parameter file into pp dict
    """
    print('spathy_4 read_setup inifile:', inifile)
    cfg = configparser.ConfigParser()
    cfg.read(str(inifile))

    pp = {}
    for s in cfg.sections():
        print(s)
        #section = s.encode('ascii', 'ignore')
        section=s
        print(section)
        pp[section] = {}
        for k, v in cfg.items(section):
            print(k,v)
            #key = k.encode('ascii', 'ignore')
            #val = v.encode('ascii', 'ignore')
            key=k
            val=v
            print(key,val)
            if section == 'General':  # 'general' section
                pp[section][key] = val
            else:
                pp[section][key] = float(val)
    pp['General']['dt'] = float(pp['General']['dt'])

    pgen = pp['General']
    pcpy = pp['CanopyGrid']
    pbu = pp['BucketGrid']
    ptop = pp['Topmodel']

    return pgen, pcpy, pbu, ptop

def initialize_netCDF(ID, gis, forc, roff=None, fpath='results', fname=None):
    """
    SpatHy netCDF4 format output file initialization
    IN:
        ID -catchment id as int or str
        gis - GisData dict
        forc - forcing data (pd.dataframe)
        roff - measured runoff (pd.Series)
        fpath - path for saving results
        fname - filename
    OUT:
        ncf - netCDF file handle. Initializes data
        ff - netCDF filename incl. path
    LAST EDIT 1.11.
    """
    from netCDF4 import Dataset, date2num  # , num2date
    from datetime import datetime

    # dimensions
    dlat, dlon = np.shape(gis['cmask'])
    dtime = None

    if fname:
        ff = pathlib.Path(spathypaths.pathy_path, fpath, fname)
        print(ff)
    else:
        ff = pathlib.Path(spathypath.spathy_path, fpath, 'Spathy_ch' + str(ID) + '.nc')
        print('initialize_netCDF:',fpath,str(ID),fname,ff)
    # create dataset & dimensions
    dirpath = pathlib.Path(spathypath.spathy_path, fpath)
    #The 'exist_ok' named parameter/keyword appeared on python 3.5 
    #The generic solution is to check if the path exists before creating it.
    if not dirpath.exists():
      dirpath.mkdir(parents=True)
    ncf = Dataset(str(ff), 'w')
    ncf.description = 'SpatHy results. Catchment : ' + str(ID)
    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'SpatHy -model v.0.99'

    ncf.createDimension('dtime', dtime)
    ncf.createDimension('dlon', dlon)
    ncf.createDimension('dlat', dlat)

    # create variables into base and groups 'forc','eval','cpy','bu','top'
    # call as createVariable(varname,type,(dimensions))
    time = ncf.createVariable('time', 'f8', ('dtime',))
    time.units = "days since 0001-01-01 00:00:00.0"
    time.calendar = 'standard'

    lat = ncf.createVariable('lat', 'f4', ('dlat',))
    lat.units = 'ETRS-TM35FIN'
    lon = ncf.createVariable('lon', 'f4', ('dlon',))
    lon.units = 'ETRS-TM35FIN'

    tvec = [k.to_datetime() for k in forc.index]
    time[:] = date2num(tvec, units=time.units, calendar=time.calendar)
    lon[:] = gis['lon0']
    lat[:] = gis['lat0']

    # evaluation data
    Roff = ncf.createVariable('/eval/Roff', 'f4', ('dtime',))
    Roff.units = 'meas. streamflow [mm]'
    if roff:
        Roff[:] = roff.values

    # CanopyGrid outputs
    W = ncf.createVariable('/cpy/W', 'f4', ('dtime', 'dlat', 'dlon',))
    W.units = 'canopy storage [mm]'
    SWE = ncf.createVariable('/cpy/SWE', 'f4', ('dtime', 'dlat', 'dlon',))
    SWE.units = 'snow water equiv. [mm]'
    Trfall = ncf.createVariable('/cpy/Trfall', 'f4', ('dtime', 'dlat', 'dlon',))
    Trfall.units = 'throughfall [mm]'
    Inter = ncf.createVariable('/cpy/Inter', 'f4', ('dtime', 'dlat', 'dlon',))
    Inter.units = 'interception [mm]'
    Potinf = ncf.createVariable('/cpy/Potinf', 'f4', ('dtime', 'dlat', 'dlon',))
    Potinf.units = 'pot. infiltration [mm]'
    ET = ncf.createVariable('/cpy/ET', 'f4', ('dtime', 'dlat', 'dlon',))
    ET.units = 'dry-canopy et. [mm]'
    Transpi = ncf.createVariable('/cpy/Transpi', 'f4', ('dtime', 'dlat', 'dlon',))
    Transpi.units = 'transpiration [mm]'
    Efloor = ncf.createVariable('/cpy/Efloor', 'f4', ('dtime', 'dlat', 'dlon',))
    Efloor.units = 'forest floor evap. [mm]'
    Evap = ncf.createVariable('/cpy/Evap', 'f4', ('dtime', 'dlat', 'dlon',))
    Evap.units = 'interception evap. [mm]'
    Mbe = ncf.createVariable('/cpy/Mbe', 'f4', ('dtime', 'dlat', 'dlon',))
    Mbe.units = 'mass-balance error [mm]'

    # BucketGrid outputs
    Wliq = ncf.createVariable('/bu/Wliq', 'f4', ('dtime', 'dlat', 'dlon',))
    Wliq.units = 'vol. water cont. [m3m-3]'
    Wsto = ncf.createVariable('/bu/Wsto', 'f4', ('dtime', 'dlat', 'dlon',))
    Wsto.units = 'water storage [mm]'
    SurSto = ncf.createVariable('/bu/SurSto', 'f4', ('dtime', 'dlat', 'dlon',))
    SurSto.units = 'pond storage [mm]'
    Infil = ncf.createVariable('/bu/Infil', 'f4', ('dtime', 'dlat', 'dlon',))
    Infil.units = 'infiltration [mm]'
    Drain = ncf.createVariable('/bu/Drain', 'f4', ('dtime', 'dlat', 'dlon',))
    Drain.units = 'drainage [mm]'
    Mbe = ncf.createVariable('/bu/Mbe', 'f4', ('dtime', 'dlat', 'dlon',))
    Mbe.units = 'mass-balance error [mm]'

    # topmodel outputs
    Qt = ncf.createVariable('/top/Qt', 'f4', ('dtime',))
    Qt.units = 'streamflow[m]'
    Qb = ncf.createVariable('/top/Qb', 'f4', ('dtime',))
    Qb.units = 'baseflow [m]'
    Qr = ncf.createVariable('/top/Qr', 'f4', ('dtime',))
    Qr.units = 'returnflow [m]'
    Qs = ncf.createVariable('/top/Qs', 'f4', ('dtime',))
    Qs.units = 'surface runoff [m]'
    R = ncf.createVariable('/top/R', 'f4', ('dtime',))
    R.units = 'average recharge [m]'
    S = ncf.createVariable('/top/S', 'f4', ('dtime',))
    S.units = 'average sat. deficit [m]'
    fsat = ncf.createVariable('/top/fsat', 'f4', ('dtime',))
    fsat.units = 'saturated area fraction [-]'

    return ncf, ff
