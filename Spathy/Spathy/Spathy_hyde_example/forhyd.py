# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:13:45 2017

@author: slauniai

THIS COMBINES CANOPYGRID & BUCKETGRID FOR SOLVING SINGLE SITE WATER BALANCE AND WATER FLUXES.
USES HYYTIÄLÄ DATA AS FORCING AND COMPARES SOME OUTPUTS TO MEASUREMENTS
"""
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import timeit
from scipy import stats

from canopygrid_3 import CanopyGrid
from bucketgrid import BucketGrid

from iotools import read_FMI_weather, read_HydeDaily

eps = np.finfo(float).eps  # machine epsilon
spathy_path = os.path.join('c:', 'c:\pyspace\spathy')


def Hyde_eval():
    setupfile = 'c:\\pyspace\\spathy\\ini\\spathy_hyde.ini'
    pgen, pcpy, pbu = read_setup(setupfile)
    
    """ read forcing data and evaluation data """
    fname = 'c:\\pyspace\\spathy\\data\\HydeDaily2000-2010.txt'   
    dat, FORC = read_HydeDaily(fname)
    FORC['Prec'] = FORC['Prec'] / pgen['dt']  # mms-1
    FORC['T'] = FORC['Ta'].copy()
    cmask = cmask=np.ones(1)
    
    # run model for different parameter combinations, save results into dataframe
    # +/- 15%, 15%, 20%, 30%
    p = [[2.2e-3, 0.8, 2.4, 6.0],
         [1.9e-3, 0.68, 1.95, 4.8 ],
         [3.0e-3, 0.92, 2.85, 7.7]]

    out = []
    for k in range(3):
        a = p[k]
        pcpy['gsref_conif'] = a[0]
        pcpy['gsref_decid'] = 1.6 * a[0]
        pcpy['f'] = a[1]
        pcpy['wmax'] = a[2]
        pcpy['wmaxsnow'] = a[3]
        
        model = ForHyd(pgen, pcpy, pbu, cmask, cpy_outputs=True, bu_outputs=True, FORC=FORC)
        nsteps=len(FORC)
        model._run(0, nsteps)
        
        out.append(model)
        del model
        
    # best model
    Wliq_mod = np.ravel(out[0].bu.results['Wliq'])
    Wliq_low = np.ravel(out[1].bu.results['Wliq'])
    Wliq_high = np.ravel(out[2].bu.results['Wliq'])
    ET_mod = np.ravel(out[0].cpy.results['ET']) 
    ET_low = np.ravel(out[1].cpy.results['ET'])
    ET_high = np.ravel(out[2].cpy.results['ET']) 
    
    E_mod = np.ravel(out[0].cpy.results['Evap']) 
    E_low = np.ravel(out[1].cpy.results['Evap'])
    E_high = np.ravel(out[2].cpy.results['Evap'])

#    SWC_mod = np.ravel(out[0].cpy.results['ET']) 
#    SWC_low = np.ravel(out[1].cpy.results['ET'])
#    SWC_high = np.ravel(out[2].cpy.results['ET'])) 
    
    SWCa = dat['SWCa']
    SWCb = dat['SWCb']
    SWCc = dat['SWCc']
    tvec = dat.index
    et_dry = dat['ET']
    et_dry[dat['Prec']>0.1] = np.NaN
    
    sns.set_style('whitegrid')
    with sns.color_palette('muted'):
        plt.figure()
        
        plt.subplot(2,3,(1,2))
                
        plt.plot(tvec, et_dry, 'o', markersize=4, alpha=0.3, label='meas')
        plt.fill_between(tvec, ET_low, ET_high, facecolor='grey', alpha=0.6, label='range')
        plt.plot(tvec, ET_mod, 'k-', alpha=0.4, lw=0.5, label='mod')
        #plt.xlim([pd.datetime(2003, 10, 1), pd.datetime(2011,1,1)])
        plt.legend(loc=2, fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.ylabel('ET$_{dry}$ (mm d$^{-1}$)', fontsize=8)
        plt.ylim([-0.05, 5.0])
        plt.xlim([pd.datetime(2002, 1, 1), pd.datetime(2011,1,1)])

        # sns.despine()        
        plt.subplot(2,3,(4,5))

        plt.plot(tvec, SWCa, 'o', markersize=4, alpha=0.3,label='meas')
        plt.fill_between(tvec, Wliq_low, Wliq_high, facecolor='grey', alpha=0.6, label='range')
        plt.plot(tvec, Wliq_mod, 'k-',alpha=0.4, lw=0.5, label='mod')        
        #plt.xlim([pd.datetime(2003, 10, 1), pd.datetime(2011,1,1)])
        plt.legend(loc=2, fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.ylabel('$\\theta$ (m$^3$ m$^{-3}$)', fontsize=8)
        plt.xlim([pd.datetime(2002, 1, 1), pd.datetime(2011,1,1)])

        # scatterplot
        plt.subplot(2,3,6)

        meas = np.array(SWCa.values.tolist())
        slope, intercept, r_value, p_value, std_err = stats.linregress(meas, Wliq_mod)
        #print slope, intercept
        xx = np.array([min(meas), max(meas)])
        plt.plot(meas, Wliq_mod, 'o', markersize=5, alpha=0.3)
        plt.plot(xx, slope*xx + intercept, 'k-')
        plt.plot([0.05, 0.45], [0.05, 0.45], 'k--')
        plt.text( 0.07, 0.42, 'y = %.2f x + %.2f' %(slope, intercept), fontsize=8)
        plt.xlim([0.05, 0.45]); plt.ylim([0.05, 0.45])
        
        ax = plt.gca()
        ax.set_yticks([0.1, 0.2, 0.3, 0.4])
        ax.set_xticks([0.1, 0.2, 0.3, 0.4])
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        plt.ylabel('$\\theta$ mod (m$^3$ m$^{-3}$)', fontsize=8)
        
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.xlabel('$\\theta$ meas (m$^3$ m$^{-3}$)', fontsize=8)
        
        # scatterplot
        plt.subplot(2,3,3)

        meas = np.array(et_dry.values.tolist())
        ix = np.where(np.isfinite(meas))
        meas=meas[ix].copy()
        mod = ET_mod[ix].copy()
        slope, intercept, r_value, p_value, std_err = stats.linregress(meas, mod)
        xx = np.array([min(meas), max(meas)])
        plt.plot(meas, mod, 'o', markersize=4, alpha=0.3)

        plt.plot(xx, slope*xx + intercept, 'k-')
        plt.plot([0, 5], [0, 5], 'k--')
        plt.text(0.3, 4.2, 'y = %.2f x + %.2f' %(slope, intercept), fontsize=8)

        plt.xlim([-0.01, 5]); plt.ylim([-0.01, 5])
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        ax = plt.gca()
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        plt.ylabel('ET$_{dry}$ mod (mm d$^{-1}$)', fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.xlabel('ET$_{dry}$ meas (mm d$^{-1}$)', fontsize=8)
        
        plt.savefig('Hyde_validate.pdf')
        plt.savefig('Hyde_validate.png')
        
        # snowpack and throughfall
        
        plt.figure()
        SWE_mod = np.ravel(out[0].cpy.results['SWE'])        
        SWE_low = np.ravel(out[1].cpy.results['SWE'])
        SWE_hi = np.ravel(out[2].cpy.results['SWE'])
        swe_meas = dat['SWE']

        plt.plot(tvec, swe_meas, 'o', markersize=10, alpha=0.3, label='meas')
        plt.fill_between(tvec, SWE_low, SWE_hi, facecolor='grey', alpha=0.6, label='range')
        plt.plot(tvec, SWE_mod, 'k-', alpha=0.4, lw=0.5, label='mod')
        plt.title('SWE'); plt.ylabel('SWE mm')
        
        
        
    return out, dat, FORC


class ForHyd():
    """
    ForHyd model class
    """
    def __init__(self, pgen, pcpy, pbu, cmask, lai_conif=None, lai_decid=None, cf=None, hc=None, 
                 sdata=None, cpy_outputs=False, bu_outputs=False, FORC=None):

        self.dt = pgen['dt']  # s
        self.id = pgen['catchment_id']
        self.spinup_end = pgen['spinup_end']
        self.pgen = pgen

        """ read forcing data """
        if FORC is None:
            FORC = read_FMI_weather(pgen['catchment_id'],
                                    pgen['start_date'],
                                    pgen['end_date'],
                                    sourcefile=pgen['forcing_file'])
            FORC['Prec'] = FORC['Prec'] / self.dt  # mms-1
    
        self.FORC = FORC 
        self.Nsteps = len(self.FORC)

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
        for k in range(fstep, fstep + Nsteps):
            print 'k=' + str(k)
            # forcing
            doy = self.FORC['doy'].iloc[k]; ta = self.FORC['T'].iloc[k]
            vpd = self.FORC['VPD'].iloc[k]; rg = self.FORC['Rg'].iloc[k]
            par = self.FORC['Par'].iloc[k]; prec = self.FORC['Prec'].iloc[k]
            u = self.FORC['U'].iloc[k];  # u = 2.0;

            Rew0 = 1.0
            Rew0 = self.bu.Rew
            beta = self.bu.Wliq / self.bu.poros
            
            # run CanopyGrid
            potinf, trfall, interc, evap, et, transpi, efloor, mbe = \
                self.cpy.run_timestep(doy, dt, ta, prec, rg, par, vpd, U=u, beta=beta, Rew=Rew0, P=101300.0)

            # run BucketGrid water balance
            # first compute potential drainage from buckets
            dra = self.bu.hydrCond()*dt  # m in dt

            infi, infi_ex, drain, _, mbes = self.bu.watbal(dt=dt, rr=1e-3*potinf, et=1e-3*et, drain=dra, retflow=0.0)

            # catchment average [m per unit area] infiltration excess runoff--> goes to stream as surface runoff
            Inex = np.nansum(infi_ex) / np.size(self.bu.poros)

            # catchment average ground water recharge [m per unit area]
            RR = np.nansum(drain) / np.size(self.bu.poros)
            # 'streamflow'
            qs = Inex + RR

            """ update outputs """
            if ncf:
                # writes to netCDF -file at every timestep; bit slow - should accumulate into temporatry variables and save every 10 days?                
                # canopygrid       
                ncf['cpy']['W'][k,:,:] = self.cpy.W
                ncf['cpy']['SWE'][k,:,:] = self.cpy.SWE
                ncf['cpy']['Trfall'][k,:,:] = trfall 
                ncf['cpy']['Potinf'][k,:,:] = potinf
                ncf['cpy']['ET'][k,:,:] = et
                ncf['cpy']['Evap'][k,:,:] = evap
                ncf['cpy']['Transpi'][k,:,:] = transpi
                ncf['cpy']['Efloor'][k,:,:] = efloor
                ncf['cpy']['Inter'][k,:,:] = interc
                ncf['cpy']['Mbe'][k,:,:] = mbe              
                
                # bucketgrid
                ncf['bu']['Drain'][k,:,:] = drain
                ncf['bu']['Infil'][k,:,:] = infi
                ncf['bu']['Wliq'][k,:,:] = self.bu.Wliq
                ncf['bu']['Wsto'][k,:,:] = self.bu.WatSto
                ncf['bu']['SurSto'][k,:,:] = self.bu.SurfSto
                ncf['bu']['Mbe'][k,:,:] = mbes
                ncf['bu']['Infiex'][k] = Inex
                ncf['bu']['Roff'][k] = qs        
                
        # end of time loop
        if ncf:
            ncf.close()  # close netCDF-file



def read_setup(inifile):
    """
    reads Spathy.ini parameter file into pp dict
    """
    import configparser
    inifile = os.path.join(spathy_path, inifile)
    
    cfg = configparser.ConfigParser()
    cfg.read(inifile)

    pp = {}
    for s in cfg.sections():
        section = s.encode('ascii', 'ignore')
        pp[section] = {}
        for k, v in cfg.items(section):
            key = k.encode('ascii', 'ignore')
            val = v.encode('ascii', 'ignore')
            if section == 'General':  # 'general' section
                pp[section][key] = val
            else:
                pp[section][key] = float(val)
    pp['General']['dt'] = float(pp['General']['dt'])

    pgen = pp['General']
    pcpy = pp['CanopyGrid']
    pbu = pp['BucketGrid']

    return pgen, pcpy, pbu
    

def initialize_netCDF(ID, gis, forc, fpath='results', fname=None):
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
    """
    from netCDF4 import Dataset, date2num  # , num2date
    from datetime import datetime

    # dimensions
    dlat, dlon = np.shape(gis['cmask'])
    dtime = None

    if fname:
        ff = os.path.join(spathy_path, fpath, fname)
        print ff
    else:
        ff = os.path.join(spathy_path, fpath, 'Spathy_ch' + str(ID) + '.nc')

    # create dataset & dimensions
    ncf = Dataset(ff, 'w')
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
    Infiex = ncf.createVariable('/bu/Infiex', 'f4', ('dtime',))
    Infiex.units = 'infiltration excess runoff [mm]'
    Roff = ncf.createVariable('/bu/Roff', 'f4', ('dtime',))
    Roff.units = 'total runoff[mm]'

    return ncf, ff


def ECsite_eval(site, setupfile):
    # setupfile = 'c:\\pyspace\\spathy\\ini\\spathy_default.ini'
    pgen, pcpy, pbu = read_setup(setupfile)
    
    """ read forcing data and evaluation data """
     
    dat, FORC = read_daily_ECdata(site)
    FORC['Prec'] = FORC['Prec'] / pgen['dt']  # mms-1
    FORC['T'] = FORC['Ta'].copy()
    cmask = cmask=np.ones(1)
    
    # run model for different parameter combinations, save results into dataframe
    # +/- 15%, 15%, 20%, 30%
    p = [[2.2e-3, 0.8, 2.4, 4.5],
         [1.9e-3, 0.68, 1.95, 5.25 ],
         [2.5e-3, 0.92, 2.85, 9.75]]

    out = []
    for k in range(3):
        a = p[k]
        pcpy['gsref_conif'] = a[0]
        pcpy['gsref_decid'] = 1.4 * a[0]
        pcpy['f'] = a[1]
        pcpy['wmax'] = a[2]
        pcpy['wmaxsnow'] = a[3]
        
        model = ForHyd(pgen, pcpy, pbu, cmask, cpy_outputs=True, bu_outputs=True, FORC=FORC)
        nsteps=len(FORC)
        model._run(0, nsteps)
        
        out.append(model)
        del model
        
    # best model
    Wliq_mod = np.ravel(out[0].bu.results['Wliq'])
    Wliq_low = np.ravel(out[1].bu.results['Wliq'])
    Wliq_high = np.ravel(out[2].bu.results['Wliq'])
    ET_mod = np.ravel(out[0].cpy.results['ET']) 
    ET_low = np.ravel(out[1].cpy.results['ET'])
    ET_high = np.ravel(out[2].cpy.results['ET']) 
    
    E_mod = np.ravel(out[0].cpy.results['Evap']) 
    E_low = np.ravel(out[1].cpy.results['Evap'])
    E_high = np.ravel(out[2].cpy.results['Evap'])

    
    tvec = dat.index
    et_dry = dat['ET']
    et_dry[dat['Prec']>0.1] = np.NaN
    
    sns.set_style('whitegrid')
    with sns.color_palette('muted'):
        plt.figure()
        
        plt.subplot(2,3,(1,2))
                
        plt.plot(tvec, et_dry, 'o', markersize=4, alpha=0.3, label='meas')
        plt.fill_between(tvec, ET_low, ET_high, facecolor='grey', alpha=0.6, label='range')
        plt.plot(tvec, ET_mod, 'k-', alpha=0.4, lw=0.5, label='mod')
        #plt.xlim([pd.datetime(2003, 10, 1), pd.datetime(2011,1,1)])
        plt.legend(loc=2, fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.ylabel('ET$_{dry}$ (mm d$^{-1}$)', fontsize=8)
        plt.ylim([-0.05, 5.0])
        #plt.xlim([pd.datetime(2002, 1, 1), pd.datetime(2011,1,1)])
        plt.title(site)
        # sns.despine()        
        
        # scatterplot
        plt.subplot(2,3,3)

        meas = np.array(et_dry.values.tolist())
        ix = np.where(np.isfinite(meas))
        meas=meas[ix].copy()
        mod = ET_mod[ix].copy()
        slope, intercept, r_value, p_value, std_err = stats.linregress(meas, mod)
        xx = np.array([min(meas), max(meas)])
        plt.plot(meas, mod, 'o', markersize=4, alpha=0.3)

        plt.plot(xx, slope*xx + intercept, 'k-')
        plt.plot([0, 5], [0, 5], 'k--')
        plt.text(0.3, 4.2, 'y = %.2f x + %.2f' %(slope, intercept), fontsize=8)

        plt.xlim([-0.01, 5]); plt.ylim([-0.01, 5])
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        ax = plt.gca()
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        plt.ylabel('ET$_{dry}$ mod (mm d$^{-1}$)', fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.xlabel('ET$_{dry}$ meas (mm d$^{-1}$)', fontsize=8)
        
#        plt.subplot(2,3,(4,5))
#        #plt.plot(tvec, SWCa, 'o', markersize=4, alpha=0.3,label='meas')
#        plt.fill_between(tvec, Wliq_low, Wliq_high, facecolor='grey', alpha=0.6, label='range')
#        plt.plot(tvec, Wliq_mod, 'k-',alpha=0.4, lw=0.5, label='mod')        
#        #plt.xlim([pd.datetime(2003, 10, 1), pd.datetime(2011,1,1)])
#        plt.legend(loc=2, fontsize=8)
#        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
#        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
#        plt.ylabel('$\\theta$ (m$^3$ m$^{-3}$)', fontsize=8)
#        # plt.xlim([pd.datetime(2002, 1, 1), pd.datetime(2011,1,1)])


    return out, dat, FORC


def read_daily_ECdata(site):

    if site=='FISod':
        fpath = 'E:\\Data\\CEIP\\DailyData_Peltoniemi\\FMI_Sodankyla'
        yrs = np.arange(2001, 2010)
        fnames = [ 'SodDaily_%4d.dat' %(k) for k in yrs]
        cols=['doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','Prec','U','Pamb',
              'SWC1','SWC2','Tsoil1', 'Tsoil2', 'Rnetflag', 'Snowdepth']        
    
    if site=='FIKal':
        fpath = 'E:\\Data\\CEIP\\DailyData_Peltoniemi\\FMI_Kalevansuo\\'
        yrs = [2005, 2006, 2007, 2008]
        fnames = ['KalevansuoDaily_%4d.dat' %(k) for k in yrs]
        cols=['doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','empty1','Prec','U','Pamb',
              'WTD', 'Snowdepth', 'Rnetflag']        

    if site=='SEKno':
        fpath = 'E:\\Data\\CEIP\\DailyData_Peltoniemi\\Knottasen\\'
        yrs =[2007]
        fnames = ['KnoDaily_%4d.dat' %(k) for k in yrs]
        cols=['doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','Prec','U','Pamb',
              'SWC1','SWC2','Tsoil1', 'Tsoil2', 'Rnetflag', 'Snowdepth']   
              
    if site=='SENor':
        fpath = 'E:\\Data\\CEIP\\DailyData_Peltoniemi\\Norunda\\'
        yrs = [1996, 1997, 1999, 2003]
        fnames = ['NorDaily_%4d.dat' %(k) for k in yrs]
        cols=['doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','Prec','U','Pamb',
              'SWC1','SWC2','Tsoil1', 'Tsoil2', 'Rnetflag', 'Snowdepth']  

    if site=='SESky2':
        fpath = 'E:\\Data\\CEIP\\DailyData_Peltoniemi\\Skyttorp2\\'
        yrs = [2005]
        fnames = ['Sky2Daily_2005.dat']
        cols=['doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','Prec','U','Pamb',
              'SWC1','SWC2','Tsoil1', 'Tsoil2', 'Rnetflag', 'Snowdepth']  
    
    if site=='FILet':
        fpath = 'E:\\Data\\CEIP\\DailyData_Peltoniemi\\FMI_Lettosuo\\'
        yrs = [2010, 2011, 2012]
        fnames = ['FILet_Daily_%4d.dat' %(k) for k in yrs]
        cols=['year', 'month', 'day', 'doy','NEE','GPP','TER','ET','H','Rnet','Rg', 'Par','Prec_site', 'Prec', 'Ta', 'RH',
              'VPD','CO2','U','Pamb', 'WTD','WTDwest','WTDsouth', 'WTDeast', 'WTDnorth', 'SWC1', 'SWC2', 'empty', 'Ts1', 'Ts2',
              'Ts3', 'Ts4', 'NEEflag', 'ETflag', 'Hflag','Rnetflag']  
    

    dat = pd.DataFrame()
    
    for k in range(len(yrs)):
        fname = os.path.join(fpath, fnames[k])
        tmp = pd.read_csv(fname,sep='\s+',header=None, names=cols)
        tmp['year'] = yrs[k]
        
        dat = dat.append(tmp)
    
    dat['doy'] = dat['doy'].astype(int)
    tvec = pd.to_datetime(dat['year'] * 1000 + dat['doy'], format='%Y%j')
    #dat.drop('year')
    dat.index = tvec
    
    # forcing data
  
    forc=dat[['doy','Ta','VPD','Prec','Par','U']]; forc['Par']= 1/4.6*forc['Par']; forc['Rg']=2.0*forc['Par']
    forc['VPD'][forc['VPD']<=0]=eps
    forc = forc.interpolate()  # fills missing values by linear interpolation  
    #relatively extractable water, from soil moisture
    forc['Rew'] = 1.0
    if site=='SEKno':
        forc['Rg'] = 1.4*forc['Rg']
#    fc=0.30
#    wp=0.10
#    Wliq=dat['SWCa']
#    Rew=np.maximum( 0.0, np.minimum( (Wliq-wp)/(fc - wp + eps), 1.0) )
#    forc['Rew']=Rew     
    
    return dat, forc