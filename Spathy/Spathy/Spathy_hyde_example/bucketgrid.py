# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 10:52:25 2017

@author: slauniai
"""
import numpy as np
eps = np.finfo(float).eps

class BucketGrid(object):
    """
    Single-layer soil water bucket model for gridded use.
    
    17.11.2016: this version reads soil type-dependent properties from 
    'soiltypefile' and then assigns properties to each
    grid cell according to values in 'sdata' matrix and 'soil_id' in soiltypefile

    """
    def __init__(self, spara, soiltypefile=None, sdata=None, cmask=None, D=None,
                 MaxPond=None, SatRatio=None, SurfSto=None, outputs=False):
        """
        Initializes BucketGrid:
        REQUIRED:
            spara - 'Bucket' parameters from ini-file
        OPTIONAL: if spatial soil properties   
            soiltypefile - filepath to soiltypefile
            sdata - soil map (topsoil) code raster, n x m array 
            cmask - catchment mask (n x m ) array. !* NEEDED IF only spara is given as input, otherwise computed from sdata
            D - bucket depth (m), (n x m ) array or scalar
            MaxPond - maximum ponding depth above bucket [m] scalar or (n x m) array                
            WatSto - initial storage [m] scalar or (n x m) array
            SurfSto - initial surface storage [m], scalar or (n x m) array
            outputs - True saves outputgrids to list elements for each timestep

        CHANGES:
            14.6.2017 Samuli: added outputs as option
        """

        if sdata is not None:  # when spatial_soil = yes

            # print '******BucketGrid: spatial soil****'
            GridShape = np.shape(sdata)
            # print GridShape
            gridmask = np.ones(GridShape)
            gridmask[np.isnan(sdata)] = np.NaN
            gridmask[sdata == -1.0] = np.NaN  # waterbodies
            poros = np.empty(GridShape)*gridmask
            Fc = poros.copy()
            Wp = poros.copy()
            Ksat = poros.copy()
            beta = poros.copy()

            """ read soil hydrol. properties from parameter file """
            pp = read_soil_properties(soiltypefile)
            # print pp

            """ create hydrologic parameter grids by matching sdata code with
            that of soil properties dict pp
            """
            for key in pp.keys():
                # print key
                soilcode = pp[key]['soil_id']
                ix = np.where(sdata == soilcode)
                poros[ix] = pp[key]['poros']
                Fc[ix] = pp[key]['fc']
                Wp[ix] = pp[key]['wp']
                Ksat[ix] = pp[key]['ksat']
                beta[ix] = pp[key]['beta']
        else:  # use constant soil properties throughout
            # print '******BucketGrid: constant soil properties****'
            GridShape = np.shape(cmask)
            gridmask = np.ones(GridShape)*cmask
            poros = spara['poros']*gridmask
            Fc = spara['fc']*gridmask
            Wp = spara['wp']*gridmask
            Ksat = spara['ksat']*gridmask
            beta = spara['beta']*gridmask

        """ set object properties. All will be (n x m) grids """
        if D is None:
            self.D = spara['depth']*gridmask
        else:
            self.D = D*gridmask

        self.poros = poros
        self.Fc = Fc
        self.Wp = Wp

#        """ SAMULI ADD 24.3.2017 """
#        # self.Ksat=Ksat  # saturated hydraulic conductivity [ms-1]
#        self.Ksat0 = Ksat
#        self.Ksat = Ksat*spara['kmult']
#        self.spara = spara
#        """ ------"""
        self.Ksat = Ksat
        self.beta = beta  # hydr. cond. exponent
        self.spara = spara
        self.MaxWatSto = self.D*self.poros  # maximum soil water storage [m]

        if MaxPond is None:
            MaxPond = spara['maxpond']
        self.MaxPond = MaxPond*gridmask

        """
        set buckets initial state:
        SatRatio and/or SurfSto are given either as (n x m) array or scalar, or
        taken from spara['satratio']
        """
        if SatRatio is None:
            SatRatio = spara['satratio']
        self.WatSto = np.minimum(SatRatio*self.D*self.poros, self.D*self.poros)

        if SurfSto is None:
            SurfSto = spara['surfsto']
        self.SurfSto = np.minimum(SurfSto*gridmask, self.MaxPond)

        self.Wliq = self.poros*self.WatSto/self.MaxWatSto
        self.Wair = self.poros - self.Wliq
        self.Sat = self.Wliq/self.poros
        self.Rew = np.minimum((self.Wliq - self.Wp) / (self.Fc - self.Wp + eps), 1.0)

        # create dictionary of empty lists for saving results
        if outputs:
            self.results = {'Infil': [], 'Drain': [], 'Roff': [],
                            'ET': [], 'MBE': [], 'Retflow': [], 'Wliq': [], 'SurfSto': []}

    def watbal(self, dt=1, rr=0.0, et=0.0, drain=0.0, retflow=0.0):
        """
        Computes bucket model water balance for one timestep dt
        IN:
            dt [unit]
            rr = potential infiltration [m]
            et [m]
            drain [m],
            retflow [m]
        OUT:
            infil [m] - infiltration [m]
            drain [m] - percolation /leakage
            roff [m] - surface runoff
            et [m] - evapotranspiration
            mbe [m] - mass balance error
        """

        SurfSto0 = self.SurfSto
        # total inflow for infiltration, drainage and net inflow
        Qin = (rr + retflow + SurfSto0)  # m, pot. inflow
        drain = np.minimum(drain, np.maximum(0.0, (self.Wliq - self.Fc))*self.D)
        #drain = np.minimum(drain, np.maximum(0.0, (self.Wliq - self.Wp))*self.D)
        # infiltr is restricted by availability, Ksat or available pore space
        infil = np.minimum(Qin,
                np.minimum(self.Ksat*dt, self.MaxWatSto - self.WatSto + drain))

        # change in soil & surface water store
        dSto = (infil - drain - et)
        dSur = (rr + retflow)-infil

        # update state variable grids
        self.setState(dSto, dSur)
        roff = np.maximum(0.0, dSur - (self.SurfSto - SurfSto0))

        # mass balance error [m]
        mbe = dSto + (self.SurfSto - SurfSto0) - (rr + retflow - et - drain - roff)

        # append results to lists; use only for testing small grids!
        if hasattr(self, 'results'):
            self.results['Infil'].append(infil)
            self.results['Roff'].append(roff)
            self.results['Drain'].append(drain)
            self.results['ET'].append(et)
            self.results['MBE'].append(mbe)
            self.results['Wliq'].append(self.Wliq)
            self.results['SurfSto'].append(self.SurfSto)

        return infil, roff, drain, et, mbe

    def setState(self, dWat, dSur=0):
        """ updates state variables by computed dSto [m] and dSurf [m] """
        self.WatSto = np.minimum(self.MaxWatSto,
                                 np.maximum(self.WatSto + dWat, eps))
        self.SurfSto = np.minimum(self.MaxPond, self.SurfSto + dSur)

        self.Wliq = self.poros*self.WatSto/self.MaxWatSto
        self.Wair = self.poros - self.Wliq
        self.Sat = self.Wliq/self.poros
        self.Rew = np.maximum(0.0,
              np.minimum((self.Wliq - self.Wp)/(self.Fc - self.Wp + eps), 1.0))

    def hydrCond(self):
        """
        returns hydraulic conductivity [ms-1] based on Campbell -formulation
        """
        k = self.Ksat*self.Sat**(2*self.beta + 3.0)
        return k


def read_soil_properties(soilparamfile):
    """
    reads soiltype-dependent hydrol. properties from text file
    IN: soilparamfile - path to file
    OUT: pp - dictionary of soil parameters
    """
    import configparser

    cfg = configparser.ConfigParser()
    cfg.read(soilparamfile)
    # print cfg
    pp = {}
    for s in cfg.sections():
        section = s.encode('ascii', 'ignore')
        pp[section] = {}
        for k, v in cfg.items(section):
            key = k.encode('ascii', 'ignore')
            val = v.encode('ascii', 'ignore')

            if key == 'gtk_code':  # and len(val)>1:
                val = map(int, val.split(','))
                pp[section][key] = val
            else:
                pp[section][key] = float(val)
    del cfg

    return pp
