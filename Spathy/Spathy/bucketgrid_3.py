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
    def __init__(self, spara, soiltypefile=None, sdata=None, cmask=None, D=None,  OrgSto=None, 
                 MaxPond=None, SatRatio=None, PondSto=None, outputs=False):
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
            PondSto - initial pond storage [m], scalar or (n x m) array
            outputs - True saves outputgrids to list elements for each timestep

        CHANGES:
            14.6.2017 Samuli: added outputs as option
        """
        self.SurfSto=None
        if sdata is not None:  # when spatial_soil = yes

            # print '******BucketGrid: spatial soil****'
            GridShape = np.shape(sdata)
            # print GridShape
            gridmask = np.ones(GridShape)
            gridmask[np.isnan(sdata)] = np.NaN
            gridmask[sdata == -1.0] = np.NaN  # waterbodies
            poros = np.empty(GridShape)*gridmask
            Fc = np.empty(GridShape)*gridmask
            Wp = np.empty(GridShape)*gridmask
            Ksat = np.empty(GridShape)*gridmask
            beta = np.empty(GridShape)*gridmask

            """ read soil hydrol. properties from parameter file """
            print('BucketGrid __ini__',soiltypefile)
            soiltypefile=str(soiltypefile)
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

        if OrgSto is None:
            self.MaxStoTop = spara['orgsto'] * gridmask
        else:
            self.MaxStoTop = OrgSto

        self.WatStoTop = self.MaxStoTop.copy()
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

        if PondSto is None:
            PondSto = spara['pondsto']
        self.PondSto = np.minimum(PondSto*gridmask, self.MaxPond)

        self.Wliq = self.poros*self.WatSto / self.MaxWatSto
        self.Wair = self.poros - self.Wliq
        self.Sat = self.Wliq/self.poros
        self.Rew = np.minimum((self.Wliq - self.Wp) / (self.Fc - self.Wp + eps), 1.0)

        # create dictionary of empty lists for saving results
        if outputs:
            self.results = {'Infil': [], 'Retflow': [], 'Drain': [], 'Roff': [], 'ET': [],
            'MBE': [], 'Retflow': [], 'Wliq': [], 'PondSto': [], 'WatStoTop': []}


    def watbal(self, dt=1, rr=0.0, tr=0.0, evap=0.0, retflow=0.0):
        """
        Computes 2-layer bucket model water balance for one timestep dt
        IN:
            dt [unit]
            rr = potential infiltration [m]
            tr = transpiration from root zone [m]
            evap = evaporation from top layer [m]
            retflow = return flow [m]
        OUT:
            infil [m] - infiltration [m]
            drain [m] - percolation /leakage
            roff [m] - surface runoff
            et [m] - evapotranspiration
            mbe [m] - mass balance error
        """
        gridshape = np.shape(self.Wliq)  # rows, cols
    
        if np.shape(retflow) != gridshape:
            retflow = retflow * np.ones(gridshape)
        if np.shape(rr) != gridshape:
            rr = rr * np.ones(gridshape)
        
        rr0 = rr.copy()
       
       # add current Pond storage to rr
        PondSto0 = self.PondSto
        rr += self.PondSto
        self.PondSto = np.zeros(gridshape)
        
        WatSto0 = self.WatSto
        WatStoTop0 = self.WatStoTop
        

        # compute interception to top layer, route remaining to root zone
        interc = np.maximum(0.0, (self.MaxStoTop - self.WatStoTop))\
                    * (1.0 - np.exp(-(rr / self.MaxStoTop)))

        rr = rr - interc
        print('rr', np.nanmax(rr), np.nanmin(rr))     
        #  update state
        dSto = interc - evap
        self.WatStoTop = np.maximum(0.0, self.WatStoTop + dSto)        
        del dSto
        
        # ********* compute bottom layer (root zone) water balance ***********
        
        # drainage. gridcells where retflow >0, set drain to zero. This delays drying of
        # cells which receive water from surroundings
        drain = np.minimum(self.Ksat * dt, np.maximum(0.0, (self.Wliq - self.Fc))*self.D)
        drain[retflow > 0.0] = 0.0

        # inflow is restricted by potential inflow or available pore space
        Qin = (retflow + rr)  # m, pot. inflow
        inflow = np.minimum(Qin, self.MaxWatSto - self.WatSto + drain + tr)

        dSto = (inflow - drain - tr)
        self.WatSto = np.minimum(self.MaxWatSto, np.maximum(self.WatSto + dSto, eps))
        
        # ***** route inflow excess ***********             
        # if inflow excess, update first top layer storage
        exfil = Qin - inflow + eps
        to_top_layer = np.minimum(exfil, self.MaxStoTop - self.WatStoTop + eps)
        # self.WatStoTop = self.WatStoTop + to_top_layer
        self.WatStoTop += to_top_layer

        # ... and pond storage ...
        to_pond = np.minimum(exfil - to_top_layer, self.MaxPond - self.PondSto + eps)
        self.PondSto += to_pond
 
        # ... and remaining routed as surface runoff
        roff = exfil - to_top_layer - to_pond
        
        # update state variables:
        self.setState()
        
        # mass balance error [m]
        mbe = (self.WatSto - WatSto0)  + (self.WatStoTop - WatStoTop0) + (self.PondSto - PondSto0) \
            - (rr0 + retflow - tr - evap - drain - roff)
        # print('mbe', np.nanmax(mbe))
        
        # append results to lists; use only for testing small grids!
        if hasattr(self, 'results'):
            self.results['Infil'].append(inflow - retflow)   # infiltration through bottom boundary
            self.results['Retflow'].append(retflow)         # return flow from below 
            self.results['Roff'].append(roff)
            self.results['Drain'].append(drain)
            self.results['ET'].append(tr + evap)
            self.results['MBE'].append(mbe)
            self.results['Wliq'].append(self.Wliq)
            self.results['PondSto'].append(self.PondSto)
            self.results['WatStoTop'].append(self.WatStoTop)
        
        return inflow, roff, drain, tr+evap, mbe

    def setState(self):
        """ updates state variables"""
#        self.WatSto = np.minimum(self.MaxWatSto,
#                                 np.maximum(self.WatSto + dWat, eps))
#This was commented!!!
#        self.SurfSto = np.minimum(self.MaxPond, self.SurfSto + dSur)

        self.Wliq = self.poros*self.WatSto / self.MaxWatSto
        self.Wair = self.poros - self.Wliq
        self.Sat = self.Wliq / self.poros
        self.Rew = np.maximum(0.0,
              np.minimum((self.Wliq - self.Wp) / (self.Fc - self.Wp + eps), 1.0))

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
    cfg.read(soilparamfile,encoding = "ISO-8859-1")
    # print cfg
    pp = {}
    for s in cfg.sections():
        section = s
        pp[section] = {}
        for k, v in cfg.items(section):
            key = k
            val = v
            if key == 'gtk_code':  # and len(val)>1:
                val = map(int, val.split(','))
                pp[section][key] = val
            else:
                pp[section][key] = float(val)
    del cfg

    return pp
