"""
Copyright: IIASA (International Institute for Applied Systems Analysis), 2016-2020; CEA (Commissariat a L'Energie Atomique) & UVSQ (Universite de Versailles et Saint-Quentin), 2016
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at), Yann Quilcaille

This software is a computer program whose purpose is to simulate the behavior of the Earth system, with a specific but not exclusive focus on anthropogenic climate change.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""

##################################################
##################################################

import numpy as np
import xarray as xr

from time import perf_counter

from core_fct.fct_loadP import load_all_param
from core_fct.fct_loadD import load_all_hist
from core_fct.fct_genMC import generate_config, generate_drivers
from core_fct.fct_genD import create_hist_drivers, create_scen_drivers


##################################################
## 1. WRAPPER FUNCTION TO RUN SIMULATIONS
##################################################

def run_model(model, inds, Par0=None, For0=None, Ini0=None, For1=None, Ini1=None, mod_region='', nMC=1, LCC='gross', save=None, output=False, **model_args):
    '''
    Wrapper function to run OSCAR. Note that specific cases it may be beneficial to write a dedicated script/function.
    
    Input:
    ------
    model (Model)       model to run (of the class Model in 'cls_main')
    inds (tuple)        starting (ind0), preindustrial (indPI), junction (indH) and ending (indF) years for simulation

    Output:
    -------
    (dict)              if requested, outputs a dict of:
                        'Par' = parameter dataset
                        'For_hist' = historical drivers dataset
                        'Ini_hist' = historical initial conditions dataset
                        'Out_hist' = historical output dataset
                        'For_scen' = scenario drivers dataset
                        'Ini_scen' = scenario initial conditions dataset
                        'Out_scen' = scenario output dataset

    Options:
    --------
    Par0 (str           if str, will load parameters dataset from 'results/' + Par0 + '_Par.nc';
    OR xr.Dataset)      if dataset, will generate MC configuration from Par0 taken as primary parameters;
                        if None, will load_all_param from 'fct_loadP' and then generate MC config;
                        default = None
    For0 (str           if str, will load historical drivers dataset from 'results/' + For0 + '_For_hist.nc';
    OR xr.Dataset)      if dataset, will generate MC drivers from For0 taken as primary historical drivers;
                        if None, will load_all_hist from 'fct_loadD' and then generate MC drivers;
                        default = None
    Ini0 (str           if str, will load historical initial conditions dataset from 'results/' + Ini0 + '_Ini_hist.nc';
    OR xr.Dataset)      if dataset, will use dataset as initial conditions;
                        if None, will set all initial conditions at zero;
                        default = None
    For1 (str           if str, will load provided dataset from 'results/' + For1 + '_For_scen.nc';
    OR xr.Dataset)      if dataset, will consistently connect For_hist and For1;
                        if None, will not run any scenarios;
                        default = None
    Ini1 (str)          if str, will load scenario initial conditions dataset from 'results/' + Ini1 + '_Ini_scen.nc';
                        also if str, will not run historical period;
                        if None, will use last historical time-step as scenario initial conditions;
                        default = None
    mod_region (str)    regional aggregation for the simulation
                        not needed if Par0 and For0 are not None;
                        default = ''
    nMC (int)           number of MC elements;
                        not needed if Par0 and For0 are not None;
                        default = 1
    LCC (str)           load either 'gross' or 'net' land-cover change data;
                        default = 'gross'
    save (str)          name/address of simulation for saving on disk in 'results' (no saving if None);
                        default = None
    output (bool)       if True, return outputs
                        default = False
    **model_args        all other arguments are passed on to model
    '''

    ## getting time counter
    print(model.name + ' loading (historical)')
    t0 = perf_counter()

    ## 0. PARAMETERS
    ## load primary
    if Par0 is None:
        Par0 = load_all_param(mod_region)
    ## and generate MC
    if type(Par0) == xr.Dataset:
        Par = generate_config(Par0, nMC=nMC)
    ## or just read file
    elif type(Par0) == str:
        with xr.open_dataset('results/' + Par0 + '_Par.nc') as TMP: Par = TMP.load()

    ##-----

    ## 1a. HISTORICAL DRIVERS
    ## load primary
    if For0 is None:
        For0 = load_all_hist(mod_region, LCC=LCC)
    ## and generate MC
    if type(For0) == xr.Dataset:
        For_hist = generate_drivers(create_hist_drivers(For0, inds=inds), nMC=nMC)
    ## or just read file
    elif type(Par0) == str:
        with xr.open_dataset('results/' + For0 + '_For_hist.nc') as TMP: For_hist = TMP.load()

    ##-----

    ## move parameters from For to Par (normally, Aland_0)
    Par = xr.merge([Par, For_hist.drop([VAR for VAR in For_hist if 'year' in For_hist[VAR].dims])])
    For_hist = For_hist.drop([VAR for VAR in For_hist if 'year' not in For_hist[VAR].dims])
    
    ## save (Par & For_hist)
    if save is not None:
        Par.to_netcdf('results/' + save + '_Par.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Par})
        For_hist.to_netcdf('results/' + save + '_For_hist.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in For_hist})

    ## printing time counter
    if type(Ini1) == str:
        print('total loading time (historical): {:.0f} seconds'.format(perf_counter() - t0))

    ##-----

    ## 1b. HISTORICAL RUN
    Ini_hist = Out_hist = None
    if not type(Ini1) == str:

        ## create initial steady state
        if Ini0 is None:
            Ini_hist = xr.Dataset()
            for VAR in list(model.var_prog):
                if len(model[VAR].core_dims) == 0: Ini_hist[VAR] = xr.DataArray(0.)
                else: Ini_hist[VAR] = sum([xr.zeros_like(Par[dim], dtype=float) if dim in Par.coords else xr.zeros_like(For_hist[dim], dtype=float) for dim in model[VAR].core_dims])
        ## or just read file
        elif type(Ini0) == str:
            with xr.open_dataset('results/' + Ini0 + '_Ini_hist.nc') as TMP: Ini_hist = TMP.load()
        
        ## save (Ini_hist)
        if save is not None:
            Ini_hist.to_netcdf('results/' + save + '_Ini_hist.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Ini_hist})
        
        ## printing time counter
        print('total loading time (historical): {:.0f} seconds'.format(perf_counter() - t0))

        ##-----

        ## run historical
        Out_hist = model(Ini_hist, Par, For_hist, **model_args)
        
        ## save (Out_hist)
        if save is not None:
            Out_hist.to_netcdf('results/' + save + '_Out_hist.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Out_hist})

    ##-----

    ## 2. SCENARIOS
    For_scen = Ini_scen = Out_scen = None
    if For1 is not None:

        ## getting time counter
        print(model.name + ' loading (scenarios)')
        t0 = perf_counter()

        ## connect to historical
        if type(For1) == xr.Dataset:
            For_scen = create_scen_drivers(For_hist, For1, inds=inds)
        ## or just read file
        elif type(For1) == str:
            with xr.open_dataset('results/' + For1 + '_For_scen.nc') as TMP: For_scen = TMP.load()
        
        ## save (For_scen)
        if save is not None:
            For_scen.to_netcdf('results/' + save + '_For_scen.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in For_scen})

        ##-----

        ## create initial state from historical
        if Ini1 is None:
            Ini_scen = Out_hist.isel(year=-1, drop=True).drop([VAR for VAR in Out_hist if VAR not in list(model.var_prog)])
        ## or just read file
        elif type(Ini1) == str:
            with xr.open_dataset('results/' + Ini1 + '_Ini_scen.nc') as TMP: Ini_scen = TMP.load()

        ## save (Ini_scen)
        if save is not None:
            Ini_scen.to_netcdf('results/' + save + '_Ini_scen.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Ini_scen})

        ## printing time counter
        print('total loading time (scenarios): {:.0f} seconds'.format(perf_counter() - t0))

        ##-----

        ## run (scenarios)
        Out_scen = model(Ini_scen, Par, For_scen, **model_args)
        
        ## save (Out_scen)
        if save is not None:
            Out_scen.to_netcdf('results/' + save + '_Out_scen.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Out_scen})
        
        ##-----

    ## return if requested
    if output:
        return {'Par':Par, 'For_hist':For_hist, 'Ini_hist':Ini_hist, 'Out_hist':Out_hist, 'For_scen':For_scen, 'Ini_scen':Ini_scen, 'Out_scen':Out_scen}

