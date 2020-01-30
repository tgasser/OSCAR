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

import random
import warnings
import numpy as np
import xarray as xr

from scipy.stats import theilslopes

from core_fct.fct_misc import extend_timeseries


##################################################
## 1. GENERATE MONTE CARLO PARAMETERS
##################################################

## generate all Monte Carlo configurations 
def generate_config(Par, nMC, ignored=['mean_', 'old_', 'off_'], fixed={}, return_details=False):
    '''
    Function to generate Monte Carlo configuration (= parameters) for OSCAR.
    
    Input:
    ------
    Par (xr.Dataset)        dataset containing primary parameters
    nMC (int)               number of MC elements
    
    Output:
    -------
    Par_mc (xr.Dataset)     dataset containing MC parameters
    mod_mc (xr.Dataset)     (optional) dataset containing the name of the chosen options for each configuration

    Options:
    --------
    ignored (list)          configurations containing any of the list's elements in their names will be ignored;
                            default = ['mean_', 'old_', 'off_']
    fixed (dict)            dict of options to be fixed (i.e. ignored by MC and taken as provided);
                            keys are configuration options (starting with 'mod_') and values are desired configuration;
                            default = {}
    return_details (bool)   whether details on each configuration should be returned
                            default = False
    '''

    print('generating MC configurations')

    ## get list of mod_ options and other dimensions
    mod_list = [var for var in Par.coords if var[:4] == 'mod_']

    ## draw configurations (loop)
    Par_mc = []
    mod_mc = []
    for n in range(nMC):

        ## pick randomly
        mod_dic = {mod:random.choice([var for var in Par[mod].values if all([test not in var for test in ignored])]) for mod in mod_list}
        
        ## apply fixed values and select config
        for mod in fixed.keys(): mod_dic[mod] = fixed[mod]
        Par_mc.append(Par.sel(mod_dic, drop=True).assign_coords(config=n).expand_dims('config', -1))
    
        ## save config details
        mod_mc.append(xr.Dataset(mod_dic).assign_coords(config=n).expand_dims('config', -1))

    ## concatenate configurations
    Par_mc = xr.concat(Par_mc, dim='config')
    mod_mc = xr.concat(mod_mc, dim='config')

    ## remove useless configuration axis
    for VAR in Par_mc:
        if np.all(Par_mc[VAR].sel(config=0, drop=True) == Par_mc[VAR]):
            Par_mc[VAR] = Par_mc[VAR].sel(config=0, drop=True)

    ## return final dataset
    if not return_details: return Par_mc
    else: return Par_mc, mod_mc


##################################################
## 2. GENERATE MONTE CARLO FORCINGS
##################################################

## generate all Monte Carlo drivers
def generate_drivers(For, nMC, fixed={}):
    '''
    Function to geenrate Monte Carlo (historical) drivers for OSCAR.
    
    Input:
    ------
    For (xr.Dataset)        dataset containing secondary historical drivers
    nMC (int)               number of MC elements
    
    Output:
    -------
    For_mc (xr.Dataset)     dataset containing MC historical drivers
    Options:
    --------
    fixed (dict)            dict of data options to be fixed (i.e. ignored by MC and taken as provided);
                            keys are data options (starting with 'data_') and values are desired configuration;
                            default = {}
    '''

    print('generating MC drivers')

    ## get list of data_ options and other dimensions
    data_list = [var for var in For.coords if var[:5] == 'data_']

    ## draw configurations (loop)
    For_mc = []
    for n in range(nMC):
        
        ## pick randomly
        data_dic = {}
        for data in data_list:
            ## masked data (typically, because of subspecies with different sets of data like with halogenated)
            if data.replace('data_', 'mask_') in For.coords:
                list_dims = [var for var in For[data.replace('data_', 'mask_')].dims if var != data]
                if len(list_dims) == 1: other_dim = list_dims[0]
                else: raise RuntimeError("mask of '{0}' for MC generation should have only one additional dimension: {1}".format(data, list_dims))
                data_random = [random.choice(For[data].where(For[data.replace('data_', 'mask_')]).sel(**{other_dim:val}).dropna(data).values) for val in For[other_dim]]
                data_dic[data] = xr.DataArray(data_random, coords={other_dim:For[other_dim]}, dims=other_dim)
            ## regular data
            else:
                data_dic[data] = random.choice([var for var in For[data].values])
        
        ## apply fixed values and select config
        for data in fixed.keys(): data_dic[data] = fixed[data]
        For_mc.append(For.sel(data_dic, drop=True).expand_dims('config', -1).assign_coords(config=[n]))
    
    ## concatenate configurations
    For_mc = xr.concat(For_mc, dim='config')

    ## drop masked coordinates
    For_mc = For_mc.drop([var for var in For_mc.coords if var[:5] == 'mask_'])

    ## remove useless configuration axis
    for VAR in For_mc:
        if np.all(For_mc[VAR].sel(config=0, drop=True) == For_mc[VAR]):
            For_mc[VAR] = For_mc[VAR].sel(config=0, drop=True)

    ## return final dataset
    return For_mc

