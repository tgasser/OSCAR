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
## 1. GENERATE ONE TIME-SERIES
##################################################

## create one full timeseries of historical drivers
def make_one_timeseries(Old, data_list, inds, data_connect, missing_data, preindustrial, ref_length=None):
    '''
    Function to put together one time-series of forcing data.
    
    Input:
    ------
    Old (xr.Dataset)        dataset containing all data sources of data_list
    data_list (list)        list of data sources in chronological order;
                            name of the master (central) data source must be surrounded with '*'
    inds (tuple)            starting (ind0), preindustrial (indPI) and ending (indH) years of the output time-series
    data_connect (str)      how data sources are connected together:
                            'raw' = data sources are simply juxtaposed;
                            'rel_change' = data sources are extended following relative variations
    missing_data (str)      how the time-series is extended forward for missing data:
                            'cst' = constant extension;
                            'trend' = extension using linear trend;
                            'safe_trend' = extension using linear trend, but negative values are set to zero
    preindustrial (str)     how preindustrial (e.g. emissions) are defined:
                            'zero' = absolute values of the forcing are kept (thus preindustrial value is zero by definition);
                            'offset' = the whole time-series is offset by its value at indPI (this offsetting value is the preindsutrial value);
                            'safe_offset' = the whole time-series is offset by its value at indPI, but negative values are set to zero
        
    Output:
    -------
    New (xr.Dataset)        dataset which consistent time-series covering ind0 to indH

    Options:
    --------
    ref_length (int)        length of time period used for extension of missing_data;
                            default = None
    '''

    ## check arguments are a possibility
    assert data_connect in ['raw', 'rel_change']
    assert missing_data in ['cst', 'trend', 'safe_trend']
    assert preindustrial in ['zero', 'offset', 'safe_offset']

    ## unpack indexes
    ind0, indPI, indH = inds[:3]

    ## get data info
    data_axis = [coord for coord in Old.dims if 'data_' in coord][0]

    ## get initial dataset
    data_0 = data_list[[var.count('*') for var in data_list].index(2)]
    New = Old.sel({data_axis:data_0.replace('*', '')}, drop=True).dropna('year', how='all')

    ## complement forward
    for data in data_list[data_list.index(data_0)+1:]:
        if data_connect == 'rel_change' and 'reg_land' in New.coords:
            New = xr.merge([extend_timeseries(New[VAR], Old[VAR].sel({data_axis:data}, drop=True), 'forward', scale_to_global=True, dump_loss_in={'reg_land':0}).to_dataset(name=VAR) for VAR in New])
        elif data_connect == 'rel_change':
            New = xr.merge([extend_timeseries(New[VAR], Old[VAR].sel({data_axis:data}, drop=True), 'forward', scale_to_global=False).to_dataset(name=VAR) for VAR in New])
        elif data_connect == 'raw':
            New = xr.merge([extend_timeseries(New[VAR], Old[VAR].sel({data_axis:data}, drop=True), 'forward', juxtaposition=True).to_dataset(name=VAR) for VAR in New])

    ## randomly complement backward
    for data in data_list[:data_list.index(data_0)][::-1]:
        if data_connect == 'rel_change' and 'reg_land' in New.coords:
            New = xr.merge([extend_timeseries(New[VAR], Old[VAR].sel({data_axis:data}, drop=True), 'backward', scale_to_global=True, dump_loss_in={'reg_land':0}).to_dataset(name=VAR) for VAR in New])
        elif data_connect == 'rel_change':
            New = xr.merge([extend_timeseries(New[VAR], Old[VAR].sel({data_axis:data}, drop=True), 'backward', scale_to_global=False).to_dataset(name=VAR) for VAR in New])
        elif data_connect == 'raw':
            New = xr.merge([extend_timeseries(New[VAR], Old[VAR].sel({data_axis:data}, drop=True), 'backward', juxtaposition=True).to_dataset(name=VAR) for VAR in New])

    ## extend further forward
    if New.year.max() < indH:
        years = xr.DataArray(np.arange(New.year.max()-ref_length+1, indH+1), coords={'year':np.arange(New.year.max()-ref_length+1, indH+1)}, dims=['year'])
        for VAR in New:        
        
            ## catching warnings from theilslopes and empty slices averaging
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')

                ## extend following latest trend
                if missing_data in ['trend', 'safe_trend']:

                    ## reference year
                    VAR_ref = New[VAR].loc[New.year.max()-ref_length+1:New.year.max()].mean('year')
                    year_ref = New[VAR].loc[New.year.max()-ref_length+1:New.year.max()].year.mean('year')
                    
                    ## local and global trends
                    trend_reg = xr.apply_ufunc(lambda X: theilslopes(X)[0], New[VAR].loc[New.year.max()-ref_length+1:New.year.max()], input_core_dims=[['year']], vectorize=True)
                    if 'reg_land' in New.coords:
                        trend_glob = theilslopes(New[VAR].loc[New.year.max()-ref_length+1:New.year.max()].sum('reg_land', min_count=1))[0]
                    
                    ## extension (remove last line if not scale_to_global)
                    Ext = (years - year_ref) * trend_reg + VAR_ref
                    if 'reg_land' in New.coords:
                        Ext.loc[{'reg_land':Ext.reg_land==0}] = np.nan_to_num(Ext.loc[{'reg_land':0}]) + (years - year_ref) * trend_glob + VAR_ref.sum('reg_land', min_count=1) - Ext.sum('reg_land', min_count=1)
                    
                    ## ensure no negative emissions (if requested)
                    if missing_data == 'safe_trend': 
                        Ext = Ext.where(Ext >= 0, other=0)
                
                ## or with constant values
                elif missing_data == 'cst':
                    Ext = 0.*years + New[VAR].loc[New.year.max()-ref_length+1:New.year.max()].mean('year')

                ## apply extension to existing data
                New = New.combine_first(Ext.to_dataset(name=VAR))

    ## extend further backward linearly to zero
    if New.year.min() > indPI:
        New = New.combine_first(xr.zeros_like(New.isel(year=0)).assign_coords(year=indPI).expand_dims('year', -1))
        New = New.combine_first(xr.merge([New[VAR].interp({'year':np.arange(indPI, indH+1, 1)}).to_dataset(name=VAR) for VAR in New]))

    ## slice over proper time period
    if New.year.min() <= indPI:
        New = New.sel(year=slice(indPI,indH))
        
        ## offset by preindustrial value (if requested)
        if preindustrial in ['offset', 'safe_offset']:
            New = New - New.sel(year=indPI, drop=True).fillna(0)
            
            ## ensure no negative emissions (if requested)
            if preindustrial == 'safe_offset':
                New = New.where(New >= 0, other=0)

    ## extend with zeros over preindustrial period
    if ind0 < indPI:
        New = New.combine_first(0.*xr.DataArray(np.arange(ind0, indPI), coords={'year':np.arange(ind0, indPI)}, dims=['year']) + xr.zeros_like(New.sum('year')))

    ## return
    return New.fillna(0.).assign_coords(**{data_axis:' '.join(data_list)}).expand_dims(data_axis, -1)


##################################################
## 2. GENERATE ALL TIME-SERIES
##################################################

##=====================
## 2.1. Check all cases
##=====================

## check all possible combinations
## /!\ WARNING: this is done randomly! it is inefficient and unsuited for too many combinations (>1000)
def check_combinations(For, inds, ignored=[], n_safe=2000, n_break=1000):
    '''
    Function to check all possible combinations of data sources to create consistent time-series of forcing data.
    
    Input:
    ------
    For (xr.Dataset)    dataset containing all desired data sources
    inds (tuple)        starting (ind0), preindustrial (indPI) and ending (indH) years of the output time-series

    Output:
    -------
    combi (list)        list of all possible combinations found;
                        each combination is given as a chronological list of data sources, with central one surrounded by '*'

    Options:
    --------
    ignored (list)      name of data sources to be ignored as master (central) data sources;
                        default = []
    n_safe (int)        random loop will assume convergence after [number of combinations found] * n_safe iterations;
                        default = 2000
    n_break (int)       random loop will break after n_safe * n_break iterations;
                        default = 1000
    '''

    ## unpack indexes
    indPI, indH = inds[1:3]

    ## get data info
    data_axis = [coord for coord in For.dims if 'data_' in coord][0]
    data_values = For.dropna(data_axis, how='all').coords[data_axis].values

    ## get years covered by emission datasets
    year_min, year_max = {}, {}
    for data in data_values:
        year_min[data] = int(For.sel({data_axis:data}).dropna('year', how='all').year.min())
        year_max[data] = int(For.sel({data_axis:data}).dropna('year', how='all').year.max())

    ## look for all combinations randomly!
    combi = set()
    n_loop = 0
    while n_loop <= len(combi) * n_safe:
        n_loop += 1

        ## set initial dataset
        data = random.choice([var for var in data_values if var not in ignored])
        data_list = ['*'+data+'*']
        years = [year_min[data], year_max[data]]

        ## randomly complement forward
        while max(years) < indH and any([year_max[var] > max(years) and year_min[var] < max(years) for var in data_values]):
            data = random.choice([var for var in data_values if year_max[var] > max(years) and year_min[var] < max(years)])
            data_list.append(data)
            years[-1] = year_max[data]

        ## randomly complement backward
        while min(years) > indPI and any([year_min[var] < min(years) and year_max[var] > min(years) for var in data_values]):
            data = random.choice([var for var in data_values if year_min[var] < min(years) and year_max[var] > min(years)])
            data_list.insert(0, data)
            years[0] = year_min[data]

        ## add to pool of combinations
        combi = combi | set([' '.join(data_list)])

        ## safe breaking of while loop
        if n_loop > n_safe*n_break:
            break

    ## return
    return combi


##======================
## 2.2. Create all cases
##======================

## create all possible combinations for each driver
def create_hist_drivers(For, inds, ignored=['EDGAR-HYDEv13', 'EDGAR-HYDEv14'], 
    Xhalo_PI={'CF4':1950, 'C2F6':1950, 'HFC-23':1950, 'CH3Br':1950, 'CH3Cl':1950},
    LCC_vars=['d_Acover', 'd_Anet', 'd_Agross']):

    '''
    Function to create historical drivers for OSCAR.
    
    Input:
    ------
    For (xr.Dataset)    dataset containing primary drivers
    inds (tuple)        starting (ind0), preindustrial (indPI) and ending (indH) years for the drivers

    Output:
    -------
    For0 (xr.Dataset)   dataset containing secondary drivers

    Options:
    --------
    ignored (list)      name of data sources to be ignored as master (central) data sources;
                        default = ['EDGAR-HYDEv13', 'EDGAR-HYDEv14']
    Xhalo_PI (dict)     keys are halogenated species and values are years to be considered preindustrial;
                        default = {'CF4':1950, 'C2F6':1950, 'HFC-23':1950, 'CH3Br':1950, 'CH3Cl':1950}
    LCC_vars (list)     list of variables that describe land-cover change;
                        this will be used to calculate preindustrial land-cover;
                        default = ['d_Acover', 'd_Anet', 'd_Agross']
    '''

    print('creating secondary historical drivers')

    ## get units
    units = {VAR:For[VAR].units for VAR in For if 'units' in For[VAR].attrs}

    ## loop on data axes
    For0 = []
    for data_axis in [var for var in For.dims if 'data_' in var]:

        ## select subset of variables defined on data_axis
        For1 = For.drop([VAR for VAR in For if data_axis not in For[VAR].dims])
        For1 = For1.drop([var for var in For1.coords if all([var not in For1[VAR].coords for VAR in For1])])

        ## actually create combinations
        ## 1. emissions (fossil fuel)
        if 'data_Eff' in data_axis:
            options = {'inds':inds, 'data_connect':'rel_change', 'missing_data':'trend', 'preindustrial':'zero', 'ref_length':5}
            data_combi = [var.split(' ') for var in check_combinations(For1, inds=inds, ignored=ignored)]
            For2 = [make_one_timeseries(For1, data_list, **options) for data_list in data_combi]
            For0.append(xr.concat(For2, dim=data_axis))

        ## 2. emissions (halogenated)
        elif 'data_E_Xhalo' in data_axis:
            For_tmp = []
            for spc in For1.spc_halo.values:
                spc_inds = (inds[0], max(inds[1], Xhalo_PI[spc]) if str(spc) in Xhalo_PI.keys() else inds[1], inds[2])
                options = {'inds':spc_inds, 'data_connect':'rel_change', 'missing_data':'cst', 'preindustrial':'offset', 'ref_length':1}
                data_combi = [var.split(' ') for var in check_combinations(For1.sel(spc_halo=spc, drop=True), inds=inds, ignored=ignored)]
                For2 = [make_one_timeseries(For1.sel(spc_halo=spc, drop=True), data_list, **options) for data_list in data_combi]
                For_tmp.append(xr.concat(For2, dim=data_axis).expand_dims('spc_halo', -2).assign_coords(spc_halo=[spc]))
            For_tmp = xr.concat(For_tmp, dim='spc_halo')
            For_tmp.coords[data_axis.replace('data_','mask_')] = np.logical_not(For_tmp['E_Xhalo'].isnull().all([var for var in For_tmp.dims if var not in ['spc_halo', data_axis]]))
            For0.append(For_tmp)            

        ## 3. emissions (rest)
        elif 'data_E' in data_axis:
            options = {'inds':inds, 'data_connect':'rel_change', 'missing_data':'trend', 'preindustrial':'offset', 'ref_length':5}
            data_combi = [var.split(' ') for var in check_combinations(For1, inds=inds, ignored=ignored)]
            For2 = [make_one_timeseries(For1, data_list, **options) for data_list in data_combi]
            For0.append(xr.concat(For2, dim=data_axis))

        ## 4. land-use
        ## note: options other than data_connect='raw' and missing_data='cst' will likely ruin data consistency!
        elif 'data_LULCC' in data_axis:
            landuse_inds = (inds[0], inds[0], inds[2])
            options = {'inds':landuse_inds, 'data_connect':'raw', 'missing_data':'cst', 'preindustrial':'zero', 'ref_length':5}
            if len(LCC_vars) > 0:
                VAR = [var for var in For1 if var in LCC_vars][0]
                For1['Aland'] = For1[VAR].sum('bio_from', min_count=1).rename({'bio_to':'bio_land'}).cumsum('year') - For1[VAR].sum('bio_to', min_count=1).cumsum('year').rename({'bio_from':'bio_land'})
                For1['Aland'] += For1.Aland_0 - For1['Aland'].sel(year=For1.Aland_0.year,drop=True)
            data_combi = [var.split(' ') for var in check_combinations(For1, inds=inds, ignored=ignored)]
            For2 = xr.concat([make_one_timeseries(For1.drop('Aland_0'), data_list, **options) for data_list in data_combi], dim=data_axis)
            if len(LCC_vars) > 0:
                For2['Aland_0'] = For2.Aland.sel(year=inds[0])
                For2 = For2.drop('Aland')
            For0.append(For2)

        ## 5. RF drivers (natural)
        elif 'data_RF_volc' in data_axis or 'data_RF_solar' in data_axis:
            options = {'inds':inds, 'data_connect':'rel_change', 'missing_data':'cst', 'preindustrial':'zero', 'ref_length':11}
            data_combi = [var.split(' ') for var in check_combinations(For1, inds=inds, ignored=ignored)]
            For2 = [make_one_timeseries(For1, data_list, **options) for data_list in data_combi]
            For0.append(xr.concat(For2, dim=data_axis))

        ## 6. RF drivers (anthropogenic)
        elif 'data_RF' in data_axis:
            options = {'inds':inds, 'data_connect':'rel_change', 'missing_data':'trend', 'preindustrial':'zero', 'ref_length':5}
            data_combi = [var.split(' ') for var in check_combinations(For1, inds=inds, ignored=ignored)]
            For2 = [make_one_timeseries(For1, data_list, **options) for data_list in data_combi]
            For0.append(xr.concat(For2, dim=data_axis))

    ## merge all variables
    For0 = xr.merge(For0)

    ## reapply units
    for VAR in For0:
        if VAR in units.keys():
            For0[VAR].attrs['units'] = units[VAR]

    ## return final dataset
    return For0


##################################################
## 3. CONNECT HISTORICAL AND SCENARIOS
##################################################

## combine historical and scenarios timeseries
def create_scen_drivers(Hist, Scen, inds, data_connect='transition', trans_length=20,
    raw_vars=['d_Acover', 'd_Anet', 'd_Agross', 'd_Ashift', 'd_Hwood'],
    rel_vars={'RF_contr':1, 'RF_solar':11, 'RF_volc':1},
    copy_last_year=True):
    
    '''
    Function to create future drivers for OSCAR that are consistent with historical ones.
    
    Input:
    ------
    Hist (xr.Dataset)       dataset containing secondary or MC historical drivers
    Scen (xr.Dataset)       dataset containing primary scenario drivers
    inds (tuple)            starting (ind0), preindustrial (indPI), junction (indH) and ending (indF) years for the drivers

    Output:
    -------
    For (xr.Dataset)        dataset containing secondary scenario drivers

    Options:
    --------
    data_connect (str)      how historical and scenarios are connected together:
                            'raw' = they are simply juxtaposed;
                            'rel_change' = scenarios used to extend historical following their relative variations;
                            'transition' = move linearly from 'rel_change' to 'raw' in trans_length time steps;
                            default = 'transition'
    trans_length (int)      number of time steps used to transition in data_connect;
                            default = 20
    raw_vars (list)         list of variables that are forced to use 'raw' data_connect;
                            default = ['d_Acover', 'd_Anet', 'd_Agross', 'd_Ashift', 'd_Hwood']
    rel_vars (dict)         keys are variables that are forced to use 'rel_change' data_connect;
                            values are the length of reference period for rescaling;
                            rel_vars = {'RF_contr':1, 'RF_solar':11, 'RF_volc':1}
    copy_last_year (bool)   if True, if only one year is missing at the end to match indF, the second to last year is replicated
    '''

    print('creating secondary scenario drivers')

    ## check arguments are a possibility
    assert data_connect in ['raw', 'rel_change', 'transition']

    ## unpack indexes
    indH, indF = inds[-2:]

    ## get units
    units = {VAR:Hist[VAR].units for VAR in Hist if 'units' in Hist[VAR].attrs}

    ## create combined variables
    For = xr.Dataset()
    for VAR in Hist:
        if VAR in Scen:

            ## stable axis for extension
            stable_axis = ['config', 'scen'] + ['spc_halo']*('spc_halo' in Scen[VAR].coords)

            ## raw juxtaposition
            if VAR in raw_vars or data_connect == 'raw':
                For[VAR] = extend_timeseries(Hist[VAR], Scen[VAR], 'forward', juxtaposition=True, stable_axis=stable_axis)

            ## use scenario for relative change (= rescaled)
            elif VAR in rel_vars.keys() or data_connect == 'rel_change':
                if VAR not in rel_vars.keys(): ref_length = 1
                else: ref_length = rel_vars[VAR]
                if 'reg_land' in Scen[VAR].coords:
                    For[VAR] = extend_timeseries(Hist[VAR], Scen[VAR], 'forward', scale_to_global=True, ref_length=ref_length, dump_loss_in={'reg_land':0}, stable_axis=stable_axis)
                else: 
                    For[VAR] = extend_timeseries(Hist[VAR], Scen[VAR], 'forward', scale_to_global=False, ref_length=ref_length, stable_axis=stable_axis)

            ## transition from rescaled to raw data
            elif data_connect == 'transition':
                ## get raw and rel_change data
                raw = extend_timeseries(Hist[VAR], Scen[VAR], 'forward', juxtaposition=True, stable_axis=stable_axis)
                if 'reg_land' in Scen[VAR].coords: 
                    rel = extend_timeseries(Hist[VAR], Scen[VAR], 'forward', scale_to_global=True, dump_loss_in={'reg_land':0}, stable_axis=stable_axis)
                else: 
                    rel = extend_timeseries(Hist[VAR], Scen[VAR], 'forward', scale_to_global=False, stable_axis=stable_axis)                    
                ## transition over trans_length
                tt = xr.concat([1 + 0.*raw.year.where(raw.year <= indH, drop=True), 0.*Scen.year.where(raw.year >= indH + trans_length, drop=True)], dim='year')
                tt = tt.combine_first(tt.interp({'year':np.arange(indH, indH + trans_length, 1)}))
                For[VAR] = raw * (1-tt) + tt * rel

        ## if missing variable
        else: For[VAR] = 0.*Scen.year

    ## slice over right period
    For = For.sel(year=slice(indH, indF, None))

    ## set last year equal to previous year if missing
    if copy_last_year:
        for VAR in For:
            if np.all(For[VAR].isel(year=-1).isnull()):
                For[VAR][dict(year=-1)] = For[VAR].isel(year=-2)

    ## reapply units
    for VAR in For:
        if VAR in units.keys():
            For[VAR].attrs['units'] = units[VAR]

    ## return final dataset
    return For.fillna(0.)

