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

"""
CONTENT
-------
1. FORCINGS
    1.1. Emissions
        load_emissions_hist
        load_emissions_scen
    1.2. Land-use
        load_landuse_hist
        load_landuse_scen
    1.3. Others
        load_RFdrivers_hist
        load_RFdrivers_scen
2. OBSERVATIONS
    2.1. Climate
        load_GMST TODO
        load_climOcean TODO
        load_climLand TODO
    2.2. Radiative forcing
        load_RF TODO
    2.3. Concentrations
        load_atmoConc TODO
    2.4. Carbon cycle
        load_GCB TODO
Z. WRAPPERS
    load_all_hist
    load_all_scen
"""

##################################################
##################################################

import os
import warnings
import numpy as np
import xarray as xr

from core_fct.fct_misc import aggreg_region, group_scenarios


##################################################
##   1. FORCINGS
##################################################

##===============
## 1.1. Emissions
##===============

## historical anthropogenic emissions
def load_emissions_hist(mod_region, 
    datasets=['ACCMIP', 'CDIAC', 'CEDS', 'EDGAR-HYDEv13', 'EDGAR-HYDEv14', 'EDGARv42', 'EDGARv42-FT2010', 'EDGARv432', 'EDGARv432-FT2016', 'EPA', 'Meinshausen_2011', 'PRIMAP'], 
    dropped_species=['PM10', 'PM2p5'], 
    **useless):
    '''
    Function to load and format primary historical emissions datasets, taken from the 'input_data' folder.
    
    Input:
    ------
    mod_region (str)            regional aggregation name

    Output:
    -------
    For (xr.Dataset)            dataset that contains the loaded datasets aggregated over 'mod_region'

    Options:
    --------
    datasets (list)             names of primary datasets to be loaded;
                                default = ['ACCMIP', 'CDIAC', 'CEDS', 'EDGAR-HYDEv13', 'EDGAR-HYDEv14', 
                                           'EDGARv42', 'EDGARv42-FT2010', 'EDGARv432', 'EDGARv432-FT2016', 
                                           'EPA', 'Meinshausen_2011', 'PRIMAP']
    dropped_species (list)      species to be excluded from the loaded datasets;
                                default = ['PM10', 'PM2p5']
    '''

    ## list of missing halogenated species taken from Meinshausen_2011
    missing_halo = ['CFC-11', 'CFC-12', 'CFC-113', 'CFC-114', 'CFC-115', 'CCl4', 'CH3CCl3', 'HCFC-22', 'Halon-1211', 'Halon-1202', 'Halon-1301', 'Halon-2402', 'CH3Br', 'CH3Cl']

    ## dictionaries for ignoring sectors
    ## CO2 emissions
    CO2_ignored_sectors = {
        'ACCMIP':['agr', 'awb', 'for', 'gra', 'wst'],
        'CEDS':['3B_Manure-management', '3D_Rice-Cultivation', '3D_Soil-emissions', '3E_Enteric-fermentation', '3I_Agriculture-other', '5A_Solid-waste-disposal', '5C_Waste-combustion', '5D_Wastewater-handling', '5E_Other-waste-handling'],
        'EDGAR-HYDEv13':['agr', 'bfc', 'liv', 'sav', 'def', 'awb', 'lfl'],
        'EDGAR-HYDEv14':['BF3', 'AGL', 'ANM', 'SAV', 'DEF', 'AGR', 'LAN'],
        'EDGARv42':['4A', '4B', '4C', '4D1', '4D2', '4D3', '4D4', '4E', '4F', '5A', '5C', '5D', '5F', '5F1', '5F2', '5FL', '5FL1', '6A', '6B', '6C', '6D'],
        'EDGARv42-FT2010':['5A', '5C', '5D', '5F2'],
        'EDGARv432':['4A', '4B', '4C', '4D1', '4D2', '4D3', '4D4', '4F', '6A', '6B', '6C', '6D'],
        'EDGARv432-FT2016':[],
        'EPA':['agr1', 'agr2', 'agr3', 'agr4', 'agr5', 'was1', 'was2', 'was3', 'was4'],
        'PRIMAP':['4', '5', '6']}
    ## non-CO2 emissions
    nonCO2_ignored_sectors = {
        'ACCMIP':['for', 'gra'],
        'CEDS':[],
        'EDGAR-HYDEv13':['sav', 'def'],
        'EDGAR-HYDEv14':['SAV', 'DEF'],
        'EDGARv42':['4E', '5A', '5C', '5D', '5F', '5F1', '5F2', '5FL', '5FL1'],
        'EDGARv42-FT2010':['4E', '5A', '5C', '5D', '5F2'],
        'EDGARv432':[],
        'EDGARv432-FT2016':[],
        'EPA':['agr5'], # slightly inconsistent as mixing e.g. agricultural waste burning & savannah burning
        'PRIMAP':['5']}

    ## main loading loop
    For0 = []
    units = {}
    for data in datasets:

        ## load data if available
        if os.path.isfile('input_data/drivers/emissions_' + data + '.nc'):
            with xr.open_dataset('input_data/drivers/emissions_' + data + '.nc') as TMP: 
                For1 = TMP.load()

        ## display message otherwise
        else: raise IOError('{0} not available'.format(data))

        ## get and check units
        for VAR in For1:
            if 'units' in For1[VAR].attrs:
                if VAR not in units.keys():
                    units[VAR] = For1[VAR].units
                else:
                    if units[VAR] != For1[VAR].units:
                        raise RuntimeWarning('inconsistent units: {0} (internal dic) vs. {1} ({2} in {3})'.format(units[VAR], For1[VAR].units, "'"+VAR+"'", "'"+data+"'"))

        ##---

        ## take only historical period
        if data in ['ACCMIP', 'Meinshausen_2011']:
            For1 = For1.sel(scen='historical', drop=True).dropna('year', how='all')

        ## take only national-based data
        if data in ['CDIAC']:
            For1 = For1.sel(data='national', drop=True)
        
        ## aggregate over fuels
        if data in ['CDIAC']:
            For1 = For1.sum('fuel', min_count=1)

        ## aggregate over sectors (ignoring some)
        if data in ['ACCMIP', 'CEDS', 'EDGAR-HYDEv13', 'EDGAR-HYDEv14', 'EDGARv42', 'EDGARv42-FT2010', 'EDGARv432', 'EDGARv432-FT2016', 'EPA', 'PRIMAP']:
            ## selection and aggregation
            for VAR in For1:
                if VAR == 'E_CO2': 
                    For1['E_CO2'] = For1['E_CO2'].sel(sector=[sec for sec in For1.sector.values if sec not in CO2_ignored_sectors[data]]).sum('sector', min_count=1)
                else: 
                    For1[VAR] = For1[VAR].sel(sector=[sec for sec in For1.sector.values if sec not in nonCO2_ignored_sectors[data]]).sum('sector', min_count=1)
            ## dropping useless coords
            For1 = For1.drop('sector')
            if 'sector_long_name' in For1.coords: 
                For1 = For1.drop('sector_long_name')

        ## take only missing halogenated species
        if data in ['Meinshausen_2011']:
            For1 = For1.drop([VAR for VAR in For1 if VAR != 'E_Xhalo'])
            For1 = For1.sel(spc_halo=missing_halo)

        ## put global data on regional axis
        if data in ['Meinshausen_2011']:
            For1 = For1.expand_dims('reg_iso', -1).assign_coords(reg_iso=[999])

        ##---

        ## interpolate to yearly data
        if not all(np.diff(For1.year.values) == 1):
            For1 = For1.interp({'year':np.arange(int(For1.year[0]), int(For1.year[-1])+1, 1)})

        ## rename CO2 emissions
        if 'E_CO2' in For1:
            For1 = For1.rename({'E_CO2':'Eff'})

        ## aggregate to model regions
        if 'reg_iso' in For1.coords:
            For1 = aggreg_region(For1, mod_region)
        else:
            For1 = aggreg_region(For1, mod_region, old_axis='region', dataset=data)

        ## append to final list (with new dimension)
        For0.append(For1.expand_dims('data', -1).assign_coords(data=[data]))
        del For1

    ## merge into one xarray
    For0 = xr.merge(For0)

    ## create one data axis per driver
    For = xr.Dataset()
    for VAR in For0:
        TMP = [For0[VAR].sel(data=data).rename({'data':'data_'+VAR}) for data in For0.data.values if not np.isnan(For0[VAR].sel(data=data).sum(min_count=1))]
        For[VAR] = xr.concat(TMP, dim='data_'+VAR)
        del TMP
    
    ## order dimensions
    For = For.transpose(*(['year', 'reg_land', 'spc_halo'] + [var for var in For.coords if 'data' in var]))

    ## drop requested species
    for vars_then_coords in [For, For.coords]:
        for VAR in vars_then_coords:
            if any([spc in VAR for spc in dropped_species]):
                For = For.drop(VAR)

    ## reapply units
    for VAR in For:
        if VAR in units.keys():
            For[VAR].attrs['units'] = units[VAR]

    ## return
    return For


#------
#------

## scenarios for anthropogenic emissions
def load_emissions_scen(mod_region, 
    datasets=['Meinshausen_2011', 'RCPdb', 'SRES', 'ScenarioMIP'],
    all_SRES=False, all_SSPdb=False, 
    dropped_species=['CCS'],
    Xhalo_offset={'CF4':{}, 'C2F6':{}, 'HFC-23':{}, 'CH3Br':{},
        'CH3Cl':{'RCPdb':3100.211, 'Meinshausen_2011':3100.211}},
    **useless):
    '''
    Function to load and format primary scenario emissions datasets, taken from the 'input_data' folder.
    
    Input:
    ------
    mod_region (str)            regional aggregation name

    Output:
    -------
    For (xr.Dataset)            dataset that contains the loaded datasets tentatively aggregated over 'mod_region'

    Options:
    --------
    datasets (list)             names of primary datasets to be loaded;
                                default = ['Meinshausen_2011', 'RCPdb', 'SRES', 'ScenarioMIP']
    all_SRES (bool)             whether to take all SRES scenarios (if loaded) or just markers;
                                default = False
    all_SSPdb (bool)            whether to take all SSP database scenarios (if loaded) or just markers;
                                default = False
    dropped_species (list)      species to be excluded from the loaded datasets;
                                default = ['CCS']
    Xhalo_offset (dict)         how the offset by Xhalo preindustrial emissions is handled;
                                keys are species whose emissions must be offset;
                                values are another dict being:
                                    either empty, in which case the offset is made on RCP2.6;
                                    or whose keys are dataset names and values are floats, for offset by the specified values;
                                default = {'CF4':{}, 'C2F6':{}, 'HFC-23':{}, 'CH3Br':{},
                                           'CH3Cl':{'RCPdb':3100.211, 'Meinshausen_2011':3100.211}}
    '''

    ## dictionaries for ignoring sectors
    ## non-CO2 emissions
    nonCO2_ignored_sectors = {
        'RCPdb':['for', 'gra'], 
        'ScenarioMIP':['Forest Burning', 'Grassland Burning', 'Peat Burning']}

    ## main loading loop
    For0 = []
    units = {}
    for data in datasets:

        ## load data if available
        if os.path.isfile('input_data/drivers/emissions_' + data + '.nc'):
            with xr.open_dataset('input_data/drivers/emissions_' + data + '.nc') as TMP: 
                For1 = TMP.load()

        ## display message otherwise
        else: raise IOError('{0} not available'.format(data))

        ## get and check units
        for VAR in For1:
            if 'units' in For1[VAR].attrs:
                if VAR not in units.keys():
                    units[VAR] = For1[VAR].units
                else:
                    if units[VAR] != For1[VAR].units:
                        raise RuntimeWarning('inconsistent units: {0} (internal dic) vs. {1} ({2} in {3})'.format(units[VAR], For1[VAR].units, "'"+VAR+"'", "'"+data+"'"))

        ##---

        ## take only halogenated species over right scenario
        ## note: all RCPs are supposed to be the same, however some inconsistent values appear but not in RCP2.6
        if data in ['Meinshausen_2011']:
            For1 = For1['E_Xhalo'].to_dataset(name='E_Xhalo')
            For1 = For1.sel(scen='RCP2.6').dropna('year', how='all')
            For1 = For1.assign_coords(scen='CMIP5').expand_dims('scen', -1)

        ## put global data on regional axis
        if data in ['Meinshausen_2011']:
            For1 = For1.expand_dims('reg_iso', -1).assign_coords(reg_iso=[999])

        ## offset with preindustrial emissions level
        if data in ['Meinshausen_2011', 'RCPdb']:
            for VAR in Xhalo_offset.keys():
                if data in Xhalo_offset[VAR].keys():
                    For1['E_Xhalo'] = xr.where(For1.spc_halo == VAR, For1['E_Xhalo'].sel(spc_halo=VAR) - Xhalo_offset[VAR][data], For1['E_Xhalo'])
                else:
                    if 'RCP2.6' in For1.scen:
                        For1['E_Xhalo'] = xr.where(For1.spc_halo == VAR, For1['E_Xhalo'].sel(spc_halo=VAR) - For1['E_Xhalo'].sel(spc_halo=VAR, scen='RCP2.6').dropna('year', how='all').isel(year=-1), For1['E_Xhalo'])
                    else:
                        For1['E_Xhalo'] = xr.where(For1.spc_halo == VAR, For1['E_Xhalo'].sel(spc_halo=VAR) - For1['E_Xhalo'].sel(spc_halo=VAR).dropna('year', how='all').isel(year=-1), For1['E_Xhalo'])

        ## aggregate over sectors (ignoring some)
        if data in ['RCPdb', 'ScenarioMIP']:
            ## selection and aggregation
            for VAR in For1:
                if 'sector' in For1[VAR].coords:
                    For1[VAR] = For1[VAR].sel(sector=[sec for sec in For1.sector.values if sec not in nonCO2_ignored_sectors[data]]).sum('sector', min_count=1)
            ## dropping useless coords
            For1 = For1.drop('sector')
            if 'sector_long_name' in For1.coords: 
                For1 = For1.drop('sector_long_name')

        ## select SRES scenarios
        if data in ['SRES']:
            ## take all scenarios (flatten array)
            if all_SRES:
                For1 = For1.stack(new_scen=('scen', 'model'))
                For1['new_scen'] = ['SRES-'+var1+' ('+var2+')' for var1, var2 in For1.new_scen.values]
                For1 = For1.dropna('new_scen').rename({'new_scen':'scen'})
            ## take only markers
            else: 
                For1 = For1.where(For1.is_marker, drop=True).sum('model', min_count=1)
                For1['scen'] = [data+'-'+var+' (marker)' for var in For1.scen.values]

        ## select SSP scenarios
        if data in ['SSPdb']:
            ## take all scenarios (flatten array)
            if all_SSPdb:
                For1 = For1.stack(new_scen=('scen_ssp', 'scen_rcp', 'model'))
                For1['new_scen'] = [var1+'-'+var2+' ('+var3+')' for var1, var2, var3 in For1.new_scen.values]
                For1 = For1.dropna('new_scen').rename({'new_scen':'scen'})
            ## take only markers (also flattened)
            else: 
                For1 = For1.where(For1.is_marker, drop=True).sum('model', min_count=1)
                For1 = For1.stack(new_scen=('scen_ssp', 'scen_rcp'))
                For1['new_scen'] = [var1+'-'+var2+' (marker)' for var1, var2 in For1.new_scen.values]
                For1 = For1.dropna('new_scen').rename({'new_scen':'scen'})

        ##---

        ## interpolate to yearly data
        if not all(np.diff(For1.year.values) == 1):
            For1 = For1.interp({'year':np.arange(int(For1.year[0]), int(For1.year[-1])+1, 1)})

        ## aggregate to model regions
        if 'reg_iso' in For1.coords:
            For1 = aggreg_region(For1, mod_region)
        else:
            For1 = aggreg_region(For1, mod_region, old_axis='region', dataset=data)

        ## append to final list
        For0.append(For1)
        del For1

    ## merge into one xarray
    For0 = xr.merge(For0)

    ## create one data axis per driver
    For = xr.Dataset()
    for VAR in For0:
        TMP = [For0[VAR].sel(scen=scen).rename({'scen':'scen_'+VAR}) for scen in For0.scen.values if not np.isnan(For0[VAR].sel(scen=scen).sum(min_count=1))]
        For[VAR] = xr.concat(TMP, dim='scen_'+VAR)
        del TMP
    
    ## order dimensions
    For = For.transpose(*(['year', 'reg_land'] + ['spc_halo']*('spc_halo' in For.coords) + [var for var in For.coords if 'scen' in var]))

    ## drop requested species
    for vars_then_coords in [For, For.coords]:
        for VAR in vars_then_coords:
            if any([spc in VAR for spc in dropped_species]):
                For = For.drop(VAR)

    ## reapply units
    for VAR in For:
        if VAR in units.keys():
            For[VAR].attrs['units'] = units[VAR]

    ## return
    return For


##==============
## 1.2. Land-use
##==============

## historical anthropogenic land-use and land-cover change
def load_landuse_hist(mod_region, 
    datasets=['LUH1', 'LUH1-TRENDYv4', 'LUH2', 'LUH2-TRENDYv8'],
    LCC='all',
    **useless):
    '''
    Function to load and format primary historical land-use datasets, taken from the 'input_data' folder.
    
    Input:
    ------
    mod_region (str)        regional aggregation name

    Output:
    -------
    For (xr.Dataset)        dataset that contains the loaded datasets aggregated over 'mod_region'

    Options:
    --------
    datasets (list)         names of primary datasets to be loaded;
                            default = ['LUH1', 'LUH1-TRENDYv4', 'LUH2', 'LUH2-TRENDYv8']
    LCC (str)               which of 'gross' or 'net' land-cover transitions should be kept;
                            unless both are kept ('all'), the driver is renamed to 'd_Acover';
                            default = 'all'
    '''

    ## main loading loop
    For0 = []
    units = {}
    for data in datasets:

        ## load data if available
        if os.path.isfile('input_data/drivers/land-use_' + data + '.nc'):
            with xr.open_dataset('input_data/drivers/land-use_' + data + '.nc') as TMP: 
                For1 = TMP.load()

        ## display message otherwise
        else: raise IOError('{0} not available'.format(data))

        ## get and check units
        for VAR in For1:
            if 'units' in For1[VAR].attrs:
                if VAR not in units.keys():
                    units[VAR] = For1[VAR].units
                else:
                    if units[VAR] != For1[VAR].units:
                        raise RuntimeWarning('inconsistent units: {0} (internal dic) vs. {1} ({2} in {3})'.format(units[VAR], For1[VAR].units, "'"+VAR+"'", "'"+data+"'"))

        ##---

        ## just assign data dimension
        if data in ['LUH1-TRENDYv4', 'LUH2-TRENDYv8']:
            For1 = For1.expand_dims('data_LULCC', -1).assign_coords(data_LULCC=[data])

        ## take historical only (+ data dim)
        if data in ['LUH1']:
            For1 = For1.sel(scen='historical', drop=True).dropna('year', how='all')
            For1 = For1.expand_dims('data_LULCC', -1).assign_coords(data_LULCC=[data])
        
        ## take historical variants (+ data dim)
        if data in ['LUH2']:
            For1 = For1.sel(scen=['historical', 'historical_high', 'historical_low']).dropna('year', how='all')
            For1 = For1.rename({'scen':'data_LULCC'}).assign_coords(data_LULCC=['LUH2', 'LUH2-High', 'LUH2-Low'])

        ##---

        ## interpolate to yearly data
        if not all(np.diff(For1.year.values) == 1):
            For1 = For1.interp({'year':np.arange(int(For1.year[0]), int(For1.year[-1])+1, 1)})

        ## aggregate to model regions
        if 'reg_iso' in For1.coords:
            For1 = aggreg_region(For1, mod_region)
        else:
            For1 = aggreg_region(For1, mod_region, old_axis='region', dataset=data)

        ## append to final list
        For0.append(For1)
        del For1

    ## merge into one xarray
    For0 = xr.merge(For0)

    ## order dimensions
    For0 = For0.transpose(*(['year', 'reg_land', 'bio_land', 'bio_from', 'bio_to'] + [var for var in For0.coords if 'data' in var]))

    ## reapply units
    for VAR in For0:
        if VAR in units.keys():
            For0[VAR].attrs['units'] = units[VAR]

    ## return (with selected net/gross LCC, if requested)
    if LCC == 'net': return For0.rename({'d_Anet':'d_Acover'}).drop('d_Agross')
    elif LCC == 'gross': return For0.rename({'d_Agross':'d_Acover'}).drop('d_Anet')
    else: return For0


#------
#------

## scenarios of anthropogenic land-use and land-cover change
def load_landuse_scen(mod_region, 
    datasets=['LUH1', 'LUH2'], 
    LCC='all',
    **useless):
    '''
    Function to load and format primary scenario land-use datasets, taken from the 'input_data' folder.
    
    Input:
    ------
    mod_region (str)        regional aggregation name

    Output:
    -------
    For (xr.Dataset)        dataset that contains the loaded datasets aggregated over 'mod_region'

    Options:
    --------
    datasets (list)         names of primary datasets to be loaded;
                            default = ['LUH1', 'LUH2']
    LCC (str)               which of 'gross' or 'net' land-cover transitions should be kept;
                            unless both are kept ('all'), the driver is renamed to 'd_Acover';
                            default = 'all'
    '''

    ## main loading loop
    For0 = []
    units = {}
    for data in datasets:

        ## load data if available
        if os.path.isfile('input_data/drivers/land-use_' + data + '.nc'):
            with xr.open_dataset('input_data/drivers/land-use_' + data + '.nc') as TMP: 
                For1 = TMP.load()

        ## display message otherwise
        else: raise IOError('{0} not available'.format(data))

        ## get and check units
        for VAR in For1:
            if 'units' in For1[VAR].attrs:
                if VAR not in units.keys():
                    units[VAR] = For1[VAR].units
                else:
                    if units[VAR] != For1[VAR].units:
                        raise RuntimeWarning('inconsistent units: {0} (internal dic) vs. {1} ({2} in {3})'.format(units[VAR], For1[VAR].units, "'"+VAR+"'", "'"+data+"'"))

        ##---

        ## take scenarios
        if data in ['LUH1', 'LUH2']:
            For1 = For1.sel(scen=[sc for sc in For1.scen.values if 'historical' not in sc]).dropna('year', how='all')
        
        ##---

        ## interpolate to yearly data
        if not all(np.diff(For1.year.values) == 1):
            For1 = For1.interp({'year':np.arange(int(For1.year[0]), int(For1.year[-1])+1, 1)})

        ## aggregate to model regions
        if 'reg_iso' in For1.coords:
            For1 = aggreg_region(For1, mod_region)
        else:
            For1 = aggreg_region(For1, mod_region, old_axis='region', dataset=data)

        ## append to final list
        For0.append(For1)
        del For1

    ## merge into one xarray (and rename scenario dimension)
    For0 = xr.merge(For0)
    For0 = For0.rename({'scen':'scen_LULCC'})

    ## order dimensions
    For0 = For0.transpose(*(['year', 'reg_land', 'bio_land', 'bio_from', 'bio_to'] + [var for var in For0.coords if 'scen' in var]))

    ## reapply units
    for VAR in For0:
        if VAR in units.keys():
            For0[VAR].attrs['units'] = units[VAR]

    ## return (with selected net/gross LCC, if requested)
    For0 = For0.drop('Aland_0')
    if LCC == 'net': return For0.rename({'d_Anet':'d_Acover'}).drop('d_Agross')
    elif LCC == 'gross': return For0.rename({'d_Agross':'d_Acover'}).drop('d_Anet')
    else: return For0


##============
## 1.3. Others
##============

## historical complementary RF drivers
def load_RFdrivers_hist(
    datasets=['radiative-forcing_AR5', 'volcanic-activity_CMIP6', 'solar-activity_CMIP6', 'aviation_ICAO'],
    extension={'RF_volc':'AOD550', 'RF_solar':'TSI', 'RF_contr':'dist_flown'}, 
    offset_volc=True,
    **useless):
    '''
    Function to load and format primary historical RF drivers datasets, taken from the 'input_data' folder.
    
    Input:
    ------
    None

    Output:
    -------
    For (xr.Dataset)        dataset that contains the loaded datasets

    Options:
    --------
    datasets (list)         names of primary datasets to be loaded;
                            default = ['radiative-forcing_AR5', 
                                       'volcanic-activity_CMIP6', 'solar-activity_CMIP6', 'aviation_ICAO']
    extension (dict)        keys are the name of RF drivers that will be extended;
                            values are the name of non-RF driversto be used for extension;
                            default = {'RF_volc':'AOD550', 'RF_solar':'TSI', 'RF_contr':'dist_flown'}
    offset_volc (bool)      whether volcanoes forcing should be offset by the average of the whole time-series;
                            default = True
    '''

    ## main loading loop
    For0 = []
    units = {}
    for data in datasets:

        ## load data if available
        if os.path.isfile('input_data/drivers/' + data + '.nc'):
            with xr.open_dataset('input_data/drivers/' + data + '.nc') as TMP: 
                For1 = TMP.load()
        elif os.path.isfile('input_data/observations/' + data + '.nc'):
            with xr.open_dataset('input_data/observations/' + data + '.nc') as TMP: 
                For1 = TMP.load()

        ## display message otherwise
        else: raise IOError('{0} not available'.format(data))

        ## get and check units
        for VAR in For1:
            if 'units' in For1[VAR].attrs:
                if VAR not in units.keys():
                    units[VAR] = For1[VAR].units
                else:
                    if units[VAR] != For1[VAR].units:
                        raise RuntimeWarning('inconsistent units: {0} (internal dic) vs. {1} ({2} in {3})'.format(units[VAR], For1[VAR].units, "'"+VAR+"'", "'"+data+"'"))

        ##---

        ## take only required drivers
        if data in ['radiative-forcing_AR5']:
            For1 = For1.drop([var for var in For1 if var not in extension.keys()])

        ## take detrended data
        if data in ['volcanic-activity_CMIP6']:
            For1 = For1.sel(data='detrended', drop=True)

        ## take historical data
        if data in ['solar-activity_CMIP6']:
            For1 = For1.sel(scen='historical', drop=True).dropna('year', how='all')

        ## offset by preindustrial value
        if data in ['solar-activity_CMIP6']:
            for VAR in For1:
                if VAR+'_0' in For1:
                    For1[VAR] -= For1[VAR+'_0']
                    For1 = For1.drop(VAR+'_0')

        ##---

        ## interpolate to yearly data
        if not all(np.diff(For1.year.values) == 1):
            For1 = For1.interp({'year':np.arange(int(For1.year[0]), int(For1.year[-1])+1, 1)})

        ## append to final list (with new dimension)
        For0.append(For1.expand_dims('data', -1).assign_coords(data=[data[::-1][:data[::-1].find('_')][::-1]]))
        del For1

    ## merge into one xarray
    For0 = xr.merge(For0)

    ## rescale extended driver to main dataset
    for VAR in extension.keys():
        for data in For0.data:

            ## catching warnings from theilslopes and empty slices averaging
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')

                if not np.isnan(For0[extension[VAR]].sel(data=data).sum(min_count=1)):
                    ## get values and overlapping years
                    ref = For0[VAR].mean('data')
                    ext = For0[extension[VAR]].sel(data=data)
                    years = ref.dropna('year').year & ext.dropna('year').year
                    ## apply rescale
                    For0[VAR] = xr.where(For0.data == data, ext * ref.sel(year=years).mean('year') / ext.sel(year=years).mean('year'), For0[VAR])

    ## offset volcanoes by average of total period (to simulate global average volcano)
    ## catching warnings from empty slices averaging
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        if offset_volc:
            For0['RF_volc'] -= For0['RF_volc'].mean('year')

    ## drop old variables
    For0 = For0.drop(list(extension.values()))

    ## create one data axis per driver
    For = xr.Dataset()
    for VAR in For0:
        TMP = [For0[VAR].sel(data=data).rename({'data':'data_'+VAR}) for data in For0.data.values if not np.isnan(For0[VAR].sel(data=data).sum(min_count=1))]
        For[VAR] = xr.concat(TMP, dim='data_'+VAR)
        del TMP

    ## order dimensions
    For = For.transpose(*(['year'] + [var for var in For.coords if 'data' in var]))

    ## reapply units
    for VAR in For:
        if VAR in units.keys():
            For[VAR].attrs['units'] = units[VAR]

    ## return
    return For

#------
#------

## scenarios of complementary RF drivers
## note: all defined in relative change; most of them hand-made!
def load_RFdrivers_scen(
    datasets=['solar-activity_CMIP6'],
    add_data=['volcanic-activity_CMIP5', 'volcanic-activity_CMIP6', 'solar-activity_CMIP5', 'contrails_CMIP5'],
    **useless):
    '''
    Function to load and format primary scenario RF drivers datasets, taken from the 'input_data' folder.
    WARNING: these are used to extend historical drivers, and are defined as relative variations!
    
    Input:
    ------
    None

    Output:
    -------
    For (xr.Dataset)        dataset that contains the loaded datasets

    Options:
    --------
    datasets (list)         names of primary datasets to be loaded;
                            default = ['solar-activity_CMIP6']
    add_data (list)         names of datasets that are hard-coded in the function!
                            default = ['volcanic-activity_CMIP5', 'volcanic-activity_CMIP6', 'solar-activity_CMIP5', 'contrails_CMIP5']
    '''

    ## main loading loop
    For0 = []
    for data in datasets:

        ## load data if available
        if os.path.isfile('input_data/drivers/' + data + '.nc'):
            with xr.open_dataset('input_data/drivers/' + data + '.nc') as TMP: 
                For1 = TMP.load()

        ## display message otherwise
        else: raise IOError('{0} not available'.format(data))

        ##---

        ## take historical data
        if data in ['solar-activity_CMIP6']:
            For1 = For1.sel(scen=['historical', 'proj_CMIP6']).mean('scen').sel(year=slice(2014-11+1, None, None)).dropna('year', how='all')

        ## offset by preindustrial value
        if data in ['solar-activity_CMIP6']:
            for VAR in For1:
                if VAR+'_0' in For1:
                    For1[VAR] -= For1[VAR+'_0']
                    For1 = For1.drop(VAR+'_0')
        
        ## rescale to reference (last 11 years) to be used as relative variation
        if data in ['solar-activity_CMIP6']:
            For1 /= For1.sel(year=slice(2014-11+1, 2014)).mean('year')
            For1 = For1.rename({'TSI':'RF_solar'})
            For1['RF_solar'].attrs['units'] = '1'

        ##---

        ## interpolate to yearly data
        if not all(np.diff(For1.year.values) == 1):
            For1 = For1.interp({'year':np.arange(int(For1.year[0]), int(For1.year[-1])+1, 1)})

        ## append to final list (with new dimension)
        For0.append(For1.expand_dims('scen', -1).assign_coords(scen=[data[::-1][:data[::-1].find('_')][::-1]]))
        del For1

    ## secondary loading loop for hand-made data
    for data in add_data:

        ## nil extension
        if data in ['volcanic-activity_CMIP5']:
            For1 = xr.DataArray(np.zeros([2500-2005+1]), coords={'year':np.arange(2005, 2500+1)}, dims='year', attrs={'units':'1'}).to_dataset(name='RF_volc')

        ## constant extension
        if data in ['solar-activity_CMIP5']:
            For1 = xr.DataArray(np.ones([2500-2005+11]), coords={'year':np.arange(2005-11+1, 2500+1)}, dims='year', attrs={'units':'1'}).to_dataset(name='RF_contr')

        if data in ['contrails_CMIP5']:
            For1 = xr.DataArray(np.ones([2500-2005+1]), coords={'year':np.arange(2005, 2500+1)}, dims='year', attrs={'units':'1'}).to_dataset(name='RF_contr')

        ## 10-year linear decrease
        if data in ['volcanic-activity_CMIP6']:
            For1 = xr.DataArray(np.append(np.arange(1, 0, -0.1), np.zeros([2300-2024+1])), coords={'year':np.arange(2014, 2300+1)}, dims='year', attrs={'units':'1'}).to_dataset(name='RF_volc')

        ## append to final list (with new dimension)
        For0.append(For1.expand_dims('scen', -1).assign_coords(scen=[data[::-1][:data[::-1].find('_')][::-1]]))
        del For1

    ## merge into one xarray
    For0 = xr.merge(For0)

    ## create one data axis per driver
    For = xr.Dataset()
    for VAR in For0:
        TMP = [For0[VAR].sel(scen=scen).rename({'scen':'scen_'+VAR}) for scen in For0.scen.values if not np.isnan(For0[VAR].sel(scen=scen).sum(min_count=1))]
        For[VAR] = xr.concat(TMP, dim='scen_'+VAR)
        del TMP

    ## order dimensions
    For = For.transpose(*(['year'] + [var for var in For.coords if 'scen' in var]))

    ## return
    return For


##################################################
##   2. OBSERVATIONS
##################################################

##=============
## 2.1. Climate
##=============



##################################################
##   Z. WRAPPERS
##################################################

## wrapping function for historical
def load_all_hist(mod_region, **kwargs):
    '''
    Wrapper function to load all primary historical drivers.
    
    Input:
    ------
    mod_region (str)        regional aggregation name       

    Output:
    -------
    For (xr.Dataset)        merged dataset

    Options:
    --------
    **kwargs                arguments to be passed on to individual loading functions
    '''

    print('loading primary historical drivers')

    ## list of loading fuctions
    load_list = [load_emissions_hist, load_landuse_hist, load_RFdrivers_hist]
    
    ## return all
    return xr.merge([load(mod_region=mod_region, **kwargs) for load in load_list])


## wrapping function for scenarios
def load_all_scen(mod_region, group_scen=None, default_scen=['CMIP6', 'CMIP5'], **kwargs):
    '''
    Wrapper function to load all primary scenario drivers.
    
    Input:
    ------
    mod_region (str)        regional aggregation name       

    Output:
    -------
    For (xr.Dataset)        merged dataset

    Options:
    --------
    group_scen (list)       scenarios that must be grouped under the same 'scen' dimension;
                            default = None
    default_scen (list)     fallback scenarios taken in order when using 'group_scen' option;
                            default = ['CMIP6', 'CMIP5']
    **kwargs                arguments to be passed on to individual loading functions
    '''
    
    print('loading primary scenario drivers')

    ## list of loading fuctions
    load_list = [load_emissions_scen, load_landuse_scen, load_RFdrivers_scen]
    
    ## get merged data
    ds = xr.merge([load(mod_region=mod_region, **kwargs) for load in load_list])

    ## group scenarios if requested
    if group_scen is not None: 
        ds = group_scenarios(ds, group_scen=group_scen, default_scen=default_scen)
        ds = ds.sel(scen=[var for var in group_scen if var in ds.scen])

    ## return dataset
    return ds

