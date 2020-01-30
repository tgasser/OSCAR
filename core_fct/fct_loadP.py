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
1. CARBON CYCLE
    1.1. Ocean
        load_ocean_struct
        load_ocean_chem
        load_ocean_CMIP5
    1.2. Land
        load_land_misc
        load_land_TRENDYv7
        load_land_CMIP5
        load_land_wooduse
        load_land_GFED3
    1.3. Permafrost
        load_permafrost_all
    1.4. Wetlands
        load_wetlands_WETCHIMP
2. ATMO CHEMISTRY
    2.1. Atmosphere
        load_atmosphere_misc
        load_atmosphere_CCMVal2
        load_atmochem_adhoc
    2.2. Methane
        load_CH4_lifetime
        load_OH_response
    2.3. Nitrous oxide
        load_N2O_lifetime
        load_hv_response
    2.4. Halogenated compounds
        load_halo_lifetime
    2.5. Ozone tropospheric
        load_regions_HTAP
        load_O3t_regional
        load_O3t_response
    2.6. Ozone stratospheric
        load_ODS_all
        load_O3s_response
    2.7. Aerosols
        load_AER_regional
        load_AER_atmoload
        load_AER_solub
3. RADIATIVE FORCING
    3.1. Greenhouse gases
        load_RFghg_all
    3.2. Short-lived species
        load_RFozo_all
        load_RFaer_direct
        load_RFaer_semidirect
        load_RFaer_indirect
    3.3. Black carbon on snow
        load_regions_Reddy_2007
        load_RFbcsnow_all
    3.4. Land-cover change
        load_RFlcc_all
    3.5. Specific RFs
        load_RF_warmeff
        load_RF_atmfrac
4. CLIMATE
    4.1. Temperature
        load_temp_CMIP5
    4.2. Precipitation
        load_prec_CMIP5
    4.3. Ocean heat content
        load_OHC_all
5. IMPACTS
    5.1. Acidification
        load_pH_all
Z. WRAPPER
    load_all_param
"""

##################################################
##################################################

import os
import numpy as np
import xarray as xr

from core_fct.fct_calib import calib_land_TRENDYv7


##################################################
##   1. CARBON DIOXIDE
##################################################

##===========
## 1.1. Ocean
##===========

## ocean structure and impulse response function
## (Joos et al., 1996; doi:10.3402/tellusb.v48i3.15921) (Appendix)
def load_ocean_struct(**useless):
    ## TODO v3.2: take latest GMD paper?

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_Focean_struct'] = ['box-diffusion', 'HILDA', 'Princeton-2D', 'Princeton-3D']
    Par.coords['box_osurf'] = np.arange(8)

    ## DIC conversion factor
    Par['a_dic'] = xr.DataArray(1.722E17, attrs={'units':'umol m3 ppm-1 kgSW-1'})

    ## preindustrial mixing layer depth
    Par['mld_0'] = xr.DataArray([75., 75., 50.0, 50.9], dims='mod_Focean_struct', attrs={'units':'m'})

    ## global ocean area
    Par['A_ocean'] = xr.DataArray(1E14 * np.array([3.62, 3.62, 3.54, 3.55]), dims='mod_Focean_struct', attrs={'units':'m2'})

    ## preindustrial global ocean temperature
    Par['To_0'] = xr.DataArray(273.15 + np.array([17.7, 18.2, 18.3, 17.7]), dims='mod_Focean_struct', attrs={'units':'K'})

    ## gaseous exchange speed at ocean surface
    Par['v_fg'] = xr.DataArray(1 / np.array([7.80, 9.06, 7.46, 7.66]), dims='mod_Focean_struct', attrs={'units':'yr-1'})

    ## fraction of ocean surface boxes
    ## note: slightly altered from original to have continuous responses
    Par['p_circ'] = xr.DataArray(
        [[0.019737, 0.031528, 0.010469, 0.050469, 0.076817, 0.11803, 0.16851] + [0.52444], # last value as sum of all fast timescales
        [np.nan] + [0.022936, 0.035549, 0.037820, 0.089318, 0.13963, 0.24278] + [0.431967], # last value as sum of all fast timescales
        [np.nan] + [0.013691, 0.012456, 0.026933, 0.026994, 0.036608, 0.067380] + [0.815938], # last value as sum of all fast timescales
        [np.nan, np.nan] + [0.014819, 0.019439, 0.038344, 0.066485, 0.24966] + [0.70367-0.092417]], # offset last value to sum to 1
        dims=['mod_Focean_struct', 'box_osurf'], attrs={'units':'1'})

    ## time-scales of oceanic surface-to-deep transport
    ## note: fastest time-scales also altered to have continuous responses
    ##-------for v3.2------- (to keep)
    Par['t_circ'] = xr.DataArray(
        [[1.e18] + [215.71, 148.77, 43.506, 14.172, 4.8702, 1.6388] + [1.6388/5], # arbitrary fastest time-scale as 20% of previous one
        [np.nan] + [1.e18] + [232.30, 68.736, 18.601, 5.2528, 1.2679] + [1.2679/5], # arbitrary fastest time-scale as 20% of previous one
        [np.nan] + [1.e18] + [331.54, 107.57, 38.946, 11.677, 10.515] + [10.515/5], # arbitrary fastest time-scale as 20% of previous one
        [np.nan, np.nan] + [1.e18] + [347.55, 65.359, 15.281, 2.3488, 0.70177]],
        dims=['mod_Focean_struct', 'box_osurf'], attrs={'units':'yr'})

    ##-------for v3.0------- (to be deleted)
    Par['t_circ'] = xr.DataArray(
        [[1.e18] + [215.71, 148.77, 43.506, 14.172, 4.8702, 1.6388] + [1/3.], # arbitrary fastest time-scale as 1/3 yr
        [np.nan] + [1.e18] + [232.30, 68.736, 18.601, 5.2528, 1.2679] + [1/3.], # arbitrary fastest time-scale as 1/3 yr
        [np.nan] + [1.e18] + [331.54, 107.57, 38.946, 11.677, 10.515] + [1/2.], # arbitrary fastest time-scale as 1/2 yr
        [np.nan, np.nan] + [1.e18] + [347.55, 65.359, 15.281, 2.3488, 0.70177]],
        dims=['mod_Focean_struct', 'box_osurf'], attrs={'units':'yr'})

    ## return
    return Par


## carbonate chemistry emulation hard-coded as process, and based on:
## (Harmann et al. 2011; ISBN: 978-0-643-10745-8)
def load_ocean_chem(**useless):
    ## TODO v3.2: add old function by joos 2001?

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_Focean_chem'] = ['CO2Sys-Pade', 'CO2Sys-Power']

    ## option to choose emulation functional form
    Par['pCO2_is_Pade'] = xr.DataArray([True, False], dims='mod_Focean_chem')

    ## return
    return Par


## transient response of oceanic carbon-cycle
## calibrated on CMIP5 models
def load_ocean_CMIP5(recalibrate=False, **useless):

    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/ocean_CMIP5.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/ocean_CMIP5.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par = calib_ocean_CMIP5()

    ## return
    return Par


##==========
## 1.2. Land
##==========

## land carbon-cycle general parameters
def load_land_misc(**useless):
    ## TODO v3.2: reformulate log/hyp into one extended log formula || add t_shift from LUH2

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_Fland_fert'] = ['log', 'hyp']

    ## option to choose NPP functional form
    Par['fert_is_Log'] = xr.DataArray([True, False], dims='mod_Fland_fert')

    ## characteristic time for shifting cultivation
    ## (Hurtt et al., 2006; doi:10.1111/j.1365-2486.2006.01150.x)
    Par['t_shift'] = xr.DataArray(15., attrs={'units':'yr'})

    ## return
    return Par


## preindustrial land carbon-cycle
## calibrated on TRENDYv7 models
def load_land_TRENDYv7(mod_region, recalibrate=False, path_in='input_data/parameters/', **useless):

    ## load from existing file
    if os.path.isfile(path_in + 'land_TRENDYv7__' + mod_region + '.nc') and not recalibrate:
        with xr.open_dataset(path_in + 'land_TRENDYv7__' + mod_region + '.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        Par = calib_land_TRENDYv7(mod_region=mod_region)

    ## return
    return Par


## transient response of land carbon-cycle
## calibrated on CMIP5 models
def load_land_CMIP5(mod_region, recalibrate=False, **useless):

    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/land_CMIP5__' + mod_region + '.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/land_CMIP5__' + mod_region + '.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par = calib_land_CMIP5(mod_region=mod_region)

    ## return
    return Par


## wood-use parameters
def load_land_wooduse(mod_region, recalibrate=False, **useless):
    ## TODO v3.2: add old Houghton elemental C and partitioning; remove/change BB mod (in pre-processed data)
    
    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_Ehwp_tau'] = ['Houghton_2001', 'Earles_2012']
    Par.coords['mod_Ehwp_speed'] = ['normal', 'fast', 'slow']
    Par.coords['box_hwp'] = np.arange(3)

    ## turnover of harwested wood product pools
    ## (Houghton & Hackler, 2001; doi:10.3334/CDIAC/lue.ndp050)
    ## (Earles et al., 2012; doi:10.1038/nclimate1535)
    Par['t_hwp'] = xr.DataArray([[1., 10., 100.], [0.5, 2., 30.]], 
        dims=['mod_Ehwp_tau', 'box_hwp'], attrs={'units':'yr'})

    ## adjustment factor for HWP decay times
    ## arbitrary, but to compensate the use of exponential decay only
    Par['w_t_hwp'] = xr.DataArray([1., 0.8/-np.log(0.2), 1.5/-np.log(0.3)], 
        dims='mod_Ehwp_speed', attrs={'units':'1'})

    ## burnt fraction of HWP pools
    ## note: arbitrary
    Par['p_hwp_bb'] = xr.DataArray([0.5, 0., 0.], dims='box_hwp')

    ## other parameters
    ## (Earles et al., 2012; doi:10.1038/nclimate1535)
    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/land_Earles_2012__' + mod_region + '.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/land_Earles_2012__' + mod_region + '.nc') as TMP: Par2 = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par2 = calib_land_Earles_2012(mod_region=mod_region)

    ## return
    return xr.merge([Par, Par2])


## biomass burning factors
## calibrated on GFED3
## (Randerson et al., 2013; doi:10.3334/ORNLDAAC/1191)
def load_land_GFED3(mod_region, recalibrate=False, **useless):
    
    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/land_GFED3__' + mod_region + '.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/land_GFED3__' + mod_region + '.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par = calib_land_GFED3(mod_region=mod_region)

    ## return
    return Par


##================
## 1.3. Permafrost
##================

## permafrost carbon parameters
## (Gasser et al., 2018; doi:10.1038/s41561-018-0227-0) (Table S4)
def load_permafrost_all(**useless):

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_Epf_main'] = ['off_', 'JSBACH', 'ORCHIDEE-MICT', 'JULES-DeepResp', 'JULES-SuppressResp']
    Par.coords['reg_pf'] = ['Eurasia', 'North America']
    Par.coords['box_thaw'] = np.arange(3)

    ## preindustrial frozen permafrost carbon pool
    Par['Cfroz_0'] = xr.DataArray([[0., 0.], [414., 277.], [271., 118.], [481., 176.], [373., 121.]], 
        dims=['mod_Epf_main', 'reg_pf'], attrs={'units':'PgC'})

    ## polar amplification of local air temperature
    Par['w_clim_pf'] = xr.DataArray([[1., 1.], [1.86, 1.95], [1.87, 1.78], [1.96, 2.03], [1.96, 2.02]], 
        dims=['mod_Epf_main', 'reg_pf'], attrs={'units':'1'})

    ## permafrost respiration sensitivities to temperature
    Par['g_respT_pf'] = xr.DataArray([[0., 0.], [0.0969, 0.101], [0.114, 0.110], [0.168, 0.144], [0.135, 0.112]], 
        dims=['mod_Epf_main', 'reg_pf'], attrs={'units':'K-1'})
    Par['g_respT2_pf'] = xr.DataArray([[0., 0.], [0.00206, 0.00214], [0.00366, 0.00345], [0.00531, 0.00380], [0.00499, 0.00339]], 
        dims=['mod_Epf_main', 'reg_pf'], attrs={'units':'K-2'})

    ## permafrost respiration buffer factor
    Par['k_resp_pf'] = xr.DataArray([[1., 1.], [0.938, 0.735], [3.38, 3.71], [1.37, 1.58], [0.662, 2.19]], 
        dims=['mod_Epf_main', 'reg_pf'], attrs={'units':'1'})

    ## minimum (negative) thawed fraction
    Par['pthaw_min'] = xr.DataArray([[0., 0.], [0.209, 0.195], [0.0807, 0.875], [0.772, 0.959], [0.870, 1.16]], 
        dims=['mod_Epf_main', 'reg_pf'], attrs={'units':'1'})

    ## thawed fraction sensitivity to temperature
    Par['g_pthaw'] = xr.DataArray([[0., 0.], [0.176, 0.286], [1.62, 0.0973], [0.153, 0.145], [0.180, 0.190]], 
        dims=['mod_Epf_main', 'reg_pf'], attrs={'units':'K-1'})

    ## shape parameter for thawed fraction function
    Par['k_pthaw'] = xr.DataArray([[1., 1.], [2.30, 1.42], [0.357, 8.44], [2.05, 1.91], [1.63, 1.68]], 
        dims=['mod_Epf_main', 'reg_pf'], attrs={'units':'1'})

    ## speed of thawing and freezing
    Par['v_thaw'] = xr.DataArray([[1., 1.], [4.89, 0.500], [0.287, 0.500], [0.266, 0.599], [0.373, 0.492]], 
        dims=['mod_Epf_main', 'reg_pf'], attrs={'units':'yr-1'})
    Par['v_froz'] = xr.DataArray([[1., 1.], [0., 0.], [0.0529, 0.0490], [0.0410, 0.0330], [0.0613, 0.0429]], 
        dims=['mod_Epf_main', 'reg_pf'], attrs={'units':'yr-1'})
 
    ## partitioning coefficient of thawed carbon
    Par['p_pf_thaw'] = xr.DataArray(
        [[[0., 0.], [0.070, 0.100], [0.244, 0.288], [0.034, 0.651], [1., 1.]],
        [[0., 0.], [0.264, 0.234], [0.756, 0.712], [0.966, 0.349], [0., 0.]],
        [[0., 0.], [0.666, 0.666], [0., 0.], [0., 0.], [0., 0.]]],
        dims=['box_thaw', 'mod_Epf_main', 'reg_pf'], attrs={'units':'1'})
 
    ## time-scales of emission of thawed carbon
    Par['t_pf_thaw'] = xr.DataArray(
        [[[1., 1.], [7.96, 11.6], [1330., 2950.], [156., 1970.], [9430., 4320.]],
        [[1.e18, 1.e18], [88.5, 103.], [15000., 23400.], [3160., 5370.], [1.e18, 1.e18]],
        [[1.e18, 1.e18], [1360, 1240], [1.e18, 1.e18], [1.e18, 1.e18], [1.e18, 1.e18]]],
        dims=['box_thaw', 'mod_Epf_main', 'reg_pf'], attrs={'units':'yr'})

    ## fraction of instantaneous permafrost emissions
    Par['p_pf_inst'] =  xr.DataArray(0., attrs={'units':'1'}) 

    ## fraction of permafrost emitted as methane
    Par['p_pf_CH4'] =  xr.DataArray([0., 0.023, 0.046], 
        coords={'mod_Epf_CH4':['zero', 'best', 'twice']}, dims='mod_Epf_CH4',  attrs={'units':'1'})

    ## return
    return Par

    
##==============
## 1.4. Wetlands
##==============

## wetlands parameters
## calibrated on WETCHIMP
def load_wetlands_WETCHIMP(mod_region, recalibrate=False, **useless):
    ## TODO v3.2: properly uncouple from land cover data

    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/wetlands_WETCHIMP__' + mod_region + '.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/wetlands_WETCHIMP__' + mod_region + '.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par = calib_wetlands_WETCHIMP(mod_region=mod_region)

    ## return
    return Par


##################################################
##   2. ATMO CHEMISTRY
##################################################

##================
## 2.1. Atmosphere
##================

## atmospheric conversion factors and preindustrial concentrations
def load_atmosphere_misc(**useless):
    ## TODO v3.2: move p_CH4geo to drivers

    ## initialization
    Par = xr.Dataset()
    Par.coords['spc_halo'] = ['HFC-23', 'HFC-32', 'HFC-125', 'HFC-134a', 'HFC-143a', 'HFC-152a', 'HFC-227ea', 'HFC-236fa', 'HFC-245fa', 'HFC-365mfc', 'HFC-43-10mee',
        'SF6', 'NF3', 'CF4', 'C2F6', 'C3F8', 'c-C4F8', 'C4F10', 'C5F12', 'C6F14', 'C7F16',
        'CFC-11', 'CFC-12', 'CFC-113', 'CFC-114', 'CFC-115', 'CCl4', 'CH3CCl3', 'HCFC-22', 'HCFC-141b', 'HCFC-142b', 'Halon-1211', 'Halon-1202', 'Halon-1301', 'Halon-2402', 'CH3Br', 'CH3Cl']

    ## required values
    ## (Prather et al., 2012; doi:10.1029/2012GL051440) (Table A1)
    ## mass of the atmosphere
    M_atm = 5.113E9 # Tg
    ## mean molar mass of air
    m_air = 28.97 # g mol-1

    ## atmospheric conversion factors
    ## greenhouse gases
    Par['a_CO2'] = xr.DataArray(M_atm / m_air * 1E-9 * 12.0, attrs={'units':'PgC ppm-1'})
    Par['a_CH4'] = xr.DataArray(M_atm / m_air * 1E-9 * 12.0, attrs={'units':'TgC ppb-1'})
    Par['a_N2O'] = xr.DataArray(M_atm / m_air * 1E-9 * 28.0, attrs={'units':'TgN ppb-1'})
    Par['a_Xhalo'] = xr.DataArray(M_atm / m_air * 1E-9 * np.array(
        [70.0, 52.0, 120.0, 102.0, 84.0, 66.0, 170.0, 152.0, 134.0, 148.1, 252.1,
        146.1, 71.0, 88.0, 138.0, 188.0, 200.0, 238.0, 288.0, 338.0, 388.1,
        137.4, 120.9, 187.4, 170.9, 154.5, 153.8, 133.4, 86.5, 117.0, 100.5, 165.4, 209.8, 148.9, 259.8, 94.9, 50.5]),
        dims=['spc_halo'], attrs={'units':'Gg ppt-1'})
    ## aerosols
    Par['a_SO4'] = xr.DataArray(96/32., attrs={'units':'Tg TgS-1'})
    Par['a_POM'] = xr.DataArray([1.4, 1.6, 1.3], coords={'mod_POA_conv':['default', 'GFDL', 'CSIRO']}, dims='mod_POA_conv', attrs={'units':'Tg TgC-1'})
    Par['a_NO3'] = xr.DataArray(62/14., attrs={'units':'Tg TgN-1'})
    
    ## preindustrial concentrations
    ## (Ciais et al., 2013; doi:10.1017/CBO9781107415324.015) (text)
    ## (Prather et al., 2013; doi:10.1017/CBO9781107415324.030) (Table AII.1.1a)
    Par['CO2_0'] = xr.DataArray(278., attrs={'units':'ppm', 'unc (90%)':5})
    Par['CH4_0'] = xr.DataArray(722., attrs={'units':'ppb', 'unc (90%)':25})
    Par['N2O_0'] = xr.DataArray(270., attrs={'units':'ppb', 'unc (90%)':7})
    Par['Xhalo_0'] = xr.DataArray(np.zeros(len(Par.coords['spc_halo'])), dims='spc_halo', attrs={'units':'ppt'})
    Par.Xhalo_0.loc['CF4'] = 35.
    ## (Meinshausen et al., 2011; doi:10.1007/s10584-011-0156-z) (Table 1)
    Par.Xhalo_0.loc['CH3Br'] = 5.8
    Par.Xhalo_0.loc['CH3Cl'] = 480.

    ## fraction of geological CH4 emitted by anthropogenic activities
    ## note: place-holder
    Par['p_CH4geo'] = xr.DataArray(0., attrs={'units':'1'})

    ## return
    return Par


## atmospheric response to climate change
## calibrated on CCMVal2
def load_atmosphere_CCMVal2(recalibrate=False, **useless):
    
    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/atmosphere_CCMVal2.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/atmosphere_CCMVal2.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par = calib_atmosphere_CCMVal2()

    ## return
    return Par


## ad hoc adjustment factors
def load_atmochem_adhoc(**useless):

    ## initialization
    Par = xr.Dataset()

    ## adjustment factor for OH lifetimes
    Par['w_t_OH'] = xr.DataArray(0.80, attrs={'units':'1'})
    ## adjustment factor for hv lifetimes
    Par['w_t_hv'] = xr.DataArray(1.06, attrs={'units':'1'})

    ## return
    return Par


##=============
## 2.2. Methane
##=============

## preindustrial lifetime of CH4
def load_CH4_lifetime(**useless):
    ## TODO v3.2: take ACCMIP as best guess? (and change ad hoc factor)

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_Foh_tau'] = ['Prather_2012', 'CESM-CAM-superfast', 'CICERO-OsloCTM2', 'CMAM', 'EMAC', 'GEOSCCM', 'GDFL-AM3', 
        'GISS-E2-R', 'GISS-E2-R-TOMAS', 'HadGEM2', 'LMDzORINCA', 'MIROC-CHEM', 'MOCAGE', 'NCAR-CAM3.5', 'STOC-HadAM3', 'TM5', 'UM-CAM']
    
    ## OH lifetime
    ## (Prather et al., 2012; doi:10.1029/2012GL051440) (Table A1)
    ## (Naik et al., 2013; 10.5194/acp-13-5277-2013) (Table 1)
    ## note: Prather_2012 is best guess, and other values rescaled to it (models' mean = 9.7 yr)
    Par['t_OH_CH4'] = xr.DataArray(11.2/9.7 * np.array([9.7, 8.4, 10.0, 9.4, 9.1, 9.6, 9.4, 10.6, 9.2, 11.6, 10.5, 8.7, 7.1, 9.2, 9.1, 9.9, 14.0]), 
        dims='mod_Foh_tau', attrs={'units':'yr'})

    ## other lifetimes
    ## (Prather et al., 2012; doi:10.1029/2012GL051440) (Table A1)
    Par['t_hv_CH4'] = xr.DataArray(120., attrs={'units':'yr'})
    Par['t_soil_CH4'] = xr.DataArray(150., attrs={'units':'yr'})
    Par['t_ocean_CH4'] = xr.DataArray(200., attrs={'units':'yr'})

    ## return
    return Par


## OH lifetime sensitivities
def load_OH_response(**useless):
    ## TODO v3.2: switch to only one formulation, and add old OxComp models

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_Foh_fct'] = ['log', 'lin']
    Par.coords['mod_Foh_trans'] = ['Holmes_2013', 'UCI-CTM', 'Oslo-CTM3', 'GEOS-Chem', 'mean_OxComp']

    ## sensitivities of OH sink
    ## (Holmes et al.; doi:10.5194/acp-13-285-2013) (Table 2)
    ## (Ehhalt et al., 2001; IPCC AR3 WG1 Chapter 4) (Table 4.11)
    Par['x_OH_Ta'] = xr.DataArray([3.0, 3.9, 2.8, 2.2, 0.], dims='mod_Foh_trans', attrs={'units':'1'})
    Par['x_OH_Qa'] = xr.DataArray([0.32, 0.32, 0.29, 0.34, 0.], dims='mod_Foh_trans', attrs={'units':'1'})
    Par['x_OH_O3s'] = xr.DataArray([-0.55, -0.66, -0.43, -0.61, 0.], dims='mod_Foh_trans', attrs={'units':'1'})
    Par['x_OH_CH4'] = xr.DataArray([-0.31, -0.363, -0.307, -0.274, -0.32], dims='mod_Foh_trans', attrs={'units':'1'})
    ## log formulation
    Par['x_OH_NOX'] = xr.DataArray([-0.14, -0.15, -0.10, -0.16, -0.137], dims='mod_Foh_trans', attrs={'units':'1'})
    Par['x_OH_CO'] = xr.DataArray([0.06, 0.066, 0.05, 0.065, 0.11], dims='mod_Foh_trans', attrs={'units':'1'})
    Par['x_OH_VOC'] = xr.DataArray([0.04, 0.04, 0.04, 0.04, 0.047], dims='mod_Foh_trans', attrs={'units':'1'})
    ## linear formulation
    Par['x2_OH_NOX'] = xr.DataArray(0.0042/-0.137 * np.array([-0.14, -0.15, -0.10, -0.16, -0.137]), dims='mod_Foh_trans', attrs={'units':'yr TgN-1'})
    Par['x2_OH_CO'] = xr.DataArray(-1.05E-4/0.11 * 28/12. * np.array([0.06, 0.066, 0.05, 0.065, 0.11]), dims='mod_Foh_trans', attrs={'units':'yr TgC-1'})
    Par['x2_OH_VOC'] = xr.DataArray(-3.14E-4/0.047 * np.array([0.04, 0.04, 0.04, 0.04, 0.047]), dims='mod_Foh_trans', attrs={'units':'yr Tg-1'})
   
    ## parameters for global atmospheric temperature and relative humidity changes
    ## (Holmes et al.; doi:10.5194/acp-13-285-2013) (Supp.)
    Par['w_clim_Ta'] = xr.DataArray(0.94, attrs={'units':'1', 'unc (1std)':'0.1'})
    Par['k_Qa'] = xr.DataArray(1.5, attrs={'units':'1', 'unc (1std)':'0.1'})
    Par['Ta_0'] = xr.DataArray(251., attrs={'units':'K', 'unc (1std)':'1'})

    ## parameters for saturating vapor pressure
    ## (Jacobson, 2005; doi:10.1017/cbo9781139165389) (Eq. 2.62)
    Par['k_svp'] = xr.DataArray(17.67, attrs={'units':'1'})
    Par['T_svp'] = xr.DataArray(243.5 - 273.15, attrs={'units':'K'}) # negative sign is not a mistake!

    ## preindustrial stratospheric ozone burden
    ## (Cionni et al., 2011; doi:10.5194/acp-11-11267-2011) (Fig. 7; guesstimated)
    Par['O3s_0'] = xr.DataArray(280., attrs={'units':'DU'})

    ## preindustrial ozone precursors emission
    ## (Skeie et al., 2011; doi:) (Table 1)
    Par['Enat_NOX'] = xr.DataArray(13., attrs={'units':'TgN yr-1'})
    Par['Enat_CO'] = xr.DataArray(180. * 12/28., attrs={'units':'TgC yr-1'})
    Par['Enat_VOC'] = xr.DataArray(39. + 220 + 175, attrs={'units':'Tg yr-1'})

    ## option to choose kOH functional form
    Par['kOH_is_Log'] = xr.DataArray([True, False], dims='mod_Foh_fct')

    ## return
    return Par


##===================
## 2.3. Nitrous Oxide
##===================

## preindustrial lifetime of N2O
def load_N2O_lifetime(**useless):

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_Fhv_tau'] = ['Prather_2015', 'GMI', 'GEOSCCM', 'G2d-M', 'G2d', 'Oslo-c29', 'Oslo-c36', 'UCI-c29', 'UCI-c36']
    
    ## strato lifetime
    ## (Prather et al., 2015; doi:10.1002/2015JD023267) (Table 2 and text)
    ## note: Prather_2015 is best guess, and other values rescaled to it (models' mean = 132.5 yr)
    Par['t_hv_N2O'] = xr.DataArray(123./132.5 * np.array([132.5, 137.4, 120.2, 127.0, 129.5, 126.1, 146.7, 126.2, 146.2]), dims='mod_Fhv_tau', attrs={'units':'yr'})

    ## return
    return Par


## hv lifetime sensitivities
def load_hv_response(**useless):
    ## TODO v3.2: check/add old MAGICC model? set old Prather_2012

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_Fhv_trans'] = ['Prather_2015', 'G2d', 'Oslo-c29', 'UCI-c29', 'Prather_2012']

    ## required values
    ## age of air
    ## (Fleming et al., 2011; 10.5194/acp-11-8515-2011) (Fig. 12; guesstimated)
    age_pd, age_pi = 4.0, 4.5 # yr
    
    ## equivalent effective stratospheric chlorine
    ## (Newman et al., 2007; doi:10.5194/acp-7-4537-2007) (Table 1)
    ## note: actual values in models unknown, replaced by rounded estimates based on 1750 and 2005 gas concentrations
    EESC_pd, EESC_pi = 1730., 420. # ppt

    ## sensitivities of hv sink
    ## (Prather et al., 2015; doi:10.1002/2015JD023267) (Table 3)
    Par['x_hv_N2O'] = xr.DataArray([0.065, 0.018/np.log(321/270.), 0.010/np.log(321/270.), 0.012/np.log(321/270.), 0.08], dims='mod_Fhv_trans', attrs={'units':'1'})
    Par['x_hv_EESC'] = xr.DataArray(1/np.log(EESC_pd/EESC_pi) * np.array([0.04, 0.048, 0.033, 0.029, 0.]), dims='mod_Fhv_trans', attrs={'units':'1'})
    Par['x_hv_ageair'] = xr.DataArray([0., 0.011/np.log(age_pd/age_pi), 0., 0., 0.], dims='mod_Fhv_trans', attrs={'units':'1'})

    ## return
    return Par


##===========================
## 2.4. Halogenated compounds
##===========================

## preindustrial lifetime of halo comp.
# (WMO, 2011; ISBN: 9966-7319-6-2) (Table 1-3)
def load_halo_lifetime(**useless):
    ## TODO v3.2: change loading and update to WMO 2014 and 2018

    ## initialization
    Par = xr.Dataset()
    Par.coords['spc_halo'] = ['HFC-23', 'HFC-32', 'HFC-125', 'HFC-134a', 'HFC-143a', 'HFC-152a', 'HFC-227ea', 'HFC-236fa', 'HFC-245fa', 'HFC-365mfc', 'HFC-43-10mee',
        'SF6', 'NF3', 'CF4', 'C2F6', 'C3F8', 'c-C4F8', 'C4F10', 'C5F12', 'C6F14', 'C7F16',
        'CFC-11', 'CFC-12', 'CFC-113', 'CFC-114', 'CFC-115', 'CCl4', 'CH3CCl3', 'HCFC-22', 'HCFC-141b', 'HCFC-142b', 'Halon-1211', 'Halon-1202', 'Halon-1301', 'Halon-2402', 'CH3Br', 'CH3Cl']

    ## OH lifetimes
    Par['t_OH_Xhalo'] = xr.DataArray(
        [245., 5.5, 32., 14.3, 55, 1.6, 44.5, 253., 8.2, 9.3, 17.9,
        1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18,
        1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 6.1, 12.8, 10.7, 19.3, 1.e18, 1.e18, 1.e18, 1.e18, 1.5, 1.9],
        dims=['spc_halo'], attrs={'units':'yr'})

    ## hv lifetimes
    Par['t_hv_Xhalo'] = xr.DataArray(
        [2347., 89., 246., 232., 327., 45.4, 310., 5676., 116., 125., 157.,
        3200., 500., 50000., 10000., 2600., 3200., 2600., 4100., 3100., 3000.,
        45., 100., 85., 190., 1020., 35., 39., 186., 64.9, 160., 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18],
        dims=['spc_halo'], attrs={'units':'yr'})

    ## other sinks lifetimes
    Par['t_other_Xhalo'] = xr.DataArray(
        [1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18,
        1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 1.e18,
        1.e18, 1.e18, 1.e18, 1.e18, 1.e18, 101., 89., 1.e18, 1.e18, 1.e18, 16., 2.9, 65., 20., 3., 1.4],
        dims=['spc_halo'], attrs={'units':'yr'})

    ## return
    return Par


##========================
## 2.5. Tropospheric ozone
##========================

## HTAP regions split fractions
def load_regions_HTAP(mod_region, recalibrate=False, **useless):
    ## TODO v3.2: maybe? correct fractions to account only for land...
    
    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/regions_HTAP__' + mod_region + '.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/regions_HTAP__' + mod_region + '.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par = calib_regions_HTAP(mod_region=mod_region)

    ## return
    return Par


## regional weighting for ozone precursors effect
## taken from HTAP experiments
def load_O3t_regional(**useless):

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_O3t_regsat'] = ['mean_HTAP', 'CAMCHEM', 'FRSGCUCI', 'GISS-PUCCINI-modelE', 'GMI', 'INCA', 'LLNL-IMPACT', 'MOZARTGFDL', 'MOZECH', 'STOC-HadAM3', 'TM5-JRC', 'UM-CAM']
    Par.coords['reg_slcf'] = ['Globe', 'EA', 'EU', 'NA', 'SA']
    
    ## temporary: tropospheric O3 perturbations
    ## (Fry et al., 2012; doi:10.1029/2011JD017134) (Table S5)
    ## raw data
    Par['D_O3t_NOX'] = xr.DataArray([[np.nan for _ in Par['mod_O3t_regsat']],
        [np.nan, -0.91, -0.55, -0.62, -0.71, -0.75, -0.67, -0.50, -0.81, -0.50, -0.46, -0.98],
        [np.nan, -0.67, -0.44, -0.22, -0.55, -0.54, -0.56, -0.36, -0.57, -0.58, -0.27, -0.72],
        [np.nan, -1.32, -0.62, -0.50, -1.14, -0.96, -1.13, -0.82, -1.18, -0.81, -0.59, -1.32],
        [np.nan, -0.41, -0.30, -0.24, -0.44, -0.36, -0.72, -0.31, -0.44, -0.35, -0.11, -0.78]],
        dims=['reg_slcf', 'mod_O3t_regsat'])
    Par['D_O3t_CO'] = xr.DataArray([[np.nan for _ in Par['mod_O3t_regsat']],
        [np.nan, -0.58, -0.49, -0.90, -0.87, -0.33, -0.40, -0.36, -0.48, -0.40, -0.42, -0.44],
        [np.nan, -0.22, -0.28, -0.36, -0.38, -0.17, -0.24, -0.35, -0.31, -0.25, -0.16, -0.28],
        [np.nan, -0.36, -0.47, -0.61, -0.52, -0.18, -0.31, -0.34, -0.46, -0.39, -0.54, -0.41],
        [np.nan, -0.33, -0.30, -0.36, -0.42, -0.23, -0.46, -0.24, -0.36, -0.25, -0.22, -0.29]],
        dims=['reg_slcf', 'mod_O3t_regsat'])
    Par['D_O3t_VOC'] = xr.DataArray([[np.nan for _ in Par['mod_O3t_regsat']],
        [np.nan, -0.21, -0.68, -0.14, -0.38, -0.30, -0.04, -0.20, -0.45, -0.45, -0.58, -0.58],
        [np.nan, -0.59, -0.78, -0.08, -0.15, -0.35, -0.10, -0.28, -0.55, -0.64, -0.49, -0.64],
        [np.nan, -0.29, -0.60, -0.08, -0.31, -0.40, -0.05, -0.19, -0.69, -0.40, -0.56, -0.38],
        [np.nan, -0.16, -0.36, -0.11, -0.24, 0., -0.03, -0.17, -0.25, -0.27, -0.28, -0.31]],
        dims=['reg_slcf', 'mod_O3t_regsat'])
    ## processing
    for var in ['D_O3t_NOX', 'D_O3t_CO', 'D_O3t_VOC']:
        ## get global
        Par[var].loc['Globe',:] = Par[var].sum('reg_slcf')
        ## get mean
        Par[var].loc[:,'mean_HTAP'] = Par[var].mean('mod_O3t_regsat')

    ## temporary: ozone precursors emission perturbations
    ## (Fiore et al., 2009; doi:10.1029/2008JD010816) (Table S2)
    ## raw_data
    Par['D_E_NOX'] = xr.DataArray(-0.2 * np.array(
        [[np.nan for _ in Par['reg_slcf']], [np.nan, 9.2, 8.2, 6.8, 3.4], [np.nan, 8.4, 8.8, 6.9, 3.0], 
        [np.nan, 8.6, 8.9, 10.8, 3.1], [np.nan, 7.5, 8.0, 5.8, 3.1], [np.nan, 7.3, 8.6, 7.2, 3.8],
        [np.nan, 9.7, 8.8, 7.3, 4.2], [np.nan, 8.6, 9.4, 5.2, 2.6], [np.nan, 7.4, 8.9, 7.0, 3.6], 
        [np.nan, 8.6, 9.0, 7.1, 3.2], [np.nan, 8.5, 8.7, 6.4, 3.0], [np.nan, 8.7, 8.9, 7.0, 3.3]]),
        dims=['mod_O3t_regsat', 'reg_slcf'])
    Par['D_E_CO'] = xr.DataArray(-0.2 * np.array(
        [[np.nan for _ in Par['reg_slcf']], [np.nan, 106., 154., 221., 145.], [np.nan, 74., 131., 129., 74.],
        [np.nan, 69., 111., 170., 67.], [np.nan, 101., 134., 194., 99.], [np.nan, 67, 74., 127., 108.],
        [np.nan, 111., 130., 153., 124.], [np.nan, 130., 124., 134., 105.], [np.nan, 85., 107., 154., 124.],
        [np.nan, 74., 129., 127., 76.], [np.nan, 75., 127., 123., 74.], [np.nan, 81., 137., 131., 79.]]),
        dims=['mod_O3t_regsat', 'reg_slcf'])
    Par['D_E_VOC'] = xr.DataArray(-0.2 * np.array(
        [[np.nan for _ in Par['reg_slcf']], [np.nan, 60.6, 92.4, 55.5, 36.7], [np.nan, 33.5, 57.3, 50.5, 33.8],
        [np.nan, 26.9, 29.0, 31.7, 24.1], [np.nan, 32.0, 60.1, 42.0, 31.6], [np.nan, 56.7, 100.1, 77.0, 45.3],
        [np.nan, 23.3, 53.2, 36.2, 30.7], [np.nan, 32.4, 68.7, 48.9, 36.8], [np.nan, 44.7, 107.0, 62.0, 36.9],
        [np.nan, 52.2, 70.9, 50.0, 43.8], [np.nan, 34.6, 50.3, 47.0, 32.4], [np.nan, 32.7, 56.4, 47.4, 34.6]]),
        dims=['mod_O3t_regsat', 'reg_slcf'])
    ## processing
    for var in ['D_E_NOX', 'D_E_CO', 'D_E_VOC']:
        ## get global
        Par[var].loc[:,'Globe'] = Par[var].sum('reg_slcf')
        ## get mean
        Par[var].loc['mean_HTAP':,] = Par[var].mean('mod_O3t_regsat')

    ## final: regional normalized sensitivities
    for var in ['NOX', 'CO', 'VOC']:
        ## per unit emission
        Par['w_reg_'+var] = Par['D_O3t_'+var] / Par['D_E_'+var]
        ## normalization
        Par['w_reg_'+var] /= Par['w_reg_'+var].sel(reg_slcf='Globe', drop=True)
        Par['w_reg_'+var].attrs['units'] = '1'
        ## remove temporary parameters
        del Par['D_O3t_'+var], Par['D_E_'+var]

    ## add option to turn off
    Par2 = xr.Dataset()
    for var in ['NOX', 'CO', 'VOC']:
        Par2['w_reg_'+var] = xr.DataArray([[1. for _ in Par['reg_slcf']]], 
            coords={'mod_O3t_regsat':['off_'], 'reg_slcf':Par['reg_slcf']}, dims=['mod_O3t_regsat', 'reg_slcf'], attrs={'units':'1'})

    ## return
    return xr.merge([Par, Par2])

    
## tropospheric ozone sensitivities
## partly calibrated on ACCMIP
def load_O3t_response(recalibrate=False, **useless):
    ## TODO v3.2: add old OxComp
    
    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/ozochem_ACCMIP.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/ozochem_ACCMIP.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par = calib_ozochem_ACCMIP()

    ## extra (old) parameterizations
    ## (Ehhalt et al., 2001; IPCC AR3 WG1 Chapter 4) (Table 4.11)
    Par2 = xr.Dataset()
    Par2['x_O3t_CH4'] = xr.DataArray([5.0], coords={'mod_O3t_emis':['mean_OxComp']}, dims='mod_O3t_emis')
    Par2['x_O3t_NOX'] = xr.DataArray([0.125], coords={'mod_O3t_emis':['mean_OxComp']}, dims='mod_O3t_emis', attrs={'units':'DU yr TgN-1'})
    Par2['x_O3t_CO'] = xr.DataArray([28/12. * 0.0011], coords={'mod_O3t_emis':['mean_OxComp']}, dims='mod_O3t_emis', attrs={'units':'DU yr TgC-1'})
    Par2['x_O3t_VOC'] = xr.DataArray([0.0033], coords={'mod_O3t_emis':['mean_OxComp']}, dims='mod_O3t_emis', attrs={'units':'DU yr Tg-1'})

    ## return
    return xr.merge([Par, Par2])


##=========================
## 2.6. Stratospheric ozone
##=========================

## ozone depleting substances parameters
def load_ODS_all(**useless):
    ## TODO v3.2: update fracrel on e.g. WMO, and replace Daniel 2010 with saturating formula
    
    ## initialization
    Par = xr.Dataset()
    Par.coords['spc_halo'] = ['CFC-11', 'CFC-12', 'CFC-113', 'CFC-114', 'CFC-115', 'CCl4', 'CH3CCl3', 'HCFC-22', 'HCFC-141b', 'HCFC-142b', 'Halon-1211', 'Halon-1202', 'Halon-1301', 'Halon-2402', 'CH3Br', 'CH3Cl']
    Par.coords['mod_O3s_fracrel'] = ['Newman_2007', 'Laube_2013']

    ## time lag for lagged concentrations
    ## note: 3 yr is typical value used by WMO, and fractional release factors depend on it
    Par['t_lag'] = xr.DataArray(3., attrs={'units':'yr'})

    ## fractional release factors
    ## (Newman et al., 2007; doi:10.5194/acp-7-4537-2007) (Table 1)
    ## (Laube et al., 2013; doi:10.5194/acp-13-2779-2013) (Table 3)
    Par['p_fracrel'] = xr.DataArray(
        [[0.47, 0.23, 0.29, 0.12, 0.04, 0.56, 0.67, 0.13, 0.08, 0.01, 0.62, 0.62, 0.28, 0.65, 0.60, 0.44], 
        [0.35, 0.19, 0.22, np.nan, np.nan, 0.42, 0.61, 0.07, 0.17, 0.05, 0.52, np.nan, 0.26, np.nan, np.nan, np.nan]],
        dims=['mod_O3s_fracrel', 'spc_halo'], attrs={'units':'1'})
    ## replace missing values by existing ones
    Par.p_fracrel.loc['Laube_2013',:] = Par['p_fracrel'].sel(mod_O3s_fracrel='Laube_2013').fillna(Par['p_fracrel'].sel(mod_O3s_fracrel='Newman_2007'))

    ## relative strength of bromine vs. chlorine
    ## (WMO, 2007; ISBN: 978-92-807-2756-2) (Chapter 8)
    Par['k_Br_Cl'] = xr.DataArray(60., attrs={'units':'1'})

    ## number of Cl and Br atoms
    Par['n_Cl'] = xr.DataArray([3, 2, 3, 2, 1, 4, 3, 1, 2, 1, 1, 0, 0, 0, 0, 1], dims='spc_halo', attrs={'units':'1'})
    Par['n_Br'] = xr.DataArray([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 2, 1, 0], dims='spc_halo', attrs={'units':'1'})

    ## sensitivity factors for N2O
    ## (Daniel et al., 2010; doi:10.5194/acp-10-7697-2010)
    ## note: taken from values in text, and based on ODP estimated for N2O
    ##-------for v3.2------- (to keep)
    ## saturating factor
    Par['EESC_x'] = xr.DataArray((240 * 1.53 - 1640 * 1) / (1 - 1.53), attrs={'units':'ppt'})
    ## conversion factor
    Par['k_EESC_N2O'] = xr.DataArray(np.array([0, 6.4]) * (1 + 1640 / Par['EESC_x'].values), coords={'mod_O3s_nitrous':['off_', 'Daniel_2010']}, dims='mod_O3s_nitrous')
    Par['k_EESC_N2O'] = Par['k_EESC_N2O'] * Par['p_fracrel'].sel({'spc_halo':'CFC-11'}, drop=True)
    Par['k_EESC_N2O'].attrs['units'] = 'ppt ppb-1'

    ##-------for v3.0------- (to be deleted)
    ## saturating factor
    Par['EESC_x'] = xr.DataArray((240 * 1 - 1640 * 1.53) / (1 - 1.53), attrs={'units':'ppt'})
    ## conversion factor
    Par['k_EESC_N2O'] = xr.DataArray(np.array([0, 6.4]) / (1 - 1640 / Par['EESC_x'].values), coords={'mod_O3s_nitrous':['off_', 'Daniel_2010']}, dims='mod_O3s_nitrous')
    Par['k_EESC_N2O'] = Par['k_EESC_N2O'] * Par['p_fracrel'].sel({'spc_halo':'CFC-11'}, drop=True)
    Par['k_EESC_N2O'].attrs['units'] = 'ppt ppb-1'

    ## return
    return Par


## stratospheric ozone sensitivities
## taken from CCMVal2
def load_O3s_response(recalibrate=False, **useless):

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_O3s_trans'] = ['mean_CCMVal2', 'AMTRAC', 'CCSR-NIES', 'CMAM', 'CNRM-ACM', 'LMDZrepro', 'MRI', 'Niwa-SOCOL', 'SOCOL', 'ULAQ', 'UMSLIMCAT', 'UMUKCA-UCAM']

    ## sensitivities to EESC and climate
    ## (Douglass et al.; 2014; doi:10.1002/2013JD021159) (Fig. 2 & extra; data provided by author)
    Par['x_O3s_EESC'] = xr.DataArray(1E-3 * np.array([-12.5, -13.7, -7.9, -8.1, -23.2, -10.6, -20.3, -5.2, -10.5, -15.4, -10.6, -11.5]),
        dims='mod_O3s_trans', attrs={'units':'DU ppt-1'})
    Par['G_O3s'] = xr.DataArray(
        np.array([0.012, -0.015, 0.013, -0.15, 0.0066, 0.042, -0.030, 0.087, 0.067, 0.062, 0.0071, 0.040]) / # trend in strato O3
        np.array([0.0210, 0.0202, 0.0188, 0.0309, 0.0200, 0.0136, 0.0157, 0.0127, 0.0220, 0.0278, 0.0237, 0.0234]), # trend in surface T
        dims='mod_O3s_trans', attrs={'units':'DU K-1'})

    ## return
    return Par


##==============
## 2.7. Aerosols
##==============

## regional weighting for (main) aerosols
## taken from HTAP experiments
## (Yu et al., 2013; doi:10.1029/2012JD018148) (Table 6 detailed; data provided by author)
def load_AER_regional(**useless):

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_SO4_regsat'] = ['mean_HTAP', 'CAMCHEM', 'GISS-PUCCINI', 'GMI', 'GOCART', 'INCA', 'LLNL-IMPACT', 'SPRINTARS']
    Par.coords['mod_POA_regsat'] = ['mean_HTAP', 'CAMCHEM', 'GISS-PUCCINI', 'GMI', 'GOCART', 'INCA', 'LLNL-IMPACT', 'SPRINTARS']
    Par.coords['mod_BC_regsat'] = ['mean_HTAP', 'CAMCHEM', 'GISS-PUCCINI', 'GMI', 'GOCART', 'INCA', 'LLNL-IMPACT', 'SPRINTARS']
    Par.coords['reg_slcf'] = ['Globe', 'NA', 'EU', 'EA', 'SA']

    ## sulfates
    Par['w_reg_SO2'] = xr.DataArray(
        [[-3.51, -3.87, -3.92, -2.85, -3.92],
        [-4.26, -5.01, -4.86, -3.29, -4.55],
        [-2.73, -3.88, -2.92, -1.78, -3.44],
        [-3.61, -3.55, -3.96, -3.30, -3.59],
        [-4.05, -3.86, -4.68, -3.79, -3.06],
        [-4.03, -4.52, -4.27, -3.33, -5.45],
        [-3.70, -3.74, -4.27, -2.84, -4.82],
        [-2.16, -2.51, -2.53, -1.64, -2.56]],
        dims=['mod_SO4_regsat', 'reg_slcf'], attrs={'units':'W m-2'})

    ## primary organic aerosols
    Par['w_reg_OC'] = xr.DataArray(
        [[-4.00, -4.39, -4.30, -3.68, -4.09],
        [-3.30, -3.90, -3.23, -2.80, -3.71],
        [-6.26, -5.86, -6.07, -6.44, -6.38],
        [-4.62, -5.30, -6.13, -4.22, -4.16],
        [-3.57, -3.65, -3.39, -3.28, -4.00],
        [-3.07, -4.27, -3.33, -2.88, -2.66],
        [-1.31, -1.41, -1.97, -0.99, -1.29],
        [-5.86, -6.32, -5.97, -5.12, -6.45]],
        dims=['mod_POA_regsat', 'reg_slcf'], attrs={'units':'W m-2'})

    ## black carbon
    Par['w_reg_BC'] = xr.DataArray(
        [[29.51, 27.31, 37.36, 28.36, 25.31],
        [27.56, 28.00, 35.71, 25.24, 24.08],
        [60.41, 51.67, 69.53, 65.06, 45.46],
        [26.68, 25.80, 42.13, 24.81, 15.00],
        [46.20, 42.30, 52.21, 45.87, 43.71],
        [17.32, 16.88, 20.37, 14.14, 23.69],
        [7.25, 6.63, 12.99, 5.66, 5.56],
        [21.16, 19.93, 28.58, 17.75, 19.69]],
        dims=['mod_BC_regsat', 'reg_slcf'], attrs={'units':'W m-2'})

    ## processing
    for var in ['w_reg_SO2', 'w_reg_OC', 'w_reg_BC']:
        Par[var] /= Par[var].sel({'reg_slcf':'Globe'}, drop=True)
        Par[var].attrs['units'] = '1'

    ## add option to turn off
    dic = {'SO2':'SO4', 'OC':'POA', 'BC':'BC'}
    Par2 = xr.Dataset()
    for var in ['SO2', 'OC', 'BC']:
        Par2['w_reg_'+var] = xr.DataArray([[1. for _ in Par['reg_slcf']]], 
            coords={'mod_'+dic[var]+'_regsat':['off_'], 'reg_slcf':Par['reg_slcf']}, dims=['mod_'+dic[var]+'_regsat', 'reg_slcf'], attrs={'units':'1'})

    ## return
    return xr.merge([Par, Par2])


## aerosols tropospheric load
## calibrated on ACCMIP and CMIP5
def load_AER_atmoload(recalibrate=False, **useless):
    ## TODO v3.2: add dust and salt
    
    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/aerchem_ACCMIP.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/aerchem_ACCMIP.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par = calib_aerchem_ACCMIP()

    ## add option to turn off nitrate, dust and salt aerosols
    Par2 = xr.Dataset()
    Par2['t_NH3'] = Par2['t_NOX'] = xr.DataArray([0.], coords={'mod_NO3_load':['off_']}, dims='mod_NO3_load', attrs={'units':'yr'})
    Par2['G_NO3'] = xr.DataArray([0.], coords={'mod_NO3_load':['off_']}, dims='mod_NO3_load', attrs={'units':'Tg K-1'})
    for var in ['dust', 'salt']:
        Par2['t_'+var] = xr.DataArray([0.], coords={'mod_M'+var+'_load':['zero']}, dims='mod_M'+var+'_load', attrs={'units':'yr'})
        Par2['G_'+var] = xr.DataArray([0.], coords={'mod_M'+var+'_load':['zero']}, dims='mod_M'+var+'_load', attrs={'units':'Tg K-1'})

    ## note: parameters for nitrate calibrated on data extracted from papers
    """
    ## (Bellouin et al., 2011; doi:10.1029/2011JD016074) (and additional CMIP5 and RCP data)
    ENOX = np.array([5.7, 37.4, 18.4, 16.3, 16.6, 23.8, 18.4, 16.3, 16.6, 23.8])
    ENH3 = np.array([16.6, 41.4, 67.2, 49.2, 63.0, 70.0, 67.2, 49.2, 63.0, 70.0])
    tas = np.array([13.55, 14.07, 15.39, 16.46, 17.07, 18.66, 13.55, 13.55, 13.55, 13.55])
    NO3 = np.array([0.05, 0.34, 0.56, 0.29, 0.41, 0.52, 0.63, 0.36, 0.50, 0.68])

    ## (Hauglustaine et al., 2014; doi:10.5194/acp-14-11031-2014) (Tables 1 & 5)
    ENOX = np.array([10., 36., 29., 26., 14., 32., 26., 14., 30., 27., 13., 38., 30., 21., 21., 14.])
    ENH3 = np.array([21., 29., 41., 46., 58., 35., 36., 33., 36., 43., 51., 42., 48., 57., 33., 57.])
    tas = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
    #NO3 = np.array([0.23, 0.48, 0.46, 0.47, 0.37, 0.48, 0.46, 0.42, 0.47, 0.48, 0.40, 0.54, 0.52, 0.52, 0.47, 0.43]) # HNO3 and NO3-
    NO3 = np.array([0.09, 0.18, 0.21, 0.23, 0.21, 0.20, 0.20, 0.18, 0.19, 0.21, 0.21, 0.23, 0.24, 0.25, 0.20, 0.22]) # NO3- only
    """

    ##-------for v3.0------- (to be deleted)
    ## remove other dust and salt options
    Par = Par.drop(['t_dust', 'G_dust', 't_salt', 'G_salt', 'mod_Mdust_load', 'mod_Msalt_load'])

    ## return
    return xr.merge([Par, Par2])


## soluble aerosols fraction
def load_AER_solub(**useless):
    ## TODO v3.2: check p_sol_salt in Lamarque_2011! WARNING: this whould change Phi_0 as well!

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_RFcloud_solub'] = ['Hansen_2005', 'Lamarque_2011']

    ## solubility of aerosols for the aerosol-cloud interaction
    ## (Hansen et al., 2005; doi:10.1029/2005JD005776) (Sect. 3.3.1)
    ## (Lamarque et al., 2011; doi:10.1007/s10584-011-0155-0) (data from RCP database and provided by author)
    ## note: in Hansen_2005, p_sol_BC taken as average of FF and BB, and p_sol_SOA taken as p_sol_POA
    Par['p_sol_SO4'] = xr.DataArray([1., 1.], dims='mod_RFcloud_solub', attrs={'units':'1'})
    Par['p_sol_POA'] = xr.DataArray([0.8, 0.86], dims='mod_RFcloud_solub', attrs={'units':'1'})
    Par['p_sol_BC'] = xr.DataArray([0.7, 0.80], dims='mod_RFcloud_solub', attrs={'units':'1'})
    Par['p_sol_NO3'] = xr.DataArray([1., 1.], dims='mod_RFcloud_solub', attrs={'units':'1'})
    Par['p_sol_SOA'] = xr.DataArray([0.8, 1.], dims='mod_RFcloud_solub', attrs={'units':'1'})
    Par['p_sol_dust'] = xr.DataArray([0., 0.12], dims='mod_RFcloud_solub', attrs={'units':'1'})
    Par['p_sol_salt'] = xr.DataArray([1., 0.05], dims='mod_RFcloud_solub', attrs={'units':'1'})

    ## return
    return Par


##################################################
##   3. RADIATIVE FORCING
##################################################

##======================
## 3.1. Greenhouse gases
##======================

## radiative forcing factors of WMGHGs
## based on IPCC AR5
def load_RFghg_all(**useless):
    ## TODO v3.2: extend to Etminan formulas

    ## initialization
    Par = xr.Dataset()
    Par.coords['spc_halo'] = ['HFC-23', 'HFC-32', 'HFC-125', 'HFC-134a', 'HFC-143a', 'HFC-152a', 'HFC-227ea', 'HFC-236fa', 'HFC-245fa', 'HFC-365mfc', 'HFC-43-10mee',
        'SF6', 'NF3', 'CF4', 'C2F6', 'C3F8', 'c-C4F8', 'C4F10', 'C5F12', 'C6F14', 'C7F16',
        'CFC-11', 'CFC-12', 'CFC-113', 'CFC-114', 'CFC-115', 'CCl4', 'CH3CCl3', 'HCFC-22', 'HCFC-141b', 'HCFC-142b', 'Halon-1211', 'Halon-1202', 'Halon-1301', 'Halon-2402', 'CH3Br', 'CH3Cl']

    ## radiative forcing factors
    ## (Myrhe et al., 2013; doi:10.1017/CBO9781107415324.018) (Table 8.SM.1)
    Par['rf_CO2'] = xr.DataArray(5.35, attrs={'units':'W m-2'})
    Par['rf_CH4'] = xr.DataArray(0.036, attrs={'units':'W m-2 ppb-0.5'})
    Par['rf_N2O'] = xr.DataArray(0.12, attrs={'units':'W m-2 ppb-0.5'})
    ## (Myrhe et al., 2013; doi:10.1017/CBO9781107415324.018) (text)
    Par['k_rf_H2Os'] = xr.DataArray(0.15, attrs={'units':'1'})

    ## radiative efficiency of halogenated compounds
    ## (Myrhe et al., 2013; doi:10.1017/CBO9781107415324.018) (Table 8.A.1)
    Par['rf_Xhalo'] = xr.DataArray(1E-3 * np.array(
        [0.18, 0.11, 0.23, 0.16, 0.16, 0.10, 0.26, 0.24, 0.24, 0.22, 0.42,
        0.57, 0.20, 0.09, 0.25, 0.28, 0.32, 0.36, 0.41, 0.44, 0.50,
        0.26, 0.32, 0.30, 0.31, 0.20, 0.17, 0.07, 0.21, 0.16, 0.19, 0.29, 0.27, 0.30, 0.31, 0.004, 0.01]),
        dims=['spc_halo'], attrs={'units':'W m-2 ppt-1'})

    ## return
    return Par


##=========================
## 3.2. Short-lived species
##=========================

# radiative efficiency of (tropo and strato) ozone
## taken from ACCMIP and ACCENT
def load_RFozo_all(**useless):

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_O3t_radeff'] = ['Myhre_2013', 'Forster_2007', 'mean_ACCMIP', 
        'CESM-CAM-superfast', 'CICERO-OsloCTM2', 'CMAM', 'EMAC', 'GEOSCCM', 'GFDL-AM3', 'GISS-E2-R', 'GISS-E2-R-TOMAS', 'HadGEM2', 'LMDzORINCA', 'MIROC-CHEM', 'MOCAGE', 'NCAR-CAM3.5', 'STOC-HadAM3', 'UM-CAM', 'TM5']
    Par.coords['mod_O3s_radeff'] = ['Forster_2007', 'mean_ACCENT', 'ULAQ', 'DLR-E39C', 'NCAR-MACCM', 'CHASER']

    ## radiative efficiency of tropospheric O3
    ## 1st value from IPCC AR5 (Myhre et al., 2013; doi:10.1017/CBO9781107415324.018)
    ## 2nd value from IPCC AR4 (Forster et al., 2007; ISBN: 978-0521-70596-7) (Chapter 2)
    ## all others from ACCMIP (Stevenson et al., 2013) (Table 3; MASK150)
    Par['rf_O3t'] = xr.DataArray(1E-3 * 
        np.array([42., 32., 377, 446, 401, 322, 460, 387, 423 , 314, 333, 303, 351, 402, 219, 433, 437, 376, 422]) /
        np.array([1., 1., 8.9, 10.0, 9.3, 7.6, 10.8, 8.7, 10.3, 8.3, 8.7, 7.3, 8.2, 9.2, 4.8, 10.2, 10.5, 8.7, 10.0]),
        dims='mod_O3t_radeff', attrs={'units':'W m-2 DU-1'})

    ## radiative efficiency of stratospheric O3
    ## 1st value from IPCC AR4 (Forster et al., 2007; ISBN: 978-0521-70596-7) (Chapter 2)
    ## all others from ACCENT (Gauss et al., 2006; doi:10.5194/acp-6-575-2006) (Tables 4 & 6)
    Par['rf_O3s'] = xr.DataArray(
        np.array([0.004, -0.058, -0.059, -0.027, -0.019, -0.126]) /
        np.array([1., -13.9, -12.6, -16.1, -12.7, -14.1]),
        dims='mod_O3s_radeff', attrs={'units':'W m-2 DU-1'})

    ## return
    return Par


## radiative efficiency for direct effect of aerosols
## taken from AeroCom2
## (Myhre et al., 2013; doi:10.5194/acp-13-1853-2013) (Tables 4-8)
def load_RFaer_direct(**useless):
    
    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_SO4_radeff'] = ['BCC', 'CAM4−Oslo', 'CAM5.1', 'GEOSCHEM', 'GISS-MATRIX', 'GISS-modelE', 'GMI', 'GOCART', 'HadGEM2', 'IMPACT-Umich', 'INCA', 'MPIHAM', 'NCAR-CAM3.5', 'OsloCTM2', 'SPRINTARS', 'mean_AeroCom2']
    Par.coords['mod_POA_radeff'] = ['BCC', 'CAM4-Oslo', 'CAM5.1', 'GEOSCHEM', 'GISS-MATRIX', 'GISS-modelE', 'GMI', 'GOCART', 'HadGEM2', 'IMPACT-Umich', 'INCA', 'MPIHAM', 'NCAR-CAM3.5', 'OsloCTM2', 'SPRINTARS', 'mean_AeroCom2']
    Par.coords['mod_BC_radeff'] = ['BCC', 'CAM4-Oslo', 'CAM5.1', 'GEOSCHEM', 'GISS-MATRIX', 'GISS-modelE', 'GMI', 'GOCART', 'HadGEM2', 'IMPACT-Umich', 'INCA', 'MPIHAM', 'NCAR-CAM3.5', 'OsloCTM2', 'SPRINTARS', 'mean_AeroCom2']
    Par.coords['mod_NO3_radeff'] = ['GEOSCHEM', 'GISS-MATRIX', 'GMI', 'HadGEM2', 'IMPACT-Umich', 'INCA', 'NCAR-CAM3.5', 'OsloCTM2', 'mean_AeroCom2']
    Par.coords['mod_SOA_radeff'] = ['CAM5.1', 'GEOSCHEM', 'IMPACT-Umich', 'MPIHAM', 'OsloCTM2', 'mean_AeroCom2']

    ## required value
    ## surface area of Earth
    A_Earth = 510072E9 # m2

    ## radiative efficiencies (as RF normalized to burden divived by Earth's area)
    ## sulfates
    Par['rf_SO4'] = xr.DataArray(1E12 / A_Earth * np.array(
        [-108, -173, -104, -123, -196, -307, -195, -238, -193, -113, -180, -125, -354, -192, -172, -185]),
        dims='mod_SO4_radeff', attrs={'units':'W m-2 Tg-1'})

    ## primary organic aerosols
    Par['rf_POA'] = xr.DataArray(1E12 / A_Earth * np.array(
        [-97, -118, -69, -95, -129, -76, -189, -144, -145, -141, -76, -41, -48, -165, -102, -113]),
        dims='mod_POA_radeff', attrs={'units':'W m-2 Tg-1'})
    
    ## black carbon
    Par['rf_BC'] = xr.DataArray(1E12 / A_Earth * np.array(
        [650, 1763, 2661, 1067, 2484, 1253, 1208, 874, 612, 1467, 1160, 1453, 1364, 2161, 1322, 1438]),
        dims='mod_BC_radeff', attrs={'units':'W m-2 Tg-1'})

    ## nitrate
    Par['rf_NO3'] = xr.DataArray(1E12 / A_Earth * np.array(
        [-136, -240, -103, -249, -155, -110, -91, -173, -166]),
        dims='mod_NO3_radeff', attrs={'units':'W m-2 Tg-1'})

    ## secondary organic aerosols
    Par['rf_SOA'] = xr.DataArray(1E12 / A_Earth * np.array(
        [-45, -45, -218, -139, -161, -122]),
        dims='mod_SOA_radeff', attrs={'units':'W m-2 Tg-1'})
    
    ## dust and salt aerosols set to zero
    Par['rf_dust'] = xr.DataArray(0., attrs={'units':'W m-2 Tg-1'})
    Par['rf_salt'] = xr.DataArray(0., attrs={'units':'W m-2 Tg-1'})

    ## return
    return Par


## adjustment coefficient for semi-direct effect
## taken from ACCMIP
def load_RFaer_semidirect(**useless):

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_BC_adjust'] = ['Boucher_2013', 'CSIRO', 'GISS', 'HadGEM2', 'ECHAM5', 'ECMWF']
    
    ## adjustment factor
    ## 1st value from IPCC AR5 (Boucher et al., 2013; doi:10.1017/CBO9781107415324.016)
    ## all others from ACCMIP (Lohmann et al., 2010; doi:10.5194/acp-10-3235-2010) (Figure 2; data provided by author)
    ## note: Boucher_2013 is best guess, and other values rescaled to it (models' mean = 0.111)
    Par['k_adj_BC'] = xr.DataArray((-0.1/0.6)/-0.111 * np.array(
        [-0.111, -0.37, -0.225, -0.13, 0.05, 0.12]),
        dims='mod_BC_adjust', attrs={'units':'1'})

    ## return
    return Par


## parameters for indirect effect
## taken from ACCMIP
def load_RFaer_indirect(**useless):

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_RFcloud_solub'] = ['Hansen_2005', 'Lamarque_2011']
    Par.coords['mod_RFcloud_erf'] = ['Boucher_2013', 'CSIRO-Mk3.6.0', 'GFDL-AM3', 'GISS-E2-R', 'HadGEM2', 'LMDzORINCA', 'MIROC-CHEM', 'NCAR-CAM5.1']
    Par.coords['mod_RFcloud_preind'] = ['low', 'median', 'high']

    ## temporary: ERF over 1850-2000 for the aerosol-cloud interaction in ACCMIP
    ## 1st value from IPCC AR5 (Boucher et al., 2013; doi:10.1017/CBO9781107415324.016)
    ## all others from ACCMIP (Shindell et al., 2013; doi:10.5194/acp-13-2939-2013) (Table 7; Aer ERF - Aer direct plus missing)
    ## note: Boucher_2013 is best guess, and other values rescaled to it (models' mean = -0.84)
    Par['ERF_cloud'] = xr.DataArray(-0.45/-0.84 * np.array(
        [-0.84, -0.99, -0.82, -0.61, -0.89, -0.21, -1.12, -1.22]),
        dims='mod_RFcloud_erf', attrs={'units':'W m-2'})

    ## temporary: load of soluble aerosols in 1850 and 2000 in ACCMIP
    ## note: values directly extracted from OSCAR v2 and copy-pasted here!
    ## /!\ WARNING: those values depend on the solubility factors p_sol_ defined elsewhere above!
    Par['AERsol_1850'] = xr.DataArray(
        [[14.426657, 16.418844, 8.222893, 9.29759, 27.672163, 18.48567, 5.0919113, 10.659031],
        [6.493315, 8.630909, 3.7277446, 7.1631308, 7.775085, 3.8228593, 3.495258, 5.2824764]],
        dims=['mod_RFcloud_solub', 'mod_RFcloud_erf'], attrs={'units':'Tg'})
    Par['AERsol_2000'] = xr.DataArray(
        [[17.000614, 18.901499, 10.227486, 12.710339, 29.372377, 20.307552, 6.6259317, 12.145359],
        [9.124497, 11.215797, 5.859241, 10.560671, 9.570654, 5.533892, 5.2857947, 6.660064]],
        dims=['mod_RFcloud_solub', 'mod_RFcloud_erf'], attrs={'units':'Tg'})

    ## final: intensity of cloud-aerosols interaction
    ## based on log formula
    Par['Phi_0'] = Par['ERF_cloud'] / np.log(Par['AERsol_2000'] / Par['AERsol_1850'])
    Par['Phi_0'].attrs['units'] = 'W m-2'
    del Par['ERF_cloud'], Par['AERsol_2000']

    ## final: preindustrial soluble aerosols load
    ## (Carslaw et al., 2013; doi:doi:10.1038/nature12674)
    factor = xr.DataArray([2., 1., 0.], dims='mod_RFcloud_preind')
    Par['AERsol_0'] = Par['AERsol_1850'] * np.exp(factor * (1.42 - 1.30) / Par['Phi_0'])
    del Par['AERsol_1850']

    ## return
    return Par


##==========================
## 3.3. Black carbon on snow
##==========================

## Reddy_2007 regions split fractions
def load_regions_Reddy_2007(mod_region, recalibrate=False, **useless):
    ## TODO v3.2: maybe? correct fractions to account only for land...
    
    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/regions_Reddy_2007__' + mod_region + '.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/regions_Reddy_2007__' + mod_region + '.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par = calib_regions_Reddy_2007(mod_region=mod_region)

    ## return
    return Par


## parameters for RF from BC deposition on snow
## partly based on ACCMIP
def load_RFbcsnow_all(**useless):

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_RFbcsnow_reg'] = ['Reddy_2007']
    Par.coords['reg_bcsnow'] = ['Globe', 'SAM', 'NAM', 'AFR', 'EUR', 'WCA', 'SAS', 'EAS', 'AUP', 'OCE']
    Par.coords['mod_RFbcsnow_radeff'] = ['Boucher_2013', 'CICERO-OsloCTM2', 'GFDL-AM3', 'GISS-E2-R', 'GISS-E2-R-TOMAS', 'HadGEM2', 'MIROC-CHEM', 'NCAR-CAM3.5', 'NCAR-CAM5.1']

    ## regional weighting coefficients
    ## (Reddy et al., 2007; doi:10.1029/2006GL028904) (Table 1)
    Par['w_reg_bcsnow'] = xr.DataArray(
        [[1., 1/6., 11/11., 1/10., 63/12., 2/3., 2/13., 17/43., 1/1., 2/1.]],
        dims=['mod_RFbcsnow_reg', 'reg_bcsnow'], attrs={'units':'1'})

    ## radiative forcing normalized to emissions
    ## 1st value from IPCC AR5 (Boucher et al., 2013; doi:10.1017/CBO9781107415324.016) (text & Table 7.1a)
    ## all others from ACCMIP (Lee et al., 2013; doi:10.5194/acp-13-2607-2013) (Table 3 & Fig. 15; data provided by author)
    ## note: Boucher_2013 is best guess, and other values rescaled to it (models' mean = 0.0146/(7.9-3.2))
    Par['rf_bcsnow'] = xr.DataArray((0.04/4.8) / (0.0146/(7.9-3.2)) * 
        np.array([0.0146, 0.0131, 0.0130, 0.0142, 0.0175, 0.0133, 0.0173, 0.0143, 0.0141]) /
        np.array([7.9-3.2, 7.8-3.1, 7.8-3.1, 8.8-4.0, 7.8-3.1, 7.8-3.1, 7.7-3.0, 7.8-3.1, 7.8-3.1]),
        dims='mod_RFbcsnow_radeff', attrs={'units':'W yr m-2 TgC-1'})

    ## return
    return Par


##=======================
## 3.4. Land cover change
##=======================

## land-cover change albedo parameters
def load_RFlcc_all(mod_region, recalibrate=False, **useless):
    ## TODO v3.2: redo? remove LUH1 pixels limitation
    
    ## initialization
    Par = xr.Dataset()

    ## average upward transmittance of the atmosphere
    ## (Lenton & Vaughan, 2009; doi:10.5194/acp-9-5539-2009)
    Par['p_trans'] = xr.DataArray(0.854, attrs={'units':'1'})

    ## other parameters
    ## /!\ compiled from various land cover, albedo and radiative flux climatologies
    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/albedo_all__' + mod_region + '.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/albedo_all__' + mod_region + '.nc') as TMP: Par2 = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par2 = calib_albedo_all(mod_region=mod_region)

    ## return
    return xr.merge([Par, Par2])


##==================
## 3.5. Specific RFs
##==================

## warming efficacies
def load_RF_warmeff(**useless):

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_RFbcsnow_warmeff'] = ['median', 'low', 'high']
    Par.coords['mod_RFlcc_warmeff'] = ['Hansen_2005', 'Davin_2007', 'Davin_2010', 'Jones_2013']

    ## volcano forcing
    ## (Gregory et al., 2016; doi:10.1007/s00382-016-3055-1)
    Par['w_warm_volc'] = xr.DataArray(0.6, attrs={'units':'1'})

    ## black carbon on snow
    ## (Boucher et al., 2013; doi:10.1017/CBO9781107415324.016) (Sect. 7.5.2.3)
    Par['w_warm_bcsnow'] = xr.DataArray([3., 2., 4.], dims='mod_RFbcsnow_warmeff', attrs={'units':'1'})

    ## land-cover change
    ## (Bright et al., 2015; doi:10.1111/gcb.12951) (Table 7)
    Par['w_warm_lcc'] = xr.DataArray([1.02, 0.50, 0.78, 0.79], dims='mod_RFlcc_warmeff', attrs={'units':'1'})

    ## return
    return Par


## atmospheric RF fractions
def load_RF_atmfrac(**useless):

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_Pg_radfact'] = ['Andrews_2010', 'Kvalevag_2013']

    ## atmospheric fractions factors from
    ## (Andrews et al., 2010; doi:10.1029/2010GL043991) (Table 3)
    ## (Kvalevag et al., 2013; doi:10.1002/grl.50318) (Table 2; highest perturbation & zero if not given)
    Par['p_atm_CO2'] = xr.DataArray([0.8, 2.0/3.6], dims='mod_Pg_radfact', attrs={'units':'1'})
    Par['p_atm_nonCO2'] = xr.DataArray([0.5, 0.5/1.8], dims='mod_Pg_radfact', attrs={'units':'1'})
    Par['p_atm_O3t'] = xr.DataArray([-0.3, 0.0], dims='mod_Pg_radfact', attrs={'units':'1'})
    Par['p_atm_strat'] = xr.DataArray([0.0, 0.0], dims='mod_Pg_radfact', attrs={'units':'1'}) # assumed
    Par['p_atm_scatter'] = xr.DataArray([0.0, 0.6/-1.4], dims='mod_Pg_radfact', attrs={'units':'1'})
    Par['p_atm_absorb'] = xr.DataArray([2.5, 3.1/0.5], dims='mod_Pg_radfact', attrs={'units':'1'})
    Par['p_atm_cloud'] = xr.DataArray([0.0, 0.0], dims='mod_Pg_radfact', attrs={'units':'1'}) # assumed
    Par['p_atm_alb'] = xr.DataArray([0.0, 0.0], dims='mod_Pg_radfact', attrs={'units':'1'})
    Par['p_atm_solar'] = xr.DataArray([0.2, 0.4/3.5], dims='mod_Pg_radfact', attrs={'units':'1'})

    ## return
    return Par


##################################################
##   4. CLIMATE
##################################################

##=================
## 4.1. Temperature
##=================

## temperature response
## calibrated on CMIP5 models
def load_temp_CMIP5(mod_region, recalibrate=False, **useless):

    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/temp_CMIP5__' + mod_region + '.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/temp_CMIP5__' + mod_region + '.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par = calib_temp_CMIP5(mod_region=mod_region)

    ## return
    return Par


##===================
## 4.2. Precipitation
##===================

## precipitation response
## calibrated on CMIP5 models
def load_prec_CMIP5(mod_region, recalibrate=False, **useless):

    ## load from existing file
    if os.path.isfile('input_data/parameters/from_OSCARv2/prec_CMIP5__' + mod_region + '.nc') and not recalibrate:
        with xr.open_dataset('input_data/parameters/from_OSCARv2/prec_CMIP5__' + mod_region + '.nc') as TMP: Par = TMP.load()
    ## otherwise, launch calibration
    else:
        raise RuntimeError('embedded calibration not available yet')
        #Par = calib_prec_CMIP5(mod_region=mod_region)

    ## return
    return Par


##========================
## 4.3. Ocean heat content
##========================

## OHC parameters
def load_OHC_all(**useless):
    ## TODO v3.2: use alternative formula

    ## initialization
    Par = xr.Dataset()
    
    ## fraction of energy not used to heat land, atmosphere or melt ice
    ## (Otto et al., 2013; doi:10.1038/ngeo1836) (SI)
    Par['p_ohc'] = xr.DataArray(0.94, attrs={'units':'1'})

    ## return
    return Par


##################################################
##   5. IMPACTS
##################################################

##===================
## 5.1. Acidification
##===================

## surface acidification parameters
def load_pH_all(**useless):
    ## TODO v3.2: remove log formula?

    ## initialization
    Par = xr.Dataset()
    Par.coords['mod_pH_fct'] = ['log', 'poly']

    ## option to choose pH functional form
    Par['pH_is_Log'] = xr.DataArray([True, False], dims='mod_pH_fct')

    ## return
    return Par


##################################################
##   Z. WRAPPER
##################################################

## wrapping function
def load_all_param(mod_region, recalibrate=False):
    '''
    Wrapper function to load all primary parameters.
    
    Input:
    ------
    mod_region (str)        regional aggregation name       

    Output:
    -------
    Par (xr.Dataset)        merged dataset

    Options:
    --------
    recalibrate (bool)      whether to recalibrate all possible parameters;
                            WARNING: currently not working;
                            default = False
    '''

    print('loading primary parameters')

    ## list of loading fuctions
    load_list = [load_ocean_struct, load_ocean_chem, load_ocean_CMIP5,
        load_land_misc, load_land_TRENDYv7, load_land_CMIP5, load_land_wooduse, load_land_GFED3,
        load_permafrost_all, 
        load_wetlands_WETCHIMP,
        load_atmosphere_misc, load_atmosphere_CCMVal2, load_atmochem_adhoc,
        load_CH4_lifetime, load_OH_response,
        load_N2O_lifetime, load_hv_response,
        load_halo_lifetime,
        load_regions_HTAP, load_O3t_regional, load_O3t_response,
        load_ODS_all, load_O3s_response,
        load_AER_regional, load_AER_atmoload, load_AER_solub,
        load_RFghg_all, 
        load_RFozo_all, load_RFaer_direct, load_RFaer_semidirect, load_RFaer_indirect,
        load_regions_Reddy_2007, load_RFbcsnow_all,
        load_RFlcc_all,
        load_RF_warmeff, load_RF_atmfrac,
        load_temp_CMIP5,
        load_prec_CMIP5,
        load_OHC_all,
        load_pH_all]
    
    ## return all
    return xr.merge([load(mod_region=mod_region, recalibrate=recalibrate) for load in load_list])

