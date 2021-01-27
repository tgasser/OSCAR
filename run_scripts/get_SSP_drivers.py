#################################################
## GENERATE SSP DRIVERS
#################################################

## the following code will create historical and scenario drivers to run the 8 marker SSPs of the AR6

## general imports
import os
import numpy as np
import xarray as xr

## specific imports
from core_fct.fct_misc import aggreg_region, load_data
from core_fct.fct_loadD import load_landuse_hist, load_RFdrivers_hist, load_landuse_scen, load_RFdrivers_scen
from core_fct.fct_loadP import load_all_param
from core_fct.fct_genMC import generate_config


##################################################
##################################################

## OPTIONS
folder_out = 'results/SSP_drivers/'
mod_region = 'RCP_5reg'

## create output folder
if not os.path.isdir(folder_out):
    os.mkdir(folder_out)

## get primary parameters
Par0 = load_all_param(mod_region)


##################################################
## 1. HISTORICAL
##################################################

## load emissions
For_hist = aggreg_region(load_data('emissions_CEDS'), mod_region).sum('sector')
For_hist = For_hist.rename({'E_CO2':'Eff'})

## load LULCC
TMP = load_landuse_hist(mod_region, ['LUH2'], LCC='gross').sel(data_LULCC='LUH2', drop=True)
For_hist = xr.merge([For_hist, TMP])

## load RF drivers
TMP = load_RFdrivers_hist().sel(data_RF_contr='ICAO', data_RF_solar='CMIP6', data_RF_volc='CMIP6', drop=True)
For_hist = xr.merge([For_hist, TMP])
For_hist['RF_contr'] *= 0

## load atmospheric concentrations
TMP = load_data('concentrations_CMIP6').sel(region='Globe', drop=True)
for VAR in ['CO2', 'CH4', 'N2O', 'Xhalo']:
    For_hist['D_'+VAR] = TMP[VAR] - Par0[VAR+'_0']

## missing drivers
For_hist['Eluc'] = 0*For_hist['Eff']
For_hist['E_N2O'] = 0*For_hist['Eff']
For_hist['E_Xhalo'] = 0*For_hist['D_Xhalo']

## format and save
For_hist = For_hist.sel(year=slice(1750, 2014)).fillna(0.)
if __name__ == '__main__':
    For_hist.to_netcdf(folder_out + 'For_hist_SSP.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in For_hist})


##################################################
## 2. SCENARIOS
##################################################

## load emissions
For_scen = aggreg_region(load_data('emissions_ScenarioMIP'), mod_region, old_axis='region', dataset='ScenarioMIP')
For_scen = xr.concat([For_scen, For_scen.sel(year=2100).assign_coords(year=2300)], dim='year')

## correction: non-CO2 emissions from FF&I decrease to 0, while those from LUC stay constant
sectors_LU = ['Forest Burning', 'Grassland Burning', 'Peat Burning'] + ['Agricultural Waste Burning', 'Agriculture']
for var in For_scen.variables:
    if (var not in For_scen.coords)  and  (var not in ['E_N2O', 'Eluc']):
        For_scen[var].loc[{'sector':[sec for sec in For_scen.sector.values if sec not in sectors_LU], 'year':2300}] = 0.
For_scen = For_scen.interp(year=np.arange(2015, 2300+1)).sum('sector')

## specific treatment for Eff. Does not follow O'Neill et al, 2016, but updates of extensions from the poster of Meinshausen 2019 (ScenarioForum).
dico_dates_end_plateau = {'SSP1-1.9':2140, 'SSP1-2.6':2140, 'SSP2-4.5':2100, 'SSP3-7.0':2100, 'SSP3-7.0-LowNTCF':2100, 'SSP4-3.4':2140, 'SSP4-6.0':2100, 'SSP5-3.4-OS':2140, 'SSP5-8.5':2100}
dico_dates_end_ramp = { 'SSP1-1.9':2190, 'SSP1-2.6':2190, 'SSP2-4.5':2250, 'SSP3-7.0':2250, 'SSP3-7.0-LowNTCF':2250, 'SSP4-3.4':2190, 'SSP4-6.0':2250, 'SSP5-3.4-OS':2170, 'SSP5-8.5':2250 }
for scen in For_scen.scen.values:
    For_scen['Eff'].loc[dict(scen=scen,year=range(2100,dico_dates_end_plateau[scen]+1))] = For_scen['Eff'].sel(scen=scen,year=2100) # constant over 2100-2140
    For_scen['Eff'].loc[dict(scen=scen,year=np.arange(dico_dates_end_plateau[scen],dico_dates_end_ramp[scen]+1))] = np.linspace( For_scen['Eff'].sel(scen=scen,year=2100) , 0.*For_scen['Eff'].sel(scen=scen,year=2100) , dico_dates_end_ramp[scen]-dico_dates_end_plateau[scen]+1)  # linear ramp over 2140-2190
    For_scen['Eff'].loc[dict(scen=scen,year=np.arange(dico_dates_end_ramp[scen]+1,2300+1))] = 0. # 0 afterwards

## specific treatment for Eluc. Does not follow O'Neill et al, 2016, but updates of extensions from the poster of Meinshausen 2019 (ScenarioForum).
For_scen['Eluc'].loc[dict(year=np.arange(2100,2150+1))] = np.linspace( For_scen['Eluc'].sel(year=2100) , 0.*For_scen['Eluc'].sel(year=2100) , 2150-2100+1)  # linear ramp over 2100-2150
For_scen['Eluc'].loc[dict(year=np.arange(2150+1,2300+1))] = 0. # 0 afterwards

## load LULCC
TMP = load_landuse_scen(mod_region, ['LUH2'], LCC='gross').rename({'scen_LULCC':'scen'})
TMP.coords['scen'] = [scen + '-OS' * (scen == 'SSP5-3.4') for scen in TMP['scen'].values]
TMP = TMP.where(TMP.year < 2100).dropna('year')
TMP = xr.concat([TMP, TMP.sel(year=2099).assign_coords(year=2300)], dim='year')
TMP = TMP.interp(year=np.arange(850, 2300+1))
TMP['d_Acover'] = TMP['d_Acover'].where(TMP.year <= 2099, 0)
For_scen = xr.merge([For_scen, TMP])

## load RF drivers
TMP = load_RFdrivers_scen().sel(scen_RF_contr='CMIP5', scen_RF_solar='CMIP6', scen_RF_volc='CMIP6', drop=True)
for VAR in TMP:
    years = (For_hist[VAR].dropna('year').year + TMP[VAR].dropna('year').year).year
    For_scen[VAR] = For_hist[VAR].sel(year=years).mean('year') / TMP[VAR].sel(year=years).mean('year') * TMP[VAR]

## load atmospheric concentrations
TMP = load_data('concentrations_ScenarioMIP').drop('Xhalo_eq')
for VAR in ['CO2', 'CH4', 'N2O', 'Xhalo']:
    For_scen['D_'+VAR] = TMP[VAR] - Par0[VAR+'_0']

## missing drivers
For_scen['E_Xhalo'] = 0*For_scen['D_Xhalo']

## format and save
For_scen = For_scen.sel(year=slice(2015, 2300))
For_scen = xr.concat([For_hist.sel(year=2014).drop('Aland_0'), For_scen], dim='year')
if __name__ == '__main__':
    For_scen.to_netcdf(folder_out + 'For_scen_SSP.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in For_scen})

