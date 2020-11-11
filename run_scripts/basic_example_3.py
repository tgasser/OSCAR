#################################################
## EXAMPLE 3
#################################################

## GOAL: simulate historical land-use CO2 emissions

## imports (generic)
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

## imports (specific)
from core_fct.fct_loadP import load_all_param
from core_fct.fct_genMC import generate_config
from core_fct.fct_loadD import load_landuse_hist
from core_fct.fct_misc import load_data, aggreg_region

## import submodel that contains only the land carbon cycle
from core_fct.mod_process import OSCAR_landC

## set run options
mod_region = 'AR6_10reg'
ind0, ind1 = 1750, 2018

## load primary parameters and generate an averaged configuration 
## note: some 'mod_' options are prescribed to mimic ORCHIDEE-CNP/IPSL-CM5A-LR
Par0 = load_all_param(mod_region)
Par_mc = generate_config(Par0, 1000, fixed={'mod_Fland_preind':'ORCHIDEE-CNP', 'mod_Efire_preind':'off_', 
    'mod_Eharv_preind':'ORCHIDEE-CNP', 'mod_Egraz_preind':'off_','mod_Eluc_agb':'ORCHIDEE-CNP',
    'mod_Fland_trans':'IPSL-CM5A-LR', 'mod_Efire_trans':'IPSL-CM5A-LR'})
Par = Par_mc.mean('config')

## load raw drivers
For_luc = load_landuse_hist(mod_region, LCC='gross')
For_clim = aggreg_region(load_data('CRU-TS'), mod_region, weight_var={'Tl':'area', 'Pl':'area'}).drop('area')
For_co2 = xr.concat([load_data('concentrations_AR5').CO2.interp({'year':np.arange(1750, 1959)}), 
    load_data('NOAA-ESRL').CO2.sel(site='Mauna Loa', drop=True)], dim='year').to_dataset(name='CO2')

## format drivers
For = xr.merge([For_luc, For_clim, For_co2]).sel(year=slice(ind0, ind1))
For['D_CO2'] = For.CO2 - Par.CO2_0
For['D_Tl'] = For.Tl - For.Tl.sel(year=slice(1901, 1930)).mean('year')
For['D_Pl'] = For.Pl - For.Pl.sel(year=slice(1901, 1930)).mean('year')
For = For.drop(['CO2', 'Tl', 'Pl'])
For = For.fillna(0).sel(year=slice(1900)).combine_first(For)

## swap preindustrial land-cover
Par['Aland_0'] = For.Aland_0
For = For.drop('Aland_0')

## run OSCAR land carbon cycle
OUT = OSCAR_landC(Ini=None, Par=Par, For=For)

## DAMNED! forgot to specify which diagnostic variables to keep
## SOLUTION: recursively call the corresponding process of the model using the computed prognostic variables as input
OUT2 = xr.merge([OSCAR_landC[VAR](OUT, Par, For, recursive=True).to_dataset(name=VAR) for VAR in ['D_Fland', 'D_Eluc']])

## simple plot
plt.figure()
for data in OUT2.data_LULCC.values:
    plt.plot(OUT2.year, OUT2.D_Eluc.sel(data_LULCC=data), label=data, alpha=0.7)
plt.legend(loc=0)
plt.title(r'$\Delta E_\mathrm{luc}$ (' + OSCAR_landC['D_Eluc'].units.replace('-1', '$^{-1}$') +')')

