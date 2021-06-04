#################################################
## QUICK CHECK OF LUC STRUCTURE
#################################################

## GOAL: test all possible structures of the LUC module

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
from core_fct.mod_process import OSCAR_landC as OSCAR

## import functions for alternative LUC structures and apply
from core_fct.fct_process_alt import split_LUC, full_LUC, lite_LUC, cut_LUC
OSCAR_split = split_LUC(OSCAR)
OSCAR_full = full_LUC(OSCAR)
OSCAR_lite = lite_LUC(OSCAR)
OSCAR_cut = cut_LUC(OSCAR)

## set run options
mod_region = 'Houghton_2017'
ind0, ind1 = 1750, 2018

## load primary parameters and generate MC ensemble
Par0 = load_all_param(mod_region)
Par = generate_config(Par0, 1000)

## load raw drivers
For_luc = load_landuse_hist(mod_region, LCC='gross').sel(data_LULCC='LUH2', drop=True)
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

## create drivers for lite version
For_alt = For.copy(deep=True)
For_alt['d_Ashift'] = For_alt['d_Ashift'].sum('bio_from', min_count=1).rename({'bio_to':'bio_land'})
For_alt['d_Again'] = For_alt['d_Acover'].sum('bio_from', min_count=1).rename({'bio_to':'bio_land'})
For_alt['d_Aloss'] = For_alt['d_Acover'].sum('bio_to', min_count=1).rename({'bio_from':'bio_land'})
For_alt = For_alt.drop(['d_Acover', 'bio_from', 'bio_to'])

## run OSCAR land carbon cycle
OUT = OSCAR(Ini=None, Par=Par, For=For, var_keep=['D_Fland', 'D_Eluc'])
OUT_split = OSCAR_split(Ini=None, Par=Par, For=For, var_keep=['D_Fland', 'D_Eluc'])
OUT_full = OSCAR_full(Ini=None, Par=Par, For=For, var_keep=['D_Fland', 'D_Eluc'])
OUT_lite = OSCAR_lite(Ini=None, Par=Par, For=For_alt, var_keep=['D_Fland', 'D_Eluc'])

## create exogenous Eluc forcing for cut version and run it
For_Eluc = OUT.D_Eluc.assign_coords(reg_land=0).expand_dims('reg_land', 1).to_dataset(name='Eluc')
OUT_cut = OSCAR_cut(Ini=None, Par=Par, For=xr.merge([For, For_Eluc]), var_keep=['D_Fland', 'D_Eluc'])

## simple plot
plt.figure()
plt.subplot(1,2,1)
plt.plot(OUT.year, OUT.D_Eluc.mean('config'), label='base', alpha=0.7)
plt.plot(OUT_split.year, OUT_split.D_Eluc.mean('config'), label='split', alpha=0.7, lw=0.5, ls='-.',  marker='x', ms=5)
plt.plot(OUT_full.year, OUT_full.D_Eluc.mean('config'), label='full', alpha=0.7, lw=0.5, marker='+', ms=5)
plt.plot(OUT_lite.year, OUT_lite.D_Eluc.mean('config'), label='lite', alpha=0.7, ls='--')
plt.plot(OUT_cut.year, OUT_cut.D_Eluc.mean('config'), label='cut', alpha=0.7, ls=':', marker='o', mfc='none', ms=3)
plt.title(r'$\Delta E_\mathrm{luc}$ (' + OSCAR['D_Eluc'].units.replace('-1', '$^{-1}$') +')')
plt.subplot(1,2,2)
plt.plot(OUT.year, OUT.D_Fland.mean('config'), label='base', alpha=0.7)
plt.plot(OUT_split.year, OUT_split.D_Fland.mean('config'), label='split', alpha=0.7, lw=0.5, ls='-.', marker='x', ms=5)
plt.plot(OUT_full.year, OUT_full.D_Fland.mean('config'), label='full', alpha=0.7, lw=0.5, marker='+', ms=5)
plt.plot(OUT_lite.year, OUT_lite.D_Fland.mean('config'), label='lite', alpha=0.7, ls='--')
plt.plot(OUT_cut.year, OUT_cut.D_Fland.mean('config'), label='cut', alpha=0.7, lw=0.5, ls=':', marker='o', mfc='none', ms=3)
plt.title(r'$\Delta F_\mathrm{land}$ (' + OSCAR['D_Fland'].units.replace('-1', '$^{-1}$') +')')
plt.legend(loc=0)

