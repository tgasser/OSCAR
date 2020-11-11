##################################################
##   EXAMPLE 2
##################################################

## GOAL: look into some feedback loops, using an averaged parameterisation and the 'OSCAR' object

## imports
import xarray as xr
import matplotlib.pyplot as plt
from core_fct.mod_process import OSCAR
from core_fct.fct_loadD import load_all_hist
from core_fct.fct_loadP import load_all_param
from core_fct.fct_genMC import generate_config

## choose a regional aggregation
mod_region = 'AR6_5reg'

## load primary parameters and drivers
Par0 = load_all_param(mod_region)
For0 = load_all_hist(mod_region)

## generate 1000 configurations and average them
Par_mc = generate_config(Par0, 1000)
Par_av = Par_mc.mean('config')

## choose best guess drivers
data_bg = ['CEDS', 'PRIMAP', 'LUH2', 'AR5']
For_bg = For0.copy(deep=True)
for data_dim in [dim for dim in For0.dims if 'data_' in dim]:
    data_val = [data for data in data_bg if data in For0[data_dim]][0]
    For_bg = For_bg.sel({data_dim:data_val}, drop=True)
For_bg['E_CH4'] = For0['E_CH4'].sel(data_E_CH4='PRIMAP', drop=True)

## select years and extrapolate missing drivers
For_bg = For_bg.sel(year=slice(1750, 2011))
for VAR in ['E_CH4', 'E_N2O']:
    For_bg[VAR].loc[{'year':1750}] = 0
    For_bg[VAR] = For_bg[VAR].dropna('year').interp({'year':range(1750, 2011+1)})

## final formating of drivers
For_bg = For_bg.fillna(0)
For_bg = For_bg.rename({'d_Agross':'d_Acover'}).drop('d_Anet')

## swap preindustrial land-cover from For to Par
Par_av['Aland_0'] = For_bg.Aland_0
For_bg = For_bg.drop('Aland_0')

## create initial state
## note: can be set to None to ask the model to create zero-valued xarrays by itself
Ini = None

## run control simulation
## note: OSCARv3 is comparatively slow when running with only one set of drivers and parameters, because of all the xarray packing/unpacking; better use a MC setup!
Out_ctrl = OSCAR(Ini, Par_av, For_bg, var_keep=['RF'])

## run simulation in which atmospheric CO2 / CH4 /climate remains at preindustrial level
## note: For_bg.RF_contr is used because it already has the right dimensions, but the extra driver can be created from scratch
Out_noCO2 = OSCAR(Ini, Par_av, xr.merge([For_bg, 0*For_bg.RF_contr.to_dataset(name='D_CO2')]), var_keep=['RF'])
Out_noCH4 = OSCAR(Ini, Par_av, xr.merge([For_bg, 0*For_bg.RF_contr.to_dataset(name='D_CH4')]), var_keep=['RF'])
Out_noClim = OSCAR(Ini, Par_av, xr.merge([For_bg, 0*For_bg.RF_contr.to_dataset(name='D_Tg'), 0*For_bg.RF_contr.to_dataset(name='D_Pg')]), var_keep=['RF'])

## plot results in terms of RF
plt.figure()
plt.subplot(1, 2, 1)
plt.plot(Out_ctrl.year, Out_ctrl.RF, color='k', label='ctrl')
plt.plot(Out_noCO2.year, Out_noCO2.RF, label='noCO2')
plt.plot(Out_noCH4.year, Out_noCH4.RF, label='noCH4')
plt.plot(Out_noClim.year, Out_noClim.RF, label='noClim')
plt.legend(loc=0)
plt.title('RF (' + Out_ctrl.RF.units + ')')
plt.subplot(1, 2, 2)
plt.plot(Out_noCO2.year, Out_ctrl.RF - Out_noCO2.RF)
plt.plot(Out_noCH4.year, Out_ctrl.RF - Out_noCH4.RF)
plt.plot(Out_noClim.year, Out_ctrl.RF - Out_noClim.RF)
plt.title('RF_ctrl - RF (' + Out_ctrl.RF.units + ')')

