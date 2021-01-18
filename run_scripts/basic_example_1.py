##################################################
##   EXAMPLE 1
##################################################

## GOAL: run a historical simulation, followed by the SSPs, in a Monte Carlo setup

## being  in the OSCAR folder, import the model
from core_fct.mod_process import OSCAR

## load primary parameters
## note: requires choosing a 'mod_region' 
from core_fct.fct_loadP import load_all_param
Par0 = load_all_param(mod_region='RCP_5reg')

## generate Monte Carlo set of parameters
## note: requires choosing 'nMC' (here, 200 elements)
from core_fct.fct_genMC import generate_config
Par = generate_config(Par0, nMC=200)

## import the historical drivers from the provided script
## note: SSPs are concentration-driven scenarios, so atmospheric GHGs are prescribed (and some emissions are missing)
from run_scripts.get_SSP_drivers import For_hist

## move parameters from For to Par
## note: some parameters may come from driving datasets (normally, Aland_0)
import xarray as xr
Par = xr.merge([Par, For_hist.drop([VAR for VAR in For_hist if 'year' in For_hist[VAR].dims])])
For_hist = For_hist.drop([VAR for VAR in For_hist if 'year' not in For_hist[VAR].dims])

## run the model
## note: Ini is set to None to ask the model to create a zero-valued initial state
Out_hist = OSCAR(Ini=None, Par=Par, For=For_hist)

## save output (if asked)
if False: Out_hist.to_netcdf('results/Out_hist_test.nc')

## import matplotlib
import matplotlib.pyplot as plt

## plot global temperature change
plt.plot(Out_hist.year, Out_hist.D_Tg.mean('config'), color='k', lw=2)
plt.fill_between(Out_hist.year, Out_hist.D_Tg.mean('config') - Out_hist.D_Tg.std('config'), Out_hist.D_Tg.mean('config') + Out_hist.D_Tg.std('config'), color='k', alpha=0.5)

##--------------------

## import the scenario drivers from the provided script (but take only up to 2100)
from run_scripts.get_SSP_drivers import For_scen
For_scen = For_scen.sel(year=slice(None, 2100))

## create initial state
Ini = Out_hist.isel(year=-1, drop=True)

## define extra which variables to be kept (in addition to state variables)
var_keep = ['D_Eluc', 'D_Focean', 'D_Fland', 'D_Epf'] + ['tau_CH4', 'tau_N2O'] + [var for var in list(OSCAR._processes.keys()) if 'RF' in var]

## run the model and save output
## note: increase number of sub-timesteps 'nt' for better stability (also increases run time)
Out_scen = OSCAR(Ini=Ini, Par=Par, For=For_scen, var_keep=var_keep, nt=4)
if False: Out_scen.to_netcdf('results/Out_scen_test.nc')

## simple plot function
## note: works only with timeseries of global (i.e. one-dimensional) variables
import matplotlib.pyplot as plt
def f_plot(VAR):
    plt.figure()
    if VAR in Out_hist:
        plt.plot(Out_hist.year, Out_hist[VAR].mean('config'), color='k', lw=2, label='hist')
        plt.fill_between(Out_hist.year, Out_hist[VAR].mean('config') - Out_hist[VAR].std('config'), 
            Out_hist[VAR].mean('config') + Out_hist[VAR].std('config'), color='k', alpha=0.5)
    for scen in Out_scen.scen.values:
        plt.plot(Out_scen.year, Out_scen[VAR].sel(scen=scen).mean('config'), lw=2, label=scen)
        plt.fill_between(Out_scen.year, Out_scen[VAR].sel(scen=scen).mean('config') - Out_scen[VAR].sel(scen=scen).std('config'), 
            Out_scen[VAR].sel(scen=scen).mean('config') + Out_scen[VAR].sel(scen=scen).std('config'), alpha=0.5)
    plt.title(VAR + ' (' + Out_scen[VAR].units + ')')

## plot some results
f_plot('D_Tg')
f_plot('D_Focean')

