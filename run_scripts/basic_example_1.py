##################################################
##   EXAMPLE 1a
##################################################

## GOAL: use the wrapper function to run historical simulations in a Monte Carlo setup

## being in the OSCAR folder, import the wrapper function
from core_fct.fct_wrap import run_model

## import the OSCAR model
from core_fct.fct_process import OSCAR

## run the model and save output
## note: the creation of secondary drivers is rather long!
## note: requires defining a 'mod_region' and the number of MC elements 'nMC'
RUN_a = run_model(OSCAR, (1750, 1750, 2014), mod_region='Houghton_2001', nMC=200, output=True)
OUT = RUN_a['Out_hist']

## import matplotlib
import matplotlib.pyplot as plt

## plot global temperature change
plt.plot(OUT.year, OUT.D_Tg.mean('config'), color='k', lw=2)
plt.fill_between(OUT.year, OUT.D_Tg.mean('config') - OUT.D_Tg.std('config'), OUT.D_Tg.mean('config') + OUT.D_Tg.std('config'), color='k', alpha=0.5)


##################################################
##   EXAMPLE 1b
##################################################

## GOAL: use the wrapper function to run the RCP scenarios in a Monte Carlo setup

## same imports as before
from core_fct.fct_wrap import run_model
from core_fct.fct_process import OSCAR

## import function to load scenarios provided with OSCAR
from core_fct.fct_loadD import load_all_scen

## load RCP-specific drivers
For_rcp = load_all_scen(mod_region='RCP_5reg', group_scen=['RCP2.6', 'RCP4.5', 'RCP6.0', 'RCP8.5'], LCC='gross')

## define which variables will be kept (in addition to state variables)
var_keep = ['D_Eluc', 'D_Focean', 'D_Fland', 'D_Epf'] + ['tau_CH4', 'tau_N2O'] + [var for var in list(OSCAR._processes.keys()) if 'RF' in var]

## run the model and save output
## note: increase number of sub-timesteps 'nt' for better stability (also increases run time)
RUN_b = run_model(OSCAR, (1750, 1750, 2014, 2100), For1=For_rcp, mod_region='RCP_5reg', nMC=200, output=True, var_keep=var_keep, nt=4)

## simple plot function
## note: works only with timeseries of global (i.e. one-dimensional) variables
import matplotlib.pyplot as plt
def f_plot(VAR, RUN=RUN_b):
    plt.figure()
    plt.plot(RUN['Out_hist'].year, RUN['Out_hist'][VAR].mean('config'), color='k', lw=2, label='hist')
    plt.fill_between(RUN['Out_hist'].year, RUN['Out_hist'][VAR].mean('config') - RUN['Out_hist'][VAR].std('config'), 
        RUN['Out_hist'][VAR].mean('config') + RUN['Out_hist'][VAR].std('config'), color='k', alpha=0.5)
    for scen in RUN['Out_scen'].scen.values:
        plt.plot(RUN['Out_scen'].year, RUN['Out_scen'][VAR].sel(scen=scen).mean('config'), lw=2, label=scen)
        plt.fill_between(RUN['Out_scen'].year, RUN['Out_scen'][VAR].sel(scen=scen).mean('config') - RUN['Out_scen'][VAR].sel(scen=scen).std('config'), 
            RUN['Out_scen'][VAR].sel(scen=scen).mean('config') + RUN['Out_scen'][VAR].sel(scen=scen).std('config'), alpha=0.5)
    plt.title(VAR + ' (' + RUN['Out_scen'][VAR].units + ')')

## plot some results
f_plot('D_Tg')
f_plot('D_CO2')

