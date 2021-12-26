##################################################
##   EXAMPLE 1a
##################################################

## GOAL: use the wrapper function to run historical simulations in a Monte Carlo setup

## being in the OSCAR folder, import the wrapper function
from core_fct.fct_wrap import run_model

## import the OSCAR model
from core_fct.mod_process import OSCAR

## run the model and save output
## note: the creation of secondary drivers is rather long!
## note: requires defining a 'mod_region' and the number of MC elements 'nMC'
RUN_a = run_model(OSCAR, (1750, 1750, 2014), mod_region='AR6_22reg', nMC=20, output=True,
                  nature_aerosols=True)
OUT = RUN_a['Out_hist']

## import matplotlib
import matplotlib.pyplot as plt

## plot global temperature change
plt.plot(OUT.year, OUT.D_Tg.mean('config'), color='k', lw=2)
plt.fill_between(OUT.year, OUT.D_Tg.mean('config') - OUT.D_Tg.std('config'), OUT.D_Tg.mean('config') + OUT.D_Tg.std('config'), color='k', alpha=0.5)

