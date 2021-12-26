##################################################
##   EXAMPLE
##################################################


## GOAL: run historical simulations to detect the diffference between OSCAR with and without nature emission module

## being in the OSCAR folder, import the wrapper function
from core_fct.fct_wrap import run_model

## import the OSCAR model
from core_fct.mod_process import OSCAR

## import matplotlib
import matplotlib.pyplot as plt

nMC = 500

## run the model and save output
## note: the creation of secondary drivers is rather long!
## note: requires defining a 'mod_region' and the number of MC elements 'nMC'
RUN_a = run_model(OSCAR, (1750, 1750, 2014), mod_region='AR6_22reg', nMC=nMC, output=True,
                  nature_aerosols=True)
OUT_a = RUN_a['Out_hist']
## plot global temperature change
plt.plot(OUT_a.year, OUT_a.D_Tg.mean('config'), color='b', lw=2)
plt.fill_between(OUT_a.year, OUT_a.D_Tg.mean('config') - OUT_a.D_Tg.std('config'), OUT_a.D_Tg.mean('config') + OUT_a.D_Tg.std('config'), color='b', alpha=0.5)

# same as above, but without nature emission module

RUN_b = run_model(OSCAR, (1750, 1750, 2014), mod_region='AR6_22reg', nMC=nMC, output=True,
                  nature_aerosols=False)
OUT_b = RUN_b['Out_hist']
## plot global temperature change
plt.plot(OUT_b.year, OUT_b.D_Tg.mean('config'), color='k', lw=2)
plt.fill_between(OUT_b.year, OUT_b.D_Tg.mean('config') - OUT_b.D_Tg.std('config'), OUT_b.D_Tg.mean('config') + OUT_b.D_Tg.std('config'), color='k', alpha=0.5)


