#################################################
## GENERATE AR6-LIKE PARAMETERS
#################################################

## the following code will generate a set of prior parameters similar to the one used for the AR6

from core_fct.fct_loadP import load_all_param
from core_fct.fct_genMC import generate_config, adjust_config

## shift & noise
shift = {'e_ohu': 0.28}
noise = {'e_ohu': ('normal', 0.25),
    'rf_CO2': ('normal', 0.10), 'rf_CH4': ('normal', 0.15), 'rf_N2O': ('normal', 0.10), 'k_rf_H2Os': ('normal', 0.40), 'rf_Xhalo': ('normal', 0.10), 
    'p_trans': ('normal', 0.50), 'p_ohc': ('normal', 0.01),
    'lambda_0': ('lognormal', 0.8), 'Th_g': ('lognormal', 0.8), 'Th_d': ('lognormal', 0.8), 'th_0': ('lognormal', 0.8),
    'rf_bcsnow': ('lognormal', 0.5), 'Phi_0': ('lognormal', 0.5),  'k_adj_BC': ('lognormal', 0.5), 'rf_O3t': ('lognormal', 0.5), 'rf_O3s': ('lognormal', 0.5), 
    'rf_SO4': ('lognormal', 0.5), 'rf_POA': ('lognormal', 0.5),  'rf_BC': ('lognormal', 0.5), 'rf_NO3': ('lognormal', 0.5), 'rf_SOA': ('lognormal', 0.5)}


## create parameters
Par0 = load_all_param('Houghton_2017')
Par_mc = generate_config(Par0, 10000)
Par_ar6 = adjust_config(Par_mc, shift=shift, noise=noise)

