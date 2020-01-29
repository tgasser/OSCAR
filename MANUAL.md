# OSCAR v3.0
_Manual_ 


## Basic simulation

Essentially, to make a basic simulation one must: 
1. import the `OSCAR` object from `core_fct.fct_process`; 
2. define the `Ini` (initial state), `For` (forcing data) and `Par` (parameters) arguments; 
3. call `OSCAR` with these arguments (and possibly other optional arguments). 

The `run_model` function found in `core_fct.fct_main` is a wrapper that does all that and more, using internal data for forcings and parameters. However, this function may not be adequate for advanced usage of OSCAR, in which case it should be used as inspiration for defining one's own run scripts. To further help with this, the `run_scripts` folder contains a few basic examples.


## Core structure

Here is a quick overview of the files contained in the `core_fct` folder and their content.

| File | Content |
| --- | --- |
| `cls_main` | definition of the `Model` and `Process` classes upon which OSCAR v3 is based |
| `fct_ancillary` | a bunch of useful functions, notably including the solving schemes, a generic loading function called `load_data`, and a function to regionally aggregate datasets called `aggreg_region` |
| `fct_genD` | functions to generate consistent timeseries of drivers |
| `fct_genMC` | functions to generate the Monte Carlo setup |
| `fct_loadD` | functions to load the primary drivers |
| `fct_loadP` | functions to load the primary parameters, some of them being loaded from files and others manually written there |
| `fct_main` | wrapper function to run the model in a not-so-flexible standard mode |
| `fct_process` | equations for the physical processes constituting OSCAR; also contains `OSCAR` and submodels |


## Dimensions, drivers, variables and parameters

### Dimensions

Here is a table summerizing the various dimensions over which OSCAR's input, internal and output data may be defined. Additional dimensions can be added freely to the `Ini`, `For` and/or `Par` arguments, in which case they will be conserved throughout the run, which allows easily parallelizing experiments (e.g. scenarios). This can be heavy on the memory, however.

| Dims | Description |
| --- | --- |
| `year` | time axis |
| `config` | Monte Carlo elements |
| `spc_halo` | species of halogenated compounds |
| `box_osurf` | pools for the surface ocean carbon cycling |
| `reg_land` | land carbon-cycle regions |
| `bio_land` | land carbon-cycle biomes |
| `bio_from` | origine biomes of the land-use perturbations |
| `bio_to` | destination biomes of the land-use perturbations |
| `box_hwp` | pools of harvested wood products  |
| `reg_pf` | regions specific to the permafrost module |
| `box_thaw` | pools of thawed permafrost |
| `spc_bb` | species from biomass burning |
| `reg_slcf` | regions specific to SLCF regional saturation effects |
| `reg_snow` | regions specific to BC deposition on snow |


### Drivers

Drivers are the forcing data that need to be prescribed to the model for it to be able to run. They must be prescribed using the `For` argument when calling a `Model` object. The model automatically connects the various processes it is made of, and deduces what input data are required, so that it will display an error message if some drivers are missing in `For`. Assuming `OSCAR` has been imported, a list of the model's drivers can be displayed with `OSCAR.var_in`. More information on the drivers is available in `core_fct.fct_loadD`.

| In code | In papers | Units | Dims |
| --- | --- | --- | --- |
| `Eff` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{FF}" /> | PgC yr<sup>-1</sup> | `year, reg_land` |
| `E_CH4` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{CH_4}" /> | TgC yr<sup>-1</sup> | `year, reg_land` |
| `E_N2O` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{N_2O}" /> | TgN yr<sup>-1</sup> | `year, reg_land` |
| `E_Xhalo` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{X}" /> | Gg yr<sup>-1</sup> | `year, reg_land, spc_halo` |
| `E_NOX` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{NO_x}" /> | TgN yr<sup>-1</sup> | `year, reg_land` |
| `E_CO` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{CO}" /> | TgC yr<sup>-1</sup> | `year, reg_land` |
| `E_VOC` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{VOC}" /> | Tg yr<sup>-1</sup> | `year, reg_land` |
| `E_SO2` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{SO_2}" /> | TgS yr<sup>-1</sup> | `year, reg_land` |
| `E_NH3` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{NH_3}" /> | TgN yr<sup>-1</sup> | `year, reg_land` |
| `E_OC` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{OC}" /> | TgC yr<sup>-1</sup> | `year, reg_land` |
| `E_BC` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{BC}" /> | TgC yr<sup>-1</sup> | `year, reg_land` |
||||
| `d_Acover` | <img src="https://latex.codecogs.com/gif.latex?\delta\!\,A" /> | Mha yr<sup>-1</sup> | `year, reg_land, bio_from, bio_to` |
| `d_Hwood` | <img src="https://latex.codecogs.com/gif.latex?\delta\!\,H" /> | PgC yr<sup>-1</sup> | `year, reg_land, bio_land` |
| `d_Ashift` | <img src="https://latex.codecogs.com/gif.latex?\delta\!\,S" /> | Mha yr<sup>-1</sup> | `year, reg_land, bio_from, bio_to` |
||||
| `RF_contr` | <img src="https://latex.codecogs.com/gif.latex?\mathrm{RF}_\mathrm{con}" /> | W m<sup>-2</sup> | `year` |
| `RF_volc` | <img src="https://latex.codecogs.com/gif.latex?\mathrm{RF}_\mathrm{volc}" /> | W m<sup>-2</sup> | `year` |
| `RF_solar` | <img src="https://latex.codecogs.com/gif.latex?\mathrm{RF}_\mathrm{solar}" /> | W m<sup>-2</sup> | `year` |

### Variables

Each of the model's variable is defined through a `Process` object; and a `Model` object is essentially a collection of connected processes. Prognostic variables (i.e. state variables) are those defined through a time-differential equation, while diagnostic variables are defined at any time _t_ as a function of prognostic variables and/or other diagnostic variables at that same time _t_. When solving, at every single timestep, the model first solves all prognostic variables, and only then calculates the diagnostic variables. Assuming `OSCAR` has been imported, a list of the model's variables can be displayed with `OSCAR.proc_all`, or somewhat equivalently with `OSCAR.var_mid | OSCAR.var_out`. Prognostic and diagnostic variables can be displayed with `OSCAR.var_prog` and `OSCAR.var_diag`, respectively. More information on each variable/process is available in `core_fct.fct_process`.

| In code | In papers | Units | Dims | Prog? |
| --- | --- | --- | --- | --- |
| `D_pCO2` | <img src="https://latex.codecogs.com/gif.latex?\mathcal{F}_\mathrm{pCO_2}" /> | ppm | - ||
 `D_mld` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,h_\mathrm{mld}" /> | m | - ||
| `D_dic` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{dic}" /> | µmol kg<sup>-1</sup> | - ||
| `D_Fin` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{in}" /> | PgC yr<sup>-1</sup> | `box_osurf` ||
| `D_Fout` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{out}" /> | PgC yr<sup>-1</sup> | `box_osurf` ||
| `D_Fcirc` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{circ}" /> | PgC yr<sup>-1</sup> | `box_osurf` ||
| `D_Focean` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,ocean}" /> | PgC yr<sup>-1</sup> | - ||
| `D_Cosurf` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{surf}" /> | PgC | - | **yes** |
||||||
| `f_fert` | <img src="https://latex.codecogs.com/gif.latex?\mathcal{F}_\mathrm{fert}" /> | 1 | `reg_land, bio_land` ||
| `D_npp` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{npp}" />| PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `f_igni` | <img src="https://latex.codecogs.com/gif.latex?\mathcal{F}_\mathrm{igni}" /> | 1 | `reg_land, bio_land` ||
| `D_efire` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,e_\mathrm{fire}" /> | PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_fmort` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,f_\mathrm{mort}" /> | PgC Mha<sup>-1</sup> yr<sup>-1</sup> |`reg_land, bio_land` ||
| `f_resp` | <img src="https://latex.codecogs.com/gif.latex?\mathcal{F}_\mathrm{resp}" /> | 1 | `reg_land, bio_land` ||
| `D_rh1` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{rh}_\mathrm{litt}" /> | PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_fmet` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,f_\mathrm{met}" /> |  PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_rh2` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{rh}_\mathrm{soil}" /> | PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_nbp` | - |  PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_cveg` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,c_\mathrm{veg}" /> | PgC Mha<sup>-1</sup> | `reg_land, bio_land` | **yes** |
| `D_csoil1` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,c_\mathrm{litt}" /> | PgC Mha<sup>-1</sup> | `reg_land, bio_land` | **yes** |
| `D_csoil2` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,c_\mathrm{soil}" /> | PgC Mha<sup>-1</sup> | `reg_land, bio_land` | **yes** |
||||||
| `D_Fveg_bk` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to` ||
| `D_Fsoil1_bk` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to` ||
| `D_Fsoil2_bk` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to` ||
| `D_Fslash` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to` ||
| `D_Fhwp` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to, box_hwp` ||
| `D_NPP_bk` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to` ||
| `D_Efire_bk` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to` ||
| `D_Fmort_bk` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to` ||
| `D_RH1_bk` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to` ||
| `D_Fmet_bk` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to` ||
| `D_RH2_bk` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to` ||
| `D_Ehwp` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to` ||
| `D_NBP_bk` | - | PgC yr<sup>-1</sup> | `reg_land, bio_from, bio_to, box_hwp` ||
| `D_Eluc` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{LUC}" /> | PgC yr<sup>-1</sup> | - ||
| `D_Fland` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,land}" /> | PgC yr<sup>-1</sup> | - ||
| `D_Aland` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,A" />  | Mha | `reg_land, bio_land` | **yes** |
| `D_Cveg_bk` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{veg,luc}" /> | PgC | `reg_land, bio_from, bio_to` | **yes** |
| `D_Csoil1_bk` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{litt,luc}" /> | PgC | `reg_land, bio_from, bio_to` | **yes** |
| `D_Csoil_bk` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{soil,luc}" /> | PgC | `reg_land, bio_from, bio_to` | **yes** |
| `D_Chwp` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{hwp,luc}" /> | PgC | `reg_land, bio_from, bio_to, box_hwp` | **yes** |
||||||
| `f_resp_pf` | - | 1 | `reg_pf` ||
| `D_pthaw_bar` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\bar{p}_\mathrm{thaw}" /> | 1 | `reg_pf` ||
| `d_pthaw` | <img src="https://latex.codecogs.com/gif.latex?\textstyle{\frac{\mathrm{d}}{\mathrm{d}t}}p_\mathrm{thaw}" /> | yr<sup>-1</sup> | `reg_pf` ||
| `D_pthaw` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,p_\mathrm{thaw}" /> | 1 | `reg_pf` | **yes** |
| `D_Fthaw` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{thaw}" /> | PgC yr<sup>-1</sup> | `reg_pf` ||
| `D_Ethaw` | - | PgC yr<sup>-1</sup> | `reg_pf, box_thaw` ||
| `D_Epf` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{pf}" /> | PgC yr<sup>-1</sup> | `reg_pf` ||
| `D_Epf_CO2` | - | PgC yr<sup>-1</sup> | `reg_pf` ||
| `D_Epf_CH4` | - | TgC yr<sup>-1</sup> | `reg_pf` ||
| `D_Cfroz` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{froz}" /> | PgC | `reg_pf` | **yes** |
| `D_Cthaw` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{thaw}" /> | PgC | `reg_pf, box_thaw` | **yes** |
||||||
| `D_CO2` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{CO_2}" /> | ppm | - | **yes** |
| `AF` | - | 1 | - ||
| `kS` | - | yr<sup>-1</sup> | - ||
| `RF_CO2` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{CO_2}" /> | W m<sup>-2</sup> | - ||
||||||
| `D_Efire` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{fire}" /> | PgC yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_Ebb_nat` | - | Tg<i>X</i> yr<sup>-1</sup> | `reg_land, bio_land, spc_bb` ||
| `D_Ebb_ant` | - | Tg<i>X</i> yr<sup>-1</sup> | `reg_land, bio_land, spc_bb` ||
| `D_Ebb` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{bb }" /> | Tg<i>X</i> yr<sup>-1</sup> | `reg_land, bio_land, spc_bb` ||
||||||
| `D_CH4_lag` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{CH_4}_\mathrm{lag}" /> | ppb | - | **yes** |
| `D_N2O_lag` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{N_2O}_\mathrm{lag}" /> | ppb | - | **yes** |
| `D_Xhalo_lag` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,X_\mathrm{lag}" /> | ppt | `spc_halo` | **yes** |
||||||
| `D_Ta` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,T_A" /> | K | - ||
| `D_f_Qa` | <img src="https://latex.codecogs.com/gif.latex?\frac{\Delta\!\,Q_A}{Q_{A,0}}" /> | 1 | - ||
| `f_kOH` | <img src="https://latex.codecogs.com/gif.latex?\mathcal{F}_\mathrm{OH}" /> | 1 | - ||
| `D_Foh_CH4` | - | TgC yr<sup>-1</sup> | - ||
| `D_Fhv_CH4` | -  | TgC yr<sup>-1</sup> | - ||
| `D_Fsoil_CH4` | - | TgC yr<sup>-1</sup> | - ||
| `D_Focean_CH4` | - | TgC yr<sup>-1</sup> | - ||
| `D_Fsink_CH4` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow}^\mathrm{CH_4}" /> | TgC yr<sup>-1</sup> | - ||
| `D_Foxi_CH4` | - | PgC yr<sup>-1</sup> | - ||
||||||
| `D_ewet` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,e_\mathrm{wet}" /> | TgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land` ||
| `D_Awet` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,A_\mathrm{wet}" /> | Mha | `reg_land` ||
| `D_Ewet` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{wet}" /> | TgC yr<sup>-1</sup> | `reg_land` ||
||||||
| `D_CH4` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{CH_4}" /> | ppb | - | **yes** |
| `tau_CH4` | - | yr | - ||
| `RF_CH4` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{CH_4}" /> | W m<sup>-2</sup> | - ||
| `RF_H2Os` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{H_2Os}" /> | W m<sup>-2</sup> | - ||
||||||
| `D_f_ageair` | - | 1 | - ||
| `f_hv` | <img src="https://latex.codecogs.com/gif.latex?\mathcal{F}_\mathrm{h\nu}" /> | 1 | - ||
| `D_Fhv_N2O` | -  | TgN yr<sup>-1</sup> | - ||
| `D_Fsink_N2O` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow}^\mathrm{N_2O}" /> | TgN yr<sup>-1</sup> | - ||
||||||
| `D_N2O` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{N_2O}" /> | ppb | - | **yes** |
| `tau_N2O` | - | yr | - ||
| `RF_N2O` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{N_2O}" /> | W m<sup>-2</sup> | - ||
||||||
| `D_Foh_Xhalo` | - | Gg yr<sup>-1</sup> | `spc_halo` ||
| `D_Fhv_CH4` | -  | Gg yr<sup>-1</sup> | `spc_halo` ||
| `D_Fother_CH4` | - | Gg yr<sup>-1</sup> | `spc_halo` ||
| `D_Fsink_CH4` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow}^X" /> | Gg yr<sup>-1</sup> | `spc_halo` ||
||||||
| `D_Xhalo` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,X" /> | ppt | `spc_halo` | **yes** |
| `RF_Xhalo` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^X" /> | W m<sup>-2</sup> | `spc_halo` ||
| `RF_halo` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{halo}" /> | W m<sup>-2</sup> | - ||
||||||
| `D_O3t` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{O_3t}" /> | DU | - ||
| `RF_O3t` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{O_3t}" /> | W m<sup>-2</sup> | - ||
||||||
| `D_EESC` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{EESC}" /> | ppt | - ||
| `D_O3s` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{O_3s}" /> | DU | - ||
| `RF_O3s` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{O_3s}" /> | W m<sup>-2</sup> | - ||
||||||
| `D_Edms` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{DMS}" /> | TgS yr<sup>-1</sup> | - ||
| `D_Ebvoc` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{BVOC}" /> | Tg yr<sup>-1</sup> | - ||
| `D_Edust` | - | Tg yr<sup>-1</sup> | - ||
| `D_Esalt` | - | Tg yr<sup>-1</sup> | - ||
||||||
| `D_SO4` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{SO_4}" /> | Tg | - ||
| `D_POA` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{POA}" /> | Tg | - ||
| `D_BC` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{BC}" /> | Tg | - ||
| `D_NO3` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{NO_3}" /> | Tg | - ||
| `D_SOA` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{SOA}" /> | Tg | - ||
| `D_Mdust` | - | Tg | - ||
| `D_Msalt` | - | Tg | - ||
| `RF_SO4` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{SO_4}" /> | W m<sup>-2</sup> | - ||
| `RF_POA` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{POA}" /> | W m<sup>-2</sup> | - ||
| `RF_BC` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{BC}" /> | W m<sup>-2</sup> | - ||
| `RF_NO3` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{NO_3}" /> | W m<sup>-2</sup> | - ||
| `RF_SOA` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{SOA}" /> | W m<sup>-2</sup> | - ||
| `RF_dust` | - | W m<sup>-2</sup> | - ||
| `RF_salt` | - | W m<sup>-2</sup> | - ||
||||||
| `D_AERsol` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{AER}_\mathrm{sol}" /> | Tg | - ||
| `RF_cloud1` | - | W m<sup>-2</sup> | - ||
| `RF_cloud2` | - | W m<sup>-2</sup> | - ||
| `RF_cloud` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{cloud}" /> | W m<sup>-2</sup> | - ||
||||||
| `RF_snow` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{BCsnow}" /> | W m<sup>-2</sup> | - ||
||||||
| `RF_lcc` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{LCC}" /> | W m<sup>-2</sup> | - ||
||||||
| `RF_nonCO2` | - | W m<sup>-2</sup> | - ||
| `RF_wmghg` | - | W m<sup>-2</sup> | - ||
| `RF_strat` | - | W m<sup>-2</sup> | - ||
| `RF_scatter` | - | W m<sup>-2</sup> | - ||
| `RF_absorb` | - | W m<sup>-2</sup> | - ||
| `RF_AERtot` | - | W m<sup>-2</sup> | - ||
| `RF_slcf` | - | W m<sup>-2</sup> | - ||
| `RF_alb` | - | W m<sup>-2</sup> | - ||
| `RF` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}" /> | W m<sup>-2</sup> | - ||
| `RF_warm` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}_\mathrm{warm}" /> | W m<sup>-2</sup> | - ||
| `RF_atm` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}_\mathrm{atm}" /> | W m<sup>-2</sup> | - ||
||||||
| `D_Tg` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,T_G" /> | K | - | **yes** |
| `D_Td` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,T_D" /> | K | - | **yes** |
| `d_Tg` | <img src="https://latex.codecogs.com/gif.latex?\textstyle{\frac{\mathrm{d}}{\mathrm{d}t}}T_G" /> | K yr<sup>-1</sup> | - |
| `D_Tl` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,T_L" /> | K | `reg_land` ||
| `D_To` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,T_S" /> | K | - ||
||||||
| `D_Pg` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,P_G" /> | mm yr<sup>-1</sup> | - | **yes** |
| `D_Pl` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,P_L" /> | mm yr<sup>-1</sup> | `reg_land` ||
||||||
| `D_OHC` | <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{OHC}" /> | ZJ | - | **yes** |
||||||
| `D_pH` | - | 1 | - ||


### Parameters

Parameters are implicitly defined when creating a model's processes. When OSCAR is run, it does not check whether the needed parameters are actually provided in the `Par` argument. Primary parameters can be loaded with the `load_all_param` function defined in `core_fct.fct_loadP`. Many parameters have several possible values, and these different configurations are defined along the various `mod_` dimensions of the dataset containing the primary parameters. Sets of randomly drawn parameters for Monte Carlo runs can be generated using the `generate_config` function defined in `core_fct.fct_genMC`. More information on each parameter is available in `core_fct.fct_loadP`.

| In code | In papers | Units | Dims | Mods
| --- | --- | --- | --- | --- |
| `a_dic` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{sol}" /> | µmol kg<sup>-1</sup> [ppm m<sup>-3</sup>]<sup>-1</sup> | - | - |
| `mld_0` | <img src="https://latex.codecogs.com/gif.latex?h_{\mathrm{mld},0}" /> | m | - | `mod_Focean_struct` |
| `A_ocean` | <img src="https://latex.codecogs.com/gif.latex?A_\mathrm{ocean}" /> | m<sup>2</sup> | - | `mod_Focean_struct` |
| `To_0` | <img src="https://latex.codecogs.com/gif.latex?T_{S,0}" /> | K | - | `mod_Focean_struct` |
| `v_fg` | <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{fg}" /> | yr<sup>-1</sup> | - | `mod_Focean_struct` |
| `p_circ` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{circ}" /> | 1 | `box_osurf` | `mod_Focean_struct` |
| `t_circ` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{circ}" /> | yr | `box_osurf` | `mod_Focean_struct` |
||||||
| `pCO2_is_Pade` | - | `bool` | - | `mod_Focean_chem` |
||||||
| `p_mld` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{mld}" /> | 1 | - | `mod_Focean_trans` |
| `g_mld` | <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{mld}" /> | K<sup>-1</sup> | - | `mod_Focean_trans` |
||||||
| `fert_is_Log` | - | `bool` | - | `mod_Fland_fert` |
| `k_met` | <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{met}" /> | 1 | - | - |
| `t_shift` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{shift}" /> | yr | - | - |
||||||
| `npp_0` | <img src="https://latex.codecogs.com/gif.latex?\eta" /> | PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind` |
| `igni_0` | <img src="https://latex.codecogs.com/gif.latex?\iota" /> | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Efire_preind` |
| `cveg_0` | <img src="https://latex.codecogs.com/gif.latex?c_{\mathrm{veg},0}" /> | PgC Mha<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind, mod_Efire_preind` |
| `mu_0` | <img src="https://latex.codecogs.com/gif.latex?\mu" /> | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind` |
| `rho1_0` | <img src="https://latex.codecogs.com/gif.latex?\rho_\mathrm{litt}" /> | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind` |
| `csoil1_0` | <img src="https://latex.codecogs.com/gif.latex?c_{\mathrm{litt},0}" /> | PgC Mha<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind, mod_Efire_preind` |
| `rho2_0` | <img src="https://latex.codecogs.com/gif.latex?\rho_\mathrm{soil}" /> | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind` |
| `csoil2_0` | <img src="https://latex.codecogs.com/gif.latex?c_{\mathrm{soil},0}" /> | PgC Mha<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind, mod_Efire_preind` |
| `p_agb` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{agb}" /> | 1 | `reg_land, bio_land` | `mod_Eluc_agb` |
||||||
| `b_npp` | <img src="https://latex.codecogs.com/gif.latex?\beta_\mathrm{npp}" /> | 1 | `reg_land, bio_land` | `mod_Fland_trans` |
| `b2_npp` | <img src="https://latex.codecogs.com/gif.latex?\tilde{\beta}_\mathrm{npp}" /> | 1 | `reg_land, bio_land` | `mod_Fland_trans` |
| `CO2_cp` | <img src="https://latex.codecogs.com/gif.latex?\mathrm{CO_2}_\mathrm{cp}" /> | ppm | `reg_land, bio_land` | `mod_Fland_trans` |
| `g_nppT` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{npp},T}" /> | K<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_trans, mod_Fland_fert` |
| `g_nppP` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{npp},P}" /> | [mm yr<sup>-1</sup>]<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_trans, mod_Fland_fert` |
| `g_rhoT` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{resp},T}" /> | K<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_trans, mod_Fland_resp` |
| `g_rhoT1` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{resp},T_1}" /> | K<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_trans, mod_Fland_resp` |
| `g_rhoT2` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{resp},T_2}" /> | K<sup>-2</sup> | `reg_land, bio_land` | `mod_Fland_trans, mod_Fland_resp` |
| `g_rhoP` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{resp},P}" /> | [mm yr<sup>-1</sup>]<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_trans` |
| `g_igniC` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{igni},C}" /> | ppm<sup>-1</sup> | `reg_land, bio_land` | `mod_Efire_trans` |
| `g_igniT` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{igni},T}" /> | K<sup>-1</sup> | `reg_land, bio_land` | `mod_Efire_trans` |
| `g_igniP` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{igni},P}" /> | [mm yr<sup>-1</sup>]<sup>-1</sup> | `reg_land, bio_land` | `mod_Efire_trans` |
||||||
| `t_hwp` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{hwp}" /> | yr | `box_hwp` | `mod_Ehwp_tau` |
| `w_t_hwp` | - | 1 | - | `mod_Ehwp_speed` |
| `p_hwp_bb` | - | 1 | `box_hwp` | - |
| `p_hwp` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{hwp}" /> | 1 | `reg_land, bio_land, box_hwp` | `mod_Ehwp_bb` |
||||||
| `a_bb` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{bb}" /> | Tg<i>X</i> PgC<sup>-1</sup> | `reg_land, bio_land, spc_bb` | - |
||||||
| `Cfroz_0` | <img src="https://latex.codecogs.com/gif.latex?C_{\mathrm{froz},0}}" /> | PgC | `reg_pf` | `mod_Epf_main` |
| `w_clim_pf` | <img src="https://latex.codecogs.com/gif.latex?\omega_{T_\mathrm{pf}}}" /> | 1 | `reg_pf` | `mod_Epf_main` |
| `g_respT_pf` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{pf},T_1}}" /> | K<sup>-1</sup> | `reg_pf` | `mod_Epf_main` |
| `g_respT2_pf` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{pf},T_2}}" /> | K<sup>-2</sup> | `reg_pf` | `mod_Epf_main` |
| `k_resp_pf` | <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{resp,pf}}" /> | 1 | `reg_pf` | `mod_Epf_main` |
| `pthaw_min` | <img src="https://latex.codecogs.com/gif.latex?p_\mathrm{thaw,min}}" /> | 1 | `reg_pf` | `mod_Epf_main` |
| `g_pthaw` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{p_\mathrm{thaw}}}" /> | K<sup>-1</sup> | `reg_pf` | `mod_Epf_main` |
| `k_pthaw` | <img src="https://latex.codecogs.com/gif.latex?\kappa_{p_\mathrm{thaw}}}" /> | 1 | `reg_pf` | `mod_Epf_main` |
| `v_thaw` | <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{thaw}}" /> | yr<sup>-1</sup> | `reg_pf` | `mod_Epf_main` |
| `v_froz` | <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{froz}}" /> | yr<sup>-1</sup> | `reg_pf` | `mod_Epf_main` |
| `p_pf_thaw` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{thaw}}" /> | 1 | `reg_pf, box_thaw` | `mod_Epf_main` |
| `t_pf_thaw` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{thaw}}" /> | yr | `reg_pf, box_thaw` | `mod_Epf_main` |
| `p_pf_inst` | - | 1 | - | - |
| `p_pf_CH4` | - | 1 | - | `mod_Epf_CH4` |
||||||
| `ewet_0` | <img src="https://latex.codecogs.com/gif.latex?e_{\mathrm{wet},0}" /> | TgC yr<sup>-1</sup> | `reg_land` | `mod_Ewet_preind` |
| `Awet_0` | <img src="https://latex.codecogs.com/gif.latex?A_{\mathrm{wet},0}" /> | Mha | `reg_land` | `mod_Ewet_preind` |
| `p_wet` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{wet}" /> | 1 | `reg_land, bio_land` | `mod_Ewet_preind` |
| `g_wetC` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{wet},C}" /> | ppm<sup>-1</sup> | `reg_land` | `mod_Awet_trans` |
| `g_wetT` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{wet},T}" /> | K<sup>-1</sup> | `reg_land` | `mod_Awet_trans` |
| `g_wetP` | <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{wet},P}" /> | [mm yr<sup>-1</sup>]<sup>-1</sup> | `reg_land` | `mod_Awet_trans` |
||||||
| `a_CO2` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{atm}^\mathrm{CO_2}" /> | PgC ppm<sup>-1</sup> | - | - |
| `a_CH4` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{atm}^\mathrm{CH_4}" /> | TgC ppb<sup>-1</sup> | - | - |
| `a_N2O` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{atm}^\mathrm{N_2O}" /> | TgC ppb<sup>-1</sup> | - | - |
| `a_Xhalo` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{atm}^X" /> | Gg ppt<sup>-1</sup> | `spc_halo` | - |
| `a_SO4` | - | Tg TgS<sup>-1</sup> | - | - |
| `a_POM` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{OM}^\mathrm{OC}" /> | Tg TgC<sup>-1</sup> | - | `mod_POA_conv` |
| `a_NO3` | - | Tg TgN<sup>-1</sup> | - | - |
| `CO2_0` | <img src="https://latex.codecogs.com/gif.latex?\mathrm{CO_2}_0" /> | ppm | - | - |
| `CH4_0` | <img src="https://latex.codecogs.com/gif.latex?\mathrm{CH_4}_0" /> | ppb | - | - |
| `N2O_0` | <img src="https://latex.codecogs.com/gif.latex?\mathrm{N_2O}_0" /> | ppb | - | - |
| `Xhalo_0` | <img src="https://latex.codecogs.com/gif.latex?X_0" />  | ppt | `spc_halo` | - |
| `p_CH4geo` | - | 1 | - | - |
||||||
| `g_ageair` | <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{age}" /> | K<sup>-1</sup> | - | `mod_Fhv_ageair` |
||||||
| `w_t_OH` | - | 1 | - | - |
| `w_t_hv` | - | 1 | - | - |
||||||
| `t_OH_CH4` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{OH}^\mathrm{CH_4}" /> | yr | - | `mod_Foh_tau` |
| `t_hv_CH4` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{h\nu}^\mathrm{CH_4}" /> | yr | - | - |
| `t_soil_CH4` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{soil}^\mathrm{CH_4}" /> | yr | - | - |
| `t_ocean_CH4` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{ocean}^\mathrm{CH_4}" /> | yr | - | - |
||||||
| `x_OH_Ta` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{T_A}^\mathrm{OH}" /> | 1 | - | `mod_Foh_trans` |
| `x_OH_Qa` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{Q_A}^\mathrm{OH}" /> | 1 | - | `mod_Foh_trans` |
| `x_OH_O3s` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{O_3s}^\mathrm{OH}" /> | 1 | - | `mod_Foh_trans` |
| `x_OH_CH4` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{CH4}^\mathrm{OH}" /> | 1 | - | `mod_Foh_trans` |
| `x_OH_NOX` | <img src="https://latex.codecogs.com/gif.latex?\tilde{\chi}_\mathrm{NO_x}^\mathrm{OH}" /> | 1 | - | `mod_Foh_trans` |
| `x_OH_CO` | <img src="https://latex.codecogs.com/gif.latex?\tilde{\chi}_\mathrm{CO}^\mathrm{OH}" /> | 1 | - | `mod_Foh_trans` |
| `x_OH_VOC` | <img src="https://latex.codecogs.com/gif.latex?\tilde{\chi}_\mathrm{VOC}^\mathrm{OH}" /> | 1 | - | `mod_Foh_trans` |
| `x2_OH_NOX` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{NO_x}^\mathrm{OH}" /> | [TgN yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_Foh_trans` |
| `x2_OH_CO` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{CO}^\mathrm{OH}" /> | [TgC yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_Foh_trans` |
| `x2_OH_VOC` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{VOC}^\mathrm{OH}" /> | [Tg yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_Foh_trans` |
| `w_clim_Ta` | <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{T_A}" /> | 1 | - | - |
| `k_Qa` | <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{Q_A}" /> | 1 | - | - |
| `Ta_0` | <img src="https://latex.codecogs.com/gif.latex?T_{A,0}" /> | K | - | - |
| `k_svp` | <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{svp}" /> | 1 | - | - |
| `T_svp` | <img src="https://latex.codecogs.com/gif.latex?T_\mathrm{svp}" /> | K | - | - |
| `O3s_0` | <img src="https://latex.codecogs.com/gif.latex?\mathrm{O_3s}_0" /> | DU | - | - |
| `Enat_NOX` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{nat}^\mathrm{NO_x}" /> | TgN yr<sup>-1</sup> | - | - |
| `Enat_CO` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{nat}^\mathrm{CO}" /> | TgC yr<sup>-1</sup> | - | - |
| `Enat_VOC` | <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{nat}^\mathrm{VOC}" /> | Tg yr<sup>-1</sup> | - | - |
| `kOH_is_Log` | - | `bool` | - | `mod_Foh_fct` |
||||||
| `t_hv_N2O` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{h\nu}^\mathrm{N_2O}" /> | yr | - | `mod_Fhv_tau` |
||||||
| `x_hv_N2O` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{N2O}^\mathrm{h\nu}" /> | 1 | - | `mod_Fhv_trans` |
| `x_hv_EESC` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{EESC}^\mathrm{h\nu}" /> | 1 | - | `mod_Fhv_trans` |
| `x_hv_age` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{age}^\mathrm{h\nu}" /> | 1 | - | `mod_Fhv_trans` |
||||||
| `t_OH_Xhalo` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{OH}^X}" /> | yr | `spc_halo` | - |
| `t_hv_Xhalo` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{h\nu}^X}" /> | yr | `spc_halo` | - |
| `t_other_Xhalo` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{othr}^X}" /> | yr | `spc_halo` | - |
||||||
| `p_reg_slcf` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{reg}" /> | 1 | `reg_land, reg_slcf` | - |
||||||
| `w_reg_NOX` | <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{NO_x}" /> | 1 | `reg_slcf` | `mod_O3t_regsat` |
| `w_reg_CO` | <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{CO}" /> | 1 | `reg_slcf` | `mod_O3t_regsat` |
| `w_reg_VOC` | <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{VOC}" /> | 1 | `reg_slcf` | `mod_O3t_regsat` |
||||||
| `x_O3t_CH4` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{CH_4}^\mathrm{O_3t}" /> | DU | - | `mod_O3t_emis` |
| `x_O3t_NOX` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{NO_x}^\mathrm{O_3t}" /> | DU [TgN yr<sup>-1</sup>]<sup>-1</sup>  | - | `mod_O3t_emis` |
| `x_O3t_CO` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{CO}^\mathrm{O_3t}" /> | DU [TgC yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_O3t_emis` |
| `x_O3t_VOC` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{VOC}^\mathrm{O_3t}" /> | DU [Tg yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_O3t_emis` |
| `G_O3t` | <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{O_3t}" /> | DU K<sup>-1</sup> | - | `mod_O3t_clim` |
||||||
| `t_lag` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{lag}" /> | yr | - | - |
| `p_fracrel` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{rel}^X" /> | 1 | `spc_halo` | `mod_O3s_fracrel` |
| `k_Br_Cl` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{Cl}^\mathrm{Br}" /> | 1 | - | - |
| `n_Cl` | <img src="https://latex.codecogs.com/gif.latex?n_\mathrm{Cl}^X}" /> | 1 | `spc_halo` | - |
| `n_Br` | <img src="https://latex.codecogs.com/gif.latex?n_\mathrm{Br}^X}" /> | 1 | `spc_halo` | - |
| `EESC_x` | <img src="https://latex.codecogs.com/gif.latex?\mathrm{EESC}_\times" /> | ppt | - | - |
| `k_EESC_N2O` | <img src="https://latex.codecogs.com/gif.latex?\frac{\chi_\mathrm{N2O}^\mathrm{O_3s}}{\chi_\mathrm{EESC}^\mathrm{O_3s}}" /> | ppt ppb<sup>-1</sup> | - | `mod_O3s_nitrous, mod_O3s_fracrel` |
||||||
| `x_O3s_EESC` | <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{EESC}^\mathrm{O_3s}" /> | DU ppt<sup>-1</sup> | - | `mod_O3s_trans` |
| `G_O3s` | <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{O_3s}" /> | DU K<sup>-1</sup> | - | `mod_O3s_trans` |
||||||
| `w_reg_SO2` | <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{SO_2}" /> | 1 | `reg_slcf` | `mod_SO4_regsat` |
| `w_reg_OC` | <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{OC}" /> | 1 | `reg_slcf` | `mod_POA_regsat` |
| `w_reg_BC` | <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{BC}" /> | 1 | `reg_slcf` | `mod_BC_regsat` |
||||||
| `t_SO2` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{SO_2}" /> | yr | - | `mod_SO4_load` |
| `t_DMS` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{DMS}" /> | yr | - | `mod_SO4_load` |
| `G_SO4` | <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{SO_4}" /> | Tg K<sup>-1</sup> | - | `mod_SO4_load` |
| `t_OMff` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{OM,ff}" /> | yr | - | `mod_POA_load` |
| `t_OMbb` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{OM,bb}" /> | yr | - | `mod_POA_load` |
| `G_POA` | <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{POA}" /> | Tg K<sup>-1</sup> | - | `mod_POA_load` |
| `t_BCff` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{BC,ff}" /> | yr | - | `mod_BC_load` |
| `t_BCbb` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{BC,bb}" /> | yr | - | `mod_BC_load` |
| `G_BC` | <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{BC}" /> | Tg K<sup>-1</sup> | - | `mod_BC_load` |
| `t_NOX` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{NO_x}" /> | yr | - | `mod_NO3_load` |
| `t_NH3` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{NH_3}" /> | yr | - | `mod_NO3_load` |
| `G_NO3` | <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{NO_3}" /> | Tg K<sup>-1</sup> | - | `mod_NO3_load` |
| `t_VOC` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{VOC}" /> | yr | - | `mod_SOA_load` |
| `t_BVOC` | <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{BVOC}" /> | yr | - | `mod_SOA_load` |
| `G_SOA` | <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{SOA}" /> | Tg K<sup>-1</sup> | - | `mod_SOA_load` |
| `t_dust` | - | yr | - | `mod_Mdust_load` |
| `G_dust` | - | Tg K<sup>-1</sup> | - | `mod_Mdust_load` |
| `t_salt` | - | yr | - | `mod_Msalt_load` |
| `G_salt` | - | Tg K<sup>-1</sup> | - | `mod_Msalt_load` |
||||||
| `p_sol_SO4` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{sol}^\mathrm{SO_4}" /> | 1 | - | `mod_RFcloud_solub` |
| `p_sol_POA` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{sol}^\mathrm{POA}" /> | 1 | - | `mod_RFcloud_solub` |
| `p_sol_BC` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{sol}^\mathrm{BC}" /> | 1 | - | `mod_RFcloud_solub` |
| `p_sol_NO3` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{sol}^\mathrm{NO_3}" /> | 1 | - | `mod_RFcloud_solub` |
| `p_sol_SOA` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{sol}^\mathrm{SOA}" /> | 1 | - | `mod_RFcloud_solub` |
| `p_sol_dust` | - | 1 | - | `mod_RFcloud_solub` |
| `p_sol_salt` | - | 1 | - | `mod_RFcloud_solub` |
||||||
| `rf_CO2` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{CO_2}" /> | W m<sup>-2</sup> | - | - |
| `rf_CH4` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{CH_4}" /> | W m<sup>-2</sup> ppb<sup>-0.5</sup> | - | - |
| `rf_N2O` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{N_2O}" /> | W m<sup>-2</sup> ppb<sup>-0.5</sup> | - | - |
| `k_rf_H2Os` | <img src="https://latex.codecogs.com/gif.latex?\frac{\alpha_\mathrm{rf}^\mathrm{H_2Os}}{\alpha_\mathrm{rf}^\mathrm{CH_4}}" /> | 1 | - | - |
| `rf_Xhalo` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^X" /> | W m<sup>-2</sup> ppt<sup>-1</sup> | `spc_halo` | - |
||||||
| `rf_O3t` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{O_3t}" /> | W m<sup>-2</sup> DU<sup>-1</sup> | - | `mod_O3t_radeff` |
| `rf_O3s` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{O_3s}" /> | W m<sup>-2</sup> DU<sup>-1</sup> | - | `mod_O3s_radeff` |
||||||
| `rf_SO4` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{SO_4}" /> | W m<sup>-2</sup> Tg<sup>-1</sup> | - | `mod_SO4_radeff` |
| `rf_POA` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{POA}" /> | W m<sup>-2</sup> Tg<sup>-1</sup> | - | `mod_POA_radeff` |
| `rf_BC` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{BC}" /> | W m<sup>-2</sup> Tg<sup>-1</sup> | - | `mod_BC_radeff` |
| `rf_NO3` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{NO_3}" /> | W m<sup>-2</sup> Tg<sup>-1</sup> | - | `mod_NO3_radeff` |
| `rf_SOA` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{SOA}" /> | W m<sup>-2</sup> Tg<sup>-1</sup> | - | `mod_SOA_radeff` |
| `rf_dust` | - | W m<sup>-2</sup> Tg<sup>-1</sup> | - | - |
| `rf_salt` | - | W m<sup>-2</sup> Tg<sup>-1</sup> | - | - |
||||||
| `k_adj_BC` | <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{adj}^\mathrm{BC}" /> | 1 | - | `mod_BC_adjust` |
||||||
| `Phi_0` | <img src="https://latex.codecogs.com/gif.latex?\Phi" /> | W m<sup>-2</sup> | - | `mod_RFcloud_erf, mod_RFcloud_solub` |
| `AERsol_0` | <img src="https://latex.codecogs.com/gif.latex?\mathrm{AER}_{\mathrm{sol},0}" /> | Tg | - | `mod_RFcloud_solub, mod_RFcloud_erf, mod_RFcloud_preind` |
||||||
| `p_reg_snow` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{reg}" /> | 1 | `reg_land, reg_snow` | - |
||||||
| `w_reg_snow` | <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{BCsnow}" /> | 1 | `reg_snow` | `mod_RFsnow_reg` |
| `rf_snow` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{BCsnow}" /> | W m<sup>-2</sup> [TgC yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_RFsnow_rf` |
||||||
| `p_trans` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{trans}" /> | 1 | - | - |
| `alpha_alb` | <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{alb}" /> | 1 | `reg_land, bio_land` | `mod_RFlcc_alb, mod_RFlcc_flux, mod_RFlcc_cover` |
| `F_rsds` | <img src="https://latex.codecogs.com/gif.latex?\phi_\mathrm{rsds}" /> | W m<sup>-2</sup> | `reg_land` | `mod_RFlcc_flux` |
||||||
| `w_warm_volc` | <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{warm}^\mathrm{volc}" /> | 1 | - | - |
| `w_warm_snow` | <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{warm}^\mathrm{BCsnow}" /> | 1 | - | `mod_RFsnow_warmeff` |
| `w_warm_lcc` | <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{warm}^\mathrm{LCC}" /> | 1 | - | `mod_RFlcc_warmeff` |
|||||||
| `p_atm_CO2` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{CO_2}" /> | 1 | - | `mod_Pg_radfact` |
| `p_atm_nonCO2` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{noCO_2}" /> | 1 | - | `mod_Pg_radfact` |
| `p_atm_O3t` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{O_3t}" /> | 1 | - | `mod_Pg_radfact` |
| `p_atm_strat` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{strat}" /> | 1 | - | `mod_Pg_radfact` |
| `p_atm_scatter` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{scatter}" /> | 1 | - | `mod_Pg_radfact` |
| `p_atm_absorb` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{absorb}" /> | 1 | - | `mod_Pg_radfact` |
| `p_atm_cloud` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{cloud}" /> | 1 | - | `mod_Pg_radfact` |
| `p_atm_alb` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{alb}" /> | 1 | - | `mod_Pg_radfact` |
| `p_atm_solar` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{solar}" /> | 1 | - | `mod_Pg_radfact` |
||||||
| `ecs_0` | <img src="https://latex.codecogs.com/gif.latex?\lambda" /> | K [W m<sup>-2</sup>]<sup>-1</sup> | - | `mod_Tg_resp` |
| `Th_g` | <img src="https://latex.codecogs.com/gif.latex?\frac{\tau_{T_G}}{\lambda}" /> | yr W m<sup>-2</sup> K<sup>-1</sup> | - | `mod_Tg_resp` |
| `Th_d` | <img src="https://latex.codecogs.com/gif.latex?\frac{\tau_{T_D}}{\lambda}" /> | yr W m<sup>-2</sup> K<sup>-1</sup> | - | `mod_Tg_resp` |
| `th_0` | <img src="https://latex.codecogs.com/gif.latex?\frac{\theta}{\lambda}" /> | W m<sup>-2</sup> K<sup>-1</sup> | - | `mod_Tg_resp` |
| `w_clim_Tl` | <img src="https://latex.codecogs.com/gif.latex?\omega_{T_L}" /> | 1 | `reg_land` | `mod_Tl_pattern, mod_Tg_resp` |
| `w_clim_To` | <img src="https://latex.codecogs.com/gif.latex?\omega_{T_S}" /> | 1 | - | `mod_Tl_pattern, mod_Tg_resp` |
||||||
| `a_prec` | <img src="https://latex.codecogs.com/gif.latex?\alpha_{P_G}" /> | mm yr<sup>-1</sup> K<sup>-1</sup> | - | `mod_Pg_resp` |
| `b_prec` | <img src="https://latex.codecogs.com/gif.latex?\beta_{P_G}" /> | mm yr<sup>-1</sup> [W m<sup>-2</sup>]<sup>-1</sup> | - | `mod_Pg_resp` |
| `w_clim_Pl` | <img src="https://latex.codecogs.com/gif.latex?\omega_{P_L}" /> | 1 | `reg_land` | `mod_Pl_pattern, mod_Pg_resp` |
||||||
| `p_ohc` | <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{ohc}" /> | 1 | - | - |
||||||
| `pH_is_Log` | - | `bool` | - | `mod_pH_fct` |
