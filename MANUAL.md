# OSCAR
_Manual_ 


## Basic simulation

Essentially, to make a basic simulation one must: 
1. import the `OSCAR` object from `core_fct.mod_process`; 
2. define the `For` (forcing data) and `Par` (parameters) arguments;
3. call `OSCAR` with these arguments (and possibly other optional arguments, like `Ini` (initial state)). 

The `run_scripts` folder contains a few extra basic examples.


## Core structure

Here is a quick overview of the files contained in the `core_fct` folder and their content.

| File | Content |
| --- | --- |
| `cls_main` | definition of the `Model` and `Process` classes upon which OSCAR v3 is based |
| `fct_calib` | functions to calibrate some of the model's parameters |
| `fct_genD` | functions to generate consistent timeseries of drivers |
| `fct_genMC` | functions to generate the Monte Carlo setup |
| `fct_loadD` | functions to load the primary drivers |
| `fct_loadP` | functions to load the primary parameters, some being loaded from files and others manually written therein |
| `fct_process_alt` | functions to replace some processes with alternative formulations |
| `fct_misc` | a bunch of useful functions, notably including the solving schemes, a generic loading function called `load_data`, and a function to regionally aggregate datasets called `aggreg_region` |
| `mod_process` | equations for the physical processes constituting OSCAR; also contains `OSCAR` and submodels |


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
| `reg_bcsnow` | regions specific to BC deposition on snow |


### Drivers

Drivers are the forcing data that need to be prescribed to the model for it to be able to run. They must be prescribed using the `For` argument when calling a `Model` object. The model automatically connects the various processes it is made of, and deduces what input data are required, so that it will display an error message if some drivers are missing in `For`. Assuming `OSCAR` has been imported, a list of the model's drivers can be displayed with `OSCAR.var_in`. More information on the drivers is available in `core_fct.fct_loadD`.

| In code | In papers | Units | Dims |
| --- | --- | --- | --- |
| `Eff` | $E_\mathrm{FF}$ | PgC yr<sup>-1</sup> | `year, reg_land` |
| `E_CH4` | $E_\mathrm{CH_4}$ | TgC yr<sup>-1</sup> | `year, reg_land` |
| `E_N2O` | $E_\mathrm{N_2O}$ | TgN yr<sup>-1</sup> | `year, reg_land` |
| `E_Xhalo` | $E_\mathrm{X}$ | Gg yr<sup>-1</sup> | `year, reg_land, spc_halo` |
| `E_NOX` | $E_\mathrm{NO_x}$ | TgN yr<sup>-1</sup> | `year, reg_land` |
| `E_CO` | $E_\mathrm{CO}$ | TgC yr<sup>-1</sup> | `year, reg_land` |
| `E_VOC` | $E_\mathrm{VOC}$ | Tg yr<sup>-1</sup> | `year, reg_land` |
| `E_SO2` | $E_\mathrm{SO_2}$ | TgS yr<sup>-1</sup> | `year, reg_land` |
| `E_NH3` | $E_\mathrm{NH_3}$ | TgN yr<sup>-1</sup> | `year, reg_land` |
| `E_OC` | $E_\mathrm{OC}$ | TgC yr<sup>-1</sup> | `year, reg_land` |
| `E_BC` | $E_\mathrm{BC}$ | TgC yr<sup>-1</sup> | `year, reg_land` |
||||
| `d_Acover` | $\delta A$ | Mha yr<sup>-1</sup> | `year, reg_land, bio_from, bio_to` |
| `d_Hwood` | $\delta H$ | PgC yr<sup>-1</sup> | `year, reg_land, bio_land` |
| `d_Ashift` | $\delta S$ | Mha yr<sup>-1</sup> | `year, reg_land, bio_from, bio_to` |
||||
| `RF_contr` | $\mathrm{RF}_\mathrm{con}$ | W m<sup>-2</sup> | `year` |
| `RF_volc` | $\mathrm{RF}_\mathrm{volc}$ | W m<sup>-2</sup> | `year` |
| `RF_solar` | $\mathrm{RF}_\mathrm{solar}$ | W m<sup>-2</sup> | `year` |

### Variables

Each of the model's variable is defined through a `Process` object; and a `Model` object is essentially a collection of connected processes. Prognostic variables (i.e. state variables) are those defined through a time-differential equation, while diagnostic variables are defined at any time _t_ as a function of prognostic variables and/or other diagnostic variables at that same time _t_. When solving, at every single timestep, the model first solves all prognostic variables, and only then calculates the diagnostic variables. Assuming `OSCAR` has been imported, a list of the model's variables can be displayed with `OSCAR.proc_all`, or somewhat equivalently with `OSCAR.var_mid | OSCAR.var_out`. Prognostic and diagnostic variables can be displayed with `OSCAR.var_prog` and `OSCAR.var_diag`, respectively. More information on each variable/process is available in `core_fct.fct_process`.

| In code | In papers | Units | Dims | Prog? |
| --- | --- | --- | --- | --- |
| `D_pCO2` | $\mathcal{F}_\mathrm{pCO_2}$ | ppm | - ||
 `D_mld` | $\Delta h_\mathrm{mld}$ | m | - ||
| `D_dic` | $\Delta \mathrm{dic}$ | µmol kg<sup>-1</sup> | - ||
| `D_Fin` | $\Delta F_\mathrm{in}$ | PgC yr<sup>-1</sup> | `box_osurf` ||
| `D_Fout` | $\Delta F_\mathrm{out}$ | PgC yr<sup>-1</sup> | `box_osurf` ||
| `D_Fcirc` | $\Delta F_\mathrm{circ}$ | PgC yr<sup>-1</sup> | `box_osurf` ||
| `D_Focean` | $\Delta F_\mathrm{\downarrow\,ocean}$ | PgC yr<sup>-1</sup> | - ||
| `D_Cosurf` | $\Delta C_\mathrm{surf}$ | PgC | - | **yes** |
||||||
| `f_fert` | $\mathcal{F}_\mathrm{fert}$ | 1 | `reg_land, bio_land` ||
| `D_npp` | $\Delta \mathrm{npp}$| PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `f_igni` | $\mathcal{F}_\mathrm{igni}$ | 1 | `reg_land, bio_land` ||
| `D_efire` | $\Delta e_\mathrm{fire}$ | PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_eharv` | $\Delta e_\mathrm{harv}$ | PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_egraz` | $\Delta e_\mathrm{graz}$ | PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_fmort1` | $\Delta f_\mathrm{mort1}$ | PgC Mha<sup>-1</sup> yr<sup>-1</sup> |`reg_land, bio_land` ||
| `D_fmort2` | $\Delta f_\mathrm{mort2}$ | PgC Mha<sup>-1</sup> yr<sup>-1</sup> |`reg_land, bio_land` ||
| `f_resp` | $\mathcal{F}_\mathrm{resp}$ | 1 | `reg_land, bio_land` ||
| `D_rh1` | $\Delta \mathrm{rh}_1$ | PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_fmet` | $\Delta f_\mathrm{met}$ |  PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_rh2` | $\Delta \mathrm{rh}_2$ | PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_nbp` | - |  PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_cveg` | $\Delta c_\mathrm{veg}$ | PgC Mha<sup>-1</sup> | `reg_land, bio_land` | **yes** |
| `D_csoil1` | $\Delta c_\mathrm{soil1}$ | PgC Mha<sup>-1</sup> | `reg_land, bio_land` | **yes** |
| `D_csoil2` | $\Delta c_\mathrm{soil2}$ | PgC Mha<sup>-1</sup> | `reg_land, bio_land` | **yes** |
||||||
| `D_Fveg_bk` | $\delta C_\mathrm{veg,bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Fsoil1_bk` | $\delta C_\mathrm{soil1,bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Fsoil2_bk` | $\delta C_\mathrm{soil2,bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Fslash1` | $\Delta F_\mathrm{slash1}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Fslash2` | $\Delta F_\mathrm{slash2}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Fhwp` | $\Delta F_\mathrm{hwp}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land, box_hwp`&#8224; ||
| `D_NPP_bk` | $\Delta \mathrm{NPP}_\mathrm{bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Efire_bk` | $\Delta E_\mathrm{fire,bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Eharv_bk` | $\Delta E_\mathrm{harv,bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Egraz_bk` | $\Delta E_\mathrm{graz,bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Fmort1_bk` | $\Delta F_\mathrm{mort1,bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Fmort2_bk` | $\Delta F_\mathrm{mort2,bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Rh1_bk` | $\Delta \mathrm{Rh}_\mathrm{1,bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Fmet_bk` | $\Delta F_\mathrm{met,bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Rh2_bk` | $\Delta \mathrm{Rh}_\mathrm{2,bk}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Ehwp` | $\Delta E_\mathrm{hwp}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land, box_hwp`&#8224; ||
| `D_NBP_bk` | - | PgC yr<sup>-1</sup> | `reg_land, bio_land`&#8224; ||
| `D_Eluc` | $\Delta E_\mathrm{LUC}$ | PgC yr<sup>-1</sup> | - ||
| `D_Fland` | $\Delta F_\mathrm{\downarrow\,land}$ | PgC yr<sup>-1</sup> | - ||
| `D_Flasc` | $\Delta F_\mathrm{LASC}$ | PgC yr<sup>-1</sup> | - ||
| `D_Aland` | $\Delta A$  | Mha | `reg_land, bio_land` | **yes** |
| `D_Cveg_bk` | $\Delta C_\mathrm{veg,bk}$ | PgC | `reg_land, bio_land`&#8224; | **yes** |
| `D_Csoil1_bk` | $\Delta C_\mathrm{soil1,bk}$ | PgC | `reg_land, bio_land`&#8224; | **yes** |
| `D_Csoil2_bk` | $\Delta C_\mathrm{soil2,bk}$ | PgC | `reg_land, bio_land`&#8224; | **yes** |
| `D_Chwp` | $\Delta C_\mathrm{hwp}$ | PgC | `reg_land, bio_land, box_hwp`&#8224; | **yes** |
||||||
| `f_resp_pf` | - | 1 | `reg_pf` ||
| `D_pthaw_bar` | $\Delta \bar{p}_\mathrm{thaw}$ | 1 | `reg_pf` ||
| `d_pthaw` | $\textstyle{\frac{\mathrm{d}}{\mathrm{d}t}}p_\mathrm{thaw}$ | yr<sup>-1</sup> | `reg_pf` ||
| `D_pthaw` | $\Delta p_\mathrm{thaw}$ | 1 | `reg_pf` | **yes** |
| `D_Fthaw` | $\Delta F_\mathrm{thaw}$ | PgC yr<sup>-1</sup> | `reg_pf` ||
| `D_Ethaw` | - | PgC yr<sup>-1</sup> | `reg_pf, box_thaw` ||
| `D_Epf` | $\Delta E_\mathrm{pf}$ | PgC yr<sup>-1</sup> | `reg_pf` ||
| `D_Epf_CO2` | - | PgC yr<sup>-1</sup> | `reg_pf` ||
| `D_Epf_CH4` | - | TgC yr<sup>-1</sup> | `reg_pf` ||
| `D_Cfroz` | $\Delta C_\mathrm{froz}$ | PgC | `reg_pf` | **yes** |
| `D_Cthaw` | $\Delta C_\mathrm{thaw}$ | PgC | `reg_pf, box_thaw` | **yes** |
||||||
| `D_CO2` | $\Delta \mathrm{CO_2}$ | ppm | - | **yes** |
| `d_CO2` | $\textstyle{\frac{\mathrm{d}}{\mathrm{d}t}}\mathrm{CO_2}$ | ppm yr<sup>-1</sup> | - ||
| `AF` | - | 1 | - ||
| `kS` | - | yr<sup>-1</sup> | - ||
| `RF_CO2` | $\Delta \mathrm{RF}^\mathrm{CO_2}$ | W m<sup>-2</sup> | - ||
||||||
| `D_Efire` | $\Delta E_\mathrm{fire}$ | PgC yr<sup>-1</sup> | `reg_land, bio_land` ||
| `D_Ebb_nat` | - | Tg<i>X</i> yr<sup>-1</sup> | `reg_land, bio_land, spc_bb` ||
| `D_Ebb_ant` | - | Tg<i>X</i> yr<sup>-1</sup> | `reg_land, bio_land, spc_bb` ||
| `D_Ebb` | $\Delta E_\mathrm{bb }$ | Tg<i>X</i> yr<sup>-1</sup> | `reg_land, bio_land, spc_bb` ||
||||||
| `D_CH4_lag` | $\Delta \mathrm{CH_4}_\mathrm{lag}$ | ppb | - | **yes** |
| `D_N2O_lag` | $\Delta \mathrm{N_2O}_\mathrm{lag}$ | ppb | - | **yes** |
| `D_Xhalo_lag` | $\Delta X_\mathrm{lag}$ | ppt | `spc_halo` | **yes** |
||||||
| `D_Ta` | $\Delta T_A$ | K | - ||
| `D_f_Qa` | $\frac{\Delta Q_A}{Q_{A,0}}$ | 1 | - ||
| `f_kOH` | $\mathcal{F}_\mathrm{OH}$ | 1 | - ||
| `D_Foh_CH4` | - | TgC yr<sup>-1</sup> | - ||
| `D_Fhv_CH4` | -  | TgC yr<sup>-1</sup> | - ||
| `D_Fsoil_CH4` | - | TgC yr<sup>-1</sup> | - ||
| `D_Focean_CH4` | - | TgC yr<sup>-1</sup> | - ||
| `D_Fsink_CH4` | $\Delta F_\mathrm{\downarrow}^\mathrm{CH_4}$ | TgC yr<sup>-1</sup> | - ||
| `D_Foxi_CH4` | - | PgC yr<sup>-1</sup> | - ||
||||||
| `D_ewet` | $\Delta e_\mathrm{wet}$ | TgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land` ||
| `D_Awet` | $\Delta A_\mathrm{wet}$ | Mha | `reg_land` ||
| `D_Ewet` | $\Delta E_\mathrm{wet}$ | TgC yr<sup>-1</sup> | `reg_land` ||
||||||
| `D_CH4` | $\Delta \mathrm{CH_4}$ | ppb | - | **yes** |
| `tau_CH4` | - | yr | - ||
| `RF_CH4` | $\Delta \mathrm{RF}^\mathrm{CH_4}$ | W m<sup>-2</sup> | - ||
| `RF_H2Os` | $\Delta \mathrm{RF}^\mathrm{H_2Os}$ | W m<sup>-2</sup> | - ||
||||||
| `D_f_ageair` | - | 1 | - ||
| `f_hv` | $\mathcal{F}_\mathrm{h\nu}$ | 1 | - ||
| `D_Fhv_N2O` | -  | TgN yr<sup>-1</sup> | - ||
| `D_Fsink_N2O` | $\Delta F_\mathrm{\downarrow}^\mathrm{N_2O}$ | TgN yr<sup>-1</sup> | - ||
||||||
| `D_N2O` | $\Delta \mathrm{N_2O}$ | ppb | - | **yes** |
| `tau_N2O` | - | yr | - ||
| `RF_N2O` | $\Delta \mathrm{RF}^\mathrm{N_2O}$ | W m<sup>-2</sup> | - ||
||||||
| `D_Foh_Xhalo` | - | Gg yr<sup>-1</sup> | `spc_halo` ||
| `D_Fhv_CH4` | -  | Gg yr<sup>-1</sup> | `spc_halo` ||
| `D_Fother_CH4` | - | Gg yr<sup>-1</sup> | `spc_halo` ||
| `D_Fsink_CH4` | $\Delta F_\mathrm{\downarrow}^X$ | Gg yr<sup>-1</sup> | `spc_halo` ||
||||||
| `D_Xhalo` | $\Delta X$ | ppt | `spc_halo` | **yes** |
| `RF_Xhalo` | $\Delta \mathrm{RF}^X$ | W m<sup>-2</sup> | `spc_halo` ||
| `RF_halo` | $\Delta \mathrm{RF}^\mathrm{halo}$ | W m<sup>-2</sup> | - ||
||||||
| `D_O3t` | $\Delta \mathrm{O_3t}$ | DU | - ||
| `RF_O3t` | $\Delta \mathrm{RF}^\mathrm{O_3t}$ | W m<sup>-2</sup> | - ||
||||||
| `D_EESC` | $\Delta \mathrm{EESC}$ | ppt | - ||
| `D_O3s` | $\Delta \mathrm{O_3s}$ | DU | - ||
| `RF_O3s` | $\Delta \mathrm{RF}^\mathrm{O_3s}$ | W m<sup>-2</sup> | - ||
||||||
| `D_Edms` | $\Delta E_\mathrm{DMS}$ | TgS yr<sup>-1</sup> | - ||
| `D_Ebvoc` | $\Delta E_\mathrm{BVOC}$ | Tg yr<sup>-1</sup> | - ||
| `D_Edust` | - | Tg yr<sup>-1</sup> | - ||
| `D_Esalt` | - | Tg yr<sup>-1</sup> | - ||
||||||
| `D_SO4` | $\Delta \mathrm{SO_4}$ | Tg | - ||
| `D_POA` | $\Delta \mathrm{POA}$ | Tg | - ||
| `D_BC` | $\Delta \mathrm{BC}$ | Tg | - ||
| `D_NO3` | $\Delta \mathrm{NO_3}$ | Tg | - ||
| `D_SOA` | $\Delta \mathrm{SOA}$ | Tg | - ||
| `D_Mdust` | - | Tg | - ||
| `D_Msalt` | - | Tg | - ||
| `RF_SO4` | $\Delta \mathrm{RF}^\mathrm{SO_4}$ | W m<sup>-2</sup> | - ||
| `RF_POA` | $\Delta \mathrm{RF}^\mathrm{POA}$ | W m<sup>-2</sup> | - ||
| `RF_BC` | $\Delta \mathrm{RF}^\mathrm{BC}$ | W m<sup>-2</sup> | - ||
| `RF_NO3` | $\Delta \mathrm{RF}^\mathrm{NO_3}$ | W m<sup>-2</sup> | - ||
| `RF_SOA` | $\Delta \mathrm{RF}^\mathrm{SOA}$ | W m<sup>-2</sup> | - ||
| `RF_dust` | - | W m<sup>-2</sup> | - ||
| `RF_salt` | - | W m<sup>-2</sup> | - ||
||||||
| `D_AERsol` | $\Delta \mathrm{AER}_\mathrm{sol}$ | Tg | - ||
| `RF_cloud1` | - | W m<sup>-2</sup> | - ||
| `RF_cloud2` | - | W m<sup>-2</sup> | - ||
| `RF_cloud` | $\Delta \mathrm{RF}^\mathrm{cloud}$ | W m<sup>-2</sup> | - ||
||||||
| `RF_BCsnow` | $\Delta \mathrm{RF}^\mathrm{BCsnow}$ | W m<sup>-2</sup> | - ||
||||||
| `RF_lcc` | $\Delta \mathrm{RF}^\mathrm{LCC}$ | W m<sup>-2</sup> | - ||
||||||
| `RF_nonCO2` | - | W m<sup>-2</sup> | - ||
| `RF_wmghg` | - | W m<sup>-2</sup> | - ||
| `RF_strat` | - | W m<sup>-2</sup> | - ||
| `RF_scatter` | - | W m<sup>-2</sup> | - ||
| `RF_absorb` | - | W m<sup>-2</sup> | - ||
| `RF_AERtot` | - | W m<sup>-2</sup> | - ||
| `RF_slcf` | - | W m<sup>-2</sup> | - ||
| `RF_alb` | - | W m<sup>-2</sup> | - ||
| `RF` | $\Delta \mathrm{RF}$ | W m<sup>-2</sup> | - ||
| `RF_warm` | $\Delta \mathrm{RF}_\mathrm{warm}$ | W m<sup>-2</sup> | - ||
| `RF_atm` | $\Delta \mathrm{RF}_\mathrm{atm}$ | W m<sup>-2</sup> | - ||
||||||
| `D_Tg` | $\Delta T_G$ | K | - | **yes** |
| `D_Td` | $\Delta T_D$ | K | - | **yes** |
| `d_Tg` | $\textstyle{\frac{\mathrm{d}}{\mathrm{d}t}}T_G$ | K yr<sup>-1</sup> | - |
| `CFF` | - | W m<sup>-2</sup> K<sup>-1</sup> | - ||
| `D_Tl` | $\Delta T_L$ | K | `reg_land` ||
| `D_To` | $\Delta T_S$ | K | - ||
||||||
| `D_Pg` | $\Delta P_G$ | mm yr<sup>-1</sup> | - | **yes** |
| `D_Pl` | $\Delta P_L$ | mm yr<sup>-1</sup> | `reg_land` ||
||||||
| `D_OHC` | $\Delta \mathrm{OHC}$ | ZJ | - | **yes** |
| `d_OHC` | $\textstyle{\frac{\mathrm{d}}{\mathrm{d}t}}\mathrm{OHC}$ | ZJ yr<sup>-1</sup> | - ||
||||||
| `D_pH` | - | 1 | - ||

&#8224; default; can be altered through the `fct_process_alt` functions.

### Parameters

Parameters are implicitly defined when creating a model's processes. When OSCAR is run, it does not check whether the needed parameters are actually provided in the `Par` argument. Primary parameters can be loaded with the `load_all_param` function defined in `core_fct.fct_loadP`. Many parameters have several possible values, and these different configurations are defined along the various `mod_` dimensions of the dataset containing the primary parameters. Sets of randomly drawn parameters for Monte Carlo runs can be generated using the `generate_config` function defined in `core_fct.fct_genMC`. More information on each parameter is available in `core_fct.fct_loadP`.

| In code | In papers | Units | Dims | Mods
| --- | --- | --- | --- | --- |
| `a_dic` | $\alpha_\mathrm{sol}$ | µmol kg<sup>-1</sup> [ppm m<sup>-3</sup>]<sup>-1</sup> | - | - |
| `mld_0` | $h_{\mathrm{mld},0}$ | m | - | `mod_Focean_struct` |
| `A_ocean` | $A_\mathrm{ocean}$ | m<sup>2</sup> | - | `mod_Focean_struct` |
| `To_0` | $T_{S,0}$ | K | - | `mod_Focean_struct` |
| `v_fg` | $\nu_\mathrm{fg}$ | yr<sup>-1</sup> | - | `mod_Focean_struct` |
| `p_circ` | $\pi_\mathrm{circ}$ | 1 | `box_osurf` | `mod_Focean_struct` |
| `t_circ` | $\tau_\mathrm{circ}$ | yr | `box_osurf` | `mod_Focean_struct` |
||||||
| `pCO2_is_Pade` | - | `bool` | - | `mod_Focean_chem` |
||||||
| `p_mld` | $\pi_\mathrm{mld}$ | 1 | - | `mod_Focean_trans` |
| `g_mld` | $\gamma_\mathrm{mld}$ | K<sup>-1</sup> | - | `mod_Focean_trans` |
||||||
| `fert_is_Log` | - | `bool` | - | `mod_Fland_fert` |
| `t_shift` | $\tau_\mathrm{shift}$ | yr | - | - |
||||||
| `npp_0` | $\eta$ | PgC Mha<sup>-1</sup> yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind` |
| `igni_0` | $\iota$ | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Efire_preind` |
| `harv_0` | $\epsilon_\mathrm{harv}$ | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Eharv_preind` |
| `graz_0` | $\epsilon_\mathrm{graz}$ | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Egraz_preind` |
| `mu1_0` | $\mu_1$ | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind` |
| `mu2_0` | $\mu_2$ | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind` |
| `muM_0` | $\mu_\mathrm{met}$ | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind` |
| `rho1_0` | $\rho_1$ | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind` |
| `rho2_0` | $\rho_2$ | yr<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_preind` |
| `p_agb` | $\pi_\mathrm{agb}$ | 1 | `reg_land, bio_land` | `mod_Eluc_agb` |
||||||
| `b_npp` | $\beta_\mathrm{npp}$ | 1 | `reg_land, bio_land` | `mod_Fland_trans` |
| `b2_npp` | $\tilde{\beta}_\mathrm{npp}$ | 1 | `reg_land, bio_land` | `mod_Fland_trans` |
| `CO2_cp` | $\mathrm{CO_2}_\mathrm{cp}$ | ppm | `reg_land, bio_land` | `mod_Fland_trans` |
| `g_nppT` | $\gamma_{\mathrm{npp},T}$ | K<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_trans, mod_Fland_fert` |
| `g_nppP` | $\gamma_{\mathrm{npp},P}$ | [mm yr<sup>-1</sup>]<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_trans, mod_Fland_fert` |
| `g_rhoT` | $\gamma_{\mathrm{resp},T}$ | K<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_trans, mod_Fland_resp` |
| `g_rhoT1` | $\gamma_{\mathrm{resp},T_1}$ | K<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_trans, mod_Fland_resp` |
| `g_rhoT2` | $\gamma_{\mathrm{resp},T_2}$ | K<sup>-2</sup> | `reg_land, bio_land` | `mod_Fland_trans, mod_Fland_resp` |
| `g_rhoP` | $\gamma_{\mathrm{resp},P}$ | [mm yr<sup>-1</sup>]<sup>-1</sup> | `reg_land, bio_land` | `mod_Fland_trans` |
| `g_igniC` | $\gamma_{\mathrm{igni},C}$ | ppm<sup>-1</sup> | `reg_land, bio_land` | `mod_Efire_trans` |
| `g_igniT` | $\gamma_{\mathrm{igni},T}$ | K<sup>-1</sup> | `reg_land, bio_land` | `mod_Efire_trans` |
| `g_igniP` | $\gamma_{\mathrm{igni},P}$ | [mm yr<sup>-1</sup>]<sup>-1</sup> | `reg_land, bio_land` | `mod_Efire_trans` |
||||||
| `t_hwp` | $\tau_\mathrm{hwp}$ | yr | `box_hwp` | `mod_Ehwp_tau` |
| `w_t_hwp` | - | 1 | - | `mod_Ehwp_speed` |
| `p_hwp_bb` | - | 1 | `box_hwp` | - |
| `p_hwp` | $\pi_\mathrm{hwp}$ | 1 | `reg_land, bio_land, box_hwp` | `mod_Ehwp_bb` |
||||||
| `a_bb` | $\alpha_\mathrm{bb}$ | Tg<i>X</i> PgC<sup>-1</sup> | `reg_land, bio_land, spc_bb` | - |
||||||
| `Cfroz_0` | $C_{\mathrm{froz},0}$ | PgC | `reg_pf` | `mod_Epf_main` |
| `w_clim_pf` | $\omega_{T_\mathrm{pf}}$ | 1 | `reg_pf` | `mod_Epf_main` |
| `g_respT_pf` | $\gamma_{\mathrm{pf},T_1}$ | K<sup>-1</sup> | `reg_pf` | `mod_Epf_main` |
| `g_respT2_pf` | $\gamma_{\mathrm{pf},T_2}$ | K<sup>-2</sup> | `reg_pf` | `mod_Epf_main` |
| `k_resp_pf` | $\kappa_\mathrm{resp,pf}$ | 1 | `reg_pf` | `mod_Epf_main` |
| `pthaw_min` | $p_\mathrm{thaw,min}$ | 1 | `reg_pf` | `mod_Epf_main` |
| `g_pthaw` | $\gamma_{p_\mathrm{thaw}}$ | K<sup>-1</sup> | `reg_pf` | `mod_Epf_main` |
| `k_pthaw` | $\kappa_{p_\mathrm{thaw}}$ | 1 | `reg_pf` | `mod_Epf_main` |
| `v_thaw` | $\nu_\mathrm{thaw}$ | yr<sup>-1</sup> | `reg_pf` | `mod_Epf_main` |
| `v_froz` | $\nu_\mathrm{froz}$ | yr<sup>-1</sup> | `reg_pf` | `mod_Epf_main` |
| `p_pf_thaw` | $\pi_\mathrm{thaw}$ | 1 | `reg_pf, box_thaw` | `mod_Epf_main` |
| `t_pf_thaw` | $\tau_\mathrm{thaw}$ | yr | `reg_pf, box_thaw` | `mod_Epf_main` |
| `p_pf_inst` | - | 1 | - | - |
| `p_pf_CH4` | - | 1 | - | `mod_Epf_CH4` |
||||||
| `ewet_0` | $e_{\mathrm{wet},0}$ | TgC yr<sup>-1</sup> | `reg_land` | `mod_Ewet_preind` |
| `Awet_0` | $A_{\mathrm{wet},0}$ | Mha | `reg_land` | `mod_Ewet_preind` |
| `p_wet` | $\pi_\mathrm{wet}$ | 1 | `reg_land, bio_land` | `mod_Ewet_preind` |
| `g_wetC` | $\gamma_{\mathrm{wet},C}$ | ppm<sup>-1</sup> | `reg_land` | `mod_Awet_trans` |
| `g_wetT` | $\gamma_{\mathrm{wet},T}$ | K<sup>-1</sup> | `reg_land` | `mod_Awet_trans` |
| `g_wetP` | $\gamma_{\mathrm{wet},P}$ | [mm yr<sup>-1</sup>]<sup>-1</sup> | `reg_land` | `mod_Awet_trans` |
||||||
| `a_CO2` | $\alpha_\mathrm{atm}^\mathrm{CO_2}$ | PgC ppm<sup>-1</sup> | - | - |
| `a_CH4` | $\alpha_\mathrm{atm}^\mathrm{CH_4}$ | TgC ppb<sup>-1</sup> | - | - |
| `a_N2O` | $\alpha_\mathrm{atm}^\mathrm{N_2O}$ | TgC ppb<sup>-1</sup> | - | - |
| `a_Xhalo` | $\alpha_\mathrm{atm}^X$ | Gg ppt<sup>-1</sup> | `spc_halo` | - |
| `a_SO4` | - | Tg TgS<sup>-1</sup> | - | - |
| `a_POM` | $\alpha_\mathrm{OM}^\mathrm{OC}$ | Tg TgC<sup>-1</sup> | - | `mod_POA_conv` |
| `a_NO3` | - | Tg TgN<sup>-1</sup> | - | - |
| `CO2_0` | $\mathrm{CO_2}_0$ | ppm | - | - |
| `CH4_0` | $\mathrm{CH_4}_0$ | ppb | - | - |
| `N2O_0` | $\mathrm{N_2O}_0$ | ppb | - | - |
| `Xhalo_0` | $X_0$  | ppt | `spc_halo` | - |
| `p_CH4geo` | - | 1 | - | - |
||||||
| `g_ageair` | $\gamma_\mathrm{age}$ | K<sup>-1</sup> | - | `mod_Fhv_ageair` |
||||||
| `w_t_OH` | - | 1 | - | - |
| `w_t_hv` | - | 1 | - | - |
||||||
| `t_OH_CH4` | $\tau_\mathrm{OH}^\mathrm{CH_4}$ | yr | - | `mod_Foh_tau` |
| `t_hv_CH4` | $\tau_\mathrm{h\nu}^\mathrm{CH_4}$ | yr | - | - |
| `t_soil_CH4` | $\tau_\mathrm{soil}^\mathrm{CH_4}$ | yr | - | - |
| `t_ocean_CH4` | $\tau_\mathrm{ocean}^\mathrm{CH_4}$ | yr | - | - |
||||||
| `x_OH_Ta` | $\chi_\mathrm{T_A}^\mathrm{OH}$ | 1 | - | `mod_Foh_trans` |
| `x_OH_Qa` | $\chi_\mathrm{Q_A}^\mathrm{OH}$ | 1 | - | `mod_Foh_trans` |
| `x_OH_O3s` | $\chi_\mathrm{O_3s}^\mathrm{OH}$ | 1 | - | `mod_Foh_trans` |
| `x_OH_CH4` | $\chi_\mathrm{CH4}^\mathrm{OH}$ | 1 | - | `mod_Foh_trans` |
| `x_OH_NOX` | $\tilde{\chi}_\mathrm{NO_x}^\mathrm{OH}$ | 1 | - | `mod_Foh_trans` |
| `x_OH_CO` | $\tilde{\chi}_\mathrm{CO}^\mathrm{OH}$ | 1 | - | `mod_Foh_trans` |
| `x_OH_VOC` | $\tilde{\chi}_\mathrm{VOC}^\mathrm{OH}$ | 1 | - | `mod_Foh_trans` |
| `x2_OH_NOX` | $\chi_\mathrm{NO_x}^\mathrm{OH}$ | [TgN yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_Foh_trans` |
| `x2_OH_CO` | $\chi_\mathrm{CO}^\mathrm{OH}$ | [TgC yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_Foh_trans` |
| `x2_OH_VOC` | $\chi_\mathrm{VOC}^\mathrm{OH}$ | [Tg yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_Foh_trans` |
| `w_clim_Ta` | $\kappa_\mathrm{T_A}$ | 1 | - | - |
| `k_Qa` | $\kappa_\mathrm{Q_A}$ | 1 | - | - |
| `Ta_0` | $T_{A,0}$ | K | - | - |
| `k_svp` | $\kappa_\mathrm{svp}$ | 1 | - | - |
| `T_svp` | $T_\mathrm{svp}$ | K | - | - |
| `O3s_0` | $\mathrm{O_3s}_0$ | DU | - | - |
| `Enat_NOX` | $E_\mathrm{nat}^\mathrm{NO_x}$ | TgN yr<sup>-1</sup> | - | - |
| `Enat_CO` | $E_\mathrm{nat}^\mathrm{CO}$ | TgC yr<sup>-1</sup> | - | - |
| `Enat_VOC` | $E_\mathrm{nat}^\mathrm{VOC}$ | Tg yr<sup>-1</sup> | - | - |
| `kOH_is_Log` | - | `bool` | - | `mod_Foh_fct` |
||||||
| `t_hv_N2O` | $\tau_\mathrm{h\nu}^\mathrm{N_2O}$ | yr | - | `mod_Fhv_tau` |
||||||
| `x_hv_N2O` | $\chi_\mathrm{N2O}^\mathrm{h\nu}$ | 1 | - | `mod_Fhv_trans` |
| `x_hv_EESC` | $\chi_\mathrm{EESC}^\mathrm{h\nu}$ | 1 | - | `mod_Fhv_trans` |
| `x_hv_age` | $\chi_\mathrm{age}^\mathrm{h\nu}$ | 1 | - | `mod_Fhv_trans` |
||||||
| `t_OH_Xhalo` | $\tau_\mathrm{OH}^X$ | yr | `spc_halo` | - |
| `t_hv_Xhalo` | $\tau_\mathrm{h\nu}^X$ | yr | `spc_halo` | - |
| `t_other_Xhalo` | $\tau_\mathrm{othr}^X$ | yr | `spc_halo` | - |
||||||
| `p_reg_slcf` | $\pi_\mathrm{reg}$ | 1 | `reg_land, reg_slcf` | - |
||||||
| `w_reg_NOX` | $\omega_\mathrm{NO_x}$ | 1 | `reg_slcf` | `mod_O3t_regsat` |
| `w_reg_CO` | $\omega_\mathrm{CO}$ | 1 | `reg_slcf` | `mod_O3t_regsat` |
| `w_reg_VOC` | $\omega_\mathrm{VOC}$ | 1 | `reg_slcf` | `mod_O3t_regsat` |
||||||
| `x_O3t_CH4` | $\chi_\mathrm{CH_4}^\mathrm{O_3t}$ | DU | - | `mod_O3t_emis` |
| `x_O3t_NOX` | $\chi_\mathrm{NO_x}^\mathrm{O_3t}$ | DU [TgN yr<sup>-1</sup>]<sup>-1</sup>  | - | `mod_O3t_emis` |
| `x_O3t_CO` | $\chi_\mathrm{CO}^\mathrm{O_3t}$ | DU [TgC yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_O3t_emis` |
| `x_O3t_VOC` | $\chi_\mathrm{VOC}^\mathrm{O_3t}$ | DU [Tg yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_O3t_emis` |
| `G_O3t` | $\Gamma_\mathrm{O_3t}$ | DU K<sup>-1</sup> | - | `mod_O3t_clim` |
||||||
| `t_lag` | $\tau_\mathrm{lag}$ | yr | - | - |
| `p_fracrel` | $\pi_\mathrm{rel}^X$ | 1 | `spc_halo` | `mod_O3s_fracrel` |
| `k_Br_Cl` | $\alpha_\mathrm{Cl}^\mathrm{Br}$ | 1 | - | - |
| `n_Cl` | $n_\mathrm{Cl}^X$ | 1 | `spc_halo` | - |
| `n_Br` | $n_\mathrm{Br}^X$ | 1 | `spc_halo` | - |
| `EESC_x` | $\mathrm{EESC}_\times$ | ppt | - | - |
| `k_EESC_N2O` | $\frac{\chi_\mathrm{N2O}^\mathrm{O_3s}}{\chi_\mathrm{EESC}^\mathrm{O_3s}}$ | ppt ppb<sup>-1</sup> | - | `mod_O3s_nitrous, mod_O3s_fracrel` |
||||||
| `x_O3s_EESC` | $\chi_\mathrm{EESC}^\mathrm{O_3s}$ | DU ppt<sup>-1</sup> | - | `mod_O3s_trans` |
| `G_O3s` | $\Gamma_\mathrm{O_3s}$ | DU K<sup>-1</sup> | - | `mod_O3s_trans` |
||||||
| `w_reg_SO2` | $\omega_\mathrm{SO_2}$ | 1 | `reg_slcf` | `mod_SO4_regsat` |
| `w_reg_OC` | $\omega_\mathrm{OC}$ | 1 | `reg_slcf` | `mod_POA_regsat` |
| `w_reg_BC` | $\omega_\mathrm{BC}$ | 1 | `reg_slcf` | `mod_BC_regsat` |
||||||
| `t_SO2` | $\tau_\mathrm{SO_2}$ | yr | - | `mod_SO4_load` |
| `t_DMS` | $\tau_\mathrm{DMS}$ | yr | - | `mod_SO4_load` |
| `G_SO4` | $\Gamma_\mathrm{SO_4}$ | Tg K<sup>-1</sup> | - | `mod_SO4_load` |
| `t_OMff` | $\tau_\mathrm{OM,ff}$ | yr | - | `mod_POA_load` |
| `t_OMbb` | $\tau_\mathrm{OM,bb}$ | yr | - | `mod_POA_load` |
| `G_POA` | $\Gamma_\mathrm{POA}$ | Tg K<sup>-1</sup> | - | `mod_POA_load` |
| `t_BCff` | $\tau_\mathrm{BC,ff}$ | yr | - | `mod_BC_load` |
| `t_BCbb` | $\tau_\mathrm{BC,bb}$ | yr | - | `mod_BC_load` |
| `G_BC` | $\Gamma_\mathrm{BC}$ | Tg K<sup>-1</sup> | - | `mod_BC_load` |
| `t_NOX` | $\tau_\mathrm{NO_x}$ | yr | - | `mod_NO3_load` |
| `t_NH3` | $\tau_\mathrm{NH_3}$ | yr | - | `mod_NO3_load` |
| `G_NO3` | $\Gamma_\mathrm{NO_3}$ | Tg K<sup>-1</sup> | - | `mod_NO3_load` |
| `t_VOC` | $\tau_\mathrm{VOC}$ | yr | - | `mod_SOA_load` |
| `t_BVOC` | $\tau_\mathrm{BVOC}$ | yr | - | `mod_SOA_load` |
| `G_SOA` | $\Gamma_\mathrm{SOA}$ | Tg K<sup>-1</sup> | - | `mod_SOA_load` |
| `t_dust` | - | yr | - | `mod_Mdust_load` |
| `G_dust` | - | Tg K<sup>-1</sup> | - | `mod_Mdust_load` |
| `t_salt` | - | yr | - | `mod_Msalt_load` |
| `G_salt` | - | Tg K<sup>-1</sup> | - | `mod_Msalt_load` |
||||||
| `p_sol_SO4` | $\pi_\mathrm{sol}^\mathrm{SO_4}$ | 1 | - | `mod_RFcloud_solub` |
| `p_sol_POA` | $\pi_\mathrm{sol}^\mathrm{POA}$ | 1 | - | `mod_RFcloud_solub` |
| `p_sol_BC` | $\pi_\mathrm{sol}^\mathrm{BC}$ | 1 | - | `mod_RFcloud_solub` |
| `p_sol_NO3` | $\pi_\mathrm{sol}^\mathrm{NO_3}$ | 1 | - | `mod_RFcloud_solub` |
| `p_sol_SOA` | $\pi_\mathrm{sol}^\mathrm{SOA}$ | 1 | - | `mod_RFcloud_solub` |
| `p_sol_dust` | - | 1 | - | `mod_RFcloud_solub` |
| `p_sol_salt` | - | 1 | - | `mod_RFcloud_solub` |
||||||
| `rf_CO2` | $\alpha_\mathrm{rf}^\mathrm{CO_2}$ | W m<sup>-2</sup> | - | - |
| `rf_CH4` | $\alpha_\mathrm{rf}^\mathrm{CH_4}$ | W m<sup>-2</sup> ppb<sup>-0.5</sup> | - | - |
| `rf_N2O` | $\alpha_\mathrm{rf}^\mathrm{N_2O}$ | W m<sup>-2</sup> ppb<sup>-0.5</sup> | - | - |
| `k_rf_H2Os` | $\frac{\alpha_\mathrm{rf}^\mathrm{H_2Os}}{\alpha_\mathrm{rf}^\mathrm{CH_4}}$ | 1 | - | - |
| `rf_Xhalo` | $\alpha_\mathrm{rf}^X$ | W m<sup>-2</sup> ppt<sup>-1</sup> | `spc_halo` | - |
||||||
| `rf_O3t` | $\alpha_\mathrm{rf}^\mathrm{O_3t}$ | W m<sup>-2</sup> DU<sup>-1</sup> | - | `mod_O3t_radeff` |
| `rf_O3s` | $\alpha_\mathrm{rf}^\mathrm{O_3s}$ | W m<sup>-2</sup> DU<sup>-1</sup> | - | `mod_O3s_radeff` |
||||||
| `rf_SO4` | $\alpha_\mathrm{rf}^\mathrm{SO_4}$ | W m<sup>-2</sup> Tg<sup>-1</sup> | - | `mod_SO4_radeff` |
| `rf_POA` | $\alpha_\mathrm{rf}^\mathrm{POA}$ | W m<sup>-2</sup> Tg<sup>-1</sup> | - | `mod_POA_radeff` |
| `rf_BC` | $\alpha_\mathrm{rf}^\mathrm{BC}$ | W m<sup>-2</sup> Tg<sup>-1</sup> | - | `mod_BC_radeff` |
| `rf_NO3` | $\alpha_\mathrm{rf}^\mathrm{NO_3}$ | W m<sup>-2</sup> Tg<sup>-1</sup> | - | `mod_NO3_radeff` |
| `rf_SOA` | $\alpha_\mathrm{rf}^\mathrm{SOA}$ | W m<sup>-2</sup> Tg<sup>-1</sup> | - | `mod_SOA_radeff` |
| `rf_dust` | - | W m<sup>-2</sup> Tg<sup>-1</sup> | - | - |
| `rf_salt` | - | W m<sup>-2</sup> Tg<sup>-1</sup> | - | - |
||||||
| `k_adj_BC` | $\kappa_\mathrm{adj}^\mathrm{BC}$ | 1 | - | `mod_BC_adjust` |
||||||
| `Phi_0` | $\Phi$ | W m<sup>-2</sup> | - | `mod_RFcloud_erf, mod_RFcloud_solub` |
| `AERsol_0` | $\mathrm{AER}_{\mathrm{sol},0}$ | Tg | - | `mod_RFcloud_solub, mod_RFcloud_erf, mod_RFcloud_preind` |
||||||
| `p_reg_bcsnow` | $\pi_\mathrm{reg}$ | 1 | `reg_land, reg_bcsnow` | - |
||||||
| `w_reg_bcsnow` | $\omega_\mathrm{BCsnow}$ | 1 | `reg_bcsnow` | `mod_RFbcsnow_reg` |
| `rf_bcsnow` | $\alpha_\mathrm{rf}^\mathrm{BCsnow}$ | W m<sup>-2</sup> [TgC yr<sup>-1</sup>]<sup>-1</sup> | - | `mod_RFbcsnow_rf` |
||||||
| `p_trans` | $\pi_\mathrm{trans}$ | 1 | - | - |
| `alpha_alb` | $\alpha_\mathrm{alb}$ | 1 | `reg_land, bio_land` | `mod_RFlcc_alb, mod_RFlcc_flux, mod_RFlcc_cover` |
| `F_rsds` | $\phi_\mathrm{rsds}$ | W m<sup>-2</sup> | `reg_land` | `mod_RFlcc_flux` |
||||||
| `w_warm_volc` | $\kappa_\mathrm{warm}^\mathrm{volc}$ | 1 | - | - |
| `w_warm_bcsnow` | $\kappa_\mathrm{warm}^\mathrm{BCsnow}$ | 1 | - | `mod_RFbcsnow_warmeff` |
| `w_warm_lcc` | $\kappa_\mathrm{warm}^\mathrm{LCC}$ | 1 | - | `mod_RFlcc_warmeff` |
|||||||
| `p_atm_CO2` | $\pi_\mathrm{atm}^\mathrm{CO_2}$ | 1 | - | `mod_Pg_radfact` |
| `p_atm_nonCO2` | $\pi_\mathrm{atm}^\mathrm{noCO_2}$ | 1 | - | `mod_Pg_radfact` |
| `p_atm_O3t` | $\pi_\mathrm{atm}^\mathrm{O_3t}$ | 1 | - | `mod_Pg_radfact` |
| `p_atm_strat` | $\pi_\mathrm{atm}^\mathrm{strat}$ | 1 | - | `mod_Pg_radfact` |
| `p_atm_scatter` | $\pi_\mathrm{atm}^\mathrm{scatter}$ | 1 | - | `mod_Pg_radfact` |
| `p_atm_absorb` | $\pi_\mathrm{atm}^\mathrm{absorb}$ | 1 | - | `mod_Pg_radfact` |
| `p_atm_cloud` | $\pi_\mathrm{atm}^\mathrm{cloud}$ | 1 | - | `mod_Pg_radfact` |
| `p_atm_alb` | $\pi_\mathrm{atm}^\mathrm{alb}$ | 1 | - | `mod_Pg_radfact` |
| `p_atm_solar` | $\pi_\mathrm{atm}^\mathrm{solar}$ | 1 | - | `mod_Pg_radfact` |
||||||
| `lambda_0` | $\lambda$ | K [W m<sup>-2</sup>]<sup>-1</sup> | - | `mod_Tg_resp` |
| `Th_g` | $\frac{\tau_{T_G}}{\lambda}$ | yr W m<sup>-2</sup> K<sup>-1</sup> | - | `mod_Tg_resp` |
| `Th_d` | $\frac{\tau_{T_D}}{\lambda}$ | yr W m<sup>-2</sup> K<sup>-1</sup> | - | `mod_Tg_resp` |
| `th_0` | $\frac{\theta}{\lambda}$ | W m<sup>-2</sup> K<sup>-1</sup> | - | `mod_Tg_resp` |
| `e_ohu` | - | 1 | - | - |
| `w_clim_Tl` | $\omega_{T_L}$ | 1 | `reg_land` | `mod_Tl_pattern, mod_Tg_resp` |
| `w_clim_To` | $\omega_{T_S}$ | 1 | - | `mod_Tl_pattern, mod_Tg_resp` |
||||||
| `a_prec` | $\alpha_{P_G}$ | mm yr<sup>-1</sup> K<sup>-1</sup> | - | `mod_Pg_resp` |
| `b_prec` | $\beta_{P_G}$ | mm yr<sup>-1</sup> [W m<sup>-2</sup>]<sup>-1</sup> | - | `mod_Pg_resp` |
| `w_clim_Pl` | $\omega_{P_L}$ | 1 | `reg_land` | `mod_Pl_pattern, mod_Pg_resp` |
||||||
| `p_ohc` | $\pi_\mathrm{ohc}$ | 1 | - | - |
||||||
| `pH_is_Log` | - | `bool` | - | `mod_pH_fct` |
