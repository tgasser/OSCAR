# OSCAR _v2.3.1_
_Manual_ 


## Basic simulation

To make a basic simulation, one needs to execute the main file `OSCAR.py` after what the `OSCAR_lite` function will be available. Execution of the function without any argument will launch a simulation with the options specified in `OSCAR.py`. One can manually change the drivers of OSCAR by specifying their values as arguments of the main function. More generally, a look at the function definition in `OSCAR-fct.py` is strongly advised, as one can see the various optional arguments (and their default value), such as the choice of output variables or the possibility to use OSCAR in a concentration-driven fashion. Some automatic plots can be activated with the plot argument, but they may be outdated. A code example to use OSCAR in a probabilistic setup is provided in the `test_overall.py` file.


## Simulation options

Any simulation with OSCAR should start with executing the main file `OSCAR.py` in which all the options for the simulation are set. The provided file is set to the "default" parameterization of the model -- i.e. a sort of "average" parameterization that simulates a reasonable climate change over the historical period. Even in a probabilistic setup the `OSCAR.py` should be executed once, so as to load everything. Only after should the probabilistic loading/simulating occur (as in `test_overall.py`). The following table describes the various options available:

| Option | Description |
| --- | --- |
| `p` | Number of time-steps within one year of simulation. |
| `fC` | (NOT IMPLEMENTED) |
| `fT` | Turn on or off the climate feedbacks. |
| `dty` | Precision of numerical computation. |
| `PI_1750` | Boolean. If `False`, the years 1700 to 1750 are simulated. This has very little effect on the numerical outputs and it can be used to test the stability of the model. |
| `ind_final` | End year of simulation. |
| `ind_attrib` | (USELESS) |
| `attrib_DRIVERS` | (USELESS) |
| `attrib_FEEDBACKS` | (USELESS) |
| `attrib_ELUCdelta` | (USELESS) |
| `attrib_ELUCampli` | (USELESS) |
| `mod_regionI` | Regional aggregation. New regions can be defined manually as aggregates of GTAP7 regions, but it requires coding and extension of the `#DATA.Regions_GTAP.csv` file. |
| `mod_regionJ` | (USELESS) |
| `mod_sector` | (USELESS) |
| `mod_kindFF` | (USELESS) |
| `mod_kindLUC` | (USELESS) |
| `mod_kindGHG` | (USELESS) |
| `mod_kindCHI` | (USELESS) |
| `mod_kingAER` | (USELESS) |
| `mod_kindRF` | (USELESS) |
| `mod_kindGE` | (USELESS) |
| `mod_biomeSHR` | How the “shrubland” biome is accounted for. |
| `mod_biomeURB` | How the “urban” biome is accounted for. |
| `data_EFF` | Dataset for fossil-fuel emissions. |
| `data_LULCC` | Dataset of land-use and land-cover change drivers. |
| `data_ECH4` | Dataset for methane emissions. |
| `data_EN2O` | Dataset for nitrous oxide emissions. |
| `data_Ehalo` | Dataset for halogenated compounds emissions. |
| `data_ENOX` | Dataset for nitrogen oxides emissions. |
| `data_ECO` | Dataset for carbon monoxide emissions. |
| `data_EVOC` | Dataset for volatile organic compounds emissions |
| `data_ESO2` | Dataset for sulfur dioxide emissions. |
| `data_ENH3` | Dataset for ammonia emissions. |
| `data_EOC` | Dataset for organic carbon emissions. |
| `data_EBC` | Dataset for black carbon emissions. |
| `data_RFant` | Dataset for other anthropogenic RF. |
| `data_RFnat` | Dataset for natural RF. |
| `mod_DATAscen` | How the transition between historical emissions and scenarios is handled. |
| `scen_ALL` | Scenario for all drivers. Has priority over other `scen_` options. |
| `scen_EFF` | Scenario for fossil-fuel emissions. |
| `scen_LULCC` | Scenario of land-use and land-cover change drivers. |
| `scen_ECH4` | Scenario for methane emissions. |
| `scen_EN2O` | Scenario for nitrous oxide emissions. |
| `scen_Ehalo` | Scenario for halogenated compounds emissions. |
| `scen_ENOX` | Scenario for nitrogen oxides emissions. |
| `scen_ECO` | Scenario for carbon monoxide emissions. |
| `scen_EVOC` | Scenario for volatile organic compounds emissions |
| `scen_ESO2` | Scenario for sulfur dioxide emissions. |
| `scen_ENH3` | Scenario for ammonia emissions. |
| `scen_EOC` | Scenario for organic carbon emissions. |
| `scen_EBC` | Scenario for black carbon emissions. |
| `scen_RFant` | Scenario for other anthropogenic RF. |
| `scen_RFnat` | Scenario for natural RF. |
| `mod_OSNKstruct` | Structure and impulse response function of the ocean C-cycle model. |
| `mod_OSNKchem` | Function emulating the carbonate chemistry. |
| `mod_OSNKtrans` | Transient response of the MLD. |
| `mod_LSNKnpp` | Functional form of fertilization effect. |
| `mod_LSNKrho` | Functional form of heterotrophic respiration rate. |
| `mod_LSNKpreind` | Preindustrial land carbon stocks and fluxes. |
| `mod_LSNKtrans` | Transient response of NPP and HR. |
| `mod_LSNKcover` | Preindustrial natural land-cover. |
| `mod_EFIREpreind` | Preindustrial wildfire fluxes. |
| `mod_EFIREtrans` | Transient response of wildfires. |
| `mod_EPFmain` | Permafrost carbon model. |
| `mod_EPFmethane` | Fraction of permafrost carbon emitted as methane. |
| `mod_ELUCagb` | Above-ground biomass fractions. |
| `mod_EHWPbb` | Biomass burning of non-commercial harvested wood products. |
| `mod_EHWPtau` | Time constants for the oxidation of HWPs. |
| `mod_EHWPfct` | Functional form for the oxidation of HWPs. |
| `mod_OHSNKtau` | Lifetime of methane with respect to the OH sink. |
| `mod_OHSNKfct` | Functional form of the OH sink’s response to ozone precursors. |
| `mod_OHSNKtrans` | Transient response of the OH sink. |
| `mod_EWETpreind` | Preindustrial wetlands emissions. |
| `mod_AWETtrans` | Transient response of wetlands area extent. |
| `mod_HVSNKtau` | Lifetime of nitrous oxide with respect to the stratospheric sink. |
| `mod_HVSNKtrans` | Transient response of the stratospheric sink. |
| `mod_HVSNKcirc` | Transient response of the stratospheric age-of-air. |
| `mod_O3Tregsat` | Regionalization of the tropospheric O3 chemistry. |
| `mod_O3Temis` | Transient response of tropospheric O3 to precursors and methane. |
| `mod_O3Tclim` | Transient response of tropospheric O3 to climate. |
| `mod_O3Tradeff` | Radiative efficiency of tropospheric O3. |
| `mod_O3Sfracrel` | Fractional release factors for EESC. |
| `mod_O3Strans` | Transient response of stratospheric O3 to EESC and climate. |
| `mod_O3Snitrous` | Transient response of stratospheric O3 to N2O. |
| `mod_O3Sradeff` | Radiative efficiency of stratospheric O3. |
| `mod_SO4regsat` | Regionalization of SO4 burden. |
| `mod_SO4load` | Transient response of SO4 burden. |
| `mod_SO4radeff` | Radiative efficiency of SO4. |
| `mod_POAconv` | Conversion factor for OC/OM. |
| `mod_POAregsat` | Regionalization of POA burden. |
| `mod_POAload` | Transient response of POA burden. |
| `mod_POAradeff` | Radiative efficiency of POA. |
| `mod_BCregsat` | Regionalization of BC burden. |
| `mod_BCload` | Transient response of BC burden. |
| `mod_BCradeff` | Radiative efficiency of BC. |
| `mod_BCadjust` | Semi-direct effect of BC. |
| `mod_NO3load` | Transient response of NO3 burden. |
| `mod_NO3radeff` | Radiative efficiency of NO3. |
| `mod_SOAload` | Transient response of SOA burden. |
| `mod_SOAradeff` | Radiative efficiency of SOA. |
| `mod_DUSTload` | (NOT DESCRIBED IN PAPER) Transient response of mineral dust burden. |
| `mod_DUSTradeff` | (NOT IMPLEMENTED) |
| `mod_SALTload` | (NOT DESCRIBED IN PAPER) Transient response of sea salt burden. |
| `mod_SALTradeff` | (NOT IMPLEMENTED) |
| `mod_CLOUDsolub` | Solubility fractions for the indirect effect. |
| `mod_CLOUDerf` | Present-day RF used as reference for the indirect effect. |
| `mod_CLOUDpreind` | Adjustment of preindustrial burden of hydrophilic aerosols. |
| `mod_ALBBCreg` | Regionalization of BC deposition on snow. |
| `mod_ALBBCrf` | Radiative efficiency with respect to emissions of BC deposition on snow. |
| `mod_ALBBCwarm` | Warming efficacy of BC deposition on snow. |
| `mod_ALBLCalb` | Albedo climatology for LCC albedo effect. |
| `mod_ALBLCcover` | Land-cover climatology for LCC albedo effect. |
| `mod_ALBLCwarm` | Warming efficacy of LCC albedo effect. |
| `mod_TEMPresp` | Global surface temperature response. |
| `mod_TEMPpattern` | Experiment for pattern scaling of temperature. |
| `mod_PRECresp` | Global precipitations response. |
| `mod_PRECradfact` | Atmospheric fractions of RF for global precipitations. |
| `mod_PRECpattern` | Experiment for pattern scaling of precipitations. |
| `mod_ACIDsurf` | (NOT DESCRIBED IN PAPER) Functional form for surface ocean acidification. |
| `mod_SLR` | (NOT IMPLEMENTED) |


## Load order
After reading the chosen simulation options, the main file execute four other files: `OSCAR-loadD.py` used to load the drivers, `OSCAR-loadP.py` used to load the parameters, `OSCAR-format.py` used to format the drivers, and `OSCAR-fct.py` used to load the model’s main function. When using the model in a probabilistic setup, the two load files can be used separately, albeit with two important limitations: `OSCAR-loadD.py` must be followed by `OSCAR-format.py`, and the two options `data_LULCC` and `mod_LSNKcover` affect both the drivers and the parameters of the model. In a probabilistic setup, it is recommended to reload `OSCAR-fct.py` each time, or to be careful with the drivers of `OSCAR_lite` since it already has default arguments from the initial loading.
The two following sections give details of the structure of the two load files.

### `OSCAR-loadD.py`

&nbsp;A. VECTORS\
&nbsp;&nbsp;&nbsp;&nbsp;A.1. Regions\
&nbsp;&nbsp;&nbsp;&nbsp;A.2. Sectors\
&nbsp;&nbsp;&nbsp;&nbsp;A.3. Kinds\
&nbsp;&nbsp;&nbsp;&nbsp;A.4. Biomes\
&nbsp;&nbsp;&nbsp;&nbsp;A.5. Halo\
&nbsp;1. GREENHOUSE GASES\
&nbsp;&nbsp;&nbsp;&nbsp;1.1. CDIAC\
&nbsp;&nbsp;&nbsp;&nbsp;1.2. EPA\
&nbsp;&nbsp;&nbsp;&nbsp;1.3. EDGAR\
&nbsp;&nbsp;&nbsp;&nbsp;1.4. EDGAR-FT\
&nbsp;&nbsp;&nbsp;&nbsp;1.5. ACCMIP\
&nbsp;&nbsp;&nbsp;&nbsp;1.6. EDGAR-HYDE\
&nbsp;&nbsp;&nbsp;&nbsp;1.7. Stern1998\
&nbsp;&nbsp;&nbsp;&nbsp;1.8. Davidson2009\
&nbsp;&nbsp;&nbsp;&nbsp;1.9. SRES\
&nbsp;&nbsp;&nbsp;&nbsp;1.10. RCP\
&nbsp;&nbsp;&nbsp;&nbsp;1.A. Past Dataset\
&nbsp;&nbsp;&nbsp;&nbsp;1.B. Final Dataset\
&nbsp;&nbsp;&nbsp;&nbsp;1.11. Peters2011\
&nbsp;2. LAND-USE CHANGE\
&nbsp;&nbsp;&nbsp;&nbsp;2.1. LUH1\
&nbsp;&nbsp;&nbsp;&nbsp;2.2. RCP\
&nbsp;&nbsp;&nbsp;&nbsp;2.A. Final Dataset\
&nbsp;3. HALOGENATED COMPOUNDS\
&nbsp;&nbsp;&nbsp;&nbsp;3.1. EDGAR\
&nbsp;&nbsp;&nbsp;&nbsp;3.2. EDGAR-FT\
&nbsp;&nbsp;&nbsp;&nbsp;3.3. CMIP5\
&nbsp;&nbsp;&nbsp;&nbsp;3.4. RCP\
&nbsp;&nbsp;&nbsp;&nbsp;3.A. Past Dataset\
&nbsp;&nbsp;&nbsp;&nbsp;3.B. Final Dataset\
&nbsp;4. SHORT-LIVED SPECIES\
&nbsp;&nbsp;&nbsp;&nbsp;4.1. EDGAR\
&nbsp;&nbsp;&nbsp;&nbsp;4.2. EDGAR-HTAP\
&nbsp;&nbsp;&nbsp;&nbsp;4.3. ACCMIP\
&nbsp;&nbsp;&nbsp;&nbsp;4.4. SRES\
&nbsp;&nbsp;&nbsp;&nbsp;4.5. RCP\
&nbsp;&nbsp;&nbsp;&nbsp;4.A. Past Dataset\
&nbsp;&nbsp;&nbsp;&nbsp;4.B. Final Dataset\
&nbsp;5. RADIATIVE FORCINGS\
&nbsp;&nbsp;&nbsp;&nbsp;5.1. Anthropogenic\
&nbsp;&nbsp;&nbsp;&nbsp;5.2. Natural\
&nbsp;B. FINAL DRIVERS\
&nbsp;&nbsp;&nbsp;&nbsp;B.1. Preindustrial\
&nbsp;&nbsp;&nbsp;&nbsp;B.2. Attribution

### `OSCAR-loadP.py`

&nbsp;1. CARBON DIOXIDE\
&nbsp;&nbsp;&nbsp;&nbsp;1.1. ATMOSPHERE\
&nbsp;&nbsp;&nbsp;&nbsp;1.2. OCEAN\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2.1. Structure\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2.2. Chemistry\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2.3. CMIP5\
&nbsp;&nbsp;&nbsp;&nbsp;1.3. LAND\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3.1. Functions\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3.2. TRENDY\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3.3. Land-Cover\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3.4. CMIP5\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3.5. Permafrost\
&nbsp;&nbsp;&nbsp;&nbsp;1.4. LAND-USE\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.4.1. TRENDY\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.4.2. Wood-Use\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.4.3. GFED\
&nbsp;2. METHANE\
&nbsp;&nbsp;&nbsp;&nbsp;2.1. ATMOSPHERE\
&nbsp;&nbsp;&nbsp;&nbsp;2.2. CHEMISTRY\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.2.1. Lifetimes\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.2.2. Holmes2013\
&nbsp;&nbsp;&nbsp;&nbsp;2.3. LAND\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.3.1. WETCHIMP\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.3.2. Permafrost\
&nbsp;3. NITROUS OXIDE\
&nbsp;&nbsp;&nbsp;&nbsp;3.1. ATMOSPHERE\
&nbsp;&nbsp;&nbsp;&nbsp;3.2. CHEMISTRY\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.2.1. Lifetimes\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.2.2. Prather2015\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.2.3. CCMVal2\
&nbsp;4. HALOGENATED COMPOUNDS\
&nbsp;&nbsp;&nbsp;&nbsp;4.1. ATMOSPHERE\
&nbsp;&nbsp;&nbsp;&nbsp;4.2. CHEMISTRY\
&nbsp;5. OZONE\
&nbsp;&nbsp;&nbsp;&nbsp;5.1. TROPOSPHERE\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5.1.2. HTAP\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5.1.2. ACCMIP\
&nbsp;&nbsp;&nbsp;&nbsp;5.2. STRATOSPHERE\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5.2.1. Chlorine\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5.2.2. CCMVal2\
&nbsp;6. AEROSOLS\
&nbsp;&nbsp;&nbsp;&nbsp;6.1. ATMOSPHERE\
&nbsp;&nbsp;&nbsp;&nbsp;6.2. CHEMISTRY\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6.2.1. HTAP\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6.2.2. ACCMIP\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6.2.3. Nitrates\
&nbsp;7. RADIATIVE FORCING\
&nbsp;&nbsp;&nbsp;&nbsp;7.A. Reconstructions\
&nbsp;&nbsp;&nbsp;&nbsp;7.1. GREENHOUSE GASES\
&nbsp;&nbsp;&nbsp;&nbsp;7.2. OZONE\
&nbsp;&nbsp;&nbsp;&nbsp;7.3. AEROSOLS\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;7.3.1. Direct\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;7.3.2. Indirect\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;7.3.3. Volcanoes\
&nbsp;&nbsp;&nbsp;&nbsp;7.4. ALBEDO\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;7.4.1. Black Carbon\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;7.4.2. Land-Cover\
&nbsp;8. CLIMATE\
&nbsp;&nbsp;&nbsp;&nbsp;8.1. GLOBE\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8.1.A. Reconstructions\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8.1.1. CMIP5\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8.1.2. Precipitations\
&nbsp;&nbsp;&nbsp;&nbsp;8.2. OCEAN\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8.2.A. Reconstructions\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8.2.1. CMIP5\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8.2.2. Heat content\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8.2.3. Acidification\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8.2.4. Sea-level\
&nbsp;&nbsp;&nbsp;&nbsp;8.3. LAND\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8.3.A. Reconstructions\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8.3.1. CMIP5


## Drivers, variables and parameters names and units

Here we provide the correspondence between the notations used in the description paper and the names used in the model’s code, as well as the (implicit) units used in the model. In some cases there is no direct correspondence, because of the way the model is actually coded. Names should remain self-explanatory.

### Drivers

| In papers | In code | Unit |
| --- | --- | --- |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{FF}" /> | `EFF` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{CH_4}" /> | `ECH4` | MtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{N_2O}" /> | `EN2O` | MtN yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{\{HFC\}}" /> | `EHFC` | kt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{\{PFC\}}" /> | `EPFC` | kt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{\{ODS\}}" /> | `EODS` | kt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{NO_x}" /> | `ENOX` | MtN yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{CO}" /> | `ECO` | MtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{VOC}" /> | `EVOC` | Mt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{SO_2}" /> | `ESO2` | MtS yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{NH_3}" /> | `ENH3` | MtN yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{OC}" /> | `EOC` | MtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{BC}" /> | `EBC` | MtC yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\delta\!\,A" /> | `LUC` | Mha yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\delta\!\,H" /> | `HARV` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\delta\!\,S" /> | `SHIFT` | Mha yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\mathrm{RF}_\mathrm{con}" /> | `RFant` | W m<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\mathrm{RF}_\mathrm{volc}" /> | `RFvolc` | W m<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\mathrm{RF}_\mathrm{solar}" /> | `RFsolar` | W m<sup>-2</sup> |

### Variables

| In papers | In code | Unit |
| --- | --- | --- |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{surf}" /> | `D_CSURF` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{dic}" /> | `D_dic` | µmol kg<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,h_\mathrm{mld}" /> | `D_mld` | m |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{in}" /> | `D_FIN` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{out}" /> | `D_FOUT` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{circ}" /> | `D_FCIRC` | GtC yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,c_\mathrm{veg}" /> | `D_cveg` | GtC Mha<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,c_\mathrm{litt}" /> | `D_csoil1` | GtC Mha<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,c_\mathrm{soil}" /> | `D_csoil2` | GtC Mha<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{npp}" /> | `D_npp` | GtC Mha<sup>-1</sup> yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,e_\mathrm{fire}" /> | `D_efire` | GtC Mha<sup>-1</sup> yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,f_\mathrm{mort}" /> | `D_fmort` | GtC Mha<sup>-1</sup> yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{rh}_\mathrm{litt}" /> | `D_rh1` | GtC Mha<sup>-1</sup> yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,f_\mathrm{met}" /> | `D_fmet` | GtC Mha<sup>-1</sup> yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{rh}_\mathrm{soil}" /> | `D_rh2` | GtC Mha<sup>-1</sup> yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,p_\mathrm{thaw}" /> | `pthaw` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\bar{p}_\mathrm{thaw}" /> | `pthaw_bar` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_{\mathrm{thaw},1}" /> | `CTHAW1` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_{\mathrm{thaw},2}" /> | `CTHAW2` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_{\mathrm{thaw},3}" /> | `CTHAW3` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{froz}" /> | `D_CFROZ` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{thaw}" /> | `FTHAW` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_{\mathrm{thaw},1}" /> | `ETHAW1` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_{\mathrm{thaw},2}" /> | `ETHAW2` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_{\mathrm{thaw},3}" /> | `ETHAW3` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{pf}^\mathrm{CO_2}" /> | `EPF_CO2` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{pf}^\mathrm{CH_4}" /> | `EPF_CH4` | MtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{pf}" /> | `EPF` | GtC yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{veg,luc}" /> | `CVEG_luc` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{hwp,luc}^{w=1}" /> | `CHWP1_luc` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{hwp,luc}^{w=2}" /> | `CHWP2_luc` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{hwp,luc}^{w=3}" /> | `CHWP3_luc` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{litt,luc}" /> | `CSOIL1_luc` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,C_\mathrm{soil,luc}" /> | `CSOIL2_luc` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,A" /> | `D_AREA` | Mha |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{CO_2}" /> | `D_CO2` | ppm |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,ocean}" /> | `OSNK` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,land}" /> | `LSNK` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{LUC}" /> | `ELUC` | GtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{CO_2}" /> | `RF_CO2` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{bb}^\mathrm{CH_4}" /> | `D_EBB_CH4` | MtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{bb}^\mathrm{N_2O}" /> | `D_EBB_N2O` | MtN yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{bb}^\mathrm{NO_x}" /> | `D_EBB_NOX` | MtN yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{bb}^\mathrm{CO}" /> | `D_EBB_CO` | MtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{bb}^\mathrm{VOC}" /> | `D_EBB_VOC` | Mt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{bb}^\mathrm{SO_2}" /> | `D_EBB_SO2` | MtS yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{bb}^\mathrm{NH_3}" /> | `D_EBB_NH3` | MtN yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{bb}^\mathrm{OC}" /> | `D_EBB_OC` | MtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{bb}^\mathrm{BC}" /> | `D_EBB_BC` | MtC yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{CH_4}_\mathrm{lag}" /> | `D_CH4_lag` | ppb |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{N_2O}_\mathrm{lag}" /> | `D_N2O_lag` | ppb |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{\{HFC\}}_\mathrm{lag}" /> | `D_HFC_lag` | ppt |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{\{PFC\}}_\mathrm{lag}" /> | `D_PFC_lag` | ppt |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{\{ODS\}}_\mathrm{lag}" /> | `D_ODS_lag` | ppt |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,OH}^\mathrm{CH_4}" /> | `D_OHSNK_CH4` | MtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,h\nu}^\mathrm{CH_4}" /> | `D_HVSNK_CH4` | MtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,othr}^\mathrm{CH_4}" /> | `D_XSNK_CH4` | MtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{oxi,CH_4}}" /> | `D_FOXI_CH4` | GtC yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,A_\mathrm{wet}" /> | `D_AWET` | Mha |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,e_\mathrm{wet}" /> | `D_ewet` | MtC Mha<sup>-1</sup> yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,E_\mathrm{wet}" /> | `D_EWET` | MtC yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{CH_4}" /> | `D_CH4` | ppb |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{CH_4}" /> | `RF_CH4` | W m<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{H_2Os}" /> | `RF_H2Os` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,h\nu}^\mathrm{N_2O}" /> | `D_HVSNK_N2O` | MtN yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{N_2O}" /> | `D_N2O` | ppb |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{N_2O}" /> | `RF_N2O` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,OH}^\mathrm{\{HFC\}}" /> | `D_OHSNK_HFC` | kt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,h\nu}^\mathrm{\{HFC\}}" /> | `D_HVSNK_HFC` | kt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,othr}^\mathrm{\{HFC\}}" /> | `D_XSNK_HFC` | kt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,OH}^\mathrm{\{PFC\}}" /> | `D_OHSNK_PFC` | kt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,h\nu}^\mathrm{\{PFC\}}" /> | `D_HVSNK_PFC` | kt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,othr}^\mathrm{\{PFC\}}" /> | `D_XSNK_PFC` | kt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,OH}^\mathrm{\{ODS\}}" /> | `D_OHSNK_ODS` | kt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,h\nu}^\mathrm{\{ODS\}}" /> | `D_HVSNK_ODS` | kt yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,F_\mathrm{\downarrow\,othr}^\mathrm{\{ODS\}}" /> | `D_XSNK_ODS` | kt yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{\{HFC\}}" /> | `D_HFC` | ppt |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{\{PFC\}}" /> | `D_PFC` | ppt |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{\{ODS\}}" /> | `D_ODS` | ppt |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{halo}" /> | `RF_halo` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{O_3t}" /> | `D_O3t` | DU |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{O_3t}" /> | `RF_O3t` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{EESC}" /> | `D_EESC` | ppt |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{O_3s}" /> | `D_O3s` | DU |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{O_3s}" /> | `RF_O3s` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{SO_4}" /> | `D_SO4` | Tg |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{POA}" /> | `D_POA` | Tg |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{BC}" /> | `D_BC` | Tg |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{NO_3}" /> | `D_NO3` | Tg |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{SOA}" /> | `D_SOA` | Tg |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{SO_4}" /> | `RF_SO4` | W m<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{POA}" /> | `RF_POA` | W m<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{BC}" /> | `RF_BC` | W m<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{NO_3}" /> | `RF_NO3` | W m<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{SOA}" /> | `RF_SOA` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{AER}_\mathrm{sol}" /> | `D_AERh` | Tg |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{cloud}" /> | `RF_cloud` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{BCsnow}" /> | `RF_bcsnow` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}^\mathrm{LCC}" /> | `RF_lcc` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}" /> | `RF` | W m<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}_\mathrm{warm}" /> | `RF_warm` | W m<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{RF}_\mathrm{atm}" /> | `RF_atm` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,T_G" /> | `D_gst` | K |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,T_D" /> | `D_gst0` | K |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,T_S" /> | `D_sst` | K |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,T_L" /> | `D_lst` | K |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,P_G" /> | `D_gyp` | mm yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,P_L" /> | `D_lyp` | mm yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\Delta\!\,\mathrm{OHC}" /> | `D_OHC` | ZJ |

### Parameters

| In papers | In code | Unit |
| --- | --- | --- |
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{fg}" /> | `v_fg` | yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?A_\mathrm{ocean}" /> | `A_ocean` | m<sup>2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{circ}" /> | `p_circ` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{circ}" /> | `tau_circ` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\mathcal{F}_\mathrm{pCO_2}" /> | `f_pCO2` | ppm |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{sol}" /> | `alpha_sol` | µmol kg<sup>-1</sup> [ppm m<sup>-3</sup>]<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{atm}^\mathrm{CO_2}" /> | `alpha_CO2` | GtC ppm<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{mld}" /> | `alpha_mld` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{mld}" /> | `gamma_mld` | K<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\eta" /> | `npp_0` | GtC Mha<sup>-1</sup> yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\mu" /> | `mu_0` | yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\rho_\mathrm{litt}" /> | `rho1_0` | yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\rho_\mathrm{soil}" /> | `rho2_0` | yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{met}" /> | `p_met` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\beta_\mathrm{npp}" /> | `beta_npp` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\tilde{\beta}_\mathrm{npp}" /> | `beta_npp0` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\mathrm{CO_2}_\mathrm{cp}" /> | `CO2_comp` | ppm |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{npp},T}" /> | `gamma_nppT` | K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{npp},P}" /> | `gamma_nppP` | [mm yr<sup>-1</sup>]<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{resp},T}" /> | `gamma_rhoT` | K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{resp},T_1}" /> | `gamma_rhoT1` | K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{resp},T_2}" /> | `gamma_rhoT2` | K<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{resp},P}" /> | `gamma_rhoP` | [mm yr<sup>-1</sup>]<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\iota" /> | `igni_0` | yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{igni},C}" /> | `gamma_igniC` | ppm<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{igni},T}" /> | `gamma_igniT` | K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{igni},P}" /> | `gamma_igniP` | [mm yr<sup>-1</sup>]<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\omega_{T_\mathrm{pf}}}" /> | `w_reg_lstPF` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{pf},T_1}}" /> | `gamma_rhoPF1` | K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{pf},T_2}}" /> | `gamma_rhoPF2` | K<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{resp,pf}}" /> | `w_rhoPF` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?p_\mathrm{thaw,min}}" /> | `pthaw_min` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\kappa_{p_\mathrm{thaw}}}" /> | `k_pthaw` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{p_\mathrm{thaw}}}" /> | `gamma_pthaw` | K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\mathcal{F}_{\nu_\mathrm{pf}}}" /> | `f_v_PF` | yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{thaw}}" /> | `v_thaw` | yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{froz}}" /> | `v_froz` | yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\pi_{\mathrm{pf},1}}" /> | `p_PF1` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_{\mathrm{pf},2}}" /> | `p_PF2` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_{\mathrm{pf},3}}" /> | `p_PF3` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\tau_{\mathrm{pf},1}}" /> | `tau_PF1` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_{\mathrm{pf},2}}" /> | `tau_PF2` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_{\mathrm{pf},3}}" /> | `tau_PF3` | yr |
| <img src="https://latex.codecogs.com/gif.latex?C_{\mathrm{froz},0}}" /> | `CFROZ_0` | GtC |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{pf,CH_4}}" /> | `p_PF_CH4` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{pf,inst}}" /> | `p_PF_inst` | 1 |
||||
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{shift}" /> | `tau_shift` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{agb}" /> | `p_AGB` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{hwp}" /> | `p_HWP` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{hwp,bb}^{w=1}" /> | `p_HWP1_BB` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{hwp}" /> | `tau_HWP` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\mathcal{F}_\mathrm{hwp}" /> | `r_HWP` | 1 |
||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{bb}^\mathrm{CH_4}" /> | `alpha_BB_CH4` | MtC GtC<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{bb}^\mathrm{N_2O}" /> | `alpha_BB_N2O` | MtN GtC<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{bb}^\mathrm{NO_x}" /> | `alpha_BB_NOX` | MtN GtC<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{bb}^\mathrm{CO}" /> | `alpha_BB_CO` | MtC GtC<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{bb}^\mathrm{VOC}" /> | `alpha_BB_VOC` | Mt GtC<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{bb}^\mathrm{SO_2}" /> | `alpha_BB_SO2` | MtS GtC<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{bb}^\mathrm{NH_3}" /> | `alpha_BB_NH3` | MtN GtC<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{bb}^\mathrm{OC}" /> | `alpha_BB_OC` | MtC GtC<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{bb}^\mathrm{BC}" /> | `alpha_BB_BC` | MtC GtC<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{lag}" /> | `tau_lag` | yr |
||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{atm}^\mathrm{CH_4}" /> | `alpha_CH4` | MtC ppb<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{OH}^\mathrm{CH_4}" /> | `tau_CH4_OH` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{h\nu}^\mathrm{CH_4}" /> | `tau_CH4_hv` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{soil}^\mathrm{CH_4}" /> | `tau_CH4_soil` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{ocean}^\mathrm{CH_4}" /> | `tau_CH4_ocean` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{O_3s}^\mathrm{OH}" /> | `chi_OH_O3s` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{T_A}^\mathrm{OH}" /> | `chi_OH_Tatm` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{Q_A}^\mathrm{OH}" /> | `chi_OH_Qatm` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{NO_x}^\mathrm{OH}" /> | `chi_OH_NOX` | [MtN yr<sup>-1</sup>]<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{CO}^\mathrm{OH}" /> | `chi_OH_CO` | [MtC yr<sup>-1</sup>]<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{VOC}^\mathrm{OH}" /> | `chi_OH_VOC` | [Mt yr<sup>-1</sup>]<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\tilde{\chi}_\mathrm{NO_x}^\mathrm{OH}" /> | `chi_OH_NOX` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\tilde{\chi}_\mathrm{CO}^\mathrm{OH}" /> | `chi_OH_CO` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\tilde{\chi}_\mathrm{VOC}^\mathrm{OH}" /> | `chi_OH_VOC` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{T_A}" /> | `k_Tatm` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{Q_A}" /> | `k_Qatm` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{svp}" /> | `k_svp` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?T_\mathrm{svp}" /> | `T_svp` | K |
| <img src="https://latex.codecogs.com/gif.latex?T_{\mathrm{atm},0}" /> | `Tatm_0` | K |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{nat}^\mathrm{NO_x}" /> | `ENOC_oh` | MtN yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{nat}^\mathrm{CO}" /> | `ECO_oh` | MtC yr<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{nat}^\mathrm{VOC}" /> | `EVOC_oh` | Mt yr<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{wet}" /> | `p_wet` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{wet},C}" /> | `gamma_wetC` | ppm<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{wet},T}" /> | `gamma_wetT` | K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_{\mathrm{wet},P}" /> | `gamma_wetP` | [mm yr<sup>-1</sup>]<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\mathcal{F}_\mathrm{over}" /> | `f_overlap` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{atm}^\mathrm{N_2O}" /> | `alpha_N2O` | MtN ppb<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{h\nu}^\mathrm{N_2O}" /> | `tau_N2O_hv` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{N2O}^\mathrm{h\nu}" /> | `chi_hv_N2O` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{EESC}^\mathrm{h\nu}" /> | `chi_hv_EESC` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{age}^\mathrm{h\nu}" /> | `chi_hv_age` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{age}" /> | `gamma_age` | K<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{atm}^\mathrm{\{HFC\}}" /> | `alpha_HFC` | kt ppt<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{atm}^\mathrm{\{PFC\}}" /> | `alpha_PFC` | kt ppt<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{atm}^\mathrm{\{ODS\}}" /> | `alpha_ODS` | kt ppt<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{OH}^\mathrm{\{HFC\}}" /> | `tau_HFC_OH` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{h\nu}^\mathrm{\{HFC\}}" /> | `tau_HFC_hv` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{othr}^\mathrm{\{HFC\}}" /> | `tau_HFC_othr` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{OH}^\mathrm{\{PFC\}}" /> | `tau_PFC_OH` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{h\nu}^\mathrm{\{PFC\}}" /> | `tau_PFC_hv` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{othr}^\mathrm{\{PFC\}}" /> | `tau_PFC_othr` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{OH}^\mathrm{\{ODS\}}" /> | `tau_ODS_OH` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{h\nu}^\mathrm{\{ODS\}}" /> | `tau_ODS_hv` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{othr}^\mathrm{\{ODS\}}" /> | `tau_ODS_othr` | yr |
||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{\{HFC\}}" /> | `radeff_HFC` | W m<sup>-2</sup> ppt<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{\{PFC\}}" /> | `radeff_PFC` | W m<sup>-2</sup> ppt<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{\{ODS\}}" /> | `radeff_ODS` | W m<sup>-2</sup> ppt<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{CH_4}^\mathrm{O_3t}" /> | `chi_O3t_CH4` | DU |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{NO_x}^\mathrm{O_3t}" /> | `chi_O3t_NOX` | DU [MtN yr<sup>-1</sup>]<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{CO}^\mathrm{O_3t}" /> | `chi_O3t_CO` | DU [MtC yr<sup>-1</sup>]<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{VOC}^\mathrm{O_3t}" /> | `chi_O3t_VOC` | DU [Mt yr<sup>-1</sup>]<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{NO_x}" /> | `w_reg_NOX` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{CO}" /> | `w_reg_CO` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{VOC}" /> | `w_reg_VOC` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{reg}" /> | `p_reg` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{O_3t}" /> | `Gamma_O3t` | DU K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{O_3t}" /> | `radeff_O3t` | W m<sup>-2</sup> DU<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{rel}^\mathrm{\{ODS\}}" /> | `f_fracrel` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?n_\mathrm{Cl}^\mathrm{\{ODS\}}" /> | `n_Cl` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?n_\mathrm{Br}^\mathrm{\{ODS\}}" /> | `n_Br` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{Cl}^\mathrm{Br}" /> | `alpha_Br` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{EESC}^\mathrm{O_3s}" /> | `chi_O3s_EESC` | DU ppt<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\chi_\mathrm{N_2O}^\mathrm{O_3s}" /> | `chi_O3s_N2O` | DU ppb<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\mathrm{EESC}_\times" /> | `EESC_x` | ppt |
| <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{O_3s}" /> | `Gamma_O3s` | DU K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{O_3s}" /> | `radeff_O3s` | W m<sup>-2</sup> DU<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{SO_2}" /> | `tau_SO2` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{DMS}" /> | `tau_DMS` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{SO_4}" /> | `Gamma_SO4` | Tg K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{SO_2}" /> | `w_reg_SO2` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{OM,ff}" /> | `tau_OMff` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{OM,bb}" /> | `tau_OMbb` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{POA}" /> | `Gamma_POA` | Tg K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{OC}" /> | `w_reg_OC` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{OM}^\mathrm{OC}" /> | `alpha_POA` | Tg TgC<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{BC,ff}" /> | `tau_BCff` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{BC,bb}" /> | `tau_BCbb` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{BC}" /> | `Gamma_BC` | Tg K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{BC}" /> | `w_reg_BC` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{NO_x}" /> | `tau_NOX` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{NH_3}" /> | `tau_NH3` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{NO_3}" /> | `Gamma_NO3` | Tg K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{VOC}" /> | `tau_VOC` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{BVOC}" /> | `tau_BVOC` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\Gamma_\mathrm{SOA}" /> | `Gamma_SOA` | Tg K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{SO_4}" /> | `radeff_SO4` | W m<sup>-2</sup> Tg<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{POA}" /> | `radeff_POA` | W m<sup>-2</sup> Tg<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{BC}" /> | `radeff_BC` | W m<sup>-2</sup> Tg<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{NO_3}" /> | `radeff_NO3` | W m<sup>-2</sup> Tg<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{SOA}" /> | `radeff_SOA` | W m<sup>-2</sup> Tg<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{adj}^\mathrm{BC}" /> | `k_BC_adjust` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{sol}^\mathrm{SO_4}" /> | `solub_SO4` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{sol}^\mathrm{POA}" /> | `solub_POA` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{sol}^\mathrm{BC}" /> | `solub_BC` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{sol}^\mathrm{NO_3}" /> | `solub_NO3` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{sol}^\mathrm{SOA}" /> | `solub_SOA` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\Phi" /> | `Phi_0` | W m<sup>-2</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\omega_\mathrm{BCsnow}" /> | `w_reg_BCsnow` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{rf}^\mathrm{BCsnow}" /> | `radeff_BCsnow` | W m<sup>-2</sup> [MtC yr<sup>-1</sup>]<sup>-1</sup> |
||||
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{trans}" /> | `p_trans` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\phi_\mathrm{rsds}" /> | `rsds_alb` | W m<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{alb}" /> | `alpha_alb` | 1 |
||||
| <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{warm}^\mathrm{BCsnow}" /> | `warmeff_BCsnow` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{warm}^\mathrm{LCC}" /> | `warmeff_LCC` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{warm}^\mathrm{volc}" /> | `warmeff_volc` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{CO_2}" /> | `p_atm_CO2` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{noCO_2}" /> | `p_atm_noCO2` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{O_3t}" /> | `p_atm_O3t` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{strat}" /> | `p_atm_strat` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{scatter}" /> | `p_atm_scatter` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{absorb}" /> | `p_atm_absorb` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{cloud}" /> | `p_atm_cloud` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{alb}" /> | `p_atm_alb` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{atm}^\mathrm{solar}" /> | `p_atm_solar` | 1 |
||||
| <img src="https://latex.codecogs.com/gif.latex?\lambda" /> | `lambda_0` | K [W m<sup>-2</sup>]<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\tau_{T_G}" /> | `tau_gst` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\tau_{T_D}" /> | `tau_gst0` | yr |
| <img src="https://latex.codecogs.com/gif.latex?\theta" /> | `theta_0` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\omega_{T_S}" /> | `w_reg_sst` | 1 |
| <img src="https://latex.codecogs.com/gif.latex?\omega_{T_L}" /> | `w_reg_lst` | 1 |
||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_{P_G}" /> | `alpha_gyp` | mm yr<sup>-1</sup> K<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\beta_{P_G}" /> | `beta_gyp` | mm yr<sup>-1</sup> [W m<sup>-2</sup>]<sup>-1</sup> |
| <img src="https://latex.codecogs.com/gif.latex?\omega_{P_L}" /> | `w_reg_lyp` | 1 |
||||
| <img src="https://latex.codecogs.com/gif.latex?\pi_\mathrm{ohc}" /> | `p_OHC` | 1 |
