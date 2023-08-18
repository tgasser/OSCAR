# OSCAR
A compact Earth system model.


## How-to

Download a release. *Read The Fine [Manual](./MANUAL.md).*

OSCAR v3 is developed in Python 3.9 and run preferentially through IPython. It makes heavy use of the `xarray` package (v0.20.1), and `netCDF4` for saving data (v1.5.7). It also relies on other common scientific packages: `numpy` (v1.23.3) and `scipy` (v1.9.3). Although the latest versions of the packages used to run OSCAR v3 are given, older versions are likely to work.

The source code is provided firstly for transparency, and only secondly for dissemination. This means that it is provided as is, and no support of any kind is guaranteed (but feel free to ask). Feedbacks and contributions are welcome.


## Changelog

##### v3.3
* Updated: global temperature response parameters to those used in the IPCC AR6 WG1 [Chapter 7](https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-7), based on 35 CMIP6 Earth system models. Note that older parameters remain available.
* Updated: effective radiative forcing functions and parameters for long-lived greenhouse gases, also following the IPCC AR6 WG1 [Chapter 7](https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-7). This also implies update of the preindustrial atmospheric concentrations of these gases.
* Changed: internal solving and saving algorithm to output mid-year average values (instead of end-of-year). Now, to be able to start-up a simulation from an initial state that is the end of another's, one must set the call argument `get_final=True` to request a second output that is the final state of the called simulation.
* Fixed: internal summation of anthropogenic and biomass-burning emissions in `mod_process` that could lead to `NaN` values under some uncommon combination of regional aggregations and land carbon cycle parameterizations.

##### v3.2
* Added: regional aggregation options following the [MESSAGE](https://docs.messageix.org/projects/global/en/latest/overview/spatial.html) integrated assessment model's 11 regions.
* Added: ODR calibration functions (based on `scipy.odr`), to be used for future module calibrations.
* Improved: ocean carbon cycle, with a reworked formulation and structural parameters based on [Bern-SCM](https://doi.org/10.5194/gmd-11-1887-2018), and a new option for the `pCO2` functional form.
* Improved: time-stepping and solving. This includes an auto-initialization of the state variables, and a semi-automatic adjustment of the number of substeps based on the value of `D_CO2`. The latter can be turned off using `adapt_nt=False` when calling the model.
* Removed: the wrapper function (that was mostly confusing and had limited functionality).
* Fixed: the values of `x_OH_NOX`, `x_OH_CO` and `x_OH_VOC` when `mod_Foh_fct == log`. The sensitivity of the OH sink to ozone precursors would be of the opposite sign in half the parameterizations.
* Fixed: an error in the formulation of `p_shift` in the bookkeeping module. It would previously overwrite some minor transitions with `NaN`, leading to a slight underestimation of land-use change emissions.
* Fixed: math display in MANUAL.

###### v3.1.2
* Added: finest possible regional aggregation, given that the parameters are still calibrated following the v2 approach.
* Added: an additional option for the structure of the bookkeping module, in which the effect of land cover change and management are separated. This is added in a temporary file that will be merged with `fct_process_alt` in the next version.
* Fixed: the definition of `Houghton_2017` regions that had Colombia and Guinea swapped. This had minimal impact on the land C cycle.

###### v3.1.1
This version is exactly the one used in the IPCC AR6. A number of functionalities were added, mostly for convenience.
* WARNING: `fct_process` was renamed `mod_process` for consistency, which may affect existing run scripts.
* Added: options in `fct_process_alt` to alter the structure of the bookkeeping module. The default structure is now defined over only one `bio_land` dimension instead of two (`bio_from` and `bio_to`) to reduce the model's memory demand. Mathematically, this is equivalent to the full structure, although extremely small differences can be seen because of the numerical solving. This can be tested using the provided `check_LUC_structure` script.
* Added: an option to add a constant shift and/or some noise to the Monte Carlo parameters, through the new `adjust_config` function in `fct_genMC`.
* Added: the `e_ohu` parameter representing the ocean heat uptake efficacy in the climate module (see e.g. [here](https://doi.org/10.1175/JCLI-D-12-00196.1)). For consistency, this parameter is set to one by default for all configurations.
* Added: new metrics (`D_Flasc`, `d_CO2`, `CFF`, `d_OHC`) in `mod_process` to help diagnose the model ex-post.
* Added: new regional aggregations consistent with the IPCC AR6.
* Added: the `get_IPCC_AR6_parameters` script, to generate a set of prior parameters identical to those used in the IPCC AR6.

##### v3.1
* Improved: land carbon cycle, exactly as described by Gasser et al. (2020). This comes with a recalibration of the module's preindustrial steady-state on TRENDYv7 models.
* Changed: the formulation of `D_ewet` (wetlands areal emissions), to account for the new flexible structure of the land carbon cycle. This reduces the speed at which CH<sub>4</sub> emissions respond to a change in net primary productivity.

###### v3.0.1
* Fixed: an error in the formula of the function describing the overlap of absorption bands between CH<sub>4</sub> and N<sub>2</sub>O, causing significantly wrong RFs. This error appeared during the conversion from v2 to v3, however, and it did not affect the earlier versions of the model.

##### v3.0
The physical equations and parameter values of this version are exactly the same as in v2.4. A few notable changes are worth mentioning here.
* Added: an option to choose the solving scheme of the differential system. The default solving scheme is now an Eulerian exponential integrator (that typically requires fewer sub-timesteps to be stable).
* Removed: the possibility of recalibrating the model's parameters on-the-fly, and especially the regional aggregation of the land carbon cycle. For now, calibrated parameters are directly imported from OSCAR v2, but this feature will be progressively re-implemented has v3 is developed further and new calibrating data become available.
* Changed: the global temperature 2-box model's formulation, to correspond to a more usual formulation found in climate science. This is purely esthetic; this does not change the projected global temperature.
* Changed: the structure of the bookkeeping module for land-use emissions. Specifically, the initialisation is now coded like any other carbon cycle flux. This does not significantly alter the module's performance.

#### v3
This is a complete revamp of OSCAR, for which the code was rewritten from scratch and moved to Python 3. OSCAR v3 relies heavily on the `xarray` package (and therefore netCDF saving format), to better structure the input and output data, to more easily deal with the many dimensions of internal data, and to parallelize simulations so that the run time of a typical Monte Carlo ensemble has been significantly reduced. It also uses its own object classes (`Model` and `Process`) that make it easy to change, extend or tune the model's structure and/or equations. Because of this complete overhaul, manipulation of the model has to be learned from scratch, however.

##### v2.4
This version is the last update of OSCAR v2. It is meant to be as close to OSCAR v3.0 as possible.
* Added: the `mod_biomeV3` option to force biome aggregation to that of v3.0.
* Added: the `mod_EHWPspeed` option that introduces further variations in the decay time of harvested wood products. In addition to the `normal` configuration, the `fast` one scales the value of `tau_HWP` so that 20% of the pool remains after 80% of the initial decay time (a rescale by ~0.5), and the `slow` one scales the value of `tau_HWP` so that 30% of the pool remains after 150% of the initial decay time (a rescale by ~1.25).
* Removed: the `mod_EHWPfct` option, so that only exponential decay of harvested wood products is now possible. This slightly increases CO<sub>2</sub> emissions from land-use change in a Monte Carlo run.
* Removed: the dependency of `p_wet` on the `mod_LSNKcover` option. It is now the average of all possible configurations, which has very little effect on the simulated wetlands CH<sub>4</sub> emissions.

###### v2.3.1
* Added: a new parameter `p_HWP1_BB` quantifying how much of the harvested wood products in pool 1 are actually burnt in the open, and thus accounted for in non-CO<sub>2</sub> anthropogenic biomass burning emissions. It is set to `0.5` to roughly match present-day estimates.
* Added: a new configuration based on `GISS-E2-R-TOMAS` to the `mod_O3Tradeff` option.
* Removed: the `Laube-HL` configuration of the `mod_O3Sfracrel` option, since it makes little physical sense to have parameters specific to high latitudes in a global model like OSCAR.
* Removed: the `Daniel2010-lin` configuration of the `mod_O3Snitrous` option, since a non-saturating effect of N<sub>2</sub>O onto stratospheric O<sub>3</sub> lacked physical ground.
* Fixed: the discretization of `r_HWP` (the response function for harvested wood products). The previous discretization caused delayed and therefore too low emissions.
* Fixed: the value of `k_BC_adjust` for the option `mod_BCadjust == CSIRO`. This has very little impact on a Monte Carlo run.
* Fixed: the value of `radeff_O3t` for the option `mod_O3Tradeff == mean_ACCMIP`. This has no impact on a Monte Carlo run.
* Fixed: the default values of the `beta_npp0` and `CO2_comp` parameters when isolating the urban biome, to prevent `NaN` from appearing during the simulation.
* Fixed: `NaN` values no longer appear in the `alpha_BB` parameters when isolating the urban biome.

##### v2.3
* Added: permafrost carbon thaw and release, exactly as described by Gasser et al. (2018). This comes with a new CO<sub>2</sub> atmospheric flux accounting for the oxidation of geologic CH<sub>4</sub> released in the atmosphere.

###### v2.2.2
* Fixed: an error in the pre-processing of the `AeroChem_ACCMIP` input data for the `CSIRO-Mk360` configuration. This was causing biased atmospheric lifetimes for POA and BC aerosols under this configuration (and slightly biased ones under the average `mean-ACCMIP` configuration).

###### v2.2.1
* Fixed: an error in the `f_pCO2` functions, causing a too efficient ocean carbon sink under high warming and high atmospheric CO<sub>2</sub>.
* Fixed: an error in the `f_pH` function, causing unrealistic surface ocean pH changes.

##### v2.2
Initial release on GitHub. Exact model used by Gasser et al. (2017).


## References

**v3.1 (partial) |** Gasser, T., L. Crepin, Y. Quilcaille, R. A. Houghton, P. Ciais & M. Obersteiner. "Historical CO<sub>2</sub> emissions from land-use and land-cover change and their uncertainty." *Biogeosciences* 17: 4075–4101 (2020). [doi:10.5194/bg-17-4075-2020](https://doi.org/doi:10.5194/bg-17-4075-2020)

**v2.3 (partial) |** Gasser, T., M. Kechiar, P. Ciais, E. J. Burke, T. Kleinen, D. Zhu, Y. Huang, A. Ekici & M. Obersteiner. "Path-dependent reductions in CO<sub>2</sub> emission budgets caused by permafrost carbon release." *Nature Geoscience* 11: 830-835 (2018). [doi:10.1038/s41561-018-0227-0](https://doi.org/doi:10.1038/s41561-018-0227-0)

**v2.2 (full) |** Gasser, T., P. Ciais, O. Boucher, Y. Quilcaille, M. Tortora, L. Bopp & D. Hauglustaine. "The compact Earth system model OSCAR v2.2: description and first results." *Geoscientific Model Development* 10: 271-319 (2017). [doi:10.5194/gmd-10-271-2017](https://doi.org/doi:10.5194/gmd-10-271-2017)
