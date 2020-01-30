# OSCAR
A compact Earth system model.


## How-to

Download a release. *Read The Fine Manual.*

OSCAR v3 is developed in Python 3.7 and run preferentially through IPython. It makes heavy use of the `xarray` package (v0.14.1), and `netCDF4` for saving data (v1.4.2). It also relies on other common scientific packages: `numpy` (v1.16.5) and `scipy` (v1.3.1). Although the latest versions of the packages used to run OSCAR v3 are given, older versions are likely to work.

The source code is provided firstly for transparency, and only secondly for dissemination. This means that it is provided as is, and no support of any kind is guaranteed (but feel free to ask). Feedbacks and contributions are always welcome.


## Changelog

##### v3.0
The physical equations and parameter values of this version are exactly the same as in v2.4. A few notable changes are worth mentioning here.
* Added: an option to choose the solving scheme of the differential system. The default solving scheme is now an Eulerian exponential integrator (that typically requires fewer sub-timesteps to be stable).
* Removed: the possibility of recalibrating the model's parameters on-the-fly, and especially the regional aggregation of the land carbon cycle. For now, calibrated parameters are directly imported from OSCAR v2, but this feature will be progressively re-implemented has v3 is developed further and new calibrating data become available.
* Changed: the global temperature 2-box model's formulation, to correspond to a more usual formulation found in climate science. This is purely esthetic; this does not change the projected global temperature.
* Changed: the structure of the bookkeeping module for land-use emissions. Specifically, the initialisation is now coded like any other carbon-cycle flux. This does not significantly alter the module's performance.

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

**v2.3 (partial) |** : Gasser, T., M. Kechiar, P. Ciais, E. J. Burke, T. Kleinen, D. Zhu, Y. Huang, A. Ekici & M. Obersteiner. "Path-dependent reductions in CO<sub>2</sub> emission budgets caused by permafrost carbon release." *Nature Geoscience* 11: 830-835 (2018). [doi:10.1038/s41561-018-0227-0](https://doi.org/doi:10.1038/s41561-018-0227-0)

**v2.2 (full) |** Gasser, T., P. Ciais, O. Boucher, Y. Quilcaille, M. Tortora, L. Bopp & D. Hauglustaine. "The compact Earth system model OSCAR v2.2: description and first results." *Geoscientific Model Development* 10: 271-319 (2017). [doi:10.5194/gmd-10-271-2017](https://doi.org/doi:10.5194/gmd-10-271-2017)
