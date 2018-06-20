# OSCAR
A compact Earth system model.

## How-to

Download a release. *Read The Fine Manual.*

OSCAR has been developed and run mostly with Python 2.7.x. It relies on very common libraries: `os`, `csv`, `numpy`, `matplotlib` and a few `scipy` functions.

The source code is provided firstly for transparency, and only secondly for dissemination. This means that it is provided as is, and no support of any kind is guaranteed. Hopefully, this will change in the future.


## Changelog

###### v2.2.2
* Fixed: an error in the pre-processing of the `AeroChem_ACCMIP` input data for the `CSIRO-Mk360` configuration. This was causing biased atmospheric lifetimes for POA and BC aerosols under this configuration (and slightly biased ones under the average `mean-ACCMIP` configuration).

###### v2.2.1
* Fixed: an error in the `f_pCO2` functions, causing a too efficient ocean carbon sink under high warming and high atmospheric CO<sub>2</sub>.
* Fixed: an error in the `f_pH` function, causing unrealistic surface ocean pH changes.

##### v2.2
Initial release on GitHub.


## Reference

**v2.2 |** Gasser, T., P. Ciais, O. Boucher, Y. Quilcaille, M. Tortora, L. Bopp & D. Hauglustaine. "The compact Earth system model OSCAR v2.2: description and first results." *Geoscientific Model Development* 10: 271-319 (2017). [doi:10.5194/gmd-10-271-2017](https://doi.org/doi:10.5194/gmd-10-271-2017)
