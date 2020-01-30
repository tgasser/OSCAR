"""
Copyright: IIASA (International Institute for Applied Systems Analysis), 2016-2020; CEA (Commissariat a L'Energie Atomique) & UVSQ (Universite de Versailles et Saint-Quentin), 2016
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at), Yann Quilcaille

This software is a computer program whose purpose is to simulate the behavior of the Earth system, with a specific but not exclusive focus on anthropogenic climate change.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""

##################################################
##################################################

"""
CONTENT
-------
1. CARBON CYCLE
    1.1. Ocean
    1.2. Land
        load_land_TRENDYv7
    1.3. Permafrost
    1.4. Wetlands
2. ATMO CHEMISTRY
3. RADIATIVE FORCING
4. CLIMATE
5. IMPACTS
"""

##################################################
##################################################

import warnings
import numpy as np
import xarray as xr

from core_fct.fct_misc import aggreg_region


##################################################
## 1. CARBON CYCLE
##################################################

##==========
## 1.2. Land
##==========

## preindustrial land carbon-cycle
## calibrated on TRENDYv7 models
def calib_land_TRENDYv7(mod_region, biome_specific_process=True, 
    path_in='calib_data/', path_out='input_data/parameters/', **useless):
    
    ## load original data
    with xr.open_dataset(path_in + 'land_TRENDYv7.nc') as TMP:
        ds = TMP.sel(sim=['S0', 'S4']).sel(weight='area3', drop=True).load()

    ## aggregate over regions
    ds = aggreg_region(ds, mod_region)

    ## separate natural and anthropogenic biomes
    ds_nat = ds.sel(bio_land=['Forest', 'Non-Forest']).sel(sim='S0', drop=True).mean('year')
    ds_ant = ds.sel(bio_land=['Cropland', 'Pasture', 'Urban']).sel(sim='S4', drop=True).sel(year=slice(1990, 2010)).mean('year')
    ds2 = xr.merge([ds_nat, ds_ant]).sel(bio_land=['Forest', 'Non-Forest']+['Cropland', 'Pasture', 'Urban'])

    ## test existing variables per model
    exist = xr.Dataset()
    for var in ['area', 'npp', 'cVeg', 'cLitter', 'cSoil', 'fFire', 'fHarvest', 'fGrazing', 'fVegLitter', 'fVegSoil', 'fLitterSoil', 'rh', 'fDOC', 'cRoot', 'cWood']:
        exist[var] = ds2[var].sum('reg_land', min_count=1).sum('bio_land', min_count=1).notnull() & (ds2[var].sum('reg_land', min_count=1).sum('bio_land', min_count=1) != 0)

    ## test critical variables
    for var in ['npp', 'cVeg', 'cSoil']:
        if exist[var].sum() < len(exist.model):
            raise RuntimeError("'{0}' must be defined for all models!".format(var))

    ## tests whether 3- or 2-box model
    is_3box = exist.cLitter & exist.fVegLitter & exist.fLitterSoil
    
    ## calculate Fsoil2, Fmort and Rh depending on existing data
    fSoilIn = (ds2.fVegSoil.fillna(0) + ds2.fLitterSoil.fillna(0)).where(ds2.fVegSoil.notnull() | ds2.fLitterSoil.notnull())
    fMort = (ds2.fVegLitter.fillna(0) + ds2.fVegSoil.fillna(0)).where(ds2.fVegLitter.notnull() | ds2.fVegSoil.notnull())
    fMort = fMort.where(exist.fVegLitter | exist.fVegSoil, ds2.npp - ds2.fFire.fillna(0) - ds2.fHarvest.fillna(0) - ds2.fGrazing.fillna(0))
    Rh = ds2.rh.where(exist.rh, ds2.npp - ds2.fFire.fillna(0) - ds2.fHarvest.fillna(0) - ds2.fGrazing.fillna(0) - ds2.fDOC.fillna(0))
    cSoilTot = (ds2.cLitter.fillna(0) + ds2.cSoil.fillna(0)).where(ds2.cLitter.notnull() | ds2.cSoil.notnull())

    ## initialisation of final array
    Par = xr.Dataset()

    ## areal net primary productivity
    Par['npp_0'] = ds2.npp / ds2.area

    ## wildfire emission rate
    Par['igni_0'] = (ds2.fFire / ds2.cVeg).where(exist.fFire)

    ## harvest index
    Par['harv_0'] = (ds2.fHarvest / ds2.cVeg).where(exist.fHarvest)

    ## grazing rate
    Par['graz_0'] = (ds2.fGrazing / ds2.cVeg).where(exist.fGrazing)

    ## mortality rates
    Par['mu1_0'] = (ds2.fVegLitter / ds2.cVeg).where(is_3box, 0)
    Par['mu2_0'] = (ds2.fVegSoil.where(exist.fVegSoil, 0) / ds2.cVeg).where(is_3box, 0) + (fMort / ds2.cVeg).where(~is_3box, 0)

    ## metabolization rate
    Par['muM_0'] = (ds2.fLitterSoil / ds2.cLitter).where(is_3box, 0)

    ## respiration rates
    Par['rho1_0'] = ((ds2.fVegLitter - ds2.fLitterSoil) / ds2.cLitter).where(is_3box, 0)
    Par['rho2_0'] = (fSoilIn / ds2.cSoil).where(is_3box, 0) + (Rh / cSoilTot).where(~is_3box, 0)

    ## above-ground biomass fraction
    Par['p_agb'] = 1 - ds2.cRoot / ds2.cVeg

    ## additional processing (some conditional)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')

        ## fill missing anthropogenic biomes:
        ## using Non-Forest parameters (except Urban NPP set to zero)
        not_fill = Par.npp_0.notnull()
        for var in Par:
            Par[var] = Par[var].where(not_fill, Par[var].sel(bio_land='Non-Forest', drop=True))
        Par['npp_0'] = Par.npp_0.where((Par.bio_land != 'Urban') | not_fill, 0)

        ## fill 'Unknown' region if needed:
        ## using average of other regions
        if 'Unknown' in Par.reg_land_long_name:
            for var in Par:
                Par[var] = Par[var].where(Par[var].notnull(), Par[var].where(Par.reg_land_long_name != 'Unknown').mean('reg_land'))

        ## assume biome-specific processes: 
        ## no wildfire on Urban, harvest only on Cropland, grazing only on Pasture
        if biome_specific_process:
            Par['igni_0'] = Par.igni_0.where(Par.bio_land != 'Urban', 0).where(exist.fFire)
            Par['harv_0'] = Par.harv_0.where(Par.bio_land == 'Cropland', 0).where(exist.fHarvest)
            Par['graz_0'] = Par.graz_0.where(Par.bio_land == 'Pasture', 0).where(exist.fGrazing)

    ## add units
    Par['npp_0'].attrs['units'] = 'PgC yr-1 Mha-1'
    Par['p_agb'].attrs['units'] = '1'
    for var in ['igni_0', 'harv_0', 'graz_0', 'mu1_0', 'mu2_0', 'muM_0', 'rho1_0', 'rho2_0']:
        Par[var].attrs['units'] = 'yr-1'

    ## create model option axes
    Par['igni_0'] = Par.igni_0.dropna('model', how='all').rename({'model':'mod_Efire_preind'})
    Par['harv_0'] = Par.harv_0.dropna('model', how='all').rename({'model':'mod_Eharv_preind'})
    Par['graz_0'] = Par.graz_0.dropna('model', how='all').rename({'model':'mod_Egraz_preind'})
    Par['p_agb'] = Par.p_agb.dropna('model', how='all').rename({'model':'mod_Eluc_agb'})
    Par = Par.rename({'model':'mod_Fland_preind'})

    ## create multi-model mean values
    try:
        mmm = xr.Dataset()
        for mod in ['mod_Fland_preind', 'mod_Efire_preind', 'mod_Eharv_preind', 'mod_Egraz_preind', 'mod_Eluc_agb']:
            mmm = xr.merge([mmm, xr.Dataset({var:Par[var].mean(mod) for var in Par if mod in Par[var].dims}).assign_coords({mod:'mean_TRENDYv7'}).expand_dims(mod, -1)])
    except:
        mmm = xr.Dataset()
        mmm = xr.merge([mmm, xr.Dataset({var:Par[var].mean('mod_Fland_preind') for var in Par if 'mod_Fland_preind' in Par[var].dims}).assign_coords(mod_Fland_preind='mean_TRENDYv7').expand_dims('mod_Fland_preind', -1)])
        mmm = xr.merge([mmm, xr.Dataset({var:Par[var].mean('mod_Efire_preind') for var in Par if 'mod_Efire_preind' in Par[var].dims}).assign_coords(mod_Efire_preind='mean_TRENDYv7').expand_dims('mod_Efire_preind', -1)])
        mmm = xr.merge([mmm, xr.Dataset({var:Par[var].mean('mod_Eharv_preind') for var in Par if 'mod_Eharv_preind' in Par[var].dims}).assign_coords(mod_Eharv_preind='mean_TRENDYv7').expand_dims('mod_Eharv_preind', -1)])
        mmm = xr.merge([mmm, xr.Dataset({var:Par[var].mean('mod_Egraz_preind') for var in Par if 'mod_Egraz_preind' in Par[var].dims}).assign_coords(mod_Egraz_preind='mean_TRENDYv7').expand_dims('mod_Egraz_preind', -1)])
        mmm = xr.merge([mmm, xr.Dataset({var:Par[var].mean('mod_Eluc_agb') for var in Par if 'mod_Eluc_agb' in Par[var].dims}).assign_coords(mod_Eluc_agb='mean_TRENDYv7').expand_dims('mod_Eluc_agb', -1)])

    ## create off values
    try:
        off = xr.Dataset()
        for mod in ['mod_Efire_preind', 'mod_Eharv_preind', 'mod_Egraz_preind']:
            off = xr.merge([off, xr.Dataset({var:Par[var].mean(mod) * 0 for var in Par if mod in Par[var].dims}).assign_coords({mod:'off_'}).expand_dims(mod, -1)])
    except:
        off = xr.Dataset()
        off = xr.merge([off, xr.Dataset({var:Par[var].mean('mod_Efire_preind') * 0 for var in Par if 'mod_Efire_preind' in Par[var].dims}).assign_coords(mod_Efire_preind='off_').expand_dims('mod_Efire_preind', -1)])
        off = xr.merge([off, xr.Dataset({var:Par[var].mean('mod_Eharv_preind') * 0 for var in Par if 'mod_Eharv_preind' in Par[var].dims}).assign_coords(mod_Eharv_preind='off_').expand_dims('mod_Eharv_preind', -1)])
        off = xr.merge([off, xr.Dataset({var:Par[var].mean('mod_Egraz_preind') * 0 for var in Par if 'mod_Egraz_preind' in Par[var].dims}).assign_coords(mod_Egraz_preind='off_').expand_dims('mod_Egraz_preind', -1)])

    ## merge, save and return
    Par = Par.combine_first(mmm).combine_first(off)
    Par.to_netcdf(path_out + 'land_TRENDYv7__' + mod_region + '.nc')
    return Par

