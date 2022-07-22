"""
Copyright: IIASA (International Institute for Applied Systems Analysis), 2016-2021; CEA (Commissariat a L'Energie Atomique) & UVSQ (Universite de Versailles et Saint-Quentin), 2014-2016
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
1. LAND-USE CHANGE
    1.1. Split
    1.2. Full
    1.3. Lite
    1.4. Cut
"""

##################################################
##################################################

import numpy as np
import xarray as xr


##################################################
##   1. LAND-USE CHANGE
##################################################

##===========
## 1.1. Split
##===========

## this version keeps the 'bio_from' x 'bio_to' matrix for all variables
## it also splits land cover change and management, as well as slash and soil carbon
## it is the most detailed (and heavy) version

def split_LUC(model_in):
    model = model_in.copy(add_name='_split_LUC')

    ## list LUC processes to be removed
    proc_LUC = ['D_Fveg_bk', 'D_Fsoil1_bk', 'D_Fsoil2_bk', 'D_Fslash1', 'D_Fslash2', 'D_Fhwp', 
        'D_NPP_bk', 'D_Efire_bk', 'D_Eharv_bk', 'D_Egraz_bk', 'D_Fmort1_bk', 'D_Fmort2_bk', 'D_Rh1_bk', 'D_Fmet_bk', 'D_Rh2_bk', 'D_Ehwp',  
        'D_Cveg_bk', 'D_Csoil1_bk', 'D_Csoil2_bk', 'D_Chwp']

    ## remove listed processes
    for proc in proc_LUC:
        if proc in model: del model[proc]

    ## land-use book-keeping initialization of vegetation carbon (LCC)
    model.process('D_Fveg_lcc', ('cveg_0', 'D_cveg', 'd_Acover'), 
        lambda Var, Par: Eq__D_Fveg_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fveg_lcc(Var, Par):
        D_Fveg_bk_lcc = -(Var.cveg_0 + Var.D_cveg).rename({'bio_land':'bio_to'}) * Var.d_Acover
        return D_Fveg_bk_lcc
        
        
    ## land-use book-keeping initialization of vegetation carbon (LU)
    model.process('D_Fveg_lu', ('cveg_0', 'D_cveg', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fveg_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fveg_lu(Var, Par):
        ## ancillary
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        ## main values
        D_Fveg_bk_harv = -0.5*(Var.d_Hwood.rename({'bio_land':'bio_from'}) + Var.d_Hwood.rename({'bio_land':'bio_to'})).where(Var.d_Ashift.bio_from == Var.d_Ashift.bio_to, 0.)
        D_Fveg_bk_shift = -((Var.cveg_0 + Var.D_cveg) * p_shift).rename({'bio_land':'bio_to'}) * Var.d_Ashift
        return D_Fveg_bk_harv + D_Fveg_bk_shift


    ## land-use book-keeping initialization of litter carbon (LCC)
    model.process('D_Fsoil1_lcc', ('csoil1_0', 'D_csoil1', 'd_Acover'), 
        lambda Var, Par: Eq__D_Fsoil1_lcc(Var, Par), 
        units='PgC yr-1')
    
    def Eq__D_Fsoil1_lcc(Var, Par):
        D_Fsoil1_bk_lcc = ((Var.csoil1_0 + Var.D_csoil1).rename({'bio_land':'bio_from'}) - (Var.csoil1_0 + Var.D_csoil1).rename({'bio_land':'bio_to'})) * Var.d_Acover
        return D_Fsoil1_bk_lcc


    ## land-use book-keeping initialization of litter carbon (LU)
    model.process('D_Fsoil1_lu', (), 
        lambda Var, Par: Eq__D_Fsoil1_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fsoil1_lu(Var, Par):
        return 0.


    ## land-use book-keeping initialization of soil carbon (LCC)
    model.process('D_Fsoil2_lcc', ('csoil2_0', 'D_csoil2', 'd_Acover'), 
        lambda Var, Par: Eq__D_Fsoil2_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fsoil2_lcc(Var, Par):
        D_Fsoil2_lcc = ((Var.csoil2_0 + Var.D_csoil2).rename({'bio_land':'bio_from'}) - (Var.csoil2_0 + Var.D_csoil2).rename({'bio_land':'bio_to'})) * Var.d_Acover
        return D_Fsoil2_lcc


    ## land-use book-keeping initialization of soil carbon (LU)
    model.process('D_Fsoil2_lu', (), 
        lambda Var, Par: Eq__D_Fsoil2_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fsoil2_lu(Var, Par):
        return 0.


    ## land-use initial slash flux to litter pool (LCC)
    model.process('D_Fslash1_lcc', ('cveg_0', 'D_cveg', 'd_Acover'), 
        lambda Var, Par: Eq__D_Fslash1_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fslash1_lcc(Var, Par):
        ## ancillary
        p_mu1 = Par.mu1_0 / (Par.mu1_0 + Par.mu2_0)
        ## main values
        D_Fslash_lcc = (p_mu1 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb))).rename({'bio_land':'bio_from'}) * Var.d_Acover
        return D_Fslash_lcc


    ## land-use initial slash flux to litter pool (LU)
    model.process('D_Fslash1_lu', ('cveg_0', 'D_cveg', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fslash1_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fslash1_lu(Var, Par):
        ## ancillary
        p_mu1 = Par.mu1_0 / (Par.mu1_0 + Par.mu2_0)
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        ## main values
        D_Fslash_harv = (p_mu1 * (1 - Par.p_hwp.sum('box_hwp', min_count=1))).rename({'bio_land':'bio_from'}) * 0.5*(Var.d_Hwood.rename({'bio_land':'bio_from'}) + Var.d_Hwood.rename({'bio_land':'bio_to'})).where(Var.d_Ashift.bio_from == Var.d_Ashift.bio_to, 0.)
        D_Fslash_shift = (p_mu1 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb)) * p_shift).rename({'bio_land':'bio_from'}) * Var.d_Ashift
        return D_Fslash_harv + D_Fslash_shift


    ## land-use initial slash flux to soil pool (LCC)
    model.process('D_Fslash2_lcc', ('cveg_0', 'D_cveg', 'd_Acover'), 
        lambda Var, Par: Eq__D_Fslash2_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fslash2_lcc(Var, Par):
        ## ancillary
        p_mu2 = Par.mu2_0 / (Par.mu1_0 + Par.mu2_0)
        ## main values
        D_Fslash_lcc = (p_mu2 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb))).rename({'bio_land':'bio_from'}) * Var.d_Acover
        return D_Fslash_lcc


    ## land-use initial slash flux to soil pool (LU)
    model.process('D_Fslash2_lu', ('cveg_0', 'D_cveg', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fslash2_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fslash2_lu(Var, Par):
        ## ancillary
        p_mu2 = Par.mu2_0 / (Par.mu1_0 + Par.mu2_0)
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        ## main values
        D_Fslash_harv = (p_mu2 * (1 - Par.p_hwp.sum('box_hwp', min_count=1))).rename({'bio_land':'bio_from'}) * 0.5*(Var.d_Hwood.rename({'bio_land':'bio_from'}) + Var.d_Hwood.rename({'bio_land':'bio_to'})).where(Var.d_Ashift.bio_from == Var.d_Ashift.bio_to, 0.)
        D_Fslash_shift = (p_mu2 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb)) * p_shift).rename({'bio_land':'bio_from'}) * Var.d_Ashift
        return D_Fslash_harv + D_Fslash_shift


    ## land-use initial flux to harvested wood product pools (LCC)
    model.process('D_Fhwp_lcc', ('cveg_0', 'D_cveg', 'd_Acover'), 
        lambda Var, Par: Eq__D_Fhwp_lcc(Var, Par), 
        units='PgC yr-1')
    
    def Eq__D_Fhwp_lcc(Var, Par):
        D_Fhwp_lcc = ((Var.cveg_0 + Var.D_cveg) * Par.p_agb * Par.p_hwp).rename({'bio_land':'bio_from'}) * Var.d_Acover
        return D_Fhwp_lcc


    ## land-use initial flux to harvested wood product pools (LU)
    model.process('D_Fhwp_lu', ('cveg_0', 'D_cveg', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fhwp_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fhwp_lu(Var, Par):
        ## ancillary
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        ## main values
        D_Fhwp_harv = Par.p_hwp.rename({'bio_land':'bio_from'}) * 0.5*(Var.d_Hwood.rename({'bio_land':'bio_from'}) + Var.d_Hwood.rename({'bio_land':'bio_to'})).where(Var.d_Ashift.bio_from == Var.d_Ashift.bio_to, 0.)
        D_Fhwp_shift = ((Var.cveg_0 + Var.D_cveg) * Par.p_agb * Par.p_hwp * p_shift).rename({'bio_land':'bio_from'}) * Var.d_Ashift
        return D_Fhwp_harv + D_Fhwp_shift


    ## altered NPP (LCC)
    model.process('D_NPP_lcc', (), 
        lambda Var, Par: Eq__D_NPP_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_NPP_lcc(Var, Par):
        return 0.


    ## altered NPP (LU)
    model.process('D_NPP_lu', (), 
        lambda Var, Par: Eq__D_NPP_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_NPP_lu(Var, Par):
        return 0.

    ## wildfire emissions (LCC)
    model.process('D_Efire_lcc', ('f_igni', 'D_Cveg_lcc'), 
        lambda Var, Par: Eq__D_Efire_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Efire_lcc(Var, Par):
        return (Par.igni_0 * Var.f_igni).rename({'bio_land':'bio_to'}) * Var.D_Cveg_lcc


    ## wildfire emissions (LU)
    model.process('D_Efire_lu', ('f_igni', 'D_Cveg_lu'), 
        lambda Var, Par: Eq__D_Efire_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Efire_lu(Var, Par):
        return (Par.igni_0 * Var.f_igni).rename({'bio_land':'bio_to'}) * Var.D_Cveg_lu


    ## crop harvest emissions (LCC)
    model.process('D_Eharv_lcc', ('D_Cveg_lcc',), 
        lambda Var, Par: Eq__D_Eharv_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Eharv_lcc(Var, Par):
        return Par.harv_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_lcc


    ## crop harvest emissions (LU)
    model.process('D_Eharv_lu', ('D_Cveg_lu',), 
        lambda Var, Par: Eq__D_Eharv_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Eharv_lu(Var, Par):
        return Par.harv_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_lu


    ## pasture grazing emissions (LCC)
    model.process('D_Egraz_lcc', ('D_Cveg_lcc',), 
        lambda Var, Par: Eq__D_Egraz_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Egraz_lcc(Var, Par):
        return Par.graz_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_lcc


    ## pasture grazing emissions (LU)
    model.process('D_Egraz_lu', ('D_Cveg_lu',), 
        lambda Var, Par: Eq__D_Egraz_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Egraz_lu(Var, Par):
        return Par.graz_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_lu


    ## mortality flux to litter (LCC)
    model.process('D_Fmort1_lcc', ('D_Cveg_lcc',), 
        lambda Var, Par: Eq__D_Fmort1_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fmort1_lcc(Var, Par):
        return Par.mu1_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_lcc


    ## mortality flux to litter (LU)
    model.process('D_Fmort1_lu', ('D_Cveg_lu',), 
        lambda Var, Par: Eq__D_Fmort1_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fmort1_lu(Var, Par):
        return Par.mu1_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_lu


    ## mortality flux to soil (LCC)
    model.process('D_Fmort2_lcc', ('D_Cveg_lcc',), 
        lambda Var, Par: Eq__D_Fmort2_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fmort2_lcc(Var, Par):
        return Par.mu2_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_lcc


    ## mortality flux to soil (LU)
    model.process('D_Fmort2_lu', ('D_Cveg_lu',), 
        lambda Var, Par: Eq__D_Fmort2_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fmort2_lu(Var, Par):
        return Par.mu2_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_lu


    ## litter respiration flux (LCC)
    model.process('D_Rh1_lcc', ('f_resp', 'D_Csoil1_lcc'), 
        lambda Var, Par: Eq__D_Rh1_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Rh1_lcc(Var, Par):
        return (Par.rho1_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil1_lcc


    ## litter respiration flux (LU)
    model.process('D_Rh1_lu', ('f_resp', 'D_Csoil1_lu'), 
        lambda Var, Par: Eq__D_Rh1_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Rh1_lu(Var, Par):
        return (Par.rho1_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil1_lu


    ## metabolisation flux (LCC)
    model.process('D_Fmet_lcc', ('f_resp', 'D_Csoil1_lcc'), 
        lambda Var, Par: Eq__D_Fmet_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fmet_lcc(Var, Par):
        return (Par.muM_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil1_lcc


    ## metabolisation flux (LU)
    model.process('D_Fmet_lu', ('f_resp', 'D_Csoil1_lu'), 
        lambda Var, Par: Eq__D_Fmet_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fmet_lu(Var, Par):
        return (Par.muM_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil1_lu


    ## soil respiration flux (LCC)
    model.process('D_Rh2_lcc', ('f_resp', 'D_Csoil2_lcc'), 
        lambda Var, Par: Eq__D_Rh2_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Rh2_lcc(Var, Par):
        return (Par.rho2_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil2_lcc


    ## soil respiration flux (LU)
    model.process('D_Rh2_lu', ('f_resp', 'D_Csoil2_lu'), 
        lambda Var, Par: Eq__D_Rh2_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Rh2_lu(Var, Par):
        return (Par.rho2_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil2_lu


    ## fast slash respiration flux (LCC)
    model.process('D_Rh1S_lcc', ('f_resp', 'D_Cslash1_lcc'), 
        lambda Var, Par: Eq__D_Rh1S_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Rh1S_lcc(Var, Par):
        return (Par.rho1_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Cslash1_lcc


    ## fast slash respiration flux (LU)
    model.process('D_Rh1S_lu', ('f_resp', 'D_Cslash1_lu'), 
        lambda Var, Par: Eq__D_Rh1S_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Rh1S_lu(Var, Par):
        return (Par.rho1_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Cslash1_lu


    ## slash metabolisation flux (LCC)
    model.process('D_FmetS_lcc', ('f_resp', 'D_Cslash1_lcc'), 
        lambda Var, Par: Eq__D_FmetS_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_FmetS_lcc(Var, Par):
        return (Par.muM_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Cslash1_lcc


    ## slash metabolisation flux (LU)
    model.process('D_FmetS_lu', ('f_resp', 'D_Cslash1_lu'), 
        lambda Var, Par: Eq__D_FmetS_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_FmetS_lu(Var, Par):
        return (Par.muM_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Cslash1_lu


    ## slow slash respiration flux (LCC)
    model.process('D_Rh2S_lcc', ('f_resp', 'D_Cslash2_lcc'), 
        lambda Var, Par: Eq__D_Rh2S_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Rh2S_lcc(Var, Par):
        return (Par.rho2_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Cslash2_lcc


    ## slow slash respiration flux (LU)
    model.process('D_Rh2S_lu', ('f_resp', 'D_Cslash2_lu'), 
        lambda Var, Par: Eq__D_Rh2S_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Rh2S_lu(Var, Par):
        return (Par.rho2_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Cslash2_lu
    

    ## harvested wood product oxidation (LCC)
    model.process('D_Ehwp_lcc', ('D_Chwp_lcc',), 
        lambda Var, Par: Eq__D_Ehwp_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Ehwp_lcc(Var, Par):
        return 1 / Par.w_t_hwp / Par.t_hwp * Var.D_Chwp_lcc


    ## harvested wood product oxidation (LU)
    model.process('D_Ehwp_lu', ('D_Chwp_lu',), 
        lambda Var, Par: Eq__D_Ehwp_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Ehwp_lu(Var, Par):
        return 1 / Par.w_t_hwp / Par.t_hwp * Var.D_Chwp_lu


    ## net biome productivity (LCC)
    model.process('D_NBP_lcc', ('D_NPP_lcc', 'D_Efire_lcc', 'D_Eharv_lcc', 'D_Egraz_lcc', 'D_Rh1_lcc', 'D_Rh2_lcc', 'D_Rh1S_lcc', 'D_Rh2S_lcc', 'D_Ehwp_lcc'), 
        lambda Var, Par: Eq__D_NBP_lcc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_NBP_lcc(Var, Par):
        return Var.D_NPP_lcc - Var.D_Efire_lcc - Var.D_Eharv_lcc - Var.D_Egraz_lcc - Var.D_Rh1_lcc - Var.D_Rh2_lcc - Var.D_Rh1S_lcc - Var.D_Rh2S_lcc - Var.D_Ehwp_lcc.sum('box_hwp', min_count=1)


    ## net biome productivity (LU)
    model.process('D_NBP_lu', ('D_NPP_lu', 'D_Efire_lu', 'D_Eharv_lu', 'D_Egraz_lu', 'D_Rh1_lu', 'D_Rh2_lu', 'D_Rh1S_lu', 'D_Rh2S_lu', 'D_Ehwp_lu'), 
        lambda Var, Par: Eq__D_NBP_lu(Var, Par), 
        units='PgC yr-1')

    def Eq__D_NBP_lu(Var, Par):
        return Var.D_NPP_lu - Var.D_Efire_lu - Var.D_Eharv_lu - Var.D_Egraz_lu - Var.D_Rh1_lu - Var.D_Rh2_lu - Var.D_Rh1S_lu - Var.D_Rh2S_lu - Var.D_Ehwp_lu.sum('box_hwp', min_count=1)


    ## net biome productivity (LCC+LU)
    model.process('D_NBP_bk', ('D_NBP_lcc', 'D_NBP_lu'), 
        lambda Var, Par: Eq__D_NBP_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_NBP_bk(Var, Par):
        return Var.D_NBP_lcc + Var.D_NBP_lu


    ## land-use change emissions
    model.process('D_Eluc', ('D_NBP_bk',), 
        lambda Var, Par: Eq__D_Eluc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Eluc(Var, Par):
        return -Var.D_NBP_bk.sum('bio_from', min_count=1).sum('bio_to', min_count=1).sum('reg_land', min_count=1)


    ## PROGNOSTIC: vegetation carbon stock (LCC)
    model.process('D_Cveg_lcc', ('D_Cveg_lcc', 'D_Fveg_lcc', 'D_NPP_lcc', 'D_Efire_lcc', 'D_Eharv_lcc', 'D_Egraz_lcc', 'D_Fmort1_lcc', 'D_Fmort2_lcc'), 
        None, lambda Var, Par: DiffEq__D_Cveg_lcc(Var, Par), lambda Par: vLin__D_Cveg_lcc(Par),
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])

    def DiffEq__D_Cveg_lcc(Var, Par):
        return Var.D_Fveg_lcc + Var.D_NPP_lcc - Var.D_Efire_lcc - Var.D_Eharv_lcc - Var.D_Egraz_lcc - Var.D_Fmort1_lcc - Var.D_Fmort2_lcc

    def vLin__D_Cveg_lcc(Par):
        return (Par.igni_0 + Par.harv_0 + Par.graz_0 + Par.mu1_0 + Par.mu2_0).rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: vegetation carbon stock (LU)
    model.process('D_Cveg_lu', ('D_Cveg_lu', 'D_Fveg_lu', 'D_NPP_lu', 'D_Efire_lu', 'D_Eharv_lu', 'D_Egraz_lu', 'D_Fmort1_lu', 'D_Fmort2_lu'), 
        None, lambda Var, Par: DiffEq__D_Cveg_lu(Var, Par), lambda Par: vLin__D_Cveg_lu(Par),
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])

    def DiffEq__D_Cveg_lu(Var, Par):
        return Var.D_Fveg_lu + Var.D_NPP_lu - Var.D_Efire_lu - Var.D_Eharv_lu - Var.D_Egraz_lu - Var.D_Fmort1_lu - Var.D_Fmort2_lu

    def vLin__D_Cveg_lu(Par):
        return (Par.igni_0 + Par.harv_0 + Par.graz_0 + Par.mu1_0 + Par.mu2_0).rename({'bio_land':'bio_to'})
    

    ## PROGNOSTIC: litter carbon stock (LCC)
    model.process('D_Csoil1_lcc', ('D_Csoil1_lcc', 'D_Fsoil1_lcc', 'D_Fmort1_lcc', 'D_Rh1_lcc', 'D_Fmet_lcc'), 
        None, lambda Var, Par: DiffEq__D_Csoil1_lcc(Var, Par), lambda Par: vLin__D_Csoil1_lcc(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])
    
    def DiffEq__D_Csoil1_lcc(Var, Par):
        return Var.D_Fsoil1_lcc + Var.D_Fmort1_lcc - Var.D_Rh1_lcc - Var.D_Fmet_lcc

    def vLin__D_Csoil1_lcc(Par):
        return (Par.rho1_0 + Par.muM_0).where(Par.rho1_0 + Par.muM_0 != 0, 1E-18).rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: litter carbon stock (LU)
    model.process('D_Csoil1_lu', ('D_Csoil1_lu', 'D_Fsoil1_lu', 'D_Fmort1_lu', 'D_Rh1_lu', 'D_Fmet_lu'), 
        None, lambda Var, Par: DiffEq__D_Csoil1_lu(Var, Par), lambda Par: vLin__D_Csoil1_lu(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])

    def DiffEq__D_Csoil1_lu(Var, Par):
        return Var.D_Fsoil1_lu + Var.D_Fmort1_lu - Var.D_Rh1_lu - Var.D_Fmet_lu

    def vLin__D_Csoil1_lu(Par):
        return (Par.rho1_0 + Par.muM_0).where(Par.rho1_0 + Par.muM_0 != 0, 1E-18).rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: soil carbon stock (LCC)
    model.process('D_Csoil2_lcc', ('D_Csoil2_lcc', 'D_Fsoil2_lcc', 'D_Fmort2_lcc', 'D_Fmet_lcc', 'D_Rh2_lcc'), 
        None, lambda Var, Par: DiffEq__D_Csoil2_lcc(Var, Par), lambda Par: vLin__D_Csoil2_lcc(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])
    
    def DiffEq__D_Csoil2_lcc(Var, Par):
        return Var.D_Fsoil2_lcc + Var.D_Fmort2_lcc + Var.D_Fmet_lcc - Var.D_Rh2_lcc

    def vLin__D_Csoil2_lcc(Par):
        return Par.rho2_0.rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: soil carbon stock (LU)
    model.process('D_Csoil2_lu', ('D_Csoil2_lu', 'D_Fsoil2_lu', 'D_Fmort2_lu', 'D_Fmet_lu', 'D_Rh2_lu'), 
        None, lambda Var, Par: DiffEq__D_Csoil2_lu(Var, Par), lambda Par: vLin__D_Csoil2_lu(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])
    
    def DiffEq__D_Csoil2_lu(Var, Par):
        return Var.D_Fsoil2_lu + Var.D_Fmort2_lu + Var.D_Fmet_lu - Var.D_Rh2_lu

    def vLin__D_Csoil2_lu(Par):
        return Par.rho2_0.rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: fast slash carbon stock (LCC)
    model.process('D_Cslash1_lcc', ('D_Cslash1_lcc', 'D_Fslash1_lcc', 'D_Rh1S_lcc', 'D_FmetS_lcc'), 
        None, lambda Var, Par: DiffEq__D_Cslash1_lcc(Var, Par), lambda Par: vLin__D_Cslash1_lcc(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])
    
    def DiffEq__D_Cslash1_lcc(Var, Par):
        return Var.D_Fslash1_lcc - Var.D_Rh1S_lcc - Var.D_FmetS_lcc

    def vLin__D_Cslash1_lcc(Par):
        return (Par.rho1_0 + Par.muM_0).where(Par.rho1_0 + Par.muM_0 != 0, 1E-18).rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: fast slash carbon stock (LU)
    model.process('D_Cslash1_lu', ('D_Cslash1_lu', 'D_Fslash1_lu', 'D_Rh1S_lu', 'D_FmetS_lu'), 
        None, lambda Var, Par: DiffEq__D_Cslash1_lu(Var, Par), lambda Par: vLin__D_Cslash1_lu(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])

    def DiffEq__D_Cslash1_lu(Var, Par):
        return Var.D_Fslash1_lu - Var.D_Rh1S_lu - Var.D_FmetS_lu

    def vLin__D_Cslash1_lu(Par):
        return (Par.rho1_0 + Par.muM_0).where(Par.rho1_0 + Par.muM_0 != 0, 1E-18).rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: slow slash carbon stock (LCC)
    model.process('D_Cslash2_lcc', ('D_Cslash2_lcc', 'D_Fslash2_lcc', 'D_FmetS_lcc', 'D_Rh2S_lcc'), 
        None, lambda Var, Par: DiffEq__D_Cslash2_lcc(Var, Par), lambda Par: vLin__D_Cslash2_lcc(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])
    
    def DiffEq__D_Cslash2_lcc(Var, Par):
        return Var.D_Fslash2_lcc + Var.D_FmetS_lcc - Var.D_Rh2S_lcc

    def vLin__D_Cslash2_lcc(Par):
        return Par.rho2_0.rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: slow slash carbon stock (LU)
    model.process('D_Cslash2_lu', ('D_Cslash2_lu', 'D_Fslash2_lu', 'D_FmetS_lu', 'D_Rh2S_lu'), 
        None, lambda Var, Par: DiffEq__D_Cslash2_lu(Var, Par), lambda Par: vLin__D_Cslash2_lu(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])
    
    def DiffEq__D_Cslash2_lu(Var, Par):
        return Var.D_Fslash2_lu + Var.D_FmetS_lu - Var.D_Rh2S_lu

    def vLin__D_Cslash2_lu(Par):
        return Par.rho2_0.rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: harvested wood products stock (LCC)
    model.process('D_Chwp_lcc', ('D_Chwp_lcc', 'D_Fhwp_lcc', 'D_Ehwp_lcc'), 
        None, lambda Var, Par: DiffEq__D_Chwp_lcc(Var, Par), lambda Par: vLin__D_Chwp_lcc(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to', 'box_hwp'])

    def DiffEq__D_Chwp_lcc(Var, Par):
        return Var.D_Fhwp_lcc - Var.D_Ehwp_lcc

    def vLin__D_Chwp_lcc(Par):
        return 1 / Par.w_t_hwp / Par.t_hwp


    ## PROGNOSTIC: harvested wood products stock (LU)
    model.process('D_Chwp_lu', ('D_Chwp_lu', 'D_Fhwp_lu', 'D_Ehwp_lu'), 
        None, lambda Var, Par: DiffEq__D_Chwp_lu(Var, Par,), lambda Par: vLin__D_Chwp_lu(Par,), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to', 'box_hwp'])
    
    def DiffEq__D_Chwp_lu(Var, Par):
        return Var.D_Fhwp_lu - Var.D_Ehwp_lu

    def vLin__D_Chwp_lu(Par):
        return 1 / Par.w_t_hwp / Par.t_hwp
    

    ## CO2 emissions from natural biomass burning (= wildfire)
    model.process('D_Efire', ('cveg_0', 'D_efire', 'D_Aland', 'D_Efire_lcc', 'D_Efire_lu'), lambda Var, Par: Eq__D_Efire(Var, Par), units='PgC yr-1')
    def Eq__D_Efire(Var, Par):
        return Par.igni_0 * Var.cveg_0 * Var.D_Aland + Var.D_efire * (Par.Aland_0 + Var.D_Aland) + (Var.D_Efire_lcc + Var.D_Efire_lu).sum('bio_from', min_count=1).rename({'bio_to':'bio_land'})


    ## non-CO2 emissions from anthropogenic biomass burning
    model.process('D_Ebb_ant', ('D_Ehwp_lcc', 'D_Ehwp_lu'), lambda Var, Par: Eq__D_Ebb_ant(Var, Par), units='TgX yr-1')
    def Eq__D_Ebb_ant(Var, Par):
        return Par.a_bb * (Par.p_hwp_bb * (Var.D_Ehwp_lcc + Var.D_Ehwp_lu)).sum('box_hwp', min_count=1).sum('bio_to', min_count=1).rename({'bio_from':'bio_land'})


    ## RETURN
    return model


##==========
## 1.2. Full
##==========

## this version keeps the 'bio_from' x 'bio_to' matrix for all variables
## it is comprehensive but heavier in memory and computation

def full_LUC(model_in):
    model = model_in.copy(add_name='_full_LUC')

    ## land-use book-keeping initialization of vegetation carbon
    model.process('D_Fveg_bk', ('cveg_0', 'D_cveg', 'd_Acover', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fveg_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fveg_bk(Var, Par):
        ## ancillary
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        ## main values
        D_Fveg_bk_lcc = -(Var.cveg_0 + Var.D_cveg).rename({'bio_land':'bio_to'}) * Var.d_Acover
        D_Fveg_bk_harv = -0.5*(Var.d_Hwood.rename({'bio_land':'bio_from'}) + Var.d_Hwood.rename({'bio_land':'bio_to'})).where(Var.d_Acover.bio_from == Var.d_Acover.bio_to, 0.)
        D_Fveg_bk_shift = -((Var.cveg_0 + Var.D_cveg) * p_shift).rename({'bio_land':'bio_to'}) * Var.d_Ashift
        ## summing
        return D_Fveg_bk_lcc + D_Fveg_bk_harv + D_Fveg_bk_shift


    ## land-use book-keeping initialization of litter carbon
    model.process('D_Fsoil1_bk', ('csoil1_0', 'D_csoil1', 'd_Acover'), 
        lambda Var, Par: Eq__D_Fsoil1_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fsoil1_bk(Var, Par):
        D_Fsoil1_bk_lcc = ((Var.csoil1_0 + Var.D_csoil1).rename({'bio_land':'bio_from'}) - (Var.csoil1_0 + Var.D_csoil1).rename({'bio_land':'bio_to'})) * Var.d_Acover
        return D_Fsoil1_bk_lcc


    ## land-use book-keeping initialization of soil carbon
    model.process('D_Fsoil2_bk', ('csoil2_0', 'D_csoil2', 'd_Acover'), 
        lambda Var, Par: Eq__D_Fsoil2_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fsoil2_bk(Var, Par):
        D_Fsoil2_bk_lcc = ((Var.csoil2_0 + Var.D_csoil2).rename({'bio_land':'bio_from'}) - (Var.csoil2_0 + Var.D_csoil2).rename({'bio_land':'bio_to'})) * Var.d_Acover
        return D_Fsoil2_bk_lcc


    ## land-use initial slash flux to litter pool
    model.process('D_Fslash1', ('cveg_0', 'D_cveg', 'd_Acover', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fslash1(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fslash1(Var, Par):
        ## ancillary
        p_mu1 = Par.mu1_0 / (Par.mu1_0 + Par.mu2_0)
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        ## main values
        D_Fslash_lcc = (p_mu1 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb))).rename({'bio_land':'bio_from'}) * Var.d_Acover
        D_Fslash_harv = (p_mu1 * (1 - Par.p_hwp.sum('box_hwp', min_count=1))).rename({'bio_land':'bio_from'}) * 0.5*(Var.d_Hwood.rename({'bio_land':'bio_from'}) + Var.d_Hwood.rename({'bio_land':'bio_to'})).where(Var.d_Acover.bio_from == Var.d_Acover.bio_to, 0.)
        D_Fslash_shift = (p_mu1 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb)) * p_shift).rename({'bio_land':'bio_from'}) * Var.d_Ashift
        ## summing
        return D_Fslash_lcc + D_Fslash_harv + D_Fslash_shift


    ## land-use initial slash flux to soil pool
    model.process('D_Fslash2', ('cveg_0', 'D_cveg', 'd_Acover', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fslash2(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fslash2(Var, Par):
        ## ancillary
        p_mu2 = Par.mu2_0 / (Par.mu1_0 + Par.mu2_0)
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        ## main values
        D_Fslash_lcc = (p_mu2 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb))).rename({'bio_land':'bio_from'}) * Var.d_Acover
        D_Fslash_harv = (p_mu2 * (1 - Par.p_hwp.sum('box_hwp', min_count=1))).rename({'bio_land':'bio_from'}) * 0.5*(Var.d_Hwood.rename({'bio_land':'bio_from'}) + Var.d_Hwood.rename({'bio_land':'bio_to'})).where(Var.d_Acover.bio_from == Var.d_Acover.bio_to, 0.)
        D_Fslash_shift = (p_mu2 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb)) * p_shift).rename({'bio_land':'bio_from'}) * Var.d_Ashift
        ## summing
        return D_Fslash_lcc + D_Fslash_harv + D_Fslash_shift


    ## land-use initial flux to harvested wood product pools
    model.process('D_Fhwp', ('cveg_0', 'D_cveg', 'd_Acover', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fhwp(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fhwp(Var, Par):
        ## ancillary
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        ## main values
        D_Fhwp_lcc = ((Var.cveg_0 + Var.D_cveg) * Par.p_agb * Par.p_hwp).rename({'bio_land':'bio_from'}) * Var.d_Acover
        D_Fhwp_harv = Par.p_hwp.rename({'bio_land':'bio_from'}) * 0.5*(Var.d_Hwood.rename({'bio_land':'bio_from'}) + Var.d_Hwood.rename({'bio_land':'bio_to'})).where(Var.d_Acover.bio_from == Var.d_Acover.bio_to, 0.)
        D_Fhwp_shift = ((Var.cveg_0 + Var.D_cveg) * Par.p_agb * Par.p_hwp * p_shift).rename({'bio_land':'bio_from'}) * Var.d_Ashift
        ## summing
        return D_Fhwp_lcc + D_Fhwp_harv + D_Fhwp_shift


    ## wildfire emissions (under book-keeping)
    model.process('D_Efire_bk', ('f_igni', 'D_Cveg_bk'), 
        lambda Var, Par: Eq__D_Efire_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Efire_bk(Var, Par):
        return (Par.igni_0 * Var.f_igni).rename({'bio_land':'bio_to'}) * Var.D_Cveg_bk


    ## crop harvest emissions (under book-keeping)
    model.process('D_Eharv_bk', ('D_Cveg_bk',), 
        lambda Var, Par: Eq__D_Eharv_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Eharv_bk(Var, Par):
        return Par.harv_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_bk


    ## pasture grazing emissions (under book-keeping)
    model.process('D_Egraz_bk', ('D_Cveg_bk',), 
        lambda Var, Par: Eq__D_Egraz_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Egraz_bk(Var, Par):
        return Par.graz_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_bk


    ## mortality flux to litter (under book-keeping)
    model.process('D_Fmort1_bk', ('D_Cveg_bk',), 
        lambda Var, Par: Eq__D_Fmort1_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fmort1_bk(Var, Par):
        return Par.mu1_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_bk


    ## mortality flux to soil (under book-keeping)
    model.process('D_Fmort2_bk', ('D_Cveg_bk',), 
        lambda Var, Par: Eq__D_Fmort2_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fmort2_bk(Var, Par):
        return Par.mu2_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_bk


    ## litter respiration flux (under book-keeping)
    model.process('D_Rh1_bk', ('f_resp', 'D_Csoil1_bk'), 
        lambda Var, Par: Eq__D_Rh1_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Rh1_bk(Var, Par):
        return (Par.rho1_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil1_bk


    ## metabolisation flux (under book-keeping)
    model.process('D_Fmet_bk', ('f_resp', 'D_Csoil1_bk'), 
        lambda Var, Par: Eq__D_Fmet_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fmet_bk(Var, Par):
        return (Par.muM_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil1_bk


    ## soil respiration flux (under book-keeping)
    model.process('D_Rh2_bk', ('f_resp', 'D_Csoil2_bk'), 
        lambda Var, Par: Eq__D_Rh2_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Rh2_bk(Var, Par):
        return (Par.rho2_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil2_bk


    ## land-use change emissions
    model.process('D_Eluc', ('D_NBP_bk',), 
        lambda Var, Par: Eq__D_Eluc(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Eluc(Var, Par):
        return -Var.D_NBP_bk.sum('bio_from', min_count=1).sum('bio_to', min_count=1).sum('reg_land', min_count=1)


    ## PROGNOSTIC: vegetation carbon stock (under book-keeping)
    model.process('D_Cveg_bk', ('D_Cveg_bk', 'D_Fveg_bk', 'D_NPP_bk', 'D_Efire_bk', 'D_Eharv_bk', 'D_Egraz_bk', 'D_Fmort1_bk', 'D_Fmort2_bk'), 
        None, lambda Var, Par: DiffEq__D_Cveg_bk(Var, Par), lambda Par: vLin__D_Cveg_bk(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])
    
    def DiffEq__D_Cveg_bk(Var, Par):
        return Var.D_Fveg_bk + Var.D_NPP_bk - Var.D_Efire_bk - Var.D_Eharv_bk - Var.D_Egraz_bk - Var.D_Fmort1_bk - Var.D_Fmort2_bk
    
    def vLin__D_Cveg_bk(Par):
        return (Par.igni_0 + Par.harv_0 + Par.graz_0 + Par.mu1_0 + Par.mu2_0).rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: litter carbon stock (under book-keeping)
    model.process('D_Csoil1_bk', ('D_Csoil1_bk', 'D_Fsoil1_bk', 'D_Fslash1', 'D_Fmort1_bk', 'D_Rh1_bk', 'D_Fmet_bk'), 
        None, lambda Var, Par: DiffEq__D_Csoil1_bk(Var, Par), lambda Par: vLin__D_Csoil1_bk(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])
    
    def DiffEq__D_Csoil1_bk(Var, Par):
        return Var.D_Fsoil1_bk + Var.D_Fslash1 + Var.D_Fmort1_bk - Var.D_Rh1_bk - Var.D_Fmet_bk

    def vLin__D_Csoil1_bk(Par):
        return (Par.rho1_0 + Par.muM_0).where(Par.rho1_0 + Par.muM_0 != 0, 1E-18).rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: soil carbon stock (under book-keeping)
    model.process('D_Csoil2_bk', ('D_Csoil2_bk', 'D_Fsoil2_bk', 'D_Fslash2', 'D_Fmort2_bk', 'D_Fmet_bk', 'D_Rh2_bk'), 
        None, lambda Var, Par: DiffEq__D_Csoil2_bk(Var, Par), lambda Par: vLin__D_Csoil2_bk(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])

    def DiffEq__D_Csoil2_bk(Var, Par):
        return Var.D_Fsoil2_bk + Var.D_Fslash2 + Var.D_Fmort2_bk + Var.D_Fmet_bk - Var.D_Rh2_bk

    def vLin__D_Csoil2_bk(Par):
        return Par.rho2_0.rename({'bio_land':'bio_to'})


    ## PROGNOSTIC: harvested wood products stock
    model.process('D_Chwp', ('D_Chwp', 'D_Fhwp', 'D_Ehwp'), 
        None, lambda Var, Par: DiffEq__D_Chwp(Var, Par), lambda Par: vLin__D_Chwp(Par), 
        units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to', 'box_hwp'])

    def DiffEq__D_Chwp(Var, Par):
        return Var.D_Fhwp - Var.D_Ehwp

    def vLin__D_Chwp(Par):
        return 1 / Par.w_t_hwp / Par.t_hwp


    ## CO2 emissions from natural biomass burning (= wildfire)
    model.process('D_Efire', ('cveg_0', 'D_efire', 'D_Aland', 'D_Efire_bk'), 
        lambda Var, Par: Eq__D_Efire__full(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Efire__full(Var, Par):
        return Par.igni_0 * Var.cveg_0 * Var.D_Aland + Var.D_efire * (Par.Aland_0 + Var.D_Aland) + Var.D_Efire_bk.sum('bio_from', min_count=1).rename({'bio_to':'bio_land'})


    ## non-CO2 emissions from anthropogenic biomass burning
    model.process('D_Ebb_ant', ('D_Ehwp',), 
        lambda Var, Par: Eq__D_Ebb_ant__full(Var, Par), 
        units='TgX yr-1')

    def Eq__D_Ebb_ant__full(Var, Par):
        return Par.a_bb * (Par.p_hwp_bb * Var.D_Ehwp).sum('box_hwp', min_count=1).sum('bio_to', min_count=1).rename({'bio_from':'bio_land'})


    ## RETURN
    return model


##==========
## 1.3. Lite
##==========

## this version takes inputs with only one 'bio_land' axis
## it is lighter and faster but at the cost of approximation
## note: 'd_Acover' is split into 'd_Again' and 'd_Aloss'

def lite_LUC(model_in):
    model = model_in.copy(add_name='_lite_LUC')

    ## land-use book-keeping initialization of vegetation carbon
    model.process('D_Fveg_bk', ('cveg_0', 'D_cveg', 'd_Again', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fveg_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fveg_bk(Var, Par):
        ## ancillary
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        ## main values
        D_Fveg_bk_lcc = -(Var.cveg_0 + Var.D_cveg) * Var.d_Again
        D_Fveg_bk_harv = -Var.d_Hwood
        D_Fveg_bk_shift = -(Var.cveg_0 + Var.D_cveg) * p_shift * Var.d_Ashift
        ## summing
        return D_Fveg_bk_lcc + D_Fveg_bk_harv + D_Fveg_bk_shift


    ## land-use book-keeping initialization of litter carbon
    model.process('D_Fsoil1_bk', ('csoil1_0', 'D_csoil1', 'd_Again', 'd_Aloss'), 
        lambda Var, Par: Eq__D_Fsoil1_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fsoil1_bk(Var, Par):
        ## ancillary
        p_gain = (Var.d_Again / Var.d_Again.sum('bio_land', min_count=1)).where(Var.d_Again.sum('bio_land', min_count=1) != 0, 0)
        ## main values
        D_Fsoil1_bk_lcc = ((Var.csoil1_0 + Var.D_csoil1) * Var.d_Aloss).sum('bio_land', min_count=1) * p_gain - (Var.csoil1_0 + Var.D_csoil1) * Var.d_Again
        return D_Fsoil1_bk_lcc


    ## land-use book-keeping initialization of soil carbon
    model.process('D_Fsoil2_bk', ('csoil2_0', 'D_csoil2', 'd_Again', 'd_Aloss'), 
        lambda Var, Par: Eq__D_Fsoil2_bk(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fsoil2_bk(Var, Par):
        ## ancillary
        p_gain = (Var.d_Again / Var.d_Again.sum('bio_land', min_count=1)).where(Var.d_Again.sum('bio_land', min_count=1) != 0, 0)
        ## main values
        D_Fsoil2_bk_lcc = ((Var.csoil2_0 + Var.D_csoil2) * Var.d_Aloss).sum('bio_land', min_count=1) * p_gain - (Var.csoil2_0 + Var.D_csoil2) * Var.d_Again
        return D_Fsoil2_bk_lcc


    ## land-use initial slash flux to litter pool
    model.process('D_Fslash1', ('cveg_0', 'D_cveg', 'd_Again', 'd_Aloss', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fslash1(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fslash1(Var, Par):
        ## ancillary
        p_mu1 = Par.mu1_0 / (Par.mu1_0 + Par.mu2_0)
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        p_gain = (Var.d_Again / Var.d_Again.sum('bio_land', min_count=1)).where(Var.d_Again.sum('bio_land', min_count=1) != 0, 0)
        ## main values
        D_Fslash_lcc = (p_mu1 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb)) * Var.d_Aloss).sum('bio_land', min_count=1) * p_gain
        D_Fslash_harv = p_mu1 * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) * Var.d_Hwood
        D_Fslash_shift = p_mu1 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb)) * p_shift * Var.d_Ashift
        ## summing
        return D_Fslash_lcc + D_Fslash_harv + D_Fslash_shift


    ## land-use initial slash flux to soil pool
    model.process('D_Fslash2', ('cveg_0', 'D_cveg', 'd_Again', 'd_Aloss', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fslash2(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fslash2(Var, Par):
        ## ancillary
        p_mu2 = Par.mu2_0 / (Par.mu1_0 + Par.mu2_0)
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        p_gain = (Var.d_Again / Var.d_Again.sum('bio_land', min_count=1)).where(Var.d_Again.sum('bio_land', min_count=1) != 0, 0)
        ## main values
        D_Fslash_lcc = (p_mu2 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb)) * Var.d_Aloss).sum('bio_land', min_count=1) * p_gain
        D_Fslash_harv = p_mu2 * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) * Var.d_Hwood
        D_Fslash_shift = p_mu2 * (Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb)) * p_shift * Var.d_Ashift
        ## summing
        return D_Fslash_lcc + D_Fslash_harv + D_Fslash_shift


    ## land-use initial flux to harvested wood product pools
    model.process('D_Fhwp', ('cveg_0', 'D_cveg', 'd_Again', 'd_Aloss', 'd_Hwood', 'd_Ashift'), 
        lambda Var, Par: Eq__D_Fhwp(Var, Par), 
        units='PgC yr-1')

    def Eq__D_Fhwp(Var, Par):
        ## ancillary
        p_shift = 1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift).where(Var.cveg_0 != 0, 0)
        ## main values
        D_Fhwp_lcc = (Var.cveg_0 + Var.D_cveg) * Par.p_agb * Par.p_hwp * Var.d_Aloss
        D_Fhwp_harv = Par.p_hwp * Var.d_Hwood
        D_Fhwp_shift = (Var.cveg_0 + Var.D_cveg) * Par.p_agb * Par.p_hwp * p_shift * Var.d_Ashift
        ## summing
        return D_Fhwp_lcc + D_Fhwp_harv + D_Fhwp_shift


    ## PROGNOSTIC: biome area
    model.process('D_Aland', ('D_Aland', 'd_Again', 'd_Aloss'), 
        None, lambda Var, Par: DiffEq__D_Aland(Var, Par), lambda Par: vLin__D_Aland(Par), 
        units='Mha', core_dims=['reg_land', 'bio_land'])

    def DiffEq__D_Aland(Var, Par):
        return Var.d_Again - Var.d_Aloss

    def vLin__D_Aland(Par):
        return 1E-18


    ## RETURN
    return model


##=========
## 1.4. Cut
##=========

## this version entirely removes the bookkeeping and the endogenous biomass burning
## it is much faster but makes the carbon cycle inconsistent (i.e. non-conservative)
## note: it requires an exogenous 'Eluc' forcing (similar to 'Eff')
## note: biome area extent keep changing unless 'd_Acover' is set to zero

def cut_LUC(model_in):
    model = model_in.copy(add_name='_cut_LUC')

    ## list LUC & BB processes to be removed
    proc_LUC = ['D_Fveg_bk', 'D_Fsoil1_bk', 'D_Fsoil2_bk', 'D_Fslash1', 'D_Fslash2', 'D_Fhwp', 
        'D_NPP_bk', 'D_Efire_bk', 'D_Eharv_bk', 'D_Egraz_bk', 'D_Fmort1_bk', 'D_Fmort2_bk', 'D_Rh1_bk', 'D_Fmet_bk', 'D_Rh2_bk', 'D_Ehwp', 'D_NBP_bk', 
        'D_Cveg_bk', 'D_Csoil1_bk', 'D_Csoil2_bk', 'D_Chwp']
    proc_BB = ['D_Efire', 'D_Ebb_nat', 'D_Ebb_ant']

    ## remove listed processes
    for proc in proc_LUC + proc_BB:
        if proc in model: del model[proc]

    ## non-CO2 total emissions from biomass burning (kept for simplicity but set to zero)
    model.process('D_Ebb', (), 
        lambda Var, Par: Eq__D_Ebb(Var, Par), 
        units='TgX yr-1')

    def Eq__D_Ebb(Var, Par):
        return xr.zeros_like(Par.a_bb)

    ## land-use change emissions
    model.process('D_Eluc', ('Eluc',), 
        lambda Var, Par: Eq__D_Eluc(Var, Par), 
        units='PgC yr-1')
        
    def Eq__D_Eluc(Var, Par):
        return Var.Eluc.sum('reg_land', min_count=1)

    ## RETURN
    return model

