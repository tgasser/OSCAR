"""
Copyright: IIASA (International Institute for Applied Systems Analysis), 2016-2020; CEA (Commissariat a L'Energie Atomique) & UVSQ (Universite de Versailles et Saint-Quentin), 2014-2016
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
A. MAIN MODEL
1. CARBON DIOXIDE
    1.1. Ocean
    1.2. Land: intensive
    1.3. Land: extensive
    1.4. Permafrost
    1.5. Atmosphere
2. NON-CO2 SPECIES
    2.1. Biomass burning
    2.2. Lagged concentrations
3. METHANE
    3.1. Atmo chemistry
    3.2. Wetlands
    3.3. Atmosphere
4. NITROUS OXIDE
    4.1. Atmo chemistry
    4.2. Atmosphere
5. HALOGENATED COMPOUNDS
    5.1. Atmo chemistry
    5.2. Atmosphere
6. OZONE
    6.1. Troposphere
    6.2. Stratosphere
7. AEROSOLS
    7.1. Natural emissions
    7.2. Direct effect
    7.3. Cloud effects
8. SURFACE ALBEDO
    8.1. Black carbon on snow
    8.2. Land-cover change
9. CLIMATE
    9.1. Radiative forcing
    9.2. Temperature
    9.3. Precipitation
    9.4. Ocean heat content
10. IMPACTS
    10.1. Acidification
B. SUB-MODELS
"""

##################################################
##################################################

import numpy as np
import xarray as xr

from core_fct.cls_main import Model
from core_fct.fct_misc import Int_ExpInt as Int_dflt


##################################################
##   A. MAIN MODEL
##################################################

## initialize
OSCAR = Model('OSCAR_v3')


##################################################
##   1. CARBON DIOXIDE
##################################################

##===========
## 1.1. Ocean
##===========

## module adapted from:
## (Joos et al., 1996; doi:10.3402/tellusb.v48i3.15921)

## PARAMETER: preindustrial dissolved inorganic carbon in the surface layer
## (Harmann et al. 2011; ISBN:978-0-643-10745-8)
dic_0 = OSCAR.process('dic_0', (), lambda Var, Par: Eq__dic_0(Var, Par), units='umol kgSW-1')
def Eq__dic_0(Var, Par):
    ## CO2Sys Pade approximant
    a0_0 = 30015.6 * (1 - 0.0226536 * (Par.To_0 - 15 - 273.15) + 0.000167105 * (Par.To_0 - 15 - 273.15)**2)
    a1_0 = 13.4574 * (1 - 0.019829 * (Par.To_0 - 15 - 273.15) + 0.000113872 * (Par.To_0 - 15 - 273.15)**2)
    a2_0 = -0.243121 * (1 + 0.000443511 * (Par.To_0 - 15 - 273.15) - 0.000473227 * (Par.To_0 - 15 - 273.15)**2)
    dic_0_Pade = a0_0 * (Par.CO2_0 / 380.) / (1 + a1_0 * (Par.CO2_0 / 380. ) + a2_0 * (Par.CO2_0 / 380.)**2)
    ## CO2Sys Power law
    p0_0 = 2160.156 * (1 - 0.00347063 * (Par.To_0 - 15 - 273.15) - 0.0000250016 * (Par.To_0 - 15 - 273.15)**2)
    p1_0 = 0.0595961 * (1 + 0.0200328 * (Par.To_0 - 15 - 273.15) + 0.000192084 * (Par.To_0 - 15 - 273.15)**2)
    p2_0 = 0.318665 * (1 - 0.00151292 * (Par.To_0 - 15 - 273.15) - 0.000198978 * (Par.To_0 - 15 - 273.15)**2)
    dic_0_Power = p0_0 * ((Par.CO2_0 / 380.) - p2_0)**p1_0    
    ## choosing configuration
    return Par.pCO2_is_Pade * dic_0_Pade + (1-Par.pCO2_is_Pade) * dic_0_Power


## partial pressure of CO2 at sea surface
## (Harmann et al. 2011; ISBN:978-0-643-10745-8)
D_pCO2 = OSCAR.process('D_pCO2', ('dic_0', 'D_dic', 'D_To'), lambda Var, Par: Eq__D_pCO2(Var, Par), units='ppm')
def Eq__D_pCO2(Var, Par):
    ## common variables
    To = Par.To_0 + Var.D_To
    dic = Var.dic_0 + Var.D_dic
    ## following CO2Sys Pade approximant
    a0 = 30015.6 * (1 - 0.0226536 * (To - 15 - 273.15) + 0.000167105 * (To - 15 - 273.15)**2)
    a1 = 13.4574 * (1 - 0.019829 * (To - 15 - 273.15) + 0.000113872 * (To - 15 - 273.15)**2)
    a2 = -0.243121 * (1 + 0.000443511 * (To - 15 - 273.15) - 0.000473227 * (To - 15 - 273.15)**2)
    D_pCO2_Pade = 380 * (a0 - a1*dic - np.sqrt((a0 - a1*dic)**2 - 4*a2*dic**2)) / (2*a2*dic) - Par.CO2_0
    ## CO2Sys Power law
    p0 = 2160.156 * (1 - 0.00347063 * (To - 15 - 273.15) - 0.0000250016 * (To - 15 - 273.15)**2)
    p1 = 0.0595961 * (1 + 0.0200328 * (To - 15 - 273.15) + 0.000192084 * (To - 15 - 273.15)**2)
    p2 = 0.318665 * (1 - 0.00151292 * (To - 15 - 273.15) - 0.000198978 * (To - 15 - 273.15)**2)
    D_pCO2_Power = 380 * (p2 + (dic/p0)**(1/p1)) - Par.CO2_0
    ## choosing configuration
    return Par.pCO2_is_Pade * D_pCO2_Pade + (1-Par.pCO2_is_Pade) * D_pCO2_Power


## mixed-layer depth
D_mld = OSCAR.process('D_mld', ('D_To',), lambda Var, Par: Eq__D_mld(Var, Par), units='m') 
def Eq__D_mld(Var, Par):
    return Par.mld_0 * Par.p_mld * np.expm1(Par.g_mld * Var.D_To)


## dissolved inorganic carbon in the surface layer
D_dic = OSCAR.process('D_dic', ('D_Cosurf', 'D_mld'), lambda Var, Par: Eq__D_dic(Var, Par), units='umol kgSW-1')
def Eq__D_dic(Var, Par):
    return Par.a_dic / Par.A_ocean / Par.a_CO2 / (Par.mld_0 + Var.D_mld) * Var.D_Cosurf.sum('box_osurf', min_count=1)


## ingoing flux to the surface ocean
D_Fin = OSCAR.process('D_Fin', ('D_CO2',), lambda Var, Par: Eq__D_Fin(Var, Par), units='PgC yr-1')
def Eq__D_Fin(Var, Par):
    return Par.p_circ * Par.v_fg * Par.a_CO2 * Var.D_CO2


## outgoing flux from the surface ocean
D_Fout = OSCAR.process('D_Fout', ('D_pCO2',), lambda Var, Par: Eq__D_Fout(Var, Par), units='PgC yr-1')
def Eq__D_Fout(Var, Par):
    return Par.p_circ * Par.v_fg * Par.a_CO2 * Var.D_pCO2


## partial pressure of CO2 at sea surface
D_Fcirc = OSCAR.process('D_Fcirc', ('D_Cosurf',), lambda Var, Par: Eq__D_Fcirc(Var, Par), units='PgC yr-1')
def Eq__D_Fcirc(Var, Par):
    return 1/Par.t_circ * Var.D_Cosurf


## ocean carbon sink
D_Focean = OSCAR.process('D_Focean', ('D_Fin', 'D_Fout'), lambda Var, Par: Eq__D_Focean(Var, Par), units='PgC yr-1')
def Eq__D_Focean(Var, Par):
    return Var.D_Fin.sum('box_osurf', min_count=1) - Var.D_Fout.sum('box_osurf', min_count=1)


## PROGNOSTIC: surface ocean carbon stock
D_Cosurf = OSCAR.process('D_Cosurf', ('D_Cosurf', 'D_Fcirc', 'D_Fout', 'D_Fin'), lambda Var, Par, **Int_args: DiffEq__D_Cosurf(Var, Par, **Int_args), units='PgC', core_dims=['box_osurf'])
def DiffEq__D_Cosurf(Var, Par, Int=Int_dflt, dt=1.):
    v = 1 / Par.t_circ
    return Int(Var.D_Cosurf, v, Var.D_Fin - Var.D_Fout - Var.D_Fcirc, dt)


##=====================
## 1.2. Land: intensive
##=====================

## flexible 3-box model
## {no reference yet}

## PARAMETER: preindustrial vegetation carbon density
cveg_0 = OSCAR.process('cveg_0', (), lambda Var, Par: Eq__cveg_0(Var, Par), units='PgC Mha-1')
def Eq__cveg_0(Var, Par):
    return Par.npp_0 / (Par.mu1_0 + Par.mu2_0 + Par.igni_0 + Par.harv_0 + Par.graz_0)


## PARAMETER: preindustrial litter carbon density
csoil1_0 = OSCAR.process('csoil1_0', (), lambda Var, Par: Eq__csoil1_0(Var, Par), units='PgC Mha-1')
def Eq__csoil1_0(Var, Par):
    cveg_0 = Par.npp_0 / (Par.mu1_0 + Par.mu2_0 + Par.igni_0 + Par.harv_0 + Par.graz_0)
    return ((Par.mu1_0 * cveg_0) / (Par.rho1_0 + Par.muM_0)).where(Par.rho1_0 + Par.muM_0 != 0, 0)


## PARAMETER: preindustrial soil carbon density
csoil2_0 = OSCAR.process('csoil2_0', (), lambda Var, Par: Eq__csoil2_0(Var, Par), units='PgC Mha-1')
def Eq__csoil2_0(Var, Par):
    cveg_0 = Par.npp_0 / (Par.mu1_0 + Par.mu2_0 + Par.igni_0 + Par.harv_0 + Par.graz_0)
    csoil1_0 = ((Par.mu1_0 * cveg_0) / (Par.rho1_0 + Par.muM_0)).where(Par.rho1_0 + Par.muM_0 != 0, 0)
    return (Par.mu2_0 * cveg_0 + Par.muM_0 * csoil1_0) / Par.rho2_0


## fertilisation effect
## (Friedlingstein et al., 1995; doi:10.1029/95GB02381)
f_fert = OSCAR.process('f_fert', ('D_CO2',), lambda Var, Par: Eq__f_fert(Var, Par), units='1')
def Eq__f_fert(Var, Par):
    ## logarithmic formulation
    f_fert_Log = 1 + Par.b_npp * np.log1p(Var.D_CO2 / Par.CO2_0)
    ## hyperbolic formulation
    f_fert_Hyp = (1 + Var.D_CO2 / (Par.CO2_0 - Par.CO2_cp)) / (Var.D_CO2 / Par.CO2_0 * (1/Par.b2_npp * (2*Par.CO2_0 - Par.CO2_cp) / (Par.CO2_0 - Par.CO2_cp) - 1) + 1)
    ## choosing configuration
    return Par.fert_is_Log * f_fert_Log + (1-Par.fert_is_Log) * f_fert_Hyp


## net primary productivity (areal)
D_npp = OSCAR.process('D_npp', ('f_fert', 'D_Tl', 'D_Pl'), lambda Var, Par: Eq__D_npp(Var, Par), units='PgC yr-1 Mha-1')
def Eq__D_npp(Var, Par):
    return Par.npp_0 * (Var.f_fert * (1 + Par.g_nppT * Var.D_Tl + Par.g_nppP * Var.D_Pl) - 1)


## wildfire intensity
f_igni = OSCAR.process('f_igni', ('D_CO2', 'D_Tl', 'D_Pl'), lambda Var, Par: Eq__f_igni(Var, Par), units='1')
def Eq__f_igni(Var, Par):
    return 1 + Par.g_igniC * Var.D_CO2 + Par.g_igniT * Var.D_Tl + Par.g_igniP * Var.D_Pl


## wildfire emissions (areal)
D_efire = OSCAR.process('D_efire', ('f_igni', 'cveg_0', 'D_cveg'), lambda Var, Par: Eq__D_efire(Var, Par), units='PgC yr-1 Mha-1')
def Eq__D_efire(Var, Par):
    return Par.igni_0 * ((Var.cveg_0 + Var.D_cveg) * Var.f_igni - Var.cveg_0)


## crop harvest emissions (areal)
D_eharv = OSCAR.process('D_eharv', ('D_cveg',), lambda Var, Par: Eq__D_eharv(Var, Par), units='PgC yr-1 Mha-1')
def Eq__D_eharv(Var, Par):
    return Par.harv_0 * Var.D_cveg


## pasture grazing emissions (areal)
D_egraz = OSCAR.process('D_egraz', ('D_cveg',), lambda Var, Par: Eq__D_egraz(Var, Par), units='PgC yr-1 Mha-1')
def Eq__D_egraz(Var, Par):
    return Par.graz_0 * Var.D_cveg


## mortality flux to litter (areal)
D_fmort1 = OSCAR.process('D_fmort1', ('D_cveg',), lambda Var, Par: Eq__D_fmort1(Var, Par), units='PgC yr-1 Mha-1')
def Eq__D_fmort1(Var, Par):
    return Par.mu1_0 * Var.D_cveg


## mortality flux to soil (areal)
D_fmort2 = OSCAR.process('D_fmort2', ('D_cveg',), lambda Var, Par: Eq__D_fmort2(Var, Par), units='PgC yr-1 Mha-1')
def Eq__D_fmort2(Var, Par):
    return Par.mu2_0 * Var.D_cveg


## heterotrophic respiration factor
## (Tuomi et al., 2008; doi:10.1016/j.ecolmodel.2007.09.003)
f_resp = OSCAR.process('f_resp', ('D_Tl', 'D_Pl'), lambda Var, Par: Eq__f_resp(Var, Par), units='1')
def Eq__f_resp(Var, Par):
    return np.exp(Par.g_respT * Var.D_Tl + Par.g_respT2 * Var.D_Tl**2) * (1 + Par.g_respP * Var.D_Pl)


## litter respiration flux (areal)
D_rh1 = OSCAR.process('D_rh1', ('f_resp', 'csoil1_0', 'D_csoil1'), lambda Var, Par: Eq__D_rh1(Var, Par), units='PgC yr-1 Mha-1')
def Eq__D_rh1(Var, Par):
    return Par.rho1_0 * ((Var.csoil1_0 + Var.D_csoil1) * Var.f_resp - Var.csoil1_0)


## metabolisation flux (areal)
D_fmet = OSCAR.process('D_fmet', ('f_resp', 'csoil1_0', 'D_csoil1'), lambda Var, Par: Eq__D_fmet(Var, Par), units='PgC yr-1 Mha-1')
def Eq__D_fmet(Var, Par):
    return Par.muM_0 * ((Var.csoil1_0 + Var.D_csoil1) * Var.f_resp - Var.csoil1_0)


## soil respiration flux (areal)
D_rh2 = OSCAR.process('D_rh2', ('f_resp', 'csoil2_0', 'D_csoil2'), lambda Var, Par: Eq__D_rh2(Var, Par), units='PgC yr-1 Mha-1')
def Eq__D_rh2(Var, Par):
    return Par.rho2_0 * ((Var.csoil2_0 + Var.D_csoil2) * Var.f_resp - Var.csoil2_0)


## net biome productivity (areal)
D_nbp = OSCAR.process('D_nbp', ('D_npp', 'D_efire', 'D_eharv', 'D_egraz', 'D_rh1', 'D_rh2'), lambda Var, Par: Eq__D_nbp(Var, Par), units='PgC yr-1 Mha-1')
def Eq__D_nbp(Var, Par):
    return Var.D_npp - Var.D_efire - Var.D_eharv - Var.D_egraz - Var.D_rh1 - Var.D_rh2


## PROGNOSTIC: vegetation carbon stock (areal)
D_cveg = OSCAR.process('D_cveg', ('D_cveg', 'D_npp', 'D_efire', 'D_eharv', 'D_egraz', 'D_fmort1', 'D_fmort2'), lambda Var, Par, **Int_args: DiffEq__D_cveg(Var, Par, **Int_args), units='PgC Mha-1', core_dims=['reg_land', 'bio_land'])
def DiffEq__D_cveg(Var, Par, Int=Int_dflt, dt=1.):
    v = Par.igni_0 + Par.harv_0 + Par.graz_0 + Par.mu1_0 + Par.mu2_0
    return Int(Var.D_cveg, v, Var.D_npp - Var.D_efire - Var.D_eharv - Var.D_egraz - Var.D_fmort1 - Var.D_fmort2, dt)


## PROGNOSTIC: litter carbon stock (areal)
D_csoil1 = OSCAR.process('D_csoil1', ('D_csoil1', 'D_fmort1', 'D_rh1', 'D_fmet'), lambda Var, Par, **Int_args: DiffEq__D_csoil1(Var, Par, **Int_args), units='PgC Mha-1', core_dims=['reg_land', 'bio_land'])
def DiffEq__D_csoil1(Var, Par, Int=Int_dflt, dt=1.):
    v = (Par.rho1_0 + Par.muM_0).where(Par.rho1_0 + Par.muM_0 != 0, 1E-18)
    return Int(Var.D_csoil1, v, Var.D_fmort1 - Var.D_rh1 - Var.D_fmet, dt)


## PROGNOSTIC: soil carbon stock (areal)
D_csoil2 = OSCAR.process('D_csoil2', ('D_csoil2', 'D_fmort2', 'D_fmet', 'D_rh2'), lambda Var, Par, **Int_args: DiffEq__D_csoil2(Var, Par, **Int_args), units='PgC Mha-1', core_dims=['reg_land', 'bio_land'])
def DiffEq__D_csoil2(Var, Par, Int=Int_dflt, dt=1.):
    v = Par.rho2_0
    return Int(Var.D_csoil2, v, Var.D_fmort2 + Var.D_fmet - Var.D_rh2, dt)


##=====================
## 1.3. Land: extensive
##=====================

## book-keeping based on:
## (Gitz & Ciais, 2003; doi:10.1029/2002GB001963)
## and reformulated following definition 3 of:
## (Gasser & Ciais, 2013; doi:10.5194/esd-4-171-2013)

## land-use book-keeping initialization of vegetation carbon
D_Fveg_bk = OSCAR.process('D_Fveg_bk', ('cveg_0', 'D_cveg', 'd_Acover', 'd_Hwood', 'd_Ashift'), lambda Var, Par: Eq__D_Fveg_bk(Var, Par), units='PgC yr-1')
def Eq__D_Fveg_bk(Var, Par):
    D_Fveg_bk_lcc = -(Var.cveg_0 + Var.D_cveg).rename({'bio_land':'bio_to'}) * Var.d_Acover
    D_Fveg_bk_harv = -0.5*(Var.d_Hwood.rename({'bio_land':'bio_from'}) + Var.d_Hwood.rename({'bio_land':'bio_to'})).where(Var.d_Acover.bio_from == Var.d_Acover.bio_to, 0.)
    D_Fveg_bk_shift = -((Var.cveg_0 + Var.D_cveg) * (1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift))).rename({'bio_land':'bio_to'}) * Var.d_Ashift
    return D_Fveg_bk_lcc + D_Fveg_bk_harv + D_Fveg_bk_shift


## land-use book-keeping initialization of litter carbon
D_Fsoil1_bk = OSCAR.process('D_Fsoil1_bk', ('csoil1_0', 'D_csoil1', 'd_Acover'), lambda Var, Par: Eq__D_Fsoil1_bk(Var, Par), units='PgC yr-1')
def Eq__D_Fsoil1_bk(Var, Par):
    D_Fsoil1_bk_lcc = ((Var.csoil1_0 + Var.D_csoil1).rename({'bio_land':'bio_from'}) - (Var.csoil1_0 + Var.D_csoil1).rename({'bio_land':'bio_to'})) * Var.d_Acover
    return D_Fsoil1_bk_lcc


## land-use book-keeping initialization of soil carbon
D_Fsoil2_bk = OSCAR.process('D_Fsoil2_bk', ('csoil2_0', 'D_csoil2', 'd_Acover'), lambda Var, Par: Eq__D_Fsoil2_bk(Var, Par), units='PgC yr-1')
def Eq__D_Fsoil2_bk(Var, Par):
    D_Fsoil2_bk_lcc = ((Var.csoil2_0 + Var.D_csoil2).rename({'bio_land':'bio_from'}) - (Var.csoil2_0 + Var.D_csoil2).rename({'bio_land':'bio_to'})) * Var.d_Acover
    return D_Fsoil2_bk_lcc


## land-use initial slash flux
D_Fslash = OSCAR.process('D_Fslash', ('cveg_0', 'D_cveg', 'd_Acover', 'd_Hwood', 'd_Ashift'), lambda Var, Par: Eq__D_Fslash(Var, Par), units='PgC yr-1')
def Eq__D_Fslash(Var, Par):
    D_Fslash_lcc = ((Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb))).rename({'bio_land':'bio_from'}) * Var.d_Acover
    D_Fslash_harv = (1 - Par.p_hwp.sum('box_hwp', min_count=1)).rename({'bio_land':'bio_from'}) * 0.5*(Var.d_Hwood.rename({'bio_land':'bio_from'}) + Var.d_Hwood.rename({'bio_land':'bio_to'})).where(Var.d_Acover.bio_from == Var.d_Acover.bio_to, 0.)
    D_Fslash_shift = ((Var.cveg_0 + Var.D_cveg) * (Par.p_agb * (1 - Par.p_hwp.sum('box_hwp', min_count=1)) + (1 - Par.p_agb)) * (1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift))).rename({'bio_land':'bio_from'}) * Var.d_Ashift
    return D_Fslash_lcc + D_Fslash_harv + D_Fslash_shift


## land-use slash flux to litter pool
D_Fslash1 = OSCAR.process('D_Fslash1', ('D_Fslash',), lambda Var, Par: Eq__D_Fslash1(Var, Par), units='PgC yr-1')
def Eq__D_Fslash1(Var, Par):
    return (Par.mu1_0 / (Par.mu1_0 + Par.mu2_0)).rename({'bio_land':'bio_from'}) * Var.D_Fslash


## land-use slash flux to soil pool
D_Fslash2 = OSCAR.process('D_Fslash2', ('D_Fslash',), lambda Var, Par: Eq__D_Fslash2(Var, Par), units='PgC yr-1')
def Eq__D_Fslash2(Var, Par):
    return (Par.mu2_0 / (Par.mu1_0 + Par.mu2_0)).rename({'bio_land':'bio_from'}) * Var.D_Fslash


## land-use initial flux to harvested wood product pools
D_Fhwp = OSCAR.process('D_Fhwp', ('cveg_0', 'D_cveg', 'd_Acover', 'd_Hwood', 'd_Ashift'), lambda Var, Par: Eq__D_Fhwp(Var, Par), units='PgC yr-1')
def Eq__D_Fhwp(Var, Par):
    D_Fhwp_lcc = ((Var.cveg_0 + Var.D_cveg) * Par.p_agb * Par.p_hwp).rename({'bio_land':'bio_from'}) * Var.d_Acover
    D_Fhwp_harv = Par.p_hwp.rename({'bio_land':'bio_from'}) * 0.5*(Var.d_Hwood.rename({'bio_land':'bio_from'}) + Var.d_Hwood.rename({'bio_land':'bio_to'})).where(Var.d_Acover.bio_from == Var.d_Acover.bio_to, 0.)
    D_Fhwp_shift = ((Var.cveg_0 + Var.D_cveg) * Par.p_agb * Par.p_hwp * (1 - np.exp(-Par.npp_0 / Var.cveg_0 * Par.t_shift))).rename({'bio_land':'bio_from'}) * Var.d_Ashift
    return D_Fhwp_lcc + D_Fhwp_harv + D_Fhwp_shift


## altered NPP (under book-keeping)
D_NPP_bk = OSCAR.process('D_NPP_bk', (), lambda Var, Par: Eq__D_NPP_bk(Var, Par), units='PgC yr-1')
def Eq__D_NPP_bk(Var, Par):
    return 0.


## wildfire emissions (under book-keeping)
D_Efire_bk = OSCAR.process('D_Efire_bk', ('f_igni', 'D_Cveg_bk'), lambda Var, Par: Eq__D_Efire_bk(Var, Par), units='PgC yr-1')
def Eq__D_Efire_bk(Var, Par):
    return (Par.igni_0 * Var.f_igni).rename({'bio_land':'bio_to'}) * Var.D_Cveg_bk


## crop harvest emissions (under book-keeping)
D_Eharv_bk = OSCAR.process('D_Eharv_bk', ('D_Cveg_bk',), lambda Var, Par: Eq__D_Eharv_bk(Var, Par), units='PgC yr-1')
def Eq__D_Eharv_bk(Var, Par):
    return Par.harv_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_bk


## pasture grazing emissions (under book-keeping)
D_Egraz_bk = OSCAR.process('D_Egraz_bk', ('D_Cveg_bk',), lambda Var, Par: Eq__D_Egraz_bk(Var, Par), units='PgC yr-1')
def Eq__D_Egraz_bk(Var, Par):
    return Par.graz_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_bk


## mortality flux to litter (under book-keeping)
D_Fmort1_bk = OSCAR.process('D_Fmort1_bk', ('D_Cveg_bk',), lambda Var, Par: Eq__D_Fmort1_bk(Var, Par), units='PgC yr-1')
def Eq__D_Fmort1_bk(Var, Par):
    return Par.mu1_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_bk


## mortality flux to soil (under book-keeping)
D_Fmort2_bk = OSCAR.process('D_Fmort2_bk', ('D_Cveg_bk',), lambda Var, Par: Eq__D_Fmort2_bk(Var, Par), units='PgC yr-1')
def Eq__D_Fmort2_bk(Var, Par):
    return Par.mu2_0.rename({'bio_land':'bio_to'}) * Var.D_Cveg_bk


## litter respiration flux (under book-keeping)
D_Rh1_bk = OSCAR.process('D_Rh1_bk', ('f_resp', 'D_Csoil1_bk'), lambda Var, Par: Eq__D_Rh1_bk(Var, Par), units='PgC yr-1')
def Eq__D_Rh1_bk(Var, Par):
    return (Par.rho1_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil1_bk


## metabolisation flux (under book-keeping)
D_Fmet_bk = OSCAR.process('D_Fmet_bk', ('f_resp', 'D_Csoil1_bk'), lambda Var, Par: Eq__D_Fmet_bk(Var, Par), units='PgC yr-1')
def Eq__D_Fmet_bk(Var, Par):
    return (Par.muM_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil1_bk


## soil respiration flux (under book-keeping)
D_Rh2_bk = OSCAR.process('D_Rh2_bk', ('f_resp', 'D_Csoil2_bk'), lambda Var, Par: Eq__D_Rh2_bk(Var, Par), units='PgC yr-1')
def Eq__D_Rh2_bk(Var, Par):
    return (Par.rho2_0 * Var.f_resp).rename({'bio_land':'bio_to'}) * Var.D_Csoil2_bk


## harvested wood product oxidation
D_Ehwp = OSCAR.process('D_Ehwp', ('D_Chwp',), lambda Var, Par: Eq__D_Ehwp(Var, Par), units='PgC yr-1')
def Eq__D_Ehwp(Var, Par):
    return 1 / Par.w_t_hwp / Par.t_hwp * Var.D_Chwp


## net biome productivity (under book-keeping)
D_NBP_bk = OSCAR.process('D_NBP_bk', ('D_NPP_bk', 'D_Efire_bk', 'D_Eharv_bk', 'D_Egraz_bk', 'D_Rh1_bk', 'D_Rh2_bk', 'D_Ehwp'), lambda Var, Par: Eq__D_NBP_bk(Var, Par), units='PgC yr-1')
def Eq__D_NBP_bk(Var, Par):
    return Var.D_NPP_bk - Var.D_Efire_bk - Var.D_Eharv_bk - Var.D_Egraz_bk - Var.D_Rh1_bk - Var.D_Rh2_bk - Var.D_Ehwp.sum('box_hwp', min_count=1)


## land-use change emissions
D_Eluc = OSCAR.process('D_Eluc', ('D_NBP_bk',), lambda Var, Par: Eq__D_Eluc(Var, Par), units='PgC yr-1')
def Eq__D_Eluc(Var, Par):
    return -Var.D_NBP_bk.sum('bio_from', min_count=1).sum('bio_to', min_count=1).sum('reg_land', min_count=1)


## land carbon sink
D_Fland = OSCAR.process('D_Fland', ('D_nbp', 'D_Aland'), lambda Var, Par: Eq__D_Fland(Var, Par), units='PgC yr-1')
def Eq__D_Fland(Var, Par):
    return (Var.D_nbp * (Par.Aland_0 + Var.D_Aland)).sum('bio_land', min_count=1).sum('reg_land', min_count=1)


## PROGNOSTIC: biome area
D_Aland = OSCAR.process('D_Aland', ('D_Aland', 'd_Acover'), lambda Var, Par, **Int_args: DiffEq__D_Aland(Var, Par, **Int_args), units='Mha', core_dims=['reg_land', 'bio_land'])
def DiffEq__D_Aland(Var, Par, Int=Int_dflt, dt=1.):
    return Int(Var.D_Aland, 1E-18, Var.d_Acover.sum('bio_from', min_count=1).rename({'bio_to':'bio_land'}) - Var.d_Acover.sum('bio_to', min_count=1).rename({'bio_from':'bio_land'}), dt)


## PROGNOSTIC: vegetation carbon stock (under book-keeping)
D_Cveg_bk = OSCAR.process('D_Cveg_bk', ('D_Cveg_bk', 'D_Fveg_bk', 'D_NPP_bk', 'D_Efire_bk', 'D_Eharv_bk', 'D_Egraz_bk', 'D_Fmort1_bk', 'D_Fmort2_bk'), lambda Var, Par, **Int_args: DiffEq__D_Cveg_bk(Var, Par, **Int_args), units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])
def DiffEq__D_Cveg_bk(Var, Par, Int=Int_dflt, dt=1.):
    v = (Par.igni_0 + Par.harv_0 + Par.graz_0 + Par.mu1_0 + Par.mu2_0).rename({'bio_land':'bio_to'})
    return Int(Var.D_Cveg_bk, v, Var.D_Fveg_bk + Var.D_NPP_bk - Var.D_Efire_bk - Var.D_Eharv_bk - Var.D_Egraz_bk - Var.D_Fmort1_bk - Var.D_Fmort2_bk, dt)


## PROGNOSTIC: litter carbon stock (under book-keeping)
D_Csoil1_bk = OSCAR.process('D_Csoil1_bk', ('D_Csoil1_bk', 'D_Fsoil1_bk', 'D_Fslash1', 'D_Fmort1_bk', 'D_Rh1_bk', 'D_Fmet_bk'), lambda Var, Par, **Int_args: DiffEq__D_Csoil1_bk(Var, Par, **Int_args), units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])
def DiffEq__D_Csoil1_bk(Var, Par, Int=Int_dflt, dt=1.):
    v = (Par.rho1_0 + Par.muM_0).where(Par.rho1_0 + Par.muM_0 != 0, 1E-18).rename({'bio_land':'bio_to'})
    return Int(Var.D_Csoil1_bk, v, Var.D_Fsoil1_bk + Var.D_Fslash1 + Var.D_Fmort1_bk - Var.D_Rh1_bk - Var.D_Fmet_bk, dt)


## PROGNOSTIC: soil carbon stock (under book-keeping)
D_Csoil2_bk = OSCAR.process('D_Csoil2_bk', ('D_Csoil2_bk', 'D_Fsoil2_bk', 'D_Fslash2', 'D_Fmort2_bk', 'D_Fmet_bk', 'D_Rh2_bk'), lambda Var, Par, **Int_args: DiffEq__D_Csoil2_bk(Var, Par, **Int_args), units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to'])
def DiffEq__D_Csoil2_bk(Var, Par, Int=Int_dflt, dt=1.):
    v = Par.rho2_0.rename({'bio_land':'bio_to'})
    return Int(Var.D_Csoil2_bk, v, Var.D_Fsoil2_bk + Var.D_Fslash2 + Var.D_Fmort2_bk + Var.D_Fmet_bk - Var.D_Rh2_bk, dt)


## PROGNOSTIC: harvested wood products stock
D_Chwp = OSCAR.process('D_Chwp', ('D_Chwp', 'D_Fhwp', 'D_Ehwp'), lambda Var, Par, **Int_args: DiffEq__D_Chwp(Var, Par, **Int_args), units='PgC', core_dims=['reg_land', 'bio_from', 'bio_to', 'box_hwp'])
def DiffEq__D_Chwp(Var, Par, Int=Int_dflt, dt=1.):
    v = 1 / Par.w_t_hwp / Par.t_hwp
    return Int(Var.D_Chwp, v, Var.D_Fhwp - Var.D_Ehwp, dt)


##================
## 1.4. Permafrost
##================

## module taken from:
## (Gasser et al., 2018; doi:10.1038/s41561-018-0227-0)

## heterotrophic respiration factor for permafrost
f_resp_pf = OSCAR.process('f_resp_pf', ('D_Tg',), lambda Var, Par: Eq__f_resp_pf(Var, Par), units='1')
def Eq__f_resp_pf(Var, Par):
    return np.exp(Par.k_resp_pf * (Par.g_respT_pf * Par.w_clim_pf * Var.D_Tg + Par.g_respT2_pf * (Par.w_clim_pf * Var.D_Tg)**2))


## theoretical thawed fraction
D_pthaw_bar = OSCAR.process('D_pthaw_bar', ('D_Tg',), lambda Var, Par: Eq__D_pthaw_bar(Var, Par), units='1')
def Eq__D_pthaw_bar(Var, Par):
    return -Par.pthaw_min + (1 + Par.pthaw_min)  /(1 + ((1/Par.pthaw_min + 1)**Par.k_pthaw - 1)*np.exp(-Par.g_pthaw * Par.k_pthaw * Par.w_clim_pf * Var.D_Tg))**(1/Par.k_pthaw)


## derivative of actual thawed fraction
d_pthaw = OSCAR.process('d_pthaw', ('D_pthaw_bar', 'D_pthaw'), lambda Var, Par: Eq__d_pthaw(Var, Par), units='yr-1')
def Eq__d_pthaw(Var, Par):
    # other way to formulate (because v_thaw > v_froz):
    # 0.5 * (Par.v_thaw + Par.v_froz) * (Var.D_pthaw_bar - Var.D_pthaw) + 0.5 * np.abs((Par.v_thaw - Par.v_froz) * (Var.D_pthaw_bar - Var.D_pthaw))
    return (Par.v_thaw * (Var.D_pthaw_bar >= Var.D_pthaw) + Par.v_froz * (Var.D_pthaw_bar < Var.D_pthaw)) * (Var.D_pthaw_bar - Var.D_pthaw)


## PROGNOSTIC: actual thawed fraction
D_pthaw = OSCAR.process('D_pthaw', ('D_pthaw', 'd_pthaw'), lambda Var, Par, **Int_args: DiffEq__D_pthaw(Var, Par, **Int_args), units='1', core_dims=['reg_pf'])
def DiffEq__D_pthaw(Var, Par, Int=Int_dflt, dt=1.):
    v = 0.5 * (Par.v_thaw + Par.v_froz)
    return Int(Var.D_pthaw, v, Var.d_pthaw, dt)


## flux of thawing permafrost carbon
D_Fthaw = OSCAR.process('D_Fthaw', ('d_pthaw',), lambda Var, Par: Eq__D_Fthaw(Var, Par), units='PgC yr-1')
def Eq__D_Fthaw(Var, Par):
    return Par.Cfroz_0 * Var.d_pthaw


## emissions from thawed permafrost carbon
D_Ethaw = OSCAR.process('D_Ethaw', ('D_Cthaw', 'f_resp_pf'), lambda Var, Par: Eq__D_Ethaw(Var, Par), units='PgC yr-1')
def Eq__D_Ethaw(Var, Par):
    return 1/Par.t_pf_thaw * Var.f_resp_pf * Var.D_Cthaw


## total permafrost carbon emissions
D_Epf = OSCAR.process('D_Epf', ('D_Ethaw', 'D_Fthaw'), lambda Var, Par: Eq__D_Epf(Var, Par), units='PgC yr-1')
def Eq__D_Epf(Var, Par):
    return Var.D_Ethaw.sum('box_thaw', min_count=1) + Par.p_pf_inst * Var.D_Fthaw


## CO2 permafrost emissions
D_Epf_CO2 = OSCAR.process('D_Epf_CO2', ('D_Epf',), lambda Var, Par: Eq__D_Epf_CO2(Var, Par), units='PgC yr-1')
def Eq__D_Epf_CO2(Var, Par):
    return (1 - Par.p_pf_CH4) * Var.D_Epf


## CH4 permafrost emissions
D_Epf_CH4 = OSCAR.process('D_Epf_CH4', ('D_Epf',), lambda Var, Par: Eq__D_Epf_CH4(Var, Par), units='TgC yr-1')
def Eq__D_Epf_CH4(Var, Par):
    a_conv = 1E3 # from {PgC} to {TgC}
    return a_conv * Par.p_pf_CH4 * Var.D_Epf


## PROGNOSTIC: frozen permafrost carbon
D_Cfroz = OSCAR.process('D_Cfroz', ('D_Cfroz', 'D_Fthaw'), lambda Var, Par, **Int_args: DiffEq__D_Cfroz(Var, Par, **Int_args), units='PgC', core_dims=['reg_pf'])
def DiffEq__D_Cfroz(Var, Par, Int=Int_dflt, dt=1.):
    return Int(Var.D_Cfroz, 1E-18, -Var.D_Fthaw, dt)


## PROGNOSTIC: thawed permafrost carbon
D_Cthaw = OSCAR.process('D_Cthaw', ('D_Cthaw', 'D_Fthaw', 'D_Ethaw'), lambda Var, Par, **Int_args: DiffEq__D_Cthaw(Var, Par, **Int_args), units='PgC', core_dims=['reg_pf', 'box_thaw'])
def DiffEq__D_Cthaw(Var, Par, Int=Int_dflt, dt=1.):
    v = 1 / Par.t_pf_thaw
    return Int(Var.D_Cthaw, v, Par.p_pf_thaw * (1 - Par.p_pf_inst) * Var.D_Fthaw - Var.D_Ethaw, dt)


##================
## 1.5. Atmosphere
##================

## PROGNOSTIC: atmospheric CO2
D_CO2 = OSCAR.process('D_CO2', ('D_CO2', 'Eff', 'D_Eluc', 'D_Epf_CO2', 'D_Fland', 'D_Focean', 'D_Foxi_CH4'), lambda Var, Par, **Int_args: DiffEq__D_CO2(Var, Par, **Int_args), units='ppm')
def DiffEq__D_CO2(Var, Par, Int=Int_dflt, dt=1.):
    return Int(Var.D_CO2, 1E-18, 1/Par.a_CO2 * (Var.Eff.sum('reg_land', min_count=1) + Var.D_Eluc + Var.D_Epf_CO2.sum('reg_pf', min_count=1) - Var.D_Fland - Var.D_Focean + Var.D_Foxi_CH4), dt)


## METRIC: airborne fraction
AF = OSCAR.process('AF', ('Eff', 'D_Eluc', 'D_Epf_CO2', 'D_Fland', 'D_Focean', 'D_Foxi_CH4'), lambda Var, Par: Eq__AF(Var, Par), units='1')
def Eq__AF(Var, Par):
    return 1 + (Var.D_Epf_CO2.sum('reg_pf', min_count=1) - Var.D_Fland - Var.D_Focean + Var.D_Foxi_CH4) / (Var.Eff.sum('reg_land', min_count=1) + Var.D_Eluc)


## METRIC: carbon sinks rate
kS = OSCAR.process('kS', ('D_Epf_CO2', 'D_Fland', 'D_Focean', 'D_Foxi_CH4', 'D_CO2'), lambda Var, Par: Eq__kS(Var, Par), units='yr-1')
def Eq__kS(Var, Par):
    return 1/Par.a_CO2 * -(Var.D_Epf_CO2.sum('reg_pf', min_count=1) - Var.D_Fland - Var.D_Focean + Var.D_Foxi_CH4) / Var.D_CO2


## CO2 radiative forcing
## (Myhre et al., 1998; doi:10.1029/98GL01908)
RF_CO2 = OSCAR.process('RF_CO2', ('D_CO2',), lambda Var, Par: Eq__RF_CO2(Var, Par), units='W m-2')
def Eq__RF_CO2(Var, Par):
    return Par.rf_CO2 * np.log1p(Var.D_CO2 / Par.CO2_0)


##################################################
##   2. NON-CO2 SPECIES
##################################################

##=====================
## 2.1. Biomass burning
##=====================

## CO2 emissions from natural biomass burning (= wildfire)
D_Efire = OSCAR.process('D_Efire', ('cveg_0', 'D_efire', 'D_Aland', 'D_Efire_bk'), lambda Var, Par: Eq__D_Efire(Var, Par), units='PgC yr-1')
def Eq__D_Efire(Var, Par):
    return Par.igni_0 * Var.cveg_0 * Var.D_Aland + Var.D_efire * (Par.Aland_0 + Var.D_Aland) + Var.D_Efire_bk.sum('bio_from', min_count=1).rename({'bio_to':'bio_land'})


## non-CO2 emissions from natural biomass burning
D_Ebb_nat = OSCAR.process('D_Ebb_nat', ('D_Efire',), lambda Var, Par: Eq__D_Ebb_nat(Var, Par), units='TgX yr-1')
def Eq__D_Ebb_nat(Var, Par):
    return Par.a_bb * Var.D_Efire


## non-CO2 emissions from anthropogenic biomass burning
D_Ebb_ant = OSCAR.process('D_Ebb_ant', ('D_Ehwp',), lambda Var, Par: Eq__D_Ebb_ant(Var, Par), units='TgX yr-1')
def Eq__D_Ebb_ant(Var, Par):
    return Par.a_bb * (Par.p_hwp_bb * Var.D_Ehwp).sum('box_hwp', min_count=1).sum('bio_to', min_count=1).rename({'bio_from':'bio_land'})


## non-CO2 total emissions from biomass burning
D_Ebb = OSCAR.process('D_Ebb', ('D_Ebb_nat', 'D_Ebb_ant'), lambda Var, Par: Eq__D_Ebb(Var, Par), units='TgX yr-1')
def Eq__D_Ebb(Var, Par):
    return Var.D_Ebb_nat + Var.D_Ebb_ant


##==========================
## 2.2. Lagged concentration
##==========================

## PROGNOSTIC: lagged atmospheric CH4
D_CH4_lag = OSCAR.process('D_CH4_lag', ('D_CH4_lag', 'D_CH4'), lambda Var, Par, **Int_args: DiffEq__D_CH4_lag(Var, Par, **Int_args), units='ppb')
def DiffEq__D_CH4_lag(Var, Par, Int=Int_dflt, dt=1.):
    v = 1 / Par.t_lag
    return Int(Var.D_CH4_lag, v, 1/Par.t_lag * (Var.D_CH4 - Var.D_CH4_lag), dt)


## PROGNOSTIC: lagged atmospheric N2O
D_N2O_lag = OSCAR.process('D_N2O_lag', ('D_N2O_lag', 'D_N2O'), lambda Var, Par, **Int_args: DiffEq__D_N2O_lag(Var, Par, **Int_args), units='ppb')
def DiffEq__D_N2O_lag(Var, Par, Int=Int_dflt, dt=1.):
    v = 1 / Par.t_lag
    return Int(Var.D_N2O_lag, v, 1/Par.t_lag * (Var.D_N2O - Var.D_N2O_lag), dt)


## PROGNOSTIC: lagged atmospheric halogenated compounds
D_Xhalo_lag = OSCAR.process('D_Xhalo_lag', ('D_Xhalo_lag', 'D_Xhalo'), lambda Var, Par, **Int_args: DiffEq__D_Xhalo_lag(Var, Par, **Int_args), units='ppt', core_dims=['spc_halo'])
def DiffEq__D_Xhalo_lag(Var, Par, Int=Int_dflt, dt=1.):
    v = 1 / Par.t_lag
    return Int(Var.D_Xhalo_lag, v, 1/Par.t_lag * (Var.D_Xhalo - Var.D_Xhalo_lag), dt)


##################################################
##   3. METHANE
##################################################

##====================
## 3.1. Atmo chemistry
##====================

## module based on:
## (Holmes et al., 2013; doi:10.5194/acp-13-285-2013)

## atmospheric temperature
D_Ta = OSCAR.process('D_Ta', ('D_Tg',), lambda Var, Par: Eq__D_Ta(Var, Par), units='K')
def Eq__D_Ta(Var, Par):
    return Par.w_clim_Ta * Var.D_Tg


## atmospheric relative humidity
D_f_Qa = OSCAR.process('D_f_Qa', ('D_Ta',), lambda Var, Par: Eq__D_f_Qa(Var, Par), units='1')
def Eq__D_f_Qa(Var, Par):
    return Par.k_Qa * np.expm1(Par.k_svp * Var.D_Ta / (Par.Ta_0 + Par.T_svp))


## hydroxyl sink intensity
f_kOH = OSCAR.process('f_kOH', ('D_CH4', 'D_O3s', 'D_Ta', 'D_f_Qa', 'E_NOX', 'E_CO', 'E_VOC', 'D_Ebb'), lambda Var, Par: Eq__f_kOH(Var, Par), units='1')
def Eq__f_kOH(Var, Par):
    ## environmental factors
    f_kOH_CH4 = Par.x_OH_CH4 * np.log1p(Var.D_CH4 / Par.CH4_0)
    f_kOH_O3s = Par.x_OH_O3s * np.log1p(Var.D_O3s / Par.O3s_0)
    f_kOH_Tg = Par.x_OH_Ta * np.log1p(Var.D_Ta / Par.Ta_0) + Par.x_OH_Qa * np.log1p(Var.D_f_Qa)
    ## logarithmic anthropogenic factors
    f_kOH_NOX_log = Par.x_OH_NOX * np.log1p((Var.E_NOX + Var.D_Ebb.sel({'spc_bb':'NOX'}, drop=True).sum( 'bio_land', min_count=1)).sum('reg_land', min_count=1) / Par.Enat_NOX)
    f_kOH_CO_log = Par.x_OH_CO * np.log1p((Var.E_CO + Var.D_Ebb.sel({'spc_bb':'CO'}, drop=True).sum('bio_land', min_count=1)).sum('reg_land', min_count=1) / Par.Enat_CO)
    f_kOH_VOC_log = Par.x_OH_VOC * np.log1p((Var.E_VOC + Var.D_Ebb.sel({'spc_bb':'VOC'}, drop=True).sum('bio_land', min_count=1)).sum('reg_land', min_count=1) / Par.Enat_VOC)
    # linear anthropogenic factors
    f_kOH_NOX_lin = Par.x2_OH_NOX * (Var.E_NOX + Var.D_Ebb.sel({'spc_bb':'NOX'}, drop=True).sum('bio_land', min_count=1)).sum('reg_land', min_count=1)
    f_kOH_CO_lin = Par.x2_OH_CO * (Var.E_CO + Var.D_Ebb.sel({'spc_bb':'CO'}, drop=True).sum('bio_land', min_count=1)).sum('reg_land', min_count=1)
    f_kOH_VOC_lin = Par.x2_OH_VOC * (Var.E_VOC + Var.D_Ebb.sel({'spc_bb':'VOC'}, drop=True).sum('bio_land', min_count=1)).sum('reg_land', min_count=1)
    ## choosing configuration
    return np.exp(f_kOH_CH4 + f_kOH_O3s + f_kOH_Tg + Par.kOH_is_Log * (f_kOH_NOX_log + f_kOH_CO_log + f_kOH_VOC_log) + (1-Par.kOH_is_Log) * (f_kOH_NOX_lin + f_kOH_CO_lin + f_kOH_VOC_lin))


## CH4 hydroxyl sink
D_Foh_CH4 = OSCAR.process('D_Foh_CH4', ('f_kOH', 'D_CH4'), lambda Var, Par: Eq__D_Foh_CH4(Var, Par), units='TgC yr-1')
def Eq__D_Foh_CH4(Var, Par):
    return Par.a_CH4 / Par.w_t_OH / Par.t_OH_CH4 * ((Par.CH4_0 + Var.D_CH4) * Var.f_kOH - Par.CH4_0)


## CH4 stratospheric sink
D_Fhv_CH4 = OSCAR.process('D_Fhv_CH4', ('f_hv', 'D_CH4_lag'), lambda Var, Par: Eq__D_Fhv_CH4(Var, Par), units='TgC yr-1')
def Eq__D_Fhv_CH4(Var, Par):
    return Par.a_CH4 / Par.w_t_hv / Par.t_hv_CH4 * ((Par.CH4_0 + Var.D_CH4_lag) * Var.f_hv - Par.CH4_0)


## CH4 dry soil sink (= methanotrophy)
D_Fsoil_CH4 = OSCAR.process('D_Fsoil_CH4', ('D_CH4',), lambda Var, Par: Eq__D_Fsoil_CH4(Var, Par), units='TgC yr-1')
def Eq__D_Fsoil_CH4(Var, Par):
    return Par.a_CH4 / Par.t_soil_CH4 * Var.D_CH4


## CH4 oceanic boundary layer sink
D_Focean_CH4 = OSCAR.process('D_Focean_CH4', ('D_CH4',), lambda Var, Par: Eq__D_Focean_CH4(Var, Par), units='TgC yr-1')
def Eq__D_Focean_CH4(Var, Par):
    return Par.a_CH4 / Par.t_ocean_CH4 * Var.D_CH4


## CH4 total atmospheric sink
D_Fsink_CH4 = OSCAR.process('D_Fsink_CH4', ('D_Foh_CH4', 'D_Fhv_CH4', 'D_Fsoil_CH4', 'D_Focean_CH4'), lambda Var, Par: Eq__D_Fsink_CH4(Var, Par), units='TgC yr-1')
def Eq__D_Fsink_CH4(Var, Par):
    return Var.D_Foh_CH4 + Var.D_Fhv_CH4 + Var.D_Fsoil_CH4 + Var.D_Focean_CH4


## flux of geological CH4 oxidized in CO2 (in CO2 units)
D_Foxi_CH4 = OSCAR.process('D_Foxi_CH4', ('E_CH4', 'D_Ewet', 'D_Ebb', 'D_Fsink_CH4'), lambda Var, Par: Eq__D_Foxi_CH4(Var, Par), units='TgC yr-1')
def Eq__D_Foxi_CH4(Var, Par):
    a_conv = 1E-3 # from {TgC} to {PgC}
    E_CH4_nongeo = (1 - Par.p_CH4geo) * (Var.E_CH4 + Var.D_Ewet + Var.D_Ebb.sel({'spc_bb':'CH4'}, drop=True).sum('bio_land', min_count=1)).sum('reg_land', min_count=1)
    return a_conv * (Var.D_Fsink_CH4 - E_CH4_nongeo)


##==============
## 3.2. Wetlands
##==============

## areal wetland emissions
D_ewet = OSCAR.process('D_ewet', ('D_rh1', 'D_rh2', 'csoil1_0', 'csoil2_0'), lambda Var, Par: Eq__D_ewet(Var, Par), units='TgC yr-1 Mha-1')
def Eq__D_ewet(Var, Par):
    return Par.ewet_0 * (Par.p_wet * (Var.D_rh1 + Var.D_rh2)).sum('bio_land', min_count=1) / (Par.p_wet * (Par.rho1_0 * Var.csoil1_0 + Par.rho2_0 * Var.csoil2_0)).sum('bio_land', min_count=1)


## wetland extent
D_Awet = OSCAR.process('D_Awet', ('D_CO2', 'D_Tl', 'D_Pl'), lambda Var, Par: Eq__D_Awet(Var, Par), units='Mha')
def Eq__D_Awet(Var, Par):
    return Par.Awet_0 * (Par.g_wetC * Var.D_CO2 + Par.g_wetT * Var.D_Tl + Par.g_wetP * Var.D_Pl)


## wetland emissions
D_Ewet = OSCAR.process('D_Ewet', ('D_ewet', 'D_Awet'), lambda Var, Par: Eq__D_Ewet(Var, Par), units='TgC yr-1')
def Eq__D_Ewet(Var, Par):
    return Par.ewet_0 * Var.D_Awet + Var.D_ewet * (Par.Awet_0 + Var.D_Awet)


##================
## 3.3. Atmosphere
##================

## PROGNOSTIC: atmospheric CH4
D_CH4 = OSCAR.process('D_CH4', ('D_CH4', 'E_CH4', 'D_Ewet', 'D_Ebb', 'D_Epf_CH4', 'D_Fsink_CH4'), lambda Var, Par, **Int_args: DiffEq__D_CH4(Var, Par, **Int_args), units='ppb')
def DiffEq__D_CH4(Var, Par, Int=Int_dflt, dt=1.):
    v = 1 / Par.w_t_OH / Par.t_OH_CH4 + 1 / Par.w_t_hv / Par.t_hv_CH4 + 1 / Par.t_soil_CH4 + 1 / Par.t_ocean_CH4
    return Int(Var.D_CH4, v, 1/Par.a_CH4 * ((Var.E_CH4 + Var.D_Ewet + Var.D_Ebb.sel({'spc_bb':'CH4'}, drop=True).sum('bio_land', min_count=1)).sum('reg_land', min_count=1) + Var.D_Epf_CH4.sum('reg_pf', min_count=1) - Var.D_Fsink_CH4), dt)


## METRIC: CH4 lifetime
tau_CH4 = OSCAR.process('tau_CH4', ('D_Fsink_CH4', 'D_CH4'), lambda Var, Par: Eq__tau_CH4(Var, Par), units='yr')
def Eq__tau_CH4(Var, Par):
    v = 1 / Par.w_t_OH / Par.t_OH_CH4 + 1 / Par.w_t_hv / Par.t_hv_CH4 + 1 / Par.t_soil_CH4 + 1 / Par.t_ocean_CH4
    return Par.a_CH4 * (Par.CH4_0 + Var.D_CH4) / (v * Par.a_CH4 * Par.CH4_0 + Var.D_Fsink_CH4)


## CH4 radiative forcing
## (Myhre et al., 1998; doi:10.1029/98GL01908)
RF_CH4 = OSCAR.process('RF_CH4', ('D_CH4', 'D_N2O'), lambda Var, Par: Eq__RF_CH4(Var, Par), units='W m-2')
def Eq__RF_CH4(Var, Par):
    RF_overlap = 0.47 * np.log1p(2.01E-5 * ((Par.CH4_0 + Var.D_CH4) * Par.N2O_0)**0.75 + 5.31E-15 * (Par.CH4_0 + Var.D_CH4) * ((Par.CH4_0 + Var.D_CH4) * (Par.N2O_0))**1.52)
    RF_overlap -= 0.47 * np.log1p(2.01E-5 * (Par.CH4_0  * Par.N2O_0)**0.75 + 5.31E-15 * Par.CH4_0 * (Par.CH4_0 * Par.N2O_0)**1.52)
    return Par.rf_CH4 * (np.sqrt(Par.CH4_0 + Var.D_CH4) - np.sqrt(Par.CH4_0)) - RF_overlap


## stratospheric H2O radiative forcing
RF_H2Os = OSCAR.process('RF_H2Os', ('D_CH4_lag',), lambda Var, Par: Eq__RF_H2Os(Var, Par), units='W m-2')
def Eq__RF_H2Os(Var, Par):
    return Par.k_rf_H2Os * Par.rf_CH4 * (np.sqrt(Par.CH4_0 + Var.D_CH4_lag) - np.sqrt(Par.CH4_0))


##################################################
##   4. NITROUS OXIDE
##################################################

##====================
## 4.1. Atmo chemistry
##====================

## stratospheric mean age of air
D_f_ageair = OSCAR.process('D_f_ageair', ('D_Tg',), lambda Var, Par: Eq__D_f_ageair(Var, Par), units='1')
def Eq__D_f_ageair(Var, Par):
    return -Par.g_ageair * Var.D_Tg / (1 + Par.g_ageair * Var.D_Tg)


## stratospheric sink intensity
f_hv = OSCAR.process('f_hv', ('D_N2O_lag', 'D_EESC', 'EESC_0', 'D_f_ageair'), lambda Var, Par: Eq__f_hv(Var, Par), units='1')
def Eq__f_hv(Var, Par):
    f_hv_N2O = Par.x_hv_N2O * np.log1p(Var.D_N2O_lag / Par.N2O_0)
    f_hv_EESC = Par.x_hv_EESC * np.log1p(Var.D_EESC / Var.EESC_0)
    f_hv_Tg = Par.x_hv_ageair * np.log1p(Var.D_f_ageair)
    return np.exp(f_hv_N2O + f_hv_EESC + f_hv_Tg)


## N2O stratospheric sink
D_Fhv_N2O = OSCAR.process('D_Fhv_N2O', ('f_hv', 'D_N2O_lag'), lambda Var, Par: Eq__D_Fhv_N2O(Var, Par), units='TgN yr-1')
def Eq__D_Fhv_N2O(Var, Par):
    return Par.a_N2O / Par.w_t_hv / Par.t_hv_N2O * ((Par.N2O_0 + Var.D_N2O_lag) * Var.f_hv - Par.N2O_0)


## N2O total atmospheric sink
D_Fsink_N2O = OSCAR.process('D_Fsink_N2O', ('D_Fhv_N2O',), lambda Var, Par: Eq__D_Fsink_N2O(Var, Par), units='TgN yr-1')
def Eq__D_Fsink_N2O(Var, Par):
    return Var.D_Fhv_N2O


##================
## 4.2. Atmosphere
##================

## PROGNOSTIC: atmospheric N2O
D_N2O = OSCAR.process('D_N2O', ('D_N2O', 'E_N2O', 'D_Ebb', 'D_Fsink_N2O'), lambda Var, Par, **Int_args: DiffEq__D_N2O(Var, Par, **Int_args), units='ppb')
def DiffEq__D_N2O(Var, Par, Int=Int_dflt, dt=1.):
    v = 1 / Par.w_t_hv / Par.t_hv_N2O
    return Int(Var.D_N2O, v, 1/Par.a_N2O * ((Var.E_N2O + Var.D_Ebb.sel({'spc_bb':'N2O'}, drop=True).sum('bio_land', min_count=1)).sum('reg_land', min_count=1) - Var.D_Fsink_N2O), dt)


## METRIC: N2O lifetime
tau_N2O = OSCAR.process('tau_N2O', ('D_Fsink_N2O', 'D_N2O'), lambda Var, Par: Eq__tau_N2O(Var, Par), units='yr')
def Eq__tau_N2O(Var, Par):
    return Par.a_N2O * (Par.N2O_0 + Var.D_N2O) / (1 / Par.w_t_hv / Par.t_hv_N2O * Par.a_N2O * Par.N2O_0 + Var.D_Fsink_N2O)


## N2O radiative forcing
## (Myhre et al., 1998; doi:10.1029/98GL01908)
RF_N2O = OSCAR.process('RF_N2O', ('D_CH4', 'D_N2O'), lambda Var, Par: Eq__RF_N2O(Var, Par), units='W m-2')
def Eq__RF_N2O(Var, Par):
    RF_overlap = 0.47 * np.log1p(2.01E-5 * (Par.CH4_0 * (Par.N2O_0 + Var.D_N2O))**0.75 + 5.31E-15 * Par.CH4_0 * (Par.CH4_0 * (Par.N2O_0 + Var.D_N2O))**1.52)
    RF_overlap -= 0.47 * np.log1p(2.01E-5 * (Par.CH4_0 * Par.N2O_0)**0.75 + 5.31E-15 * Par.CH4_0 * (Par.CH4_0 * Par.N2O_0)**1.52)
    return Par.rf_N2O * (np.sqrt(Par.N2O_0 + Var.D_N2O) - np.sqrt(Par.N2O_0)) - RF_overlap


##################################################
##   5. HALOGENATED COMPOUNDS
##################################################

##====================
## 5.1. Atmo chemistry
##====================

## Xhalo hydroxyl sink
D_Foh_Xhalo = OSCAR.process('D_Foh_Xhalo', ('f_kOH', 'D_Xhalo'), lambda Var, Par: Eq__D_Foh_Xhalo(Var, Par), units='Gg yr-1')
def Eq__D_Foh_Xhalo(Var, Par):
    return Par.a_Xhalo / Par.w_t_OH / Par.t_OH_Xhalo * ((Par.Xhalo_0 + Var.D_Xhalo) * Var.f_kOH - Par.Xhalo_0)


## Xhalo stratospheric sink
D_Fhv_Xhalo = OSCAR.process('D_Fhv_Xhalo', ('f_hv', 'D_Xhalo_lag'), lambda Var, Par: Eq__D_Fhv_Xhalo(Var, Par), units='Gg yr-1')
def Eq__D_Fhv_Xhalo(Var, Par):
    return Par.a_Xhalo / Par.w_t_hv / Par.t_hv_Xhalo * ((Par.Xhalo_0 + Var.D_Xhalo_lag) * Var.f_hv - Par.Xhalo_0)


## Xhalo other (constant) sinks
D_Fother_Xhalo = OSCAR.process('D_Fother_Xhalo', ('D_Xhalo',), lambda Var, Par: Eq__D_Fother_Xhalo(Var, Par), units='Gg yr-1')
def Eq__D_Fother_Xhalo(Var, Par):
    return Par.a_Xhalo / Par.t_other_Xhalo * Var.D_Xhalo


## Xhalo total atmospheric sink
D_Fsink_Xhalo = OSCAR.process('D_Fsink_Xhalo', ('D_Foh_Xhalo', 'D_Fhv_Xhalo', 'D_Fother_Xhalo'), lambda Var, Par: Eq__D_Fsink_Xhalo(Var, Par), units='Gg yr-1')
def Eq__D_Fsink_Xhalo(Var, Par):
    return Var.D_Foh_Xhalo + Var.D_Fhv_Xhalo + Var.D_Fother_Xhalo


##================
## 5.2. Atmosphere
##================

## PROGNOSTIC: atmospheric Xhalo
D_Xhalo = OSCAR.process('D_Xhalo', ('D_Xhalo', 'E_Xhalo', 'D_Fsink_Xhalo'), lambda Var, Par, **Int_args: DiffEq__D_Xhalo(Var, Par, **Int_args), units='ppt', core_dims=['spc_halo'])
def DiffEq__D_Xhalo(Var, Par, Int=Int_dflt, dt=1.):
    v = 1 / Par.w_t_OH / Par.t_OH_Xhalo + 1 / Par.w_t_hv / Par.t_hv_Xhalo + 1/Par.t_other_Xhalo
    return Int(Var.D_Xhalo, v, 1/Par.a_Xhalo * (Var.E_Xhalo.sum('reg_land', min_count=1) - Var.D_Fsink_Xhalo), dt)


## radiative forcing of individual halogenated compounds
RF_Xhalo = OSCAR.process('RF_Xhalo', ('D_Xhalo',), lambda Var, Par: Eq__RF_Xhalo(Var, Par), units='W m-2')
def Eq__RF_Xhalo(Var, Par):
    return Par.rf_Xhalo * Var.D_Xhalo


## radiative forcing of combined halogenated compounds
RF_halo = OSCAR.process('RF_halo', ('RF_Xhalo',), lambda Var, Par: Eq__RF_halo(Var, Par), units='W m-2')
def Eq__RF_halo(Var, Par):
    return Var.RF_Xhalo.sum('spc_halo', min_count=1)


##################################################
##   6. OZONE
##################################################

##=================
## 6.1. Troposphere
##=================

## module adapted from:
## (Ehhalt et al., 2001; IPCC AR3 WG1 Chapter 4)

## O3 tropospheric burden
D_O3t = OSCAR.process('D_O3t', ('D_CH4', 'D_Tg', 'E_NOX', 'E_CO', 'E_VOC', 'D_Ebb'), lambda Var, Par: Eq__D_O3t(Var, Par), units='DU')
def Eq__D_O3t(Var, Par):
    D_O3t_CH4 = Par.x_O3t_CH4 * np.log1p(Var.D_CH4 / Par.CH4_0)
    D_O3t_Tg = Par.G_O3t * Var.D_Tg
    D_O3t_NOX = Par.x_O3t_NOX * (Par.w_reg_NOX * (Par.p_reg_slcf * (Var.E_NOX + Var.D_Ebb.sel({'spc_bb':'NOX'}, drop=True).sum('bio_land', min_count=1))).sum('reg_land', min_count=1)).sum('reg_slcf', min_count=1)
    D_O3t_CO = Par.x_O3t_CO * (Par.w_reg_CO * (Par.p_reg_slcf * (Var.E_CO + Var.D_Ebb.sel({'spc_bb':'CO'}, drop=True).sum('bio_land', min_count=1))).sum('reg_land', min_count=1)).sum('reg_slcf', min_count=1)
    D_O3t_VOC = Par.x_O3t_VOC * (Par.w_reg_VOC * (Par.p_reg_slcf * (Var.E_VOC + Var.D_Ebb.sel({'spc_bb':'VOC'}, drop=True).sum('bio_land', min_count=1))).sum('reg_land', min_count=1)).sum('reg_slcf', min_count=1)
    return D_O3t_CH4 + D_O3t_Tg + D_O3t_NOX + D_O3t_CO + D_O3t_VOC


## tropospheric O3 radiative forcing
RF_O3t = OSCAR.process('RF_O3t', ('D_O3t',), lambda Var, Par: Eq__RF_O3t(Var, Par), units='W m-2')
def Eq__RF_O3t(Var, Par):
    return Par.rf_O3t * Var.D_O3t


##==================
## 6.2. Stratosphere
##==================

## PARAMETER: preindustrial equivalent effective stratospheric chlorine
## (Newman et al., 2017; doi:10.5194/acp-7-4537-2007)
EESC_0 = OSCAR.process('EESC_0', (), lambda Var, Par: Eq__EESC_0(Var, Par), units='ppt')
def Eq__EESC_0(Var, Par):
    return (Par.p_fracrel * (Par.n_Cl + Par.k_Br_Cl * Par.n_Br) * Par.Xhalo_0).sum('spc_halo', min_count=1)


## equivalent effective stratospheric chlorine
## (Newman et al., 2017; doi:10.5194/acp-7-4537-2007)
D_EESC = OSCAR.process('D_EESC', ('D_Xhalo_lag',), lambda Var, Par: Eq__D_EESC(Var, Par), units='ppt')
def Eq__D_EESC(Var, Par):
    return (Par.p_fracrel * (Par.n_Cl + Par.k_Br_Cl * Par.n_Br) * Var.D_Xhalo_lag).sum('spc_halo', min_count=1)


## O3 stratospheric burden
D_O3s = OSCAR.process('D_O3s', ('D_EESC', 'D_N2O_lag', 'D_Tg'), lambda Var, Par: Eq__D_O3s(Var, Par), units='DU')
def Eq__D_O3s(Var, Par):
    D_O3s_EESC = Par.x_O3s_EESC * Var.D_EESC
    D_O3s_N2O = Par.x_O3s_EESC * Par.k_EESC_N2O * (1 - Var.D_EESC / Par.EESC_x) * Var.D_N2O_lag
    D_O3s_Tg = Par.G_O3s * Var.D_Tg
    return D_O3s_EESC + D_O3s_N2O + D_O3s_Tg


## stratospheric O3 radiative forcing
RF_O3s = OSCAR.process('RF_O3s', ('D_O3s',), lambda Var, Par: Eq__RF_O3s(Var, Par), units='W m-2')
def Eq__RF_O3s(Var, Par):
    return Par.rf_O3s * Var.D_O3s


##################################################
##   7. AEROSOLS
##################################################

##=======================
## 7.1. Natural emissions
##=======================

## PARAMETER: natural emissions of DMS
D_Edms = OSCAR.process('D_Edms', (), lambda Var, Par: Eq__D_Edms(Var, Par), units='TgS yr-1')
def Eq__D_Edms(Var, Par):
    return 0.


## PARAMETER: natural emissions of biogenic VOC
D_Ebvoc = OSCAR.process('D_Ebvoc', (), lambda Var, Par: Eq__D_Ebvoc(Var, Par), units='Tg yr-1')
def Eq__D_Ebvoc(Var, Par):
    return 0.


## PARAMETER: natural emissions of mineral dust
D_Edust = OSCAR.process('D_Edust', (), lambda Var, Par: Eq__D_Edust(Var, Par), units='Tg yr-1')
def Eq__D_Edust(Var, Par):
    return 0.


## PARAMETER: natural emissions of sea salt
D_Esalt = OSCAR.process('D_Esalt', (), lambda Var, Par: Eq__D_Esalt(Var, Par), units='Tg yr-1')
def Eq__D_Esalt(Var, Par):
    return 0.


##===================
## 7.2. Direct effect
##===================

## SO4 tropospheric burden
D_SO4 = OSCAR.process('D_SO4', ('E_SO2', 'D_Ebb', 'D_Edms', 'D_Tg'), lambda Var, Par: Eq__D_SO4(Var, Par), units='Tg')
def Eq__D_SO4(Var, Par):
    D_SO4_SO2 = Par.a_SO4 * Par.t_SO2 * (Par.w_reg_SO2 * (Par.p_reg_slcf * (Var.E_SO2 + Var.D_Ebb.sel({'spc_bb':'SO2'}, drop=True).sum('bio_land', min_count=1))).sum('reg_land', min_count=1)).sum('reg_slcf', min_count=1)
    D_SO4_DMS = Par.a_SO4 * Par.t_DMS * Var.D_Edms
    D_SO4_Tg = Par.G_SO4 * Var.D_Tg
    return D_SO4_SO2 + D_SO4_DMS + D_SO4_Tg


## POA tropospheric burden
D_POA = OSCAR.process('D_POA', ('E_OC', 'D_Ebb', 'D_Tg'), lambda Var, Par: Eq__D_POA(Var, Par), units='Tg')
def Eq__D_POA(Var, Par):
    D_POA_OMff = Par.a_POM * Par.t_OMff * (Par.w_reg_OC * (Par.p_reg_slcf * Var.E_OC).sum('reg_land', min_count=1)).sum('reg_slcf', min_count=1)
    D_POA_OMbb = Par.a_POM * Par.t_OMbb * Var.D_Ebb.sel({'spc_bb':'OC'}, drop=True).sum('bio_land', min_count=1).sum('reg_land', min_count=1)
    D_POA_Tg = Par.G_POA * Var.D_Tg
    return D_POA_OMff + D_POA_OMbb + D_POA_Tg


## BC tropospheric burden
D_BC = OSCAR.process('D_BC', ('E_BC', 'D_Ebb', 'D_Tg'), lambda Var, Par: Eq__D_BC(Var, Par), units='Tg')
def Eq__D_BC(Var, Par):
    D_BC_BCff = Par.t_BCff * (Par.w_reg_BC * (Par.p_reg_slcf * Var.E_BC).sum('reg_land', min_count=1)).sum('reg_slcf', min_count=1)
    D_BC_BCbb = Par.t_BCbb * Var.D_Ebb.sel({'spc_bb':'BC'}, drop=True).sum('bio_land', min_count=1).sum('reg_land', min_count=1)
    D_BC_Tg = Par.G_BC * Var.D_Tg
    return D_BC_BCff + D_BC_BCbb + D_BC_Tg


## NO3 tropospheric burden
D_NO3 = OSCAR.process('D_NO3', ('E_NOX', 'E_NH3', 'D_Ebb', 'D_Tg'), lambda Var, Par: Eq__D_NO3(Var, Par), units='Tg')
def Eq__D_NO3(Var, Par):
    D_NO3_NOX = Par.a_NO3 * Par.t_NOX * (Var.E_NOX + Var.D_Ebb.sel({'spc_bb':'NOX'}, drop=True).sum('bio_land', min_count=1)).sum('reg_land', min_count=1)
    D_NO3_NH3 = Par.a_NO3 * Par.t_NH3 * (Var.E_NH3 + Var.D_Ebb.sel({'spc_bb':'NH3'}, drop=True).sum('bio_land', min_count=1)).sum('reg_land', min_count=1)
    D_NO3_Tg = Par.G_NO3 * Var.D_Tg
    return D_NO3_NOX + D_NO3_NH3 + D_NO3_Tg


## SOA tropospheric burden
D_SOA = OSCAR.process('D_SOA', ('E_VOC', 'D_Ebvoc', 'D_Ebb', 'D_Tg'), lambda Var, Par: Eq__D_SOA(Var, Par), units='Tg')
def Eq__D_SOA(Var, Par):
    D_SOA_VOC = Par.t_VOC * (Var.E_VOC + Var.D_Ebb.sel({'spc_bb':'VOC'}, drop=True).sum('bio_land', min_count=1)).sum('reg_land', min_count=1)
    D_SOA_BVOC = Par.t_BVOC * Var.D_Ebvoc
    D_SOA_Tg = Par.G_SOA * Var.D_Tg
    return D_SOA_VOC + D_SOA_BVOC + D_SOA_Tg


## mineral dust tropospheric burden
D_Mdust = OSCAR.process('D_Mdust', ('D_Edust', 'D_Tg'), lambda Var, Par: Eq__D_Mdust(Var, Par), units='Tg')
def Eq__D_Mdust(Var, Par):
    return Par.t_dust * Var.D_Edust + Par.G_dust * Var.D_Tg


## sea salt tropospheric burden
D_Msalt = OSCAR.process('D_Msalt', ('D_Esalt', 'D_Tg'), lambda Var, Par: Eq__D_Msalt(Var, Par), units='Tg')
def Eq__D_Msalt(Var, Par):
    return Par.t_salt * Var.D_Esalt + Par.G_salt * Var.D_Tg


## SO4 radiative forcing (direct)
RF_SO4 = OSCAR.process('RF_SO4', ('D_SO4',), lambda Var, Par: Eq__RF_SO4(Var, Par), units='W m-2')
def Eq__RF_SO4(Var, Par):
    return Par.rf_SO4 * Var.D_SO4


## POA radiative forcing (direct)
RF_POA = OSCAR.process('RF_POA', ('D_POA',), lambda Var, Par: Eq__RF_POA(Var, Par), units='W m-2')
def Eq__RF_POA(Var, Par):
    return Par.rf_POA * Var.D_POA


## BC radiative forcing (direct)
RF_BC = OSCAR.process('RF_BC', ('D_BC',), lambda Var, Par: Eq__RF_BC(Var, Par), units='W m-2')
def Eq__RF_BC(Var, Par):
    return Par.rf_BC * Var.D_BC


## NO3 radiative forcing (direct)
RF_NO3 = OSCAR.process('RF_NO3', ('D_NO3',), lambda Var, Par: Eq__RF_NO3(Var, Par), units='W m-2')
def Eq__RF_NO3(Var, Par):
    return Par.rf_NO3 * Var.D_NO3


## SOA radiative forcing (direct)
RF_SOA = OSCAR.process('RF_SOA', ('D_SOA',), lambda Var, Par: Eq__RF_SOA(Var, Par), units='W m-2')
def Eq__RF_SOA(Var, Par):
    return Par.rf_SOA * Var.D_SOA


## mineral dust radiative forcing (direct)
RF_dust = OSCAR.process('RF_dust', ('D_Mdust',), lambda Var, Par: Eq__RF_dust(Var, Par), units='W m-2')
def Eq__RF_dust(Var, Par):
    return Par.rf_dust * Var.D_Mdust


## sea salt radiative forcing (direct)
RF_salt = OSCAR.process('RF_salt', ('D_Msalt',), lambda Var, Par: Eq__RF_salt(Var, Par), units='W m-2')
def Eq__RF_salt(Var, Par):
    return Par.rf_salt * Var.D_Msalt


##===================
## 7.3. Cloud effects
##===================

## hydrophilic aerosol tropospheric burden
D_AERsol = OSCAR.process('D_AERsol', ('D_SO4', 'D_POA', 'D_BC', 'D_NO3', 'D_SOA', 'D_Mdust', 'D_Msalt'), lambda Var, Par: Eq__D_AERsol(Var, Par), units='Tg')
def Eq__D_AERsol(Var, Par):
    return Par.p_sol_SO4 * Var.D_SO4 + Par.p_sol_POA * Var.D_POA + Par.p_sol_BC * Var.D_BC + Par.p_sol_NO3 * Var.D_NO3 + Par.p_sol_SOA * Var.D_SOA + Par.p_sol_dust * Var.D_Mdust + Par.p_sol_salt * Var.D_Msalt


## semi-direct effect (= adjustment of aerosol-radiation interactions)
RF_cloud1 = OSCAR.process('RF_cloud1', ('RF_BC',), lambda Var, Par: Eq__RF_cloud1(Var, Par), units='W m-2')
def Eq__RF_cloud1(Var, Par):
    return Par.k_adj_BC * Var.RF_BC


## second indirect effect (= aerosol-cloud interactions)
RF_cloud2 = OSCAR.process('RF_cloud2', ('D_AERsol',), lambda Var, Par: Eq__RF_cloud2(Var, Par), units='W m-2')
def Eq__RF_cloud2(Var, Par):
    return Par.Phi_0 * np.log1p(Var.D_AERsol / Par.AERsol_0)


## radiative forcing of cloud effects
RF_cloud = OSCAR.process('RF_cloud', ('RF_cloud1', 'RF_cloud2'), lambda Var, Par: Eq__RF_cloud(Var, Par), units='W m-2')
def Eq__RF_cloud(Var, Par):
    return Var.RF_cloud1 + Var.RF_cloud2


##################################################
##   8. SURFACE ALBEDO
##################################################

##==========================
## 8.1. Black carbon on snow
##==========================

## radiative forcing from BC deposition on snow
RF_BCsnow = OSCAR.process('RF_BCsnow', ('E_BC', 'D_Ebb'), lambda Var, Par: Eq__RF_BCsnow(Var, Par), units='W m-2')
def Eq__RF_BCsnow(Var, Par):
    return Par.rf_bcsnow * (Par.w_reg_bcsnow * (Par.p_reg_bcsnow * (Var.E_BC + Var.D_Ebb.sel({'spc_bb':'BC'}, drop=True).sum('bio_land', min_count=1))).sum('reg_land', min_count=1)).sum('reg_bcsnow', min_count=1)


##=======================
## 8.2. Land-cover change
##=======================

## radiative forcing from land-cover change
## (Bright & Kvalevag, 2013; doi:10.5194/acp-13-11169-2013)
RF_lcc = OSCAR.process('RF_lcc', ('D_Aland',), lambda Var, Par: Eq__RF_lcc(Var, Par), units='W m-2')
def Eq__RF_lcc(Var, Par):
    A_Earth = 510072E-1 # Mha
    return -Par.p_trans / A_Earth * (Par.F_rsds * (Par.a_alb * Var.D_Aland).sum('bio_land', min_count=1)).sum('reg_land', min_count=1)


##################################################
##   9. CLIMATE
##################################################

##=======================
## 9.1. Radiative forcing
##=======================

## radiative forcing of non-CO2 WMGHGs
RF_nonCO2 = OSCAR.process('RF_nonCO2', ('RF_CH4', 'RF_N2O', 'RF_halo'), lambda Var, Par: Eq__RF_nonCO2(Var, Par), units='W m-2')
def Eq__RF_nonCO2(Var, Par):
    return Var.RF_CH4 + Var.RF_N2O + Var.RF_halo


## radiative forcing of all WMGHGs
RF_wmghg = OSCAR.process('RF_wmghg', ('RF_CO2', 'RF_nonCO2'), lambda Var, Par: Eq__RF_wmghg(Var, Par), units='W m-2')
def Eq__RF_wmghg(Var, Par):
    return Var.RF_CO2 + Var.RF_nonCO2


## radiative forcing of stratospheric GHGs
RF_strat = OSCAR.process('RF_strat', ('RF_H2Os', 'RF_O3s'), lambda Var, Par: Eq__RF_strat(Var, Par), units='W m-2')
def Eq__RF_strat(Var, Par):
    return Var.RF_H2Os + Var.RF_O3s


## radiative forcing of scattering aerosols
RF_scatter = OSCAR.process('RF_scatter', ('RF_SO4', 'RF_POA', 'RF_NO3', 'RF_SOA', 'RF_dust', 'RF_salt'), lambda Var, Par: Eq__RF_scatter(Var, Par), units='W m-2')
def Eq__RF_scatter(Var, Par):
    return Var.RF_SO4 + Var.RF_POA + Var.RF_NO3 + Var.RF_SOA + Var.RF_dust + Var.RF_salt


## radiative forcing of absorbing aerosols
RF_absorb = OSCAR.process('RF_absorb', ('RF_BC',), lambda Var, Par: Eq__RF_absorb(Var, Par), units='W m-2')
def Eq__RF_absorb(Var, Par):
    return Var.RF_BC

    
## radiative forcing of all aerosols effects
RF_AERtot = OSCAR.process('RF_AERtot', ('RF_scatter', 'RF_absorb', 'RF_cloud'), lambda Var, Par: Eq__RF_AERtot(Var, Par), units='W m-2')
def Eq__RF_AERtot(Var, Par):
    return Var.RF_scatter + Var.RF_absorb + Var.RF_cloud
   

## radiative forcing of all SLCFs
RF_slcf = OSCAR.process('RF_slcf', ('RF_O3t', 'RF_strat', 'RF_AERtot'), lambda Var, Par: Eq__RF_slcf(Var, Par), units='W m-2')
def Eq__RF_slcf(Var, Par):
    return Var.RF_O3t + Var.RF_strat + Var.RF_AERtot


## radiative forcing of albedo effects
RF_alb = OSCAR.process('RF_alb', ('RF_BCsnow', 'RF_lcc'), lambda Var, Par: Eq__RF_alb(Var, Par), units='W m-2')
def Eq__RF_alb(Var, Par):
    return Var.RF_BCsnow + Var.RF_lcc


## total radiative forcing (= TOA RF)
RF = OSCAR.process('RF', ('RF_wmghg', 'RF_slcf', 'RF_alb', 'RF_volc', 'RF_solar', 'RF_contr'), lambda Var, Par: Eq__RF(Var, Par), units='W m-2')
def Eq__RF(Var, Par):
    return Var.RF_wmghg + Var.RF_slcf + Var.RF_alb + Var.RF_volc + Var.RF_solar + Var.RF_contr


## warming radiative forcing (~= ERF)
RF_warm = OSCAR.process('RF_warm', ('RF_wmghg', 'RF_slcf', 'RF_BCsnow', 'RF_lcc', 'RF_volc', 'RF_solar', 'RF_contr'), lambda Var, Par: Eq__RF_warm(Var, Par), units='W m-2')
def Eq__RF_warm(Var, Par):
    return Var.RF_wmghg + Var.RF_slcf + Par.w_warm_bcsnow * Var.RF_BCsnow + Par.w_warm_lcc * Var.RF_lcc + Par.w_warm_volc * Var.RF_volc + Var.RF_solar + Var.RF_contr


## atmospheric fraction of radiative forcing
RF_atm = OSCAR.process('RF_atm', ('RF_CO2', 'RF_nonCO2', 'RF_O3t', 'RF_strat', 'RF_scatter', 'RF_absorb', 'RF_cloud', 'RF_alb', 'RF_volc', 'RF_solar', 'RF_contr'), lambda Var, Par: Eq__RF_atm(Var, Par), units='W m-2')
def Eq__RF_atm(Var, Par):
    RF_atm_CO2 = Par.p_atm_CO2 * Var.RF_CO2
    RF_atm_nonCO2 = Par.p_atm_nonCO2 * Var.RF_nonCO2
    RF_atm_O3t = Par.p_atm_O3t * Var.RF_O3t
    RF_atm_strat = Par.p_atm_strat * Var.RF_strat
    RF_atm_scatter = Par.p_atm_scatter * (Var.RF_scatter + Var.RF_volc)
    RF_atm_absorb = Par.p_atm_absorb * Var.RF_absorb
    RF_atm_cloud = Par.p_atm_cloud * (Var.RF_cloud + Var.RF_contr)
    RF_atm_alb = Par.p_atm_alb * Var.RF_alb
    RF_atm_solar = Par.p_atm_solar * Var.RF_solar
    return RF_atm_CO2 + RF_atm_nonCO2 + RF_atm_O3t + RF_atm_strat + RF_atm_scatter + RF_atm_absorb + RF_atm_cloud + RF_atm_alb + RF_atm_solar


##=================
## 9.2. Temperature
##=================

## PROGNOSTIC: global mean surface temperature
## (Geoffroy et al., 2013; doi:10.1175/JCLI-D-12-00195.1)
D_Tg = OSCAR.process('D_Tg', ('D_Tg', 'D_Td', 'RF_warm'), lambda Var, Par, **Int_args: DiffEq__D_Tg(Var, Par, **Int_args), units='K')
def DiffEq__D_Tg(Var, Par, Int=Int_dflt, dt=1.):
    v = 1/Par.Th_g * (1/Par.lambda_0 + Par.th_0)
    return Int(Var.D_Tg, v, 1/Par.Th_g * (Var.RF_warm - Var.D_Tg / Par.lambda_0 - Par.th_0 * (Var.D_Tg - Var.D_Td)), dt)


## PROGNOSTIC: deep ocean temperature
## (Geoffroy et al., 2013; doi:10.1175/JCLI-D-12-00195.1)
D_Td = OSCAR.process('D_Td', ('D_Td', 'D_Tg'), lambda Var, Par, **Int_args: DiffEq__D_Td(Var, Par, **Int_args), units='K')
def DiffEq__D_Td(Var, Par, Int=Int_dflt, dt=1.):
    v = Par.th_0 / Par.Th_d
    return Int(Var.D_Td, v, 1/Par.Th_d * Par.th_0 * (Var.D_Tg - Var.D_Td), dt)


## METRIC: global mean surface temperature rate of change
d_Tg = OSCAR.process('d_Tg', ('D_Tg', 'D_Td', 'RF_warm'), lambda Var, Par: Eq__d_Tg(Var, Par), units='K yr-1')
def Eq__d_Tg(Var, Par):
    return 1/Par.Th_g * (Var.RF_warm - Var.D_Tg / Par.lambda_0 - Par.th_0 * (Var.D_Tg - Var.D_Td))


## land surface temperature
D_Tl = OSCAR.process('D_Tl', ('D_Tg',), lambda Var, Par: Eq__D_Tl(Var, Par), units='K')
def Eq__D_Tl(Var, Par):
    return Par.w_clim_Tl * Var.D_Tg


## sea surface temperature
D_To = OSCAR.process('D_To', ('D_Tg',), lambda Var, Par: Eq__D_To(Var, Par), units='K')
def Eq__D_To(Var, Par):
    return Par.w_clim_To * Var.D_Tg


##===================
## 9.3. Precipitation
##===================

## NODE: global precipitation
## (Allan et al., 2013; doi:10.1007/s10712-012-9213-z)
D_Pg = OSCAR.process('D_Pg', ('D_Pg', 'D_Tg', 'RF_atm'), lambda Var, Par, **Int_args: Eq__D_Pg(Var, Par), units='mm yr-1')
def Eq__D_Pg(Var, Par):
    return Par.a_prec * Var.D_Tg + Par.b_prec / Par.p_atm_CO2 * Var.RF_atm


## land precipitation
D_Pl = OSCAR.process('D_Pl', ('D_Pg',), lambda Var, Par: Eq__D_Pl(Var, Par), units='mm yr-1')
def Eq__D_Pl(Var, Par):
    return Par.w_clim_Pl * Var.D_Pg


##========================
## 9.4. Ocean heat content
##========================

## PROGNOSTIC: ocean heat content
D_OHC = OSCAR.process('D_OHC', ('D_OHC', 'D_Tg', 'RF'), lambda Var, Par, **Int_args: DiffEq__D_OHC(Var, Par, **Int_args), units='ZJ')
def DiffEq__D_OHC(Var, Par, Int=Int_dflt, dt=1.):
    a_conv = 3600*24*365.25 / 1E21 # from {W yr} to {ZJ}
    A_Earth = 510072E9 # m2
    return Int(Var.D_OHC, 1E-18, a_conv * A_Earth * Par.p_ohc * (Var.RF - Var.D_Tg / Par.lambda_0), dt)


##################################################
##   10. IMPACTS
##################################################

##====================
## 10.1. Acidification
##====================

## global surface ocean acidification
D_pH = OSCAR.process('D_pH', ('D_CO2',), lambda Var, Par: Eq__D_pH(Var, Par), units='1')
def Eq__D_pH(Var, Par):
    ## logarithmic formulation (approximative)
    ## (Tans 2009; doi:10.5670/oceanog.2009.94)
    D_pH_log = -0.85 * np.log1p(Var.D_CO2 / Par.CO2_0)
    ## polynomial formulation
    ## (Bernie et al., 2010; doi:10.1029/2010GL043181)
    D_pH_poly = -0.00173 * Var.D_CO2 + 1.3264E-6 * (2*Par.CO2_0 * Var.D_CO2 + Var.D_CO2**2) - 4.4943E-10 * (3*Var.D_CO2 * Par.CO2_0**2 + 3*Par.CO2_0 * Var.D_CO2**2 + Var.D_CO2**3)
    ## choose configuration
    return Par.pH_is_Log * D_pH_log + (1-Par.pH_is_Log) * D_pH_poly


##################################################
##   B. SUB-MODELS
##################################################

## ocean carbon cycle only
proc_oceanC = ['dic_0', 'D_pCO2', 'D_mld', 'D_dic', 'D_Fin', 'D_Fout', 'D_Fcirc', 'D_Focean', 'D_Cosurf']
OSCAR_oceanC = OSCAR.copy(add_name='_oceanC', only=proc_oceanC)

## land carbon cycle only
proc_landC = ['cveg_0', 'csoil1_0', 'csoil2_0', 'f_fert', 'D_npp', 'f_igni', 'D_efire', 'D_eharv', 'D_egraz', 'D_fmort1', 'D_fmort2', 'f_resp', 'D_rh1', 'D_fmet', 'D_rh2', 'D_nbp', 'D_cveg', 'D_csoil1', 'D_csoil2']
proc_landC += ['D_Fveg_bk', 'D_Fsoil1_bk', 'D_Fsoil2_bk', 'D_Fslash', 'D_Fslash1', 'D_Fslash2', 'D_Fhwp', 'D_NPP_bk', 'D_Efire_bk', 'D_Eharv_bk', 'D_Egraz_bk', 'D_Fmort1_bk', 'D_Fmort2_bk', 'D_Rh1_bk', 'D_Fmet_bk', 'D_Rh2_bk', 'D_Ehwp', 'D_NBP_bk', 'D_Eluc', 'D_Fland', 'D_Aland', 'D_Cveg_bk', 'D_Csoil1_bk', 'D_Csoil2_bk', 'D_Chwp']
OSCAR_landC = OSCAR.copy(add_name='_landC', only=proc_landC)

