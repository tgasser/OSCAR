## Copyright CEA (Commissariat a L'Energie Atomique) & UVSQ (Universite de Versailles et Saint-Quentin), (2016)
## Contributor(s): Thomas Gasser (tgasser@lsce.ipsl.fr)

## This software is a computer program whose purpose is to simulate the behavior of the Earth system, with a specific but not exclusive focus on anthropogenic climate change.

## This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

## As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

## In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

## The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.

##################################################
##################################################
##################################################


import os
import csv

import numpy as np

print 'LOADING: DRIVERS'


##################################################
#   A. VECTORS
##################################################

# ============
# A.1. Regions
# ============

for X in ['I','J']:

    exec('mod_region = mod_region'+X)
    
    if (mod_region == 'SRES4'):
        region       = ['OECD90','REF','ASIA','ALM']
        region_name  = region
        region_color = ['#0000FF','#660099','#006600','#FF6600']
    elif (mod_region == 'SRES11'):
        region       = ['NAM','WEU','PAO','EEU','FSU','CPA','SAS','PAS','MEA','LAM','AFR']
        region_name  = ['North America','Western Europe','Pacific OECD','Central and Eastern Europe','Former Soviet Union','Centrally Planned Asia','South Asia','Pacific Asia','Middle East and North Africa','Latin America and the Caribbean','Sub-Saharan Africa']
        region_color = ['#000099','#0000FF','#00FFFF','#660099','#FF0099','#006600','#00FF00','#33FF33','#FFFF00','#FF6600','#FF0000']
    elif (mod_region == 'RECCAP*'):
        region       = ['L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','LX']
        region_name  = ['Africa','Arctic','Oceania','China','S.E. Asia','S. Asia','Europe','N. America','Russia','S. America','Rest']
        region_color = ['','','','','','','','','','','']
    elif (mod_region == 'Raupach*'):
        region       = ['USA','EU','Japan','D1','FSU','China','India','D2','D3']
        region_name  = region
        region_color = ['#FF0000','#FF6600','#FFFF00','#00FF00','#00FFFF','#0000FF','#FF0099','#660066','#000000']
    elif (mod_region == 'Houghton'):
        region       = ['N. Am.','S. & C. Am.','Europe','N. Afr. & M. East','Trop. Afr.','FSU','China','S. & S.E. Asia','Pacific Dvp.']
        region_name  = ['North America','South & Central America','Europe','North Africa & Middle East','Tropical Africa','Former Soviet Union','China region','South & South-East Asia','Pacific Developed region']
        region_color = ['#FF0000','#FFFF33','#00FF00','#00FFFF','#0000FF','#FF00FF','#009999','#FF6600','#990066']
    elif (mod_region == 'IMACLIM'):
        region       = ['USA','Canada','Europe','OECD Pacific','CEI','China','India','Brazil','Middle East','Africa','Rest of Asia','Rest of LAM']
        region_name  = region
        region_color = ['','','','','','','','','','','','']
    elif (mod_region == 'Kyoto'):
        region       = ['B','nB']
        region_name  = ['Annex B','non- Annex B']
        region_color = ['#006600','#CCCC99']
    elif (mod_region == 'RCP5'):
        region       = ['ASIA','LAM','MAF','OECD90','REF']
        region_name  = ['Asia region','Latin America','Middle-East & Africa','OECD countries in 1990','Reforming countries']
        region_color = ['','','','','']
    elif (mod_region == 'RCP10*'):
        region       = ['China +','India +','Rest of Asia','Latin America','Middle East','Africa','Western Europe','Northern America','Pacific OECD','Reforming Economies']
        region_name  = region
        region_color = ['','','','','','','','','','']
    else:
        region       = []
        region_name  = []
        region_color = []

    region_index = {}
    TMP = np.array([line for line in csv.reader(open('data/Regions_GTAP/#DATA.Regions_GTAP.csv','r'))])[:,2:]
    for n in range(1,len(TMP)):
        if mod_region in list(TMP[0]):
            region_index[int(TMP[n,0])] = int(TMP[n,list(TMP[0]).index(mod_region)])
        else:
            region_index[int(TMP[n,0])] = 0

    exec('region'+X+' = ["n/a"]+region')
    exec('region'+X+'_name = ["n/a"]+region_name')
    exec('region'+X+'_color = ["0.5"]+region_color')
    exec('region'+X+'_index = dict([(0,0)]+region_index.items())')
    exec('nb_region'+X+' = len(region'+X+')')

    del region,region_name,region_color,region_index

# ============
# A.2. Sectors
# ============

if (mod_sector == 'Time')&(300 < ind_final <= 310):
    sector       = ['<1850','1850-1899','1900-1909','1910-1919','1920-1929','1930-1939','1940-1949','1950-1959']\
                 + [str(t)+'-'+str(t+4) for t in range(1960,1990,5)]\
                 + [str(t)+'-'+str(t+1) for t in range(1990,2000,2)]\
                 + [str(t) for t in range(2000,1700+ind_final+1)]
    sector_name  = sector
    sector_color = ['0.5' for n in range(len(sector))]
elif (mod_sector == 'TimeRCP')&(ind_final == 400):
    sector       = ['<2011']+[str(t+1)+'-'+str(t+5) for t in range(2010,2100,5)]
    sector_name  = sector
    sector_color = ['0.5' for n in range(len(sector))]
else:
    sector       = ['n/a']
    sector_name  = sector
    sector_color = ['0.5']

nb_sector = len(sector)

# ==========
# A.3. Kinds
# ==========

kind       = ['n/a']
kind_name  = ['n/a']
kind_color = ['0.5']

# fossil CO2
if (mod_kindFF == 'one'):
    kind       += ['FF']
    kind_name  += ['Fossil Fuel']
    kind_color += ['#FF0000'] #or ['#666666']
    kFF  = 1
    kLUC = kFF+1
elif (mod_kindFF == 'CDIAC'):
    kind       += ['FF Solids','FF Liquids','FF Gas','FF Cement','FF Flaring']
    kind_name  += ['FF Solids','FF Liquids','FF Gas','FF Cement','FF Flaring']
    kind_color += ['','','','','']
    kFF  = 1
    kLUC = kFF+5
else:
    kFF = 0
    kLUC = max(0)+1

if (mod_kindFF == 'CDIAC'):
    kindFF_index = {'sol':kFF,'liq':kFF+1,'gas':kFF+2,'cem':kFF+3,'fla':kFF+4}
else:
    kindFF_index = {'sol':kFF,'liq':kFF,'gas':kFF,'cem':kFF,'fla':kFF}

# land-use
if (mod_kindLUC == 'one'):
    kind       += ['LULCC']
    kind_name  += ['Land-Use and Land-Cover Change']
    kind_color += ['#993300']
    kGHG = kLUC+1
elif (mod_kindLUC == 'all'):
    kind       += ['LUC-CO2','LUC-BB']
    kind_name  += ['LUC CO2 only','LUC BB non-CO2']
    kind_color += []
    kGHG = kLUC+2
else:
    kLUC = 0
    kGHG = max(kFF)+1

if (mod_kindLUC == 'all'):
    kindLUC_index = {'CO2':kLUC,'BB':kLUC+1}
else:
    kindLUC_index = {'CO2':kLUC,'BB':kLUC}  

# other GHG
if (mod_kindGHG == 'one'):
    kind       += ['non-CO2']
    kind_name  += ['non-CO2']
    kind_color += ['#FFCC00']
    kCHI = kGHG+1
elif (mod_kindGHG == 'RCP'):
    kind       += ['CH4','N2O','HaloC']
    kind_name  += ['Methane','Nitrous Oxide','Halocarbons']
    kind_color += ['#FF6600','#FFCC00','#FF9999']
    kCHI = kGHG+3
else:
    kGHG = 0
    kCHI = max(kFF,kLUC)+1

if (mod_kindGHG == 'RCP'):
    kindGHG_index = {'CH4':kGHG,'N2O':kGHG+1,'HFC':kGHG+2,'PFC':kGHG+2,'ODS':kGHG+2}
else:
    kindGHG_index = {'CH4':kGHG,'N2O':kGHG,'HFC':kGHG,'PFC':kGHG,'ODS':kGHG}

# active species
if (mod_kindCHI == 'one'):
    kind       += ['OzPrec.']
    kind_name  += ['Ozone Precursors']
    kind_color += ['#66FF66']
    kAER = kCHI+1
elif (mod_kindCHI == 'all'):
    kind       += ['NOx','CO','NMVOC']
    kind_name  += ['Nitrogen Oxides','Carbon Monoxide','Non-Methane Volatile Organic Compounds']
    kind_color += ['#66FF66','#009999','#006600']
    kAER = kCHI+3
else:
    kCHI = 0
    kAER = max(kFF,kLUC,kGHG)+1

if (mod_kindCHI == 'all'):
    kindCHI_index = {'NOX':kCHI,'CO':kCHI+1,'VOC':kCHI+2}
else:
    kindCHI_index = {'NOX':kCHI,'CO':kCHI,'VOC':kCHI}    

# aerosols
if (mod_kindAER == 'one'):
    kind       += ['AER']
    kind_name  += ['Aerosols']
    kind_color += ['#0000FF']
    kRF = kAER+1
elif (mod_kindAER == 'all'):
    kind       += ['SO2','NH3','OC','BC']
    kind_name  += ['Sulfur Dioxide','Ammonia','Organic Carbon','Black Carbon']
    kind_color += ['#00CCFF','#0000FF','#660099','#CC0066']
    kRF = kAER+4
else:
    kAER = 0
    kRF = max(kFF,kLUC,kGHG,kCHI)+1

if (mod_kindAER == 'all'):
    kindAER_index = {'SO2':kAER,'NH3':kAER+1,'OC':kAER+2,'BC':kAER+3}
else:
    kindAER_index = {'SO2':kAER,'NH3':kAER,'OC':kAER,'BC':kAER}

# other radiative forcings
if (mod_kindRF == 'one'):
    kind       += ['RFother']
    kind_name  += ['Other RF']
    kind_color += ['#999999']
    kGE = kRF+1
elif (mod_kindRF == 'two'):
    kind       += ['RFant','RFnat']
    kind_name  += ['Anthropogenic RF','Natural RF']
    kind_color += ['','']
    kGE = kRF+2
elif (mod_kindRF == 'all'):
    kind       += ['RFcon','RFsol','RFvol']
    kind_name  += ['Contrails RF','Solar RF','Volcanoes RF']
    kind_color += ['','','']
    kGE = kRF+3    
else:
    kRF = 0
    kGE = max(kFF,kLUC,kGHG,kCHI,kAER)+1

if (mod_kindRF == 'all'):
    kindRF_index = {'RFcon':kRF,'RFsol':kRF+1,'RFvol':kRF+2}
elif (mod_kindRF == 'two'):
    kindRF_index = {'RFcon':kRF,'RFsol':kRF+1,'RFvol':kRF+1}
else:
    kindRF_index = {'RFcon':kRF,'RFsol':kRF,'RFvol':kRF}

# geoengineering
if (mod_kindGE == 'PUP'):
    kind       += ['AFO','CCS','ALB','AER']
    kind_name  += ['Aforestation','Carbon Capture and Storage','Surface Albedo','Sulfate Aerosols']
    kind_color += ['','','','']
else:
    kGE = 0

nb_kind = len(kind)

# ===========
# A.4. Biomes
# ===========

if (mod_biomeSHR == 'w/FOR')&(mod_biomeURB == 'w/DES'):
    biome = ['DES+','FOR+','GRA','CRO','PAS']
    biome_name = ['Desert & Urban','Forest & Shrubland','Grassland','Cropland','Pasture']
    biome_color = ['','','','','']
    biome_index = {'des':0,'for':1,'shr':1,'gra':2,'cro':3,'pas':4,'urb':0}
elif (mod_biomeSHR == 'w/FOR')&(mod_biomeURB == 'URB'):
    biome = ['DES','FOR+','GRA','CRO','PAS','URB']
    biome_name = ['Desert','Forest & Shrubland','Grassland','Cropland','Pasture','Urban']
    biome_color = ['','','','','','']
    biome_index = {'des':0,'for':1,'shr':1,'gra':2,'cro':3,'pas':4,'urb':5}
elif (mod_biomeSHR == 'w/GRA')&(mod_biomeURB == 'w/DES'):
    biome = ['DES+','FOR','GRA+','CRO','PAS']
    biome_name = ['Desert & Urban','Forest','Grassland & Shrubland','Cropland','Pasture']
    biome_color = ['','','','','']
    biome_index = {'des':0,'for':1,'shr':2,'gra':2,'cro':3,'pas':4,'urb':0}
elif (mod_biomeSHR == 'w/GRA')&(mod_biomeURB == 'URB'):
    biome = ['DES','FOR','GRA+','CRO','PAS','URB']
    biome_name = ['Desert','Forest','Grassland & Shrubland','Cropland','Pasture','Urban']
    biome_color = ['','','','','','']
    biome_index = {'des':0,'for':1,'shr':2,'gra':2,'cro':3,'pas':4,'urb':5}
elif (mod_biomeSHR == 'SHR')&(mod_biomeURB == 'w/DES'):
    biome = ['DES+','FOR','SHR','GRA','CRO','PAS']
    biome_name = ['Desert & Urban','Forest','Shrubland','Grassland','Cropland','Pasture']
    biome_color = ['','','','','','']
    biome_index = {'des':0,'for':1,'shr':2,'gra':3,'cro':4,'pas':5,'urb':0}
elif (mod_biomeSHR == 'SHR')&(mod_biomeURB == 'URB'):
    biome = ['DES','FOR','SHR','GRA','CRO','PAS','URB']
    biome_name = ['Desert','Forest','Shrubland','Grassland','Cropland','Pasture','Urban']
    biome_color = ['','','','','','','']
    biome_index = {'des':0,'for':1,'shr':2,'gra':3,'cro':4,'pas':5,'urb':6}
else:
    biome = ['all']
    biome_name = biome
    biome_color = []
    biome_index = {'des':0,'for':0,'shr':0,'gra':0,'cro':0,'pas':0,'urb':0}

nb_biome = len(biome)

# =========
# A.5. Halo
# =========

HFC = ['HFC23','HFC32','HFC125','HFC134a','HFC143a','HFC152a','HFC227ea','HFC236fa','HFC245fa','HFC365mfc','HFC4310mee']
PFC = ['SF6','NF3','CF4','C2F6','C3F8','cC4F8','C4F10','C5F12','C6F14','C7F16']
ODS = ['CFC11','CFC12','CFC113','CFC114','CFC115','CCl4','CH3CCl3','HCFC22','HCFC141b','HCFC142b','Halon1211','Halon1202','Halon1301','Halon2402','CH3Br','CH3Cl']

nb_HFC = len(HFC)
nb_PFC = len(PFC)
nb_ODS = len(ODS)


##################################################
#   1. GREENHOUSE GASES
##################################################

EFF = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty) # {GtC/yr}
ECH4 = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty) # {TgC/yr}
EN2O = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty) # {TgN/yr}

# ==========
# 1.1. CDIAC
# ==========

# load CDIAC emissions
# from [Boden et al., 2012]
ind_cdiac = 310
kin = ['sol','liq','gas','cem','fla']

# total emissions
EFFcdiac = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for k in range(len(kin)):
    TMP = np.array([line for line in csv.reader(open('data/EFossil_CDIAC/#DATA.EFossil_CDIAC.1751-2010_114reg0.EFF_'+kin[k]+'.csv','r'))], dtype=dty)
    for i in range(114+1):
        EFFcdiac[51:ind_cdiac+1,regionJ_index[i],0,kFF,regionI_index[i]] += TMP[:ind_cdiac-51+1,i]

# distribution among fuels
# TODO

# ========
# 1.2. EPA
# ========

# load EPA emissions
# see [EPA, 2012]
ind_epa = 310
sec_epa = ['agr_ferm','agr_manu','agr_othr','agr_rice','agr_soil','ene_burn','ene_coal','ene_comb','ene_ngos','ene_othr','ind_acid']
sec_epa1 = ['agr_ferm','agr_manu','agr_othr','agr_rice','agr_soil','ene_coal','ene_comb','ene_ngos','ene_othr','ind_acid']

# CH4
ECH4epa = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for s in range(len(sec_epa1)):
    if os.path.isfile('data/EMethane_EPA/#DATA.EMethane_EPA.1990-'+str(1700+ind_epa)+'_114reg0.ECH4_'+sec_epa1[s]+'.csv'):
        TMP = np.array([line for line in csv.reader(open('data/EMethane_EPA/#DATA.EMethane_EPA.1990-'+str(1700+ind_epa)+'_114reg0.ECH4_'+sec_epa1[s]+'.csv','r'))], dtype=dty)
        for i in range(114+1):
            ECH4epa[290:ind_epa+1,regionJ_index[i],0,kindGHG_index['CH4'],regionI_index[i]] += TMP[:ind_epa-290+1,i]

# N2O
EN2Oepa = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for s in range(len(sec_epa1)):
    if os.path.isfile('data/ENitrousOx_EPA/#DATA.ENitrousOx_EPA.1990-'+str(1700+ind_epa)+'_114reg0.EN2O_'+sec_epa1[s]+'.csv'):
        TMP = np.array([line for line in csv.reader(open('data/ENitrousOx_EPA/#DATA.ENitrousOx_EPA.1990-'+str(1700+ind_epa)+'_114reg0.EN2O_'+sec_epa1[s]+'.csv','r'))], dtype=dty)
        for i in range(114+1):
            EN2Oepa[290:ind_epa+1,regionJ_index[i],0,kindGHG_index['N2O'],regionI_index[i]] += TMP[:ind_epa-290+1,i]

# ==========
# 1.3. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]
ind_edgar = 308
sec_ehyde = ['oo','fc','fp','bc','in','al','an','aw','lf','sb','df']
sec_accmip = ['ooo','ene','ind','tra','shp','air','dom','slv','agr','awb','wst','for','gra']

# FF
EFFedgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for s in range(1,len(sec_ehyde)-2):
    if os.path.isfile('data/ECarbon_EDGAR/#DATA.ECarbon_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.ECO2_'+sec_ehyde[s]+'.csv'):
        TMP = np.array([line for line in csv.reader(open('data/ECarbon_EDGAR/#DATA.ECarbon_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.ECO2_'+sec_ehyde[s]+'.csv','r'))], dtype=dty)
        for i in range(114+1):
            EFFedgar[270:ind_edgar+1,regionJ_index[i],0,kFF,regionI_index[i]] += TMP[:ind_edgar-270+1,i]

# CH4
ECH4edgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for s in range(1,len(sec_accmip)-2):
    if os.path.isfile('data/EMethane_EDGAR/#DATA.EMethane_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.ECH4_'+sec_accmip[s]+'.csv'):
        TMP = np.array([line for line in csv.reader(open('data/EMethane_EDGAR/#DATA.EMethane_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.ECH4_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
        for i in range(114+1):
            ECH4edgar[270:ind_edgar+1,regionJ_index[i],0,kindGHG_index['CH4'],regionI_index[i]] += TMP[:ind_edgar-270+1,i]

# N2O
EN2Oedgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for s in range(1,len(sec_ehyde)-2):
    if os.path.isfile('data/ENitrousOx_EDGAR/#DATA.ENitrousOx_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.EN2O_'+sec_ehyde[s]+'.csv'):
        TMP = np.array([line for line in csv.reader(open('data/ENitrousOx_EDGAR/#DATA.ENitrousOx_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.EN2O_'+sec_ehyde[s]+'.csv','r'))], dtype=dty)
        for i in range(114+1):
            EN2Oedgar[270:ind_edgar+1,regionJ_index[i],0,kindGHG_index['N2O'],regionI_index[i]] += TMP[:ind_edgar-270+1,i]

# =============
# 1.4. EDGAR-FT
# =============

# load emissions from EDGAR-FT v4.2-FT2010
# see [JRC, 2013]

# FF
EFFeft = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for s in range(1,len(sec_ehyde)-2):
    if os.path.isfile('data/ECarbon_EDGAR-FT/#DATA.ECarbon_EDGAR-FT.2008-2010_114reg0.ECO2_'+sec_ehyde[s]+'.csv'):
        TMP = np.array([line for line in csv.reader(open('data/ECarbon_EDGAR-FT/#DATA.ECarbon_EDGAR-FT.2008-2010_114reg0.ECO2_'+sec_ehyde[s]+'.csv','r'))], dtype=dty)
        for i in range(114+1):
            EFFeft[308:310+1,regionJ_index[i],0,kFF,regionI_index[i]] += TMP[:310-308+1,i]

# CH4
ECH4eft = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for s in range(1,len(sec_accmip)-2):
    if os.path.isfile('data/EMethane_EDGAR-FT/#DATA.EMethane_EDGAR-FT.2008-2010_114reg0.ECH4_'+sec_accmip[s]+'.csv'):
        TMP = np.array([line for line in csv.reader(open('data/EMethane_EDGAR-FT/#DATA.EMethane_EDGAR-FT.2008-2010_114reg0.ECH4_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
        for i in range(114+1):
            ECH4eft[308:310+1,regionJ_index[i],0,kindGHG_index['CH4'],regionI_index[i]] += TMP[:310-308+1,i]

# N2O
EN2Oeft = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for s in range(1,len(sec_ehyde)-2):
    if os.path.isfile('data/ENitrousOx_EDGAR-FT/#DATA.ENitrousOx_EDGAR-FT.2008-2010_114reg0.EN2O_'+sec_ehyde[s]+'.csv'):
        TMP = np.array([line for line in csv.reader(open('data/ENitrousOx_EDGAR-FT/#DATA.ENitrousOx_EDGAR-FT.2008-2010_114reg0.EN2O_'+sec_ehyde[s]+'.csv','r'))], dtype=dty)
        for i in range(114+1):
            EN2Oeft[308:310+1,regionJ_index[i],0,kindGHG_index['N2O'],regionI_index[i]] += TMP[:310-308+1,i]

# ===========
# 1.5. ACCMIP
# ===========

# load emissions from ACCMIP
# from [Lamarque et al., 2010]

# CH4
ECH4accmip = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
p_ECH4_bio = np.zeros([nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for s in range(1,len(sec_accmip)-2):
    if os.path.isfile('data/EMethane_ACCMIP/#DATA.EMethane_ACCMIP.1850-2000_114reg0.ECH4_'+sec_accmip[s]+'.csv'):
        TMP = np.array([line for line in csv.reader(open('data/EMethane_ACCMIP/#DATA.EMethane_ACCMIP.1850-2000_114reg0.ECH4_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
        for i in range(114+1):
            ECH4accmip[150:300+1,regionJ_index[i],0,kindGHG_index['CH4'],regionI_index[i]] += TMP[:300-150+1,i]
            if sec_accmip[s] in ['agr','awb','wst']:
                p_ECH4_bio[regionJ_index[i],0,kindGHG_index['CH4'],regionI_index[i]] += TMP[0,i]
p_ECH4_bio /= ECH4accmip[150]
p_ECH4_bio[np.isnan(p_ECH4_bio)|np.isinf(p_ECH4_bio)] = 0

# ===============
# 1.6. EDGAR-HYDE
# ===============

# load emissions from EDGAR-HYDE v1.3
# from [van Aardenne et al., 2001]

# N2O
EN2Oehyde = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
p_EN2O_bio = np.zeros([nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for s in range(1,len(sec_ehyde)-2):
    if os.path.isfile('data/ENitrousOx_EDGAR-HYDE/#DATA.ENitrousOx_EDGAR-HYDE.1890-1990_114reg0.EN2O_'+sec_ehyde[s]+'.csv'):
        TMP = np.array([line for line in csv.reader(open('data/ENitrousOx_EDGAR-HYDE/#DATA.ENitrousOx_EDGAR-HYDE.1890-1990_114reg0.EN2O_'+sec_ehyde[s]+'.csv','r'))], dtype=dty)
        for i in range(114+1):
            EN2Oehyde[190:290+1,regionJ_index[i],0,kindGHG_index['N2O'],regionI_index[i]] += TMP[:290-190+1,i]
            if sec_ehyde[s] in ['al','an','aw','lf']:
                p_EN2O_bio[regionJ_index[i],0,kindGHG_index['N2O'],regionI_index[i]] += TMP[0,i]
p_EN2O_bio /= EN2Oehyde[190]
p_EN2O_bio[np.isnan(p_EN2O_bio)|np.isinf(p_EN2O_bio)] = 0

# ==============
# 1.7. Stern1998
# ==============

# load emissions from [Stern et al., 1998]
ECH4stern = np.zeros([ind_cdiac+1], dtype=dty)
TMP = np.array([line for line in csv.reader(open('data/EMethane_Stern1998/#DATA.EMethane_Stern1998.1860-1994_(7sec).ECH4.csv','r'))][1:], dtype=dty)
lgd = [line for line in csv.reader(open('data/EMethane_Stern1998/#DATA.EMethane_Stern1998.1860-1994_(7sec).ECH4.csv','r'))][0]
for s in range(len(lgd)):
    if not lgd[s] in ['Biomass Burning']:
        ECH4stern[160:294+1] += TMP[:294-160+1,s]

# =================
# 1.8. Davidson2009
# =================

# load emissions from [Davidson et al., 2009]
EN2Odavidson = np.zeros([ind_cdiac+1], dtype=dty)
TMP = np.array([line for line in csv.reader(open('data/ENitrousOx_Davidson2009/#DATA.ENitrousOx_Davidson2009.1860-2005_(5sec).EN2O.csv','r'))][1:], dtype=dty)
lgd = [line for line in csv.reader(open('data/ENitrousOx_Davidson2009/#DATA.ENitrousOx_Davidson2009.1860-2005_(5sec).EN2O.csv','r'))][0]
for s in range(len(lgd)):
    if lgd[s] in ['nyl_prod','ff_burn','biogen']:
        EN2Odavidson[160:305+1] += TMP[:305-160+1,s]

# global rescaling of EDGAR-HYDE emissions
EN2OehydeR = EN2Oehyde.copy()
if True:
    EN2OehydeR[190:290+1,...] *= (EN2Odavidson[190:290+1]/np.sum(np.sum(np.sum(np.sum(EN2Oehyde[190:290+1,...],4),3),2),1))[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]

# =========
# 1.9. SRES
# =========

# initialization of projected drivers
for VAR in ['FF','CH4','N2O']:
    exec('E'+VAR+'proj = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')

# projection of emissions under SRES scenarios
# from [IPCC, 2000]

# FF
if (scen_EFF[:4] == 'SRES')&(ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open('data/EFossil_SRES/#DATA.EFossil_SRES.2000-2100_4reg0.'+scen_EFF[5:]+'_EFF.csv','r'))], dtype=dty)
    for i in range(4+1):
        if (mod_regionI == 'SRES4')&(mod_regionJ == 'SRES4'):
            EFFproj[300:min(ind_final,400)+1,i,0,kFF,i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI == 'SRES4')&(mod_regionJ != 'SRES4'):
            EFFproj[300:min(ind_final,400)+1,0,0,kFF,i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI != 'SRES4')&(mod_regionJ == 'SRES4'):
            EFFproj[300:min(ind_final,400)+1,i,0,kFF,0] += TMP[:min(ind_final,400)-300+1,i]
        elif(mod_regionI != 'SRES4')&(mod_regionJ != 'SRES4'):
            EFFproj[300:min(ind_final,400)+1,0,0,kFF,0] += TMP[:min(ind_final,400)-300+1,i]

# CH4
if (scen_ECH4[:4] == 'SRES')&(ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open('data/EMethane_SRES/#DATA.EMethane_SRES.2000-2100_4reg0.'+scen_ECH4[5:]+'_ECH4.csv','r'))], dtype=dty)
    for i in range(4+1):
        if (mod_regionI == 'SRES4')&(mod_regionJ == 'SRES4'):
            ECH4proj[300:min(ind_final,400)+1,i,0,kindGHG_index['CH4'],i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI == 'SRES4')&(mod_regionJ != 'SRES4'):
            ECH4proj[300:min(ind_final,400)+1,0,0,kindGHG_index['CH4'],i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI != 'SRES4')&(mod_regionJ == 'SRES4'):
            ECH4proj[300:min(ind_final,400)+1,i,0,kindGHG_index['CH4'],0] += TMP[:min(ind_final,400)-300+1,i]
        elif(mod_regionI != 'SRES4')&(mod_regionJ != 'SRES4'):
            ECH4proj[300:min(ind_final,400)+1,0,0,kindGHG_index['CH4'],0] += TMP[:min(ind_final,400)-300+1,i]

# N2O
if (scen_EN2O[:4] == 'SRES')&(ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open('data/ENitrousOx_SRES/#DATA.ENitrousOx_SRES.2000-2100_4reg0.'+scen_EN2O[5:]+'_EN2O.csv','r'))], dtype=dty)
    for i in range(4+1):
        if (mod_regionI == 'SRES4')&(mod_regionJ == 'SRES4'):
            EN2Oproj[300:min(ind_final,400)+1,i,0,kindGHG_index['N2O'],i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI == 'SRES4')&(mod_regionJ != 'SRES4'):
            EN2Oproj[300:min(ind_final,400)+1,0,0,kindGHG_index['N2O'],i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI != 'SRES4')&(mod_regionJ == 'SRES4'):
            EN2Oproj[300:min(ind_final,400)+1,i,0,kindGHG_index['N2O'],0] += TMP[:min(ind_final,400)-300+1,i]
        elif(mod_regionI != 'SRES4')&(mod_regionJ != 'SRES4'):
            EN2Oproj[300:min(ind_final,400)+1,0,0,kindGHG_index['N2O'],0] += TMP[:min(ind_final,400)-300+1,i]

# =========
# 1.10. RCP
# =========

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# FF
if (scen_EFF[:3] == 'RCP')&(ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open('data/EFossil_RCP/#DATA.EFossil_RCP.2000-2100_5reg0.rcp'+scen_EFF[3]+scen_EFF[5]+'_EFF.csv','r'))], dtype=dty)
    for i in range(5+1):
        if (mod_regionI == 'RCP5')&(mod_regionJ == 'RCP5'):
            EFFproj[300:min(ind_final,400)+1,i,0,kFF,i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI == 'RCP5')&(mod_regionJ != 'RCP5'):
            EFFproj[300:min(ind_final,400)+1,0,0,kFF,i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI != 'RCP5')&(mod_regionJ == 'RCP5'):
            EFFproj[300:min(ind_final,400)+1,i,0,kFF,0] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI != 'RCP5')&(mod_regionJ != 'RCP5'):
            EFFproj[300:min(ind_final,400)+1,0,0,kFF,0] += TMP[:min(ind_final,400)-300+1,i]

# CH4
if (scen_ECH4[:3] == 'RCP')&(ind_final > ind_cdiac):
    for s in range(1,len(sec_accmip)-2):
        if os.path.isfile('data/EMethane_RCP/#DATA.EMethane_RCP.2000-2100_5reg0.rcp'+scen_ECH4[3]+scen_ECH4[5]+'_ECH4_'+sec_accmip[s]+'.csv'):
            TMP = np.array([line for line in csv.reader(open('data/EMethane_RCP/#DATA.EMethane_RCP.2000-2100_5reg0.rcp'+scen_ECH4[3]+scen_ECH4[5]+'_ECH4_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
            for i in range(5+1):
                if (mod_regionI == 'RCP5')&(mod_regionJ == 'RCP5'):
                    ECH4proj[300:min(ind_final,400)+1,i,0,kindGHG_index['CH4'],i] += TMP[:min(ind_final,400)-300+1,i]
                elif (mod_regionI == 'RCP5')&(mod_regionJ != 'RCP5'):
                    ECH4proj[300:min(ind_final,400)+1,0,0,kindGHG_index['CH4'],i] += TMP[:min(ind_final,400)-300+1,i]
                elif (mod_regionI != 'RCP5')&(mod_regionJ == 'RCP5'):
                    ECH4proj[300:min(ind_final,400)+1,i,0,kindGHG_index['CH4'],0] += TMP[:min(ind_final,400)-300+1,i]
                elif (mod_regionI != 'RCP5')&(mod_regionJ != 'RCP5'):
                    ECH4proj[300:min(ind_final,400)+1,0,0,kindGHG_index['CH4'],0] += TMP[:min(ind_final,400)-300+1,i]

# N2O
if (scen_EN2O[:3] == 'RCP')&(ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open('data/ENitrousOx_RCP/#DATA.ENitrousOx_RCP.2000-2100_5reg0.rcp'+scen_EN2O[3]+scen_EN2O[5]+'_EN2O.csv','r'))], dtype=dty)
    for i in range(5+1):
        if (mod_regionI == 'RCP5')&(mod_regionJ == 'RCP5'):
            EN2Oproj[300:min(ind_final,400)+1,i,0,kindGHG_index['N2O'],i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI == 'RCP5')&(mod_regionJ != 'RCP5'):
            EN2Oproj[300:min(ind_final,400)+1,0,0,kindGHG_index['N2O'],i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI != 'RCP5')&(mod_regionJ == 'RCP5'):
            EN2Oproj[300:min(ind_final,400)+1,i,0,kindGHG_index['N2O'],0] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI != 'RCP5')&(mod_regionJ != 'RCP5'):
            EN2Oproj[300:min(ind_final,400)+1,0,0,kindGHG_index['N2O'],0] += TMP[:min(ind_final,400)-300+1,i]

# =================
# 1.A. PAST DATASET
# =================

# datasets mixed following trends
for VAR in ['FF','CH4','N2O']:

    if VAR in ['FF']:
        
        # with CDIAC as reference
        if (data_EFF == 'CDIAC'):
            exec('E'+VAR+'past = E'+VAR+'cdiac.copy()')

        # with EDGAR as reference
        elif (data_EFF == 'EDGAR'):
            exec('E'+VAR+'past = E'+VAR+'edgar.copy()')
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar+1,ind_cdiac+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * E'+VAR+'eft[t,...]/E'+VAR+'eft[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'eft[t,...])/np.sum(E'+VAR+'eft[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')        
            # follow CDIAC variations before 1970
            for t in range(50,270)[::-1]:
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t+1,...] * E'+VAR+'cdiac[t,...]/E'+VAR+'cdiac[t+1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t+1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'cdiac[t,...])/np.sum(E'+VAR+'cdiac[t+1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')

    elif VAR in ['CH4']:

        # with EDGAR as reference
        if (data_ECH4 == 'EDGAR'):
            exec('E'+VAR+'past = E'+VAR+'edgar.copy()')
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar+1,ind_cdiac+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * E'+VAR+'eft[t,...]/E'+VAR+'eft[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'eft[t,...])/np.sum(E'+VAR+'eft[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')        
            # follow ACCMIP variations before 1970
            for t in range(150,270)[::-1]:
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t+1,...] * E'+VAR+'accmip[t,...]/E'+VAR+'accmip[t+1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t+1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'accmip[t,...])/np.sum(E'+VAR+'accmip[t+1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')
            # linear extrapolation before 1850
            for t in range(50,150):
                exec('E'+VAR+'past[t,...] = (1-p_E'+VAR+'_bio) * E'+VAR+'past[150,...] * (t-50)/float(150-50)')
            for t in range(0,150):
                exec('E'+VAR+'past[t,...] += p_E'+VAR+'_bio * E'+VAR+'past[150,...] * (t--200)/float(150--200)')            
            exec('E'+VAR+'_0 = E'+VAR+'past[50,...]')

        # with ACCMIP as reference
        elif (data_ECH4 == 'ACCMIP'):
            exec('E'+VAR+'past = E'+VAR+'accmip.copy()')
            # follow EDGAR variations after 2000
            for t in range(300+1,ind_edgar+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * E'+VAR+'edgar[t,...]/E'+VAR+'edgar[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'edgar[t,...])/np.sum(E'+VAR+'edgar[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')        
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar+1,ind_cdiac+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * E'+VAR+'eft[t,...]/E'+VAR+'eft[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'eft[t,...])/np.sum(E'+VAR+'eft[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')        
            # linear extrapolation before 1850
            for t in range(50,150):
                exec('E'+VAR+'past[t,...] = (1-p_E'+VAR+'_bio) * E'+VAR+'past[150,...] * (t-50)/float(150-50)')
            for t in range(0,150):
                exec('E'+VAR+'past[t,...] += p_E'+VAR+'_bio * E'+VAR+'past[150,...] * (t--200)/float(150--200)')            
            exec('E'+VAR+'_0 = E'+VAR+'past[50,...]')

        # with EPA as reference
        elif (data_ECH4 == 'EPA'):
            exec('E'+VAR+'past = E'+VAR+'epa.copy()')
            # follow ACCMIP variations before 1990
            for t in range(150,290)[::-1]:
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t+1,...] * E'+VAR+'accmip[t,...]/E'+VAR+'accmip[t+1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t+1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'accmip[t,...])/np.sum(E'+VAR+'accmip[t+1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')
            # linear extrapolation before 1850
            for t in range(50,150):
                exec('E'+VAR+'past[t,...] = (1-p_E'+VAR+'_bio) * E'+VAR+'past[150,...] * (t-50)/float(150-50)')
            for t in range(0,150):
                exec('E'+VAR+'past[t,...] += p_E'+VAR+'_bio * E'+VAR+'past[150,...] * (t--200)/float(150--200)')            
            exec('E'+VAR+'_0 = E'+VAR+'past[50,...]')

    elif VAR in ['N2O']:

        # with EDGAR as reference
        if (data_EN2O == 'EDGAR'):
            exec('E'+VAR+'past = E'+VAR+'edgar.copy()')
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar+1,ind_cdiac+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * E'+VAR+'eft[t,...]/E'+VAR+'eft[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'eft[t,...])/np.sum(E'+VAR+'eft[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')        
            # follow EDGAR-HYDE variations before 1970
            for t in range(190,270)[::-1]:
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t+1,...] * E'+VAR+'ehydeR[t,...]/E'+VAR+'ehydeR[t+1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t+1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'ehydeR[t,...])/np.sum(E'+VAR+'ehydeR[t+1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')
            # linear extrapolation before 1890
            for t in range(50,190):
                exec('E'+VAR+'past[t,...] = (1-p_E'+VAR+'_bio) * E'+VAR+'past[190,...] * (t-50)/float(190-50)')
            for t in range(0,190):
                exec('E'+VAR+'past[t,...] += p_E'+VAR+'_bio * E'+VAR+'past[190,...] * (t--200)/float(190--200)')            
            exec('E'+VAR+'_0 = E'+VAR+'past[50,...]')

        # with EPA as reference
        elif (data_EN2O == 'EPA'):
            exec('E'+VAR+'past = E'+VAR+'epa.copy()')
            # follow EDGAR-HYDE variations before 1990
            for t in range(190,290)[::-1]:
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t+1,...] * E'+VAR+'ehydeR[t,...]/E'+VAR+'ehydeR[t+1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t+1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'ehydeR[t,...])/np.sum(E'+VAR+'ehydeR[t+1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')
            # linear extrapolation before 1850
            for t in range(50,190):
                exec('E'+VAR+'past[t,...] = (1-p_E'+VAR+'_bio) * E'+VAR+'past[190,...] * (t-50)/float(190-50)')
            for t in range(0,190):
                exec('E'+VAR+'past[t,...] += p_E'+VAR+'_bio * E'+VAR+'past[190,...] * (t--200)/float(190--200)')            
            exec('E'+VAR+'_0 = E'+VAR+'past[50,...]')

    # cut past dataset to right length
    exec('E'+VAR+'[:min(ind_cdiac,ind_final)+1,...] = E'+VAR+'past[:min(ind_cdiac,ind_final)+1,...]')

# ==================
# 1.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR in ['FF','CH4','N2O']:
    exec('scen = scen_E'+VAR)

    # stop emissions
    if (scen == 'stop')&(ind_final > ind_cdiac):
        if VAR in ['FF']:
            exec('E'+VAR+'[ind_cdiac+1:,...] = 0')
        else:
            exec('E'+VAR+'[ind_cdiac+1:,...] = E'+VAR+'_0[np.newaxis,...]')

    # constant emissions
    elif (scen == 'cst')&(ind_final > ind_cdiac):
        exec('E'+VAR+'[ind_cdiac+1:,...] = E'+VAR+'[ind_cdiac,...][np.newaxis,...]')        

    # RCP or SRES scenarios
    elif ((scen[:4] == 'SRES')|(scen[:3] == 'RCP'))&(ind_final > ind_cdiac):
        
        # raw discontinuity
        if (mod_DATAscen == 'raw'):
            exec('E'+VAR+'[ind_cdiac+1:,...] = E'+VAR+'proj[ind_cdiac+1:,...]')

        # offset at transition point
        elif (mod_DATAscen == 'offset'):
            exec('E'+VAR+'[ind_cdiac+1:,...] = E'+VAR+'proj[ind_cdiac+1:,...] - E'+VAR+'proj[ind_cdiac,...] + E'+VAR+'[ind_cdiac,...]')
            for t in range(ind_cdiac+1,ind_final+1):
                exec('def_regI = bool(np.sum(E'+VAR+'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E'+VAR+'proj[t,1:,...,:]))')
                if not def_regI:
                    exec('E'+VAR+'[t,:,...,0] += np.sum(E'+VAR+'[t,:,...,1:],-1)')
                    exec('E'+VAR+'[t,:,...,1:] = 0')
                if not def_regJ:
                    exec('E'+VAR+'[t,0,...,:] += np.sum(E'+VAR+'[t,1:,...,:],0)')
                    exec('E'+VAR+'[t,1:,...,:] = 0')                  

        # linear transition over N years
        elif (mod_DATAscen[:6] == 'smooth'):
            N = int(mod_DATAscen[6:])
            if (ind_final >= ind_cdiac+N):
                for t in range(ind_cdiac+1,ind_cdiac+N):
                    exec('E'+VAR+'[t,...] = (1-(t-ind_cdiac)/float(N)) * E'+VAR+'[ind_cdiac,...] + (t-ind_cdiac)/float(N) * E'+VAR+'proj[ind_cdiac+N,...]')
                    exec('def_regI = bool(np.sum(E'+VAR+'proj[t,:,...,1:]))')
                    exec('def_regJ = bool(np.sum(E'+VAR+'proj[t,1:,...,:]))')
                    if not def_regI:
                        exec('E'+VAR+'[t,:,...,0] += np.sum(E'+VAR+'[t,:,...,1:],-1)')
                        exec('E'+VAR+'[t,:,...,1:] = 0')
                    if not def_regJ:
                        exec('E'+VAR+'[t,0,...,:] += np.sum(E'+VAR+'[t,1:,...,:],0)')
                        exec('E'+VAR+'[t,1:,...,:] = 0')   
                exec('E'+VAR+'[ind_cdiac+N:,...] = E'+VAR+'proj[ind_cdiac+N:,...]')

        # follow trends of projection
        elif (mod_DATAscen == 'trends'):
            for t in range(ind_cdiac+1,ind_final+1):
                exec('def_regI = bool(np.sum(E'+VAR+'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E'+VAR+'proj[t,1:,...,:]))')
                if def_regI and def_regJ:
                    exec('E'+VAR+'[t,...] = E'+VAR+'[t-1,...] * E'+VAR+'proj[t,...]/E'+VAR+'proj[t-1,...]')
                    exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')
                    exec('E'+VAR+'[t,...] *= np.sum(E'+VAR+'[t-1,...])/np.sum(E'+VAR+'[t,...]) * np.sum(E'+VAR+'proj[t,...])/np.sum(E'+VAR+'proj[t-1,...])')
                elif not def_regI and def_regJ:
                    exec('E'+VAR+'[t,:,...,0] = np.sum(E'+VAR+'[t-1,:,...,:],-1) * np.sum(E'+VAR+'proj[t,:,...,:],-1)/np.sum(E'+VAR+'proj[t-1,:,...,:],-1)')
                    exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')
                    exec('E'+VAR+'[t,...] *= np.sum(E'+VAR+'[t-1,...])/np.sum(E'+VAR+'[t,...]) * np.sum(E'+VAR+'proj[t,...])/np.sum(E'+VAR+'proj[t-1,...])')
                elif def_regI and not def_regJ:
                    exec('E'+VAR+'[t,0,...,:] = np.sum(E'+VAR+'[t-1,:,...,:],0) * np.sum(E'+VAR+'proj[t,:,...,:],0)/np.sum(E'+VAR+'proj[t-1,:,...,:],0)')
                    exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')
                    exec('E'+VAR+'[t,...] *= np.sum(E'+VAR+'[t-1,...])/np.sum(E'+VAR+'[t,...]) * np.sum(E'+VAR+'proj[t,...])/np.sum(E'+VAR+'proj[t-1,...])')
                elif not def_regI and not def_regJ:
                    exec('E'+VAR+'[t,0,...,0] = np.sum(np.sum(E'+VAR+'[t-1,:,...,:],-1),0) * np.sum(np.sum(E'+VAR+'proj[t,:,...,:],-1),0)/np.sum(np.sum(E'+VAR+'proj[t-1,:,...,:],-1),0)')
            exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')

# delete individual datasets
for VAR in ['FF']:
    exec('del E'+VAR+'cdiac,E'+VAR+'edgar,E'+VAR+'eft,E'+VAR+'past,E'+VAR+'proj')
for VAR in ['CH4']:
    exec('del E'+VAR+'epa,E'+VAR+'edgar,E'+VAR+'eft,E'+VAR+'accmip,E'+VAR+'past,E'+VAR+'proj,E'+VAR+'stern,p_E'+VAR+'_bio')
for VAR in ['N2O']:
    exec('del E'+VAR+'epa,E'+VAR+'edgar,E'+VAR+'eft,E'+VAR+'ehyde,E'+VAR+'past,E'+VAR+'proj,E'+VAR+'davidson,p_E'+VAR+'_bio')

# ===========
# 1.11. PETERS
# ===========

# TODO


##################################################
#   2. LAND-USE CHANGE
##################################################

LUC = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_biome,nb_biome], dtype=dty) # {Mha/yr}
HARV = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_biome], dtype=dty) # {GtC/yr}
SHIFT = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_biome,nb_biome], dtype=dty) # {Mha/yr}

# =========
# 2.1. LUH1
# =========

# load land-use data from LUH1
# from [Hurtt et al., 2011] and updated
bio = ['des','for','shr','gra','cro','pas','urb']

# LUC
if (data_LULCC[:3] == 'LUH'):
    for b1 in range(len(bio)):
        for b2 in range(len(bio)):
            if os.path.isfile('data/LandUse_'+data_LULCC+'/#DATA.LandUse_'+data_LULCC+'_'+mod_LSNKcover+'.1501-2015_114reg1.LUC_'+bio[b1]+'2'+bio[b2]+'.csv'):
                TMP = np.array([line for line in csv.reader(open('data/LandUse_LUH1/#DATA.LandUse_'+data_LULCC+'_'+mod_LSNKcover+'.1501-2015_114reg1.LUC_'+bio[b1]+'2'+bio[b2]+'.csv','r'))], dtype=dty)
                for i in range(1,114+1):
                    LUC[1:min(ind_final,ind_cdiac)+1,regionJ_index[i],0,kLUC,regionI_index[i],biome_index[bio[b1]],biome_index[bio[b2]]] += TMP[200:200+min(ind_final,ind_cdiac)-1+1,i-1]

# HARV
if (data_LULCC[:3] == 'LUH'):
    for b in range(len(bio)):
        if os.path.isfile('data/LandUse_'+data_LULCC+'/#DATA.LandUse_'+data_LULCC+'_'+mod_LSNKcover+'.1501-2015_114reg1.HARV_'+bio[b]+'.csv'):
            TMP = np.array([line for line in csv.reader(open('data/LandUse_LUH1/#DATA.LandUse_'+data_LULCC+'_'+mod_LSNKcover+'.1501-2015_114reg1.HARV_'+bio[b]+'.csv','r'))], dtype=dty)
            for i in range(1,114+1):
                HARV[1:min(ind_final,ind_cdiac)+1,regionJ_index[i],0,kLUC,regionI_index[i],biome_index[bio[b]]] += TMP[200:200+min(ind_final,ind_cdiac)-1+1,i-1]

# SHIFT
if (data_LULCC[:3] == 'LUH'):
    for b1 in range(len(bio)):
        for b2 in range(b1,len(bio)):
            if os.path.isfile('data/LandUse_'+data_LULCC+'/#DATA.LandUse_'+data_LULCC+'_'+mod_LSNKcover+'.1501-2015_114reg1.SHIFT_'+bio[b1]+'2'+bio[b2]+'.csv'):
                TMP = np.array([line for line in csv.reader(open('data/LandUse_LUH1/#DATA.LandUse_'+data_LULCC+'_'+mod_LSNKcover+'.1501-2015_114reg1.SHIFT_'+bio[b1]+'2'+bio[b2]+'.csv','r'))], dtype=dty)
                for i in range(1,114+1):
                    SHIFT[1:min(ind_final,ind_cdiac)+1,regionJ_index[i],0,kLUC,regionI_index[i],biome_index[bio[b1]],biome_index[bio[b2]]] += TMP[200:200+min(ind_final,ind_cdiac)-1+1,i-1]

# ========
# 2.2. RCP
# ========

# initialization of projected drivers
LUCproj = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_biome,nb_biome], dtype=dty)
HARVproj = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_biome], dtype=dty)
SHIFTproj = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_biome,nb_biome], dtype=dty)

# projection of land-use data under RCP scenarios
# from [Hurtt et al., 2011] and [Meinshausen et al., 2011]

# LUC
if (scen_LULCC[:3] == 'RCP')&(ind_final > ind_cdiac):
    for b1 in range(len(bio)):
        for b2 in range(len(bio)):
            if os.path.isfile('data/LandUse_RCP/#DATA.LandUse_RCP_'+mod_LSNKcover+'.2006-2100_114reg1.rcp'+scen_LULCC[3]+scen_LULCC[5]+'_LUC_'+bio[b1]+'2'+bio[b2]+'.csv'):
                TMP = np.array([line for line in csv.reader(open('data/LandUse_RCP/#DATA.LandUse_RCP_'+mod_LSNKcover+'.2006-2100_114reg1.rcp'+scen_LULCC[3]+scen_LULCC[5]+'_LUC_'+bio[b1]+'2'+bio[b2]+'.csv','r'))], dtype=dty)
                for i in range(1,114+1):
                    LUCproj[306:min(ind_final,400)+1,regionJ_index[i],0,kLUC,regionI_index[i],biome_index[bio[b1]],biome_index[bio[b2]]] += TMP[:min(ind_final,400)-306+1,i-1]

# HARV
if (scen_LULCC[:3] == 'RCP')&(ind_final > ind_cdiac):
    for b in range(len(bio)):
        if os.path.isfile('data/LandUse_RCP/#DATA.LandUse_RCP_'+mod_LSNKcover+'.2006-2100_114reg1.rcp'+scen_LULCC[3]+scen_LULCC[5]+'_HARV_'+bio[b]+'.csv'):
            TMP = np.array([line for line in csv.reader(open('data/LandUse_RCP/#DATA.LandUse_RCP_'+mod_LSNKcover+'.2006-2100_114reg1.rcp'+scen_LULCC[3]+scen_LULCC[5]+'_HARV_'+bio[b]+'.csv','r'))], dtype=dty)
            for i in range(1,114+1):
                HARVproj[306:min(ind_final,400)+1,regionJ_index[i],0,kLUC,regionI_index[i],biome_index[bio[b]]] += TMP[:min(ind_final,400)-306+1,i-1]

# SHIFT
if (scen_LULCC[:3] == 'RCP')&(ind_final > ind_cdiac):
    for b1 in range(len(bio)):
        for b2 in range(b1,len(bio)):
            if os.path.isfile('data/LandUse_RCP/#DATA.LandUse_RCP_'+mod_LSNKcover+'.2006-2100_114reg1.rcp'+scen_LULCC[3]+scen_LULCC[5]+'_SHIFT_'+bio[b1]+'2'+bio[b2]+'.csv'):
                TMP = np.array([line for line in csv.reader(open('data/LandUse_RCP/#DATA.LandUse_RCP_'+mod_LSNKcover+'.2006-2100_114reg1.rcp'+scen_LULCC[3]+scen_LULCC[5]+'_SHIFT_'+bio[b1]+'2'+bio[b2]+'.csv','r'))], dtype=dty)
                for i in range(1,114+1):
                    SHIFTproj[306:min(ind_final,400)+1,regionJ_index[i],0,kLUC,regionI_index[i],biome_index[bio[b1]],biome_index[bio[b2]]] += TMP[:min(ind_final,400)-306+1,i-1]    

# ==================
# 2.A. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR in ['LUC','HARV','SHIFT']:

    # stop emissions
    if (scen_LULCC == 'stop')&(ind_final > ind_cdiac):
        exec(VAR+'[ind_cdiac+1:,...] = 0')

    # constant emissions
    elif (scen_LULCC == 'cst')&(ind_final > ind_cdiac):
        exec(VAR+'[ind_cdiac+1:,...] = '+VAR+'[ind_cdiac,...][np.newaxis,...]')        

    # RCP scenarios
    # always raw discontinuity
    elif (scen_LULCC[:3] == 'RCP')&(ind_final > ind_cdiac):
        exec(VAR+'[ind_cdiac+1:,...] = '+VAR+'proj[ind_cdiac+1:,...]')

# delete individual datasets
for VAR in ['LUC','HARV','SHIFT']:
    exec('del '+VAR+'proj')


##################################################
#   3. HALOGENATED COMPOUNDS
##################################################

EHFC = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_HFC], dtype=dty) # {kt/yr}
EPFC = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_PFC], dtype=dty) # {kt/yr}
EODS = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_ODS], dtype=dty) # {kt/yr}

# ==========
# 3.1. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]

# HFCs
EHFCedgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_HFC], dtype=dty)
for VAR in HFC:
    TMP = np.array([line for line in csv.reader(open('data/EHaloComp_EDGAR/#DATA.EHaloComp_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.E'+VAR+'.csv','r'))], dtype=dty)
    for i in range(114+1):
        EHFCedgar[270:ind_edgar+1,regionJ_index[i],0,kindGHG_index['HFC'],regionI_index[i],HFC.index(VAR)] += TMP[:ind_edgar-270+1,i]

# PFCs
EPFCedgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_PFC], dtype=dty)
for VAR in PFC:
    TMP = np.array([line for line in csv.reader(open('data/EHaloComp_EDGAR/#DATA.EHaloComp_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.E'+VAR+'.csv','r'))], dtype=dty)
    for i in range(114+1):
        EPFCedgar[270:ind_edgar+1,regionJ_index[i],0,kindGHG_index['PFC'],regionI_index[i],PFC.index(VAR)] += TMP[:ind_edgar-270+1,i]

# =============
# 3.2. EDGAR-FT
# =============

# load emissions from EDGAR-FT v4.2-FT2010
# see [JRC, 2013]

# HFCs
EHFCeft = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_HFC], dtype=dty)
for VAR in HFC:
    TMP = np.array([line for line in csv.reader(open('data/EHaloComp_EDGAR-FT/#DATA.EHaloComp_EDGAR-FT.2008-2010_114reg0.E'+VAR+'.csv','r'))], dtype=dty)
    for i in range(114+1):
        EHFCeft[308:310+1,regionJ_index[i],0,kindGHG_index['HFC'],regionI_index[i],HFC.index(VAR)] += TMP[:310-308+1,i]

# PFCs
EPFCeft = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_PFC], dtype=dty)
for VAR in PFC:
    TMP = np.array([line for line in csv.reader(open('data/EHaloComp_EDGAR-FT/#DATA.EHaloComp_EDGAR-FT.2008-2010_114reg0.E'+VAR+'.csv','r'))], dtype=dty)
    for i in range(114+1):
        EPFCeft[308:310+1,regionJ_index[i],0,kindGHG_index['PFC'],regionI_index[i],PFC.index(VAR)] += TMP[:310-308+1,i]

# ==========
# 3.3. CMIP5
# ==========

# ODSs
EODScmip5 = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_ODS], dtype=dty)
TMP = np.array([line for line in csv.reader(open('data/EHaloComp_CMIP5/#DATA.EHaloComp_CMIP5.1765-2005_(16ghg).EODS.csv','r'))][1:], dtype=dty)
lgd = [line for line in csv.reader(open('data/EHaloComp_CMIP5/#DATA.EHaloComp_CMIP5.1765-2005_(16ghg).EODS.csv','r'))][0]
for x in range(len(lgd)):
    EODScmip5[65:305+1,0,0,kindGHG_index['ODS'],0,ODS.index(lgd[x])] = TMP[:min(ind_final,305)-65+1,x]

# extend dataset following RCP unique projection
TMP = np.array([line for line in csv.reader(open('data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp85_EODS.csv','r'))][1:], dtype=dty)
lgd = [line for line in csv.reader(open('data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp85_EODS.csv','r'))][0]
for x in range(len(lgd)):
    EODScmip5[306:ind_cdiac+1,0,0,kindGHG_index['ODS'],0,ODS.index(lgd[x])] = TMP[6:min(ind_final,ind_cdiac)-300+1,x]

# ========
# 3.4. RCP
# ========

# initialization of projected drivers
for VAR in ['HFC','PFC','ODS']:
    exec('E'+VAR+'proj = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_'+VAR+'], dtype=dty)')

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# HFCs
if (scen_Ehalo[:3] == 'RCP')&(ind_final > ind_cdiac):
    for VAR in HFC:
        if (os.path.isfile('data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp'+scen_Ehalo[3]+scen_Ehalo[5]+'_E'+VAR+'.csv')):
            TMP = np.array([line for line in csv.reader(open('data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp'+scen_Ehalo[3]+scen_Ehalo[5]+'_E'+VAR+'.csv','r'))], dtype=dty)
            for i in range(4+1):
                if (mod_regionI == 'RCP5')&(mod_regionJ == 'RCP5'):
                    EHFCproj[300:min(ind_final,400)+1,i,0,kindGHG_index['HFC'],i,HFC.index(VAR)] += TMP[:min(ind_final,400)-300+1,i]
                elif (mod_regionI == 'RCP5')&(mod_regionJ != 'RCP5'):
                    EHFCproj[300:min(ind_final,400)+1,0,0,kindGHG_index['HFC'],i,HFC.index(VAR)] += TMP[:min(ind_final,400)-300+1,i]
                elif (mod_regionI != 'RCP5')&(mod_regionJ == 'RCP5'):
                    EHFCproj[300:min(ind_final,400)+1,i,0,kindGHG_index['HFC'],0,HFC.index(VAR)] += TMP[:min(ind_final,400)-300+1,i]
                else:
                    EHFCproj[300:min(ind_final,400)+1,0,0,kindGHG_index['HFC'],0,HFC.index(VAR)] += TMP[:min(ind_final,400)-300+1,i]

# PFCs
if (scen_Ehalo[:3] == 'RCP')&(ind_final > ind_cdiac):
    for VAR in PFC:
        if (os.path.isfile('data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp'+scen_Ehalo[3]+scen_Ehalo[5]+'_E'+VAR+'.csv')):
            TMP = np.array([line for line in csv.reader(open('data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp'+scen_Ehalo[3]+scen_Ehalo[5]+'_E'+VAR+'.csv','r'))], dtype=dty)
            for i in range(4+1):
                if (mod_regionI == 'RCP5')&(mod_regionJ == 'RCP5'):
                    EPFCproj[300:min(ind_final,400)+1,i,0,kindGHG_index['PFC'],i,PFC.index(VAR)] += TMP[:min(ind_final,400)-300+1,i]
                elif (mod_regionI == 'RCP5')&(mod_regionJ != 'RCP5'):
                    EPFCproj[300:min(ind_final,400)+1,0,0,kindGHG_index['PFC'],i,PFC.index(VAR)] += TMP[:min(ind_final,400)-300+1,i]
                elif (mod_regionI != 'RCP5')&(mod_regionJ == 'RCP5'):
                    EPFCproj[300:min(ind_final,400)+1,i,0,kindGHG_index['PFC'],0,PFC.index(VAR)] += TMP[:min(ind_final,400)-300+1,i]
                else:
                    EPFCproj[300:min(ind_final,400)+1,0,0,kindGHG_index['PFC'],0,PFC.index(VAR)] += TMP[:min(ind_final,400)-300+1,i]

# ODSs
if (scen_Ehalo[:3] == 'RCP')&(ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open('data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp'+scen_Ehalo[3]+scen_Ehalo[5]+'_EODS.csv','r'))][1:], dtype=dty)
    lgd = [line for line in csv.reader(open('data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp'+scen_Ehalo[3]+scen_Ehalo[5]+'_EODS.csv','r'))][0]
    for x in range(len(lgd)):
        EODSproj[300:min(ind_final,400)+1,0,0,kindGHG_index['ODS'],0,ODS.index(lgd[x])] = TMP[:min(ind_final,400)-300+1,x]

# =================
# 3.A. PAST DATASET
# =================

# datasets mixed following trends
for VAR in ['HFC','PFC','ODS']:

    if VAR in ['HFC']:
        
        # with EDGAR as reference
        if (data_Ehalo == 'EDGAR'):
            exec('E'+VAR+'past = E'+VAR+'edgar.copy()')
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar+1,ind_cdiac+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * E'+VAR+'eft[t,...]/E'+VAR+'eft[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'eft[t,...])/np.sum(E'+VAR+'eft[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')        
            # quadratic extrapolation before 1970
            # starting year of emission based on [Meinshausen et al., 2011]
            for x in range(nb_HFC):
                if (HFC[x] == 'HFC23'):
                    for t in range(230,270):
                        exec('E'+VAR+'past[t,...,x] = E'+VAR+'past[270,...,x] * ((t-230)/(270-230.))**2')

    elif VAR in ['PFC']:

        # with EDGAR as reference
        if (data_Ehalo == 'EDGAR'):
            exec('E'+VAR+'past = E'+VAR+'edgar.copy()')
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar+1,ind_cdiac+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * E'+VAR+'eft[t,...]/E'+VAR+'eft[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'eft[t,...])/np.sum(E'+VAR+'eft[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')        
            # quadratic extrapolation before 1970
            # starting year of emission based on [Meinshausen et al., 2011]
            for x in range(nb_PFC):
                if (PFC[x] == 'SF6'):
                    for t in range(250,270):
                        exec('E'+VAR+'past[t,...,x] = E'+VAR+'past[270,...,x] * ((t-250.)/(270-250.))**2')
                elif (PFC[x] == 'CF4'):
                    for t in range(222,270):
                        exec('E'+VAR+'past[t,...,x] = E'+VAR+'past[270,...,x] * ((t-222.)/(270-222.))**2')
                elif (PFC[x] == 'C2F6'):
                    for t in range(189,270):
                        exec('E'+VAR+'past[t,...,x] = E'+VAR+'past[270,...,x] * ((t-189)/(270-189.))**2')

    elif VAR in ['ODS']:

        # with CMIP5 as reference
        if (data_Ehalo == data_Ehalo):
            exec('E'+VAR+'past = E'+VAR+'cmip5.copy()')
            # linear extrapolation before 1765
            for t in range(50,65):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[65,...] * (t-50)/float(65-50)')

    # cut past dataset to right length
    exec('E'+VAR+'[:min(ind_cdiac,ind_final)+1,...] = E'+VAR+'past[:min(ind_cdiac,ind_final)+1,...]')

# ==================
# 3.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR in ['HFC','PFC','ODS']:

    # stop emissions
    if (scen_Ehalo == 'stop')&(ind_final > ind_cdiac):
        exec('E'+VAR+'[ind_cdiac+1:,...] = 0')

    # constant emissions
    elif (scen_Ehalo == 'cst')&(ind_final > ind_cdiac):
        exec('E'+VAR+'[ind_cdiac+1:,...] = E'+VAR+'[ind_cdiac,...][np.newaxis,...]')        

    # RCP scenarios
    elif (scen_Ehalo[:3] == 'RCP')&(ind_final > ind_cdiac):
        
        # raw discontinuity
        if (mod_DATAscen == 'raw'):
            exec('E'+VAR+'[ind_cdiac+1:,...] = E'+VAR+'proj[ind_cdiac+1:,...]')

        # offset at transition point
        elif (mod_DATAscen == 'offset'):
            exec('E'+VAR+'[ind_cdiac+1:,...] = E'+VAR+'proj[ind_cdiac+1:,...] - E'+VAR+'proj[ind_cdiac,...] + E'+VAR+'[ind_cdiac,...]')
            for t in range(ind_cdiac+1,ind_final+1):
                exec('def_regI = bool(np.sum(E'+VAR+'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E'+VAR+'proj[t,1:,...,:]))')
                if not def_regI:
                    exec('E'+VAR+'[t,:,...,0] += np.sum(E'+VAR+'[t,:,...,1:],-1)')
                    exec('E'+VAR+'[t,:,...,1:] = 0')
                if not def_regJ:
                    exec('E'+VAR+'[t,0,...,:] += np.sum(E'+VAR+'[t,1:,...,:],0)')
                    exec('E'+VAR+'[t,1:,...,:] = 0')                  

        # linear transition over N years
        elif (mod_DATAscen[:6] == 'smooth'):
            N = int(mod_DATAscen[6:])
            if (ind_final >= ind_cdiac+N):
                for t in range(ind_cdiac+1,ind_cdiac+N):
                    exec('E'+VAR+'[t,...] = (1-(t-ind_cdiac)/float(N)) * E'+VAR+'[ind_cdiac,...] + (t-ind_cdiac)/float(N) * E'+VAR+'proj[ind_cdiac+N,...]')
                    exec('def_regI = bool(np.sum(E'+VAR+'proj[t,:,...,1:]))')
                    exec('def_regJ = bool(np.sum(E'+VAR+'proj[t,1:,...,:]))')
                    if not def_regI:
                        exec('E'+VAR+'[t,:,...,0] += np.sum(E'+VAR+'[t,:,...,1:],-1)')
                        exec('E'+VAR+'[t,:,...,1:] = 0')
                    if not def_regJ:
                        exec('E'+VAR+'[t,0,...,:] += np.sum(E'+VAR+'[t,1:,...,:],0)')
                        exec('E'+VAR+'[t,1:,...,:] = 0')   
                exec('E'+VAR+'[ind_cdiac+N:,...] = E'+VAR+'proj[ind_cdiac+N:,...]')

        # follow trends of projection
        elif (mod_DATAscen == 'trends'):
            for t in range(ind_cdiac+1,ind_final+1):
                exec('def_regI = bool(np.sum(E'+VAR+'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E'+VAR+'proj[t,1:,...,:]))')
                if def_regI and def_regJ:
                    exec('E'+VAR+'[t,...] = E'+VAR+'[t-1,...] * E'+VAR+'proj[t,...]/E'+VAR+'proj[t-1,...]')
                    exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')
                    exec('E'+VAR+'[t,...] *= np.sum(E'+VAR+'[t-1,...])/np.sum(E'+VAR+'[t,...]) * np.sum(E'+VAR+'proj[t,...])/np.sum(E'+VAR+'proj[t-1,...])')
                elif not def_regI and def_regJ:
                    exec('E'+VAR+'[t,:,...,0] = np.sum(E'+VAR+'[t-1,:,...,:],0) * np.sum(E'+VAR+'proj[t,:,...,:],-1)/np.sum(E'+VAR+'proj[t-1,:,...,:],-1)')
                    exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')
                    exec('E'+VAR+'[t,...] *= np.sum(E'+VAR+'[t-1,...])/np.sum(E'+VAR+'[t,...]) * np.sum(E'+VAR+'proj[t,...])/np.sum(E'+VAR+'proj[t-1,...])')
                elif def_regI and not def_regJ:
                    exec('E'+VAR+'[t,0,...,:] = np.sum(E'+VAR+'[t-1,:,...,:],0) * np.sum(E'+VAR+'proj[t,:,...,:],0)/np.sum(E'+VAR+'proj[t-1,:,...,:],0)')
                    exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')
                    exec('E'+VAR+'[t,...] *= np.sum(E'+VAR+'[t-1,...])/np.sum(E'+VAR+'[t,...]) * np.sum(E'+VAR+'proj[t,...])/np.sum(E'+VAR+'proj[t-1,...])')
                elif not def_regI and not def_regJ:
                    exec('E'+VAR+'[t,0,...,0] = np.sum(np.sum(E'+VAR+'[t-1,:,...,:],-1),0) * np.sum(np.sum(E'+VAR+'proj[t,:,...,:],-1),0)/np.sum(np.sum(E'+VAR+'proj[t-1,:,...,:],-1),0)')
            exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')

# delete individual datasets
for VAR in ['HFC','PFC']:
    exec('del E'+VAR+'edgar,E'+VAR+'eft,E'+VAR+'past,E'+VAR+'proj')
for VAR in ['ODS']:
    exec('del E'+VAR+'cmip5,E'+VAR+'past,E'+VAR+'proj')


##################################################
#   4. SHORT-LIVED SPECIES
##################################################

ENOX = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty) # {TgN/yr}
ECO = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty) # {TgC/yr}
EVOC = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty) # {Tg/yr}
ESO2 = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty) # {TgS/yr}
ENH3 = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty) # {TgN/yr}
EOC = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty) # {Tg/yr}
EBC = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty) # {Tg/yr}

# ==========
# 4.1. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]

# OzoPrec
for VAR in ['NOX','CO','VOC']:
    exec('E'+VAR+'edgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1,len(sec_accmip)-2):
        if os.path.isfile('data/EOzoPrec_EDGAR/#DATA.EOzoPrec_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv'):
            TMP = np.array([line for line in csv.reader(open('data/EOzoPrec_EDGAR/#DATA.EOzoPrec_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
            for i in range(114+1):
                exec('E'+VAR+'edgar[270:ind_edgar+1,regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]')

# AeroPrec
for VAR in ['SO2','NH3']:
    exec('E'+VAR+'edgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1,len(sec_accmip)-2):
        if os.path.isfile('data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv'):
            TMP = np.array([line for line in csv.reader(open('data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
            for i in range(114+1):
                exec('E'+VAR+'edgar[270:ind_edgar+1,regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]')

# PM10 (proxy for OC/BC)
EPM10edgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)
for s in range(1,len(sec_accmip)-2):
    if os.path.isfile('data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.EPM10_'+sec_accmip[s]+'.csv'):
        TMP = np.array([line for line in csv.reader(open('data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-'+str(1700+ind_edgar)+'_114reg0.EPM10_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
        for i in range(114+1):
            EPM10edgar[270:ind_edgar+1,regionJ_index[i],0,kindAER_index['OC'],regionI_index[i]] += TMP[:ind_edgar-270+1,i]/2
            EPM10edgar[270:ind_edgar+1,regionJ_index[i],0,kindAER_index['BC'],regionI_index[i]] += TMP[:ind_edgar-270+1,i]/2

# ===============
# 4.2. EDGAR-HTAP
# ===============

# load emissions from EDGAR-HTAP v2
# see [JRC, 2013]

# OzoPrec
for VAR in ['NOX','CO','VOC']:
    exec('E'+VAR+'ehtap = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1,len(sec_accmip)-2):
        if os.path.isfile('data/EOzoPrec_EDGAR-HTAP/#DATA.EOzoPrec_EDGAR-HTAP.2008-2010_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv'):
            TMP = np.array([line for line in csv.reader(open('data/EOzoPrec_EDGAR-HTAP/#DATA.EOzoPrec_EDGAR-HTAP.2008-2010_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
            for i in range(114+1):
                exec('E'+VAR+'ehtap[308:310+1,regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]')

# AeroPrec
for VAR in ['SO2','NH3','OC','BC']:
    exec('E'+VAR+'ehtap = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1,len(sec_accmip)-2):
        if os.path.isfile('data/EAeroPrec_EDGAR-HTAP/#DATA.EAeroPrec_EDGAR-HTAP.2008-2010_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv'):
            TMP = np.array([line for line in csv.reader(open('data/EAeroPrec_EDGAR-HTAP/#DATA.EAeroPrec_EDGAR-HTAP.2008-2010_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
            for i in range(114+1):
                exec('E'+VAR+'ehtap[308:310+1,regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]')

# ===========
# 4.3. ACCMIP
# ===========

# load emissions from ACCMIP
# from [Lamarque et al., 2010]

# OzoPrec
for VAR in ['NOX','CO','VOC']:
    exec('E'+VAR+'accmip = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    exec('p_E'+VAR+'_bio = np.zeros([nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1,len(sec_accmip)-2):
        if os.path.isfile('data/EOzoPrec_ACCMIP/#DATA.EOzoPrec_ACCMIP.1850-2000_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv'):
            TMP = np.array([line for line in csv.reader(open('data/EOzoPrec_ACCMIP/#DATA.EOzoPrec_ACCMIP.1850-2000_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
            for i in range(114+1):
                exec('E'+VAR+'accmip[150:300+1,regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[:300-150+1,i]')
                if sec_accmip[s] in ['agr','awb','wst']:
                    exec('p_E'+VAR+'_bio[regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[0,i]')
    exec('p_E'+VAR+'_bio /= E'+VAR+'accmip[150]')
    exec('p_E'+VAR+'_bio[np.isnan(p_E'+VAR+'_bio)|np.isinf(p_E'+VAR+'_bio)] = 0')


# AeroPrec
for VAR in ['SO2','NH3','OC','BC']:
    exec('E'+VAR+'accmip = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    exec('p_E'+VAR+'_bio = np.zeros([nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1,len(sec_accmip)-2):
        if os.path.isfile('data/EAeroPrec_ACCMIP/#DATA.EAeroPrec_ACCMIP.1850-2000_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv'):
            TMP = np.array([line for line in csv.reader(open('data/EAeroPrec_ACCMIP/#DATA.EAeroPrec_ACCMIP.1850-2000_114reg0.E'+VAR+'_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
            for i in range(114+1):
                exec('E'+VAR+'accmip[150:300+1,regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[:300-150+1,i]')
                if sec_accmip[s] in ['agr','awb','wst']:
                    exec('p_E'+VAR+'_bio[regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[0,i]')
    exec('p_E'+VAR+'_bio /= E'+VAR+'accmip[150]')
    exec('p_E'+VAR+'_bio[np.isnan(p_E'+VAR+'_bio)|np.isinf(p_E'+VAR+'_bio)] = 0')


# =========
# 4.4. SRES
# =========

# initialization of projected drivers
for VAR in ['NOX','CO','VOC']+['SO2','NH3','OC','BC']:
    exec('E'+VAR+'proj = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')

# projection of emissions under SRES scenarios
# from [IPCC, 2000]

# OzoPrec
for VAR in ['NOX','CO','VOC']:
    exec('scen = scen_E'+VAR)
    if (scen[:4] == 'SRES')&(ind_final > ind_cdiac):
        TMP = np.array([line for line in csv.reader(open('data/EOzoPrec_SRES/#DATA.EOzoPrec_SRES.2000-2100_4reg0.'+scen[5:]+'_E'+VAR+'.csv','r'))], dtype=dty)
        for i in range(4+1):
            if (mod_regionI == 'SRES4')&(mod_regionJ == 'SRES4'):
                exec('E'+VAR+'proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
            elif (mod_regionI == 'SRES4')&(mod_regionJ != 'SRES4'):
                exec('E'+VAR+'proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
            elif (mod_regionI != 'SRES4')&(mod_regionJ == 'SRES4'):
                exec('E'+VAR+'proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')
            elif(mod_regionI != 'SRES4')&(mod_regionJ != 'SRES4'):
                exec('E'+VAR+'proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')

# AeroPrec (SO2 only)
if (scen_ESO2[:4] == 'SRES')&(ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open('data/EAeroPrec_SRES/#DATA.EAeroPrec_SRES.2000-2100_4reg0.'+scen_ESO2[5:]+'_ESO2.csv','r'))], dtype=dty)
    for i in range(4+1):
        if (mod_regionI == 'SRES4')&(mod_regionJ == 'SRES4'):
            ESO2proj[300:min(ind_final,400)+1,i,0,kindAER_index['SO2'],i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI == 'SRES4')&(mod_regionJ != 'SRES4'):
            ESO2proj[300:min(ind_final,400)+1,0,0,kindAER_index['SO2'],i] += TMP[:min(ind_final,400)-300+1,i]
        elif (mod_regionI != 'SRES4')&(mod_regionJ == 'SRES4'):
            ESO2proj[300:min(ind_final,400)+1,i,0,kindAER_index['SO2'],0] += TMP[:min(ind_final,400)-300+1,i]
        elif(mod_regionI != 'SRES4')&(mod_regionJ != 'SRES4'):
            ESO2proj[300:min(ind_final,400)+1,0,0,kindAER_index['SO2'],0] += TMP[:min(ind_final,400)-300+1,i]

# ========
# 4.5. RCP
# ========

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# OzoPrec
for VAR in ['NOX','CO','VOC']:
    exec('scen = scen_E'+VAR)
    if (scen[:3] == 'RCP')&(ind_final > ind_cdiac):
        for s in range(1,len(sec_accmip)-2):
            if os.path.isfile('data/EOzoPrec_RCP/#DATA.EOzoPrec_RCP.2000-2100_5reg0.rcp'+scen[3]+scen[5]+'_E'+VAR+'_'+sec_accmip[s]+'.csv'):
                TMP = np.array([line for line in csv.reader(open('data/EOzoPrec_RCP/#DATA.EOzoPrec_RCP.2000-2100_5reg0.rcp'+scen[3]+scen[5]+'_E'+VAR+'_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
                for i in range(5+1):
                    if (mod_regionI == 'RCP5')&(mod_regionJ == 'RCP5'):
                        exec('E'+VAR+'proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI == 'RCP5')&(mod_regionJ != 'RCP5'):
                        exec('E'+VAR+'proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI != 'RCP5')&(mod_regionJ == 'RCP5'):
                        exec('E'+VAR+'proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI != 'RCP5')&(mod_regionJ != 'RCP5'):
                        exec('E'+VAR+'proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')

# AeroPrec
for VAR in ['SO2','NH3','OC','BC']:
    exec('scen = scen_E'+VAR)
    if (scen[:3] == 'RCP')&(ind_final > ind_cdiac):
        for s in range(1,len(sec_accmip)-2):
            if os.path.isfile('data/EAeroPrec_RCP/#DATA.EAeroPrec_RCP.2000-2100_5reg0.rcp'+scen[3]+scen[5]+'_E'+VAR+'_'+sec_accmip[s]+'.csv'):
                TMP = np.array([line for line in csv.reader(open('data/EAeroPrec_RCP/#DATA.EAeroPrec_RCP.2000-2100_5reg0.rcp'+scen[3]+scen[5]+'_E'+VAR+'_'+sec_accmip[s]+'.csv','r'))], dtype=dty)
                for i in range(5+1):
                    if (mod_regionI == 'RCP5')&(mod_regionJ == 'RCP5'):
                        exec('E'+VAR+'proj[300:min(ind_final,400)+1,i,0,kindAER_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI == 'RCP5')&(mod_regionJ != 'RCP5'):
                        exec('E'+VAR+'proj[300:min(ind_final,400)+1,0,0,kindAER_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI != 'RCP5')&(mod_regionJ == 'RCP5'):
                        exec('E'+VAR+'proj[300:min(ind_final,400)+1,i,0,kindAER_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI != 'RCP5')&(mod_regionJ != 'RCP5'):
                        exec('E'+VAR+'proj[300:min(ind_final,400)+1,0,0,kindAER_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')

# =================
# 4.A. PAST DATASET
# =================

# datasets mixed following trends
for VAR in ['NOX','CO','VOC']+['SO2','NH3','OC','BC']:

    if VAR in ['NOX','CO','VOC']+['SO2','NH3']:
        exec('data = data_E'+VAR)
        
        # with EDGAR as reference
        if (data == 'EDGAR'):
            exec('E'+VAR+'past = E'+VAR+'edgar.copy()')
            # follow EDGAR-HTAP variations after 2008
            for t in range(ind_edgar+1,ind_cdiac+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * E'+VAR+'ehtap[t,...]/E'+VAR+'ehtap[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'ehtap[t,...])/np.sum(E'+VAR+'ehtap[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')        
            # follow ACCMIP variations before 1970
            for t in range(150,270)[::-1]:
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t+1,...] * E'+VAR+'accmip[t,...]/E'+VAR+'accmip[t+1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t+1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'accmip[t,...])/np.sum(E'+VAR+'accmip[t+1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')
            # linear extrapolation before 1850
            for t in range(50,150):
                exec('E'+VAR+'past[t,...] = (1-p_E'+VAR+'_bio) * E'+VAR+'past[150,...] * (t-50)/float(150-50)')
            for t in range(0,150):
                exec('E'+VAR+'past[t,...] += p_E'+VAR+'_bio * E'+VAR+'past[150,...] * (t--200)/float(150--200)')            
            exec('E'+VAR+'_0 = E'+VAR+'past[50,...]')

        # with ACCMIP as reference
        elif (data == 'ACCMIP'):
            exec('E'+VAR+'past = E'+VAR+'accmip.copy()')
            # follow EDGAR variations after 2000
            for t in range(300+1,ind_edgar+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * E'+VAR+'edgar[t,...]/E'+VAR+'edgar[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'edgar[t,...])/np.sum(E'+VAR+'edgar[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')        
            # follow EDGAR-HTAP variations after 2008
            for t in range(ind_edgar+1,ind_cdiac+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * E'+VAR+'ehtap[t,...]/E'+VAR+'ehtap[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'ehtap[t,...])/np.sum(E'+VAR+'ehtap[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')        
            # linear extrapolation before 1850
            for t in range(50,150):
                exec('E'+VAR+'past[t,...] = (1-p_E'+VAR+'_bio) * E'+VAR+'past[150,...] * (t-50)/float(150-50)')
            for t in range(0,150):
                exec('E'+VAR+'past[t,...] += p_E'+VAR+'_bio * E'+VAR+'past[150,...] * (t--200)/float(150--200)')            
            exec('E'+VAR+'_0 = E'+VAR+'past[50,...]')

    elif VAR in ['OC','BC']:
        exec('data = data_E'+VAR)

        # with ACCMIP as reference
        if (data == 'ACCMIP'):
            exec('E'+VAR+'past = E'+VAR+'accmip.copy()')
            # follow EDGAR (PM10) variations after 2000
            for t in range(300+1,ind_edgar+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * EPM10edgar[t,...]/EPM10edgar[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(EPM10edgar[t,...])/np.sum(EPM10edgar[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')    
            # follow EDGAR-HTAP variations after 2008
            for t in range(ind_edgar+1,ind_cdiac+1):
                exec('E'+VAR+'past[t,...] = E'+VAR+'past[t-1,...] * E'+VAR+'ehtap[t,...]/E'+VAR+'ehtap[t-1,...]')
                exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0') 
                exec('E'+VAR+'past[t,...] *= np.sum(E'+VAR+'past[t-1,...])/np.sum(E'+VAR+'past[t,...]) * np.sum(E'+VAR+'ehtap[t,...])/np.sum(E'+VAR+'ehtap[t-1,...])')
            exec('E'+VAR+'past[np.isnan(E'+VAR+'past)|np.isinf(E'+VAR+'past)] = 0')        
            # linear extrapolation before 1850
            for t in range(50,150):
                exec('E'+VAR+'past[t,...] = (1-p_E'+VAR+'_bio) * E'+VAR+'past[150,...] * (t-50)/float(150-50)')
            for t in range(0,150):
                exec('E'+VAR+'past[t,...] += p_E'+VAR+'_bio * E'+VAR+'past[150,...] * (t--200)/float(150--200)')            
            exec('E'+VAR+'_0 = E'+VAR+'past[50,...]')

    # cut past dataset to right length
    exec('E'+VAR+'[:min(ind_cdiac,ind_final)+1,...] = E'+VAR+'past[:min(ind_cdiac,ind_final)+1,...]')

# ==================
# 4.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR in ['NOX','CO','VOC']+['SO2','NH3','OC','BC']:
    exec('scen = scen_E'+VAR)

    # stop emissions
    if (scen == 'stop')&(ind_final > ind_cdiac):
        exec('E'+VAR+'[ind_cdiac+1:,...] = E'+VAR+'_0[np.newaxis,...]')

    # constant emissions
    elif (scen == 'cst')&(ind_final > ind_cdiac):
        exec('E'+VAR+'[ind_cdiac+1:,...] = E'+VAR+'[ind_cdiac,...][np.newaxis,...]')        

    # RCP or SRES scenarios
    elif ((scen[:4] == 'SRES')|(scen[:3] == 'RCP'))&(ind_final > ind_cdiac):
        
        # raw discontinuity
        if (mod_DATAscen == 'raw'):
            exec('E'+VAR+'[ind_cdiac+1:,...] = E'+VAR+'proj[ind_cdiac+1:,...]')

        # offset at transition point
        elif (mod_DATAscen == 'offset'):
            exec('E'+VAR+'[ind_cdiac+1:,...] = E'+VAR+'proj[ind_cdiac+1:,...] - E'+VAR+'proj[ind_cdiac,...] + E'+VAR+'[ind_cdiac,...]')
            for t in range(ind_cdiac+1,ind_final+1):
                exec('def_regI = bool(np.sum(E'+VAR+'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E'+VAR+'proj[t,1:,...,:]))')
                if not def_regI:
                    exec('E'+VAR+'[t,:,...,0] += np.sum(E'+VAR+'[t,:,...,1:],-1)')
                    exec('E'+VAR+'[t,:,...,1:] = 0')
                if not def_regJ:
                    exec('E'+VAR+'[t,0,...,:] += np.sum(E'+VAR+'[t,1:,...,:],0)')
                    exec('E'+VAR+'[t,1:,...,:] = 0')                  

        # linear transition over N years
        elif (mod_DATAscen[:6] == 'smooth'):
            N = int(mod_DATAscen[6:])
            if (ind_final >= ind_cdiac+N):
                for t in range(ind_cdiac+1,ind_cdiac+N):
                    exec('E'+VAR+'[t,...] = (1-(t-ind_cdiac)/float(N)) * E'+VAR+'[ind_cdiac,...] + (t-ind_cdiac)/float(N) * E'+VAR+'proj[ind_cdiac+N,...]')
                    exec('def_regI = bool(np.sum(E'+VAR+'proj[t,:,...,1:]))')
                    exec('def_regJ = bool(np.sum(E'+VAR+'proj[t,1:,...,:]))')
                    if not def_regI:
                        exec('E'+VAR+'[t,:,...,0] += np.sum(E'+VAR+'[t,:,...,1:],-1)')
                        exec('E'+VAR+'[t,:,...,1:] = 0')
                    if not def_regJ:
                        exec('E'+VAR+'[t,0,...,:] += np.sum(E'+VAR+'[t,1:,...,:],0)')
                        exec('E'+VAR+'[t,1:,...,:] = 0')   
                exec('E'+VAR+'[ind_cdiac+N:,...] = E'+VAR+'proj[ind_cdiac+N:,...]')

        # follow trends of projection
        elif (mod_DATAscen == 'trends'):
            for t in range(ind_cdiac+1,ind_final+1):
                exec('def_regI = bool(np.sum(E'+VAR+'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E'+VAR+'proj[t,1:,...,:]))')
                if def_regI and def_regJ:
                    exec('E'+VAR+'[t,...] = E'+VAR+'[t-1,...] * E'+VAR+'proj[t,...]/E'+VAR+'proj[t-1,...]')
                    exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')
                    exec('E'+VAR+'[t,...] *= np.sum(E'+VAR+'[t-1,...])/np.sum(E'+VAR+'[t,...]) * np.sum(E'+VAR+'proj[t,...])/np.sum(E'+VAR+'proj[t-1,...])')
                elif not def_regI and def_regJ:
                    exec('E'+VAR+'[t,:,...,0] = np.sum(E'+VAR+'[t-1,:,...,:],-1) * np.sum(E'+VAR+'proj[t,:,...,:],-1)/np.sum(E'+VAR+'proj[t-1,:,...,:],-1)')
                    exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')
                    exec('E'+VAR+'[t,...] *= np.sum(E'+VAR+'[t-1,...])/np.sum(E'+VAR+'[t,...]) * np.sum(E'+VAR+'proj[t,...])/np.sum(E'+VAR+'proj[t-1,...])')
                elif def_regI and not def_regJ:
                    exec('E'+VAR+'[t,0,...,:] = np.sum(E'+VAR+'[t-1,:,...,:],0) * np.sum(E'+VAR+'proj[t,:,...,:],0)/np.sum(E'+VAR+'proj[t-1,:,...,:],0)')
                    exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')
                    exec('E'+VAR+'[t,...] *= np.sum(E'+VAR+'[t-1,...])/np.sum(E'+VAR+'[t,...]) * np.sum(E'+VAR+'proj[t,...])/np.sum(E'+VAR+'proj[t-1,...])')
                elif not def_regI and not def_regJ:
                    exec('E'+VAR+'[t,0,...,0] = np.sum(np.sum(E'+VAR+'[t-1,:,...,:],-1),0) * np.sum(np.sum(E'+VAR+'proj[t,:,...,:],-1),0)/np.sum(np.sum(E'+VAR+'proj[t-1,:,...,:],-1),0)')
            exec('E'+VAR+'[np.isnan(E'+VAR+')|np.isinf(E'+VAR+')] = 0')

# delete individual datasets
for VAR in ['NOX','CO','VOC']+['SO2','NH3']:
    exec('del E'+VAR+'edgar,E'+VAR+'ehtap,E'+VAR+'accmip,E'+VAR+'past,E'+VAR+'proj,p_E'+VAR+'_bio')
for VAR in ['OC','BC']:
    exec('del E'+VAR+'ehtap,E'+VAR+'accmip,E'+VAR+'past,E'+VAR+'proj,p_E'+VAR+'_bio')
for VAR in ['PM10']:
    exec('del E'+VAR+'edgar')


##################################################
#   5. RADIATIVE FORCING
##################################################

RFcon = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind], dtype=dty) # {W/m2}
RFvolc = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind], dtype=dty) # {W/m2}
RFsolar = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind], dtype=dty) # {W/m2}

# ==================
# 5.1. Anthropogenic
# ==================

# load RF estimates from IPCC AR5
# from [IPCC, 2013] (annexe 2)
if (data_RFant == 'IPCC-AR5'):
    TMP = np.array([line for line in csv.reader(open('data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv','r'))][1:], dtype=dty)
    lgd = [line for line in csv.reader(open('data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv','r'))][0]
    # past estimates
    for x in range(len(lgd)):
        if (lgd[x] == 'Contrails'):
            RFcon[51:min(ind_final,ind_cdiac)+1,0,0,kindRF_index['RFcon']] += TMP[1:1+min(ind_final,ind_cdiac)-51+1,x]
    # arbitrary extension
    if (scen_RFant == 'cst')&(ind_final > ind_cdiac):
        for x in range(len(lgd)):
            if (lgd[x] == 'Contrails'):
                RFcon[ind_cdiac+1:,0,0,kindRF_index['RFcon']] = TMP[-1,x]

# ============
# 5.2. Natural
# ============

# load RF estimates from IPCC AR5
# from [IPCC, 2013] (annexe 2)
if (data_RFnat == 'IPCC-AR5'):
    TMP = np.array([line for line in csv.reader(open('data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv','r'))][1:], dtype=dty)
    lgd = [line for line in csv.reader(open('data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv','r'))][0]
    # past estimates
    for x in range(len(lgd)):
        if (lgd[x] == 'Volcano'):
            RFvolc[51:min(ind_final,ind_cdiac)+1,0,0,kindRF_index['RFvol']] += TMP[1:1+min(ind_final,ind_cdiac)-51+1,x] - np.mean(TMP[1:,x])          
        elif (lgd[x] == 'Solar'):
            RFsolar[51:min(ind_final,ind_cdiac)+1,0,0,kindRF_index['RFsol']] += TMP[1:1+min(ind_final,ind_cdiac)-51+1,x]
    # arbitrary extension
    if (scen_RFnat == 'cst')&(ind_final > ind_cdiac):
        for x in range(len(lgd)):
            if (lgd[x] == 'Volcano'):
                RFvolc[ind_cdiac+1:,0,0,kindRF_index['RFvol']] = 0.
            elif (lgd[x] == 'Solar'):
                RFsolar[ind_cdiac+1:,0,0,kindRF_index['RFsol']] = np.mean(TMP[-11:,x])


##################################################
#   B. FINAL DRIVERS
##################################################

# ==================
# B.1. PREINDUSTRIAL
# ==================

# reference emissions for preindustrial
for VAR in ['ECH4','EN2O']+['ENOX','ECO','EVOC','ESO2','ENH3','EOC','EBC']:
    exec(VAR+'[:,...] -= '+VAR+'_0[np.newaxis,...]')

# force true 1750 preindustrial (drivers)
if PI_1750:
    for VAR in ['EFF','ECH4','EN2O']+['LUC','HARV','SHIFT']+['EHFC','EPFC','EODS']+['ENOX','ECO','EVOC','ESO2','ENH3','EOC','EBC']+['RFcon','RFvolc','RFsolar']:
        exec(VAR+'[:50+1] *= 0')

# ================
# B.2. ATTRIBUTION
# ================

# starting year of regional attribution
for VAR in ['EFF','ECH4','EN2O']+['LUC','HARV','SHIFT']+['EHFC','EPFC','EODS']+['ENOX','ECO','EVOC','ESO2','ENH3','EOC','EBC']+['RFcon','RFvolc','RFsolar']:
    exec(VAR+'[:ind_attrib,0,:,:,...] = np.sum('+VAR+'[:ind_attrib,:,:,:,...],1)')
    exec(VAR+'[:ind_attrib,1:,:,:,...] = 0')

# timeframed sectoral attribution
if (mod_sector == 'Time')&(300 < ind_final <= 310):
    for VAR in ['EFF','ECH4','EN2O']+['LUC','HARV','SHIFT']+['EHFC','EPFC','EODS']+['ENOX','ECO','EVOC','ESO2','ENH3','EOC','EBC']+['RFcon','RFvolc','RFsolar']:
        exec(VAR+'[150:200,:,:,1,:,...] = '+VAR+'[150:200,:,:,0,:,...]')
        exec(VAR+'[200:210,:,:,2,:,...] = '+VAR+'[200:210,:,:,0,:,...]')
        exec(VAR+'[210:220,:,:,3,:,...] = '+VAR+'[210:220,:,:,0,:,...]')
        exec(VAR+'[220:230,:,:,4,:,...] = '+VAR+'[220:230,:,:,0,:,...]')
        exec(VAR+'[230:240,:,:,5,:,...] = '+VAR+'[230:240,:,:,0,:,...]')
        exec(VAR+'[240:250,:,:,6,:,...] = '+VAR+'[240:250,:,:,0,:,...]')
        exec(VAR+'[250:260,:,:,7,:,...] = '+VAR+'[250:260,:,:,0,:,...]')
        exec(VAR+'[260:265,:,:,8,:,...] = '+VAR+'[260:265,:,:,0,:,...]')
        exec(VAR+'[265:270,:,:,9,:,...] = '+VAR+'[265:270,:,:,0,:,...]')
        exec(VAR+'[270:275,:,:,10,:,...] = '+VAR+'[270:275,:,:,0,:,...]')
        exec(VAR+'[275:280,:,:,11,:,...] = '+VAR+'[275:280,:,:,0,:,...]')
        exec(VAR+'[280:285,:,:,12,:,...] = '+VAR+'[280:285,:,:,0,:,...]')
        exec(VAR+'[285:290,:,:,13,:,...] = '+VAR+'[285:290,:,:,0,:,...]')
        exec(VAR+'[290:292,:,:,14,:,...] = '+VAR+'[290:292,:,:,0,:,...]')
        exec(VAR+'[292:294,:,:,15,:,...] = '+VAR+'[292:294,:,:,0,:,...]')
        exec(VAR+'[294:296,:,:,16,:,...] = '+VAR+'[294:296,:,:,0,:,...]')
        exec(VAR+'[296:298,:,:,17,:,...] = '+VAR+'[296:298,:,:,0,:,...]')
        exec(VAR+'[298:300,:,:,18,:,...] = '+VAR+'[298:300,:,:,0,:,...]')
        for t in range(300,ind_final+1):
            exec(VAR+'[t,:,:,19+t-300,:,...] = '+VAR+'[t,:,:,0,:,...]')
        exec(VAR+'[150:,:,:,0,:,...] = 0')
elif (mod_sector == 'TimeRCP')&(ind_final == 400):
    for VAR in ['EFF','ECH4','EN2O']+['LUC','HARV','SHIFT']+['EHFC','EPFC','EODS']+['ENOX','ECO','EVOC','ESO2','ENH3','EOC','EBC']+['RFcon','RFvolc','RFsolar']:
        exec(VAR+'[311:316,:,1,:,...] = '+VAR+'[311:316,:,0,:,...]')
        exec(VAR+'[316:321,:,2,:,...] = '+VAR+'[316:321,:,0,:,...]')
        exec(VAR+'[321:326,:,3,:,...] = '+VAR+'[321:326,:,0,:,...]')
        exec(VAR+'[326:331,:,4,:,...] = '+VAR+'[326:331,:,0,:,...]')
        exec(VAR+'[331:336,:,5,:,...] = '+VAR+'[331:336,:,0,:,...]')
        exec(VAR+'[336:341,:,6,:,...] = '+VAR+'[336:341,:,0,:,...]')
        exec(VAR+'[341:346,:,7,:,...] = '+VAR+'[341:346,:,0,:,...]')
        exec(VAR+'[346:351,:,8,:,...] = '+VAR+'[346:351,:,0,:,...]')
        exec(VAR+'[351:356,:,9,:,...] = '+VAR+'[351:356,:,0,:,...]')
        exec(VAR+'[356:361,:,10,:,...] = '+VAR+'[356:361,:,0,:,...]')
        exec(VAR+'[361:366,:,11,:,...] = '+VAR+'[361:366,:,0,:,...]')
        exec(VAR+'[366:371,:,12,:,...] = '+VAR+'[366:371,:,0,:,...]')
        exec(VAR+'[371:376,:,13,:,...] = '+VAR+'[371:376,:,0,:,...]')
        exec(VAR+'[376:381,:,14,:,...] = '+VAR+'[376:381,:,0,:,...]')
        exec(VAR+'[381:386,:,15,:,...] = '+VAR+'[381:386,:,0,:,...]')
        exec(VAR+'[386:391,:,16,:,...] = '+VAR+'[386:391,:,0,:,...]')
        exec(VAR+'[391:396,:,17,:,...] = '+VAR+'[391:396,:,0,:,...]')
        exec(VAR+'[396:401,:,18,:,...] = '+VAR+'[396:401,:,0,:,...]')
        exec(VAR+'[311:,:,0,:,...] = 0')

