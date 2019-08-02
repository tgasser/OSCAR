import csv

import numpy as np


##################################################################
##################################################################
##################################################################

execfile('OSCAR.py')

nb_run = 500
break_if_error = False

##################################################################

tag_break = False
drive_list = []
param_list = []
# runs
for nrun in range(nb_run):
    print str(nrun+1)+'/'+str(nb_run)

    # drivers
    drive = ['CDIAC','EDGAR']   
    data_EFF = drive[np.random.random_integers(len(drive))-1]
    drive = ['LUH1']   
    data_LULCC = drive[np.random.random_integers(len(drive))-1]
    drive = ['EDGAR','ACCMIP','EPA']     
    data_ECH4 = drive[np.random.random_integers(len(drive))-1]
    drive = ['EDGAR','EPA']   
    data_EN2O = drive[np.random.random_integers(len(drive))-1]
    drive = ['EDGAR']   
    data_Ehalo = drive[np.random.random_integers(len(drive))-1]
    drive = ['EDGAR','ACCMIP']    
    data_ENOX = drive[np.random.random_integers(len(drive))-1]
    drive = ['EDGAR','ACCMIP']    
    data_ECO = drive[np.random.random_integers(len(drive))-1]    
    drive = ['EDGAR','ACCMIP']   
    data_EVOC = drive[np.random.random_integers(len(drive))-1]
    drive = ['EDGAR','ACCMIP']    
    data_ESO2 = drive[np.random.random_integers(len(drive))-1]
    drive = ['EDGAR','ACCMIP']   
    data_ENH3 = drive[np.random.random_integers(len(drive))-1]
    drive = ['ACCMIP']   
    data_EOC = drive[np.random.random_integers(len(drive))-1]   
    drive = ['ACCMIP']   
    data_EBC = drive[np.random.random_integers(len(drive))-1]       
    drive = ['IPCC-AR5']   
    data_RFant = drive[np.random.random_integers(len(drive))-1]   
    drive = ['IPCC-AR5']   
    data_RFnat = drive[np.random.random_integers(len(drive))-1]    
    
    # ocean
    param = ['HILDA','BD-model','2D-model','3D-model']
    mod_OSNKstruct = param[np.random.random_integers(len(param))-1]
    param = ['CO2SysPade','CO2SysPower']
    mod_OSNKchem = param[np.random.random_integers(len(param))-1]
    param = ['CESM1-BGC','IPSL-CM5A-LR','MPI-ESM-LR']
    mod_OSNKtrans = param[np.random.random_integers(len(param))-1]
    
    # biosphere
    param = ['log','hyp']
    mod_LSNKnpp = param[np.random.random_integers(len(param))-1]
    param = ['exp','gauss']    
    mod_LSNKrho = param[np.random.random_integers(len(param))-1]
    param = ['CLM-45','JSBACH','JULES','LPJ','LPJ-GUESS','LPX-Bern','OCN','ORCHIDEE','VISIT']
    mod_LSNKpreind = param[np.random.random_integers(len(param))-1]
    param = ['BCC-CSM-11','CESM1-BGC','CanESM2','HadGEM2-ES','IPSL-CM5A-LR','MPI-ESM-LR','NorESM1-ME']
    mod_LSNKtrans = param[np.random.random_integers(len(param))-1]
    param = ['ESA-CCI','MODIS','Ramankutty1999','Levavasseur2012','CLM-45','JSBACH','JULES','LPJ','LPJ-GUESS','LPX-Bern','OCN','ORCHIDEE','VISIT']
    mod_LSNKcover = param[np.random.random_integers(len(param))-1]

    # wildfires
    param = ['','CLM-45','JSBACH','LPJ','LPJ-GUESS','ORCHIDEE','VISIT']
    mod_EFIREpreind = param[np.random.random_integers(len(param))-1]
    param = ['','CESM1-BGC','IPSL-CM5A-LR','MPI-ESM-LR','NorESM1-ME']
    mod_EFIREtrans = param[np.random.random_integers(len(param))-1]

    # permafrost
    param = ['','JSBACH','ORCHIDEE-MICT','JULES-DeepResp','JULES-SuppressResp']
    mod_EPFmain = param[np.random.random_integers(len(param))-1]
    param = ['zero','best','twice']
    mod_EPFmethane = param[np.random.random_integers(len(param))-1]

    # land-use
    param = ['CLM-45','LPJ-GUESS','ORCHIDEE']
    mod_ELUCagb = param[np.random.random_integers(len(param))-1]
    param = ['high','low']
    mod_EHWPbb = param[np.random.random_integers(len(param))-1]
    param = ['Houghton2001','Earles2012']
    mod_EHWPtau = param[np.random.random_integers(len(param))-1]
    param = ['normal','fast','slow']
    mod_EHWPspeed = param[np.random.random_integers(len(param))-1]

    # hydroxyl
    param = ['Prather2012','CESM-CAM-superfast','CICERO-OsloCTM2','CMAM','EMAC','GEOSCCM','GDFL-AM3','GISS-E2-R','GISS-E2-R-TOMAS','HadGEM2','LMDzORINCA','MIROC-CHEM','MOCAGE','NCAR-CAM-35','STOC-HadAM3','TM5','UM-CAM']
    mod_OHSNKtau = param[np.random.random_integers(len(param))-1]
    param = ['lin','log']
    mod_OHSNKfct = param[np.random.random_integers(len(param))-1]
    param = ['mean-OxComp','Holmes2013','GEOS-Chem','Oslo-CTM3','UCI-CTM']
    mod_OHSNKtrans = param[np.random.random_integers(len(param))-1]

    # wetlands
    param = ['','CLM-4Me','DLEM','IAP-RAS','LPJ-Bern','LPJ-WSL','ORCHIDEE','SDGVM']
    mod_EWETpreind = param[np.random.random_integers(len(param))-1]
    param = ['','CLM-4Me','DLEM','LPJ-Bern','ORCHIDEE','SDGVM','UVic-ESCM']
    mod_AWETtrans = param[np.random.random_integers(len(param))-1]

    # photolysis
    param = ['Prather2015','GMI','GEOSCCM','G2d-M','G2d','Oslo-c29','Oslo-c36','UCI-c29','UCI-c36']
    mod_HVSNKtau = param[np.random.random_integers(len(param))-1]
    param = ['Prather2012','Prather2015','G2d','Oslo-c29','UCI-c29']
    mod_HVSNKtrans = param[np.random.random_integers(len(param))-1]
    param = ['CAM-35','CMAM','Niwa-SOCOL','SOCOL','ULAQ','UMUKCA-UCAM']
    mod_HVSNKcirc = param[np.random.random_integers(len(param))-1]

    # ozone tropo
    param = ['','CAMCHEM','FRSGCUCI','GISS-modelE','GMI','INCA','LLNL-IMPACT','MOZART-GFDL','MOZECH','STOC-HadAM3','TM5-JRC','UM-CAM']
    mod_O3Tregsat = param[np.random.random_integers(len(param))-1]
    param = ['mean-OxComp','CICERO-OsloCTM2','NCAR-CAM-35','STOC-HadAM3','UM-CAM']
    mod_O3Temis = param[np.random.random_integers(len(param))-1]
    param = ['','CESM-CAM-superfast','GFDL-AM3','GISS-E2-R','MIROC-CHEM','MOCAGE','NCAR-CAM-35','STOC-HadAM3','UM-CAM']
    mod_O3Tclim = param[np.random.random_integers(len(param))-1]
    param = ['IPCC-AR5','IPCC-AR4','CESM-CAM-superfast','CICERO-OsloCTM2','CMAM','EMAC','GEOSCCM','GFDL-AM3','GISS-E2-R','GISS-E2-R-TOMAS','HadGEM2','LMDzORINCA','MIROC-CHEM','MOCAGE','NCAR-CAM-35','STOC-HadAM3','UM-CAM','TM5']
    mod_O3Tradeff = param[np.random.random_integers(len(param))-1]

    # ozone strato
    param = ['Newman2006','Laube2013']
    mod_O3Sfracrel = param[np.random.random_integers(len(param))-1]
    param = ['AMTRAC','CCSR-NIES','CMAM','CNRM-ACM','LMDZrepro','MRI','Niwa-SOCOL','SOCOL','ULAQ','UMSLIMCAT','UMUKCA-UCAM']
    mod_O3Strans = param[np.random.random_integers(len(param))-1]
    param = ['','Daniel2010']
    mod_O3Snitrous = param[np.random.random_integers(len(param))-1]
    param = ['IPCC-AR4','ULAQ','DLR-E39C','NCAR-MACCM','CHASER']
    mod_O3Sradeff = param[np.random.random_integers(len(param))-1]

    # sulfate
    param = ['','CAMCHEM','GISS-PUCCINI','GMI','GOCART','INCA2','LLNL-IMPACT','SPRINTARS']
    mod_SO4regsat = param[np.random.random_integers(len(param))-1]
    param = ['CSIRO-Mk360','GFDL-AM3','GISS-E2-R','MIROC-CHEM']
    mod_SO4load = param[np.random.random_integers(len(param))-1]
    param = ['BCC','CAM4-Oslo','CAM-51','GEOS-CHEM','GISS-MATRIX','GISS-modelE','GMI','GOCART','HadGEM2','IMPACT-Umich','INCA','MPIHAM','NCAR-CAM-35','OsloCTM2','SPRINTARS']
    mod_SO4radeff = param[np.random.random_integers(len(param))-1]
    
    # poa
    param = ['default','GFDL','CSIRO']
    mod_POAconv = param[np.random.random_integers(len(param))-1]
    param = ['','CAMCHEM','GISS-PUCCINI','GMI','GOCART','INCA2','LLNL-IMPACT','SPRINTARS']
    mod_POAregsat = param[np.random.random_integers(len(param))-1]
    param = ['CSIRO-Mk360','GFDL-AM3','GISS-E2-R','MIROC-CHEM']
    mod_POAload = param[np.random.random_integers(len(param))-1]
    param = ['BCC','CAM4-Oslo','CAM-51','GEOS-CHEM','GISS-MATRIX','GISS-modelE','GMI','GOCART','HadGEM2','IMPACT-Umich','INCA','MPIHAM','NCAR-CAM-35','OsloCTM2','SPRINTARS']
    mod_POAradeff = param[np.random.random_integers(len(param))-1]

    # bc
    param = ['','CAMCHEM','GISS-PUCCINI','GMI','GOCART','INCA2','LLNL-IMPACT','SPRINTARS']
    mod_BCregsat = param[np.random.random_integers(len(param))-1]
    param = ['CSIRO-Mk360','GFDL-AM3','GISS-E2-R','MIROC-CHEM']
    mod_BCload = param[np.random.random_integers(len(param))-1]
    param = ['BCC','CAM4-Oslo','CAM-51','GEOS-CHEM','GISS-MATRIX','GISS-modelE','GMI','GOCART','HadGEM2','IMPACT-Umich','INCA','MPIHAM','NCAR-CAM-35','OsloCTM2','SPRINTARS']
    mod_BCradeff = param[np.random.random_integers(len(param))-1]
    param = ['Boucher2013','CSIRO','GISS','HadGEM2','ECHAM5','ECMWF']
    mod_BCadjust = param[np.random.random_integers(len(param))-1]
     
    # nitrate
    param = ['Bellouin2011','Hauglustaine2014']
    mod_NO3load = param[np.random.random_integers(len(param))-1]
    param = ['GEOS-CHEM','GISS-MATRIX','GMI','HadGEM2','IMPACT-Umich','INCA','NCAR-CAM-35','OsloCTM2']
    mod_NO3radeff = param[np.random.random_integers(len(param))-1]

    # soa
    param = ['','GFDL-AM3','GISS-E2-R']
    mod_SOAload = param[np.random.random_integers(len(param))-1]
    param = ['CAM-51','GEOS-CHEM','IMPACT-Umich','MPIHAM','OsloCTM2']
    mod_SOAradeff = param[np.random.random_integers(len(param))-1]

    # dust
    param = ['CSIRO-Mk360','GFDL-AM3','GISS-E2-R','MIROC-CHEM']
    mod_DUSTload = param[np.random.random_integers(len(param))-1]
    param = [''] # [no value for now]
    mod_DUSTradeff = param[np.random.random_integers(len(param))-1]
   
    # salt
    param = ['GFDL-AM3','GISS-E2-R','MIROC-CHEM']
    mod_SALTload = param[np.random.random_integers(len(param))-1]
    param = [''] # [no value for now]
    mod_SALTradeff = param[np.random.random_integers(len(param))-1]
   
    # cloud
    param = ['Hansen2005','Lamarque2011']
    mod_CLOUDsolub = param[np.random.random_integers(len(param))-1]
    param = ['mean-ACCMIP','CSIRO-Mk360','GFDL-AM3','GISS-E2-R','HadGEM2','LMDzORINCA','MIROC-CHEM','NCAR-CAM-51']
    mod_CLOUDerf = param[np.random.random_integers(len(param))-1]    
    param = ['low','median','high']
    mod_CLOUDpreind = param[np.random.random_integers(len(param))-1]    

    # albedo bc
    param = ['Reddy2007']
    mod_ALBBCreg = param[np.random.random_integers(len(param))-1]
    param = ['CICERO-OsloCTM2','GFDL-AM3','GISS-E2-R','GISS-E2-R-TOMAS','HadGEM2','MIROC-CHEM','NCAR-CAM-35','NCAR-CAM-51']
    mod_ALBBCrf = param[np.random.random_integers(len(param))-1]    
    param = ['low','median','high']
    mod_ALBBCwarm = param[np.random.random_integers(len(param))-1]    

    # albedo lc
    param = ['CERES','GEWEX','MERRA']
    mod_ALBLCflux = param[np.random.random_integers(len(param))-1]    
    param = ['GlobAlbedo','MODIS']
    mod_ALBLCalb = param[np.random.random_integers(len(param))-1]
    param = ['ESA-CCI','MODIS']
    mod_ALBLCcover = param[np.random.random_integers(len(param))-1]    
    param = ['Hansen2005','Davin2007','Davin2010','Jones2013']
    mod_ALBLCwarm = param[np.random.random_integers(len(param))-1]    

    # temperature
    param = ['ACCESS-10','ACCESS-13','BCC-CSM-11','BCC-CSM-11m','CanESM2','CCSM4','CNRM-CM5','CNRM-CM5-2','CSIRO-Mk360','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H','GISS-E2-R','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MPI-ESM-P','MRI-CGCM3','NorESM1-M']
    mod_TEMPresp = param[np.random.random_integers(len(param))-1]
    param = ['4xCO2','hist&RCPs']
    mod_TEMPpattern = param[np.random.random_integers(len(param))-1]    
    
    # precipitations
    param = ['ACCESS-10','ACCESS-13','BCC-CSM-11','BCC-CSM-11m','CanESM2','CCSM4','CNRM-CM5','CNRM-CM5-2','CSIRO-Mk360','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H','GISS-E2-R','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MPI-ESM-P','MRI-CGCM3','NorESM1-M']
    mod_PRECresp = param[np.random.random_integers(len(param))-1]
    param = ['Andrews2010','Kvalevag2013']
    mod_PRECradfact = param[np.random.random_integers(len(param))-1]    
    param = ['4xCO2','hist&RCPs']
    mod_PRECpattern = param[np.random.random_integers(len(param))-1]    

    # acidification
    param = ['Tans2009','Bernie2010']
    mod_ACIDsurf = param[np.random.random_integers(len(param))-1]    
    
    # slr
    param = [''] # [no value for now]
    mod_SLR = param[np.random.random_integers(len(param))-1]    

    # SAVE DRIVERS
    drive_list.append([data_EFF,data_LULCC,\
                       data_ECH4,data_EN2O,data_Ehalo,\
                       data_ENOX,data_ECO,data_EVOC,\
                       data_ESO2,data_ENH3,data_EOC,data_EBC,\
                       data_RFant,data_RFnat])
    
    # SAVE PARAMETERS
    param_list.append([mod_OSNKstruct,mod_OSNKchem,mod_OSNKtrans,\
                       mod_LSNKnpp,mod_LSNKrho,mod_LSNKpreind,mod_LSNKtrans,mod_LSNKcover,\
                       mod_EFIREpreind,mod_EFIREtrans,\
                       mod_ELUCagb,mod_EHWPbb,mod_EHWPtau,mod_EHWPfct,\
                       mod_OHSNKtau,mod_OHSNKfct,mod_OHSNKtrans,\
                       mod_EWETpreind,mod_AWETtrans,\
                       mod_HVSNKtau,mod_HVSNKtrans,mod_HVSNKcirc,\
                       mod_O3Tregsat,mod_O3Temis,mod_O3Tclim,mod_O3Tradeff,\
                       mod_O3Sfracrel,mod_O3Strans,mod_O3Snitrous,mod_O3Sradeff,\
                       mod_SO4regsat,mod_SO4load,mod_SO4radeff,\
                       mod_POAconv,mod_POAregsat,mod_POAload,mod_POAradeff,\
                       mod_BCregsat,mod_BCload,mod_BCradeff,mod_BCadjust,\
                       mod_NO3load,mod_NO3radeff,\
                       mod_SOAload,mod_SOAradeff,\
                       mod_DUSTload,mod_DUSTradeff,\
                       mod_SALTload,mod_SALTradeff,\
                       mod_CLOUDsolub,mod_CLOUDerf,mod_CLOUDpreind,\
                       mod_ALBBCreg,mod_ALBBCrf,mod_ALBBCwarm,\
                       mod_ALBLCflux,mod_ALBLCalb,mod_ALBLCcover,mod_ALBLCwarm,\
                       mod_TEMPresp,mod_TEMPpattern,\
                       mod_PRECresp,mod_PRECradfact,mod_PRECpattern,\
                       mod_ACIDsurf,\
                       mod_SLR])

    # RELOAD
    execfile('OSCAR-loadD.py')
    execfile('OSCAR-loadP.py')
    execfile('OSCAR-format.py')
    execfile('OSCAR-fct.py')

    # RUN
    OUT = OSCAR_lite(var_output=['D_CO2','D_CH4','D_N2O','RF_halo','D_O3t','D_O3s','D_SO4','D_POA','D_BC','D_NO3','D_SOA','D_AERh','RF','D_gst'])

    # CHECK
    for n in range(len(OUT)):
        if np.any(np.isnan(OUT[n])):
            tag_break = tag_break|break_if_error
        else:
            tag_break = tag_break|False

    # BREAK
    if tag_break:
        break

# WRITE
writer = csv.writer(open('results/drive_test.csv','wb'))
writer.writerows(drive_list)
writer = csv.writer(open('results/param_test.csv','wb'))
writer.writerows(param_list)
writer = csv.writer(open('results/#OUT.USELESS.empty','wb'))

