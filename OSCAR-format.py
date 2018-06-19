"""
Copyright: IIASA (International Institute for Applied Systems Analysis), 2016-2018; CEA (Commissariat a L'Energie Atomique) & UVSQ (Universite de Versailles et Saint-Quentin), 2016
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at)

This software is a computer program whose purpose is to simulate the behavior of the Earth system, with a specific but not exclusive focus on anthropogenic climate change.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""


##################################################
##################################################
##################################################


import numpy as np


##################################################
#   A. FORMAT DRIVERS
##################################################
print 'FORMATING'

# remove attribution axis
# drivers
for VAR in ['EFF','ECH4','EN2O']+['LUC','HARV','SHIFT']+['EHFC','EPFC','EODS']+['ENOX','ECO','EVOC','ESO2','ENH3','EOC','EBC']+['RFcon','RFvolc','RFsolar']:
    exec(VAR+' = np.sum(np.sum(np.sum('+VAR+',3),2),1)')
# parameters
for VAR in ['ECH4','EN2O']+['ENOX','ECO','EVOC','ESO2','ENH3','EOC','EBC']:
    exec(VAR+'_0 = np.sum(np.sum(np.sum('+VAR+'_0,2),1),0)')

