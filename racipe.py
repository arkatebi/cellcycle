#/usr/bin/env python
''''
    racipe: this program implements the RACIPE method to generate an 
ensemble of network models from a given gene network topology. 
    How to run this tool?
    Mode 1: topology mode-supply only topology file-
       > python racipe.py -M='T' \
                          -I1=TS.tpo \
     Mode 2: parameter mode-supply parameter file as well-
       > python racipe.py -M='C' \
                          -I1=TS.tpo \
                          -I2=TS.prs \
'''
import os
import sys

CUR_DIR=os.getcwd()
PROG_PATH=os.path.dirname(os.path.abspath(sys.argv[0]))
os.chdir(PROG_PATH)
import networkManager as nm 
import modelGenerator as mg 
import modelGenerator_Stochastic as mgs 
import modelGenerator_Stochastic_Annealing as mgsa 

import dataAnalyzer as da 

from collections import defaultdict 
from collections import OrderedDict
import time
from os.path import basename 

def are_seeds_valid(USER_SEED, SEED):
    if (int(USER_SEED)<=0 and int(SEED)<=0):
        return False 
    return True 

#======================================================================#
def main():
    import ctypes

    #parsed_dict=ap.parse_args()
 
    # create an instance of Network class: 
    network=nm.Network(CUR_DIR) 
    #print("here:")
    #sys.exit(0)

    config_dict=network.get_config_dict()

    SEED=int(config_dict['SEED'])
    USER_SEED=int(config_dict['USER_SEED'])

    #print(USER_SEED, "\t", SEED)
    if (not are_seeds_valid(USER_SEED,SEED)):
        print("both user seed(US) and internal seed(S) cannot be non-positive.")
        print("exiting ...")
        sys.exit(0)

    clib=ctypes.cdll.LoadLibrary('./simulation_clib.so')
    clib.randu.argtypes=(ctypes.c_double,ctypes.c_double)
    clib.randu.restype=ctypes.c_double
    #clib.set_seed(ctypes.c_int(SEED))
    #set seed for estimation of parameter ranges:
    #clib.set_seed(ctypes.c_int(SEED),ctypes.c_int(USER_SEED))
    #clib.set_seed(ctypes.c_int(USER_SEED), ctypes.c_int(SEED))
    #load network topology and parameter ranges:
    mode=network.process_network(USER_SEED, SEED)

    #generate models:
    if(mode=='A' or mode=='P'):
        print('Generating models ...')
        #set seed for generating models:
        #clib.set_seed(ctypes.c_int(SEED*2))
        #user supplied value for USER_SEED supercedes 
        #default value in racipe.cfg:
        #clib.set_seed(ctypes.c_int(SEED*2),ctypes.c_int(USER_SEED))
        clib.set_seed(ctypes.c_int(USER_SEED), ctypes.c_int(SEED*2))
        mg.generate_models(network)

        #if(mode=='A'): 
        #    network.delete_intermediate_files()
        #else:    
        #    network.delete_intermediate_files_2()

    #generate models - stochastic:
    # simulation when both I1 and I2 are supplied:
    if(mode=='SCP'):
        print('Generating models ...')
        #set seed for generating models:
        #clib.set_seed(ctypes.c_int(SEED*2))
        #user supplied value for USER_SEED supercedes 
        #default value in racipe.cfg:
        #clib.set_seed(ctypes.c_int(SEED*2),ctypes.c_int(USER_SEED))
        clib.set_seed(ctypes.c_int(USER_SEED), ctypes.c_int(SEED*2))
        mgs.generate_models(network)
        #network.delete_intermediate_files()
 
    # simulation when only I1 is supplied:
    if(mode=='SC'):
        print('Generating models ...')
        #set seed for generating models:
        #clib.set_seed(ctypes.c_int(SEED*2))
        #user supplied value for USER_SEED supercedes 
        #default value in racipe.cfg:
        #clib.set_seed(ctypes.c_int(SEED*2),ctypes.c_int(USER_SEED))
        clib.set_seed(ctypes.c_int(USER_SEED), ctypes.c_int(SEED*2))
        mgs.generate_models(network)
        #network.delete_intermediate_files()
        
    # simulated annealing when only I1 is supplied:
    if(mode=='SA'):
        print('Generating models ...')
        #set seed for generating models:
        #clib.set_seed(ctypes.c_int(SEED*2))
        #user supplied value for USER_SEED supercedes 
        #default value in racipe.cfg:
        #clib.set_seed(ctypes.c_int(SEED*2),ctypes.c_int(USER_SEED))
        clib.set_seed(ctypes.c_int(USER_SEED), ctypes.c_int(SEED*2))
        mgsa.generate_models(network)
        #network.delete_intermediate_files()
    
    # simulated annealing when both I1 and I2 are supplied:
    if(mode=='SAP'):
        print('Generating models ...')
        #set seed for generating models:
        #clib.set_seed(ctypes.c_int(SEED*2))
        #user supplied value for USER_SEED supercedes 
        #default value in racipe.cfg:
        #clib.set_seed(ctypes.c_int(SEED*2),ctypes.c_int(USER_SEED))
        clib.set_seed(ctypes.c_int(USER_SEED), ctypes.c_int(SEED*2))
        mgsa.generate_models(network)
        #network.delete_intermediate_files()

    return None 


#**********************************************************************#
if __name__=='__main__':
   start_time=time.time()
   main()
   os.chdir(CUR_DIR)
   print("--- %s seconds ---" % (time.time() - start_time))
