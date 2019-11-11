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
import networkManager_exp as nm 

from collections import defaultdict 
from collections import OrderedDict
import time
from os.path import basename 

print(os.getcwd())
#import exprCalculator as epcr 
import modelGenerator as mg


def are_seeds_valid(USER_SEED, SEED):
    if (int(USER_SEED)<=0 and int(SEED)<=0):
        return False 
    return True 

#======================================================================#
def main():
    import ctypes

    # create an instance of Network class: 
    network=nm.Network(CUR_DIR) 

    config_dict=network.get_config_dict()

    SEED=int(config_dict['SEED'])
    USER_SEED=int(config_dict['USER_SEED'])

    #print(USER_SEED, "\t", SEED)
    if (not are_seeds_valid(USER_SEED,SEED)):
        print("both user seed(US) and internal seed(S) cannot be non-positive.")
        print("exiting ...")
        sys.exit(0)

    #clib=ctypes.cdll.LoadLibrary('./expr_sim.so')
    clib=ctypes.cdll.LoadLibrary('./simulation_clib.so')
    clib.randu.argtypes=(ctypes.c_double,ctypes.c_double)
    clib.randu.restype=ctypes.c_double

    mode=network.process_network(USER_SEED, SEED)

    print('Generating TF activities ...')

    clib.set_seed(ctypes.c_int(USER_SEED), ctypes.c_int(SEED*2))
    #epcr.generate_expressions(network)
    mg.cal_expressions(network)

    return None 


#**********************************************************************#
if __name__=='__main__':
   start_time=time.time()
   main()
   os.chdir(CUR_DIR)
   print("--- %s seconds ---" % (time.time() - start_time))
