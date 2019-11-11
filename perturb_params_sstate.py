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
#import networkManager_ensemble as nm 
import networkManager_pertParams as nm 


#import modelGenerator as mg
import modelGenerator_pertParams as mg

#import modelAnalyzer_annealing as ma
#import ensembleAnalyzer_annealing as ma
import ensembleAnalyzer_pertParams as ma

from collections import defaultdict 
from collections import OrderedDict
import time
from os.path import basename 

import numpy as np

def are_seeds_valid(USER_SEED, SEED):
    if (int(USER_SEED)<=0 and int(SEED)<=0):
        return False 
    return True 

#======================================================================#
def main():
    # create an instance of Network class: 
    network=nm.Network(CUR_DIR) 

    config_dict=network.get_config_dict()

    SEED=int(config_dict['SEED'])
    USER_SEED=int(config_dict['USER_SEED'])

    #load network topology and parameter ranges:
    mode=network.process_network(USER_SEED, SEED)

    fname_params_prefix = CUR_DIR + '/cellcycle.'
    #fname_params_prefix = CUR_DIR + '/repressilator.'
    fname_models = CUR_DIR + '/models.txt'

    # name of the file containing steady states and cluster number:
    fname_sstate = fname_params_prefix + 'states.unnormed.txt'

    #import model ids from the input file:
    with open(fname_models, 'r') as fh_models:
        content = fh_models.read()
    models = content.strip().split("\n")

    MODELS_TO_INSPECT = list(map(int, models))

    #import perturbed params from the params file as an ordered dictionary:
    #(parameters are perturbed in an R function)
    model_params_dict = ma.import_model_params(fname_params_prefix, \
                                            network, MODELS_TO_INSPECT)

    #import steady state expressions of the current cluster:
    sstate_dict = ma.import_steady_states(fname_sstate, \
                                          network, \
                                          MODELS_TO_INSPECT)

    # generate models using perturbed parameters and random ICs:
    # (number of random ICs can be set in racipe.cfg file)
    #mg.generate_models_random_ICs(network, model_params_dict)

    # generate models using perturbed parameters and steady state as IC:
    mg.generate_models_sstate(network, model_params_dict, sstate_dict)

    return None


#**********************************************************************#
if __name__=='__main__':
   start_time=time.time()
   main()
   os.chdir(CUR_DIR)
   print("--- %s seconds ---" % (time.time() - start_time))
