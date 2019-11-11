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

    #fname_params_prefix = CUR_DIR + '/cellcycle.'
    fname_params_prefix = CUR_DIR + '/emt.'
    fname_models = CUR_DIR + '/models.txt'

    #name of the file containing steady states and cluster number:
    #fname_sstate = fname_params_prefix + 'states.unnormed.txt'
    fname_sstate = fname_params_prefix + 'states.clustered.txt'


    #import model ids from the input file:
    with open(fname_models, 'r') as fh_models:
        content = fh_models.read()
    models = content.strip().split("\n")

    #print(models)

    MODELS_TO_INSPECT = list(map(int, models))

    #import perturbed params from the params file as an ordered dictionary:
    #(parameters are perturbed in an R function)
    model_params_dict = ma.import_model_params(fname_params_prefix, \
                                            network, MODELS_TO_INSPECT)

    # Import starting state
    # ---------------------
    CLUSTER_NO=1
    sstate_dict = ma.import_states_byCluster(fname_sstate, \
                                             network, \
                                             MODELS_TO_INSPECT, 
                                             CLUSTER_NO)

    # Import starting state
    # ---------------------
    CLUSTER_NO=2
    estate_dict = ma.import_states_byCluster(fname_sstate, \
                                             network, \
                                             MODELS_TO_INSPECT, 
                                             CLUSTER_NO)

    # generate trajectory of the models from E to M
    # --------------------------------------------
    #mg.generate_e2m_traj_mpr_signal(network, model_params_dict, sstate_dict)
    #FCHANGE = 1.5 #2**1 #2**2 #2**3 #2**4 #2**5 #2**6 #2**7 #2**8 - phase.20
    #FCHANGE = 2**856 #2**136 #2**128 #2**64 #2**24 #2**16 #2**8 - phase 21

    SIGNALING_NODE_ID = 8
    FCHANGE = 2**8 #6.84 # phase 32
    mg.generate_traj_incr_sig_gradually(network, model_params_dict, sstate_dict, \
                                        SIGNALING_NODE_ID, FCHANGE)
    return None

#**********************************************************************#
if __name__=='__main__':
   start_time=time.time()
   main()
   os.chdir(CUR_DIR)
   print("--- %s seconds ---" % (time.time() - start_time))
