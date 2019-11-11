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
import dataAnalyzer as da 
import modelAnalyzer_annealing as ma

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

    fname_params = CUR_DIR + '/cellcycle.params.txt'
    #fname_models = CUR_DIR + '/model_list.txt'
    fname_models = CUR_DIR + '/model_list.10.10.txt'
    fname_state_med = CUR_DIR + '/cellcycle.states.med.csv'

    #fname_params = CUR_DIR + '/EMT22.params.txt'
    #fname_models = CUR_DIR + '/model_list.txt'
    #fname_state_med = CUR_DIR + '/EMT22.states.med.csv'

    #import median expressions of the states: 
    fh_med = open(fname_state_med, 'r') 
    next(fh_med) 
    line = fh_med.readline()
    fields = line.strip().split(',') 
    fields.pop(0)
    node_exp_med = list(map(float, fields))
    
    #convert median expressions from log scale to original expression levels:
    node_exp_med = np.power(2,node_exp_med) 

    #import model ids from the input file:
    with open(fname_models, 'r') as fh_models: 
        content = fh_models.read()
    models = content.strip().split("\n")

    MODELS_TO_INSPECT = list(map(int, models))

    #import params from the params file as an ordered dictionary:
    model_params_dict = ma.import_model_params(open(fname_params, "r"), \
                                            network, MODELS_TO_INSPECT)

    if (len(set(MODELS_TO_INSPECT)) == len(set(model_params_dict.keys()))):
        print("verified model ids ...")
        
    print("running simulated annealing ...")
    for model_id in list(model_params_dict.keys()):
        #print(model_id)
        params_list = model_params_dict[model_id]
        #ma.emulate_a_model(network, model_id, params_list)
        #ma.anneal_network_stochastic(network, model_id, params_list)
        ma.anneal_network_stochastic_gnw(network, model_id, params_list)
        #ma.anneal_network_stochastic(network, model_id, 
        #        params_list, node_exp_med)
        #break
    return None 


#**********************************************************************#
if __name__=='__main__':
   start_time=time.time()
   main()
   os.chdir(CUR_DIR)
   print("--- %s seconds ---" % (time.time() - start_time))
