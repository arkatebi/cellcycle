#/usr/bin/env python
''''
    racipe: this program implements the RACIPE method to generate an 
ensemble of network models from a given gene network topology. 
    How to run this tool?
    Mode 1: topology mode-supply only topology file-
       > python probeModel.py -M='T' \
                          -I1=TS.tpo \
     Mode 2: parameter mode-supply parameter file as well-
       > python probeModel.py -M='C' \
                          -I1=TS.tpo \
                          -I2=TS.prs \
'''
import os
import sys

#after saving the current directory, move 
#to the program directory:

CUR_DIR=os.getcwd()
PROG_PATH=os.path.dirname(os.path.abspath(sys.argv[0]))
os.chdir(PROG_PATH)

import argParserCal as ap
import networkManager as nm 
import modelGenerator as mg 
import modelAnalyzer as ma
import configHandler as ch 

from collections import defaultdict 
from collections import OrderedDict
import time
from os.path import basename 

#------------------------------------------------------------------#
def build_network(fname_tpo):
    node_dict,edge_dict,source_dict,target_dict,master_dict=\
        nm.upload_topology(open(fname_tpo,'r'),config_dict,
                        node_dict,source_dict,
                        target_dict,master_dict) 

    sys.exit(0)
    node_id_dict,edge_source_dict,edge_target_dict,edge_type_dict=\
        nm.build_map_dict(node_dict,master_dict,node_id_dict,
                       edge_source_dict,edge_target_dict,
                       edge_type_dict)
    return None 




#======================================================================#
def main():
    import ctypes
    clib=ctypes.cdll.LoadLibrary('./analyzer_clib.so')

    model_to_inspect=2

    parsed_dict=ap.parse_args()
    print("here:")
    print("tpo file name: ",parsed_dict['fname_tpo'])
    print("EXP file name: ",parsed_dict['fname_exp'])
    print("TSH file name: ",parsed_dict['fname_params'])

    tpo_fname=parsed_dict['fname_tpo']

    tpo_path=os.path.dirname(parsed_dict['fname_tpo'])

    print(tpo_path)
    #fn_prefix=(basename(tpo_fname)).strip().\
    #    split(self.config_dict['TOPOLOGY_FNAME_EXTENSION'])[0]

    # Get .tpo filename with absolute path:
    print(CUR_DIR)
    if (not tpo_path or tpo_path=="."): 
        tpo_path=CUR_DIR
    tpo_fname=tpo_path+'/'+tpo_fname

    tmp_fname=parsed_dict['fname_exp']
    exp_fname=CUR_DIR+'/'+parsed_dict['fname_exp']
    params_fname=CUR_DIR+'/'+parsed_dict['fname_params'] 

    CONFIG_FNAME=CUR_DIR+'/'+'racipe.cfg'
    config_dict = ch.read_config(CONFIG_FNAME) 
    build_network(tpo_fname,config_dict)

    sys.exit(0)
    #tpo_fname=self.parsed_dict['tpo_fname']
 
    #create an instance of Network class: 
    network=nm.Network(CUR_DIR)


    #load network topology and parameter ranges:
    network.build_network()
    
    #obtain racipe config:
    config_dict=network.get_config_dict()
    SEED=int(config_dict['SEED'])
    user_seed=int(config_dict['USER_SEED'])
    #set the seed for model generation:
    clib.set_seed(ctypes.c_int(SEED*2),ctypes.c_int(user_seed))

    ma.inspect_model(network,model_to_inspect)
    sys.exit(0)

    SEED=int(config_dict['SEED'])
    user_seed=int(config_dict['USER_SEED'])

    clib=ctypes.cdll.LoadLibrary('./simulation_clib.so')
    clib.randu.argtypes=(ctypes.c_double,ctypes.c_double)
    clib.randu.restype=ctypes.c_double
    #clib.set_seed(ctypes.c_int(SEED))
    #set seed for estimation of parameter ranges:
    clib.set_seed(ctypes.c_int(SEED),ctypes.c_int(user_seed))

    mode=network.process_network()

    #generate models:
    if(mode=='A' or mode=='P'):
        print('Generating models ...')
        #set seed for generating models:
        #clib.set_seed(ctypes.c_int(SEED*2))
        clib.set_seed(ctypes.c_int(SEED*2),ctypes.c_int(user_seed))
        #mg.generate_models_usingPython(network) 
        mg.generate_models_usingC(network)
    return None 


#**********************************************************************#
if __name__=='__main__':
   start_time=time.time()
   main()
   os.chdir(CUR_DIR)
   print("--- %s seconds ---" % (time.time() - start_time))
