#!/usr/bin/env python
'''
    This module has the following methods:
'''
import sys
import os
import numpy as np
from collections import OrderedDict
from collections import defaultdict
from os.path import basename
import copy
import time
import ctypes
np.random.seed(1)

clib=ctypes.cdll.LoadLibrary('./simulation_clib.so')

#clib.randu.argtypes=(ctypes.c_double,ctypes.c_double)
#clib.randu.restype=ctypes.c_double

SLOW_EDGE_TYPES=[1,2,3,4]
FAST_EDGE_TYPES=[5,6]

EXCITATION_TYPES=[1,4,5] #4: degrdation inhibition => activation of degradation
INHIBITION_TYPES=[2,3,6] #3: dedgration activation => inhibition of degradation

#-----------------------------------------------------------------------------#
def print_dict_log2(d): 
    for k,v in sorted(d.items()): 
        print(k, ': ', np.log2(v))
    print('\n')
    return None

#-----------------------------------------------------------------------------#
def print_dict(d): 
    for k,v in d.items():
        print(k, ': ', v)
    print('\n')
    return None

#-----------------------------------------------------------------------------#
def print_nested_defaultdict(node_dict,d): 
    for k in node_dict.keys(): 
        v= d[k]
        print(k,':')
        for idx,e in v.items(): 
            print('  ',idx,': ', e)
    print('\n')
    return None

#-----------------------------------------------------------------------------#
def print_nested_defaultdict_log2(d): 
    #print(d)
    #sys.exit(0)
    for idx,v in d.items(): 
        print(idx,':')
        for X,e in v.items(): 
            print('  ',X,': ', np.log2(e))
    print('\n')
    return None

#-----------------------------------------------------------------------------#
def save_parameters(fh_params,model_no,nodeParam_dict,
                      master_dict,solution_dict):
    #model no + number of states:
    outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    #skipping node type information nodeParam_dict[X][0]
    #save MPRs:
    for X in nodeParam_dict.keys():
        #outstr+='\t'+str('%10.6f'%nodeParam_dict[X][0])
        outstr+='\t'+str('%10.6f'%nodeParam_dict[X][1])
    outstr+='\n'
    fh_params.write(outstr)
    #save all DNRs:
    outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    for X in nodeParam_dict.keys():
        #outstr+='\t'+str('%10.6f'%nodeParam_dict[X][1])
        outstr+='\t'+str('%10.6f'%nodeParam_dict[X][2])
    outstr+='\n'
    fh_params.write(outstr)

    #save all TSHs:
    outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    for idx in master_dict.keys():
        outstr+='\t'+str('%10.6f'%master_dict[idx][3])
    outstr+='\n'
    fh_params.write(outstr)
    #save all HCOs:
    outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    for idx in master_dict.keys():
        outstr+='\t'+str(master_dict[idx][4])
    outstr+='\n'
    fh_params.write(outstr)
    #save all FCHs:
    outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    for idx in master_dict.keys():
        outstr+='\t'+str('%10.6f'%master_dict[idx][5])
    outstr+='\n'
    fh_params.write(outstr)

    fh_params.flush()
    return None

def save_parameters_old(fh_params,model_no,nodeParam_dict,
                      master_dict,solution_dict):
    #model no + number of states:
    outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    #skipping node type information nodeParam_dict[X][0]
    #write MPR to outstr:
    for X in nodeParam_dict.keys():
        #outstr+='\t'+str('%10.6f'%nodeParam_dict[X][0])
        outstr+='\t'+str('%10.6f'%nodeParam_dict[X][1])
    #write DNR to outstr:
    for X in nodeParam_dict.keys():
        #outstr+='\t'+str('%10.6f'%nodeParam_dict[X][1])
        outstr+='\t'+str('%10.6f'%nodeParam_dict[X][2])
    #write edge parameters (TSH,HCO,FCH) to outstr:
    for idx in master_dict.keys():
        outstr+='\t'+str('%10.6f'%master_dict[idx][3])
        outstr+='\t'+str(master_dict[idx][4])
        outstr+='\t'+str('%10.6f'%master_dict[idx][5])
    outstr+='\n'
    fh_params.write(outstr)
    fh_params.flush()
    return None

#----------------------------------------------------------------------#
def read_params(fh_params,network,MODEL_TO_INSPECT):
    '''
    This method loads the network topology and places the topology
    information into the specific data structures. 
    fh_tpo: file handle to a network topology file
    It updates the following data structures: 
    node_dict: a dictionary of genes/nodes in the network
    target_dict: 
    source_dict: 
    master_dict: 
    '''
    #
    node_id_dict=network.get_node_id_dict()
    edge_source_dict=network.get_edge_source_dict()
    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

    MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #import node and edge information:
    node_dict=network.get_node_dict()
    master_dict=network.get_master_dict()
    node_count=len(node_dict.keys())
    edge_count=len(master_dict.keys())

    for line in fh_params:
        if not line.strip():
            continue 
        #read line 1:MPR
        #fields[0]:model index,fields[1]:number of states
        #fields[2] ~ fields[.]:node_count number of MPRs 
        fields=line.strip().split()
        if (int(fields[0]) is not int(MODEL_TO_INSPECT)): 
            #for the current model, skip all subsequent lines:
            next(fh_params) #skip DNR line
            next(fh_params) #skip TSH line
            next(fh_params) #skip HCO line
            next(fh_params) #skip FCH line
            continue

        for i in range(2,len(fields)):
            MPR_arr[i-2]=fields[i]

        if(node_count is not len(fields[2:])): 
            print('mismatch in node count and MPR count')
            print('program exiting ...')
            sys.exit(0)

        #read line 2:DNR
        line2=fh_params.readline()
        fields=line2.strip().split()
        for i in range(2,len(fields)):
            DNR_arr[i-2]=fields[i]

        #read line 3:TSH
        line3=fh_params.readline()
        fields=line3.strip().split()
        for i in range(2,len(fields)):
            TSH_arr[i-2]=fields[i]

        if(edge_count is not len(fields[2:])):
            print('mismatch in edge count and TSH count')
            print('program exiting ...')
            sys.exit(0)
 
        #read line 4:HCO 
        line4=fh_params.readline()
        fields=line4.strip().split()
        for i in range(2,len(fields)):
            HCO_arr[i-2]=fields[i]

        #read line 5:FCH 
        line5=fh_params.readline()
        fields=line5.strip().split()
        for i in range(2,len(fields)):
            FCH_arr[i-2]=fields[i]
    return MPR_arr,DNR_arr,TSH_arr,HCO_arr,FCH_arr


#----------------------------------------------------------------------#
def import_model_params(fh_params, network, MODEL_LIST):
    '''
    This method imports the parameters for each model in 
    the MODEL_LIST and places them in a dictionary  model_params_dict

    model_params_dict: 
    key: model no 
    values: a list of five lists for each model: 
            MPR list: a list of MPR parameters
            DNR list: a list of DNR parameters
            TSH list: a list of TSH parameters 
            HCO list: a list of HCO parameters 
            FCH list: a list of FCH parameters
    '''
    # define a dictionary to store the parameters for the model set:
    model_params_dict = OrderedDict()

    # convert the list elements to int if they are not in int format already:
    MODEL_LIST = list(map(int, MODEL_LIST))

    node_id_dict=network.get_node_id_dict()
    edge_source_dict=network.get_edge_source_dict()
    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

    #MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    #DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    #NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    #TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    #HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    #FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #import node and edge information: 
    #need this information for verification of number of fields in different lines:
    node_dict=network.get_node_dict()
    master_dict=network.get_master_dict()
    node_count=len(node_dict.keys())
    edge_count=len(master_dict.keys())

    for line in fh_params:
        if not line.strip():
            continue 
        #read line 1:MPR
        fields=line.strip().split()

        #fields[0]:model index,fields[1]:number of states
        #fields[2] ~ fields[.]:node_count number of MPRs 

        if (int(fields[0]) not in MODEL_LIST): 
            #for the current model, skip all subsequent lines:
            next(fh_params) #skip DNR line
            next(fh_params) #skip TSH line
            next(fh_params) #skip HCO line
            next(fh_params) #skip FCH line
            continue
        
        # initialize the lists for all parameter types:
        MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
        DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
        NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

        TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
        HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
        FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

        for i in range(2,len(fields)):
            MPR_arr[i-2]=fields[i]

        if(node_count is not len(fields[2:])): 
            print('mismatch in node count and MPR count')
            print('program exiting ...')
            sys.exit(0)

        #read line 2:DNR
        line2=fh_params.readline()
        fields=line2.strip().split()
        for i in range(2,len(fields)):
            DNR_arr[i-2]=fields[i]

        #read line 3:TSH
        line3=fh_params.readline()
        fields=line3.strip().split()
        for i in range(2,len(fields)):
            TSH_arr[i-2]=fields[i]

        if(edge_count is not len(fields[2:])):
            print('mismatch in edge count and TSH count')
            print('program exiting ...')
            sys.exit(0)
 
        #read line 4:HCO 
        line4=fh_params.readline()
        fields=line4.strip().split()
        for i in range(2,len(fields)):
            HCO_arr[i-2]=fields[i]

        #read line 5:FCH 
        line5=fh_params.readline()
        fields=line5.strip().split()
        for i in range(2,len(fields)):
            FCH_arr[i-2]=fields[i]

        # insert the params in the dictionary: 
        model_params_dict[int(fields[0])] = [MPR_arr, DNR_arr, TSH_arr, HCO_arr, FCH_arr]

    return model_params_dict

#----------------------------------------------------------------------#
def emulate_a_model(network, model_no, params_list): 
    '''
    This method generates RACIPE models. 
    '''

    global clib

    import numpy
    import numpy as np
    import ctypes
    import sys

    import modelManager as mm

    #sys.exit(0)
    # add appropriate noise: 
    NOISE=1.0 
    #NOISE=10.0 # 30.0 

    NOISE_SHOT=0

    node_dict = network.get_node_dict()
    node_id_dict=network.get_node_id_dict()

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

    #create arrays for storing noise level for each node:
    NOISE_strength=np.zeros(NUM_NODES,dtype=np.double)
    NOISE_strength_shot=np.zeros(NUM_NODES,dtype=np.double)


    #create array for storing the source nodes of the edges:
    edge_source_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_source_dict.keys():
        edge_source_arr[i]=edge_source_dict[idx]
        i+=1

    #create array for storing the target nodes of the edges:
    edge_target_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_target_dict.keys():
        edge_target_arr[i]=edge_target_dict[idx]
        i+=1
    #create character array for work directory string:
    WORK_DIR=(network.get_work_dir()).encode('utf-8')

    #create character arrays for file names to be used in C:
    fname_dict_simu = network.get_fname_dict_simu()
    FNAME_STATES = fname_dict_simu['FNAME_STATES'].encode('utf-8')
    FNAME_LIMITCYCLES = fname_dict_simu['FNAME_LIMITCYCLES'].encode('utf-8')
    FNAME_SUMMARY = fname_dict_simu['FNAME_SUMMARY'].encode('utf-8')

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    #create array for storing the source nodes of the edges:
    edge_source_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    i=0
    for idx in edge_source_dict.keys():
        edge_source_arr[i]=edge_source_dict[idx]
        i+=1

    #create array for storing the target nodes of the edges:
    edge_target_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_target_dict.keys():
        edge_target_arr[i]=edge_target_dict[idx]
        i+=1

    #create array for storing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    # create arrays for storing the parameters:
    MPR_arr = params_list[0]
    DNR_arr = params_list[1]

    TSH_arr = params_list[2]
    HCO_arr = params_list[3]
    FCH_arr = params_list[4]

    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)
    node_dict=network.get_node_dict()
    i=0
    for node in node_dict.keys():
        (node_type,mpr_range,dnr_range)=node_dict[node]
        NODE_TYPE_arr[i]=node_type

        NOISE_strength[i]=1
        NOISE_strength_shot[i]=1
        i=i+1

    #import config_dict:
    config_dict=network.get_config_dict()

    EXP_dict_arr=np.zeros(NUM_NODES*int(config_dict['NUM_RANDOM_ICS']),
                          dtype=np.double) 
    clib.find_solutions_stochastic(\
         ctypes.create_string_buffer(WORK_DIR),
         #ctypes.create_string_buffer(FNAME_STATES),
         #ctypes.create_string_buffer(FNAME_LIMITCYCLES),
         #ctypes.create_string_buffer(FNAME_SUMMARY),
         ctypes.c_int(model_no),
         ctypes.c_int(NUM_NODES),
         ctypes.c_int(NUM_EDGES),
         ctypes.c_int(int(config_dict['NUM_RANDOM_ICS'])),
         ctypes.c_int(int(config_dict['ITER_FOR_ODE'])),
         ctypes.c_double(float(config_dict['EULER_SIM_TIME'])),
         ctypes.c_double(float(config_dict['EULER_SIM_STEP_SIZE'])),
         ctypes.c_longdouble(float(config_dict['CONVERGENCE_PROXIMITY'])),
         ctypes.c_double(float(config_dict['TRANS_RATE_FACTOR'])),
         ctypes.c_void_p(MPR_arr.ctypes.data), 
         ctypes.c_void_p(DNR_arr.ctypes.data), 
         ctypes.c_void_p(NODE_TYPE_arr.ctypes.data),
         ctypes.c_void_p(edge_source_arr.ctypes.data), 
         ctypes.c_void_p(edge_target_arr.ctypes.data), 
         ctypes.c_void_p(edge_type_arr.ctypes.data), 
         ctypes.c_void_p(TSH_arr.ctypes.data), 
         ctypes.c_void_p(HCO_arr.ctypes.data), 
         ctypes.c_void_p(FCH_arr.ctypes.data), 
         ctypes.c_void_p(EXP_dict_arr.ctypes.data),
         ctypes.c_void_p(NOISE_strength.ctypes.data), 
         ctypes.c_void_p(NOISE_strength_shot.ctypes.data), 
         ctypes.c_double(NOISE),
         ctypes.c_double(NOISE_SHOT)
         )
    return None

#----------------------------------------------------------------------#
def simulate_network_cnoise(network, model_no, params_list): 
    '''
    This method generates RACIPE models. 
    '''

    global clib

    import numpy
    import numpy as np
    import ctypes
    import sys
    import modelManager as mm

    node_dict = network.get_node_dict()
    node_id_dict=network.get_node_id_dict()

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

    #create arrays for storing noise level for each node:
    NOISE_strength=np.zeros(NUM_NODES,dtype=np.double)
    NOISE_strength_shot=np.zeros(NUM_NODES,dtype=np.double)

    #create array for storing the source nodes of the edges:
    edge_source_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_source_dict.keys():
        edge_source_arr[i]=edge_source_dict[idx]
        i+=1

    #create array for storing the target nodes of the edges:
    edge_target_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_target_dict.keys():
        edge_target_arr[i]=edge_target_dict[idx]
        i+=1
    #create character array for work directory string:
    WORK_DIR=(network.get_work_dir()).encode('utf-8')

    #create character arrays for file names to be used in C:
    fname_dict_simu = network.get_fname_dict_simu()
    #FNAME_STATES = fname_dict_simu['FNAME_STATES'].encode('utf-8')
    FNAME_STATES = fname_dict_simu['FNAME_STATES_CNOISE'].encode('utf-8')
    FNAME_LIMITCYCLES = fname_dict_simu['FNAME_LIMITCYCLES'].encode('utf-8')
    FNAME_SUMMARY = fname_dict_simu['FNAME_SUMMARY'].encode('utf-8')

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    #create array for storing the source nodes of the edges:
    edge_source_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    i=0
    for idx in edge_source_dict.keys():
        edge_source_arr[i]=edge_source_dict[idx]
        i+=1

    #create array for storing the target nodes of the edges:
    edge_target_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_target_dict.keys():
        edge_target_arr[i]=edge_target_dict[idx]
        i+=1

    #create array for storing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    # create arrays for storing the parameters:
    MPR_arr = params_list[0]
    DNR_arr = params_list[1]

    TSH_arr = params_list[2]
    HCO_arr = params_list[3]
    FCH_arr = params_list[4]

    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)
    node_dict=network.get_node_dict()
    i=0
    for node in node_dict.keys():
        (node_type,mpr_range,dnr_range)=node_dict[node]
        NODE_TYPE_arr[i]=node_type
        i=i+1

    #import config_dict:
    config_dict=network.get_config_dict()

    EXP_dict_arr=np.zeros(NUM_NODES*int(config_dict['NUM_RANDOM_ICS']),
                          dtype=np.double) 

    clib.simulate_network_cnoise_gnw(\
         ctypes.create_string_buffer(WORK_DIR),
         ctypes.create_string_buffer(FNAME_STATES),
         ctypes.create_string_buffer(FNAME_LIMITCYCLES),
         ctypes.create_string_buffer(FNAME_SUMMARY),
         ctypes.c_int(model_no),
         ctypes.c_int(NUM_NODES),
         ctypes.c_int(NUM_EDGES),
         ctypes.c_int(int(config_dict['NUM_RANDOM_ICS'])),
         ctypes.c_int(int(config_dict['ITER_FOR_ODE'])),
         ctypes.c_double(float(config_dict['EULER_SIM_TIME'])),
         ctypes.c_double(float(config_dict['EULER_SIM_STEP_SIZE'])),
         ctypes.c_longdouble(float(config_dict['CONVERGENCE_PROXIMITY'])),
         ctypes.c_double(float(config_dict['TRANS_RATE_FACTOR'])),
         ctypes.c_void_p(MPR_arr.ctypes.data), 
         ctypes.c_void_p(DNR_arr.ctypes.data), 
         ctypes.c_void_p(NODE_TYPE_arr.ctypes.data),
         ctypes.c_void_p(edge_source_arr.ctypes.data), 
         ctypes.c_void_p(edge_target_arr.ctypes.data), 
         ctypes.c_void_p(edge_type_arr.ctypes.data), 
         ctypes.c_void_p(TSH_arr.ctypes.data), 
         ctypes.c_void_p(HCO_arr.ctypes.data), 
         ctypes.c_void_p(FCH_arr.ctypes.data), 
         ctypes.c_void_p(EXP_dict_arr.ctypes.data)
         )
    return None


#-----------------------------------------------------------------------------#
if __name__ == '__main__':
    print (sys.argv[0] + ':')
    print(__doc__)
    sys.exit(0)
