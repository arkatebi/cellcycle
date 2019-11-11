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
import ctypes

#clib=ctypes.cdll.LoadLibrary('./simulation_clib.so')
clib=ctypes.cdll.LoadLibrary('./pertParamSim.so')
clib.randu.argtypes=(ctypes.c_double,ctypes.c_double)
clib.randu.restype=ctypes.c_double

SLOW_EDGE_TYPES=[1,2,3,4]
FAST_EDGE_TYPES=[5,6]

EXCITATION_TYPES=[1,4,5] #4: degrdation inhibition => activation of degradation
INHIBITION_TYPES=[2,3,6] #3: dedradation activation => inhibition of degradation

ZERO_PRODUCTION_RATE = 0.000001

MPR_MIN_KD = 1.0 
MPR_MAX_KD = 5.0

DNR_MIN_KD = 0.1 
DNR_MAX_KD = 1.0

#-----------------------------------------------------------------------------#
def set_parameters(network):
    global clib
    #nodeParam_dict: node to (MPR,DNR) mapping:
    #key: node, value: tuple (MPR,DNR)
    nodeParam_dict=OrderedDict()

    arg_dict = network.get_parsed_dict()

    if (arg_dict['knockdownlist']):  
        knock_down_list = arg_dict['knockdownlist'][0].split(",")
    else: 
        knock_down_list = arg_dict['knockdownlist']

    #import config_dict, node_dict, source_dict, target_dict:
    config_dict=network.get_config_dict()
    node_dict=network.get_node_dict()
    source_dict=network.get_source_dict()
    target_dict=network.get_target_dict()
    master_dict=network.get_master_dict()

    #update nodeParam_dict, source_dict, master_dict:
    #with edge parameters:
    for node in node_dict.keys():
        (node_type,mpr_range,dnr_range)=node_dict[node]

        if (node not in knock_down_list):
            nodeParam_dict[node]=(node_type,
                                  clib.randu(ctypes.c_double(float(mpr_range[0])),
                                             ctypes.c_double(float(mpr_range[1]))),
                                  clib.randu(ctypes.c_double(float(dnr_range[0])),
                                             ctypes.c_double(float(dnr_range[1]))))
            #print(node)
            #print(nodeParam_dict[node])
        else: 
            #nodeParam_dict[node]=(node_type, ZERO_PRODUCTION_RATE, ZERO_PRODUCTION_RATE)
            #nodeParam_dict[node]=(node_type, MPR_MIN_KD, MPR_MAX_KD)
            nodeParam_dict[node]=(node_type, 
                                  clib.randu(ctypes.c_double(float(MPR_MIN_KD)), 
                                             ctypes.c_double(float(MPR_MAX_KD))), 
                                  clib.randu(ctypes.c_double(float(dnr_range[0])), 
                                             ctypes.c_double(float(dnr_range[1]))))
            #print(node)
            #print(nodeParam_dict[node])

        #print(nodeParam_dict[node])
        for (idx, e) in target_dict[node].items():
            target,reg_type,tsh_range,hco_range,fch_range=e
            tsh=clib.randu(ctypes.c_double(tsh_range[0]),\
                           ctypes.c_double(tsh_range[1]))
            hco=int(clib.randu(ctypes.c_double(hco_range[0]),
                               ctypes.c_double(hco_range[1]+1)))

            if int(reg_type) in EXCITATION_TYPES: #any of the EXCITATORY links
                fch=clib.randu(ctypes.c_double(fch_range[0]),
                               ctypes.c_double(fch_range[1]))
            elif int(reg_type) in INHIBITION_TYPES: #any of the inhibitory links
                fch=1/clib.randu(ctypes.c_double(fch_range[0]),\
                                 ctypes.c_double(fch_range[1]))
            #store the edge parameters in source_dict:
            source_dict[target][idx]=source_dict[target][idx]+[tsh,hco,fch]
            #store the edge parameters in master_dict:
            master_dict[idx]=(node,target,reg_type,tsh,hco,fch)
    return (nodeParam_dict,source_dict,master_dict)

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
    for idx,v in d.items(): 
        print(idx,':')
        for X,e in v.items(): 
            print('  ',X,': ', np.log2(e))
    print('\n')
    return None



#----------------------------------------------------------------------#
def generate_models_random_ICs(network, model_params_dict): 
    '''
    This method generates RACIPE models. 
    Two arguments: 
        network: an object of type Network 
        model_params_dict: a dictionary whose 
            keys: models
            values: a dictionary whose 
                keys: parameter names mpr, dnr, tsh, hco, fch
                values: numpy array of values of parameters for each type
    '''
    import numpy
    import numpy as np
    import ctypes
    import sys
    global clib

    print("Generating models ... ")

    node_dict = network.get_node_dict()
     
    node_id_dict=network.get_node_id_dict()

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

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

    #create arrays for placing node and edge parameters 
    #to be passed to C function:
    MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #create array for placing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    #import config_dict:
    config_dict=network.get_config_dict()

    # place node types onto numpy array
    # facilitates sending them to C:
    i=0
    for node in node_dict.keys():
        (node_type, mpr_range, dnr_range) = node_dict[node]
        NODE_TYPE_arr[i] = node_type
        i=i+1

    # calculate steady states for each model:
    for model_no, param_dict in model_params_dict.items():
        #print(model_no)

        # place the params onto numpy array  
        # facilitates sending them to C:
        MPR_arr =  param_dict['mpr']
        DNR_arr =  param_dict['dnr']
        TSH_arr =  param_dict['tsh']
        HCO_arr =  param_dict['hco']
        FCH_arr =  param_dict['fch']

        # numpy array to fetch solutions from C to python:
        EXP_dict_arr=np.zeros(NUM_NODES*int(config_dict['NUM_RANDOM_ICS']),
                              dtype=np.double) 

        clib.find_solutions(\
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
        #end_time=time.time()
    return 0


#----------------------------------------------------------------------#
#def generate_models_random_ICs(network, model_params_dict):
def generate_models_signaled_random_ICs(network, model_params_dict):
    '''
    This method generates RACIPE models. 
    Two arguments: 
        network: an object of type Network 
        model_params_dict: a dictionary whose 
            keys: models
            values: a dictionary whose 
                keys: parameter names mpr, dnr, tsh, hco, fch
                values: numpy array of values of parameters for each type
    '''
    import numpy
    import numpy as np
    import ctypes
    import sys
    global clib

    print("Generating models ... ")

    node_dict = network.get_node_dict()
     
    node_id_dict=network.get_node_id_dict()

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

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

    #create arrays for placing node and edge parameters 
    #to be passed to C function:
    MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #create array for placing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    #import config_dict:
    config_dict=network.get_config_dict()

    # place node types onto numpy array
    # facilitates sending them to C:
    i=0
    for node in node_dict.keys():
        (node_type, mpr_range, dnr_range) = node_dict[node]
        NODE_TYPE_arr[i] = node_type
        i=i+1

    #print("number of nodes")
    #print(len(MPR_arr))

    SIGNALING_NODE_ID=8

    #print(SIGNALING_NODE_ID)

    FCHANGE = 2**29 #r15:2**28 #r14:2**27 #r13:2**26 #r12:2**25 #r11:2**24 
    #r4:2**8 #r6:2**15 #r5:2**10 #r7:2**20 #r8:2**21 
    #r9:2**22 #r10:2**23 #r3:2**7 #r2:2**6 #r1:2**5 
    #r24:2**(-9) #r23:2**(-8) #r22:2**(-7) #r21:2**(-6)

    #FCHANGE = 2**23 #r12:2**(-6) #r11:0.03125 
    #r8:2**22 #7:2**21 #r6:2**20 #r5:2**10  #r4:256 #r3:128 #r2:64  r1:32

    # calculate steady states for each model:
    for model_no, param_dict in model_params_dict.items():
        print(model_no)

        # place the params onto numpy array  
        # facilitates sending them to C:
        MPR_arr =  param_dict['mpr']
        DNR_arr =  param_dict['dnr']
        TSH_arr =  param_dict['tsh']
        HCO_arr =  param_dict['hco']
        FCH_arr =  param_dict['fch']

        MPR_arr[SIGNALING_NODE_ID] = MPR_arr[SIGNALING_NODE_ID] * FCHANGE

        # numpy array to fetch solutions from C to python:
        EXP_dict_arr=np.zeros(NUM_NODES*int(config_dict['NUM_RANDOM_ICS']),
                              dtype=np.double) 

        clib.find_solutions(\
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
        #end_time=time.time()
    return 0


#----------------------------------------------------------------------#
def generate_models_signaled_sstate(network, model_params_dict, sstate_dict):
    '''
    This method generates RACIPE models. 
    Thre arguments: 
        network: an object of type Network 
        model_params_dict: a dictionary whose 
            keys: models
            values: a dictionary whose 
                keys: parameter names mpr, dnr, tsh, hco, fch
                values: numpy array of values of parameters for each type
        sstate_dict: a dictionary whose 
            keys: monostabel models
            values: steady state expressions for the model
    '''
    import numpy
    import numpy as np
    import ctypes
    import sys
    global clib

    print("Generating models ... ")

    node_dict = network.get_node_dict()
     
    node_id_dict=network.get_node_id_dict()

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

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

    #create arrays for placing node and edge parameters 
    #to be passed to C function:
    MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #create array for placing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    #import config_dict:
    config_dict=network.get_config_dict()

    # place node types onto numpy array
    # facilitates sending them to C:
    i=0
    for node in node_dict.keys():
        (node_type, mpr_range, dnr_range) = node_dict[node]
        NODE_TYPE_arr[i] = node_type
        i=i+1

    SIGNALING_NODE_ID = 8
    FCHANGE = 2**66 
 
    # calculate steady states for each model:
    for model_no, param_dict in model_params_dict.items():
        #print(model_no)
        sstate_exp = sstate_dict[model_no]

        #print(sstate_exp) 
        #sys.exit(0)

        # place the params onto numpy array  
        # facilitates sending them to C:
        MPR_arr =  param_dict['mpr']
        DNR_arr =  param_dict['dnr']
        TSH_arr =  param_dict['tsh']
        HCO_arr =  param_dict['hco']
        FCH_arr =  param_dict['fch']

        MPR_arr[SIGNALING_NODE_ID] = MPR_arr[SIGNALING_NODE_ID] * FCHANGE

        # numpy array to fetch solutions from C to python:
        EXP_dict_arr=np.zeros(NUM_NODES*int(config_dict['NUM_RANDOM_ICS']),
                              dtype=np.double) 

        clib.find_solutions_sstate(\
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
             ctypes.c_void_p(sstate_exp.ctypes.data),
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
        #end_time=time.time()
        #break
    return 0


#----------------------------------------------------------------------#
def generate_models_sstate(network, model_params_dict, sstate_dict):
    '''
    This method generates RACIPE models. 
    Thre arguments: 
        network: an object of type Network 
        model_params_dict: a dictionary whose 
            keys: models
            values: a dictionary whose 
                keys: parameter names mpr, dnr, tsh, hco, fch
                values: numpy array of values of parameters for each type
        sstate_dict: a dictionary whose 
            keys: monostabel models
            values: steady state expressions for the model
    '''
    import numpy
    import numpy as np
    import ctypes
    import sys
    global clib

    print("Generating models ... ")

    node_dict = network.get_node_dict()
     
    node_id_dict=network.get_node_id_dict()

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

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

    #create arrays for placing node and edge parameters 
    #to be passed to C function:
    MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #create array for placing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    #import config_dict:
    config_dict=network.get_config_dict()

    # place node types onto numpy array
    # facilitates sending them to C:
    i=0
    for node in node_dict.keys():
        (node_type, mpr_range, dnr_range) = node_dict[node]
        NODE_TYPE_arr[i] = node_type
        i=i+1

    # calculate steady states for each model:
    for model_no, param_dict in model_params_dict.items():
        #print(model_no)
        sstate_exp = sstate_dict[model_no]

        #print(sstate_exp) 
        #sys.exit(0)

        # place the params onto numpy array  
        # facilitates sending them to C:
        MPR_arr =  param_dict['mpr']
        DNR_arr =  param_dict['dnr']
        TSH_arr =  param_dict['tsh']
        HCO_arr =  param_dict['hco']
        FCH_arr =  param_dict['fch']

        # numpy array to fetch solutions from C to python:
        EXP_dict_arr=np.zeros(NUM_NODES*int(config_dict['NUM_RANDOM_ICS']),
                              dtype=np.double) 

        clib.find_solutions_sstate(\
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
             ctypes.c_void_p(sstate_exp.ctypes.data),
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
        #end_time=time.time()
    return 0

#----------------------------------------------------------------------#
def generate_e2m_traj_mpr_signal(network, model_params_dict, sstate_dict):
    '''
    This method generates trajectories for RACIPE models. 
    Three arguments: 
        network: an object of type Network 
        model_params_dict: a dictionary whose 
            keys: models
            values: a dictionary whose 
                keys: parameter names mpr, dnr, tsh, hco, fch
                values: numpy array of values of parameters for each type
        sstate_dict: a dictionary whose 
            keys: monostabel models
            values: steady state expressions for the model
    '''
    import numpy
    import numpy as np
    import ctypes
    import sys
    global clib

    print("Generating trajectories ... ")

    node_dict = network.get_node_dict()
     
    node_id_dict=network.get_node_id_dict()

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

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

    #create arrays for placing node and edge parameters 
    #to be passed to C function:
    MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #create array for placing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    #import config_dict:
    config_dict=network.get_config_dict()

    # place node types onto numpy array
    # facilitates sending them to C:
    i=0
    for node in node_dict.keys():
        (node_type, mpr_range, dnr_range) = node_dict[node]
        NODE_TYPE_arr[i] = node_type
        i=i+1

    SIGNALING_NODE_ID=8 
        # TGFbeta: 9th node - counting starts from 0 in python becomes 8th node

    # MPR signal fold change:
    FCHANGE = 2**16 

    # calculate steady states for each model:
    for model_no, param_dict in model_params_dict.items():
        sstate_exp = sstate_dict[model_no]
        #sstate_exp = estate_dict[model_no]

        # place the params onto numpy array  
        # facilitates sending them to C:
        MPR_arr =  param_dict['mpr']
        DNR_arr =  param_dict['dnr']
        TSH_arr =  param_dict['tsh']
        HCO_arr =  param_dict['hco']
        FCH_arr =  param_dict['fch']

        # Increase MPR signal:
        MPR_arr[SIGNALING_NODE_ID] = MPR_arr[SIGNALING_NODE_ID] * FCHANGE

        # numpy array to fetch solutions from C to python:
        EXP_dict_arr=np.zeros(NUM_NODES*int(config_dict['NUM_RANDOM_ICS']),
                              dtype=np.double) 

        clib.find_traj_sstate_signaled(\
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
             ctypes.c_void_p(sstate_exp.ctypes.data),
             ctypes.c_void_p(MPR_arr.ctypes.data), 
             ctypes.c_void_p(DNR_arr.ctypes.data), 
             ctypes.c_void_p(NODE_TYPE_arr.ctypes.data),
             ctypes.c_void_p(edge_source_arr.ctypes.data), 
             ctypes.c_void_p(edge_target_arr.ctypes.data), 
             ctypes.c_void_p(edge_type_arr.ctypes.data), 
             ctypes.c_void_p(TSH_arr.ctypes.data), 
             ctypes.c_void_p(HCO_arr.ctypes.data), 
             ctypes.c_void_p(FCH_arr.ctypes.data)
             )
        #end_time=time.time()
        #sys.exit(0)
    return 0


#----------------------------------------------------------------------#
def generate_traj_mpr_signal(network, model_params_dict, sstate_dict, FCHANGE):
    '''
    This method generates trajectories for RACIPE models. 
    Three arguments: 
        network: an object of type Network 
        model_params_dict: a dictionary whose 
            keys: models
            values: a dictionary whose 
                keys: parameter names mpr, dnr, tsh, hco, fch
                values: numpy array of values of parameters for each type
        sstate_dict: a dictionary whose 
            keys: monostabel models
            values: steady state expressions for the model
    '''
    import numpy
    import numpy as np
    import ctypes
    import sys
    global clib

    print("Generating trajectories ... ")

    node_dict = network.get_node_dict()
     
    node_id_dict=network.get_node_id_dict()

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

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

    #create arrays for placing node and edge parameters 
    #to be passed to C function:
    MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #create array for placing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    #import config_dict:
    config_dict=network.get_config_dict()

    # place node types onto numpy array
    # facilitates sending them to C:
    i=0
    for node in node_dict.keys():
        (node_type, mpr_range, dnr_range) = node_dict[node]
        NODE_TYPE_arr[i] = node_type
        i=i+1

    SIGNALING_NODE_ID=8 
        # TGFbeta: 9th node - counting starts from 0 in python becomes 8th node

    # Calculate steady states for each model:
    for model_no, param_dict in model_params_dict.items():
        sstate_exp = sstate_dict[model_no]

        # place the params onto numpy array  
        # facilitates sending them to C:
        MPR_arr =  param_dict['mpr']
        DNR_arr =  param_dict['dnr']
        TSH_arr =  param_dict['tsh']
        HCO_arr =  param_dict['hco']
        FCH_arr =  param_dict['fch']

        # Update MPR signal:
        MPR_arr[SIGNALING_NODE_ID] = MPR_arr[SIGNALING_NODE_ID] * FCHANGE

        # numpy array to fetch solutions from C to python:
        EXP_dict_arr=np.zeros(NUM_NODES*int(config_dict['NUM_RANDOM_ICS']),
                              dtype=np.double) 

        clib.find_traj_sstate_signaled(\
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
             ctypes.c_void_p(sstate_exp.ctypes.data),
             ctypes.c_void_p(MPR_arr.ctypes.data), 
             ctypes.c_void_p(DNR_arr.ctypes.data), 
             ctypes.c_void_p(NODE_TYPE_arr.ctypes.data),
             ctypes.c_void_p(edge_source_arr.ctypes.data), 
             ctypes.c_void_p(edge_target_arr.ctypes.data), 
             ctypes.c_void_p(edge_type_arr.ctypes.data), 
             ctypes.c_void_p(TSH_arr.ctypes.data), 
             ctypes.c_void_p(HCO_arr.ctypes.data), 
             ctypes.c_void_p(FCH_arr.ctypes.data)
             )
        #end_time=time.time()
        #sys.exit(0)
    return 0


#----------------------------------------------------------------------#
def generate_traj_incr_sig_gradually(network, model_params_dict, sstate_dict, \
                                     SIGNALING_NODE_ID, FCHANGE):
    '''
    This method generates trajectories for RACIPE models. 
    Three arguments: 
        network: an object of type Network 
        model_params_dict: a dictionary whose 
            keys: models
            values: a dictionary whose 
                keys: parameter names mpr, dnr, tsh, hco, fch
                values: numpy array of values of parameters for each type
        sstate_dict: a dictionary whose 
            keys: monostable models
            values: steady state expressions for the model
    '''
    import numpy
    import numpy as np
    import ctypes
    import sys
    global clib

    print("Generating trajectories ... ")

    node_dict = network.get_node_dict()
     
    node_id_dict=network.get_node_id_dict()

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

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

    #create arrays for placing node and edge parameters 
    #to be passed to C function:
    MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #create array for placing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    #print("here")
    #sys.exit(0)
    #import config_dict:
    config_dict=network.get_config_dict()

    # place node types onto numpy array
    # facilitates sending them to C:
    i=0
    for node in node_dict.keys():
        (node_type, mpr_range, dnr_range) = node_dict[node]
        NODE_TYPE_arr[i] = node_type
        i=i+1

    #SIGNALING_NODE_ID=8 # TGFbeta: 9th node - counting starts from 0 in python becomes 8th node
    #SIGNALING_NODE_ID=2 # node C: 3rd node - counting starts from 0 in python becomes 2nd node

    # Calculate steady states for each model:
    for model_no, param_dict in model_params_dict.items():
        sstate_exp = sstate_dict[model_no]

        # place the params onto numpy array  
        # facilitates sending them to C:
        MPR_arr =  param_dict['mpr']
        DNR_arr =  param_dict['dnr']
        TSH_arr =  param_dict['tsh']
        HCO_arr =  param_dict['hco']
        FCH_arr =  param_dict['fch']

        clib.find_traj_sstate_incr_sig_gradually(\
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
             ctypes.c_void_p(sstate_exp.ctypes.data),
             ctypes.c_void_p(MPR_arr.ctypes.data), 
             ctypes.c_void_p(DNR_arr.ctypes.data), 
             ctypes.c_void_p(NODE_TYPE_arr.ctypes.data),
             ctypes.c_void_p(edge_source_arr.ctypes.data), 
             ctypes.c_void_p(edge_target_arr.ctypes.data), 
             ctypes.c_void_p(edge_type_arr.ctypes.data), 
             ctypes.c_void_p(TSH_arr.ctypes.data), 
             ctypes.c_void_p(HCO_arr.ctypes.data), 
             ctypes.c_void_p(FCH_arr.ctypes.data),
             ctypes.c_int(SIGNALING_NODE_ID),
             ctypes.c_double(FCHANGE)
             )
        #end_time=time.time()
        #sys.exit(0)
    return 0

#----------------------------------------------------------------------#
def generate_traj_decr_sig_gradually(network, model_params_dict, sstate_dict, \
                                     SIGNALING_NODE_ID, FCHANGE):
    '''
    This method generates trajectories for RACIPE models. 
    Three arguments: 
        network: an object of type Network 
        model_params_dict: a dictionary whose 
            keys: models
            values: a dictionary whose 
                keys: parameter names mpr, dnr, tsh, hco, fch
                values: numpy array of values of parameters for each type
        sstate_dict: a dictionary whose 
            keys: monostabel models
            values: steady state expressions for the model
    '''
    import numpy
    import numpy as np
    import ctypes
    import sys
    global clib

    print("Generating trajectories ... ")

    node_dict = network.get_node_dict()
     
    node_id_dict=network.get_node_id_dict()

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

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

    #create arrays for placing node and edge parameters 
    #to be passed to C function:
    MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #create array for placing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    #import config_dict:
    config_dict=network.get_config_dict()

    # place node types onto numpy array
    # facilitates sending them to C:
    i=0
    for node in node_dict.keys():
        (node_type, mpr_range, dnr_range) = node_dict[node]
        NODE_TYPE_arr[i] = node_type
        i=i+1

    #SIGNALING_NODE_ID=8 # TGFbeta: 9th node - counting starts from 0 in python becomes 8th node
    #SIGNALING_NODE_ID=2 # node C: 3rd node - counting starts from 0 in python becomes 2nd node

    # Calculate steady states for each model:
    for model_no, param_dict in model_params_dict.items():
        sstate_exp = sstate_dict[model_no]

        # place the params onto numpy array  
        # facilitates sending them to C:
        MPR_arr =  param_dict['mpr']
        DNR_arr =  param_dict['dnr']
        TSH_arr =  param_dict['tsh']
        HCO_arr =  param_dict['hco']
        FCH_arr =  param_dict['fch']

        clib.find_traj_sstate_decr_sig_gradually(\
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
             ctypes.c_void_p(sstate_exp.ctypes.data),
             ctypes.c_void_p(MPR_arr.ctypes.data), 
             ctypes.c_void_p(DNR_arr.ctypes.data), 
             ctypes.c_void_p(NODE_TYPE_arr.ctypes.data),
             ctypes.c_void_p(edge_source_arr.ctypes.data), 
             ctypes.c_void_p(edge_target_arr.ctypes.data), 
             ctypes.c_void_p(edge_type_arr.ctypes.data), 
             ctypes.c_void_p(TSH_arr.ctypes.data), 
             ctypes.c_void_p(HCO_arr.ctypes.data), 
             ctypes.c_void_p(FCH_arr.ctypes.data),
             ctypes.c_int(SIGNALING_NODE_ID),
             ctypes.c_double(FCHANGE)
             )
        #end_time=time.time()
        #sys.exit(0)
    return 0

#----------------------------------------------------------------------#
def generate_trajectory_sstate(network, model_params_dict, sstate_dict, estate_dict):
    '''
    This method generates trajectories for RACIPE models. 
    Three arguments: 
        network: an object of type Network 
        model_params_dict: a dictionary whose 
            keys: models
            values: a dictionary whose 
                keys: parameter names mpr, dnr, tsh, hco, fch
                values: numpy array of values of parameters for each type
        sstate_dict: a dictionary whose 
            keys: monostabel models
            values: steady state expressions for the model
    '''
    import numpy
    import numpy as np
    import ctypes
    import sys
    global clib

    print("Generating trajectories ... ")

    node_dict = network.get_node_dict()
     
    node_id_dict=network.get_node_id_dict()

    edge_source_dict=network.get_edge_source_dict()
    edge_target_dict=network.get_edge_target_dict()
    edge_type_dict=network.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

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

    #create arrays for placing node and edge parameters 
    #to be passed to C function:
    MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #create array for placing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    #import config_dict:
    config_dict=network.get_config_dict()

    # place node types onto numpy array
    # facilitates sending them to C:
    i=0
    for node in node_dict.keys():
        (node_type, mpr_range, dnr_range) = node_dict[node]
        NODE_TYPE_arr[i] = node_type
        i=i+1

    print(node_dict.keys())
    print(NODE_TYPE_arr)
    SIGNALING_NODE_ID=8 
    print(SIGNALING_NODE_ID)

    # calculate steady states for each model:
    for model_no, param_dict in model_params_dict.items():
        print(model_no)
        sstate_exp = sstate_dict[model_no]
        #sstate_exp = estate_dict[model_no]

        #print(sstate_exp) 
        #sys.exit(0)

        # place the params onto numpy array  
        # facilitates sending them to C:
        MPR_arr =  param_dict['mpr']
        DNR_arr =  param_dict['dnr']
        TSH_arr =  param_dict['tsh']
        HCO_arr =  param_dict['hco']
        FCH_arr =  param_dict['fch']

        # numpy array to fetch solutions from C to python:
        EXP_dict_arr=np.zeros(NUM_NODES*int(config_dict['NUM_RANDOM_ICS']),
                              dtype=np.double) 

        clib.find_trajectory_sstate(\
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
             ctypes.c_void_p(sstate_exp.ctypes.data),
             ctypes.c_void_p(MPR_arr.ctypes.data), 
             ctypes.c_void_p(DNR_arr.ctypes.data), 
             ctypes.c_void_p(NODE_TYPE_arr.ctypes.data),
             ctypes.c_void_p(edge_source_arr.ctypes.data), 
             ctypes.c_void_p(edge_target_arr.ctypes.data), 
             ctypes.c_void_p(edge_type_arr.ctypes.data), 
             ctypes.c_void_p(TSH_arr.ctypes.data), 
             ctypes.c_void_p(HCO_arr.ctypes.data), 
             ctypes.c_void_p(FCH_arr.ctypes.data), 
             ctypes.c_int(SIGNALING_NODE_ID),
             ctypes.c_void_p(EXP_dict_arr.ctypes.data)
             )
        #end_time=time.time()
        #sys.exit(0)
    return 0


#-----------------------------------------------------------------------------#
if __name__ == '__main__':
    print (sys.argv[0] + ':')
    print(__doc__)
    sys.exit(0)
