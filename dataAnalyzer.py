#!/usr/bin/env python

# dataAnalyzer.py

'''
    This module has methods for the following:
     (1) read parameter file 
     (2) calculate probabilities each edge being activated 
         (half-functional rule)
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

#------------------------------------------------------------------------------#
def print_nested_defaultdict_log2(d): 
    #print(d)
    #sys.exit(0)
    for idx,v in d.items(): 
        print(idx,':')
        for X,e in v.items(): 
            print('  ',X,': ', np.log2(e))
    print('\n')
    return None

#------------------------------------------------------------------------------#
#                            read_line_parameters
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
#------------------------------------------------------------------------------#

def read_line_parameters(fh_params, line1, 
                         NUM_NODES, NUM_EDGES):
    MPR_arr = np.zeros(NUM_NODES, dtype=np.double)
    DNR_arr = np.zeros(NUM_NODES, dtype=np.double)
    NODE_TYPE_arr = np.zeros(NUM_NODES, dtype=np.intc)

    TSH_arr = np.zeros(NUM_EDGES, dtype=np.double)
    HCO_arr = np.zeros(NUM_EDGES, dtype=np.intc)
    FCH_arr = np.zeros(NUM_EDGES, dtype=np.double)

    #read line 1:MPR
    fields = line1.strip().split()
    for i in range(2, len(fields)):
        MPR_arr[i-2] = fields[i]

    if(NUM_NODES is not len(fields[2:])): 
        print('mismatch in node count and MPR count')
        print('program exiting ...')
        sys.exit(0)

    #read line 2:DNR
    line2 = fh_params.readline()
    fields = line2.strip().split()
    
    for i in range(2, len(fields)):
        DNR_arr[i-2] = fields[i]

    #read line 3:TSH
    line3 = fh_params.readline()
    fields = line3.strip().split()
    for i in range(2, len(fields)):
        TSH_arr[i-2]=fields[i]

    if(NUM_EDGES is not len(fields[2:])):
        print('mismatch in edge count and TSH count')
        print('program exiting ...')
        sys.exit(0)

    #read line 4:HCO 
    line4 = fh_params.readline()
    fields = line4.strip().split()
    for i in range(2, len(fields)):
        HCO_arr[i-2] = fields[i]

    #read line 5:FCH 
    line5 = fh_params.readline()
    fields = line5.strip().split()
    for i in range(2, len(fields)):
        FCH_arr[i-2]=fields[i]

    return MPR_arr,DNR_arr,TSH_arr,HCO_arr,FCH_arr

#------------------------------------------------------------------------------#
#           calProb(node_dict, master_dict, exp_fname, 
#                   params_input_fname, edge_stat_fname)
    '''
    node_dict: 
    master_dict: 
    exp_fname: 
    params_input_fname: 
    edge_stat_fname:

    This function calculates the probability of each edge to 
    be activated and write it to the output file. 

    Step 1: construct TSH_dict 
        Fetch the TSH of each edge from the parameters file. 
        Do this for all models. 

    Step 2: construct edge_count_arr
        For each model: 
            (1) construct expr_arr: 
                extract expression of all nodes for the current 
                model.
            (2) construct TSH_arr: 
                extract TSH of all edges from TSH_dict for the 
                current model.
            (3) construct edge_count_arr:
                for each edge, 
                    calculate whether the expression 
                    level of the source node exceeds the TSH 
                    level of the edge from the source to the 
                    target. 
                    If it does, add this edge to the edge count. 
        edge_count_arr contains the number of models where the 
            source node expression exceeds the TSH of the edge
            from the source to the target. 

    Step 3: calculate edge probability
        (1) probability of an edge is measured as the ratio between 
            number of states where the source node expression 
            exceeds TSH of the edge and the number of models in 
            that cluster. 
            Equivalently, divide each count value in edge_count_arr 
            by the total count of models. This gives the probability 
            of that edge being activated. 

        (2) write the calculated probability to the output file. 
    '''
#------------------------------------------------------------------------------#

def calProb(node_dict, master_dict, exp_fname, 
            params_input_fname, edge_stat_fname): 
    import numpy
    import numpy as np
    import ctypes
    import sys

    NUM_NODES = len(node_dict.keys()) 
    NUM_EDGES = len(master_dict.keys()) 
    fh_input_params = open(params_input_fname,'r')

    node_id_dict=OrderedDict() 
    idx=0
    for k in node_dict.keys(): 
        node_id_dict[k] = idx
        idx = idx+1

    TSH_dict=OrderedDict() 
    idx=1
    while True:
        line = fh_input_params.readline()
        if line =="\n": # check whether a new line
           continue
        if not line: # check whether end of file is reached
            break 
        MPR_arr, DNR_arr, TSH_arr, HCO_arr, FCH_arr = \
            read_line_parameters(fh_input_params, line, 
                                 NUM_NODES, NUM_EDGES) 
        TSH_dict[idx] = np.log2(TSH_arr)
        idx = idx+1
    fh_input_params.close()

    exp_arr = np.zeros(NUM_NODES, dtype=np.double)
    edge_count_arr = np.zeros(NUM_EDGES, dtype=np.intc)
    model_count = 0

    #sys.exit(0)
    fh_exp = open(exp_fname,'r')
    next(fh_exp) # skip header line
    while True:
        line = fh_exp.readline()
        if line == "\n": # check whether a new line
           continue
        if not line: # check whether end of file is reached
            break
        fields = line.strip().split() 
        # line format: 
        # field[0]: model number 
        # field[1]: number of states
        # filed[2]: state number 
        # field[3] and forward: node expression level
        model_no = int(fields[0]) 
        for i in range(3,len(fields)):
            exp_arr[i-3] = fields[i]

        TSH_arr = TSH_dict[model_no]
        for k in master_dict.keys():
            source = master_dict[k][0]
            node_id = node_id_dict[source]

            if(exp_arr[node_id]>TSH_arr[k]):
                edge_count_arr[k] = edge_count_arr[k]+1
        model_count = model_count+1
    fh_exp.close()

    # calculate the probability of each edge: 
    fh_edge_stat = open(edge_stat_fname,'w')
    for k in master_dict.keys():
        fh_edge_stat.write(str(edge_count_arr[k]/model_count))
        fh_edge_stat.write("\n")
    fh_edge_stat.close()
    print(edge_count_arr/model_count)
    return None 


#-----------------------------------------------------------------------------#
if __name__ == '__main__':
    print (sys.argv[0] + ':')
    print(__doc__)
    sys.exit(0)
