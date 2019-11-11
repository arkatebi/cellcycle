#!/usr/bin/env python
'''
    This module has the following methods:
'''
import sys
import os
import numpy as np
from collections import OrderedDict
from collections import defaultdict
#import matplotlib.pyplot as plt
from os.path import basename
import copy
import time
import ctypes
 
TOPOLOGY_FNAME_EXTENSION='.tpo'
SOLUTION_FNAME_INFIX='_solution_'
SOLUTION_FNAME_EXTENSION='.dat'

#NUM_MODELS=1000
NUM_MODELS=100

MAX_STABLE_STATES=10

TS_path='./workspace/TS'
TS1SA_path='./workspace/TS1SA'
TS2SA_path='./workspace/TS2SA'

#-----------------------------------------------------------------------------#
def open_solution_files(rcp): 
    ''' 
    This method opens output files for writing solutions. 
    It creats MAX_STABLE_STATES number of file handles in the write mode, 
    places in a dictionary, and returns the dictionary. 
    '''
    tpo_fname=rcp.get_tpo_fname()
    fn_prefix=basename(tpo_fname.strip()).split(TOPOLOGY_FNAME_EXTENSION)[0]

    #create file names and save them in a dictionary:
    fname_dict_solutions=OrderedDict()
    for num_states in range(1,MAX_STABLE_STATES+1):
        out_fname=fn_prefix+SOLUTION_FNAME_INFIX+\
                  str(num_states)+SOLUTION_FNAME_EXTENSION
        fname_dict_solutions[num_states]=out_fname
    #for idx in fname_dict_solutions.keys(): 
    #    print(fname_dict_solutions[idx])
    
    #open files and save the file handles in a dictionary: 
    fh_dict_solutions=OrderedDict()
    for idx in fname_dict_solutions.keys(): 
        #fh=open(fname_dict_solutions[idx],'w')
        fh_dict_solutions[idx]=open(fname_dict_solutions[idx],'w')
    #for idx in fname_dict_solutions.keys(): 
    #    print(str(idx),' ',fname_dict_solutions[idx])
    return fname_dict_solutions,fh_dict_solutions 

#----------------------------------------------------------------------#
def save_solutions(fh_dict_solutions,model_no,solution_dict): 
    ''' 
    This method writes the solutions to the files. 
    It will generate one output file for each solution. 
    '''
    num_states=len(solution_dict.keys())
    #print(num_states)
    outstr=str(model_no)+"\t"+str(num_states)
    for state_no in range(1,len(solution_dict.keys())+1):
#        current_solution_dict=solution_dict[state_no]
#        for X in current_solution_dict.keys(): 
#            outstr+="\t"+str(current_solution_dict[X])
        for X in solution_dict[state_no].keys(): 
            #outstr+="\t"+str(solution_dict[state_no][X])
            #outstr+="\t"+str(np.log2(solution_dict[state_no][X]))
            outstr+="\t"+str('%10.6f'%(np.log2(solution_dict[state_no][X])))
    outstr+="\n"
    fh_dict_solutions[num_states].write(outstr)
    return None

#----------------------------------------------------------------------#
def close_solution_files(fname_dict_solutions,fh_dict_solutions): 
    for idx in fname_dict_solutions.keys(): 
        fh_dict_solutions[idx].close()
    for idx in fname_dict_solutions.keys(): 
        fsize=os.path.getsize(fname_dict_solutions[idx])
        if (not fsize):
            os.remove(fname_dict_solutions[idx])
    return None 

#----------------------------------------------------------------------#
def analyze_data(): 
    #fname_dict_solutions,fh_dict_solutions=open_solution_files(rcp)
    #close_solution_files(fname_dict_solutions,fh_dict_solutions)

    fh_TS1_solution_1=open(TS1_path+'\'+'TS_solution_1.dat','r') 
    expression_dict=OrderedDict
    while line in fh_TS1_solution_1: 
           
    return None 
#-----------------------------------------------------------------------------#
if __name__ == '__main__':
    print (sys.argv[0] + ':')
    print(__doc__)
    sys.exit(0)
