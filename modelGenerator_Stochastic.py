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


clib=ctypes.cdll.LoadLibrary('./simulation_clib.so')
clib.randu.argtypes=(ctypes.c_double,ctypes.c_double)
clib.randu.restype=ctypes.c_double

SLOW_EDGE_TYPES=[1,2,3,4]
FAST_EDGE_TYPES=[5,6]

EXCITATION_TYPES=[1,4,5] #4: degrdation inhibition => activation of degradation
INHIBITION_TYPES=[2,3,6] #3: dedradation activation => inhibition of degradation

#-----------------------------------------------------------------------------#
def set_parameters(rcp):
    global clib
    #nodeParam_dict: node to (MPR,DNR) mapping:
    #key: node, value: tuple (MPR,DNR)
    nodeParam_dict=OrderedDict()

    #import config_dict, node_dict, source_dict, target_dict:
    config_dict=rcp.get_config_dict()
    node_dict=rcp.get_node_dict()
    source_dict=rcp.get_source_dict()
    target_dict=rcp.get_target_dict()
    master_dict=rcp.get_master_dict()

    #update nodeParam_dict, source_dict, master_dict:
    #with edge parameters:
    for node in node_dict.keys():
        (node_type,mpr_range,dnr_range)=node_dict[node]
        nodeParam_dict[node]=(node_type,
                              clib.randu(ctypes.c_double(float(mpr_range[0])),
                                         ctypes.c_double(float(mpr_range[1]))),
                              clib.randu(ctypes.c_double(float(dnr_range[0])),
                                         ctypes.c_double(float(dnr_range[1]))))
        for (idx,e) in target_dict[node].items():
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
def sum_delta(expression_dict,expression_dict_prev): 
    sumDelta=0
    for X in expression_dict.keys(): 
        sumDelta+=np.square(expression_dict[X]-expression_dict_prev[X])
    return sumDelta

#-----------------------------------------------------------------------------#
def find_solutions(nodeParam_dict,config_dict,source_dict):
    count_iteration=0
    expression_dict_ICs=OrderedDict()
    while count_iteration<int(config_dict['NUM_RANDOM_ICS']): 
        expression_dict=set_ICs(config_dict,nodeParam_dict,source_dict)

        expression_dict_ICs[count_iteration]=(estimate_stable_expression(\
                                              config_dict,
                                              nodeParam_dict,source_dict, 
                                              expression_dict.copy())).copy()
        count_iteration+=1
    return expression_dict_ICs

#-----------------------------------------------------------------------------#
def cal_average(solution_count,solution_dict,expression_dict):
    updated_solution_dict=OrderedDict()
    for X in solution_dict.keys():
        updated_solution_dict[X]=(solution_count*solution_dict[X]+\
                                  expression_dict[X])/(solution_count+1)
    return updated_solution_dict

#-----------------------------------------------------------------------------#
def print_limitCycle(expression_dict_ICs): 
    IC_indices=sorted(expression_dict_ICs.keys())
    #for k in expression_dict_ICs[0].keys():
        #print(expression_dict_ICs[0][k])
    for IC_no in range(IC_indices[0],IC_indices[-1]+1): 
        print('\t'.join(str(v) for k,v in expression_dict_ICs[IC_no].items()))
    return None

def write_limitCycle_trace(expression_dict_ICs,fh_LCtrace): 
    IC_indices=sorted(expression_dict_ICs.keys())
    for IC_no in range(IC_indices[0],IC_indices[-1]+1): 
        outstr ='\t'.join(str('%10.6f'%(np.log2(v))) \
                 for k,v in expression_dict_ICs[IC_no].items())
        outstr+='\n'
        fh_LCtrace.write(outstr) 
    return None

def write_limitCycle_trace_old(expression_dict_ICs,fh_LCtrace): 
    IC_indices=sorted(expression_dict_ICs.keys())
    for IC_no in range(IC_indices[0],IC_indices[-1]+1): 
        outstr=""
        for k,v in expression_dict_ICs[IC_no].items():
            outstr+=str(v)+'\t'
        outstr+='\n'
        fh_LCtrace.write(outstr) 
    return None

#-----------------------------------------------------------------------------#
def count_states(config_dict,expression_dict_ICs,fh_LCtrace):
    solution_dict=defaultdict(lambda:OrderedDict())
    solution_count_dict=defaultdict(int)

    state_count=1
    solution_dict[state_count]=expression_dict_ICs[0]
    solution_count_dict[state_count]=1 

    IC_indices=sorted(expression_dict_ICs.keys())
    for IC_no in range(IC_indices[0],IC_indices[-1]+1): 
        found_match=False
        for sol_no in range(1,len(solution_dict.keys())+1):
            testDelta=sum_delta(solution_dict[sol_no],
                                expression_dict_ICs[IC_no])
            if testDelta<=float(config_dict['CONVERGENCE_PROXIMITY']):
                solution_dict[sol_no]= \
                    cal_average(solution_count_dict[sol_no],
                                                  solution_dict[sol_no],
                                                  expression_dict_ICs[IC_no])
                solution_count_dict[sol_no]+=1
                found_match=True
                break
        if(not found_match):
            state_count+=1
            solution_dict[state_count]=expression_dict_ICs[IC_no]
            solution_count_dict[state_count]=1
        if state_count>=int(config_dict['MAX_STABLE_STATES']):
            write_limitCycle_trace(expression_dict_ICs,fh_LCtrace)
            fh_LCtrace.write("\n") 
            break
    return (solution_dict,solution_count_dict)

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
    for idx,v in d.items(): 
        print(idx,':')
        for X,e in v.items(): 
            print('  ',X,': ', np.log2(e))
    print('\n')
    return None

#-----------------------------------------------------------------------------#
def save_parameters(fh_params,model_no,nodeParam_dict,
                    master_dict,solution_dict):

    #print(type(fh_params))
    #model no + number of states:
    outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    #skip node type information nodeParam_dict[X][0]
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
    outstr='\n'
    fh_params.write(outstr)

    fh_params.flush()
    return None

def save_parameters_5(fh_dict_nodeparams, fh_dict_edgeparams,
                      model_no,nodeParam_dict,master_dict,solution_dict):

    #for (k,fh_params) in fh_dict_nodeparams.items(): 
    #    print(k) 
    #    print(type(fh_params))

    (fh_mpr,fh_dnr)=fh_dict_nodeparams.values() 
    #print(fh_mpr) 
    #print(fh_dnr) 
    #print(type(fh_mpr))
    #model no + number of states:
    #outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    outstr=str(model_no)
    #skip node type information nodeParam_dict[X][0]
    #save MPRs:
    for X in nodeParam_dict.keys():
        #outstr+='\t'+str('%10.6f'%nodeParam_dict[X][0])
        outstr+='\t'+str('%10.6f'%nodeParam_dict[X][1])
    outstr+='\n'
    #fh_params.write(outstr)
    fh_mpr.write(outstr)

    #save all DNRs:
    #outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    outstr=str(model_no)
    for X in nodeParam_dict.keys():
        #outstr+='\t'+str('%10.6f'%nodeParam_dict[X][1])
        outstr+='\t'+str('%10.6f'%nodeParam_dict[X][2])
    outstr+='\n'
    #fh_params.write(outstr)
    fh_dnr.write(outstr)

    (fh_tsh,fh_hco,fh_fch)=fh_dict_edgeparams.values() 
    #save all TSHs:
    #outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    outstr=str(model_no)
    for idx in master_dict.keys():
        outstr+='\t'+str('%10.6f'%master_dict[idx][3])
    outstr+='\n'
    #fh_params.write(outstr)
    fh_tsh.write(outstr)
    #save all HCOs:
    #outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    outstr=str(model_no)
    for idx in master_dict.keys():
        outstr+='\t'+str(master_dict[idx][4])
    outstr+='\n'
    #fh_params.write(outstr)
    fh_hco.write(outstr)

    #save all FCHs:
    #outstr=str(model_no)+'\t'+str(len(solution_dict.keys()))
    outstr=str(model_no)
    for idx in master_dict.keys():
        outstr+='\t'+str('%10.6f'%master_dict[idx][5])
    outstr+='\n'
    #fh_params.write(outstr)
    fh_fch.write(outstr)

    #fh_params.flush()
    fh_mpr.flush()
    fh_dnr.flush()
    fh_tsh.flush()
    fh_hco.flush()
    fh_fch.flush()
    return None

#-----------------------------------------------------------------------------#
def open_parameter_files(rcp): 
    ''' 
    This method opens output files for writing solutions. 
    It creats MAX_STABLE_STATES number of file handles in the write mode, 
    places in a dictionary, and returns the dictionary. 
    '''

    fh_dict_nodeparams=OrderedDict()
    fh_dict_edgeparams=OrderedDict()

    fname_dict_nodeparams=rcp.get_fname_nodeparams()
    fname_dict_edgeparams=rcp.get_fname_edgeparams()

    for (k,fname) in fname_dict_nodeparams.items(): 
        fh_dict_nodeparams[k]=open(fname,'a')
        #print(type(fh_dict_nodeparams[k]))

    for (k,fname) in fname_dict_edgeparams.items(): 
        fh_dict_edgeparams[k]=open(fname,'a')
        #print(type(fh_dict_edgeparams[k]))
    return fh_dict_nodeparams,fh_dict_edgeparams

#----------------------------------------------------------------------#
def close_parameter_files(fh_dict_nodeparams,fh_dict_edgeparams):
    for idx in fh_dict_nodeparams.keys():
        fh_dict_nodeparams[idx].close()

    for idx in fh_dict_edgeparams.keys():
        fh_dict_edgeparams[idx].close()
    return None 

#-----------------------------------------------------------------------------#
def open_solution_files(rcp): 
    ''' 
    This method opens output files for writing solutions. 
    It creats MAX_STABLE_STATES number of file handles in the write mode, 
    places in a dictionary, and returns the dictionary. 
    '''
    save_path=rcp.get_work_dir()
    #print('Saving files in the folder: '+save_path)
    config_dict=rcp.get_config_dict()
    tpo_fname=rcp.get_tpo_fname()
    fn_prefix=basename(tpo_fname.strip()).\
              split(config_dict['TOPOLOGY_FNAME_EXTENSION'])[0]

    #create file names and save them in a dictionary:
    fname_dict_solutions=OrderedDict()
    for num_states in range(1,int(config_dict['MAX_STABLE_STATES'])+1):
        out_fname=save_path+'/'+\
                  fn_prefix+\
                  config_dict['SOLUTION_FNAME_INFIX']+\
                  str(num_states)+config_dict['SOLUTION_FNAME_EXTENSION']
        fname_dict_solutions[num_states]=out_fname
    #open files and save the file handles in a dictionary: 
    fh_dict_solutions=OrderedDict()
    for idx in fname_dict_solutions.keys(): 
        fh_dict_solutions[idx]=open(fname_dict_solutions[idx],'w')
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
        for X in solution_dict[state_no].keys(): 
            outstr+="\t"+str('%10.6f'%(np.log2(solution_dict[state_no][X])))
    outstr+="\n"
    fh_dict_solutions[num_states].write(outstr)
    return None

#----------------------------------------------------------------------#
def flush_solutions(fh_dict_solutions):
    for idx in fh_dict_solutions.keys():
        fh_dict_solutions[idx].flush()
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
def generate_models(rcp): 
    '''
    This method generates RACIPE models. 
    '''
    import numpy
    import numpy as np
    import ctypes
    import sys
    global clib

    node_id_dict=rcp.get_node_id_dict()

    edge_source_dict=rcp.get_edge_source_dict()
    edge_target_dict=rcp.get_edge_target_dict()
    edge_type_dict=rcp.get_edge_type_dict()

    NUM_NODES=len(node_id_dict.keys()) #size of node_id_arr
    NUM_EDGES=len(edge_source_dict.keys()) #size of edge_id_arr

    # update these to add noise:
    NOISE=1.0
    NOISE_SHOT=0

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
    WORK_DIR=(rcp.get_work_dir()).encode('utf-8')

    #create arrays for storing node and edge parameters:
    MPR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NOISE_strength=np.zeros(NUM_NODES,dtype=np.double)
    NOISE_strength_shot=np.zeros(NUM_NODES,dtype=np.double)
    DNR_arr=np.zeros(NUM_NODES,dtype=np.double)
    NODE_TYPE_arr=np.zeros(NUM_NODES,dtype=np.intc)

    TSH_arr=np.zeros(NUM_EDGES,dtype=np.double)
    HCO_arr=np.zeros(NUM_EDGES,dtype=np.intc)
    FCH_arr=np.zeros(NUM_EDGES,dtype=np.double)

    #create array for storing the types of the edges:
    edge_type_arr=np.zeros(NUM_EDGES,dtype=np.intc) 
    i=0
    for idx in edge_type_dict.keys():
        edge_type_arr[i]=edge_type_dict[idx]
        i+=1

    #open file to write parameters:
    params_fname=rcp.get_params_fname()
    fh_params=open(params_fname,'w')

    #open files to write parameters:
    fh_dict_nodeparams=OrderedDict()
    fh_dict_edgeparams=OrderedDict()
    fh_dict_nodeparams,fh_dict_edgeparams=open_parameter_files(rcp)

    #open files to write solutions:
    fname_dict_solutions,fh_dict_solutions=open_solution_files(rcp)
    #open file to write limit cycle trace:
    LCtrace_fname=rcp.get_LCtrace_fname()
    fh_LCtrace=open(LCtrace_fname,'w')
   
    #import config_dict:
    config_dict=rcp.get_config_dict()

    model_no=1
    while model_no<=int(config_dict['NUM_MODELS']):
        (nodeParam_dict,source_dict,master_dict)=set_parameters(rcp)

        i=0
        for X in nodeParam_dict.keys():  
            NODE_TYPE_arr[i]=nodeParam_dict[X][0]
            MPR_arr[i]=nodeParam_dict[X][1]
            NOISE_strength[i]=1
            NOISE_strength_shot[i]=1
            DNR_arr[i]=nodeParam_dict[X][2]
            i+=1

        for idx,e in master_dict.items():  
            TSH_arr[idx]=master_dict[idx][3]
            HCO_arr[idx]=master_dict[idx][4]
            FCH_arr[idx]=master_dict[idx][5]
        EXP_dict_arr=np.zeros(NUM_NODES*int(config_dict['NUM_RANDOM_ICS']),
                              dtype=np.double) 

        #start_time=time.time()
        clib.find_solutions_stochastic(ctypes.create_string_buffer(WORK_DIR),
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
        #end_time=time.time()

        #place the solutions in a dictionary:
        expression_dict=OrderedDict()
        expression_dict_ICs=OrderedDict()
        count_iteration=0
        for i in range(0,NUM_NODES*int(config_dict['NUM_RANDOM_ICS']),
                         NUM_NODES):
            for node,node_id in node_id_dict.items():
                expression_dict[node]=EXP_dict_arr[i+node_id]
            expression_dict_ICs[count_iteration]=expression_dict.copy()
            count_iteration+=1
        (solution_dict,solution_count_dict)=count_states(config_dict,
                                                         expression_dict_ICs,
                                                         fh_LCtrace)
        save_parameters(fh_params,model_no,nodeParam_dict,
                        master_dict,solution_dict)
        save_parameters_5(fh_dict_nodeparams, fh_dict_edgeparams,model_no,
                        nodeParam_dict,master_dict,solution_dict)

        save_solutions(fh_dict_solutions,model_no,solution_dict)

        if (not model_no%100):
            flush_solutions(fh_dict_solutions)
        model_no+=1
    fh_params.close()
    fh_LCtrace.close()
    close_solution_files(fname_dict_solutions,fh_dict_solutions)
    close_parameter_files(fh_dict_nodeparams,fh_dict_edgeparams)
    return None 

#-----------------------------------------------------------------------------#
if __name__ == '__main__':
    print (sys.argv[0] + ':')
    print(__doc__)
    sys.exit(0)
