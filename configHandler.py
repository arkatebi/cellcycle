#!/usr/bin/env python
'''
    This module has the following two methods:
    read_config: 
        This method tries to locate the configuration file in 
        the current directory or any of the subdirectories. If 
        it does not find the configuration file, it creates one
        by invoking create_config. At the end, it returns the
        content of the configuration file as an ordered dictionary.
    create_config: 
        This method creates a configureation file in the current 
        directory. 
'''
import os
import sys
import configparser as cp
from collections import OrderedDict

configDict=OrderedDict([('MAX_RUN_TIME',23.50),
                        ('SEED',139),
                        ('USER_SEED',0),
                        ('DIST_NAME','uniform'),
                        ('NUM_MODELS',10),
                        ('NUM_RANDOM_ICS',1000),
                        ('CONVERGENCE_PROXIMITY',1.0e-12),
                        ('SAMEPOINT_PROXIMITY',0.1),
                        ('ITER_FOR_ODE',200),
                        ('ITER_FOR_RELAXATION',10),
                        ('ITER_FOR_LIMIT_CYCLE',20),
                        ('EULER_SIM_TIME',10),
                        ('EULER_SIM_TIME_LAST_ANN_STEP',2000),
                        ('STARTING_NOISE_ANN',20),
                        ('STARTING_NOISE_CONSTANT',20),
                        ('NOISE_SCALING_FACTOR_ANN',0.90),
                        ('LIMIT_CYCLE_SIM_TIME',10),
                        ('EULER_SIM_STEP_SIZE',0.01),
                        ('LIMIT_CYCLE_SIM_STEP_SIZE',0.01),
                        ('MAX_STABLE_STATES',20),
                        ('MAX_LIMIT_CYCLES',10),
                        ('MAX_ALLOWED_PERIODS',100),
                        ('NO_OF_SAMPLED_PERIODS',3),
                        ('ALLOWED_ERROR_IN_PERIODS',3),
                        ('MPR_MIN',1.00),
                        ('MPR_MAX',100.00),
                        ('DNR_MIN',0.10),
                        ('DNR_MAX',1.00),
                        ('HCO_MIN',1),
                        ('HCO_MAX',6),
                        ('FCH_MIN',1.00),
                        ('FCH_MAX',100.00),
                        ('FCH_MAX_DEG',10.00),
                        ('TRANS_RATE_FACTOR',10.00),
                        ('TSH_SCALE_FACTOR_MIN',0.02),
                        ('TSH_SCALE_FACTOR_MAX',1.98),
                        ('NUM_SIM_THRESHOLD',100000),
                        ('REG_EXCITATORY',1),
                        ('REG_INHIBITORY',2),
                        ('NODE_ATTRIBUTES','MPR,DNR'),
                        ('EDGE_ATTRIBUTES','TSH,HCO,FCH'),
                        ('PRECISION_TOTAL',10),
                        ('PRECISION_AFTER_DECIMAL',6),
                        ('TOPOLOGY_FNAME_EXTENSION','.tpo'),
                        ('CONFIG_FNAME_EXTENSION','.cfg'),
                        ('PARAMETER_RANGE_FNAME_EXTENSION','.prs'),
                        ('PARAMETER_FNAME_EXTENSION','.params'),
                        ('SOLUTION_FNAME_INFIX','_solution_'),
                        ('SOLUTION_FNAME_EXTENSION','.dat'),
                        ('TOPOLOGY_FILE_HEADER','Source,Target,Type'),
                        ('PARAMETER_FILE_HEADER','Index,Parameter_name,'\
                               'Minimum_value,Maximum_value,Regulation_type'),
                        ('SUMMARY_FILE_HEADER','MODEL_NO,NO_STATES,'\
                               'NO_LIMITCYCLES,CAL_TIME_FOR_LCS')
                        ])

#----------------------------------------------------------------------#
def read_config(fname_config):
    """
    This method reads the conig file supplied by fname_config and returns
    the configuration as an ordered dictionary.
    If the config file is not found in the current direcotry or
    any subdirectory, it creates one by invoking create_config method.
    """
    #if config file does not exist, create one: 
    if(not os.path.isfile(fname_config)):
        #create_config(open(fname_config,'w'))
        create_config(fname_config)

    # Read config file:
    Config_handle = cp.ConfigParser()
    Config_handle.read(fname_config)
    configParam = OrderedDict()
    #ConfigParam['workdir'] = Config_handle.get('WORKDIR', 'DEFAULT_PATH')
    for k in configDict.keys():
        configParam[k] = Config_handle.get('DEFAULTS', k)
    return configParam

#----------------------------------------------------------------------#
def create_config_old(fh_config):
    block_name='[DEFAULTS]'+'\n'
    fh_config.write(block_name)
    for k,v in configDict.items():
        #print(k, ':', v)
        outstr=k+':'+str(v)+'\n'
        fh_config.write(outstr)
    fh_config.close()
    return None

def create_config(fname_config):
    #fh_config=open(fname_config,"w",encoding='UTF-8')
    fh_config=open(fname_config,"w")
    block_name='[DEFAULTS]'+'\n'
    fh_config.write(block_name)
    for k,v in configDict.items():
        #print(k, ':', v)
        outstr=k+':'+str(v)+'\n'
        fh_config.write(outstr)
    fh_config.close()
    return None


#**********************************************************************#
if __name__ == '__main__':
    print (sys.argv[0] + ':')
    print (__doc__)
    sys.exit(0)
