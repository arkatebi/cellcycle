#/usr/bin/env python
"""
# -*- coding: utf-8 -*-
"""
from collections import defaultdict 
from collections import OrderedDict
import sys
import os
from os.path import basename 
import time

import argParser as ap
import configHandler as ch 
import thresholdEstimator as te
import dataAnalyzer as da 

import ctypes;

SLOW_NODE=1
FAST_NODE=2
NODE_TYPES=[SLOW_NODE,FAST_NODE]

SLOW_EDGE_TYPES=[1,2,3,4]
FAST_EDGE_TYPES=[5,6]

LIMIT_CYLCE_FNAME_SUFFIX='_LCtrace.dat'
CONFIG_FNAME='racipe.cfg'

HEADER_1_PARAMS=["MODEL_NO"]
HEADER_1_STATES =["MODEL_NO","NO_STATES","STATE_NO"]
HEADER_1_LIMITCYCLES =["MODEL_NO","LIMITCYCLE_NO","PERIOD"]


HEADER_ANNEALED_STATES =["MODEL_NO","NOISE_LEVEL","STATE_NO"]

#**********************************************************************#
class Network: 
    def __init__(self,WORK_DIR):
        #ch.create_config(open(CONFIG_FNAME,'w')) 
        #Collect config file entries:
        #self.config_dict = ch.read_config(CONFIG_FNAME) 
        #print(self.config_dict)

        #self.work_dir = self.ConfigParam['workdir']
        #self.work_dir='workspace'
        self.work_dir=WORK_DIR

        # Look for workspace, and if none exists create one:
        #if not os.path.exists(self.work_dir):
        #    os.makedirs(self.work_dir) # Create work space 

        # Collect user arguments into a dictionary:
        self.parsed_dict=ap.parse_args()

        #self.user_seed=int(self.parsed_dict['user_seed']) 

        cfg_fname=self.parsed_dict['cfg_fname']
        if (not cfg_fname):
            cfg_fname=CONFIG_FNAME
        cfg_path=os.path.dirname(cfg_fname)
        if (not cfg_path or cfg_path=="."):
            cfg_path=self.work_dir
        self.cfg_fname=cfg_path+'/'+basename(cfg_fname)
        #self.cfg_fname=self.work_dir+'/'+cfg_fname
        #self.config_dict = ch.read_config(CONFIG_FNAME)
        self.config_dict=ch.read_config(self.cfg_fname)
        #self.cfg_fname=self.work_dir+'/'+basename(cfg_fname)

        #Update the internal SEED, if user supplies this:
        if ((int(self.parsed_dict['SEED'])) != (int(self.config_dict['SEED']))):
            self.config_dict['SEED']=self.parsed_dict['SEED']

        #update the config_dict with USER_SEED from commandline:
        if ((int(self.parsed_dict['user_seed'])) != (int(self.config_dict['USER_SEED']))):
            self.config_dict['USER_SEED']=self.parsed_dict['user_seed']
        #Get .tpo filename:
        tpo_fname=self.parsed_dict['tpo_fname']
        #tpo_path=os.path.dirname(self.parsed_dict['tpo_fname'])
        #
        #fn_prefix=(basename(tpo_fname)).strip().\
        #          split(self.config_dict['TOPOLOGY_FNAME_EXTENSION'])[0]
        fn_prefix=tpo_fname.strip().\
                  split(self.config_dict['TOPOLOGY_FNAME_EXTENSION'])[0]

        # Get .tpo filename with absolute path:
        #if (not tpo_path or tpo_path=="."): 
        #    tpo_path=self.work_dir
        #self.tpo_fname=tpo_path+'/'+basename(tpo_fname)
        self.tpo_fname=self.work_dir+'/'+tpo_fname

        # Get .prs filename with path:
        prs_fname=self.parsed_dict['prs_fname']
        # if .prs file is not supplied, create file name for it:
        if not prs_fname:
            self.prs_fname=self.work_dir+'/'+\
                           fn_prefix+\
                           self.config_dict['PARAMETER_RANGE_FNAME_EXTENSION']
        else: #otherwise, append working directory with it:
            self.prs_fname=self.work_dir+'/'+prs_fname
                           

        # Get .exp filename with path:
        self.exp_fname=self.work_dir+'/'+self.parsed_dict['exp_fname']

        # Get .params filename with path:
        self.params_input_fname=self.work_dir+'/'+self.parsed_dict['params_input_fname']
        # threshold test results
        self.edge_stat_fname=self.work_dir+'/'+ self.parsed_dict['outfile']+'.edge.stat.txt'

        # Create other output filenames based on .tpo file name prefix:
        self.cfg_fname=self.work_dir+fn_prefix+\
                       self.config_dict['CONFIG_FNAME_EXTENSION']
        self.params_fname=self.work_dir+'/'+\
                          fn_prefix+\
                          self.config_dict['PARAMETER_FNAME_EXTENSION']
        self.LCtrace_fname=self.work_dir+'/'+\
                           fn_prefix+\
                           LIMIT_CYLCE_FNAME_SUFFIX
        self.prs_fname_bak=self.work_dir+'/'+\
                           fn_prefix+\
                           self.config_dict[\
                           'PARAMETER_RANGE_FNAME_EXTENSION']+".bak"

        #output files for C: 
        #to write stable states:
        self.fname_states=self.work_dir+'/'+fn_prefix+".states.txt"
        self.fname_limitcycles=self.work_dir+'/'+fn_prefix+".limitcycles.txt"
        self.fname_summary=self.work_dir+'/'+fn_prefix+".summary.txt"

        #print(self.config_dict['STARTING_NOISE_CONSTANT'])
        #print(type(self.config_dict['STARTING_NOISE_CONSTANT']))
        #print("here")

        fn_suffix = ".states.cnoise." + self.config_dict['STARTING_NOISE_CONSTANT'] + ".txt"
        self.fname_states_cnoise=self.work_dir+'/'+fn_prefix+fn_suffix


        fn_suffix = ".states.annealed.txt"
        self.fname_states_annealed=self.work_dir+'/'+fn_prefix+fn_suffix

        #print(self.fname_states_cnoise)

        #create dictionary for file names to be used in C simulation:
        self.fname_dict_simu = OrderedDict()
        self.fname_dict_simu['FNAME_STATES'] = self.fname_states
        self.fname_dict_simu['FNAME_LIMITCYCLES'] = self.fname_limitcycles
        self.fname_dict_simu['FNAME_SUMMARY'] = self.fname_summary

        self.fname_dict_simu['FNAME_STATES_CNOISE'] = self.fname_states_cnoise
        self.fname_dict_simu['FNAME_STATES_ANNEALED'] = self.fname_states_annealed

        #file names for parameter ranges:MPR, DNR, TSH, HCO, and FCH:
        self.fname_dict_nodeprs=OrderedDict()
        self.fname_dict_nodeprs['mpr']=self.work_dir+'/'+fn_prefix+".mpr.prs"
        self.fname_dict_nodeprs['dnr']=self.work_dir+'/'+fn_prefix+".dnr.prs"

        self.fname_dict_edgeprs=OrderedDict()
        self.fname_dict_edgeprs['tsh']=self.work_dir+'/'+fn_prefix+".tsh.prs"
        self.fname_dict_edgeprs['hco']=self.work_dir+'/'+fn_prefix+".hco.prs"
        self.fname_dict_edgeprs['fch']=self.work_dir+'/'+fn_prefix+".fch.prs"


        #output file names for parameters:MPR, DNR, TSH, HCO, and FCH:
        self.fname_dict_nodeparams=OrderedDict()
        self.fname_dict_nodeparams['mpr']=self.work_dir+'/'+fn_prefix+".mpr.txt"
        self.fname_dict_nodeparams['dnr']=self.work_dir+'/'+fn_prefix+".dnr.txt"

        self.fname_dict_edgeparams=OrderedDict()
        self.fname_dict_edgeparams['tsh']=self.work_dir+'/'+fn_prefix+".tsh.txt"
        self.fname_dict_edgeparams['hco']=self.work_dir+'/'+fn_prefix+".hco.txt"
        self.fname_dict_edgeparams['fch']=self.work_dir+'/'+fn_prefix+".fch.txt"

        #node_dict - gene/node names:
        self.node_dict=OrderedDict()

        #edge_dict - gene/node names:
        self.edge_dict=OrderedDict()

        #source_dict - maps target to source(s):
        self.source_dict=defaultdict(lambda: defaultdict(list))
        #self.source_dict=defaultdict(lambda: OrderedDict())

        #target_dict - maps source to target(s):
        self.target_dict=defaultdict(lambda: defaultdict(list))
        #self.target_dict=defaultdict(lambda: OrderedDict())

        #master_dict - maps each equation to a unique index:
        #Key:index
        #Value:(source,target,regulation type,TSH,HLC,FCH)
        #self.master_dict=defaultdict(list)
        #self.master_dict=defaultdict(tuple)
        self.master_dict=OrderedDict()
        
        #special data structures for C-interface:
        #node_id_dict - maps each gene/node to a unique id: 
        self.node_id_dict=OrderedDict()

        #edge_id_dict - maps each equation to a unique id: 
        self.edge_id_dict=OrderedDict()
      
        #edge_source_dict - maps each edge index to its source id: 
        self.edge_source_dict=OrderedDict()

        #edge_target_dict - maps each edge index to its target id: 
        self.edge_target_dict=OrderedDict()
 
        #edge_type_dict - maps each edge index to its type id: 
        self.edge_type_dict=OrderedDict()
 
        return None 

    #------------------------------------------------------------------#
    def __print_prolog(self):
        print ("*************************************************")
        print ("Running RACIPE !!!!!")
        print ('Following is a list of user supplied inputs:')
        for arg in self.parsed_dict:
            print (arg + ': ' + str(self.parsed_dict[arg]))
        print ('*********************************************\n')
        return None

    #---------------------------------------------------------------------#
    def __create_params_outfiles(self,node_dict): 
        ''' 
        This method creates FIVE output files in the home directory: 
        '''

        #create a dictionary of edges: source-target
        edge_name_dict=OrderedDict()
        for idx in self.master_dict.keys(): 
            edge_name_dict[idx]=self.master_dict[idx][0]+"-"+self.master_dict[idx][1]  

        #for idx in edge_name_dict.keys(): 
        #    print(self.master_dict[idx]) 
        #    print(edge_name_dict[idx]) 
        #sys.exit(0)

        #write header line in TWO files for node parameters:
        for fname in self.fname_dict_nodeparams.values(): 
            #print(fname)
            #write generic information type
            fh_params=open(fname,'w')
            for header_info in HEADER_1_PARAMS:
                fh_params.write(header_info+'\t')
            #write node names
            for source in node_dict.keys():
                fh_params.write(source+'\t')
            #write a new line    
            fh_params.write('\n')
            fh_params.close()

        #write header line in THREE files for edge parameters:
        for fname in self.fname_dict_edgeparams.values(): 
            #print(fname)
            #write generic information type
            fh_params=open(fname,'w')
            for header_info in HEADER_1_PARAMS:
                fh_params.write(header_info+'\t')
            #write edge name:
            #for source in node_dict.keys():
            for idx in edge_name_dict.keys(): 
                fh_params.write(edge_name_dict[idx]+'\t')
            #write a new line    
            fh_params.write('\n')
            fh_params.close()
        return None 

    #---------------------------------------------------------------------#
    def __create_outfiles_forC(self,node_dict): 
        ''' 
        This method creates THREE output files in the home directory: 
        '''

        #write header line in the .states.txt file:
        fh_ss=open(self.fname_states,'w')
        #write generic information type
        for header_info in HEADER_1_STATES:
            fh_ss.write(header_info+'\t')
        #write node names
        for source in node_dict.keys():
            fh_ss.write(source+'\t')
        #write a new line
        fh_ss.write('\n')
        fh_ss.close()


        #write header line in the .states.cnoise.txt file:
        #fh_ss=open(self.fname_states_cnoise,'w')
        #write generic information type
        #for header_info in HEADER_1_STATES:
        #    fh_ss.write(header_info+'\t')
        #write node names
        #for source in node_dict.keys():
        #    fh_ss.write(source+'\t')
        #write a new line    
        #fh_ss.write('\n')
        #fh_ss.close()


        #write header line in the .states.annealed.txt file:
        #fh_ss=open(self.fname_states_annealed,'w')
        #write generic information type
        #for header_info in HEADER_ANNEALED_STATES:
        #    fh_ss.write(header_info+'\t')
        #write node names
        #for source in node_dict.keys():
        #    fh_ss.write(source+'\t')
        #write a new line    
        #fh_ss.write('\n')
        #fh_ss.close()

        #write header line in the .limitcycles.txt file:
        fh_lc=open(self.fname_limitcycles,'w')
        for header_info in HEADER_1_LIMITCYCLES:
            fh_lc.write(header_info+'\t')
        for source in node_dict.keys():
            fh_lc.write(source+'\t')
        fh_lc.write('\n')
        fh_lc.close()

        #write header line in the .summary.txt file:
        fh_summary=open(self.fname_summary,'w')
        for k in self.config_dict.get('SUMMARY_FILE_HEADER').split(","):
            fh_summary.write(k+'\t')
        fh_summary.write('\n')
        fh_summary.close()
        return None 

    #------------------------------------------------------------------#
    def __print_epilog(self):
        if os.path.exists(self.score_fname):
            print(bcolors.OKGREEN + 'RACIPE generated models are saved ' +\
                                    'in the following output file: ' + \
                  bcolors.ENDC)
            print('         ' + basename(self.score_fname))
        else:
            print(bcolors.WARNING + 'No output file ' +\
                                    'is created' + bcolors.ENDC)
        return None

    #------------------------------------------------------------------#
    def __print_dict(self, d):
        for k,v in d.items():
            print(k, ':', v) 
        print('\n')
        return None 

    #------------------------------------------------------------------#
    def __print_nested_defaultdict(self,node_dict,d):
        for k in node_dict.keys():
            v=d[k]
            print(k,': ') 
            for (idx,e) in v.items():
                print('  ',idx,': ', e)
        print('\n')
        return None 

    #------------------------------------------------------------------#
    def get_fname_dict_simu(self):
        return self.fname_dict_simu 

    #------------------------------------------------------------------#
    def get_parsed_dict(self):
        return self.parsed_dict


    #------------------------------------------------------------------#
    def get_fname_nodeprs(self):
        return self.fname_dict_nodeprs

    #------------------------------------------------------------------#
    def get_fname_edgeprs(self):
        return self.fname_dict_edgeprs 

    #------------------------------------------------------------------#
    def get_fname_nodeparams(self):
        return self.fname_dict_nodeparams

    #------------------------------------------------------------------#
    def get_fname_edgeparams(self):
        return self.fname_dict_edgeparams 

    #------------------------------------------------------------------#
    def get_work_dir(self):
        return self.work_dir 

    #------------------------------------------------------------------#
    def get_config_dict(self):
        return self.config_dict 

    #------------------------------------------------------------------#
    def get_node_dict(self):
        return self.node_dict 

    #------------------------------------------------------------------#
    def get_node_id_dict(self):
        return self.node_id_dict 

    #------------------------------------------------------------------#
    def get_edge_source_dict(self):
        return self.edge_source_dict 

    def get_edge_target_dict(self):
        return self.edge_target_dict 

    def get_edge_type_dict(self):
        return self.edge_type_dict 
    #------------------------------------------------------------------#
    def get_target_dict(self):
        return self.target_dict 

    #------------------------------------------------------------------#
    def get_source_dict(self):
        return self.source_dict 

    #------------------------------------------------------------------#
    def get_user_seed(self):
        return self.user_seed

    #------------------------------------------------------------------#
    def get_master_dict(self):
        return self.master_dict 

    #------------------------------------------------------------------#
    def get_params_fname(self):
        return self.params_fname

    #------------------------------------------------------------------#
    def get_tpo_fname(self):
        return self.tpo_fname

    #------------------------------------------------------------------#
    def get_LCtrace_fname(self):
        return self.LCtrace_fname

    #------------------------------------------------------------------#
    def build_network(self):
        self.node_dict,self.edge_dict,self.source_dict,self.target_dict,self.master_dict=\
             upload_topology(open(self.tpo_fname,'r'),self.config_dict,
                                self.node_dict,self.source_dict,
                                self.target_dict,self.master_dict) 

        self.node_id_dict,self.edge_source_dict,self.edge_target_dict,self.edge_type_dict=\
             build_map_dict(self.node_dict,self.master_dict,self.node_id_dict,
                               self.edge_source_dict,self.edge_target_dict,
                               self.edge_type_dict)

        return None 


    #------------------------------------------------------------------#
    def delete_intermediate_files(self):
        '''
        This method deletes the intermediate files created by the program 
        '''
        for k,fname in self.fname_dict_nodeprs.items(): 
            os.remove(fname)

        for k,fname in self.fname_dict_edgeprs.items(): 
            os.remove(fname)

        for k,fname in self.fname_dict_nodeparams.items(): 
            os.remove(fname)

        for k,fname in self.fname_dict_edgeparams.items(): 
            os.remove(fname)
        return None

    #------------------------------------------------------------------#
    def delete_intermediate_files_2(self):
        '''
        This method deletes the intermediate files created by the program 
        '''
        #for k,fname in self.fname_dict_nodeprs.items(): 
        #    os.remove(fname)

        #for k,fname in self.fname_dict_edgeprs.items(): 
        #    os.remove(fname)

        for k,fname in self.fname_dict_nodeparams.items(): 
            os.remove(fname)

        for k,fname in self.fname_dict_edgeparams.items(): 
            os.remove(fname)
        return None

    #------------------------------------------------------------------#
    def delete_intermediate_prs_files(self):
        '''
        This method deletes the intermediate files created by the program 
        '''
        for k,fname in self.fname_dict_nodeprs.items(): 
            os.remove(fname)

        for k,fname in self.fname_dict_edgeprs.items(): 
            os.remove(fname)

        return None


    #------------------------------------------------------------------#
    def process_network(self, USER_SEED, SEED):
        '''
        This method invokes other methods to upload the network topology.
        '''
        
        #clib=ctypes.cdll.LoadLibrary('./simulation_clib.so')
        clib=ctypes.cdll.LoadLibrary('./pertParamSim.so')

        # Print the wellcome message:
        #self.__print_prolog()

        self.node_dict,self.edge_dict,self.source_dict,self.target_dict,self.master_dict=\
             upload_topology(open(self.tpo_fname,'r'),self.config_dict,
                                self.node_dict,self.source_dict,
                                self.target_dict,self.master_dict) 
 
        #self.__print_nested_defaultdict(self.node_dict, self.source_dict) 
        #self.__print_dict(self.node_dict) 
        #self.__print_nested_defaultdict(self.node_dict, self.target_dict) 

        #construct arrays for passing to C program: 
        self.node_id_dict,self.edge_source_dict,self.edge_target_dict,self.edge_type_dict=\
             build_map_dict(self.node_dict,self.master_dict,self.node_id_dict,
                               self.edge_source_dict,self.edge_target_dict,
                               self.edge_type_dict)

        mode=self.parsed_dict['mode']
        
        # Creates racipe.cfg file: 
        if mode=='C': 
            #create config file:
            #ch.create_config(open(self.cfg_fname,"w"))
            ch.create_config(self.cfg_fname)
            print('racipe.cfg file is generated and saved in working folder')

        elif mode=='A': 
            #set seed for random number generator:
            clib.set_seed(ctypes.c_int(USER_SEED), ctypes.c_int(SEED))
            print('Estimating threshold ranges ...')
            self.target_dict=te.estimateThreshold(self.config_dict,
                                                  self.node_dict,
                                                  self.source_dict,
                                                  self.target_dict)
            #print('after TSH estimation:')
            #self.__print_nested_defaultdict(self.node_dict, self.target_dict) 

            print('Saving parameter ranges ...')
            writeParameterRanges(open(self.prs_fname,'w'),self.config_dict,
                                    self.node_dict,self.target_dict)
            writeParameterRanges_5(self.fname_dict_nodeprs, self.fname_dict_edgeprs,
                                   self.node_dict,self.edge_dict,self.target_dict)

            #create files to save in C program
            self.__create_outfiles_forC(self.node_dict)
            #self.__create_params_outfiles(self.node_dict)
            #sys.exit(0)
        elif mode=='P':
            #upload parameter ranges from the .prs file:
            print('Uploading parameter ranges ...')
            self.target_dict=uploadParameterRanges(open(self.prs_fname,'r'),
                                                    self.config_dict,
                                                    self.node_dict,
                                                    self.target_dict)
            
            print('Saving uploaded parameter ranges to .bak file ...')
            writeParameterRanges(open(self.prs_fname_bak,'w'),self.config_dict,
                                    self.node_dict,self.target_dict)

           #create files to save in C program
            self.__create_outfiles_forC(self.node_dict)
            #self.__create_params_outfiles(self.node_dict)
            #self.delete_intermediate_files()
            #sys.exit(0)
        elif mode=='T': 
            #set seed for random number generator:
            clib.set_seed(ctypes.c_int(USER_SEED), ctypes.c_int(SEED))

            print('Estimating threshold ranges ...')
            self.target_dict=te.estimateThreshold(self.config_dict,
                                                  self.node_dict,
                                                  self.source_dict,
                                                  self.target_dict)
            print('Saving parameter ranges ...')
            writeParameterRanges(open(self.prs_fname,'w'),self.config_dict,
                                    self.node_dict,self.target_dict)
            writeParameterRanges_5(self.fname_dict_nodeprs, self.fname_dict_edgeprs,
                                   self.node_dict,self.edge_dict,self.target_dict)
            #self.delete_intermediate_prs_files()

        elif mode=='S': 
            # threshold half-functional rule test
            da.calProb(self.node_dict,self.master_dict,self.exp_fname, 
                       self.params_input_fname, self.edge_stat_fname)
            #print(self.tpo_fname) 
            #print(self.exp_fname)
            #print(self.params_input_fname)
            print(self.edge_stat_fname)

        # stochastic with both I1 and I2:
        elif mode=='SCP': 
            #upload parameter ranges from the .prs file:
            print('Uploading parameter ranges ...')
            self.target_dict=uploadParameterRanges(open(self.prs_fname,'r'),
                                                    self.config_dict,
                                                    self.node_dict,
                                                    self.target_dict)
            
            print('Saving uploaded parameter ranges to .bak file ...')
            writeParameterRanges(open(self.prs_fname_bak,'w'),self.config_dict,
                                    self.node_dict,self.target_dict)

           #create files to save in C program
            self.__create_outfiles_forC(self.node_dict)
            #self.__create_params_outfiles(self.node_dict)

        # stochastic with only I1:
        elif mode=='SC': 
            #set seed for random number generator:
            clib.set_seed(ctypes.c_int(USER_SEED), ctypes.c_int(SEED))
            print('Estimating threshold ranges ...')
            self.target_dict=te.estimateThreshold(self.config_dict,
                                                  self.node_dict,
                                                  self.source_dict,
                                                  self.target_dict)
            #print('after TSH estimation:')
            #self.__print_nested_defaultdict(self.node_dict, self.target_dict) 

            print('Saving parameter ranges ...')
            writeParameterRanges(open(self.prs_fname,'w'),self.config_dict,
                                    self.node_dict,self.target_dict)
            writeParameterRanges_5(self.fname_dict_nodeprs, self.fname_dict_edgeprs,
                                   self.node_dict,self.edge_dict,self.target_dict)

           #create files to save in C program
            self.__create_outfiles_forC(self.node_dict)
            #self.__create_params_outfiles(self.node_dict)
            #sys.exit(0)        
            
        # simulated annealing when both I1 and I2 are supplied:
        elif mode=='SAP': 
            #upload parameter ranges from the .prs file:
            print('Uploading parameter ranges ...')
            self.target_dict=uploadParameterRanges(open(self.prs_fname,'r'),
                                                    self.config_dict,
                                                    self.node_dict,
                                                    self.target_dict)
            
            print('Saving uploaded parameter ranges to .bak file ...')
            writeParameterRanges(open(self.prs_fname_bak,'w'),self.config_dict,
                                    self.node_dict,self.target_dict)

           #create files to save in C program
            self.__create_outfiles_forC(self.node_dict)
            #self.__create_params_outfiles(self.node_dict)

        # simulated annealing when only I1 is supplied:
        elif mode=='SA': 
            #set seed for random number generator:
            clib.set_seed(ctypes.c_int(USER_SEED), ctypes.c_int(SEED))
            print('Estimating threshold ranges ...')
            self.target_dict=te.estimateThreshold(self.config_dict,
                                                  self.node_dict,
                                                  self.source_dict,
                                                  self.target_dict)
            #print('after TSH estimation:')
            #self.__print_nested_defaultdict(self.node_dict, self.target_dict) 

            print('Saving parameter ranges ...')
            writeParameterRanges(open(self.prs_fname,'w'),self.config_dict,
                                    self.node_dict,self.target_dict)
            writeParameterRanges_5(self.fname_dict_nodeprs, self.fname_dict_edgeprs,
                                   self.node_dict,self.edge_dict,self.target_dict)

           #create files to save in C program
            self.__create_outfiles_forC(self.node_dict)
            #self.__create_params_outfiles(self.node_dict)
            #sys.exit(0)



        else: 
            print('Unrecognized mode. Exiting ...')
            sys.exit(0)

        #Print the summary of running this program:
        #self.__print_epilog()
        return mode 


#**********************************************************************#
def upload_topology(fh_tpo,config_dict,node_dict, 
                           source_dict,target_dict,master_dict): 
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
    edge_dict=OrderedDict()
    #skip header line:
    next(fh_tpo)
    #import network topology into relevant data structures: 
    index=0
    for line in fh_tpo:
        #tpo[0]:source,tpo[1]:target,tpo[2]:type
        tpo=line.strip().split()

        #update node_dict:
        if tpo[0] not in node_dict.keys(): 
            node_dict[tpo[0]]=[FAST_NODE,(float(config_dict['MPR_MIN']),
                                          float(config_dict['MPR_MAX'])),
                                         (float(config_dict['DNR_MIN']),
                                          float(config_dict['DNR_MAX']))]
        if tpo[1] not in node_dict.keys():
            node_dict[tpo[1]]=[FAST_NODE,(float(config_dict['MPR_MIN']),
                                          float(config_dict['MPR_MAX'])),
                                         (float(config_dict['DNR_MIN']),
                                          float(config_dict['DNR_MAX']))]

        # update node_dict for tpo[1] (target), if the edge is SLOW:
        if int(tpo[2]) in SLOW_EDGE_TYPES:
            node_dict[tpo[1]]=[SLOW_NODE,(float(config_dict['MPR_MIN']),
                                          float(config_dict['MPR_MAX'])),
                                         (float(config_dict['DNR_MIN']),
                                          float(config_dict['DNR_MAX']))]
        #update source_dict:
        source_dict[tpo[1]][index]=[tpo[0],int(tpo[2])]

        #update target_dict:
        hco_range=(int(config_dict['HCO_MIN']),
                   int(config_dict['HCO_MAX']))
        if(int(tpo[2])==3 or int(tpo[2])==4):  
            fch_range=(float(config_dict['FCH_MIN']),
                       float(config_dict['FCH_MAX_DEG']))
        else: 
            fch_range=(float(config_dict['FCH_MIN']),
                       float(config_dict['FCH_MAX']))
        target_dict[tpo[0]][index]=[tpo[1],int(tpo[2]),
                                    hco_range,fch_range]
        #update master_dict:
        master_dict[index]=[tpo[0],tpo[1],tpo[2]]
        edge_dict[index]=tpo[0]+'-'+tpo[1]
        index+=1
    for node in node_dict.keys(): 
        if not bool(source_dict[node]):
            node_dict[node]=[SLOW_NODE,(float(config_dict['MPR_MIN']),
                                        float(config_dict['MPR_MAX'])),
                                       (float(config_dict['DNR_MIN']),
                                        float(config_dict['DNR_MAX']))]
    return node_dict,edge_dict,source_dict,target_dict,master_dict 

#------------------------------------------------------------------#
def uploadParameterRanges(fh_prs,config_dict,node_dict,target_dict):
    #skip the header line:
    next(fh_prs)
    
    node_attributes=config_dict['NODE_ATTRIBUTES'].strip().split(',')
    edge_attributes=config_dict['EDGE_ATTRIBUTES'].strip().split(',')
    #upload node parameters:
    mpr_tag=node_attributes[0]+'_'
    dnr_tag=node_attributes[1]+'_'
    while True: 
        line1=fh_prs.readline()
        if line1=='\n': 
            break
        #node_idx,par_name,mpr_min,mpr_max,reg_type=line1.strip().split('\t')
        node_idx,par_name,mpr_min,mpr_max,node_type=line1.strip().split('\t')
        node=par_name.split(mpr_tag)[1]
        line2=fh_prs.readline()
        node_idx,par_name,dnr_min,dnr_max,node_type=line2.strip().split('\t')
        node=par_name.split(dnr_tag)[1]
        #node_dict[node]=[(float(mpr_min),float(mpr_max)),
        #                      (float(dnr_min),float(dnr_max))] 
        node_dict[node]=[int(node_type),
                         (float(mpr_min),float(mpr_max)),
                         (float(dnr_min),float(dnr_max))] 
    #load edge parameters:
    tsh_tag=edge_attributes[0]+'_'
    hco_tag=edge_attributes[1]+'_'
    fch_tag=edge_attributes[2]+'_'
    while True: 
        line1=fh_prs.readline()
        if line1=="\n":
            continue
        elif not line1: 
            break
        edge_idx,par_name,tsh_min,tsh_max,reg_type=line1.strip().split('\t')
        edge=par_name.split(tsh_tag)[1]
        source,target=edge.split('_TO_')

        line2=fh_prs.readline()
        edge_idx,par_name,hco_min,hco_max,reg_type=line2.strip().split('\t')

        line3=fh_prs.readline()
        edge_idx,par_name,fch_min,fch_max,reg_type=line3.strip().split('\t')
        target_dict[source][int(edge_idx)]=[target,reg_type,
                                                (float(tsh_min),float(tsh_max)),
                                                (int(hco_min),int(hco_max)),
                                                (float(fch_min),float(fch_max))]
    return target_dict

#------------------------------------------------------------------#
def build_map_dict(node_dict,master_dict,node_id_dict,edge_source_dict,
                   edge_target_dict,edge_type_dict):
    '''
    This method builds dictionary to map node and edges to unique identifiers. 
    '''

    node_id=0
    for X in node_dict.keys(): 
        node_id_dict[X]=node_id
        node_id+=1 

    for idx,e in master_dict.items(): 
        edge_source_dict[idx]=node_id_dict[e[0]]
        edge_target_dict[idx]=node_id_dict[e[1]]
        edge_type_dict[idx]=e[2]

    #return None
    return node_id_dict,edge_source_dict,edge_target_dict,edge_type_dict

#-----------------------------------------------------------------------------#
def open_prs_files(fname_dict_nodeprs, fname_dict_edgeprs): 
    ''' 
    This method opens output files for writing solutions. 
    It creats MAX_STABLE_STATES number of file handles in the write mode, 
    places in a dictionary, and returns the dictionary. 
    '''
    fh_dict_nodeprs=OrderedDict()
    fh_dict_edgeprs=OrderedDict()
    for (k,fname) in fname_dict_nodeprs.items(): 
        fh_dict_nodeprs[k]=open(fname,'w')

    for (k,fname) in fname_dict_edgeprs.items(): 
        fh_dict_edgeprs[k]=open(fname,'w')
    return fh_dict_nodeprs,fh_dict_edgeprs

#----------------------------------------------------------------------#
def close_prs_files(fh_dict_nodeprs,fh_dict_edgeprs):
    for idx in fh_dict_nodeprs.keys():
        fh_dict_nodeprs[idx].close()

    for idx in fh_dict_edgeprs.keys():
        fh_dict_edgeprs[idx].close()
    return None 

#------------------------------------------------------------------#
def writeParameterRanges_5(fname_dict_nodeprs, fname_dict_edgeprs,
                           node_dict,edge_dict, target_dict):
    #write parameter ranges in FIVE separate files:

    (fh_dict_nodeprs, fh_dict_edgeprs)=open_prs_files(fname_dict_nodeprs,fname_dict_edgeprs) 
    (fh_mpr_prs, fh_dnr_prs)=fh_dict_nodeprs.values()
    (fh_tsh_prs,fh_hco_prs,fh_fch_prs)=fh_dict_edgeprs.values()

    fh_mpr_prs.write('\n')
    fh_mpr_prs.write('\t')
    for source in node_dict.keys():
        fh_mpr_prs.write(source+'\t')

    fh_dnr_prs.write('\t')
    for source in node_dict.keys():
        fh_dnr_prs.write(source+'\t')

    outstr_mpr_min='\n'+'MIN'
    outstr_mpr_max='\n'+'MAX'
    outstr_dnr_min='\n'+'MIN'
    outstr_dnr_max='\n'+'MAX'
    for node in node_dict.keys():
        (node_type,mpr_range,dnr_range)=node_dict[node]
        outstr_mpr_min+='\t'+str('%10.6f'%(float(mpr_range[0])))
        outstr_mpr_max+='\t'+str('%10.6f'%(float(mpr_range[1])))
        outstr_dnr_min+='\t'+str('%10.6f'%(float(dnr_range[0])))
        outstr_dnr_max+='\t'+str('%10.6f'%(float(dnr_range[1])))
    fh_mpr_prs.write(outstr_mpr_min)
    fh_mpr_prs.write(outstr_mpr_max)
    fh_dnr_prs.write(outstr_dnr_min)
    fh_dnr_prs.write(outstr_dnr_max)

    #write edge parameters:
    fh_tsh_prs.write('\t')
    for edge_name in edge_dict.values():
        fh_tsh_prs.write(edge_name+'\t')
        
    fh_hco_prs.write('\t')
    for edge_name in edge_dict.values():
        fh_hco_prs.write(edge_name+'\t')

    fh_fch_prs.write('\t')
    for edge_name in edge_dict.values():
        fh_fch_prs.write(edge_name+'\t')

    outstr_tsh_min='\n'+'MIN'
    outstr_tsh_max='\n'+'MAX'
    outstr_hco_min='\n'+'MIN'
    outstr_hco_max='\n'+'MAX'
    outstr_fch_min='\n'+'MIN'
    outstr_fch_max='\n'+'MAX'
    for source in node_dict.keys():
        for (idx,e) in target_dict[source].items():
            target,reg_type,tsh_range,hco_range,fch_range=e
            #construct output strings:
            outstr_tsh_min+='\t'+str('%10.6f'%(tsh_range[0]))
            outstr_tsh_max+='\t'+str('%10.6f'%(tsh_range[1]))

            outstr_hco_min+='\t'+str('%10d'%(hco_range[0]))
            outstr_hco_max+='\t'+str('%10d'%(hco_range[1]))

            outstr_fch_min+='\t'+str('%10.6f'%(fch_range[0]))
            outstr_fch_max+='\t'+str('%10.6f'%(fch_range[1]))

    fh_tsh_prs.write(outstr_tsh_min)
    fh_tsh_prs.write(outstr_tsh_max)

    fh_hco_prs.write(outstr_hco_min)
    fh_hco_prs.write(outstr_hco_max)

    fh_fch_prs.write(outstr_fch_min)
    fh_fch_prs.write(outstr_fch_max)

    close_prs_files(fh_dict_nodeprs, fh_dict_edgeprs)
    return None

#------------------------------------------------------------------#
def writeParameterRanges(fh_prs,config_dict,node_dict,target_dict): 
    #write the header line to the file:
    tmp_list=config_dict['PARAMETER_FILE_HEADER'].strip().split(',')
    header_str=''
    for st in tmp_list: 
        header_str+=st+'\t'
    header_str+='\n' 
    fh_prs.write(header_str) 
   
    node_attributes=config_dict['NODE_ATTRIBUTES'].strip().split(',')
    edge_attributes=config_dict['EDGE_ATTRIBUTES'].strip().split(',')

    #write parameter ranges for nodes:
    idx=0
    for node in node_dict.keys():
        #mpr_range=node_dict[node][0] 
        #dnr_range=node_dict[node][1] 
        (node_type,mpr_range,dnr_range)=node_dict[node]
        #construct output strings:
        outstr_mpr=str(idx)+'\t'+node_attributes[0]+'_'+node+'\t'+\
                            str('%10.6f'%(float(mpr_range[0])))+'\t'+\
                            str('%10.6f'%(float(mpr_range[1])))+'\t'+\
                            str(node_type)+'\n'
        outstr_dnr=str(idx)+'\t'+node_attributes[1]+'_'+node+'\t'+\
                            str('%10.6f'%(float(dnr_range[0])))+'\t'+\
                            str('%10.6f'%(float(dnr_range[1])))+'\t'+\
                            str(node_type)+'\n'
        fh_prs.write(outstr_mpr)
        fh_prs.write(outstr_dnr)
        idx+=1
    #write a new line to the output file as a separator:
    fh_prs.write("\n")

    #write parameter ranges for regulations/edges:
    for source in node_dict.keys():
        for (idx,e) in target_dict[source].items():
            target,reg_type,tsh_range,hco_range,fch_range=e
            #'_'+source+'_TO_'+target:
            link='_'+source+'_TO_'+target
            #construct output strings:
            outstr_tsh=str(idx)+'\t'+edge_attributes[0]+link+'\t'+\
                                     str('%10.6f'%(tsh_range[0]))+'\t'+\
                                     str('%10.6f'%(tsh_range[1]))+'\t'+\
                                     str(reg_type)+'\n'
            outstr_hco=str(idx)+'\t'+edge_attributes[1]+link+'\t'+\
                                     str('%10d'%(hco_range[0]))+'\t'+\
                                     str('%10d'%(hco_range[1]))+'\t'+\
                                     str(reg_type)+'\n'
            outstr_fch=str(idx)+'\t'+edge_attributes[2]+link+'\t'+\
                                     str('%10.6f'%(fch_range[0]))+'\t'+\
                                     str('%10.6f'%(fch_range[1]))+'\t'+\
                                     str(reg_type)+'\n'
            fh_prs.write(outstr_tsh)
            fh_prs.write(outstr_hco)
            fh_prs.write(outstr_fch)
        #write a new line after all edges between a (source,target) pair:
        fh_prs.write("\n")
    return None


#**********************************************************************#
if __name__=='__main__':
   print("This file contains modules to read and write the network" +\
         " topology and parameters.")
   sys.exit(0)
