#/usr/bin/env python 

import sys

from collections import OrderedDict
 
from os.path import basename  


#**********************************************************************#
class Models: 
    def __init__ (self, MODEL_DIR): 

        self.work_dir = MODEL_DIR

        FNAME_PREFIX = "cellcycle."
        # define input file names:
        self.fname_models = self.work_dir + 'models.transitioning.csv'
        self.fname_states = self.work_dir + FNAME_PREFIX + 'states.csv'
        self.fname_mpr = self.work_dir + FNAME_PREFIX + 'mpr.csv'
        self.fname_dnr = self.work_dir + FNAME_PREFIX + 'dnr.csv'
        self.fname_tsh = self.work_dir + FNAME_PREFIX + 'tsh.csv'
        self.fname_fch = self.work_dir + FNAME_PREFIX + 'fch.csv'
        self.fname_hco = self.work_dir + FNAME_PREFIX + 'hco.csv'

        # define dictionaries for storing models
        #self.model_dict = OrderedDict()

        # define dictionaries for storing expressions and parameters:
        self.EXP_dict = OrderedDict()
        self.MPR_dict = OrderedDict()
        self.DNR_dict = OrderedDict()
        self.TSH_dict = OrderedDict()
        self.FCH_dict = OrderedDict()
        self.HCO_dict = OrderedDict()

        self._load_expressions()

        #print(self.EXP_dict.keys())
        #print(self.EXP_dict.values())
        

        self.MPR_dict = self._load_node_params(self.fname_mpr)
        #print(self.MPR_dict.keys())
        #print(self.MPR_dict.values())

        self.DNR_dict = self._load_node_params(self.fname_dnr)
        #print(self.DNR_dict.keys())
        #print(self.DNR_dict.values())

        self.TSH_dict = self._load_edge_params(self.fname_tsh)
        #print(self.TSH_dict.keys())
        #print(self.TSH_dict.values())

        self.FCH_dict = self._load_edge_params(self.fname_fch)
        #print(self.FCH_dict.keys())
        #print(self.FCH_dict.values())

        self.HCO_dict = self._load_edge_params(self.fname_hco)
        #print(self.HCO_dict.keys())
        #print(self.HCO_dict.values())

        return None


    #-----------------------------------------------------------------#
    def _load_expressions(self): 
        #print("loading expressions")

        NO_META_COLS = 4

        fhandle = open(self.fname_states, 'r')  
        # skip the header line: 
        next(fhandle)

        # load expression pairs into the dictionary:
        while True: 
            # first state:
            line = fhandle.readline()

            fields = line.strip().split(sep=',') 
            #if(not representInt(model_no)): break 
            if(not representInt(fields[0])): break 

            #model_no = fields[0]
            model_no = int(fields[0])
            # get cluster no of the state:
            cluster_no_state_1 = fields[3]

            # get state expressions:
            exp_state_1 = fields[(NO_META_COLS):] 

            # second state:
            line = fhandle.readline()
            fields = line.strip().split(sep=',') 
            #model_no = fields[0]

            # get cluster no of the state:
            cluster_no_state_2 = fields[3]

            # get state expressions:
            exp_state_2 = fields[(NO_META_COLS):] 

            # insert the expression pairs intot the dictionary:
            self.EXP_dict[model_no] = (cluster_no_state_1, exp_state_1, \
                    cluster_no_state_2, exp_state_2)

            if not line: break # EOF
        fhandle.close()
        return None


    #-----------------------------------------------------------------#
    def _load_node_params(self, fname_node_params):
        #print("loading node params")

        NO_META_COLS = 1
        param_dict = OrderedDict()

        fhandle = open(fname_node_params, 'r')  
        # skip the header line: 
        next(fhandle)

        # load node params into the dictionary:
        while True:
            # first state:
            line = fhandle.readline()
            fields = line.strip().split(sep=',') 

            #if(not representInt(model_no)): break 
            if(not representInt(fields[0])): break 

            #model_no = fields[0]
            model_no = int(fields[0])
            node_params = fields[(NO_META_COLS):] 

            # insert the params into the dictionary:
            param_dict[model_no] = node_params 

            if not line: break # EOF
        fhandle.close()
        return param_dict


    #-----------------------------------------------------------------#
    def _load_edge_params(self, fname_edge_params):
        #print("loading edge params")

        NO_META_COLS = 1
        param_dict = OrderedDict()

        fhandle = open(fname_edge_params, 'r')  
        # skip the header line: 
        next(fhandle)

        # load node params into the dictionary:
        while True:
            # first state:
            line = fhandle.readline()
            fields = line.strip().split(sep=',') 

            #if(not representInt(model_no)): break 
            if(not representInt(fields[0])): break 
            #model_no = fields[0]
            model_no = int(fields[0])
            node_params = fields[(NO_META_COLS):] 

            # insert the params into the dictionary:
            param_dict[model_no] = node_params 

            if not line: break # EOF
        fhandle.close()
        return param_dict

    #-----------------------------------------------------------------#
    def get_EXP_dict(self):
        return self.EXP_dict 

    #-----------------------------------------------------------------#
    def get_MPR_dict(self):
        return self.MPR_dict 

    #-----------------------------------------------------------------#
    def get_DNR_dict(self):
        return self.DNR_dict 

    #-----------------------------------------------------------------#
    def get_TSH_dict(self):
        return self.TSH_dict 

    #-----------------------------------------------------------------#
    def get_FCH_dict(self):
        return self.FCH_dict 

    #-----------------------------------------------------------------#
    def get_HCO_dict(self):
        return self.HCO_dict 


#---------------------------------------------------------------------#
def representInt(s): 
    try: 
        int(s) 
        return True 
    except ValueError: 
        return False


#**********************************************************************#
#-----------------------------------------------------------------------------#
def set_parameters_bymodel(rcp, MO, model_no):

    #nodeParam_dict: node to (MPR,DNR) mapping:
    #key: node, value: tuple (MPR,DNR)
    nodeParam_dict=OrderedDict()

    #import config_dict, node_dict, source_dict, target_dict:
    config_dict=rcp.get_config_dict()
    node_dict=rcp.get_node_dict()
    source_dict=rcp.get_source_dict()
    target_dict=rcp.get_target_dict()
    master_dict=rcp.get_master_dict()

    # extract model specific parameters
    EXP_dict = MO.get_EXP_dict()

    MPR_dict = MO.get_MPR_dict()
    DNR_dict = MO.get_DNR_dict()

    TSH_dict = MO.get_TSH_dict()
    FCH_dict = MO.get_FCH_dict()
    HCO_dict = MO.get_HCO_dict()

 
    mpr_list = MPR_dict[model_no] 
    dnr_list = DNR_dict[model_no] 

    tsh_list = list(TSH_dict[model_no])
    fch_list = list(FCH_dict[model_no])
    hco_list = list(HCO_dict[model_no])

    #update nodeParam_dict, source_dict, master_dict:
    #with node and edge parameters:
    node_no = 0
    edge_no = 0
    for node in node_dict.keys():
        (node_type,mpr_range,dnr_range)=node_dict[node]
        nodeParam_dict[node]=(node_type, float(mpr_list[node_no]), \
                                         float(dnr_list[node_no]))

        for (idx,e) in target_dict[node].items():
            target,reg_type,tsh_range,hco_range,fch_range=e

            #store the edge parameters in source_dict:
            source_dict[target][idx]=source_dict[target][idx]+\
                                     [float(tsh_list[edge_no]),\
                                     float(hco_list[edge_no]),\
                                     float(fch_list[edge_no])]

            # update master_dict with edge parameters from supplied values:
            master_dict[idx]=(node, target, reg_type, float(tsh_list[edge_no]), \
                                                      float(hco_list[edge_no]), \
                                                      float(fch_list[edge_no]))
            edge_no = edge_no + 1
        node_no = node_no + 1
    return (nodeParam_dict,source_dict,master_dict)



#**********************************************************************#
if __name__ == '__main__':
    print("This file contains modules to read raw data and parameters " +\
            "for a set of models.")
     
    '''
    MODEL_DIR = "/Users/kateba/research/cellcycle-3.stoch/data-raw/earlyG1-midG1/"
    models = Models(MODEL_DIR)
    d = models.get_EXP_dict()
    #print(d.keys())
    #print(d.values())
    for k,v in d.items(): 
        print(k) 
        print(v[0], v[1])
        print(v[2], v[3])
   ''' 
    sys.exit(0)

