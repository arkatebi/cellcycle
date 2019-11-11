#/usr/bin/env python
''''
    The program racipe.py invokes the methods defined in this file.
'''

import numpy as np
#import matplotlib.pyplot as plt

import numpy as np
import sys
import ctypes
np.random.seed(1)

clib=ctypes.cdll.LoadLibrary('./simulation_clib.so')
clib.randu.argtypes=(ctypes.c_double,ctypes.c_double)
clib.randu.restype=ctypes.c_double

#constants:
SLOW_NODE=1
FAST_NODE=2
NODE_TYPES=[SLOW_NODE,FAST_NODE]

SLOW_EDGE_TYPES=[1,2,3,4]
FAST_EDGE_TYPES=[5,6]

EXCITATION_TYPES=[1,3,5]
INHIBITION_TYPES=[4,6] 

#----------------------------------------------------------------------#
def eval_shiftedHill_fn(X,X0,nX,lamb):
    '''
    This method calculates and returns the effect of shifted Hill function.
    '''
    return lamb+(1.0-lamb)/(1.0+(X/X0)**nX)

#----------------------------------------------------------------------#
def estimateThreshold_noRegulators(config_dict): 
    '''
    config_dict: OrderedDict
    where: 
         Key: name of Configuration attribute
         Value: the default value for the attribute
    '''
    global clib
    import ctypes

    X=[]
    for i_sim in range(int(config_dict['NUM_SIM_THRESHOLD'])): 
        g=clib.randu(ctypes.c_double(float(config_dict['MPR_MIN'])),
                     ctypes.c_double(float(config_dict['MPR_MAX'])))
        k=clib.randu(ctypes.c_double(float(config_dict['DNR_MIN'])),
                     ctypes.c_double(float(config_dict['DNR_MAX']))) 
        X.append(g/k)
    #X_median=np.median(X)
    #logX=np.log2(X)
    #plt.hist(logX,bins='auto',normed=1)
    #plt.hist(logX,bins='auto')
    #plt.show()
    return np.median(X)

def estimateThreshold_noRegulators_npRandom(config_dict): 
    '''
    config_dict: OrderedDict
    where: 
         Key: name of Configuration attribute
         Value: the default value for the attribute
    '''
    global clib
    import ctypes

    X=[]
    for i_sim in range(int(config_dict['NUM_SIM_THRESHOLD'])): 
        g=np.random.uniform(float(config_dict['MPR_MIN']),
                    float(config_dict['MPR_MAX']))
        k=np.random.uniform(float(config_dict['DNR_MIN']),
                    float(config_dict['DNR_MAX'])) 
        X.append(g/k)
    return np.median(X)

#----------------------------------------------------------------------#
def estimateThreshold_withRegulators(config_dict,inwardEdges): 
    '''
    config_dict: OrderedDict
    where: 
         Key: name of Configuration attribute
         Value: the default value for the attribute
    ---
    inwardEdges: defaultdict(list)
    Key: index
    Value: a list of two entries (source, regulation type)
    '''
    global clib
    import ctypes

    TSH_singleGene=estimateThreshold_noRegulators(config_dict)
    #X for storing threshold values for all simulations:
    X=[]
    for i_sim in range(int(config_dict['NUM_SIM_THRESHOLD'])): 
        g=clib.randu(ctypes.c_double(float(config_dict['MPR_MIN'])),
                     ctypes.c_double(float(config_dict['MPR_MAX'])))
        k=clib.randu(ctypes.c_double(float(config_dict['DNR_MIN'])),
                     ctypes.c_double(float(config_dict['DNR_MAX']))) 
        tmp_tsh=g/k
        for (idx,e) in inwardEdges.items():
            gY=clib.randu(ctypes.c_double(float(config_dict['MPR_MIN'])),
                     ctypes.c_double(float(config_dict['MPR_MAX'])))
            kY=clib.randu(ctypes.c_double(float(config_dict['DNR_MIN'])),
                     ctypes.c_double(float(config_dict['DNR_MAX']))) 


            Y0=clib.randu(ctypes.c_double(float(config_dict['TSH_SCALE_FACTOR_MIN'])*\
                                               TSH_singleGene),
                           ctypes.c_double(float(config_dict['TSH_SCALE_FACTOR_MAX'])*\
                                               TSH_singleGene))
            nY=int(clib.randu(ctypes.c_double(int(config_dict['HCO_MIN'])),
                              ctypes.c_double(int(config_dict['HCO_MAX'])+1)))

            (Y,reg_type)=e
            #if int(reg_type)==int(config_dict['REG_EXCITATORY']): #excitatory link
            if int(reg_type)==1 or int(reg_type)==5: #activation by TF or protein
                lY=clib.randu(ctypes.c_double(float(config_dict['FCH_MIN'])),
                              ctypes.c_double(float(config_dict['FCH_MAX'])))
                hill_impact=eval_shiftedHill_fn(gY/kY,Y0,nY,lY)/lY
                tmp_tsh=tmp_tsh*hill_impact
            #elif int(reg_type)==int(config_dict['REG_INHIBITORY']): #inhibitory link
            elif int(reg_type)==2 or int(reg_type)==6: #inhibition by TF or protein
                lY=1.0/clib.randu(ctypes.c_double(float(config_dict['FCH_MIN'])),
                                  ctypes.c_double(float(config_dict['FCH_MAX'])))
                hill_impact=eval_shiftedHill_fn(gY/kY,Y0,nY,lY)
                tmp_tsh=tmp_tsh*hill_impact
            elif int(reg_type)==4: #degradation inhibition => activation in degradation rate
                lY=clib.randu(ctypes.c_double(float(config_dict['FCH_MIN'])),
                                  ctypes.c_double(float(config_dict['FCH_MAX_DEG'])))
                hill_impact=eval_shiftedHill_fn(gY/kY,Y0,nY,lY)
                tmp_tsh=tmp_tsh/hill_impact
            elif int(reg_type)==3: #degradation activation => inhibition in degradation rate
                lY=1.0/clib.randu(ctypes.c_double(float(config_dict['FCH_MIN'])),
                              ctypes.c_double(float(config_dict['FCH_MAX_DEG'])))
                #hill_impact=eval_shiftedHill_fn(gY/kY,Y0,nY,lY)/lY
                hill_impact=eval_shiftedHill_fn(gY/kY,Y0,nY,lY)
                tmp_tsh=tmp_tsh/hill_impact

            else:
                print('unknown regulation type found.')
                print('program exiting...')
                sys.exit(0)
            #tmp_tsh=tmp_tsh*hill_impact
        X.append(tmp_tsh)

    #X_median=np.median(X)
    #logX=np.log2(X)
    #plt.hist(logX,bins='auto',normed=1)
    #plt.hist(logX,bins='auto')
    #plt.show()
    return np.median(X)

#------------------------------------------------------------------#
def estimateThreshold(config_dict,node_dict,
                      source_dict,target_dict):
    for X in node_dict.keys():
        #print('X: ', X)  
        inwardEdge_dict=source_dict[X]
        outwardEdge_dict=target_dict[X]
        if inwardEdge_dict and outwardEdge_dict:
            medianTSH=estimateThreshold_withRegulators(config_dict,
                                                       inwardEdge_dict)
            #print('node type 1: ', X)
            for (idx,e) in outwardEdge_dict.items(): 
                #update the entry in self.target_dict:
                tsh_min=float(config_dict['TSH_SCALE_FACTOR_MIN'])*medianTSH
                tsh_max=float(config_dict['TSH_SCALE_FACTOR_MAX'])*medianTSH
                #place threshold in the edge attributes:
                #target_dict[X][idx]=e[:2]+[(tsh_min,tsh_max)]+e[2:]
                target_dict[X][idx]=e[:2]+\
                                     [(float('%10.6f'%tsh_min),
                                       float('%10.6f'%tsh_max))]+\
                                    e[2:]
        #elif len(outwardEdge_dict) and not len(inwardEdge_dict):
        else:
            #print('node type 2: ',X)
            #sys.exit(0)
            medianTSH=estimateThreshold_noRegulators(config_dict)
            for (idx,e) in outwardEdge_dict.items(): 
                #update the entry in self.target_dict:
                tsh_min=float(config_dict['TSH_SCALE_FACTOR_MIN'])*medianTSH
                tsh_max=float(config_dict['TSH_SCALE_FACTOR_MAX'])*medianTSH
                #place threshold in the edge attributes:
                #target_dict[X][idx]=e[:2]+[(tsh_min,tsh_max)]+e[2:]
                target_dict[X][idx]=e[:2]+[(float('%10.6f'%tsh_min),float('%10.6f'%tsh_max))]+e[2:]
        #elif len(inwardEdge_dict) and not len(outwardEdge_dict):
            #medianTSH=te.estimateThreshold_withRegulators(inwardEdge_dict)
            #medianTSH=str('%.2f'%(float(medianTSH)))
            #print('node type 3: ', X)
    return target_dict 


#----------------------------------------------------------------------#
if __name__=='__main__':
   print('This file contains modules for threshold estimation for' +\
         ' network nodes and edges.')
