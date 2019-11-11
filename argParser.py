#!/usr/bin/env python
'''
    The methods in this module handles the input arguments of the racipe 
software.The entry point of this module is parse_args() method which calls
other methods to collect user supplied arguments, parses and verifies them. 
Description of these methods are the following:
   
    collect_args: This method collects the user supplied arguments and 
        returns them as an aprgparse ArgumentParser object. 

    extract_args: This method puts the user supplied arguments into an 
        ordered dictionary and returns it at the end.

    check_args: This method verifies the correctness of the user supplied
        arguments and puts them into an ordered dictionary which it returns
        at the end. 

    parse_args: This method calls the above methods and returns the final 
        dictionary of the user supplied arguments to the calling point.
'''

import os
import sys
import argparse
import re
from collections import OrderedDict

#-----------------------------------------------------------------------------#    
def collect_args_old():
    """ 
    This method collects the user supplied arguments and returns them 
    as a argparse object.
    """
    parser = argparse.ArgumentParser(description='Generate network models ' + \
        'according to RACIPE method.')
    parser.add_argument('-M', '--mode', default='C', help='Specifies mode ' + \
        'of the program. Default value is A')
    parser.add_argument('-US', '--user_seed', default=0, help='Specifies seed ' + \
        'to multiply the internal seed. Default value is 0')
    parser.add_argument('-S', '--SEED', default=1, help='Changes internal seed.' + \
        'Default value is 1')
    parser.add_argument('-I1', '--tpo', help='Specifies path to the ' + \
        'circuit topology file. This option is mandatory.')
    parser.add_argument('-I2', '--prs', default='', help='Specifies path to the ' + \
        'circuit parameter file. If not supplied, the program will create one.')
    parser.add_argument('-I3', '--cfg', default='', help='Specifies path to the ' + \
        'simulation configuration file. If not supplied, the program will create one.')

    #parser.add_argument('-I4', '--exp', help='Specifies the ' + \
    #    'expression file.')
    parser.add_argument('-I4', '--exp', default='', help='Specifies the ' + \
        'expression file.')

    #parser.add_argument('-I5', '--params', help='Specifies the ' + \
    #    'simulation parameter file. ')
    parser.add_argument('-I5', '--params', default='', help='Specifies the ' + \
        'simulation parameter file. ')

    parser.add_argument('-O', '--output', default='', help='Provides user ' + \
        'an option to specify an output filename prefix. When not ' + \
        'specified, the program will create an output file name.')
    return parser


def collect_args():
    """ 
    This method collects the user supplied arguments and returns them 
    as a argparse object.
    """
    parser = argparse.ArgumentParser(description='Generate network models ' + \
        'according to RACIPE method.')
    parser.add_argument('-M', '--mode', default='C', help='Specifies mode ' + \
        'of the program. Default value is A')
    parser.add_argument('-US', '--user_seed', default=0, help='Specifies seed ' + \
        'to multiply the internal seed. Default value is 0')
    parser.add_argument('-S', '--SEED', default=1, help='Changes internal seed.' + \
        'Default value is 1')
    parser.add_argument('-I1', '--tpo', help='Specifies path to the ' + \
        'circuit topology file. This option is mandatory.')
    parser.add_argument('-I2', '--prs', default='', help='Specifies path to the ' + \
        'circuit parameter file. If not supplied, the program will create one.')
    parser.add_argument('-I3', '--cfg', default='', help='Specifies path to the ' + \
        'simulation configuration file. If not supplied, the program will create one.')

    #parser.add_argument('-I4', '--exp', help='Specifies the ' + \
    #    'expression file.')
    parser.add_argument('-I4', '--exp', default='', help='Specifies the ' + \
        'expression file.')

    #parser.add_argument('-I5', '--params', help='Specifies the ' + \
    #    'simulation parameter file. ')
    parser.add_argument('-I5', '--params', default='', help='Specifies the ' + \
        'simulation parameter file. ')

    parser.add_argument('-KD','--knockdownlist', nargs='*', default=[], help=\
        'This parameter can take the list of knock out genes.')    

    parser.add_argument('-O', '--output', default='', help='Provides user ' + \
        'an option to specify an output filename prefix. When not ' + \
        'specified, the program will create an output file name.')
    return parser

#-----------------------------------------------------------------------------#    
def extract_args_old(args):
    """
     This method builds a dictionary from the user supplied arguments
     and returns the constructed dictionary at the end.
    """
    args_dict = OrderedDict() 
    args_dict['mode'] = args.mode
    args_dict['user_seed'] = args.user_seed
    args_dict['SEED'] = args.SEED
    args_dict['tpo_fname'] = args.tpo
    args_dict['prs_fname'] = args.prs
    args_dict['cfg_fname'] = args.cfg
    args_dict['exp_fname'] = args.exp
    args_dict['params_input_fname'] = args.params   
    args_dict['outfile'] = args.output
    return args_dict

def extract_args(args):
    """
     This method builds a dictionary from the user supplied arguments
     and returns the constructed dictionary at the end.
    """
    args_dict = OrderedDict() 
    args_dict['mode'] = args.mode
    args_dict['user_seed'] = args.user_seed
    args_dict['SEED'] = args.SEED
    args_dict['tpo_fname'] = args.tpo
    args_dict['prs_fname'] = args.prs
    args_dict['cfg_fname'] = args.cfg
    args_dict['exp_fname'] = args.exp
    args_dict['params_input_fname'] = args.params   
    args_dict['knockdownlist'] = args.knockdownlist
    args_dict['outfile'] = args.output
    return args_dict

#-----------------------------------------------------------------------------#    
def check_args(args_dict,parser):
    """
    This method checks the user arguments for consistency. It builds a new 
    dictionary from these arguments and finally returns this newly created 
    dictionary. 
    """
    user_dict = OrderedDict() 

    for arg in args_dict:
        if arg == 'mode':
            if args_dict[arg] == None:
                user_dict['mode'] ='C' 
            else:
                user_dict['mode'] = args_dict[arg]
        elif arg == 'user_seed':
            if args_dict[arg] == None:
                user_dict['user_seed'] =0 
            else:
                user_dict['user_seed'] = args_dict[arg]
        elif arg == 'SEED':
            if args_dict[arg] == None:
                user_dict['SEED'] = 1 
            else:
                user_dict['SEED'] = args_dict[arg]
        elif arg == 'tpo_fname':
            if args_dict[arg] == None:
                print ('Missing topology file\n')
                print (parser.parse_args(['--help']))
            else:
                user_dict['tpo_fname'] = args_dict[arg]
        elif arg == 'prs_fname':
            if args_dict[arg] == None:
                print ('Missing parameter file\n')
                print (parser.parse_args(['--help']))
            else:
                user_dict['prs_fname'] = args_dict[arg]
        elif arg == 'cfg_fname':
            if args_dict[arg] == None:
                print ('Missing configuration file\n')
                print (parser.parse_args(['--help']))
            else:
                user_dict['cfg_fname'] = args_dict[arg]
        elif arg == 'exp_fname':
            if args_dict[arg] == None:
                print ('Missing expression file\n')
                #print("here 2")
                #sys.exit(0)
                print (parser.parse_args(['--help']))
            else:
                user_dict['exp_fname'] = args_dict[arg]
        elif arg == 'params_input_fname':
            if args_dict[arg] == None:
                print ('Missing parameters file\n')
                print (parser.parse_args(['--help']))
            else:
                user_dict['params_input_fname'] = args_dict[arg]

        elif arg == 'outfile':
            user_dict[arg] = args_dict[arg]
        elif arg == 'knockdownlist': #black list
            if args_dict[arg] == None: 
                print("Missing knock out list") 
                print (parser.parse_args(['--help']))
            else: 
                user_dict['knockdownlist'] = args_dict[arg]

    return user_dict

#-----------------------------------------------------------------------------#    
def parse_args():
    """ 
    This is the entry point for the other methods in this module. It
      1. invokes collect_args to collect the user arguments.
      2. invokes extract_args to put those arguments into an 
         ordered dictionary.
      3. checks the consistency of those arguments by invoking 
         check_args which returns an ordered dictionary of the 
         arguments.
      4. returns the dictionary of the arguments.
    """

    # Collect user arguments:
    parser = collect_args() 
    args_dict = {}
    args, unknown = parser.parse_known_args()
    if len(unknown) > 0:
        print ('\n*********************************')
        print ("Invalid Arguments")
        print ('*********************************\n')
        print (parser.parse_args(['--help']))
    # Places the user arguments into a dictionary:
    args_dict = extract_args(args) 

    # Check the arguments in the dictionary:
    user_dict = check_args(args_dict, parser) 

    return user_dict

if __name__ == '__main__':
    print (sys.argv[0] + ':')
    print (__doc__)
    sys.exit(0)
