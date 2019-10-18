#!/usr/bin/python
'''
    This methods in this script check the format of an input file. 
    It has the following methods to check the file format:

    check_gaf_format():
        It checks wheter the format of the file is in GAF. 
        If the file is in GAF format, it returns True
        Otherwise, it returns False

    check_sprot_format(fh_sprot):
        This method checks whether the format of the file
        (with file handle fh_sprot) is in UniProtKB/Swissprot 
        format.
        If the file is in UniProtKB/Swissprot format format,
            it returns True
        Otherwise,
            it returns False.

    check_benchmark_format:
        This method returns False:
            if the input file name is an empty string or
            if the file does not exist or
            if the file size is zero or
            if the file is in correct format
        Otherwise, it returns True
'''
import os
import sys
import re
from os.path import basename

from Bio import SwissProt as sp
import stat

def check_gaf_format(fh_goa):
    """
    This method checks whether the format of the file
    (with file handle fh_goa) is in GAF 1.0 or GAF 2.0.
    If the file is in GAF format, it returns True
    Otherwise, it returns False.
    """
    firstline = fh_goa.readline()
    fields = firstline.strip().split('\t')
    if re.search('^\!gaf', firstline):
        return True
    elif len(fields) == 15:
        return True
    else:
        return False

def check_sprot_format(fh_sprot):
    """
    This method checks whether the format of the file
    (with file handle fh_sprot) is in UniProtKB/Swissprot format.
    If the file is in UniProtKB/Swissprot format format,
        it returns True
    Otherwise,
       it returns False.
    """
    iter_handle = sp.parse(fh_sprot) # sp.parse method returns a generator
    try:
        for rec in iter_handle:
            break
    except:
        return False
    else:
        return True

def check_benchmark_format(fh_benchmark):
    """
    This method checks the format of a benchmark file. 
    It returns False:
        if the the file is NOT in correct 2-column format
    Otherwise, it returns True
    """
    for lines in fh_benchmark:
        cols = lines.strip().split('\t')
        if len(cols) != 2:
            return False
    return True

if __name__ == '__main__':
    print (sys.argv[0] + ':')
    print (__doc__)
