from __future__ import print_function
import os
import sys
import traceback
import numpy as np
import copy
import glob
import time
import fitsio
import pandas

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Check if information in a fit file is compatible with what we expect')
    parser.add_argument('--file', default='',
                        help='Fit file to check')
    parser.add_argument('--filexps', default='all_zones.riz',
                        help='list of exposures (in lieu of separate exps)')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    if args.filexps != '':
        print('Read file ',args.filexps)
        with open(args.filexps) as fin:
            exps = [ line.strip() for line in fin if line[0] != '#' ]
        #print('File includes %d exposures'%len(exps))
    exps = sorted(exps)
    
    data = fitsio.read(args.file)
    data = data.astype(data.dtype.newbyteorder('='))
    df = pandas.DataFrame(data)    
    
    for exp in exps:
        exp = int(exp)
       
        #print(data[[data['expnum']==exp]]['msurtemp'][0] )
       
        if(df[df['expnum'] == exp]['expnum'].count() == 0):
            print('Exposure',  exp ,  'Is mising')
        elif (df[df['expnum'] == exp]['expnum'].count() > 1):
            print('Exposure',  exp ,  'Is repeated')
        

        
        
            

        
        

if __name__ == "__main__":
    main()
