import astropy.io.fits as pyfits
import numpy
import json
import matplotlib
import matplotlib.pyplot as plt
import os
import sys


def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Create the fits of rho2 byzones, basically it reads all the .json outputs an calculates the mean rho2')

    # Drectory arguments
    parser.add_argument('--inpath', default='/home2/dfa/sobreira/alsina/catalogs/output/y3a1-v29/byexp/',
                        help='location of work directory')
    parser.add_argument('--outpath', default='./',
                        help='location of the output of the files')
    parser.add_argument('--explist', default='/home/dfa/sobreira/alsina/DESWL/psf/ally3.riz',
                        help='list of run/exposures (in lieu of separate exps, runs)')
 
    args = parser.parse_args()
    return args

def write_fit(data, file_name):
    import fitsio
    from fitsio import FITS,FITSHDR
    import os.path
    if not os.path.isfile(file_name):
        fitsio.write(file_name,  data, clobber=False)
    else:
        fits = FITS(file_name,'rw')
        fits[-1].append(data)

def load_explist(expfile):
    if expfile != '':
        print('Read file ', expfile)
        with open(expfile) as fin:
            exps = [ line.strip() for line in fin if line[0] != '#' ]
        print('File includes %d exposures'%len(exps))
        exps = sorted(exps)
    else:
        print('WARNING: Not exposure list')
    return exps
        
def overall_rho(work, expnum,  outpath):
    import numpy as np
    
    base_keys = ['g', 'r', 'i', 'z']
    keys = [ 'all_' + k for k in base_keys ]

    for key in keys:
        stat_file = os.path.join(work + str(expnum), "rho_" + key + ".json")
        if not os.path.isfile(stat_file):
            print 'File not found: ',stat_file
            continue

        # Read the json file 
        with open(stat_file,'r') as f:
            stats = json.load(f)

        #print' stats = ',stats
        if len(stats) == 1:  # I used to save a list of length 1 that in turn was a list
            stats = stats[0]

        ( meanlogr,
          rho1p,
          rho1p_im,
          rho1m,
          rho1m_im,
          var1,
          rho2p,
          rho2p_im,
          rho2m,
          rho2m_im,
          var2,
          rho3p,
          rho3p_im,
          rho3m,
          rho3m_im,
          var3,
          rho4p,
          rho4p_im,
          rho4m,
          rho4m_im,
          var4,
          rho5p,
          rho5p_im,
          rho5m,
          rho5m_im,
          var5,
        ) = stats[:26]

        #Finally this are the arrays with the data
        meanr = numpy.exp(meanlogr)
        rho1p = numpy.array(rho1p)
        rho1m = numpy.array(rho1m)
        rho2p = numpy.array(rho2p)
        rho2m = numpy.array(rho2m)
        rho3p = numpy.array(rho3p)
        rho3m = numpy.array(rho3m)
        rho4p = numpy.array(rho4p)
        rho4m = numpy.array(rho4m)
        rho5p = numpy.array(rho5p)
        rho5m = numpy.array(rho5m)
        sig_rho1 = numpy.sqrt(var1)
        sig_rho2 = numpy.sqrt(var2)
        sig_rho3 = numpy.sqrt(var3)
        sig_rho4 = numpy.sqrt(var4)
        sig_rho5 = numpy.sqrt(var5)
        sqrtn = 1

        #TAKING THE MEAN by EXP
        names =  ['expnum', 'mrho2p', 'mrho2m',  'msig_rho2']
        formats = ['i4' ,'f8' ,'f8' ,'f8' ]
        values =  [int(expnum), np.nanmean(rho2p), np.nanmean(rho2m),  np.nanmean(sig_rho2)]
        #values =  [zonenum, rho2p[25], rho2m[25],  sig_rho2[25]]    
        print(values)
        dtype = dict(names = names, formats=formats)
        outdata = np.recarray((1, ), dtype=dtype)
        for key, i in zip(names, range(len(names))):
            outdata[key] = values[i]
    
        file_name_out = "y3a1-v29_rho2byexp.fits"
        write_fit(outdata,  file_name_out) 
               

def main():
    import os
    import glob
 
    args = parse_args()

    work = os.path.expanduser(args.inpath)

    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        

    exps = load_explist(args.explist)
    for exp in exps:
        expnum = int(exp)
        overall_rho(work, expnum,  outpath)

 
 

if __name__ == "__main__":
    main()
