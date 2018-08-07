import astropy.io.fits as pyfits
import numpy
import json
import matplotlib
import matplotlib.pyplot as plt
import os
import sys


def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='Run PSFEx on a set of runs/exposures')

    # Drectory arguments
    parser.add_argument('--work', default='./',
                        help='location of work directory')
    parser.add_argument('--outpath', default='./',
                        help='location of the output of the files')
    parser.add_argument('--tag', default=None,
                        help='A version tag to add to the directory name')

    # Exposure inputs
    parser.add_argument('--exp_match', default='*_[0-9][0-9].fits.fz',
                        help='regexp to search for files in exp_dir')
    parser.add_argument('--file', default='',
                        help='list of run/exposures (in lieu of separate exps, runs)')
    parser.add_argument('--exps', default='', nargs='+',
                        help='list of exposures to run')
    parser.add_argument('--runs', default='', nargs='+',
                        help='list of runs')
    parser.add_argument('--expnum', default=None, type=int, 
                        help='')
    

    args = parser.parse_args()
    return args


def plot_overall_rho(work,  outpath, args):
    import numpy as np
    
    base_keys = ['riz']
    keys = [ 'all_' + k for k in base_keys ]

    for key in keys:
        stat_file = os.path.join(work, "rho_" + key + ".json")
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

        #ORIGINAL
        '''
        names =  ['expnum', 'meanr',  'rho2p', 'rho2m',  'sig_rho2']
        formats = ['i4', 'f8' ,'f8' ,'f8' ,'f8' ]
        values =  [int(args.expnum), meanr[25], rho2p[25], rho2m[25],  sig_rho2[25]]    
        '''

        #TAKING THE MEAN by exposure
        '''
        names =  ['expnum', 'mrho2p', 'mrho2m',  'msig_rho2']
        formats = ['i4' ,'f8' ,'f8' ,'f8' ]
        values =  [int(args.expnum), np.mean(rho2p), np.mean(rho2m),  np.mean(sig_rho2)]    
        dtype = dict(names = names, formats=formats)
        outdata = np.recarray((1, ), dtype=dtype)
        for key, i in zip(names, range(len(names))):
            outdata[key] = values[i]
    
        file_name = "y3a1-v29_rho2byexposure2.fit"
        write_fit(outdata,  file_name) 
        '''
        #TAKING THE MEAN by zone
        '''
        name = args.file
        znum = int(name[-name[::-1].index('/')+4:-4])
        names =  ['zonenum', 'mrho2p', 'mrho2m',  'msig_rho2']
        formats = ['i4' ,'f8' ,'f8' ,'f8' ]
        #values =  [znum, np.nanmean(rho2p), np.nanmean(rho2m),  np.nanmean(sig_rho2)]
        values =  [znum, rho2p[25], rho2m[25],  sig_rho2[25]]    
        print(values)
        dtype = dict(names = names, formats=formats)
        outdata = np.recarray((1, ), dtype=dtype)
        for key, i in zip(names, range(len(names))):
            outdata[key] = values[i]
    
        file_name = "y3a1-v29_rho2byzone.fits"
        write_fit(outdata,  file_name) 
        '''

        #TAKING THE MEAN by ccd
    
        ccdnum = int(args.ccdnum)
        names =  ['ccdnum', 'mrho2p', 'mrho2m',  'msig_rho2']
        formats = ['i4' ,'f8' ,'f8' ,'f8' ]
        values =  [ccdnum, np.nanmean(rho2p), np.nanmean(rho2m),  np.nanmean(sig_rho2)]
        #values =  [ccdnum, rho2p[25], rho2m[25],  sig_rho2[25]]    
        print(values)
        dtype = dict(names = names, formats=formats)
        outdata = np.recarray((1, ), dtype=dtype)
        for key, i in zip(names, range(len(names))):
            outdata[key] = values[i]
    
        file_name = "y3a1-v29_rho2byccd.fits"
        write_fit(outdata,  file_name) 
        
       
        

def main():
    import os
    import glob
 
    args = parse_args()

    work = os.path.expanduser(args.work)
    #print 'work dir = ',work

    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    #plot_single_rho(args,work)

    #for zn in range(1, 216):
    #    plot_overall_rho(work,  outpath,  args)

    for ccd in range(1, 216):
    #    plot_overall_rho(work,  outpath,  args)
 

if __name__ == "__main__":
    main()
