#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

from __future__ import print_function
import os
import numpy as np
from read_psf_cats_andres import read_data
#from read_psf_cats import read_data

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description='Run rho stats on a set of exposures')

    # Drectory arguments
    parser.add_argument('--inpath', default='/astro/u/mjarvis/work/y3_piff',
                        help='location of work directory')
    parser.add_argument('--explist', default='',
                        help='list of exposures (in lieu of separate exps)')

    args = parser.parse_args()
    return args

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
def write_fit(data, file_name):
    import fitsio
    from fitsio import FITS,FITSHDR
    import os.path
    if not os.path.isfile(file_name):
        fitsio.write(file_name,  data, clobber=False)
    else:
        fits = FITS(file_name,'rw')
        fits[-1].append(data)


def getting_data(catalogpath,  expolist,  fields,  ccdnum):
    import fitsio
    import numpy as np

    inpath = os.path.expanduser(catalogpath)
    if not os.path.exists(inpath):
        print('The path of the catalog does not exist!')
        return None
    
    exps = load_explist(expolist) 

  
    all_data = { key : [] for key in fields }
    #Getting all fields for each exposure and for a particular ccd
    print('Getting data for ccdnum =',  ccdnum)
    for exp in exps:
        #print('Start work on exp = ',exp)
        expnum = int(exp)
        #print('expnum = ',expnum)
        indir = os.path.join(inpath, exp)
        try:
            expinfo = fitsio.read(os.path.join(indir, 'exp_psf_cat_%d.fits'%expnum))
            data = expinfo.astype(expinfo.dtype.newbyteorder('='))
            flag =  (data['piff_flag']==0) & (data['ccdnum']==ccdnum)
            data =  data[flag]
            for key in fields:
                all_data[key].append(data[key])
        except (OSError, IOError):
            print('Unable to open exp_psf_cat %s.  Skipping this file.'%expinfo)

    all_data_final = { key : [] for key in fields }         
    for key in fields:
        all_data_final[key] = np.concatenate(all_data[key])

    #print(len(all_data_final['obs_e1']))
    #print(all_data_final['obs_e1'])
        
#Adding line by line
    names =  ['ccdnum', 'mean_obs_e1', 'mean_obs_e2', 'mean_obs_e', 'mean_obs_epw2', 'mean_piff_e1','mean_piff_e2', 'mean_piff_e', 'mean_piff_epw2', 'mean_de1',  'mean_de2']
    formats = ['i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4','f4', 'f4' ]
    dtype = dict(names = names, formats=formats)
    outdata = np.recarray((1, ), dtype=dtype)
    outfile_name = 'mean_ellip_byccd.fits'
    outdata['ccdnum'] = ccdnum
    #print(all_data_final)
    
    outdata['mean_obs_e1'] = np.mean(all_data_final['obs_e1'])
    outdata['mean_obs_e2'] = np.mean(all_data_final['obs_e2'])
    outdata['mean_obs_e'] = np.mean(np.sqrt(all_data_final['obs_e1']**2+all_data_final['obs_e2']**2))
    outdata['mean_obs_epw2'] = np.mean( all_data_final['obs_e1']**2+all_data_final['obs_e2']**2 )
    outdata['mean_piff_e1'] = np.mean(all_data_final['piff_e1'])
    outdata['mean_piff_e2'] = np.mean(all_data_final['piff_e2'])
    outdata['mean_piff_e'] = np.mean(np.sqrt(all_data_final['piff_e1']**2+all_data_final['piff_e2']**2))
    outdata['mean_piff_epw2'] = np.mean(all_data_final['piff_e1']**2+all_data_final['piff_e2']**2 )
    outdata['mean_de1'] = np.mean(all_data_final['obs_e1'] -  all_data_final['piff_e1'])
    outdata['mean_de2'] = np.mean(all_data_final['obs_e2'] -  all_data_final['piff_e2'])
    
    #print(outdata)
    write_fit(outdata,  outfile_name) 


def main():

    args = parse_args()
    fields = [ 'obs_e1',  'obs_e2',  'piff_e1',  'piff_e2']
    for k in range(39, 63):
        getting_data(args.inpath, args.explist,  fields,  int(k))

    '''
    names =  ['ccdnum', 'mean_obs_e1', 'mean_obs_e2',  'mean_obs_e', 'mean_piff_e1' ,'mean_piff_e2',  'mean_piff_e']
    formats = ['i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4' ]
    dtype = dict(names = names, formats=formats)
    outdata = np.recarray((1, ), dtype=dtype)
    file_name = "mean_ellip_byccd.fits"
    
    keys = [ 'obs_e1', 'obs_e2', 'piff_e1', 'piff_e2']
    for k in range(1,2):
        
        try:
            ccdnum = int(k)
            data, bands, tilings = read_data(exps, work, keys, limit_bands='riz', single_ccd=ccdnum,  prefix='piff',
                                         use_reserved=False, frac=1.0)
            print(type(data))
            outdata['ccdnum'] = ccdnum
            outdata['mean_obs_e1'] = np.mean(data['obs_e1'])
            outdata['mean_obs_e2'] = np.mean(data['obs_e2'])
            outdata['mean_obs_e'] = np.mean(np.sqrt(data['obs_e1']**2+data['obs_e2']**2))
            outdata['mean_piff_e1'] = np.mean(data['piff_e1'])
            outdata['mean_piff_e2'] = np.mean(data['piff_e2'])
            outdata['mean_piff_e'] = np.mean(np.sqrt(data['piff_e1']**2+data['piff_e2']**2))
        
            print(outdata)                                                                                                                                       
            #write_fit(outdata,  file_name)
        except (OSError, IOError):
            print('Unable to open exp_psf_cat %s.  Skipping this file.'%expinfo)
    '''


if __name__ == "__main__":
    main()
