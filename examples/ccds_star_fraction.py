#! /usr/bin/env python
# Compute rho statistics on PSFEx outputs.
# This involves creating catalogs of shapes based on the PSFEx files, and then using
# TreeCorr to compute the correlation functions.

from __future__ import print_function
import os
import sys
import numpy as np
import fitsio
from toFocal import toFocal


def parse_args():
    import argparse
    
    parser = argparse.ArgumentParser(description='This script calculate the mean number of rejected ccds base on its mean ellipticities and a threshold ')

    # Drectory arguments
    parser.add_argument('--inpath', default='/home2/dfa/sobreira/alsina/catalogs/y3a1-v29',
                        help='location of work directory')
    parser.add_argument('--explist', default='/home/dfa/sobreira/alsina/DESWL/psf/ally3.riz',
                        help='list of run/exposures (in lieu of separate exps, runs)')
 
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
def read_data(exps, work, keys, limit_bands=None, single_ccd=None , threshold=0.1 ,   prefix='piff', use_reserved=False, frac=1.):
    RESERVED = 64
    NOT_STAR = 128

    BAD_CCDS = [2, 31, 61]
    MAX_TILING = 10

    all_keys = keys

    if 'x' in keys:
        # Note: use this syntax rather than += to avoid changing input keys
        all_keys = all_keys + ['fov_x', 'fov_y']

    all_keys = all_keys + ['exp', 'ccd', 'band', 'tiling']

    all_data = { key : [] for key in all_keys }

    bands = set()   # This is the set of all bands being used
    tilings = set()   # This is the set of all tilings being used

    n_reject_mean_dt = 0
    n_reject_mean_de1 = 0
    n_reject_mean_de2 = 0
    n_reject_mean_e1 = 0
    n_reject_mean_e2 = 0
    n_reject_std_dt = 0
    n_reject_std_de1 = 0
    n_reject_std_de2 = 0
    n_reject_rho2 = 0

    nrows = 0
    n_good_obj = 0
    n_bad_obj = 0

    n_reject_ccds = 0

    for exp in sorted(exps):

        expnum = int(exp)
        old = False

        expname = os.path.join(work, exp, 'exp_psf_cat_%d.fits'%expnum)
        if not os.path.exists(expname):
            # Old file name (v21 and earlier)
            old = True
            expname = os.path.join(work, exp, 'exp_info_%d.fits'%expnum)
        if not os.path.exists(expname):
            print('Neither kind of exposure catalog exists.')
            print(os.path.join(work, exp, 'exp_psf_cat_%d.fits'%expnum))
            print(expname)
            print('Skip this exposure')
            continue

        try:
            if old:
                expcat = fitsio.read(expname)
            else:
                print(expname)
                expcat = fitsio.read(expname, ext='info')
        except Exception as e:
            print('Error reading exposure catalog')
            print('Caught ',e)
            print('Skip this exposure')
            continue

        if len(expcat) == 0:
            print('expcat for exp=%d has no entries!',expnum)
            continue

        band = expcat['band'][0]
        if (limit_bands is not None) and (band not in limit_bands):
            #print('Not doing band = %s.'%band)
            continue

        if expcat['expnum'][0] != expnum:
            print('%s != %s'%(expcat['expnum'][0], expnum))
            print('expnum in catalog does not match expected value')
            sys.exit()

        print('Start work on exp = ',exp)
        print('band = ',band)

        if 'tiling' in expcat.dtype.names:
            tiling = int(expcat['tiling'][0])
            if tiling == 0:
                # This shouldn't happen, but it did for a few exposures.  Just skip them, since this
                # might indicate some kind of problem.
                print('tiling == 0.  Skip this exposure.')
                continue
            if tiling > MAX_TILING:
                print('tiling is > %d.  Skip this exposure.'%MAX_TILING)
                continue
            print('tiling = ',tiling)
        else:
            tiling = 0

        if old:
            mask = ~np.in1d(expcat['ccdnum'], BAD_CCDS)
            mask &= expcat['flag'] == 0
            data, ccdnums = old_read_ccd_data(expcat[mask], work, expnum)
        else:
            data = fitsio.read(expname, ext='stars')
            ccdnums = data['ccdnum'].astype(int)

        if single_ccd:
            print('Enter in single_ccd mode')
            data = data[data['ccdnum']==single_ccd]

            
        flag = data[prefix+'_flag'].astype(int)
        ntot = len(data)
        nused = np.sum((flag & 1) != 0)
        nreserved = np.sum((flag & RESERVED) != 0)
        ngood = np.sum(flag == 0)
        print('ntot = ',ntot)
        print('nused = ',nused)
        print('nreserved = ',nreserved)
        print('ngood = ',ngood)

        mask = (flag & NOT_STAR) == 0
        mask &= ~np.in1d(ccdnums, BAD_CCDS)
        

        #SKIPING CANNON REGION.
        #IsInCannon = ( (data['ra'] <  0) & (data['dec'] > -10 ) )
        #mask &= ~IsInCannon

        
        if use_reserved:
            mask &= (flag & RESERVED) != 0
        used = (flag & ~RESERVED) == 0
        #print('flag = ',flag)
        #print('mask = ',mask)
        #print('used = ',used)
        #print('flag where flag == RESERVED: ',flag[flag==RESERVED+1])
        #print('mask where flag == RESERVED: ',mask[flag==RESERVED+1])
        #print('used where flag == RESERVED: ',mask[flag==RESERVED+1])
        print('nmask = ',np.sum(mask))
        print('nused = ',np.sum(used))
        
            
        T = data['obs_T']
        e1 = data['obs_e1']
        e2 = data['obs_e2']
        dT = data['obs_T'] - data[prefix + '_T']
        de1 = data['obs_e1'] - data[prefix + '_e1']
        de2 = data['obs_e2'] - data[prefix + '_e2']
        print(expnum, len(dT), band)
        #print('T = ',np.mean(T[used]),np.std(T[used]))
        #print('e1 = ',np.mean(e1[used]),np.std(e1[used]))
        #print('e2 = ',np.mean(e2[used]),np.std(e2[used]))
        #print('dT/T = ',np.mean(dT[used]/T[used]),np.std(dT[used]/T[used]))
        #print('de1 = ',np.mean(de1[used]),np.std(de1[used]))
        #print('de2 = ',np.mean(de2[used]),np.std(de2[used]))
        rho1 = (de1 - 1j*de2) * (de1 + 1j*de2)
        #print('mean rho1 = ',np.mean(rho1[used]))
        rho2 = (e1 - 1j*e2) * (de1 + 1j*de2)
        #print('mean rho2 = ',np.mean(rho2[used]))

        #Neglecting ccds with high mean ellipticity mod andres
        neglected_ccds = []
        if threshold:
            for cnum in range(1, 63):
                dataux =  data[data['ccdnum']==int(cnum)]
                eccdobs =  np.nanmean(np.sqrt(dataux['obs_e1']**2 + dataux['obs_e2']**2) )
                eccdpiff =  np.nanmean(np.sqrt(dataux[prefix +'_e1']**2 + dataux[prefix + '_e2']**2) )
                if ((eccdpiff >   threshold)|(eccdobs>threshold)):             
                    neglected_ccds.append(int(cnum))
                    #print('Blacklisting stars in ccdnum',  cnum,  'from exposure ',  expnum )
        print("List of neglected ccds", neglected_ccds)
        nccds = set(neglected_ccds).union(set(BAD_CCDS))
        n_reject_ccds += len(nccds)
        

        # Filter out egregiously bad values.  Just in case.
        n1 = np.sum(mask)
        good = (abs(dT/T) < 0.1) & (abs(de1) < 0.1) & (abs(de2) < 0.1)
        mask &= ~np.in1d(ccdnums, neglected_ccds) ##suspicius ccds modify andres
        mask = mask & good
        n2 = np.sum(mask)
        print('"good" filter removed %d/%d objects'%(n1-n2,n1))
        n_good_obj += n2
        n_bad_obj += n1-n2
        #print('mask = ',len(mask),np.sum(mask),mask)
        mask = np.where(mask)[0]
        #print('mask = ',len(mask),mask)
        if frac != 1.:
            mask = np.random.choice(mask, int(frac * len(mask)), replace=False)
        #print('mask = ',len(mask),mask)
        ngood = len(mask)
        print('ngood = ',ngood,'/',len(data))
        assert ngood == len(data[mask])
        if ngood == 0:
            print('All objects in exp %d are flagged.'%expnum)
            continue

        # Start with just the input keys, which should be columns in data.
        for key in keys:
            all_data[key].append(data[key][mask])

        # Now add the extra ones we added to all_keys
        if 'x' in keys:
            # Convert to focal position.
            x,y = toFocal(ccdnums[mask], data['x'][mask], data['y'][mask])
            # This comes back in units of mm.  Convert to arcsec.
            # 1 pixel = 15e-3 mm = 0.263 arcsec
            x *= 0.263/15e-3
            y *= 0.263/15e-3
            all_data['fov_x'].append(x)
            all_data['fov_y'].append(y)

        all_data['ccd'].append(ccdnums[mask])
        all_data['exp'].append([expnum] * ngood)
        all_data['band'].append([band] * ngood)
        all_data['tiling'].append([tiling] * ngood)
        bands.add(band)
        tilings.add(tiling)
        nrows += ngood

    print('\nFinished processing %d exposures'%len(exps))
    print('bands = ',bands)
    print('tilings = ',tilings)
    print('total "good" stars = ',n_good_obj)
    print('total "bad" stars = ',n_bad_obj)
    print('total good stars selected = ',nrows)
    print ('fraction of select good objects /total objects =', float(n_good_obj) / float(n_good_obj + n_bad_obj) )
    print ('total rejected ccds', n_reject_ccds )
    print ('fraction of rejected ccds', float(n_reject_ccds) / float(62 * len(exps) ))
    # But this rejection would be at the level of exposures
    print('Potential rejections: (not enabled):')
    print('n_reject_mean_dt = ',n_reject_mean_dt)
    print('n_reject_mean_de1 = ',n_reject_mean_de1)
    print('n_reject_mean_de2 = ',n_reject_mean_de2)
    print('n_reject_mean_e1 = ',n_reject_mean_de1)
    print('n_reject_mean_e2 = ',n_reject_mean_de2)
    print('n_reject_std_dt = ',n_reject_std_dt)
    print('n_reject_std_de1 = ',n_reject_std_de1)
    print('n_reject_std_de2 = ',n_reject_std_de2)
    print('n_reject_rho2 = ',n_reject_rho2)

 

def main():
 
    args = parse_args()

    exps = load_explist(args.explist)
    work =  args.inpath
    keys = ['ra', 'dec', 'x', 'y', 'obs_e1', 'obs_e2', 'obs_T', 'piff_e1','piff_e2', 'piff_T']
    read_data(exps, work, keys, limit_bands='grizY', single_ccd=None, threshold=0.1 ,  prefix='piff' , use_reserved=True, frac=1.)
 

if __name__ == "__main__":
    main()
