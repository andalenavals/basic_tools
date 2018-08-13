import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Read catalog program, it is assume a structure of the catalog. First the expousurelist is a column of numbers. Second each expousure have a folder. Third each folder of each expousure have a exp_info_%d.fits, which have the ccd numbers,')
    
    parser.add_argument('--explist', default='',
                        help='txt list with the number identifier of the exposure')
    parser.add_argument('--inpath', default='',
                        help='Place where input catalogs is, it is assumed that each expousure have a folder')

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


def getting_somedata(catalogpath,  expolist):
    import fitsio
    import numpy as np
    import pandas
    inpath = os.path.expanduser(catalogpath)
    if not os.path.exists(inpath):
        print('The path of the catalog does not exist!')
        return None

    exps = load_explist(expolist)
    
    names =  ['expnum', 'mean_obs_e1', 'mean_obs_e2', 'mean_obs_e', 'mean_obs_epw2', 'mean_piff_e1','mean_piff_e2', 'mean_piff_e', 'mean_piff_epw2', 'mean_de1',  'mean_de2',  'mean_de', 'mean_depw2']
    formats = ['i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4','f4', 'f4', 'f4', 'f4' ]
    dtype = dict(names = names, formats=formats)
    outdata = np.recarray((1, ), dtype=dtype)
    outfile_name = 'mean_ellip_byexp.fits'
    for exp in sorted(exps):
        expnum = int(exp)
        indir = os.path.join(inpath, exp)
        try:
            expname = os.path.join(indir, 'exp_psf_cat_%d.fits'%expnum)  
            expfile = fitsio.read(expname)
            data = expfile.astype(expfile.dtype.newbyteorder('='))
            data =  data[data['piff_flag']==0]
           
            outdata['expnum'] = expnum
            outdata['mean_obs_e1'] = np.mean(data['obs_e1'])
            outdata['mean_obs_e2'] = np.mean(data['obs_e2'])
            outdata['mean_obs_e'] = np.mean(np.sqrt(data['obs_e1']**2+data['obs_e2']**2))
            outdata['mean_obs_epw2'] = np.mean( data['obs_e1']**2+data['obs_e2']**2 )
            outdata['mean_piff_e1'] = np.mean(data['piff_e1'])
            outdata['mean_piff_e2'] = np.mean(data['piff_e2'])
            outdata['mean_piff_e'] = np.mean(np.sqrt(data['piff_e1']**2+data['piff_e2']**2))
            outdata['mean_piff_epw2'] = np.mean(data['piff_e1']**2+data['piff_e2']**2 )
            outdata['mean_de1'] = np.mean(data['obs_e1'] -  data['piff_e1'])
            outdata['mean_de2'] = np.mean(data['obs_e2'] -  data['piff_e2'])
            outdata['mean_de'] = np.mean( np.sqrt((data['obs_e1'] -  data['piff_e1'])**2 + (data['obs_e2'] -  data['piff_e2']) ** 2 ) )
            outdata['mean_depw2'] = np.mean( (data['obs_e1'] -  data['piff_e1'])**2 + (data['obs_e2'] -  data['piff_e2']) ** 2 ) 
             
            #print(outdata)
            write_fit(outdata,  outfile_name) 
        except (OSError, IOError):
            print('Unable to open exp_psf_cat %s.  Skipping this file.'%expinfo)


def main():
    args=parse_args()
    getting_somedata(args.inpath,   args.explist)
  

    
if __name__ == "__main__":
    main()
