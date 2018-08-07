import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Read catalog program, it is assume a structure of the catalog. First the expousurelist is a column of numbers. Second each expousure have a folder. Third each folder of each expousure have a exp_info_%d.fits, which have the ccd numbers,')
    
    parser.add_argument('--explist', default='',
                        help='txt list with the number identifier of the exposure')
    parser.add_argument('--explists', nargs='+',  type=str, 
                        help='lists of txt files with the number of exposure to analyse')
    parser.add_argument('--fields', nargs='+',  type=str, 
                        help='list of fields you want to read from the catalog')
    parser.add_argument('--inpath', default='',
                        help='Place where input catalogs is, it is assumed that each expousure have a folder')
    parser.add_argument('--outname', default='', type=str, 
                        help='Name of the output image example zone01.png')

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

def read_somedata(catalogpath,  expolist):
    import fitsio
    import numpy as np
    import pandas
    inpath = os.path.expanduser(catalogpath)
    if not os.path.exists(inpath):
        print('The path of the catalog does not exist!')
        return None

    exps = load_explist(expolist)
    
    names =  ['expnum', 'mean_obs_e1', 'mean_obs_e2',  'mean_obs_e',
    'mean_piff_e1' 'mean_piff_e2',  'mean_piff_e']
    formats = ['i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4' ]
    dtype = dict(names = names, formats=formats)
    outdata = np.recarray((1, ), dtype=dtype)
    
    for exp in sorted(exps):
        expnum = int(exp)
        indir = os.path.join(inpath, exp)
        try:
            expname = os.path.join(indir, 'exp_psf_cat_%d.fits'%expnum)  
            expfile = fitsio.read(expname)
            data = expfile.astype(expfile.dtype.newbyteorder('='))
            outdata['expnum'] = expnum
            outdata['mean_obs_e1'] = np.mean(data['obs_e1'])
            outdata['mean_obs_e2'] = np.mean(data['obs_e2'])
            outdata['mean_obs_e'] = np.mean(data['obs_e1'])
            outdata['mean_piff_e1'] = np.mean(data['obs_e1'])
            outdata['mean_piff_e2'] = np.mean(data['obs_e1'])
            outdata['mean_piff_e'] = np.mean(data['obs_e1'])
             
            print(outdata)
            #write_fit(outdata,  file_name)
        except (OSError, IOError):
            print('Unable to open exp_psf_cat %s.  Skipping this file.'%expinfo)

       
   
   



def write_fit(data, file_name):
    import fitsio
    from fitsio import FITS,FITSHDR
    import os.path
    if not os.path.isfile(file_name):
        fitsio.write(file_name,  data, clobber=False)
    else:
        fits = FITS(file_name,'rw')
        fits[-1].append(data)


def main():
    from astropy.io import fits
    import fitsio
    import numpy as np
    import pandas
    from os.path import basename
    args = parse_args()

    #APP
    read_somedata2(args.inpath,   args.explist)
  

    
if __name__ == "__main__":
    main()
