import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Filter exposure with certain property,')
    
    parser.add_argument('--explist', default='',
                        help='txt list with the number identifier of the exposure')
    parser.add_argument('--file', default='',
                        help='File where is the information for the filter')

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


def filter_data(filefilter,   expolist):
    import fitsio
    import numpy as np
    import pandas

    exps = load_explist(expolist)
    

    for exp in sorted(exps):
        expnum = int(exp)
        try:
            expfile = fitsio.read(filefilter)
            data = expfile.astype(expfile.dtype.newbyteorder('='))
            flag =  (data['mean_obs_e']<0.06)&(data['expnum'] == expnum)
            data =  data[flag]
            if len(data) > 0:
                print(data['expnum'][0])

            #print(outdata)

        except (OSError, IOError):
            print('Unable to open exp_psf_cat %s.  Skipping this file.'%expinfo)


def main():
    args=parse_args()
    filter_data(args.file,   args.explist)
  

    
if __name__ == "__main__":
    main()
