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

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Check if information in a fit file is compatible with what we expect')
    parser.add_argument('--file', default='',
                        help='Fit file to check')
    parser.add_argument('--explist', default='/data/git_repositories/basic_tools/all_zones.riz',
                        help='txt list with the number identifier of the exposure')
    args = parser.parse_args()
    return args

def main():
    import fitsio
    from fitsio import FITS,FITSHDR
    import pandas
    args = parse_args()


    #data = fitsio.read(args.file)
    filename =  '/data/git_repositories/basic_tools/Y3A1_extrafields.fits'
    data = fitsio.read(filename)
    data = data.astype(data.dtype.newbyteorder('='))
    df = pandas.DataFrame(data)  
 
    BAD_CCDS = [2, 31, 61]
    MAX_TILING = 10
    boolean =  (df['ccdnum'].isin(BAD_CCDS)) | (df['tiling'] >  MAX_TILING) | (df['tiling'] == 0)
    dfcl =  df.loc[ ~boolean ]
    dfcl = dfcl.reset_index(drop=True)

    dfout =  pandas.DataFrame(index=[0])
    columns =  ['magzp', 'telra', 'teldec', 'telha', 'tiling', 'airmass', 'sat', 'fwhm', 'sky', 'sigsky', 'humidity', 'pressure', 'dimmseeing', 'dT', 'outtemp', 'msurtemp', 'winddir', 'windspd']
    exps = sorted(load_explist(args.explist))
    for exp in exps:
    #for exp in [226650]:
        expnum =  int(exp)
        dfaux = dfcl.loc[(dfcl['expnum']==expnum)]
        dfaux = dfaux.reset_index(drop=True)
        print(expnum)
        dfout['expnum'] = expnum
        dfout['band'] = dfaux['band'][0]
        for col in columns:
            value = dfaux.loc[(dfaux[col] != -999.) ,  col].mean() 
            dfout[col] = value
        #print(dfout)
        file_name = "y3a1-v29_rho2byexposure2_env.fits"
        write_fit(dfout.to_records(index=False),  file_name)
        
    

       
    

    
        
        
            

        
        

if __name__ == "__main__":
    main()
