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
    args = parser.parse_args()
    return args

def main():
    import fitsio
    from fitsio import FITS,FITSHDR
    import pandas
    args = parse_args()


    #data = fitsio.read(args.file)
    filename =  '/data/git_repositories/basic_tools/y3a1-v29_rho2byexposure2_extended_final.fits'
    data = fitsio.read(filename)
    data = data.astype(data.dtype.newbyteorder('='))
    df = pandas.DataFrame(data)  
 
    columns1 =  ['musestars', 'mtotalstars']
    columns2 =  ['magzp', 'telra', 'teldec', 'telha', 'tiling', 'airmass', 'sat', 'fwhm', 'sky', 'sigsky', 'humidity', 'pressure', 'dimmseeing', 'dT', 'outtemp', 'msurtemp', 'winddir', 'windspd']
    

    dfout =  pandas.DataFrame(index=[0])
    for zn in range(1, 216):
        print(zn)
        zonefile = '/data/git_repositories/basic_tools/riz/zone' + str(zn) + '.riz'
        exps = load_explist(zonefile)
        #exps = [226650, 226651]
        #dfaux =  df.loc[, 'expnum']
        dfaux =  df.loc[df['expnum'].isin(exps)]
        dfaux = dfaux.reset_index(drop=True)
        dfout['zonenum'] =  zn
        for f1 in columns1:
            value = dfaux.loc[(dfaux[f1] != -999.) ,  f1].sum() 
            dfout[f1] = value
        for f2 in columns2:
            value = dfaux.loc[(dfaux[f2] != -999.) ,  f2].mean() 
            dfout[f2] = value

        file_name = "y3a1-v29_rho2byzone_ext.fits"
        write_fit(dfout.to_records(index=False),  file_name)

    
       
    

    
        
        
            

        
        

if __name__ == "__main__":
    main()
