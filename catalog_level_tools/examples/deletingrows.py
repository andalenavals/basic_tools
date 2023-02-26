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
    import fitsio
    from fitsio import FITS,FITSHDR
    import pandas

    #data = fitsio.read(args.file)
    filename =  '/data/git_repositories/basic_tools/Y3A1_extrafields.fits'
    data = fitsio.read(filename)
    data = data.astype(data.dtype.newbyteorder('='))
    df = pandas.DataFrame(data)  

    dfout =  df.loc[~df['expnum'].isin([233444])]
    
    
    file_name = "Y3A1_extrafields2.fits"
    write_fit(dfout.to_records(index=False),  file_name)

    
       
    

    
        
        
            

        
        

if __name__ == "__main__":
    main()
