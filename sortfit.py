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
    import pandas
    if not os.path.isfile(file_name):
        fitsio.write(file_name,  data, clobber=False)
    else:
        fits = FITS(file_name,'rw')
        fits[-1].append(data)

def main():
    
    rho_data = fitsio.read('Y3A1_atmos_pos_condition_conca.fits')
    rho_data = rho_data.astype(rho_data.dtype.newbyteorder('='))
    rho_df = pandas.DataFrame(rho_data)    
    
    newdf =  rho_df.sort_values('expnum')
    newdf = newdf.reset_index(drop=True)

    write_fit(newdf.to_records(index=False), 'Y3A1_atmos_pos_condition_conca_sorted.fits')

if __name__ == "__main__":
    main()
