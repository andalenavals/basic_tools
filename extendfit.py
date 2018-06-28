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

    #Adding columns
    
    rho_data = fitsio.read('y3a1-v29_rho2byexposure_extended.fits')
    rho_data = rho_data.astype(rho_data.dtype.newbyteorder('='))
    rho_df = pandas.DataFrame(rho_data)
    
    env_data = fitsio.read('Y3A1_atmos_pos_condition_conca_sorted.fits')
    env_data = env_data.astype(env_data.dtype.newbyteorder('='))
    env_df = pandas.DataFrame(env_data)

    result = pandas.merge(rho_df, env_df, on='expnum')

    write_fit(result.to_records(index=False), 'y3a1-v29_rho2byexposure_extended_final.fits')
    

    #Adding rows
    '''
    data0 = fitsio.read('Y3A1_atmos_pos_condition0.fits')
    data0 = data0.astype(data0.dtype.newbyteorder('='))
    df0 = pandas.DataFrame(data0)

    data1 = fitsio.read('Y3A1_atmos_pos_condition1.fits')
    data1 = data1.astype(data1.dtype.newbyteorder('='))
    df1 = pandas.DataFrame(data1)

    data2 = fitsio.read('Y3A1_atmos_pos_condition2.fits')
    data2 = data2.astype(data2.dtype.newbyteorder('='))
    df2 = pandas.DataFrame(data2)

    frames = [df0, df1, df2]
    result = pandas.concat(frames, ignore_index=True)
    
    write_fit(result.to_records(index=False), 'Y3A1_atmos_pos_condition_conca.fits')
    '''

if __name__ == "__main__":
    main()
