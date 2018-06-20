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
    
    rho_data = fitsio.read('y3a1-v29_rho2byexposure.fit')
    rho_data = rho_data.astype(rho_data.dtype.newbyteorder('='))
    rho_df = pandas.DataFrame(rho_data)
    
    env_data = fitsio.read('Y3A1_atmos_pos_condition.fit')
    env_data = env_data.astype(env_data.dtype.newbyteorder('='))
    env_df = pandas.DataFrame(env_data)

    result = pandas.merge(rho_df, env_df, on='expnum')

    write_fit(result.to_records(index=False), 'y3a1-v29_rho2byexposure_extended.fit')

if __name__ == "__main__":
    main()
