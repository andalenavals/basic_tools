def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Read catalog program format h5')
    
    parser.add_argument('--psf_cat', default='/data/catalogs/y3_master/Y3_mastercat_v1_6_20_18_subsampled.h5',
                        help='Full Path to the catalog')
 

    args = parser.parse_args()

    return args


def main():
    from astropy.io import fits
    import fitsio
    import numpy as np
    import pandas as pd
    import h5py
    
    args = parse_args()

    h5 =  pd.HDFStore(args.psf_cat, 'r')
    print(h5)
    '''
    h5 = h5py.File(args.psf_cat, 'r')
    print(h5.keys())
    cat = h5['catalog']
    print(cat.shape())
    print(cat.dtype())
    '''
    #psf_data=fitsio.read(args.psf_cat)

if __name__ == "__main__":
    main()
