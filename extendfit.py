def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple example program to merge or concat fits files')
    
    parser.add_argument('--file', default='',
                        help='Fits image to extend, by default it looks for the first hdu')
    parser.add_argument('--ref', default='',
                        help='Fits with the additional information')
    parser.add_argument('--merge', default=False, action='store_const', const=True,
                        help='boolean to print the header of the hdu')
    parser.add_argument('--concact', default=False, action='store_const', const=True,
                        help='boolean to print the columns of the hdu')
    parser.add_argument('--fields', nargs='+',  type=str, 
                        help='Print only some columns, --fields string1 string2 string3 ... ')
 
    args = parser.parse_args()

    return args


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
    import fitsio
    import pandas

    args = parse_args()
    #Adding columns

    try:
        data = fitsio.read(args.file)
        data = data.astype(data.dtype.newbyteorder('='))
        df = pandas.DataFrame(data)
    except:
        print('File could not be open!')
        return None

    try:
        ref_data = fitsio.read(args.ref)
        ref_data = ref_data.astype(ref_data.dtype.newbyteorder('='))
        ref_df = pandas.DataFrame(ref_data)
    except:
        print('Reference File could not be open!')
        return None
    
    

    result = pandas.merge(df, ref_df, on='zonenum')


    write_fit(result.to_records(index=False), 'y3a1-v29_rho2byzone_extended_final.fits')

    

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
