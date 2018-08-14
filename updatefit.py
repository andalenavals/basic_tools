def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple example program to merge or concat fits files')
    
    parser.add_argument('--file', default='',
                        help='Fits image to extend, by default it looks for the first hdu')
    parser.add_argument('--ref', default='',
                        help='Fits with the additional information')
    parser.add_argument('--out_filename', default='outupdate.fits',
                        help='Filename for the output ')
 
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
    #Updating files
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
    
    df.update( ref_df)
    write_fit(df.to_records(index=False), args.out_filename)

   

if __name__ == "__main__":
    main()
