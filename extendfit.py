def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple example program to merge or concat fits files. Or create new hdu from different files')
    
    parser.add_argument('--file', default='',
                        help='Fits image to extend, by default it looks for the first hdu. Reciver')
    parser.add_argument('--ref', default='',
                        help='Fits with the additional information. Donor')
    parser.add_argument('--hdu', type=int, 
                        help='hdu of the donor ')
    parser.add_argument('--hduconcat', default=False,
                        action='store_const', const=True,
                        help='boolean to active the concatenation of hdu mode')
    parser.add_argument('--merge', default=False, action='store_const', const=True,
                        help='boolean to merge files')
    parser.add_argument('--on', default='expnum', type=str, 
                        help='field for the camparinson and merging')
    parser.add_argument('--concat', default=False, action='store_const', const=True,
                        help='boolean to concact the files')
    parser.add_argument('--fields', nargs='+',  type=str, 
                        help='Print only some columns, --fields string1 string2 string3 ... ')
    parser.add_argument('--out_filename', default='outmc.fits',
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
        #appending data to the las hdu
        fits[-1].append(data)

def concathdu(receptorfitname, donorfitname,  hdunum):
    import fitsio
    from fitsio import FITS,FITSHDR
    from astropy.io import fits
    try:
        donorfit =  fits.open(donorfitname)
    except:
        print('File', donorfitname, 'could not be open!')
        return None
    try:
        receptorfit =  fits.open(receptorfitname)
    except:
        print('File', receptorfitname, 'could not be open!')
        return None
    print(receptorfit)
    receptorfit.insert(1, donorfit[hdunum])
    #receptorfit.append(donorfit[hdunum])
    receptorfit.writeto(receptorfitname, clobber=True)
 
def main():
    import fitsio
    import pandas

    args = parse_args()
    #Adding columns

    if (args.merge):
        try:
            data = fitsio.read(args.file)
            data = data.astype(data.dtype.newbyteorder('='))
            df = pandas.DataFrame(data)
        except:
            print('File could not be open!')
            return None
        try:
            if(args.fields):
                args.fields.append(args.on)
                ref_data = fitsio.read(args.ref,  columns=args.fields)
            else:
                ref_data = fitsio.read(args.ref)
            ref_data = ref_data.astype(ref_data.dtype.newbyteorder('='))
            ref_df = pandas.DataFrame(ref_data)
        except:
            print('Reference File could not be open!')
            return None
    
        result = pandas.merge(df, ref_df, on=args.on)
        write_fit(result.to_records(index=False), args.out_filename)

    elif (args.concat):

        data0 = fitsio.read(args.file)
        data0 = data0.astype(data0.dtype.newbyteorder('='))
        df0 = pandas.DataFrame(data0)
        
        data1 = fitsio.read(args.ref)
        data1 = data1.astype(data1.dtype.newbyteorder('='))
        df1 = pandas.DataFrame(data1)
        
        frames = [df0, df1]
        result = pandas.concat(frames, ignore_index=True)
        
        write_fit(result.to_records(index=False),  args.out_filename)

    else:
        print("Extension unespecified!!: Neither --merge, nor --concat flags were defined. ")

    if(args.hduconcat):
        concathdu(args.file, args.ref,  args.hdu)
    else:
        print('Mode of concatenation or mixturing of hdu unactivated.')

if __name__ == "__main__":
    main()
