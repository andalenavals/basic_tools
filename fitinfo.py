def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple program to read .fits files')
    
    parser.add_argument('--file', default='',
                        help='Fits image to read')
    parser.add_argument('--hdu', default='',
                        help='int number of the hdu of interes')
    parser.add_argument('--header', default=False, action='store_const', const=True,
                        help='boolean to print the header of the hdu')
    parser.add_argument('--columns', default=False, action='store_const', const=True,
                        help='boolean to print the columns of the hdu')
    parser.add_argument('--alldata', default=False, action='store_const', const=True,
                        help='boolean to print all the data of the image or tabla')
    parser.add_argument('--fields', nargs='+',  type=str, 
                        help='Print only some columns, --fields string1 string2 string3 ... ')
    parser.add_argument('--rows', nargs='+',  type=int, 
                        help='Print only some rows, --rows 0 50 (it is equivalente to the numpy operation [0:50]) ... ')
    args = parser.parse_args()

    return args

def PrettyPrint(data):
    from prettytable import PrettyTable
    x = PrettyTable(data.dtype.names)    
    for row in data:
        x.add_row(row)
        # Change some column alignments; default was 'c'
        #x.align['column_one'] = 'r'
        #x.align['column_3'] = 'l'
    print(x)
  
def fields_view(arr, fields):
    import numpy as np
    dtype2 = np.dtype({name:arr.dtype.fields[name] for name in fields})
    return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)
def main():
    import fitsio
    from astropy.io import fits
    import pandas
    
    args = parse_args()

    try:
        hdu =  fits.open(args.file)
    except:
        print('File could not be open!')
        return None
    if not (args.fields or args.rows or args.alldata):
        hdu.info()
    if (args.hdu):
        nhdu = int(args.hdu)
    if (args.hdu and args.header):
        print ('\n Reading HDU ',nhdu,' header\n',  repr(hdu[nhdu].header) )
    if (args.hdu and args.columns):       
        try:
            hdu[nhdu].columns
            print('Reading columns information of the HDU', nhdu,'\n', repr(hdu[nhdu].columns))
        except:
            print ('The',  type(hdu[nhdu]),'of the HDU',nhdu, 'does not have columns, probably it is a image instead a table')
    if (args.hdu and args.alldata and not args.fields ):
        try:
             data =  hdu[nhdu].data
        except:
            print ('The',  type(hdu[nhdu]),'of the HDU',nhdu,   'is not readable')
        if not (args.rows):
            PrettyPrint(data)
        else:
            try:
                aux = data[args.rows[0]:args.rows[1]]
                PrettyPrint(aux)
            except:
                print('Error, problem reading some of the rows, probably data[0] does not exist' )
                return None     
    if (args.hdu and args.fields ):
        try:
             data =  hdu[nhdu].data
             aux1 = fields_view(data, args.fields)
        except:
            print ('The', type(hdu[nhdu]),'of the HDU',nhdu, 'is not a BinTable. Or the header does not exist. Or the field does not exist')
            return None
        if not(args.rows):
            PrettyPrint(aux1)
        else:
            try:
                aux2 = aux1[args.rows[0]:args.rows[1]]
                PrettyPrint(aux2)
            except:
                print('Error, Rows were not properly selected' )
                return None
    
if __name__ == "__main__":
    main()
    

 
