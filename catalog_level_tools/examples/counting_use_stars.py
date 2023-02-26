def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Counting stars in catalog')
    
    parser.add_argument('--explist', default='',
                        help='txt list with the number identifier of the exposure')

    args = parser.parse_args()

    return args

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

def main():
    args = parse_args()
    import fitsio
    import pandas
    import numpy as np

    data =  fitsio.read('stars2.fits')
    data = data.astype(data.dtype.newbyteorder('='))
    df = pandas.DataFrame(data)
    
    exps = load_explist(args.explist)
    
    names =  ['expnum', 'musestars',  'mtotalstars']
    formats = ['i4', 'i4',  'i4' ]
    dtype = dict(names = names, formats=formats)
    outdata = np.recarray((len(exps), ), dtype=dtype)
    nstarslist =  []
    explist =  []
    tstarslist =  []
    
    for exp in sorted(exps):
        expnum =  int(exp)
        boolean1 = (df['exp'] == expnum)     
        tstarslist.append(len(df[boolean1]))
        boolean2 = (df['exp'] == expnum)  & (df['use']== True)
        nstarslist.append( len( df[boolean2] ) )
        explist.append(expnum)

    values =  []
    values.append(explist)
    values.append(nstarslist)
    values.append(tstarslist)
    for key, i in zip(names, range(len(names))):
        outdata[key] = values[i]

    file_name = "stars2_byexp.fits"
    write_fit(outdata,  file_name)

        
if __name__ == "__main__":
    main()
