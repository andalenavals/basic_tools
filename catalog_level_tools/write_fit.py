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
    import numpy as np
    names =  ['expnum', 'meanr',  'rho2p', 'rho2m',  'sig_rho2']
    formats = ['i4', 'f8' ,'f8' ,'f8' ,'f8' ]
    values =  [1, 2, 3, 4, 5]    
    dtype = dict(names = names, formats=formats)
    outdata = np.recarray((1, ), dtype=dtype)
    for key, i in zip(names, range(len(names))):
        outdata[key] = values[i]
        
    file_name = "y3a1-v29_rho2byexposure.fit"
    write_fit(outdata,  file_name)

        
if __name__ == "__main__":
    main()
