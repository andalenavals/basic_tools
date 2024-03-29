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


def obfuscate(s):
    mask = 'I Have a Dream'
    lmask = len(mask)
    nmask = [ord(c) for c in mask]
    return ''.join([chr(ord(c) ^ nmask[i % lmask]) for i, c in enumerate(s)])

def ps():
    return obfuscate('~\x10+\t\x1f\x15S\x07T3')
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Wget data based on erinfiles that summarizes the path where there located in the cloud the exposures, some fields are expnum, ccdnum, band, path,  magzp, telra, teldec,  telha , tiling, airmass, sat, fwhm, sky, sigsky, humidity, pressure,  dimmseeing,  dT, outtemp, msurtemp,  winddir,  windspd ')
    parser.add_argument('--outdir', default='wgetdataoutdir',
                        help='location of outputs')
    parser.add_argument('--filexps', default='../ally3.riz',
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--erinfile', default='../exposures-ccds-Y3A1_COADD.fits',
                        help='list of exposures (in lieu of separate exps)')
    args = parser.parse_args()
    return args

def remove_temp_files(wdir, root):
    """Remove wdir/root* except for any files listed in the keep_files list
    """
    files = sorted(glob.glob('%s/%s*'%(wdir,root)))
    #print('   Removing the following files from ',wdir)
    for f in files:
        #print('       ',os.path.split(f)[1])
        os.remove(f)
    #print('   Done')
def read_image_header(row, img_file):
    import galsim
    """Read some information from the image header and write into the df row.
    """
    hdu = 1

    # Note: The next line usually works, but fitsio doesn't support CONTINUE lines, which DES
    #       image headers sometimes include.
    #h = fitsio.read_header(img_file, hdu)
    # I don't  care about any of the lines the sometimes use CONITNUE (e.g. OBSERVER), so I
    # just remove them and make the header with the rest of the entries.
    f = fitsio.FITS(img_file)
    header_list = f[hdu].read_header_list()
    header_list = [ d for d in header_list if 'CONTINUE' not in d['name'] ]
    h = fitsio.FITSHDR(header_list)
    try:
        date = h['DATE-OBS']
        date, time = date.strip().split('T',1)

        filter = h['FILTER']
        filter = filter.split()[0]

        band =  h['BAND']
        sat = h['SATURATE']
        fwhm = h['FWHM']

        ccdnum = int(h['CCDNUM'])
        detpos = h['DETPOS'].strip()

        telra = h['TELRA']
        teldec = h['TELDEC']
        telha = h['HA']
        telra = galsim.Angle.from_hms(telra) / galsim.degrees
        teldec = galsim.Angle.from_dms(teldec) / galsim.degrees
        telha = galsim.Angle.from_hms(telha) / galsim.degrees
        dimmseeing =  float(h.get('DIMMSEE',-999))
        
        airmass = float(h.get('AIRMASS',-999))
        pressure = float(h.get('PRESSURE',-999))
        sky = float(h.get('SKYBRITE',-999))
        sigsky = float(h.get('SKYSIGMA',-999))
        humidity = float(h.get('HUMIDITY',-999))
        windspd = float(h.get('WINDSPD',-999))
        winddir = float(h.get('WINDDIR',-999))

        tiling = int(h.get('TILING',0))
        outtemp = float( h.get('OUTTEMP', - 999) )
        msurtemp = float(h.get('MSURTEMP',-999)) 
        dT = float(h.get('MSURTEMP',-999)) - float(h.get('OUTTEMP',0) )
        hex = int(h.get('HEX',0))


    except:
        print("Cannot read header information from " + img_file)
        raise

    row['sat'] = sat
    row['fits_filter'] = filter
    row['fwhm'] = fwhm
    row['fits_ccdnum'] = ccdnum
    row['airmass'] = airmass
    row['sky'] = sky
    row['sigsky'] = sigsky
    row['humidity'] = humidity
    row['tiling'] = tiling
    row['hex'] = hex
    row['band'] =  band
    row['telra'] = telra
    row['teldec'] = teldec
    row['telha'] = telha
    row['dimmseeing'] =  dimmseeing
    row['dT'] =  dT
    row['outtemp'] = outtemp
    row['msurtemp'] = msurtemp
    row['windspd'] = windspd
    row['winddir'] = winddir 
    
def wget(url_base, path, wdir, file):
    url = url_base + path + file
    full_file = os.path.join(wdir,file)

    import wget
    if not os.path.isfile(full_file):
        #print('Downloading ',full_file)
        # Sometimes this fails with an "http protocol error, bad status line".
        # Maybe from too many requests at once or something.  So we retry up to 5 times.
        nattempts = 5
        for attempt in range(1,nattempts+1):
            print('wget %s  (attempt %d)'%(url, attempt))
            try:
                wget.download(url, bar=None, out=full_file)
            except KeyboardInterrupt:
                raise
            except Exception as e:
                print('Caught ',e)
                if attempt < nattempts:
                    print('Try again.')
                    import time
                    time.sleep(2)
                continue
            else:
                break
    return full_file

def write_fit(data, file_name):
    import fitsio
    from fitsio import FITS,FITSHDR
    import os.path
    if not os.path.isfile(file_name):
        fitsio.write(file_name,  data, clobber=False)
    else:
        try:
            fits = FITS(file_name,'rw')
            fits[-1].append(data)
        except OSError as e:
            print("Ignore OSError from writing:",  data['expnum'])
            print(e)
            pass
 
#Getting data from particular ccd
def GetData1(args):
    wdir = os.path.expanduser(args.outdir)
    print('work dir = ',wdir)
    try:
        if not os.path.exists(wdir):
            os.makedirs(wdir)
    except OSError as e:
        print("Ignore OSError from makedirs(work):")
        print(e)
        pass

    url_base = 'https://rmjarvis:%s@desar2.cosmology.illinois.edu/DESFiles/desarchive/'%ps()
    all_exp = fitsio.read(args.erinfile)
    all_exp = all_exp.astype(all_exp.dtype.newbyteorder('='))
   
    if args.filexps != '':
        print('Read file ',args.filexps)
        with open(args.filexps) as fin:
            exps = [ line.strip() for line in fin if line[0] != '#' ]
        print('File includes %d exposures'%len(exps))

    #exps = sorted(exps)

    #for exp in exps:
    for exp in [229360, 229362]:
        exp = int(exp)
        #print(exp)
        data = all_exp[all_exp['expnum'] == exp]
        #print(data) This line is just because we are interested in
        #the quantities at the level of exposure, so grossly any ccd
        #have that enverimental quantity not changes assumed among ccds
        data =  data[data['ccdnum']==28]
        #save in mem only particular columns of erisn file
        exp_df = pandas.DataFrame(data,  columns=['expnum', 'ccdnum', 'band',  'path',  'magzp'])

        # Add some blank columns to be filled in below.
        
        for k in [ 'telra', 'teldec',  'telha' , 'tiling',  'airmass', 'sat', 'fwhm', 'sky',  'sigsky',  'humidity',  'pressure',  'dimmseeing',  'dT', 'outtemp',  'msurtemp',  'winddir',  'windspd']:
            exp_df[k] = [-999.] * len(data)

        
        try:
            row = pandas.Series(exp_df.iloc[0])
            
        except:
            print("Unxpected error from exp:",  exp)
            print(sys.exc_info()[0])
            continue
        
        path = row['path'].strip()
        #print('path = ',path)
        
        base_path, _, _, image_file_name = path.rsplit('/',3)
        root, ext = image_file_name.rsplit('_',1)
        #print('root, ext = |%s| |%s|'%(root,ext))
        image_file = wget(url_base, base_path + '/red/immask/', wdir, root + '_' + ext)
        #print('image_file = ',image_file)
        
        read_image_header(row, image_file)
        remove_temp_files(wdir,  root)
        exp_df.iloc[0] = row
     
        #file_name = os.path.join(wdir, '%d_Y3A1_atmos_pos_condition.fits'%exp)
        file_name = os.path.join(wdir, 'Y3A1_some_parameters.fits')
        write_fit(exp_df.to_records(index=False), file_name)
        #print('Done with exposure ',exp)       
        
    print('\nFinished processing all exposures')

#Getting all the data
def GetData2(args):
    wdir = os.path.expanduser(args.outdir)
    print('work dir = ',wdir)
    try:
        if not os.path.exists(wdir):
            os.makedirs(wdir)
    except OSError as e:
        print("Ignore OSError from makedirs(work):")
        print(e)
        pass

    url_base = 'https://rmjarvis:%s@desar2.cosmology.illinois.edu/DESFiles/desarchive/'%ps()
    all_exp = fitsio.read(args.erinfile)
    all_exp = all_exp.astype(all_exp.dtype.newbyteorder('='))
   
    if args.filexps != '':
        print('Read file ',args.filexps)
        with open(args.filexps) as fin:
            exps = [ line.strip() for line in fin if line[0] != '#' ]
        print('File includes %d exposures'%len(exps))

    for exp in sorted(exps):
        exp = int(exp)
        #print(exp)
        data = all_exp[all_exp['expnum'] == exp]
        #print(data)
        exp_df = pandas.DataFrame(data,  columns=['expnum', 'ccdnum', 'band',  'path',  'magzp'])

        # Add some blank columns to be filled in below.
        
        for k in [ 'telra', 'teldec',  'telha' , 'tiling',  'airmass', 'sat', 'fwhm', 'sky',  'sigsky',  'humidity',  'pressure',  'dimmseeing',  'dT', 'outtemp',  'msurtemp',  'winddir',  'windspd']:
            exp_df[k] = [-999.] * len(data)



        for k, row in exp_df.iterrows():
            ccdnum =  row['ccdnum']
            try:
                path = row['path'].strip()
                base_path, _, _, image_file_name = path.rsplit('/',3)
                root, ext = image_file_name.rsplit('_',1)
                #print('root, ext = |%s| |%s|'%(root,ext))
                print("Place where file is being looked", url_base + path+ '/red/immask/' + root + '_' + ext)
                image_file = wget(url_base, base_path + '/red/immask/', wdir, root + '_' + ext)
                #print('image_file = ',image_file)
                read_image_header(row, image_file)
                #remove_temp_files(wdir,  root)
                exp_df.iloc[k] = row
                
            except:
                print("Unxpected error from exp, ccd:",  exp,  ccdnum)
                print(sys.exc_info()[0])
                continue
        
        
        #print('path = ',path)
     
        #file_name = os.path.join(wdir, '%d_Y3A1_atmos_pos_condition.fits'%exp)
        file_name = os.path.join(wdir, 'Y3A1_extrafields.fits')
        write_fit(exp_df.to_records(index=False), file_name)
        #print('Done with exposure ',exp)       
        
    print('\nFinished processing all exposures')

def main():
    args = parse_args()
    
    GetData2(args)

if __name__ == "__main__":
    main()
