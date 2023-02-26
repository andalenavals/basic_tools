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


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Wget a file ')
    parser.add_argument('--url_file', default='https://rmjarvis:70chipsftw@desar2.cosmology.illinois.edu/DESFiles/desarchive/OPS/finalcut/Y2A1/Y1-2355/20130815/D00226650/p02/red/immask/D00226650_z_c38_r2355p02_immasked.fits.fz',
                        help='URL of the file to download')
    parser.add_argument('--outdir', default='',
                        help='Local dir for the outfile')
    args = parser.parse_args()
    return args
 
def wget(url_file, out_file):
    import wget
    if not os.path.isfile(out_file):
        #print('Downloading ',full_file)
        # Sometimes this fails with an "http protocol error, bad status line".
        # Maybe from too many requests at once or something.  So we retry up to 5 times.
        nattempts = 5
        for attempt in range(1,nattempts+1):
            print('wget %s  (attempt %d)'%(url_file, attempt))
            try:
                wget.download(url_file, bar=None, out=out_file)
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
    return out_file


def main():
    args = parse_args()
    wdir = os.path.expanduser(args.outdir)
    print('work dir = ',wdir)
    try:
        if not os.path.exists(wdir):
            os.makedirs(wdir)
    except OSError as e:
        print("Ignore OSError from makedirs(work):")
        print(e)
        pass
    url = args.url_file
    filename = url[url.rfind("/")+1:]
    print('Filename is', filename)
    image_file = wget( url, os.path.join(wdir, filename) )
    print(image_file)

if __name__ == "__main__":
    main()
