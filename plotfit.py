def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Read Fits program')
    
    parser.add_argument('--file', default='',
                        help='Fits image to read')
    parser.add_argument('--query', default='',
                        help='SQL format to query')
    parser.add_argument('--full', default=False, action='store_const', const=True,
                        help='full info flag')
    parser.add_argument('--headers', nargs='+',  type=str, 
                        help='Print only some columns')
    args = parser.parse_args()

    return args

def testinputfiles(args):
    if args.file != '':
        print ('\n File ',args.file, 'sucessfully read. OK! \n' )
    else:
        print ('No input file1')

def PrettyPrint(data):
    from prettytable import PrettyTable
    x = PrettyTable(data.dtype.names)    
    for row in data:
        x.add_row(row)
        # Change some column alignments; default was 'c'
        #x.align['column_one'] = 'r'
        #x.align['column_3'] = 'l'
    print(x)
  
def Print_to_image(data):
    import matplotlib.pyplot as plt
    from astropy.visualization import astropy_mpl_style
    from scipy.misc import imsave
    plt.style.use(astropy_mpl_style)
    plt.figure()
    plt.imshow(data,  cmap='gray')
    plt.colorbar()
    plt.show()
    imsave('image.png', data)


def unpack_file(file_name):
    import os
    """Create the unpacked file in the work directory if necessary.

    If the unpacked file already exists, then a link is made.
    Otherwise funpack is run, outputting the result into the work directory.
    """
    print('unpack ',file_name)

    img_file = os.path.splitext(file_name)[0]
    if os.path.lexists(img_file):
        print('   %s exists already.  Removing.'%img_file)
        os.remove(img_file)
    print('   unpacking fz file')
    cmd = 'funpack -O {outf} {inf}'.format(outf=img_file, inf=file_name)
    print(cmd)
    for itry in range(5):
        run_with_timeout(cmd, 120)
        if os.path.lexists(img_file): break
        print('%s was not properly made.  Retrying.'%img_file)
        time.sleep(10)

    if not os.path.lexists(img_file):
        print('Unable to create %s.  Skip this file.'%img_file)
        return None

    return img_file
 
def run_with_timeout(cmd, timeout_sec):
    # cf. https://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
    import subprocess, shlex
    from threading import Timer

    proc = subprocess.Popen(shlex.split(cmd))
    kill_proc = lambda p: p.kill()
    timer = Timer(timeout_sec, kill_proc, [proc])
    try:
        timer.start()
        proc.communicate()
    finally:
        timer.cancel()


def main():
    import fitsio
    from astropy.io import fits


    args = parse_args()
    testinputfiles(args)

    '''   
    hdu =  fits.open(args.file)
    data = hdu[1].data
    Print_to_image(data)
    '''

    '''
    import fitsio
    image_file =  '1.fz'
    bkg_file = '2.fz'
    with fitsio.FITS(image_file, 'rw') as f:
                    bkg = fitsio.read(bkg_file)
                    img = f[1].read()
                    #print(img)
                    img -= bkg
                    f[1].write(img)
    unpack_image_file = unpack_file('1.fz')
    #print('unpack_image_file =' , unpack_image_file)
    '''

    import galsim
    full_image = galsim.fits.read('1.fz', hdu=1)
    print(full_image)

    stamp_size = 48
    
    #hdu =  fits.open(args.file)
    #data = hdu[1].data
    #Print_to_image(data)

if __name__ == "__main__":
    main()
    


    
    
 
