import os
PATH="/data/PhD/observation/crab_nebula_processed/"
RIMG=os.path.join(PATH, "median_R.fit")
GIMG=os.path.join(PATH, "median_LUM.fit")
BIMG=os.path.join(PATH, "median_B.fit")

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple program to read .fits files')
    parser.add_argument('--inputfolder', default='%sBIAS'%(PATH),
                        help='Fits image to read')
    parser.add_argument('--outfile', default='%s/masterbias.fit'%(PATH),
                        help='Fits image to read')
    args = parser.parse_args()

    return args

def main():
    import numpy as np
    import matplotlib.pyplot as plt
    from tools import read_data
    from astropy.visualization import make_lupton_rgb
    import matplotlib.image
    from astropy.convolution import Gaussian2DKernel, convolve
    
    args = parse_args()

    rarray=read_data(RIMG)[0]
    #garray=read_data(GIMG)[0]
    barray=read_data(BIMG)[0]

    
    #rkernel = Gaussian2DKernel(x_stddev=0.5)
    #rarray = convolve(rarray, rkernel)
    #bkernel = Gaussian2DKernel(x_stddev=0.5)
    #barray = convolve(barray, bkernel)

    #rgbimage = make_lupton_rgb(rarray, garray, barray, Q=5, stretch=30)
    rgbimage = make_lupton_rgb(rarray, 0.5*(barray+rarray), barray, Q=5, stretch=30)

    combinationname=os.path.join(PATH, "crabnebula.png")
    matplotlib.image.imsave(combinationname, rgbimage, origin="lower", dpi=200)
    #fig = plt.figure(dpi=200)
    #ax = fig.subplots()
  

if __name__ == "__main__":
    main()
    
