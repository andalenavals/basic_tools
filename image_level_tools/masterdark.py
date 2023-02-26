PATH="/data/PhD/observation/crab_nebula_processed/"

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple program to read .fits files')
    parser.add_argument('--inputfolder', default='%sDARK'%(PATH),
                        help='Fits image to read')
    parser.add_argument('--masterbias', default='%s/masterbias.fit'%(PATH),
                        help='Fits image to read')
    parser.add_argument('--outfile', default='%s/masterdark.fit'%(PATH),
                        help='Fits image to read')
    args = parser.parse_args()

    return args

def main():
    import os, glob
    import numpy as np
    import matplotlib.pyplot as plt
    from tools import read_data, write_data, visualize_image

    args = parse_args()
    darkfiles=sorted(glob.glob(os.path.join(args.inputfolder, "*.fit")))
    dark_list=[(read_data(f)[0]-read_data(args.masterbias)[0])/read_data(f)[1] for f in darkfiles]
    darkcube=np.dstack(dark_list)
    masterdark= np.float32(np.median(darkcube, axis=2))
    write_data(masterdark, args.outfile)

    fig = plt.figure(dpi=200)
    ax = fig.subplots()
    visualize_image(ax,masterdark, scale='zscale', colorbar=False)
    #fig.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.0)
    fig.savefig(args.outfile.replace(".fit",".png"), bbox_inches='tight', pad_inches=0)

if __name__ == "__main__":
    main()
    
