PATH="/data/PhD/observation/crab_nebula_processed"
FILTER="R"

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple program to read .fits files')
    parser.add_argument('--inputfolder', default='%s/FLAT/%s'%(PATH,FILTER),
                        help='Fits image to read')
    parser.add_argument('--masterbias', default='%s/masterbias.fit'%(PATH),
                        help='Fits image to read')
    parser.add_argument('--masterdark', default='%s/masterdark.fit'%(PATH),
                        help='Fits image to read')
    parser.add_argument('--outfile', default='%s/masterflat_%s.fit'%(PATH,FILTER),
                        help='Fits image to read')
    args = parser.parse_args()

    return args

def main():
    import os, glob
    import numpy as np
    import matplotlib.pyplot as plt
    from tools import read_data, write_data, visualize_image

    args = parse_args()
    flatfiles=sorted(glob.glob(os.path.join(args.inputfolder, "*.fit")))    
    flatsprocesslist=[(read_data(f)[0]-read_data(args.masterbias)[0]-read_data(f)[1]*read_data(args.masterdark)[0])for f in flatfiles]
    flatsprocesslist=[f/np.median(f) for f in flatsprocesslist]
    flatcube=np.dstack(flatsprocesslist)
    masterflat= np.float32(np.median(flatcube, axis=2))
    write_data(masterflat, args.outfile)

    fig = plt.figure(dpi=200)
    ax = fig.subplots()
    visualize_image(ax,masterflat, scale='zscale', colorbar=False)
    #fig.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.0)
    fig.savefig(args.outfile.replace(".fit",".png"), bbox_inches='tight', pad_inches=0)

if __name__ == "__main__":
    main()
    
