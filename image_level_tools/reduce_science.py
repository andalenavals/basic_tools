REDUCE=[True,False][0]
RUNSEX=[True,False][0]
RUNSCAMP=[True,False][0]
RUNSWARP=[True,False][0]
BACK,SEG=[False,False]

#PATH="/data/PhD/observation/crab_nebula_processed"
#RA='05h34m33s'
#DEC='+22d00m52.7s'
PATH="/data/PhD/observation/owl_nebula_processed"
#RA='11h14m47.724s'
#DEC='+55d01m8.672s'
RA='11h14m47.715s'
DEC='+55d01m8.682s'

#FILTER="LUM"
#FILTER="HA"
#FILTER="G"
FILTER="B"

from astropy.coordinates import SkyCoord
coord = SkyCoord(ra=RA, dec=DEC, frame='icrs')
npix_x=3072
npix_y=2048
xref=npix_x//2
yref=npix_y//2
pixel_scale=0.413 #arsec
binnig=1
arcsec=1./3600 #deg

cd1_1=binnig*pixel_scale*arcsec
cd1_2=0.0
cd2_1=0.0
cd2_2=binnig*pixel_scale*arcsec

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A simple program to read .fits files')
    parser.add_argument('--inputfolder', default='%s/SCIENCE/%s'%(PATH,FILTER),
                        help='Fits image to read')
    parser.add_argument('--masterbias', default='%s/masterbias.fit'%(PATH),
                        help='Fits image to read')
    parser.add_argument('--masterdark', default='%s/masterdark.fit'%(PATH),
                        help='Fits image to read')
    parser.add_argument('--masterflat',
                        #default='%s/masterflat_%s.fit'%(PATH,"G"),#%(PATH,FILTER),
                        default='%s/masterflat_%s.fit'%(PATH,FILTER),
                        help='Fits image to read')
    parser.add_argument('--outfolder', default='%s/SCIENCE_RED/%s'%(PATH,FILTER),
                        help='Fits image to read')
    parser.add_argument('--sex_args',
                        default='/data/git_repositories/basic_tools/image_level_tools/sexconf/sexconf.yaml',
                        help='yaml config define sextractor inputfiles, sex_config, sex_params, sex_filter, sex_nnw')
    parser.add_argument('--scamp_args',
                        default='/data/git_repositories/basic_tools/image_level_tools/scampconf/scamp.conf',
                        help='yaml config define scamp inputfiles, scamp_bin, scamp_config, globar_header')
    parser.add_argument('--scamp_global_header',
                        default='/data/git_repositories/basic_tools/image_level_tools/scampconf/configs/default2.ahead',
                        help='yaml config define scamp global config')
    parser.add_argument('--swarp_args',
                        default='/data/git_repositories/basic_tools/image_level_tools/swarpconf/swarp.conf',
                        help='yaml config define swarp inputfiles swarp_bin, swarp_config')
    args = parser.parse_args()

    return args

def make_dir(dirname):
    import os
    try:
        if not os.path.exists(dirname):
            os.makedirs(dirname)
    except OSError:
        if not os.path.exists(dirname): raise
    
def read_yaml(filename):
    import yaml
    try:
        with open(filename) as file:
            sexargs = yaml.load(file, Loader=yaml.FullLoader)
    except OSError :
            with open(filename) as file: raise
    assert sexargs is not None
    return sexargs


def main():
    import os, glob
    import numpy as np
    import matplotlib.pyplot as plt
    from tools import read_data, write_data, visualize_image, plot_image
    from tools import make_global_header, run_sex, run_scamp, run_swarp

    args = parse_args()
        
    sexargs=read_yaml(args.sex_args)
    scampargs=read_yaml(args.scamp_args)
    swarpargs=read_yaml(args.swarp_args)

    make_global_header(args.scamp_global_header, coord, xref, yref, cd1_1, cd1_2, cd2_1, cd2_2 )
    scampargs.update({"global_header":args.scamp_global_header})
    
    masterbias=read_data(args.masterbias)[0]
    masterdark=read_data(args.masterdark)[0]
    masterflat=read_data(args.masterflat)[0]    
    files=sorted(glob.glob(os.path.join(args.inputfolder, "*.fit")))
    make_dir(args.outfolder)
    
    for f in files:
        if REDUCE:
            s_proc=(read_data(f)[0]-masterbias-read_data(f)[1]*masterdark)/masterflat
            outfilename=os.path.join(args.outfolder, os.path.basename(f))
            write_data(s_proc, outfilename )
        if RUNSEX:
            if BACK:
                chekcatname=os.path.join(args.outfolder, os.path.basename(f).replace(".fit","_back.fit"))
                sexargs.update({"check_file":chekcatname, "check_type":"BACKGROUND"})
            #chekcatname=os.path.join(args.outfolder, os.path.basename(f).replace(".fit","rms.weight.fit"))
            #sexargs.update({"check_file":chekcatname, "check_type":"MINIBACK_RMS"})
            if SEG:
                chekcatname=os.path.join(args.outfolder, os.path.basename(f).replace(".fit", "_seg.fit"))
                sexargs.update({"check_file":chekcatname, "check_type":"SEGMENTATION"})
                
            sexcatname=os.path.join(args.outfolder, os.path.basename(f).replace(".fit", ".cat"))
            run_sex(outfilename, sexcatname,  **sexargs)

    if RUNSCAMP:
        cats=glob.glob(os.path.join(args.outfolder, "*.cat"))
        #index=[0,2]
        index=list(range(len(cats)))
        cats=np.array(cats)[index]
        cats_string=" ".join(cats)
        run_scamp(cats_string, **scampargs)

    if RUNSWARP:
        #index=[0,1] #list(range(len(files))) #[0,2]
        index=list(range(len(files)))
        reducedfiles=np.array(sorted(glob.glob(os.path.join(args.outfolder, "*.fit"))))[index]
        redfiles_string=" ".join(reducedfiles)
        imgout_name= os.path.join(PATH, "median_%s.fit"%(FILTER))
        wout_name= os.path.join(PATH, "median_%s_weight.fit"%(FILTER))
        center="%s,%s"%(RA.replace("h", ":").replace("m", ":").replace("s", ""), DEC.replace("d", ":").replace("m", ":").replace("s", ""))
        img_size="3500,3500"
        swarpargs.update({"imgout_name":imgout_name, "wout_name":wout_name,"center":center, "img_size":img_size})
        run_swarp(redfiles_string, **swarpargs )

    
        fig = plt.figure(dpi=200)
        ax = fig.subplots()
        coadd=read_data(imgout_name)[0]
        visualize_image(ax,coadd, scale='zscale', colorbar=False)
        #fig.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.0)
        fig.savefig(imgout_name.replace(".fit",".png"), bbox_inches='tight', pad_inches=0)
        #plot_image(coadd)
    
    
if __name__ == "__main__":
    main()
    
