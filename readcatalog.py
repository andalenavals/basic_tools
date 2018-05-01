import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Read catalog program, it is assume a structure of the catalog. First the expousurelist is a column of numbers. Second each expousure have a folder. Third each folder of each expousure have a exp_info_%d.fits, which have the ccd numbers,')
    
    parser.add_argument('--explist', default='',
                        help='txt list with the number identifier of the exposure')
    parser.add_argument('--explists', nargs='+',  
                        help='lists of txt files with the number of exposure to analyse')
    parser.add_argument('--fields', nargs='+',  type=str, 
                        help='list of fields you want to read from the catalog')
    parser.add_argument('--inpath', default='',
                        help='Place where input catalogs is, it is assumed that each expousure have a folder')
    parser.add_argument('--outname', default='', type=str, 
                        help='Name of the output image example zone01.png')

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

#Version for V10 psf catalog
def read_alldata(args):
    import fitsio
    import numpy as np

    inpath = os.path.expanduser(args.inpath)
    if not os.path.exists(inpath):
        print('The path of the catalog does not exist!')
        return None
    
    exps = load_explist(args.explist) 

    keys = args.fields
    
    all_data = { key : [] for key in keys }
    all_keys = keys

    for exp in exps:
        #print('Start work on exp = ',exp)
        expnum = int(exp)
        #print('expnum = ',expnum)
        indir = os.path.join(inpath, exp)
        expinfo = fitsio.read(os.path.join(indir, 'exp_info_%d.fits'%expnum))

        if expnum not in expinfo['expnum']:
            print('expnum is not in expinfo!')
            print('expinfo[expnum] = ',expinfo['expnum'])
            print('Could not find information about this expnum.  Skipping ',run,exp)
            continue

        for k in range(len(expinfo)):
            ccdnum = expinfo[k]['ccdnum']
            cat_file = os.path.join(indir, 'psf_cat_%d_%d.fits'%(expnum,ccdnum))
            #print('Reading data from' , cat_file)
            try:
                dataux = fitsio.read(cat_file)
            except (OSError, IOError):
                print('Unable to open cat_file %s.  Skipping this file.'%cat_file) 
            for key in all_keys:
                all_data[key].append(dataux[key])

    all_data_final = { key : [] for key in keys }         
    for key in all_keys:
        all_data_final[key] = np.concatenate(all_data[key])         

    return all_data_final

#Version for V23 psf catalog
def read_alldata2(args):
    import fitsio
    import numpy as np

    inpath = os.path.expanduser(args.inpath)
    if not os.path.exists(inpath):
        print('The path of the catalog does not exist!')
        return None
    
    exps = load_explist(args.explist) 

    keys = args.fields
    
    all_data = { key : [] for key in keys }
    all_keys = keys

    for exp in exps:
        #print('Start work on exp = ',exp)
        expnum = int(exp)
        #print('expnum = ',expnum)
        indir = os.path.join(inpath, exp)
        try:
            expinfo = fitsio.read(os.path.join(indir, 'exp_psf_cat_%d.fits'%expnum))
            # print('File exp_psf_cat %d.  sucessfully read'%expnum) 
        except (OSError, IOError):
            print('Unable to open exp_psf_cat %s.  Skipping this file.'%expinfo) 
        for key in all_keys:
            all_data[key].append(expinfo[key])

    all_data_final = { key : [] for key in keys }         
    for key in all_keys:
        all_data_final[key] = np.concatenate(all_data[key])         

    return all_data_final
#Version BBC cartalog
def read_alldata3(args):
    import fitsio
    import numpy as np

    inpath = os.path.expanduser(args.inpath)
    if not os.path.exists(inpath):
        print('The path of the catalog does not exist!')
        return None

    nums = load_explist(args.explist)
    
    keys = args.fields
    
    all_data = { key : [] for key in keys }
    all_keys = keys

    for num in nums:
        inum = int(num)
        try:
            einfo = fitsio.read(os.path.join(inpath, 'aardvarkv1.0_des_lenscat_s2n20.%d.fit'%inum))
            # print('File exp_psf_cat %d.  sucessfully read'%expnum)
        except (OSError, IOError):
            print('Unable to open aardvarkv1.0_des_lenscat_s2n20.%s.fits. Skipping it.'%inum) 
        for key in all_keys:
            all_data[key].append(einfo[key])

    all_data_final = { key : [] for key in keys }         
    for key in all_keys:
        all_data_final[key] = np.concatenate(all_data[key])         

    return all_data_final
def plotRaDecRoot(data,  name):
    from ROOT import TCanvas, TGraph,  TH2F,  TH2,  TH1
    from ROOT import gROOT, gSystem,  Double
    from array import array
    import numpy as np

    c1 =  TCanvas('c1', '', 800, 600)
    c1.SetBottomMargin( 0.15 )
    c1.SetTopMargin( 0.05 )
    c1.SetLeftMargin( 0.15 )
    c1.SetRightMargin( 0.15 )

    gr = TGraph( len(data['dec']) - 1, data['dec'], data['ra'] )   
 
    gr.SetTitle("")
    gr.GetYaxis().SetTitle("RA")
    gr.GetYaxis().CenterTitle()
    gr.GetXaxis().SetTitle("DEC")
    gr.GetXaxis().CenterTitle()
    gr.GetXaxis().SetTitleOffset(0.9 ) 
    gr.GetYaxis().SetTitleOffset(0.9)
    gr.GetXaxis().SetTitleSize(0.065) 
    gr.GetYaxis().SetTitleSize(0.065)
    
    #gr.GetXaxis().SetLimits( - 11, - 5.4);
    #gr.GetYaxis().SetRangeUser( 21, 24.4);
    #gr.SetMarkerStyle(20)
    #gr.SetMarkerSize(0.1)

    gr.Draw('AP')
    #
    c1.Print(name)

def plotSkymapRoot(data,  name):
    from ROOT import TCanvas, TGraph,  TH2F,  TH2,  TH1
    from ROOT import gROOT, gSystem,  Double
    from ROOT import gApplication
    # gApplication.Run()
    #from rootpy.interactive import wait
    from array import array
    import numpy as np

    
    c1 =  TCanvas('c1', '', 800, 600)
    c1.SetBottomMargin( 0.15 )
    c1.SetTopMargin( 0.05 )
    c1.SetLeftMargin( 0.15 )
    c1.SetRightMargin( 0.15 )
  
    nbins = 10000
    minx = Double( np.min(data['dec']) )
    maxx = Double( np.max(data['dec']) )
    miny = Double( np.min(data['ra']) )
    maxy = Double( np.max(data['ra']) )
    
    print (nbins, minx,  maxx,  miny,  maxy)
    gr = TH2F('gr', '',  nbins ,  minx ,  maxx , nbins, miny, maxy  )
    for i in range(nbins):
        #print(i, data['dec'][i], data['ra'][i] )
        gr.Fill( Double(data['dec'][i] ), Double(data['ra'][i])  )
    #gr.Draw("COLZ")
    gr.Draw('aitoff')
    c1.Print(name)
    
def plotRaDec(data,  name): 
    import matplotlib 
    matplotlib.use('Agg')
    import matplotlib.pylab as pl
    from astropy.coordinates import SkyCoord
    from astropy import units

    #pl.figure()
    coords = SkyCoord(ra=data['ra'], dec=data['dec'], unit='degree')
    ra = coords.ra.wrap_at(180 * units.deg)
    #ra = coords.ra.radian
    dec = coords.dec
    pl.plot( ra, dec, 'ko', markersize=0.05 )
    pl.xlabel('R.A')
    pl.ylabel('DEC')
    pl.legend()
    pl.grid()
    pl.savefig(name, dpi=150)
def plotSkymap(data,  name):
    import matplotlib 
    matplotlib.use('Agg')
    import matplotlib.pylab as pl
    from astropy.coordinates import SkyCoord
    from astropy import units

    coords = SkyCoord(ra=data['ra'], dec=data['dec'], unit='degree', frame='icrs')
    ra = coords.ra.wrap_at(180 * units.deg).radian
    dec = coords.dec.radian
    
    #pl.figure()
    pl.subplot(111, projection="aitoff")  # must use radians
    #pl.subplot(111, projection="lambert")
    # Make sure to bin the regions, otherwise Memory errors appear.
    #pl.plot( ra , dec, 'ko', markersize=0.05  )
    
    color_map = pl.cm.Spectral_r
    image = pl.hexbin(ra, dec, cmap=color_map, gridsize=100, mincnt=1, bins='log')
    pl.colorbar(image, spacing='uniform', extend='max')
    
    pl.xlabel('R.A')
    pl.ylabel('DEC')
    pl.legend()
    pl.grid(True)
    pl.savefig(name, dpi=150)
def main():
    from astropy.io import fits
    import fitsio
    import numpy as np
    import pandas
    from os.path import basename
    args = parse_args()
    '''
    data =  read_alldata2(args)
    print('Data was read succesfully')
    plotRaDec(data,  args.outname)
    plotSkymap(data, 'sky' +  args.outname)
    print('plot succesfully done')
    '''
    keys =  args.explists
    data = { key : [] for key in keys }
    for key in keys:
        basename = os.path.splitext(basename(keys))[0]
        data[basename ] = readalldata2(key)
    
    
    

    
if __name__ == "__main__":
    main()
