import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Read catalog program, it is assume a structure of the catalog. First the expousurelist is a column of numbers. Second each expousure have a folder. Third each folder of each expousure have a exp_info_%d.fits, which have the ccd numbers,')
    
    parser.add_argument('--explist', default='',
                        help='txt list with the number identifier of the exposure')
    parser.add_argument('--explists', nargs='+',  type=str, 
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
def read_alldata(catalogpath,  expolist, fields):
    import fitsio
    import numpy as np

    inpath = os.path.expanduser(catalogpath)
    if not os.path.exists(inpath):
        print('The path of the catalog does not exist!')
        return None
    
    exps = load_explist(expolist) 

    keys = fields
    
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
def read_alldata2(catalogpath,  expolist,  fields):
    import fitsio
    import numpy as np

    inpath = os.path.expanduser(catalogpath)
    if not os.path.exists(inpath):
        print('The path of the catalog does not exist!')
        return None
    
    exps = load_explist(expolist) 

  
    keys = fields +  ['expnum']
    all_data = { key : [] for key in keys }

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
        for key in fields:
            all_data[key].append(expinfo[key])
        all_data['expnum'].append(len(expinfo[key]) * [expnum])

    all_data_final = { key : [] for key in keys }         
    for key in keys:
        all_data_final[key] = np.concatenate(all_data[key])         

    return all_data_final
#Version BBC cartalog
def read_alldata3(catalogpath, expolist, fields):
    import fitsio
    import numpy as np

    inpath = os.path.expanduser(catalogpath)
    if not os.path.exists(inpath):
        print('The path of the catalog does not exist!')
        return None

    nums = load_explist(expolist)
    
    keys = fields
    
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
#Usually is not a good idea save all the catalog in memory
def read_somedata(catalogpath,  expolist):
    import fitsio
    import numpy as np
    import pandas
    inpath = os.path.expanduser(catalogpath)
    if not os.path.exists(inpath):
        print('The path of the catalog does not exist!')
        return None

    exps = load_explist(expolist)
    
    names =  ['expnum', 'musestars',  'mtotalstars']
    formats = ['i4', 'i4',  'i4' ]
    dtype = dict(names = names, formats=formats)
    outdata = np.recarray((len(exps), ), dtype=dtype)
    nstarslist =  []
    explist =  []
    tstarslist =  []

    
    for exp in sorted(exps):
        #print('Start work on exp = ',exp)
        expnum = int(exp)
        #print('expnum = ',expnum)
        indir = os.path.join(inpath, exp)
        try:
            expname = os.path.join(indir, 'exp_psf_cat_%d.fits'%expnum)  
            expstars = fitsio.read(expname, ext='stars')
            data = expstars.astype(expstars.dtype.newbyteorder('='))
            ccdnums = data['ccdnum'].astype(int)
            df = pandas.DataFrame(data)
            #expcat = fitsio.read(expname, ext='info')
            # print('File exp_psf_cat %d.  sucessfully read'%expnum) 
        except (OSError, IOError):
            print('Unable to open exp_psf_cat %s.  Skipping this file.'%expinfo)

        boo = (df['ccdnum'] == 28)     
        tstarslist.append(len(df[boo]))
        boolean = (df['ccdnum'] == 28) & (df['use']== True)
        nstarslist.append( len( df[boolean] ) )
        explist.append(expnum)
        

    values =  []
    values.append(explist)
    values.append(nstarslist)
    values.append(tstarslist)
    for key, i in zip(names, range(len(names))):
        outdata[key] = values[i]

    file_name = "stars.fits"
    write_fit(outdata,  file_name)

#Implenting Mikes method
def read_somedata2(catalogpath,  expolist):
    from read_psf_cats import read_data


    exps = load_explist(expolist) 
    exps = sorted(exps)
    
    keys = ['use']
    data, bands, tilings = read_data(exps, catalogpath, keys, limit_bands='riz', prefix='piff',
                                     use_reserved=False, frac=1.)

    file_name = "stars2.fits"
    print('Finished reading data')
    write_fit(data,  file_name)

def wsc_range(data):
    import numpy as np
    ra_min =  np.min(data['ra'])
    ra_max =  np.max(data['ra'])
    dec_min =  np.min(data['dec'])
    dec_max =  np.max(data['dec'])
    wscr =  [ra_min, ra_max,  dec_min,  dec_max]
    print ( wscr )
    return wscr
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
def plotRaDecs(data,  names,  outname): 
    import matplotlib 
    matplotlib.use('Agg')
    import matplotlib.pylab as pl
    pl.rcParams['agg.path.chunksize'] = 100000
    from astropy.coordinates import SkyCoord
    from astropy import units

    #pl.figure()
    colors = ['k', 'b', 'r', 'g', 'm', 'grey' , 'y', 'c', 'olive', 'darkblue']
    for na, i  in zip(names, range(len(names)) ):
        coords = SkyCoord(ra=data[na]['ra'], dec=data[na]['dec'], unit='degree')
        ra = coords.ra.wrap_at(180 * units.deg)
        #ra = coords.ra.radian
        dec = coords.dec
        pl.plot( ra, dec, color=colors[i], marker=',' )
    pl.xlabel('R.A')
    pl.ylabel('DEC')
    pl.legend()
    pl.grid()
    pl.savefig(outname, dpi=1500)
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
    from astropy.io import fits
    import fitsio
    import numpy as np
    import pandas
    from os.path import basename
    args = parse_args()

    read_somedata2(args.inpath,   args.explist)


    #APP
    '''
    wsc_range(data)
    print('Data was read succesfully')
    plotRaDec(data,  args.outname)
    plotSkymap(data, 'sky' +  args.outname)
    print('plot succesfully done')
    '''

    #APP
    '''
    data = {}
    names = []
    for key in args.explists:
        basen = os.path.splitext(basename(key))[0]
        names.append(basen)
        data[basen ] = read_alldata2(args.inpath,  key,  args.fields)
    print('Data was read succesfully')
    
    plotRaDecs(data, names , args.outname)
    #plotSkymap(data, 'sky' +  args.outname)
    print('plot succesfully done')
    '''
    
    
    

    
if __name__ == "__main__":
    main()
