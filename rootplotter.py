def plotRaDec(data,  name): 
    import matplotlib 
    matplotlib.use('Agg')
    import matplotlib.pylab as pl
    from astropy.coordinates import SkyCoord
    from astropy import units

    #pl.figure()
    coords = SkyCoord(ra=data['telra'], dec=data['teldec'], unit='degree')
    ra = coords.ra.wrap_at(180 * units.deg)
    #ra = coords.ra.radian
    dec = coords.dec
    pl.plot( ra, dec, 'k.', markersize=0.05 )
    pl.xlabel('R.A')
    pl.ylabel('DEC')
    pl.legend()
    pl.grid()
    pl.savefig(name, dpi=150)
def plotRaDec2D(data,  name): 
    import matplotlib 
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.coordinates import SkyCoord
    from astropy import units
    import numpy as np
    from matplotlib.colors import LogNorm

    data = data[data['rho2p']>0]
    coords = SkyCoord(ra=data['telra'], dec=data['teldec'], unit='degree')
    ra = coords.ra.wrap_at(180 * units.deg)
    #ra = coords.ra.radian
    dec = coords.dec
    rho2p = data['rho2p']
 
    fig = plt.figure()
    #rho2p = (rho2p -np.nanmin(rho2p)) /(np.nanmax(rho2p) - np.nanmin(rho2p))
 
    sct = plt.scatter(ra, dec, s=1, marker='p' , cmap='gnuplot', c=rho2p, norm=LogNorm( vmin=np.nanmin(rho2p), vmax=np.nanmax(rho2p))  )
    plt.colorbar(sct)
    plt.xlabel('R.A')
    plt.ylabel('DEC')
    #pl.legend()
    #pl.grid()
    plt.savefig(name, dpi=150)
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

    gr = TGraph( len(data['teldec']) - 1, data['teldec'], data['telra'] )   
 
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

def plotTH1(data, field, Nbins , xtitle, ytitle, outfile_name):
    from ROOT import TCanvas, TGraph,  TH2F,  TH2,  TH1,  TH1F
    from ROOT import gROOT, gSystem,  Double,  TAxis,  gStyle
    import numpy as np

    ignore =  (data[field] ==-999.) |  (data[field] == None) | ( np.isnan(data[field])) 
    data = data[~ignore]
    
    gStyle.SetOptStat(0)
    c1 =  TCanvas('c1', '', 1000,1000)
    c1.SetBottomMargin( 0.15 )
    c1.SetTopMargin( 0.05 )
    c1.SetLeftMargin( 0.15 )
    c1.SetRightMargin( 0.15 )
    c1.Divide(1, 2)

    #print (len(data[field]))
    #print (list(data[field]))
    nbins = Nbins
    minx = Double( np.min(data[field]) )
    maxx = Double( np.max(data[field]) )
    #print (minx, ' ', maxx, ' ',  nbins)
   
    h0 = TH1F('h0', '',  nbins ,  minx ,  maxx )
    for k in range(len(data)):
        h0.Fill(data[field][k]) 
    
    h = TH1F('h', '',  nbins ,  minx ,  maxx )
    mh = TH1F('mh', '',  nbins ,  minx ,  maxx )
    for i in range(1, nbins + 1):
        binxlow = h.GetBinLowEdge(i)
        binxup =  h.GetXaxis().GetBinUpEdge(i)
        bool1 = ( data[field] <= binxup ) & (binxlow<data[field]) & (data['rho2p']>0)
        bool2 = ( data[field] <= binxup ) & (binxlow<data[field]) & (data['rho2p']<0)
        rho2p = data['rho2p'][bool1]
        mrho2p = data['rho2p'][bool2]
        if (len(rho2p) == 0):
            meanr2p = 0
        else:
            meanr2p =  np.mean(rho2p)
        if (len(mrho2p) == 0):
            meanmr2p = 0
        else:
            meanmr2p =  np.mean(mrho2p)

        #print( - meanmr2p)
        h.SetBinContent( i , meanr2p)
        mh.SetBinContent( i , -meanmr2p)
        
    c1.cd(1)
    h0.Draw('')
    h0.GetYaxis().CenterTitle()
    h0.GetYaxis().SetTitle("Entries")
    h0.SetTitleSize(0.045, "y");
    #h.SetTitleOffset( 0.6 , "y"); 
    c1.cd(2)
    c1.cd(2).SetLogy()
    h.Draw('')
    h.SetTitleOffset(  0.6 , "y"); 
    h.SetTitleOffset(  0.6 , "x"); 
    h.GetYaxis().CenterTitle()
    h.GetYaxis().SetTitle(ytitle)
    h.GetXaxis().CenterTitle()
    h.GetXaxis().SetTitle(xtitle)
    h.SetTitleSize(0.065, "x"); 
    h.SetTitleSize(0.065, "y");
    mh.SetLineColor(1)
    mh.SetLineStyle(2)
    mh.Draw('same')
    c1.Print(outfile_name)
    
def main():
    from astropy.io import fits
    import fitsio
    import numpy as np
    import pandas
   
    data = fitsio.read('y3a1-v29_rho2byexposure_extended_final.fit')
    data = data.astype(data.dtype.newbyteorder('='))
    #df = pandas.DataFrame(data)

    #plotRaDecRoot(data, 'footprint.pdf')
    plotRaDec2D(data, 'footprint.pdf')
    
    '''
    field =  'teldec'
    units =  '[deg]'
    plotTH1(data, field , 500, field + units, "#bar{#rho_{2}}", "rho2_vs_" + field + ".pdf")
    '''

    
    
    

    
if __name__ == "__main__":
    main()
