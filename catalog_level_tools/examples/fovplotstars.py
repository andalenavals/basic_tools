def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='This program plot the position of the stars in a exp file in FOV. I am interested into distinguish stars with high ellipticity ')
    
    parser.add_argument('--file', default='/home/dfa/sobreira/alsina/basic_tools/exp_psf_cat_282956.fits',
                        help='File where is the information for the filter')

    args = parser.parse_args()

    return args

    
def plotScatter(data, xfield, yfield, xtitle='', ytitle='', title='', xlog=False, ylog=False, grid=False, savefig=False):
    import fitsio
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.coordinates import SkyCoord
    from astropy import units
    import numpy as np
    from toFocal import toFocal


    print("the number of stars is",  len(data))
    ignorex =  (data[xfield] ==-999.) |  (data[xfield] == None) | ( np.isnan(data[xfield])) 
    data = data[~ignorex]
    ignorey =  (data[yfield] ==-999.) |  (data[yfield] == None) | ( np.isnan(data[yfield])) 
    data = data[~ignorey]
    #data = data[data['piff_flag']==0]
    print("the number of piff stars is",  len(data))

    
    plt.figure()
    #for k in range(len(data)):
    #    x, y = toFocal(data['ccdnum'][k], data[xfield][k], data[yfield][k] )
    #    plt.plot( x , y , color='blue', marker= '.',  markersize=0.5,  linewidth=0 )
    data1 = data[np.sqrt(data['obs_e1']**2+data['osb_e2']**2)<0.1]
    x1, y1 =  toFocal(data1['ccdnum'], data1[xfield], data1[yfield] )
    plt.plot( x1 , y1 , color='blue', marker= '.',  markersize=0.5,  linewidth=0 )
    data2 = data[np.sqrt(data['osb_e1']**2+data['obs_e2']**2)>=0.1]
    x2, y2 =  toFocal(data2['ccdnum'], data2[xfield], data2[yfield] )
    plt.plot( x2 , y2 , color='black', marker= '.',  markersize=0.6,  linewidth=0 )
    

    if(ylog):
        plt.yscale('log')
    if(xlog):
        plt.xscale('log')
    if(grid): 
        plt.grid()
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(title)
    if(savefig):
        plt.savefig(yfield+'_vs_'+xfield +'.pdf', dpi=150)  
    
    
def main():
    from astropy.io import fits
    import fitsio
    import numpy as np
    import pandas

    args=parse_args()
    data = fitsio.read(args.file)
    data = data.astype(data.dtype.newbyteorder('='))


    
    plotScatter(data=data, xfield='x' , yfield='y', savefig=True)


 
    
    

    
if __name__ == "__main__":
    main()
