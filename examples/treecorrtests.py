def pretty_rho(meanr, rho, sig=None,  legend=None, lfontsize=24, color='black', marker='o', ylabel='Correlation'):
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    plt.plot(meanr, rho, color=color, label=legend)
    plt.plot(meanr, -rho, color=color, ls=':')
    #rigth quadrants
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color=color, ls='', marker=marker,  capsize=2)
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color=color, ls='', marker=marker,  capsize=2)
    #leftquadrants
    plt.errorbar( -meanr, rho, yerr=sig, color=color,  marker='^',  capsize=2)
    plt.errorbar( -meanr,-rho, yerr=sig, color=color,  marker='^', ls=':', capsize=2)
 
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlabel(r'$\theta$ (rad)', fontsize=24)
    plt.ylabel(ylabel, fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()

def main():
    import treecorr
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
      
    ax =  np.array([0, 1, 2])
    ay =  np.array([0, 0, 0])
    a1 = np.array([ -1,-2, 3])
    a2 = np.array([ - 4, 5, 6])
    acat = treecorr.Catalog(x=ax, y=ay, g1=a1,  g2=a2)

    bx =  np.array([0, 1.1, 2.1])
    by =  np.array([0, 0, 0])
    b1 = np.array([ -7,8, 9])
    b2 = np.array([ -10, 11,- 12])
    bcat = treecorr.Catalog(x=bx, y=by, g1=b1,  g2=b2)

    
    gg = treecorr.GGCorrelation(min_sep=0.1, max_sep=3, bin_size=0.1)
    gg.process(acat, bcat)
    
    var = np.array((2*gg.varxi).tolist())
    npairs = np.array(gg.npairs.tolist())
    meanr = np.array(gg.meanlogr.tolist())
    xip_real = np.array(gg.xip.tolist())
    xip_im = np.array(gg.xip_im.tolist())
    xim_real = np.array(gg.xim.tolist())
    xim_im =  np.array(gg.xim_im.tolist())

    print(meanr)
    print(xip_real)
    #print(xip_im)
    #print(xim_real)
    #print(xim_im)
    
    plt.clf()
    pretty_rho(meanr, npairs, sig=np.zeros(len(meanr)),  legend=None, lfontsize=24, color='black', marker='o', ylabel='Npairs')
    #plt.xlim( [0.01,10.] )
    fname =  '/home/dfa/sobreira/alsina/basic_tools/examples/npairs.pdf'
    plt.savefig(fname)
    
    plt.clf()
    pretty_rho(meanr, xip_real, sig=np.zeros(len(meanr)) , legend='xip_real', lfontsize=24, color='black', marker='o')
    #pretty_rho(meanr, xip_real, sig=np.sqrt(var) , legend='xip_real', lfontsize=24, color='black', marker='o')
    plt.xlim( [0.01,5.] )
    fname =  '/home/dfa/sobreira/alsina/basic_tools/examples/xip_real.pdf'
    plt.savefig(fname)

    '''
    plt.clf()
    pretty_rho(meanr, xip_im, legend='xip_im', lfontsize=24, color='black', marker='o')
    fname =  '/home/dfa/sobreira/alsina/basic_tools/examples/xip_im.pdf'
    plt.savefig(fname)

    plt.clf()
    pretty_rho(meanr, xim_real, legend='xim_real', lfontsize=24, color='black', marker='o')
    fname =  '/home/dfa/sobreira/alsina/basic_tools/examples/xim_real.pdf'
    plt.savefig(fname)

    plt.clf()
    pretty_rho(meanr, xim_im, legend='xim_im', lfontsize=24, color='black', marker='o')
    fname =  '/home/dfa/sobreira/alsina/basic_tools/examples/xim_im.pdf'
    plt.savefig(fname)
    '''

if __name__ == "__main__":
    main()
