import os
#PATH="/data/PhD/observation/crab_nebula_processed/"
PATH="/data/PhD/observation/owl_nebula_processed/"

#RIMG=os.path.join(PATH, "median_R.fit")
RIMG=os.path.join(PATH, "median_HA.fit")
GIMG=os.path.join(PATH, "median_G.fit")
BIMG=os.path.join(PATH, "median_B.fit")
LUMIMG=os.path.join(PATH, "median_LUM.fit")



def main():
    import numpy as np
    import matplotlib.pyplot as plt
    from tools import read_data
    import astropy.visualization as vis
    import matplotlib.image
    from matplotlib import transforms
    import astropy.convolution as conv
    from scipy import ndimage
    import tools
    

    rarray=read_data(RIMG)[0]
    garray=read_data(GIMG)[0]
    barray=read_data(BIMG)[0]
    larray=read_data(LUMIMG)[0]

    #rarray=0.5*(rarray+larray)
    #garray=0.5*(garray+larray)
    #barray=0.4*(barray+larray)
    for a in [rarray, garray, barray, larray]:
        print(np.min(a),np.max(a), np.median(a), np.mean(a), np.std(a))

    #rarray-=np.min(rarray)
    #rarray/=np.max(rarray)

    #garray-=np.min(garray)
    #garray/=np.max(garray)

    #barray-=np.min(barray)
    #barray/=np.max(barray)
 

    if False:
        BIN=2
        operation='mean'
        newshape=list((np.array(rarray.shape[:2])/BIN).astype(int))
        rarray=tools.bin_ndarray(rarray, newshape, operation=operation)
        garray=tools.bin_ndarray(garray, newshape, operation=operation)
        barray=tools.bin_ndarray(barray, newshape, operation=operation)

    if True:
        val=0
        std=5
        barray[np.where(barray<val)]=np.nan
        barray=conv.interpolate_replace_nans(barray, conv.Gaussian2DKernel(x_stddev=std))
        rarray[np.where(rarray<val)]=np.nan
        rarray=conv.interpolate_replace_nans(rarray, conv.Gaussian2DKernel(x_stddev=std))
        garray[np.where(garray<val)]=np.nan
        garray=conv.interpolate_replace_nans(garray, conv.Gaussian2DKernel(x_stddev=1))
    
    if False:
        size=0.8
        rkernel = conv.Tophat2DKernel(radius=size)
        rarray = conv.convolve(rarray, rkernel)
        bkernel = conv.Tophat2DKernel(radius=size)
        barray = conv.convolve(barray, bkernel)
        gkernel = conv.Tophat2DKernel(radius=size)
        garray = conv.convolve(garray, gkernel)



    rgbimage = vis.make_lupton_rgb(rarray, garray, barray, Q=5, stretch=10)
    #rgbimage = make_lupton_rgb(rarray, 0.5*(barray+rarray), barray, Q=5, stretch=30)
   
    combinationname=os.path.join(PATH, "OWLnebulaFOV.jpg")
    #matplotlib.image.imsave(combinationname, rgbimage, origin="lower", dpi=5000)
    fig = plt.figure()
    ax = fig.subplots()
    rgbimage_50=ndimage.rotate(rgbimage,-48,reshape=True)
    print(rgbimage_50.shape)
    rgbimage50=rgbimage_50[1550:-1600,900:-1100,:]
    print(rgbimage50.shape)
    #newshape=list((np.array(rgbimage50.shape[:2])/2).astype(int))+[3]
    #rgbimage50=tools.bin_ndarray(rgbimage50, newshape, operation='sum')
    #rgbimage50=ndimage.zoom(rgbimage50,2)
    ax.imshow(rgbimage50, interpolation="nearest",interpolation_stage="rgba")
    ax.set_axis_off( )
    plt.savefig(combinationname, dpi=400,bbox_inches = 'tight',pad_inches = 0.0)
    

    
    '''
    fig = plt.figure(dpi=200)
    ax = fig.subplots()
    RGA_array=np.stack([rarray, garray, barray],axis=2)
    #normalizing

    
    minarray=np.min(RGA_array)
    RGA_array-=minarray
    maxarray=np.max(RGA_array)
    RGA_array/=maxarray
    #norm = vis.ImageNormalize(RGA_array, interval=vis.ZScaleInterval(),stretch=vis.LinearStretch())
    #norm = vis.ImageNormalize(RGA_array, interval=vis.ZScaleInterval(),stretch=vis.SqrtStretch())
    #norm = vis.ImageNormalize(RGA_array, interval=vis.ZScaleInterval(),stretch=vis.LogStretch())
    interval=vis.ZScaleInterval(nsamples=1000, contrast=0.25,
                           max_reject=0.5,
                           min_npixels=5, krej=2.5,
                           max_iterations=5)
    vmin,vmax=interval.get_limits(RGA_array)
        
    print(RGA_array.shape)
    combname=os.path.join(PATH, "OWLnebulaFOV_t.png")
    ax.imshow(RGA_array, vmin=vmin,vmax=vmax)
    ax.set_axis_off( )
    plt.savefig(combname, dpi=200,bbox_inches = 'tight',pad_inches = 0.01)
    '''
if __name__ == "__main__":
    main()
    
