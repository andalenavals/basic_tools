import astropy
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1
import matplotlib.axes as maxes
import astropy.visualization
import astropy.io.fits as fits

import logging
logger = logging.getLogger(__name__)

def bin_ndarray(ndarray, new_shape, operation='sum'):
    operation = operation.lower()
    if not operation in ['sum', 'mean']:
        raise ValueError("Operation not supported.")
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                  ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        op = getattr(ndarray, operation)
        ndarray = op(-1*(i+1))
    return ndarray

def plot_image(image): # remove directory and channel_list
    import astropy.visualization as vis
    fig = plt.figure(1)
    plt.tight_layout(pad=0.90, h_pad=None, w_pad=None, rect=None)
    plt.minorticks_on()
    norm = vis.ImageNormalize(image, interval=vis.ZScaleInterval(),stretch=vis.LinearStretch())
    plot = plt.imshow(image,norm=norm,origin='lower',cmap='gray')
    
    
    cbar = fig.colorbar(mappable=plot) #cax=cb_ax)
    cbar.set_label('ADU')
    plt.xlabel('pixel')
    plt.ylabel('pixel')
    plt.savefig('test.pNG',dpi=500,bbox_inches = 'tight',pad_inches = 0.01)
    return 0

def visualize_image(ax, array, cmap=plt.cm.gray, vmin=None, vmax=None, scale=None, percentile=None, 
                    lower_percentile=None, upper_percentile=None,colorbar=True, show_axis=False,
                   barpos='right'):
    if scale is not None:
        if scale=='zscale':
            vis=astropy.visualization.ZScaleInterval(nsamples=1000, contrast=0.25, max_reject=0.5, 
                                             min_npixels=5, krej=2.5, max_iterations=5)
        elif scale=='minmax': vis=astropy.visualization.MinMaxInterval()
        if percentile is not None: vis=astropy.visualization.PercentileInterval(percentile)
        if (lower_percentile is not None )& (upper_percentile is not None):
                vis=astropy.visualization.AsymmetricPercentileInterval(lower_percentile, upper_percentile )       
        vmin,vmax=vis.get_limits(array)
    print("using vmin, vmax: %.3f, %.3f"%(vmin,vmax))
    if not show_axis: ax.set_axis_off( ) #ax.axis('off')
    #ax.axes.xaxis.set_visible(False)
    #ax.axes.yaxis.set_visible(False)
    ec=ax.imshow( array , cmap = cmap , vmin = vmin , vmax =vmax , origin = 'lower')
    if colorbar:
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
        if barpos=='right':
            cax = divider.append_axes(barpos, size="5%", pad=0.1, axes_class=maxes.Axes)
            cbar=plt.gcf().colorbar(ec, cax=cax)
        if barpos=='left':
            if len(ax.get_ylabel()) == 0: pad=0.1
            else: pad=0.9
            cax = divider.append_axes(barpos, size="5%", pad=pad, axes_class=maxes.Axes)
            cbar=plt.gcf().colorbar(ec, cax=cax)
            cbar.ax.yaxis.set_label_position('left')
            cbar.ax.yaxis.set_ticks_position('left')
        if barpos=='top':
            cax = divider.append_axes(barpos, size="5%", pad=0.1, axes_class=maxes.Axes)
            cbar=plt.gcf().colorbar(ec, cax=cax, orientation='horizontal')
            cbar.ax.xaxis.set_label_position('top')
            cbar.ax.xaxis.set_ticks_position('top')
        if barpos=='bottom':
            if len(ax.get_xlabel()) ==0: pad=0.1
            else: pad=0.7
            cax = divider.append_axes(barpos, size="5%", pad=pad, axes_class=maxes.Axes)
            cbar=plt.gcf().colorbar(ec, cax=cax, orientation='horizontal')   
    #ax.grid( color = 'gray')

def read_data(filename):
    hdul = fits.open(filename)
    array = hdul[0].data
    if "EXPOSURE" in hdul[0].header:
        exptime = hdul[0].header["EXPOSURE"]
    elif "EXPTIME" in hdul[0].header:
        exptime = hdul[0].header["EXPTIME"]
    else:
        exptime=None
    hdul.close()
    #print("Done reading ", filepath, ": exptime =", exptime)
    return array, exptime

def write_data(array, filename):
    hdu = fits.PrimaryHDU(array)
    hdu.writeto(filename, overwrite=True)
    print("Wrote", filename, ": dtype =", array.dtype)

def run_swarp(img_files, swarp_bin=None, swarp_config=None, imgout_name=None, wout_name=None, center=None, img_size=None ):
    import subprocess, shlex
    from subprocess import DEVNULL, PIPE
    """Run SWarp
    """
    logger.info('running swarp')

    cmd = "{swarp_bin} {img_files} -c {config} -IMAGEOUT_NAME {imgout_name} -WEIGHTOUT_NAME {wout_name} -CENTER {center} -IMAGE_SIZE {img_size}  ".format(
        swarp_bin=swarp_bin, img_files=img_files, config=swarp_config, imgout_name=imgout_name, wout_name=wout_name, center=center, img_size=img_size)
    
    print("FULL CMD COMMAND: %s"%(cmd))
    res=subprocess.run(shlex.split(cmd))#, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    #res = subprocess.run(shlex.split(cmd), text=True, capture_output=False)
    assert res.returncode ==0


def run_scamp(cat_file, scamp_bin=None, scamp_config=None, global_header=None):
    import subprocess, shlex
    from subprocess import DEVNULL, PIPE
    """Run Scamp
    """
    logger.info('running scamp')

    cmd = "{scamp_bin}  {cat_file} -c {config} -AHEADER_GLOBAL {global_header} ".format(
        scamp_bin=scamp_bin, cat_file=cat_file, config=scamp_config,global_header=global_header )
    
    print("FULL CMD COMMAND: %s"%(cmd))
    #res=subprocess.run(shlex.split(cmd))#, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    res = subprocess.run(shlex.split(cmd), text=True, capture_output=False)
    print(res)
    assert res.returncode ==0

        
def run_sex( img_file, cat_file, check_file=None,check_type=None, sex_bin=None, sex_config=None, sex_params=None, sex_filter=None, sex_nnw=None):

    import subprocess, shlex
    from subprocess import DEVNULL, PIPE
    """Run sextractor, but only if the output file does not exist yet.
    """
    logger.info('running sextractor')
    if check_file is not None: check_file="-CHECKIMAGE_NAME %s"%(check_file)
    else: check_file=""
    if check_type is not None: check_type="-CHECKIMAGE_TYPE %s"%(check_type)
    else: check_type=""

    cmd = "{sex_bin} {img_file} -c {config} -CATALOG_NAME {cat_file} -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME {params} -FILTER_NAME {filter} -STARNNW_NAME {nnw} {check_file} {check_type} ".format(
        sex_bin=sex_bin, img_file=img_file, config=sex_config,
        cat_file=cat_file, params=sex_params, filter=sex_filter,
        nnw=sex_nnw ,check_file=check_file, check_type=check_type)
    
    print("FULL CMD COMMAND: %s"%(cmd))
    res=subprocess.run(shlex.split(cmd))#, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    assert res.returncode ==0
    
    #print(res)
    
    '''
    if res.returncode ==0:
        print(res.stdout)
    else:
        print(res.stderr)
        raise
    '''


def make_global_header(filename, ref_sky_coord, xref, yref, cd1_1, cd1_2, cd2_1, cd2_2 ):
    ra=ref_sky_coord.ra.degree
    dec=ref_sky_coord.dec.degree
    
    string="WCSAXES =          2 / no comment \nCTYPE1  = 'RA---TAN' / TAN (gnomic) projection \nCTYPE2  = 'DEC--TAN' / TAN (gnomic) projection \nEQUINOX =  2000.0    / Equatorial coordinates definition (yr) \nCRVAL1  =  {ra} / RA  of reference point \nCRVAL2  =  {dec} / DEC of reference point \nCRPIX1  =  {xref}    / X reference pixel \nCRPIX2  =  {yref}   / Y reference pixel \nCUNIT1  = 'deg' / X pixel scale units \nCUNIT2  = 'deg' / Y pixel scale units \nCD1_1   =   {cd1_1} / Transformation matrix \nCD1_2   =   {cd1_2}      / no comment \nCD2_1   =   {cd2_1}      / no comment \nCD2_2   =   {cd2_2} / no comment".format(ra="%.6f"%(ra), dec="%.6f"%(dec), xref="%i"%(xref), yref="%i"%(yref), cd1_1="%.10f"%(cd1_1), cd1_2="%.10f"%(cd1_2), cd2_1="%.10f"%(cd2_1), cd2_2="%.10f"%(cd2_2) )
    with open(filename, 'w') as f:
        f.write(string)
