'''
Stack comet diff images. Creates a fits file that can be used to measure the extent of the coma
'''
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import glob
import os

from astropy.io import fits
from astropy import visualization as aviz
from astropy.wcs import WCS
from astropy.wcs import find_all_wcs

from ccdproc import Combiner
from ccdproc import wcs_project
from astropy.nddata import CCDData
from astropy import units as u

tracklet="atlas_T08-T3753359"

# choose which data to stack
# fits_type="diff"
fits_type="reduced"

# choose which object to stack on
stack_on="stars"
# stack_on="comet"

path="/Users/jrobinson/confirm-a-comet/{}/{}-{}".format(tracklet,tracklet,fits_type)
print(path)
if fits_type=="reduced":
    file_list=glob.glob("{}/*.fits.fz".format(path))
if fits_type=="diff":
    file_list=glob.glob("{}/*.diff.fz".format(path))
save_path=tracklet
stack_save_file="{}_median_stack_{}.fits".format(tracklet,stack_on)

file_list.sort()
# file_list=np.delete(file_list,2) # drop a specific frame (e.g. comet close to star)
file_list=file_list[::-1] # reverse order
print(file_list)

img_list=[]
hdr_list=[]

# load reference tpv fits example file
hdu_tpv = fits.open("extra_files/tpv.fits")
img_tpv = hdu_tpv[0].data
hdr_tpv = hdu_tpv[0].header
hdu_tpv.close()

ref_file=0

for i,fi in enumerate(file_list):

    #load fits file
    hdu = fits.open(fi)
    # print(hdu.info())
    # print(len(hdu))

    #define the image and header data
    img = hdu[0].data
    hdr=hdu[0].header

    # print(repr(hdr))
    for key in ['DATE-OBS','OBSNAME','CRVAL1','CRVAL2','CRPIX1','CRPIX2']:
        print("{}\t{}\t{}".format(key,hdr[key],hdr.comments[key]))

    # print(hdr['CTYPE1'],hdr.comments['CTYPE1'])
    # print(hdr['CTYPE2'],hdr.comments['CTYPE2'])

    if i==ref_file:
        orig_hdr=hdr.copy()

    # remove the extra keys from the ATLAS file
    key_keep_list=list(hdr_tpv.keys())
    del_list=[]
    key_list=[]
    for key in list(hdr.keys()):
        if key not in key_keep_list:
            del_list.append(key)
        else:
            key_list.append(key)

    # drop duplicate entries, e.g. DATE
    del_list=list(dict.fromkeys(del_list))
    # print(del_list)
    del_list=[d for d in del_list if not d.startswith("PV")] # keep higher order PV terms
    # print(del_list)

    # delete extra entries
    for key in del_list:
        # print(key,hdr[key])
        del hdr[key]

    # fix the warning
    hdr['RADESYSa']=(hdr['RADECSYS'],hdr.comments['RADECSYS'])
    del hdr['RADECSYS']

    # print(repr(hdr))

    if i==ref_file:
        target_wcs=WCS(hdr)
        target_hdr=hdr
    print(target_wcs)

    hdu[0].data=img
    hdu[0].header=hdr

    # print(repr(hdu[0].header))

    hdu.writeto('{}/{}_{}_fixed.fits'.format(save_path,fi.split("/")[-1].split(".diff")[0],fits_type),overwrite=True)

    hdu.close()

fixed_file_list=glob.glob("{}/*{}_fixed.fits".format(save_path,fits_type))
# fixed_file_list=glob.glob("{}/{}/*.fits".format(save_path,tracklet))
# fixed_file_list=fixed_file_list[:-1]
print(fixed_file_list)

print("\nload images\n")
img_CCD_list=[]
# load the fixed fits files as CCDdata objects
for i,fi in enumerate(fixed_file_list):
    print("\n{}:".format(fi))
    img=CCDData.read(fi, unit=u.dimensionless_unscaled)
    # print(repr(img.header))
    # print(img.header[''])
    if stack_on=="stars":
        img = wcs_project(img, target_wcs) # transform to stack on background stars
    img_CCD_list.append(img)
    print(img_CCD_list[i].wcs)
    print()

print("\nstack images\n")
combiner = Combiner(img_CCD_list)
combiner.sigma_clipping(low_thresh=2, high_thresh=2, func=np.ma.median)
# scaling_func = lambda arr: 1/np.ma.average(arr)
# combiner.scaling = scaling_func
# stacked_image = combiner.average_combine()
stacked_image = combiner.median_combine()
# stacked_image = combiner.sum_combine()

# stacked_image.wcs=target_wcs
# stacked_image.header=orig_hdr
stacked_image.header=target_hdr
stacked_image.header['COMMENT']="this image was stacked with {} on".format(os.path.basename(__file__),stack_on)
print("\n\n")
print(stacked_image.wcs)

# stacked_image.write('median_combiner.fits',overwrite=True,output_verify='fix')
img=stacked_image.data

# # subtract the image median?
# print(' Mean = {:.2f}, Median = {:.2f}, Std. Dev.= {:.2f} '.format(np.mean(img), np.median(img), np.std(img)))
# img=img-np.median(img)

#set up the normalisation scheme for the image
norm = aviz.ImageNormalize(img,interval=aviz.ZScaleInterval())

fig = pyplot.figure()
gs = gridspec.GridSpec(1, 1)
ax1 = pyplot.subplot(gs[0,0], projection = target_wcs)
# ax1 = pyplot.subplot(gs[0,0])

#display the image, note that you need to set the origin
s1=ax1.imshow(img, norm=norm, origin='lower')

# draw a box around the centre
r_len=25
x_cent=len(img[0,:])/2
y_cent=len(img[:,0])/2
r = Rectangle((x_cent-r_len, y_cent-r_len), 2*r_len, 2*r_len, edgecolor='yellow', facecolor='none')
ax1.add_patch(r)

cbar1=fig.colorbar(s1)
ax1.set_xlabel("RA")
ax1.set_ylabel("DEC")

cbar1.set_label("counts")

# display as pythion based visulaiser

outfile = '{}/{}'.format(save_path,stack_save_file)
fits.writeto(outfile,stacked_image.data,stacked_image.header,overwrite=True,output_verify='fix')
print("save stacked file: {}".format(outfile))

pyplot.show()
