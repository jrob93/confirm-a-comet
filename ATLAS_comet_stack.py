'''
Stack comet diff images. Creates a fits file that can be used to measure the extent of the coma
'''
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
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

path="./test_data/atlas_T05-T2394163-diff"
file_list=glob.glob("{}/*.diff".format(path))
file_list.sort()
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

    if i==ref_file:
        orig_hdr=hdr.copy()

    # remove the extra keys from the ATLAS file
    key_keep_list=list(hdr_tpv.keys())
    del_list=[]
    for key in list(hdr.keys()):
        if key not in key_keep_list:
            del_list.append(key)

    # drop duplicate entries, e.g. DATE
    del_list=list(dict.fromkeys(del_list))

    # delete extra entries
    for key in del_list:
        # print(key,hdr[key])
        del hdr[key]

    # fix the warning
    hdr['RADESYSa']=(hdr['RADECSYS'],hdr.comments['RADECSYS'])
    del hdr['RADECSYS']

    # print(repr(hdr))
    print(i,fi,hdr['DATE-OBS'])
    if i==ref_file:
        target_wcs=WCS(hdr)
        target_hdr=hdr

    img_list.append(img)
    hdr_list.append(hdr)

    # print(repr(hdu[0].header))
    hdu.writeto('test_results/{}_diff_fixed.fits'.format(fi.split("/")[-1].split(".diff")[0]),overwrite=True)

    hdu.close()

fixed_file_list=glob.glob("test_results/*_fixed.fits")
print(fixed_file_list)

img_CCD_list=[]
# load the fixed fits files as CCDdata objects
for i,fi in enumerate(fixed_file_list):
    img_CCD_list.append(CCDData.read(fi, unit=u.dimensionless_unscaled))
    print("\n{}:".format(fi))
    print(img_CCD_list[-1].wcs)

# stack_img=1
#
# if stack_img==1:
combiner = Combiner(img_CCD_list)
combiner.sigma_clipping(low_thresh=2, high_thresh=2, func=np.ma.median)
# stacked_image = combiner.average_combine()
stacked_image = combiner.median_combine()

# stacked_image.wcs=target_wcs
stacked_image.header=orig_hdr
stacked_image.header['COMMENT']="this image was stacked with {}".format(os.path.basename(__file__))
print("\n\n")
print(stacked_image.wcs)
# print(repr(stacked_image.header))

# stacked_image.write('median_combiner.fits',overwrite=True,output_verify='fix')
img=stacked_image.data

# else:
#     ccd_stack=CCDData.read('median_combiner.fits', unit=u.dimensionless_unscaled)
#     img=ccd_stack.data

#set up the normalisation scheme for the image
norm = aviz.ImageNormalize(img,interval=aviz.ZScaleInterval())

fig = pyplot.figure()
gs = gridspec.GridSpec(1, 1)
ax1 = pyplot.subplot(gs[0,0], projection = target_wcs)
# ax1 = pyplot.subplot(gs[0,0])

#display the image, note that you need to set the origin
s1=ax1.imshow(img, norm=norm, origin='lower')

cbar1=fig.colorbar(s1)

outfile = 'test_results/comet_median_stack.fits'
fits.writeto(outfile,stacked_image.data,stacked_image.header,overwrite=True,output_verify='fix')
print("save stacked file: {}".format(outfile))

pyplot.show()
