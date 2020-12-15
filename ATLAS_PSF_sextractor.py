'''
Run sextractor on a reduced ATLAS image and plot the resulting FWHMs onto the candidate comet and background stars.
'''

import glob
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import subprocess
import os

from astropy.io import fits
from astropy import visualization as aviz
from astropy.io import ascii
from photutils import aperture_photometry, CircularAperture
from astropy.wcs import WCS

tracklet="atlas_T08-T3753359"

# do PSF on stars or comet
# stack_on="comet"
stack_on="stars"

hist_dist=1 # when looking at stars do a hist of fwhm

comet_pos_range=3 # radius for central object position (pixels)
pix_scale=1.86 # arcsec
mag_range = 1000 # select stars within +/- mag_range of the comet magnitude

sex_settings="./{}/{}_{}.sex".format(tracklet,tracklet,stack_on) # the sex file is set up in ATLAS_tracklet_data.py
sex_save_file="{}/{}_{}.cat".format(tracklet,tracklet,stack_on) # change this line in default.sex!

# # find all the reduced data files
# data_extension=".fits.fz"
# path="/Users/jrobinson/confirm-a-comet/{}/{}-reduced".format(tracklet,tracklet)
# file_list=glob.glob("{}/*{}".format(path,data_extension))
# file_list.sort()
# print(file_list)

# load a stacked file
file_list=["/Users/jrobinson/confirm-a-comet/{}/{}_median_stack_{}.fits".format(tracklet,tracklet,stack_on)]

for i,fi in enumerate(file_list):

    # if i!=2:
    #     continue

    #load fits file
    hdu = fits.open(fi)

    #define the image and header data
    img = hdu[0].data
    hdr=hdu[0].header

    # print(repr(hdr))

    # break

    fi_name=fi.split("/")[-1]

# we only run the sextractor command if the results file does not already exist - NOTE NEED TO RUN TWICE FOR COMET AND STAR STACK!
if not os.path.isfile(sex_save_file):
    # run the sextractor command here:
    # sex 02a58833o0083o.fits.fz -c extra_files/sextractor_files/default.set_xlim
    sex_cmd='source-extractor {} -c {}'.format(fi,sex_settings)
    print(sex_cmd)
    p=subprocess.Popen(sex_cmd,shell=True).wait()

else:
    # load the sextractor data
    data=ascii.read(sex_save_file)
    print("{}: results loaded".format(sex_save_file))

# # run the sextractor command here:
# sex_cmd='source-extractor {} -c {}'.format(fi,sex_settings)
# print(sex_cmd)
# p=subprocess.Popen(sex_cmd,shell=True).wait()
#
# # load the sextractor data
# data=ascii.read(sex_save_file)
# print("{}: results loaded".format(sex_save_file))

print(data)
print(len(data))

print(len(img))
print(len(img)/2)

# create mask for the central object
mask1=((data["X_IMAGE"]>((len(img)/2)-comet_pos_range)) & (data["X_IMAGE"]<((len(img)/2)+comet_pos_range)))
mask2=((data["Y_IMAGE"]>((len(img)/2)-comet_pos_range)) & (data["Y_IMAGE"]<((len(img)/2)+comet_pos_range)))
mask_comet=(mask1 & mask2)
mask_star=~mask_comet

if "comet" in fi_name:
    # show only the comet at image centre
    data=data[mask_comet]

if "star" in fi_name:
    # load the comet psf from the other fit to get the comet fwhm and mag
    try:
        sex_save_file_comet="{}/{}_comet.cat".format(tracklet,tracklet)
        data_comet=ascii.read(sex_save_file_comet)
        print("{}: results loaded".format(sex_save_file_comet))
    except:
        print("need to run PSF on image stacked on comet")
    mask_comet1=((data_comet["X_IMAGE"]>((len(img)/2)-comet_pos_range)) & (data_comet["X_IMAGE"]<((len(img)/2)+comet_pos_range)))
    mask_comet2=((data_comet["Y_IMAGE"]>((len(img)/2)-comet_pos_range)) & (data_comet["Y_IMAGE"]<((len(img)/2)+comet_pos_range)))
    data_comet=data_comet[mask_comet1 & mask_comet2]
    print(data_comet)
    comet_fwhm=float(data_comet["FWHM_IMAGE"]*pix_scale)
    comet_mag=float(data_comet["MAG_PSF"])
    print("comet fwhm = {}\", comet mag = {}".format(comet_fwhm,comet_mag))

    # exclude the central object
    data=data[mask_star]

    # show only stars with similar brightness
    data=data[(data["MAG_PSF"]<comet_mag+mag_range) & (data["MAG_PSF"]>comet_mag-mag_range)]

print(data)
print(len(data))

# sigma clip the FWHM distribution
fwhm=np.array(data["FWHM_IMAGE"])
std=np.std(fwhm)
med=np.median(fwhm)
low=2
high=2
clip_mask=((fwhm < (med-(std*low))) | (fwhm > (med+(std*high))))
data=data[~clip_mask]
print(data)

# extract the positions of all sources
positions = np.transpose((data["X_IMAGE"], data["Y_IMAGE"]))

# sextractor doesn't index from zero?
positions=positions-1

# define image normalisation
norm = aviz.ImageNormalize(img,interval=aviz.ZScaleInterval())

fig = pyplot.figure()

if (hist_dist==1) and ("star" in fi_name):

    gs = gridspec.GridSpec(2, 1)
    ax2 = pyplot.subplot(gs[1,0])

    #plot the distribution of FWHM
    ax2.hist(data["FWHM_IMAGE"]*pix_scale,bins="auto")
    ax2.axvline(np.median(data["FWHM_IMAGE"]*pix_scale),label="median = {:.3f}".format(np.median(data["FWHM_IMAGE"]*pix_scale)), c='r')
    ax2.axvline(np.mean(data["FWHM_IMAGE"]*pix_scale),label="mean = {:.3f}".format(np.mean(data["FWHM_IMAGE"]*pix_scale)), c='r', ls="--")
    ax2.axvline(comet_fwhm,label="comet = {:.3f}".format(comet_fwhm), c='r', ls=":")
    ax2.legend()
    ax2.set_xlabel("FWHM(\")")
    ax2.set_ylabel("N")
else:
    gs = gridspec.GridSpec(1, 1)

if "stars" in fi_name or "comet" in fi_name:
    target_wcs=WCS(hdr)
    ax1 = pyplot.subplot(gs[0,0], projection = target_wcs)
    ax1.set_xlabel("RA")
    ax1.set_ylabel("DEC")
else:
    ax1 = pyplot.subplot(gs[0,0])
    ax1.set_xlabel("x pixels")
    ax1.set_ylabel("y pixels")

# stop overlays messing with zoom
# ax1.set_xlim(0,len(img[0,:]))
# ax1.set_ylim(0,len(img[:,0]))
mid_x=len(img[0,:])/2
mid_y=len(img[:,0])/2
# x_range=100
# y_range=100
x_range=75
y_range=75
ax1.set_xlim(mid_x-x_range,mid_x+x_range)
ax1.set_ylim(mid_y-y_range,mid_y+y_range)

#display the image, note that you need to set the origin
s1=ax1.imshow(img, norm=norm, origin='lower')

# plot detected sources
for i in range(len(positions)):
    if data["FWHM_IMAGE"][i]<=0:
        continue
    # print(positions[i], data["FWHM_IMAGE"][i])
    aperture = CircularAperture(positions[i], r=data["FWHM_IMAGE"][i]/2)
    aperture.plot(color='red', lw=1.5, alpha=0.5)
    # ax1.scatter(positions[i][0],positions[i][1],marker="+",color="r")

    if hist_dist==1:
        if ((positions[i][0]>(mid_x-x_range)) and (positions[i][0]<(mid_x+x_range))) and ((positions[i][1]>(mid_y-y_range)) and (positions[i][1]<(mid_y+y_range))):
            # ax1.text(positions[i][0],positions[i][1],"{:.2f}''".format(data["FWHM_IMAGE"][i]*pix_scale),color="r",clip_on=True)
            ax1.text(positions[i][0],positions[i][1],"{:.2f},{:.2f}".format(data["FWHM_IMAGE"][i]*pix_scale,data["MAG_PSF"][i]),color="r",clip_on=True)
    else:
        # ax1.text(positions[i][0],positions[i][1],"{:.2f}''".format(data["FWHM_IMAGE"][i]*pix_scale),color="r",clip_on=True)
        ax1.text(positions[i][0],positions[i][1],"{:.2f},{:.2f}".format(data["FWHM_IMAGE"][i]*pix_scale,data["MAG_PSF"][i]),color="r",clip_on=True)

    # ax1.text(positions[i][0],positions[i][1],"{:.2f},{:.2f}".format(data["FWHM_IMAGE"][i]*pix_scale,data["MAG_PSF"][i]),color="r",clip_on=True)

# add a central rectangle
if "stars" in fi_name:
    x_len=10
else:
    x_len=30
y_len=x_len
r = Rectangle((mid_x-x_len/2, mid_y-y_len/2), x_len, y_len, edgecolor='yellow', facecolor='none')
ax1.add_patch(r)

cbar1=fig.colorbar(s1)
cbar1.set_label("counts")

pyplot.title(fi.split("/")[-1])
# pyplot.tight_layout()

fname="{}/{}_{}.png".format(tracklet,os.path.basename(__file__).split('.')[0],fi.split("/")[-1])
print(fname)
pyplot.savefig(fname, bbox_inches='tight')

pyplot.show()
