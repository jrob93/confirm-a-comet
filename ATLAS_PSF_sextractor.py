import glob
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import subprocess

from astropy.io import fits
from astropy import visualization as aviz
from astropy.io import ascii
from photutils import aperture_photometry, CircularAperture

path="./test_data/atlas_T05-T2394163-reduced/"
file_list=glob.glob("{}/*.fits.fz".format(path))
file_list.sort()
print(file_list)

for i,fi in enumerate(file_list):

    #load fits file
    hdu = fits.open(fi)

    #define the image and header data
    img = hdu[0].data
    hdr=hdu[0].header

    # print(repr(hdr))

    break

# run the sextractor command here:
# sex 02a58833o0083o.fits.fz -c extra_files/sextractor_files/default.sex
p=subprocess.Popen('source-extractor {} -c extra_files/sextractor_files/default.sex'.format(fi),shell=True)

# load the sextractor data
data=ascii.read("test_results/test.cat")
positions = np.transpose((data["X_IMAGE"], data["Y_IMAGE"]))

# sextractor doesn't index from zero?
positions=positions-1

# define image normalisation
norm = aviz.ImageNormalize(img,interval=aviz.ZScaleInterval())

fig = pyplot.figure()
gs = gridspec.GridSpec(1, 1)
# ax1 = pyplot.subplot(gs[0,0], projection = target_wcs)
ax1 = pyplot.subplot(gs[0,0])

# stop overlays messing with zoom
# ax1.set_xlim(0,len(img[0,:]))
# ax1.set_ylim(0,len(img[:,0]))
mid_x=len(img[0,:])/2
mid_y=len(img[:,0])/2
x_range=100
y_range=100
ax1.set_xlim(mid_x-x_range,mid_x+x_range)
ax1.set_ylim(mid_y-y_range,mid_y+y_range)

#display the image, note that you need to set the origin
s1=ax1.imshow(img, norm=norm, origin='lower')

# plot detected sources
for i in range(len(positions)):
    if data["FWHM_IMAGE"][i]<=0:
        continue
    # print(positions[i], data["FWHM_IMAGE"][i])
    aperture = CircularAperture(positions[i], r=data["FWHM_IMAGE"][i])
    aperture.plot(color='red', lw=1.5, alpha=0.5)
    # ax1.scatter(positions[i][0],positions[i][1],marker="+",color="r")
    ax1.text(positions[i][0],positions[i][1],data["FWHM_IMAGE"][i],color="r",clip_on=True)

# add a central rectangle
x_len=30
y_len=x_len
r = Rectangle((mid_x-x_len/2, mid_y-y_len/2), x_len, y_len, edgecolor='yellow', facecolor='none')
ax1.add_patch(r)

cbar1=fig.colorbar(s1)

pyplot.show()