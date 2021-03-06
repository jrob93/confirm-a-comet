# confirm-a-comet

ATLAS_comet_stack.py: Astropy and related packages are used to do a median combiner of the diff tracklet data.
With the stacked image, use ds9 to measure length and direction of tail. This can be done by:
- edit > region
- region > shape > vector

Draw this onto the coma and obtain the size and direction for the MPC new comet report.
ATLAS_PSF_sextractor.py: To confirm cometary activity one usually compares the PSF of the detected object to the comae of background stars.

Here we use the source-extractor, which should be installed on your machine and used via the command-line 'sex' command. This can be found on github at https://sextractor.readthedocs.io/en/latest/index.html
and the original paper (Bertin and Arnouts 1996) is at
https://aas.aanda.org/articles/aas/abs/1996/08/ds1060/ds1060.html

A very faint asteroid may appear fuzzy, due to low signal to noise.
FWHM of the candidate should be compared to FWHM of similar brightness stars.
Comets will most likely have low motion across the sky.
