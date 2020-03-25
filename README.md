# confirm-a-comet

ATLAS_comet_stack.py: Astropy and related packages are used to do a median combiner of the diff tracklet data.
With the stacked image, use ds9 to measure length and direction of tail. This can be done by:
- edit > region
- region > shape > vector

Draw this onto the coma and obtain the size and direction for the MPC new comet report.

ATLAS_PSF_sextractor.py: To confirm cometary activity one usually compares the PSF of the detected object to the comae of background stars. 
Here we use the source-extractor, which should be installed on your machine as the sex command.
