# Copyright(c) 2002-2004 Association of Universities for Research in Astronomy, Inc.
#
# GMOS example reductions script: Typical reduction of standard star observation 
#                                 with the IFU in 2-slit mode
#
# This data processing was done to aid the data quality assessment performed
# by the Gemini staff. The data processing is not designed to give the best
# possible result. Better signal-to-noise and better cleaning for cosmic-ray
# hits and bad pixels will most likely be possible. The user of these data
# is encouraged to use the provided co-added images only as guide lines and
# to re-reduce the data to obtain the best possible reduction.
#
# Gemini GMOS data reduction script
# Observation UT date: 2001sep09 
# Data processor:   Inger Jorgensen
# Data reduction date: 2002sep05 
#
# Brief data description: Observations of BD284211, GMOS-N commissioning data
#   IFU in 2-slit mode
# The default parameters for the IFU tasks are set to reduce 2-slit data
# If reducing 1-slit data, check the help files to see which parameters need
# modification
# 
set mdata=/staging2/gmos/
set gcalib=/usr/dataproc/gmos/2001sep10/Basecalib/

# Set the logfile
gmos.logfile="GN-CAL20010909_example.log"

# Set common parameter values
gfreduce.rawpath="mdata$2001sep09/"
gfreduce.bias="gcalib$N20010806S018_bias"
gfreduce.fl_fluxcal=no
gfreduce.xoffset=-18.0

# GCAL flat
gfreduce N20010908S092 fl_gscrrej- fl_wavtran- fl_skysub- fl_inter-

# Twilight flat - also reduced in GN-2001B-SV-110_example.cl
gfreduce N20010908S112 fl_wavtran- fl_skysub- fl_inter- trace- ref=ergN20010908S092 \
fl_gscrrej-

# Make response curves with twilight correction
gfresponse ergN20010908S092 ergN20010908S092_resp112 sky=ergN20010908S112 \
order=95 fl_inter- func=spline3 sample="*"

# Arc - overscan subtract instead of using the bias, since the bias does not 
# apply to fast read
gfreduce N20010908S093.fits fl_wavtran- fl_inter- ref=ergN20010908S092 \
recenter- trace- fl_skysub- fl_gscrrej- fl_bias- fl_over+ order=1 weights=none

# Establish the wavelength calibration
gswavelength ergN20010908S093 fl_inter+ nlost=10 

# Reduce the standard star observation
gfreduce N20010908S091 fl_inter- verb+  refer=ergN20010908S092 recenter- trace- \
fl_wavtran+ wavtran=ergN20010908S093 response=ergN20010908S092_resp112

# Inspect the reduced data interactively, keep the reconstructed white light image as 
# stexrgN20010908S091_2D
# The input image as one image extension after the processing with gfreduce. 
# Therefore gfdisplay has to be used with ver=1
gfdisplay stexrgN20010908S091 output="stexrgN20010908S091_2D"  ver=1

# Sum all the standard star spectra
gfapsum stexrgN20010908S091.fits fl_inter-

# Determine sensitivity function
gsstandard astexrgN20010908S091.fits std sens starname=bd284211 \
caldir=onedstds$spec50cal/ fl_inter-

# Apply flux calibration calibration
gscalibrate astexrgN20010908S091.fits 

