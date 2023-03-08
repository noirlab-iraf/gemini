# Copyright(c) 2002-2009 Association of Universities for Research in Astronomy, Inc.
#
# GMOS example reductions script: Typical reduction of science observation 
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
# Brief data description: Observations of NGC1068, GMOS-N commissioning data
#   IFU in 2-slit mode
# The default parameters for the IFU tasks are set to reduce 2-slit data
# If reducing 1-slit data, check the help files to see which parameters need
# modification
# 
set mdata=/net/sabrina/staging2/gmos/
set gcalib=/usr/dataproc/gmos/2001sep10/Basecalib/

# Set the logfile
gmos.logfile="GN-2001B-SV-110_example.log"

# Set common parameter values
gfreduce.rawpath="mdata$2001sep09/"
gfreduce.bias="gcalib$N20010806S018_bias" # bias is overscan subtracted
gfreduce.fl_fluxcal=no
gfreduce.xoffset=-18.0

# GCAL flat
gfreduce N20010908S105 fl_gscrrej- fl_wavtran- fl_skysub- fl_inter- fl_over+

# Inspect the flat field interactively
gfdisplay ergN20010908S105

# Twilight flat
gfreduce N20010908S112 fl_wavtran- fl_skysub- fl_inter- trace- \
ref=ergN20010908S105 fl_gscrrej- fl_over+ biasrows="3:64"

# Inspect the flat field interactively
gfdisplay ergN20010908S112 

# Make response curves with twilight correction
gfresponse ergN20010908S105 ergN20010908S105_resp112 sky=ergN20010908S112 \
order=95 fl_inter- func=spline3 sample="*"

# Arc - overscan subtract instead of using the bias, since the bias does not apply 
# to fast read
gfreduce N20010908S108.fits fl_wavtran- fl_inter- ref=ergN20010908S105 \
recenter- trace- fl_skysub- fl_gscrrej- fl_bias- fl_over+ order=1 \
weights=none fl_over+ biasrows="3:64"

# Establish the wavelength calibration
gswavelength ergN20010908S108 fl_inter+ nlost=10 

# Reduce the science data
gfreduce N20010908S101 fl_inter- verb+  refer=ergN20010908S105 recenter- 
trace- fl_wavtran+ wavtran=ergN20010908S108 response=ergN20010908S105_resp112 
fl_over+ biasrows="3:64"

# Inspect the reduced data, keep the 2D white light image
gfdisplay stexrgN20010908S101 output="stexrgN20010908S101_2D"  ver=1

# Make a datacube of the reduced spectra
gfcube stexrgN20010908S101

# Standard calibrate the reduced data - requires GN-CAL20010908_example.cl run first
# Apply calibration to extracted spectra
gscalibrate stexrgN20010908S101.fits 

# Make a datacube of the standard calibrated spectra
gfcube cstexrgN20010908S101.fits 

