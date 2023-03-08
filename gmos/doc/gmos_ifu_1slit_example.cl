# Copyright(c) 2002-2004 Association of Universities for Research in Astronomy, Inc.
#
# GMOS example reductions script: Typical reduction of science observation
#                                 with the IFU in 1-slit mode
#
# This data processing was done to aid the data quality assessment performed
# by the Gemini staff. The data processing is not designed to give the best
# possible result. Better signal-to-noise and better cleaning for cosmic-ray
# hits and bad pixels will most likely be possible. The user of these data
# is encouraged to use the provided co-added images only as guide lines and
# to re-reduce the data to obtain the best possible reduction.
#
# Gemini GMOS data reduction script
# Observation UT date: 2002aug13 
# Data processor:   Inger Jorgensen
# Data reduction date: 2002aug28
#
# Brief data description: GN-2002B-Q-80 science target
#   IFU in 1-slit mode, data binned 2x1
# 
set gdata=/net/sabrina/staging4/gmos/
set gcalib=/usr/dataproc/gmos/2002aug13/Basecalib/

# Set the logfile
gmos.logfile="GN-2002B-Q-80_example.log"

# Set common parameter values
gfreduce.rawpath="mdata$2002aug07/"
gfreduce.bias="gcalib$N20020803S0025_bias"
gfreduce.fl_fluxcal=no
gfreduce.xoffset=-6.1
gfreduce.xoffset=INDEF

# GCAL flat
gfreduce N20020807S0043.fits fl_gscrrej- fl_wavtran- fl_skysub- fl_inter- slits=red

# Inspect the flat field interactively
gfdisplay ergN20020807S0043 ver=1 \
ver=1 config=gmos$data/gnifu_slitr.cfg deadfib=gmos$data/gnifu_deadr.cfg

# Twilight flat
gfreduce N20020807S0085.fits fl_wavtran- fl_skysub- fl_inter- trace- \
ref=ergN20020807S0043  fl_gscrrej- slits=red

# Inspect the flat field interactively
gfdisplay ergN20020807S0085 ver=1 \
ver=1 config=gmos$data/gnifu_slitr.cfg deadfib=gmos$data/gnifu_deadr.cfg

# Make response curves with twilight correction
gfresponse ergN20020807S0043.fits ergN20020807S0043_resp085.fits  \
sky=ergN20020807S0085 order=95 fl_inter- func=spline3 sample="*"

# Arc - overscan subtract instead of using the bias, since the bias does not 
# apply to fast read
gfreduce N20020807S0038.fits fl_wavtran- fl_inter- ref=ergN20020807S0043 recenter- \
trace- fl_skysub- fl_gscrrej- fl_bias- fl_over+ order=1 weights=none slits=red

# Establish the wavelength calibration
gswavelength ergN20020807S0038 fl_inter+ nlost=10  ntarget=15 \
coordlist="gmos$data/GCALcuar.dat" aiddebug=s threshold=25  section="first line"

# science target
# This observatons has an object in the sky IFU, so the sky subtraction would be 
# quite bad with the default value of "expr". "expr" is set such that only the 
# upper part # of the small IFU field is used to define the sky spectrum
gfreduce N20020807S0040 fl_inter- verb+  refer=ergN20020807S0043 recenter- trace- \
fl_wavtran+ wavtran=ergN20020807S0038 response=ergN20020807S0043_resp085.fits \
fl_gscrrej+ slits=red expr="XINST>10 && YINST>2"

# Inspect the reduced spectra interactively
gfdisplay stexrgN20020807S0040 out=stergN20020807S0040_2D \
ver=1 config=gmos$data/gnifu_slitr.cfg deadfib=gmos$data/gnifu_deadr.cfg

# Resample the reduced spectra to a data cube
gfcube stexrgN20020807S0040.fits outimage=stergN20020807S0040_3D.fits

# Standard calibration  - this requires reduction of a spectrophotometric
# standard star as outlined in teh IFU 2-slit examples, see
#    gmosexamples IFU-specstd
#    gmosexamples IFU-2
