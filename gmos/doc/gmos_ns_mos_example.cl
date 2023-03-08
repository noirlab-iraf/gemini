# Copyright(c) 2002-2004 Association of Universities for Research in Astronomy, Inc.
#
# GMOS example reductions script: Typical reduction of N&S MOS data
#
# This data processing was done to aid the data quality assessment performed
# by the Gemini staff. The data processing is not designed to give the best
# possible result. Better signal-to-noise and better cleaning for cosmic-ray
# hits and bad pixels will most likely be possible. The user of these data
# is encouraged to use the provided co-added images only as guide lines and
# to re-reduce the data to obtain the best possible reduction.
#
# Gemini GMOS data reduction script
# Observation UT date:  2002aug10
# Data processor:  Dione Scheltus, Inger Jorgensen
# Data reduction date: 2003jun03
#
# Brief data description: N&S MOS data for the SV program GN-2002A-SV-78

gmos.logfile="GN-2002A-SV-78.log"

# Define various directories containing raw data and calibration data
set raw="../gmos_testdata/"
set gcalib="../gmos_testdata/"

gsreduce N20020810S0117,N20020810S0118,N20020810S0123,N20020810S0124,\
N20020810S0125,N20020810S0126 bias=gcalib$N20020803S0025_bias.fits \
fl_over- fl_flat- fl_gmosaic- fl_gsappwave- fl_cut- fl_title- \
mdffile=GN2002A-SV-78-7.fits rawpath=raw$ mdfdir=raw$

# gnscombine has to call gnsskysub otherwise this does not work
# Using DTAX offsets from the observing log
#  Offset in y in pixels = int( -dDTAX / 13.5 )
delete offsets.dat ver-
printf("0 0\n0 -2\n0 -4\n0 2\n0 4\n0 -4\n", >> "offsets.dat")

gnscombine gsN20020810S0117,gsN20020810S0118,gsN20020810S0123,gsN20020810S0124,\
gsN20020810S0125,gsN20020810S0126 offsets=offsets.dat outim=N20020810S0117_comb \
outch=N20020810S0117_combcr outmed=N20020810S0117_sky

gsreduce N20020810S0117_comb.fits outpref=r fl_fixpix- fl_trim- fl_bias- fl_flat- \
fl_gsappwave- fl_cut- fl_title- geointer=nearest

# After this part of the reduction, the data should be flat fielded, cut into
# slitlets and wavelength calibrated. The spectra may then be extracted.
