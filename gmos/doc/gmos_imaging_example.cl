# Copyright(c) 2002-2015 Association of Universities for Research in Astronomy, Inc.
#
# GMOS example reductions script: Typical reduction of imaging data
#
# This data processing is illustrative of data quality assessment performed
# by the Gemini staff. The processing is not designed to give the best
# possible result. Better signal-to-noise and better cleaning for cosmic-ray
# hits and bad pixels will most likely be possible. The user of these data
# is encouraged to use the provided co-added images only as guide lines and
# to re-reduce the data to obtain the best possible reduction.
#
# Gemini GMOS data reduction script
# Observation UT date:  2014jul10
# Data processor:  Inger Jorgensen / James Turner
# Data reduction date:  2015jan22
# Revisions: 2002aug26 garith changed to gemarith
#            2004oct11 changed to use gifringe and girmfringe
#            2014jan22 use new GS Hamamatsu CCD data; enable QE correction
#
# Brief data description: Imaging for the SV program GS-2014B-SV-153

# # Uncomment to clean everything up:
# imdel *.fits,*.pl
# dele *.log,*_pos,*_cen,*_trn,*.dat,tmp*,*_mag,*_coo,*_see

gmos.logfile="GS-2014B-SV-153.log"

# Location of calibration files:
set raw=/path/to/2014jul10/
set cal=/path/to/2014jul10/Basecalib/

delete glist ver-
gemlist S20140710S 34-39 > "glist"

# The bias and flat used here have been overscan- & QE-corrected

# Do the basic reduction. The flag fl_qecorr+ is only for Hamamatsu data, not
# for GMOS-N (and its use is optional, since QE differences will otherwise be
# removed by the flat). Note that Hamamatsu flats from 2014B (& possibly later)
# in the Gemini archive are not QE corrected, so you would need to re-generate
# them with giflat to follow this example with fl_qecorr+.
gireduce @glist bias=cal$S20140712S0688_bias flat1=cal$S20140720S0152_flat \
    rawpath="raw$" fl_over+ fl_qecorr-

delete rglist ver-
sections rg//@glist > rglist

# Make a fringe frame:
fields rglist 1 lines="1" | scan(s1)
s1=s1//"_frg"
gifringe("@rglist",s1)

# Subtract off the fringe frame:
girmfringe("@rglist",s1,scale=1.0)
delete frglist ver-
sections f//@rglist > frglist

# Mosaic the frames:
gmosaic @frglist 

# Coadd the frames; the default bad pixel mask is found automatically:
delete mrglist ver-
sections m//@frglist > mrglist
imcoadd @mrglist fwhm=4 threshold=100 logfile="GS-2014B-SV-153.log" 

# Measure image quality & update the headers with quality assessment info.:
gemseeing mfrgS20140710S0034_add fl_inter+ fl_over+ \
    logfile=GS-2014B-SV-153.log 
gemqa  mfrgS20140710S0034_add pass pass pass pass logfile=GS-2014B-SV-153.log

