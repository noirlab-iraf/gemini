# Copyright(c) 2009 Association of Universities for Research in Astronomy, Inc.
#
# GMOS example reductions script: Typical reduction of N&S longslit data
#
# This data processing was done to aid the data quality assessment performed
# by the Gemini staff. The data processing is not designed to give the best
# possible result. Better signal-to-noise and better cleaning for cosmic-ray
# hits and bad pixels will most likely be possible. The user of these data
# is encouraged to use the provided co-added images only as guide lines and
# to re-reduce the data to obtain the best possible reduction.
#
# Gemini GMOS data reduction script
# Observation UT date: 2005April06 
# Data processor:   Jen Holt, Inger Jorgensen, Kathy Roth
#
# Brief data description: Nod and shuffle 3 amp longslit GMOS-N
#

for (i=57;i<=59;i+=1) {
    print("N20050111S00"//i//".fits", >> "dark.list")
}
print("N20050406S0071.fits", >> "flat.list")
print("N20050406S0075.fits", >> "flat.list")

for (i=72;i<=74;i+=1) {
    print("N20050406S00"//i//".fits", >> "sci.list")
}

print("N20050508S0016.fits", >> "arc.list")

# You always want to overscan correct everything.  If the overscan level 
# is different on the science compared to the dark and you don't overscan 
# correct the dark correction might not work very well. 
gndark @dark.list gnsdark.fits rawpath="rawdir$" fl_vardq+ fl_over+

#   If there is a Nod&Shuffle dark taken close in time to the science
# (within a few months should be fine) There is not a need to do a bias 
# subtraction, because you would have to bias correct both the dark and 
# the science and so the bias correction effectively subtracts out when 
# you dark correct the science.  The exception to this might be if you 
# are using a very very old dark, in which case you might wish to bias 
# correct the dark using a bias frame constructed from biases taken close
# in time to when the dark was taken and similarly for the science. 

# Find number for "nshuffle" distance in header keyword: NODPIX note:
# NODPIX will be in the science header, not the flat field header.

gsflat @flat.list gnsflat.fits fl_double+ rawpath="rawdir$" fl_over+ \
    fl_trim+ fl_vardq+ nshuffle=38 fl_bias- fl_dark- mdfdir="rawdir$" \
    fl_detec+ bias="" fl_usegrad+ \
    fl_fixpix- combflat="flat_comb" fl_keep+ verb+

gsreduce g@sci.list fl_dark+ fl_flat+ fl_bias- flat="gnsflat" \
    dark="gnsdark" fl_over+ fl_trim+ fl_vardq+ bias="" \
    rawpath="rawdir$" mdfdir="rawdir$" fl_gscrrej+ fl_cut- fl_gmosaic- \
    fl_fixpix- fl_gsappwave-

sections gs//@sci.list > gssci.list
gnscombine "gssci.list" outimage="scicomb.fits" offsets="rawdir$offset_list"

gmosaic flat_comb
gscut mflat_comb fl_update+ gradimage=mflat_comb

gsreduce scicomb.fits fl_cut+ fl_gmosaic+ fl_fixpix+ fl_gsappwave+ \
    refimage="mflat_comb.fits" fl_flat- fl_bias- fl_trim- fl_over-

gsreduce @arc.list fl_dark- fl_flat- fl_bias- fl_over+  fl_vardq+ \
    bias="" mdfdir="rawdir$" rawpath="rawdir$" fl_trim+

gswavelength gsN20050508S0016 fl_inter-
gstransform gsN20050508S0016 wavtran=gsN20050508S0016
gstransform gsscicomb wavetran=gsN20050508S0016

gsextract tgsscicomb.fits
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
