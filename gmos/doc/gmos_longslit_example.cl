# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
#
# GMOS example script: Typical reduction of long-slit spectrum
#
# This example should only be used as a guide and can likely be improved on,
# eg. for better cosmetics and signal-to-noise ratio; the user is encouraged
# to optimize the sequence carefully for his/her particular dataset (see each
# task's help pages for available options).
#
# A subset of the data from GS-2015A-Q-99 have been used here, with the kind
# agreement of the PI, and will be made available after publication (by 2016).
#
# Define directory containing raw data. A separate directory containing
# pre-processed calibration data can also be defined here, if applicable.
set rawdir="raw/"

# Set up the log file for this reduction
gmos.logfile="GS-2015A-Q-99.log"

# Define path to the optional dependency LA Cosmic (see "help gemcrspec"), if
# using the new fl_crspec option (rather than fl_gscrrej) to clean cosmic rays
# in gsreduce and lacos_spec has not already been declared in loginuser.cl.
task lacos_spec=./lacos_spec.cl

# Most of the steps below have an fl_inter option that specifies whether to
# use an interactive plot for fitting, tracing, line identification etc. or
# whether to run through the process uninterrupted using defaults. This option
# should be enabled the first time each step is run, to inspect and possibly
# tweak the results, recording any parameter changes so they can be added to
# the command line for repeatibility; it can then be disabled if the step
# needs repeating subsequently (eg. after improving the cosmetics for some
# earlier step).
# 
# The non-default "fl_vardq" option that creates or propagates variance and
# data quality planes is also enabled here for most of the steps, in order to
# track random errors and cosmetics through the reduction process.

# # Uncomment to delete all the files from running the steps previously:
# delete *.log,*.lis,tmp*
# imdel *.fits,*.pl
# !rm -fr database

# # Edit & uncomment to copy the flux calibration from "gmosexamples
# standard_longslit" here, if not run in the same directory:
# copy /some/path/sens.fits .

# Make the combined bias, if not using a pre-processed bias from the archive.
# Exposures from additional nights could be included to reduce the noise
# further (assuming they have the same detector region, binning & read mode).
# The default overscan fitting order is 7 for the Hamamatsu detectors (see
# "help gireduce") and should be kept consistent with gsflat & gsreduce.
gemlist "S20150503S" "275-278,280" > "bias_sci.lis"
gbias @bias_sci.lis S20150503S0275_bias rawpath=rawdir$ fl_over+ fl_trim+ \
  fl_vardq+ fl_inter-

# Process the CuAr arc enough to measure a wavelength solution (without
# flat fielding). Since the arc was taken at night with slow read-out in this
# case, the bias frame for science data can be used.
gsreduce S20150503S0212 rawpath=rawdir$ fl_bias+ fl_flat- fl_fixpix- \
    bias=S20150503S0275_bias ovs_flinter-

# (For CuAr taken with the 2.0 arcsec slit, the data must be smoothed so that
# accurate line centers can be found by GSWAVELENGTH. It is also necessary to
# adjust the line width, centering radius and minimum separation parameters.)

# Establish a wavelength solution from the processed arc, using line detection
# parameters that approximately match the slit width & binning (see help). If
# needed, an edited line list may be supplied via the "coordlist" parameter.
# This step should always be run interactively first time, checking the line
# identifications and the resulting wavelength solution. 
gswavelength gsS20150503S0212 fwidth=5 gsigma=0.7 cradius=12 minsep=6 fl_inter-

# Transform (resample) the arc spectrum to linear wavelength, just to check
# that the solution is good. The task gdisplay can be used (with DS9 or ximtool
# already running) to inspect this and other results, in this case checking
# that all the rows have the same alignment (ie. the arc lines are straight).
gstransform gsS20150503S0212 wavtran=gsS20150503S0212

# Reduce and normalize the flat field.
# Scattered light correction is currently available for science data only
# (except when using the IFU), so is not included here.
gsflat S20150503S0211 S20150503S0211_flat order=7 niter=3 ovs_order=7 \
    rawpath=rawdir$ bias=S20150503S0275_bias qe_refim=gsS20150503S0212 \
    fl_qecorr+ fl_vardq+ fl_fulldq+ fl_detec+ ovs_flinter- fl_inter-

# Reduce the observations of the science target, enabling three new options:
# 1.) correction for differences in QE between the detectors as a function of
# wavelength (which is important with the Hamamatsu CCDs, to avoid continuum
# discontinuities; see gqecorr), 2.) subtraction of scattered light and
# 3.) cleaning of cosmic rays using LA Cosmic. The name of the output image
# is determined by the default prefix, outpref="gs".
#
# Steps in gsreduce 
#    subtract off the bias
#    apply QE correction
#    apply flat field correction (before or after mosaicking)
#    clean cosmic rays
#    mosaic the 3-12 detector amplifier arrays
#    interpolate across the chip gaps to aid later reduction steps
#    model & subtract scattered light
#    run gsappwave to get an approximate wavelength calibration in the 
#       header of the output image
#
# As can be seen here, GMOS-S Hamamatsu data taken prior to 26 August 2015
# suffer to varying degrees from spurious background structure on amplifier 5,
# associated with a saturated column. There is no established procedure for
# recovering information faithfully from this amplifier at present. See
# http://www.gemini.edu/sciops/instruments/gmos/status-and-availability. A
# dark band is also observed on amplifier 3 due to another saturated feature.
# Neither of these problems exists in newer datasets.
#
# The cosmic ray correction is slightly less effective than average with the
# 4x4 binning used here (as well as the amp 5 artifact) but still provides a
# considerable improvement with sufficiently conservative parameters. Below,
# the order of the object fit has been lowered and the sigma clipping threshold
# raised, to avoid spurious rejection (see the output DQ extension). Signal
# from a strong sky line still gets rejected erroneously towards the ends of
# the slit, but not coincident with the target spectra.
#
# The apfind column for scattered light correction has been adjusted from its
# default to avoid mistaking cosmic ray residuals for object spectra and the
# sample range has been reduced to avoid a bright artifact from fitting dead
# rows at the bottom of the detector. The orders have also been adjusted to
# avoid overfitting the data. The scattered light subtraction removes the mean
# sky level as well (unlike its counterpart for the IFU, gfscatsub).
#
gemlist "S20150503S" "213-215" > "sci.lis"
gsreduce @sci.lis rawpath=rawdir$ fl_crspec+ fl_qecorr+ fl_scatsub+ fl_vardq+ \
    fl_fulldq+ bias=S20150503S0275_bias flat=S20150503S0211_flat \
    qe_refim=gsS20150503S0212 cr_xorder=3 cr_sigclip=7.0 \
    sc_nfind=2 sc_column=810 sc_order1=4 sc_sample1="40:1030" sc_order2=5 \
    ovs_flinter- fl_inter-

# Interpolate over the amp. 5 saturated columns (586-591), to minimize problems
# with auto-scaling plots etc. A bad pixel mask is created from the first image.
# This should only be included when processing Hamamatsu data taken before/in
# August 2015, but can also be adapted for other bad columns/regions, adding
# instances of "imreplace badcols[dq][x1:x2,y1:y2] 1" after imexpr if needed.
imexpr "I >= 586 && I <= 591 ? 1 : 0*a" badcols[dq] gsS20150503S0213[sci] \
    outtype="ushort"
addbpm gs@sci.lis badcols
gemfix gs@sci.lis pgs@sci.lis method="fixpix" bitmask=1

# Transform (resample) the standard star spectrum to linear wavelength. The
# output image name is determined by the default prefix, outprefix="t".
gstransform pgs@sci.lis wavtran=gsS20150503S0212 fl_vardq+

# Subtract sky using spectra from clean sample regions on either side of the
# target. The sample is selected manually after inspection of the images.
# The high/low rejection limits currently don't get applied properly to the
# second interactive (gfit1d) step, but this only affects the noise estimate.
gsskysub tpgs@sci.lis long_sample="415:465,630:700" fl_vardq+ fl_inter-

# Apply flux calibration derived from a spectrophometric standard star
# observation, reduced separately as in "gmosexamples standard_longslit". The
# calibration file from that reduction, "sens.fits", is assumed to be in the
# current directory. See the help page regarding the default normalization
# (gscalibrate.fluxscale parameter). This step can also be performed after
# extracting a 1D spectrum, instead of here.
gscalibrate stpgs@sci.lis observatory="gemini-south" fl_vardq+

# Inspect the first calibrated spectrum.
imexam cstpgsS20150503S0213.fits[sci] 1

# Extract the primary target near the middle of the slit, using a 3" aperture.
# The starting column has been adjusted (using "coloffset") both to ignore a
# cosmic ray residual and to coincide with an emission line so that apall will
# pick the upper target over the lower one (which has a slightly brighter
# continuum). The intended target(s) can also be specified interactively if
# needed. At present, spectra can only be extracted for multiple targets by
# running gsextract more than once and selecting a different object each time
# (or by defining multiple apertures interactively, causing spectra to be
# stacked in the second image axis). The fl_vardq option appears not to work
# accurately here at present, but the user could attempt to extract the input
# variance in the same manner as the science array, with sci_ext="VAR".
gsextract cstpgs@sci.lis apwidth=3 coloffset=-63 fl_inter-

# Co-add the 3 exposures (this could also be done before extraction, if the
# exposures are taken at the same pointing and the conditions are similar,
# but doing it here is more generally-applicable). The new gemscombine task
# has undergone limited testing and its fl_vardq option is currently omitted
# here, having been observed to cause corruption of the science array. An
# alternative method -- only for exposures taken at the same wavelength
# setting -- is to use gemcombine instead.
gemscombine ecstpgs@sci.lis ecstpgsS20150503S0213_add

# Inspect the extracted spectrum
splot ecstpgsS20150503S0213_add.fits[sci]

