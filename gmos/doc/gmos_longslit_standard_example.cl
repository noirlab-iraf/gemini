# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
#
# GMOS example reduction script: Typical reduction of spectrophomotometric 
#                                standard star observation
#
# This example should only be used as a guide and can likely be improved on,
# eg. for better cosmetics, wavelength accuracy and signal-to-noise ratio; the
# user is encouraged to optimize the sequence carefully for his/her particular
# dataset (see each task's help pages for available options).
#
# Data taken from GS-2015A-SV-201.
#
# Define directory containing raw data. A separate directory containing
# pre-processed calibration data could also be defined here, if applicable.
set rawdir="raw/"

# Set up the logfile for this reduction
gmos.logfile="GS-2015A-SV-201.log"

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
# delete *.log,*.lis,tmp*,std
# imdel *.fits
# !rm -fr database/

# Make the combined bias, if not using a pre-processed bias from the archive.
# Exposures from additional nights could be included to reduce the noise
# further (assuming they have the same detector region, binning & read mode).
gemlist "S20150511S" "341-345" > "bias_std.lis"
gbias @bias_std.lis S20150503S0341_bias rawpath=rawdir$ fl_over+ fl_trim+ \
    fl_vardq+ fl_inter-

# Process the CuAr arc enough to measure a wavelength solution (without
# flat fielding). Leave out the full bias correction, which would require
# a separate bias frame for fast read-out.
gsreduce S20150511S0265 rawpath=rawdir$ fl_bias- fl_flat- fl_fixpix- \
    ovs_flinter-

# Establish a wavelength solution from the processed arc, using line detection
# parameters that approximately match the slit width & binning (see help).
# Since this baseline calibration exposure was taken separately in the morning,
# its wavelength zero point will likely differ from the target's by a few
# pixels of flexure, but this is not normally significant for flux calibration.
# If needed, an edited line list may be supplied via the "coordlist" parameter.
# This step should always be run interactively first time, checking the line
# identifications and the resulting wavelength solution.
gswavelength gsS20150511S0265 fwidth=7 gsigma=1.5 cradius=12 minsep=7 fl_inter-

# Transform (resample) the arc spectrum to linear wavelength, just to check
# that the solution is good. The task gdisplay can be used (with DS9 or ximtool
# already running) to inspect this and other results, in this case checking
# that all the rows have the same alignment (ie. the arc lines are straight).
gstransform gsS20150511S0265 wavtran=gsS20150511S0265

# Reduce and normalize the flat field
# Scattered light correction is currently available for science data only
# (except when using the IFU), so is not included here.
gsflat S20150511S0216 S20150511S0216_flat order=9 niter=3 rawpath=rawdir$ \
    bias=S20150503S0341_bias qe_refim=gsS20150511S0265 fl_qecorr+ fl_vardq+ \
    fl_fulldq+ fl_detec+ ovs_flinter- fl_inter-

# Reduce the observations of the star, enabling correction for differences
# in QE between the detectors as a function of wavelength (which is important
# with the Hamamatsu CCDs, to avoid continuum discontinuities; see gqecorr)
# and subtraction of scattered light. The name of the output image is
# determined by the default prefix, outpref="gs".
#
# Steps in gsreduce 
#    subtract off the bias
#    apply QE correction
#    apply flat field correction (before or after mosaicking)
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
# http://www.gemini.edu/sciops/instruments/gmos/status-and-availability.
#
# Here gsscatsub's nfind parameter is overridden because there are additional
# sources visible in the slit. Only 5 of 8 visible sources get auto-detected
# and excluded but the others are comparatively very faint and do not bias the
# low-order background fit significantly, as can be seen in the plots. Users
# requiring absolute spectrophotometric calibration should define the regions
# carefully enough to avoid subtracting signal in the wings.
#
gsreduce S20150511S0215 rawpath=rawdir$ fl_qecorr+ fl_scatsub+ fl_vardq+ \
    fl_fulldq+ bias=S20150503S0341_bias flat=S20150511S0216_flat \
    qe_refim=gsS20150511S0265 sc_nfind=6 sc_torder=5 sc_order1=4 sc_order2=5 \
    ovs_flinter- fl_inter-

# Transform (resample) the standard star spectrum to linear wavelength:
gstransform gsS20150511S0215 wavtran=gsS20150511S0265 fl_vardq+

# Sky subtract the standard star spectrum, using spectra from clean sample
# regions on either side of the target. The sample range is selected based on
# inspection of the images. The output image name is determined by the default
# outpref="s".
gsskysub tgsS20150511S0215 long_sample="110:160,320:370" fl_vardq+ fl_inter-

# Extract the 1D spectrum (central 1 arcsec by default).
gsextract stgsS20150511S0215 fl_inter-

# Extablish the sensitivity function. The standard star is EG131, for which a
# flux table is kept in gmos$calib/ (the others live under onedstds$). The
# flux band overlapping the bad column could be deleted interactively or
# excluded from a copy of the table eg131.dat if needed, but has minimal
# effect on the low-order fit in any case.
gsstandard estgsS20150511S0215 std sens starname=EG131 caldir=gmos$calib/ \
    order=5 fl_inter-

# Apply the flux calibration to the extracted 1D spectrum
# The output from gsstandard is found using the default filename for the
# sensitivity function, sfunction = "sens".
gscalibrate estgsS20150511S0215 observatory="gemini-south"

