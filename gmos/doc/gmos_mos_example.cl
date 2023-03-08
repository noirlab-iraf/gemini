# Copyright(c) 2002-2015 Association of Universities for Research in Astronomy, Inc.
#
# GMOS example reductions script: Typical reduction of MOS data
#
# This data processing was done to aid the data quality assessment performed
# by the Gemini staff. The data processing is not designed to give the best
# possible result. Better signal-to-noise and better cleaning for cosmic-ray
# hits and bad pixels will most likely be possible. The user of these data
# is encouraged to use the provided co-added images only as guide lines and
# to re-reduce the data to obtain the best possible reduction.
#
# Gemini GMOS data reduction script
# Observation UT date: 2001dec22 
# Data processor:   Inger Jorgensen
# Data reduction date: 2002feb25 (updated 2014jan16 with QE etc. by JT)
#
# Brief data description: MOS observations from the queue program GN-2001B-Q-10
#
# Update: Illustrate the use of the automatic slit edge finding feature (v1.7)
#         Illustrate use of the draft QE correction, LACosmic wrapper, script
#           for combining exposures & option not to oversize slit lengths.

# Define various directories containing raw data and calibration data
set raw=/net/sabrina/staging2/gmos/2001dec22/
set gcalib=../../2001dec25/Basecalib/

# set up the logfile for this reduction
gmos.logfile="GN-2001B-Q-10_log5.txt"

# The bias image has been made with gbias and has been overscan subtracted
# and trimmed.

# Bias subtract the flats. This can be done directly by gsflat but currently
# a break is needed if one wants to apply the new gqecorr and/or gemcrspec
# (because the first gsflat run doesn't keep the intermediate files to which
# we apply the correction).
gsreduce N20011221S139,N20011221S142,N20011221S143 fl_flat- fl_gmosaic- \
    fl_fixpix- fl_gsappwave- fl_cut- rawpath=raw$ \
    bias=gcalib$N20011217S153_bias.fits \
    mdfdir=/usr/dataproc/gmos/GN-2001B-distrib/GN-2001B-Q-10/

# The gsflat call below requires input data to be mosaicked when already bias
# subtracted, so do that here.
gmosaic gsN20011221S139,gsN20011221S142,gsN20011221S143 fl_fixpix+

# Make the initial flat field (as a reference for tracing/cutting)
# To let GSCUT find the edges automatically, fl_usegrad is set to yes.
# Then, since we will later let GSCUT find the slit edges automatically, 
# we need to keep the 'combflat' image.
gsflat mgsN20011221S139,mgsN20011221S142,mgsN20011221S143 \
    rgN20011221S139_flat.fits order=29 rawpath=./ fl_over- fl_trim- fl_bias- \
    fl_usegrad+ combflat=rgN20011221S139_comb.fits fl_keep+ fl_oversize-

# Update the _comb frame created by GSFLAT above with the location of
# the slit edges.
gscut rgN20011221S139_comb.fits fl_update+ fl_oversize- \
    gradimage=rgN20011221S139_comb

# Reduce the CuAr spectrum, using the slit positions determined from the flat.
# The CuAr is not flat fielded. This is needed before we can continue with
# gqecorr for the flats & science data below.
gsreduce N20011221S016 rawpath=raw$ \
    fl_flat- fl_dark- bias=gcalib$N20011217S153_bias.fits \
    mdfdir=/usr/dataproc/gmos/GN-2001B-distrib/GN-2001B-Q-10/ \
    refimage=rgN20011221S139_comb.fits fl_over+

# Establish the wavelength calibration (using the default GMOS line list).
# This step requires a bit of interactions since the automatic line id
# fails in some of the cases.
# For MOS data, the step parameter must be set equal to 2 so that the
# distortion in the spatial direction can be properly determined.
gswavelength gsN20011221S016.fits fwidth=13 cradius=13 minsep=5 aiddebug=s \
    order=4 match=-10 fitcxord=4 fitcyord=4 step=2

# Transform the CuAr spectrum, for checking that the transformation is ok.
# Output image name is defined by the default value outpref="t".
gstransform gsN20011221S016 wavtran=gsN20011221S016

# Inspect the transformed arc by using imexamine on each of the SCI 
# extensions
for(i=1;i<=34;i+=1) {
  imexa("tgsN20011221S016[SCI,"//str(i)//"]")
}

# Apply a correction to the flat exposures for differences in quantum
# efficiency variation with wavelength between CCDs (using the wavelength
# solution determined above) before re-generating the final flat field.
imdel rgN20011221S139_flat.fits
gqecorr gsN20011221S139,gsN20011221S142,gsN20011221S143 \
    refimages=gsN20011221S016 fl_keep+

# Mosaic the QE-corrected data for gsflat as previously.
gmosaic qgsN20011221S139,qgsN20011221S142,qgsN20011221S143 fl_fixpix+

# Re-generate the flat field after correcting for QE varations between
# detectors, re-using the slit positions determined the first time.
gsflat mqgsN20011221S139,mqgsN20011221S142,mqgsN20011221S143 \
    rgN20011221S139_flat.fits order=29 rawpath=./ fl_over- fl_trim- fl_bias- \
    fl_usegrad+ refimage=rgN20011221S139_comb.fits fl_keep-

# In this wavelength range, flats taken before late 2004 contain a pair of
# emission lines due to a fluorescent surface in GCal. Set the pixels in those
# lines to one to avoid introducing spurious features.
imreplace rgN20011221S139_flat.fits[sci,1][1430:3890,*] 1 low=1.02 upper=INDEF

# Reduce the science observations, the name of the output image is 
# defined by the default value outpref="gs"
# Steps in gsreduce 
#    subtract off the bias 
#    clean for cosmic ray hits (if using the old gscrrej method)
#    (apply QE correction, currently done separately)
#    mosaic the 3 detectors 
#    interpolate accross the chip gaps to aid later reduction steps
#    flat field correction 
#    cut the images into slitlets
#    run gsappwave to get an approximate wavelength calibration in the 
#       header of the output image

# Bias-subtract the spectra.
gsreduce N20011221S140,N20011221S141 rawpath=raw$ \
    bias=gcalib$N20011217S153_bias.fits fl_gmosaic- fl_flat- fl_fixpix- \
    fl_gsappwave- fl_cut- fl_over+ fl_gscrrej- \
    mdfdir=/usr/dataproc/gmos/GN-2001B-distrib/GN-2001B-Q-10/

# Use the optional gemcrspec for cosmic ray rejection (instead of enabling
# the older gscrrej option in gsreduce). This requires that you have installed
# lacos_spec.cl separately (see "help gemcrspec").
gemcrspec gsN20011221S140,gsN20011221S141 xgsN20011221S140,xgsN20011221S141

# QE-correct the spectra, using the same corrections derived for the flat.
gqecorr xgsN20011221S140,xgsN20011221S141 refimages=gsN20011221S016 \
    corrimages=qecorrgsN20011221S016.fits

# Finish reducing the spectra with gsreduce.
# 'refimage' contains information on the slit locations that were found
# by GSCUT.
gsreduce qxgsN20011221S140,qxgsN20011221S141 rawpath=./ \
    flat=rgN20011221S139_flat.fits refimage=rgN20011221S139_comb.fits \
    fl_over- fl_trim- fl_bias- fl_gscrrej-

# Transform the science exposures to linear wavelength.
gstransform gsqxgsN20011221S140,gsqxgsN20011221S141 wavtran=gsN20011221S016

# Sky subtract the science exposures. Output image names are defined by the
# default value outpref="s".
gsskysub tgsqxgsN20011221S140,tgsqxgsN20011221S141 mos_sample=0.8 \
    fl_oversize- fl_inter+

# Inspect the sky subtracted spectra by imexamine on each of the SCI 
# extensions
for(i=1;i<=34;i+=1) {
  imexa("stgsqxgsN20011221S140[SCI,"//str(i)//"]")
}
for(i=1;i<=34;i+=1) {
  imexa("stgsqxgsN20011221S141[SCI,"//str(i)//"]")
}

# Test gsextract works, adjust the parameters tsum and tstep to make
# it possible to trace faint objects.
# Run the task interactively to make sure all the faint objects are 
# traced correctly.
gsextract stgsqxgsN20011221S140,stgsqxgsN20011221S141 tnsum=100 tstep=100 \
    fl_inter+

# Inspect the extracted spectra
for(i=1;i<=34;i+=1) {
  splot("estgsqxgsN20011221S140[SCI,"//str(i)//"]")
}
for(i=1;i<=34;i+=1) {
  splot("estgsqxgsN20011221S141[SCI,"//str(i)//"]")
}

# Combine the extracted spectra (aligning them in wavelength if necessary).
# This can also be used on 2D spectra if they are aligned spatially, which
# may help extract faint objects.
gemscombine estgsqxgsN20011221S140,estgsqxgsN20011221S141 \
    estgsqxgsN20011221S140_comb lthreshold=-100. hthreshold=2000.

for(i=1;i<=34;i+=1) {
  splot("estgsqxgsN20011221S140_comb[SCI,"//str(i)//"]")
}

# Flux calibrate 1D spectra, this only works correctly if the sensitivity
# function covers the full wavelength range spanned but the different
# slitlets in the MOS observation. This is not normally the case.
gscalibrate estgsqxgsN20011221S140_comb sfunction=../Basecalib_spec_old/sens

# Inspect the flux calibrated spectra
for(i=1;i<=34;i+=1) {
  splot("cestgsqxgsN20011221S140_comb[SCI,"//str(i)//"]")
}

