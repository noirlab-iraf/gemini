# Copyright(c) 2004-2009 Association of Universities for Research in Astronomy, Inc.
#
# Michelle/T-ReCS example reduction script: 
#			Typical reduction of spectroscopic data
#
# This example is based on Michelle spectroscopic data but the same steps 
# apply to T-ReCS spectroscopic data.
#
# This data processing was done to aid the data quality assessment performed
# by the Gemini staff. The data processing is not designed to give the best
# possible result. The user is encouraged to carry out a more careful reduction
# to improve the resulting signal-to-noise ratio of the final images.
#
# Gemini Michelle data reduction script
# Observation UT date: 2004jul05
#
# Brief data description: Michelle engineering data.
#
#                   GN-ENG-MICHELLE-117 : Spectroscopic observation of BD+30
#                   GN-ENG-MICHELLE-116 : Spectroscopic observation of BS 7796
#
# In this script it is assumed that the raw data files are in directory
#
#	/net/sabrina/staging5/gemini_testdata/michelle_testdata
#
# and that the scripts is run is a different directory (to avoid changing the
# raw data files).  Change the directory name to whatever value is apropriate
# for your setup.
#
# One needs to load the gemini and midir packages before carrying out these
# processing steps.  Note that the midir package automatically load the gnirs
# package.

# Define directory containing raw data
set rawdir=/net/sabrina/staging5/gemini_testdata/michelle_testdata/

# Set up the logfile for this reduction
midir.logfile="michspec.log"

# Set the header keyword definitions
nsheaders michelle

# Create image lists
delete targetlist ver- >& dev$null
delete standardlist ver- >& dev$null
gemlist N20040705S "146" > targetlist
gemlist N20040705S "142" > standardlist

# Create combined difference (spectrum) and reference
# (sky lines) images for the target
imdelete r@targetlist ver- >& dev$null
imdelete c@targetlist ver- >& dev$null
mireduce @targetlist rawpath=rawdir$ frametype="dif"
mireduce @targetlist rawpath=rawdir$ frametype="src"

# Create combined difference (spectrum) and reference
# (sky lines) images for the reference star
imdelete r@standardlist ver- >& dev$null
imdelete c@standardlist ver- >& dev$null
mireduce @standardlist rawpath=rawdir$ frametype="dif"
mireduce @standardlist rawpath=rawdir$ frametype="src"

# Identify sky lines for the target and standard, every tenth
# line in the longslit spectrum of the sky reference images
imdelete wc@targetlist ver- >& dev$null
imdelete wc@standardlist ver- >& dev$null
nswavelength c@targetlist coordlist="gnirs$data/sky.dat" order=2
nswavelength c@standardlist coordlist="gnirs$data/sky.dat" order=2

# Find the wavelength transformation function across the
# longslit sky spectra
imdelete tfr@targetlist ver- >& dev$null
imdelete tfr@standardlist ver- >& dev$null
nsfitcoords r@targetlist lamp=wcN20040705S0146
nstransform fr@targetlist
nsfitcoords r@standardlist lamp=wcN20040705S0142
nstransform fr@standardlist 

# Extract a wavelength calibrated spectrum from an aperture near
# line 101 in the target and standard observations, using the 
# wavelength calibration from the reference images
imdelete xtfr@targetlist ver- >& dev$null
imdelete xtfr@standardlist ver- >& dev$null
nsextract tfr@targetlist line=101 fl_inter-
nsextract tfr@standardlist line=101 fl_inter-

# Finally, ratio the target and standard spectrum and then 
# multiply by a normalized blackbody spectrum to produce a 
# "calibrated" target spectrum.  This spectrum should have 
# the proper shape but does not have any overall normalization 
# to the magnitude of the standard in N-band.
#
# "mstelluric" calls "telluric" from the noao.onespec package.
#
# The blackbody temperature used makes little difference at 
# these wavelengths as long as it is sufficiently high (i.e. 
# more than 3000).  A default value of 10000 is recommended.
#
# For best results it may be necessary to re-run "mstelluric" 
# interactively and find the best shift/scaling from the 
# standard spectrum to the target spectrum.  The following 
# call uses the automatic processing:
#
imdelete axtfr@targetlist ver- >& dev$null
mstelluric xtfrN20040705S0146 xtrN20040705S0142 bbody=10000 xcorr- fl_inter-

# Inspect the output spectrum
splot axtfrN20040705S0146.fits[1]

# Alternative with interactive cross-correlation:
#
#   First enter the interactive cross-correlation routine when prompted.
#   Use ":shift 0.0" and ":scale 0.95" to get the best result 
#   (smallest residual ozone band and residual atmospheric lines).  
#   Then quit with "q".
#
mstelluric xtfrN20040705S0146 xtrN20040705S0142 bbody=10000 \
    xcorr+ fl_inter+ outpref="z"

# Inspect the second output spectrum
splot zxtfrN20040705S0146.fits[1]

# In this case the second spectrum is better.
