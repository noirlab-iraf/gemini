# Copyright(c) 2002-2006 Association of Universities for Research in Astronomy, Inc.
#
# NIRI example reductions script: Typical reduction of imaging data 
#   updated for 16 Sep 2004 NIRI package release
#   updated for May 2005 new scripts
#
# This data processing was done to aid the data quality assessment performed
# by the Gemini staff. The data processing is not designed to give the best
# possible result. Better signal-to-noise and better cleaning for cosmic-ray
# hits and bad pixels will most likely be possible. The user of these data
# is encouraged to use the provided co-added images only as guide lines and
# to re-reduce the data to obtain the best possible reduction.
#
# Gemini NIRI queue data reduction script
# Observation UT date:  2002 Jun 28
# Data filename prefix: N20020628S
#
# IMAGING CALIBRATIONS:
# darks				121-130
# GCAL J flat shutter closed	131-134,136-137,139-142
# GCAL J flat shutter open	144-153                 nN20020628S0145_flat
# FS135 J			160-164
# FS135 H2v=1-0			170-174
# darks				246-255




############################
#
# !!!! WARNING !!!!
#
#   It is highly recommended to run the following 'unlearn' commands
#   BEFORE you start.
#
#   ni> unlearn gmos
#   ni> unlearn gemtools
#   ni> unlearn niri
#
############################


string name,rawdir
name = "N20020628S"
rawdir = "/net/sabrina/staging2/niri/2002jun28"
niri.logfile="20020628_GN-CAL20020628.log"

####################
# J-band flat      #
# GN-CAL20020628-4 #
####################

# Imaging flat field and bad pixel mask derived from the calibration unit.
# Start by generating the file lists and updating the headers:

delete darklist ver-
for(i=123; i<=130; i+=1) {
  print(name//"0000"+i, >> "darklist")
}
delete offlist ver-
for(i=132; i<=142; i+=1) {
  if(i!=135 && i!=138) print(name//"0000"+i, >> "offlist")
  ;
}
delete onlist ver-
for(i=145; i<=153; i+=1) {
  print(name//"0000"+i, >> "onlist")
}
delete inlist ver-
concatenate darklist,onlist,offlist > inlist
nprepare("@inlist", rawpath=rawdir)

# Construct the normalized flat field, using short darks to identify
# bad pixels.  The imaging flats are derived from images taken with
# the calibration unit shutter closed ("lamps off") subtracted from
# equal exposres taken with the shutter open ("lamps on").

delete nonlist ver-
sections n@onlist > nonlist
delete nofflist ver-
sections n@offlist > nofflist
delete ndarklist ver-
sections n@darklist > ndarklist

niflat @nonlist lampsoff=@nofflist darks=@ndarklist


####################
# J-band standard  #
# FS 135           #
# GN-CAL20020628-5 #
####################

# Sky-subtract, flatten, and combine 5 images of a standard star.
# Start by updating the headers, including the BPM in the data quality
# plane of the raw images:

delete inlist ver-
for(i=160; i<=164; i+=1) {
	print(name//"0000"+i, >> "inlist")
}
nprepare("@inlist", rawpath=rawdir, bpm="nN20020628S0145_bpm.pl")

# Flag any residual images from previously saturated or non-linear
# pixels

delete ninlist ver-
sections n@inlist > ninlist
nresidual @ninlist proptime=0.5

# Construct a sky image by identifying objects in each frame, removing
# them, and averaging the remaining good pixels. 

delete bninlist ver-
sections b@ninlist > bninlist
nisky @bninlist

# Reduce the raw images by subtracting the sky and dividing by the 
# normalized flat field image:

nireduce @bninlist skyimage=bnN20020628S0160_sky.fits \
flatimage=nN20020628S0145_flat.fits

# Now combine the individual images
delete rbninlist ver-
sections r@bninlist > rbninlist

imcoadd @rbninlist logfile="20020628_GN-CAL20020628.log" geofitgeom=shift \
rotate-  fl_over+ fl_scale- \
fl_fixpix+ fl_find+ fl_map+ fl_trn+ fl_med+ fl_add+ fl_avg+ \
badpix="nN20020628S0145_bpm.pl" niter=1

# To improve the sky masks, use nisupersky and repeat reduction:
nisupersky rnN20020628S0160_add

rm rnN20020628S016*.fits
rm rnN20020628S0160_sky.fits

nisky @bninlist

nireduce @bninlist skyimage=bnN20020628S0160_sky.fits \
flatimage=nN20020628S0145_flat.fits

imcoadd @rbninlist logfile="20020628_GN-CAL20020628.log" geofitgeom=shift \
rotate-  fl_over+ fl_scale- \
fl_fixpix+ fl_find+ fl_map+ fl_trn+ fl_med+ fl_add+ fl_avg+ \
badpix="nN20020628S0145_bpm.pl" niter=1



#####################
# H_2 band standard #
# FS 135            # 
# GN-CAL20020628-7  #
#####################

# Here's an example of how a sky flat can be used instead of the
# calibration unit flat.

# Start by updating the headers:

delete inlist ver-
for(i=170; i<=174; i+=1)
	print(name//"0000"+i, >> "inlist")
nprepare("@inlist", rawpath=rawdir)

delete darklist ver-
for(i=248; i<=255; i+=1) 
	print(name//"0000"+i, >> "darklist")
nprepare("@darklist", rawpath=rawdir)

# Make the sky flat using NIFLAT with the flags set appropriately:

delete darklist ver-
for(i=248; i<=255; i+=1) 
	print("n"//name//"0000"+i, >> "darklist")
delete ninlist ver-
sections n@inlist > ninlist
niflat @ninlist lampsoff=@darklist darks=@darklist fl_rmstar+ fl_keepmasks+

# Now make the sky image to subtract:
nisky @ninlist

# Subtract the sky and divide by the sky flat:
nireduce @ninlist skyimage=nN20020628S0170_sky.fits \
flatimage=nN20020628S0170_flat.fits

# Offset and average the individual images to get the final image:
delete rninlist ver-
sections r@ninlist > rninlist

imcoadd @rninlist logfile="20020628_GN-CAL20020628.log" geofitgeom=shift \
rotate-  fl_over+ fl_scale- fwhm=9.5 \
fl_fixpix+ fl_find+ fl_map+ fl_trn+ fl_med+ fl_add+ fl_avg+ \
badpix="nN20020628S0170_bpm.pl" niter=1


####################################
# Example of WCS fix for subarrays #
# Example of WCS de-rotation       #
####################################
 
# nprepare subarray images to fix the WCS:
delete testlist ver-
print("N20040805S0102", >> "testlist")
print("N20040805S0112", >> "testlist")
print("N20040805S0122", >> "testlist")
nprepare @testlist
 
# de-rotate data taken with the CR fixed
delete ntestlist ver-
sections n@testlist > ntestlist
nirotate ntestlist
 
delete ntestlist,testlist ver-
