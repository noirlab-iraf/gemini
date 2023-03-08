# Copyright(c) 2002-2012 Association of Universities for Research in Astronomy, Inc.
#
# NIRI example reductions script: Typical reduction of spectroscopy data 
#
# WARNING: Redirection should not be used to run this script (e.g.
#       cl < GN-CAL20020628_spec_example.cl).  You should rather define the
#       script as task
#
#            ni> task $example = niri$doc/GN-CAL20020628_spec_example.cl
#
#       The IRAF task IDENTIFY will not work correctly if the redirection
#       is used (IDENTIFY used the standard input).
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
# Data processor: T. Beck, J. Jensen
# Data reduction date:  2002 Jun 28, 2002 Aug 19
# Data filename prefix: N20020628S
#
# Updated to use GNIRS package, 2004 Oct 12, KL




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
#   ni> unlearn gnirs
#
############################

string name,rawdir
name = "N20020628S"
rawdir = "/net/sabrina/staging2/niri/2002jun28"

#
# Brief data description:     files	ObsID
#   Hip69931 K-grism          49-52     2 
#   GCAL flat K-grism         67-76	3
#   Ar K-grism                46,360	1,15
#   dark 10s		      56-63,65-66	3                     

# set up the logfile for this reduction
niri.logfile="20020628_GN-CAL20020628.log"

######################
# Configure for NIRI #
######################

nsheaders niri

####################
# K-band flat      #
# GN-CAL20020628-3 #
####################

# Spectroscopic flat field and bad pixel mask from the calibration unit.
# Start by updating the headers and identifying and extracting the SPECSEC 
# we will use:

# Create a couple lists (flatlist and darklist)
delete darklist ver-
for(i=56; i<=63; i+=1) {
  print(name//"0000"+i, >> "darklist")
}
print(name//"0000"+65, >> "darklist")

delete flatlist ver-
for(i=67; i<=76; i+=1) {
  print(name//"0000"+i, >> "flatlist")
}

# Prepare the darks and flats
delete inlist ver-
concatenate darklist,flatlist > inlist
nprepare("@inlist", rawpath=rawdir, fl_vardq-)


# Cut the flats.  Do not cut the darks, let NSFLAT handle the darks.
delete nflatlist ver-
sections n@flatlist > nflatlist
nscut @nflatlist

# Now construct the normalized flat field and bad pixel mask (short
# darks are used for identifying bad pixels):

delete snflatlist,ndarklist ver-
sections sn@flatlist > snflatlist
sections n@darklist > ndarklist
nsflat @snflatlist darks=@ndarklist 

####################
# Arc lamps  K     #
# GN-CAL20020628-1 #
####################

# Argon arc lamp spectra for determining the wavelength calibration.
# Start by updating the headers and identifying and extracting the
# SPECSEC we will use.  Add approximate wavelength calibration to use
# as a starting point for NSWAVELENGTH:

delete arclist ver-
for(i=46; i<=46; i+=1) {
	print(name//"0000"+i, >> "arclist")
}
nprepare("@arclist",rawpath=rawdir,fl_vardq+,bpm="snN20020628S0067_bpm.pl")
delete narclist ver-
sections n@arclist > narclist
nscut @narclist

# Identify the lines and construct the final wavelength calibration:
# This will run non-interactively, but it is better to check it and
# make sure.  Delete spurious IDs with "d" and refit with "f".

delete snarclist ver-
sections s@narclist > snarclist
nswavelength @snarclist fl_inter+

# Do the other arc lamp spectrum, exactly the same way:

delete arclist ver-
for(i=360; i<=360; i+=1) {
	print(name//"0000"+i, >> "arclist")
}
nprepare("@arclist",rawpath=rawdir,fl_vardq+,bpm="snN20020628S0067_bpm.pl")
delete narclist ver-
sections n@arclist > narclist
nscut @narclist

delete snarclist ver-
sections s1@narclist > snarclist
nswavelength @snarclist fl_inter+


############################
# S-distortion measurement #
# 2002 Jun 30              #
# GN-ENG20020630-6         #
############################

# The S-distortion was measured by stepping a star along the slit
# at 3 arcsec intervals.  Start by updating the headers, extracting
# SPECSEC, and putting the approximate wavelength solution in the
# headers:

delete("inlist",ver-)
delete("ninlist",ver-)
delete("rinlist",ver-)
delete("sinlist",ver-)
for(i=259; i<=275; i+=1) {
  print("N20020630S0000"+i, >> "inlist")
}
sections n@inlist > ninlist
sections s@ninlist >sinlist
sections r@sinlist > rinlist

nprepare @inlist rawpath="/net/sabrina/staging2/niri/2002jun30/" \
fl_vardq-

nscut @ninlist

# Next subtract the sky for each stepped positon:

nsreduce @sinlist fl_flat- fl_vardq- fl_cut- fl_process_cut+ skyrange=1800

# Now determine the S-distortion.  NSSDIST starts by combining all the
# separate star spectra into one image.  Use NSSDIST interactively for
# now (the non-interactive mode is not yet reliable; debugging is in
# progress); use "m" to mark stars, "f" to fit, and then "q"

nssdist @rinlist fl_inter+ firstycoord=6. fwidth=16.


####################
# Hip69931 K       #
# GN-CAL20020628-2 #
####################

# This standard star was observed using an ABBA offset pattern.
# Start by fixing the headers, selecting SPECSEC1, cutting it out,
# and adding an approximate wavelength calibration:

delete stdlist ver-
for(i=49; i<=52; i+=1) {
	print(name//"0000"+i, >> "stdlist")
}
nprepare("@stdlist",rawpath=rawdir,fl_vardq+,bpm="snN20020628S0067_bpm.pl")
delete nstdlist ver- 
sections n@stdlist > nstdlist
nscut @nstdlist

# Next combine the "A" positions and subtract from the "B" positions
# and vice versa.  Divide by the normalized flat:

delete sstdlist ver- 
sections s@nstdlist > sstdlist
nsreduce @sstdlist fl_vardq+ fl_cut- fl_proc+ skyrange=120 \
    flatim=snN20020628S0067_flat.fits

# Next we combine the two "A" positions and the two "B" positions
# without offsets

delete slist ver-
sections rsn@stdlist > slist
nsstack @slist

# Now offset the resulting _stack spectra and combine them.  The 
# output _add image has ALL the flux in the positive spectrum, so there
# is no need to extract the negative spectra.

delete slist ver-
print("rsnN20020628S0049_stack", >> "slist")
print("rsnN20020628S0050_stack", >> "slist")
nscombine @slist rejtype=sigclip lsigma=3 hsigma=3

# Let's rectify and wavelength calibrate the result:
# (it runs non-interactively, but it is always safer to have a look...)

delete "tmpout" ver-
nsfitcoords rsnN20020628S0049_stack_comb outspectra="tmpout" \
outprefix = ""  lamp=wsnN20020628S0046 sdist=rsnN20020630S0259_sdist \
fl_inter+ lxorder=4 lyorder=3 sxorder=4 xyorder=3

nstransform inimages="tmpout" outspectra="trsnN20020628S0049_stack_comb"

# residual sky line removal is not needed here, but works

nsressky trsnN20020628S0049_stack_comb

# Extract a 1-d spectrum (once again, the negative spectra are not
# used; NSSTACK+NSCOMBINE already included them in the positive spectrum)

nsextract trsnN20020628S0049_stack_comb fl_apall+ fl_findneg- fl_inter- \
    fl_trace+ lower=-5 upper=5

# To look at the results,

splot xtrsnN20020628S0049_stack_comb[sci,1]

# Telluric feature removal by dividing by a standard is done using
# NSTELLURIC


