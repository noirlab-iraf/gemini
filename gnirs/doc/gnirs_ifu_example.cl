# Copyright(c) 2005-2012 Association of Universities for Research in Astronomy, Inc.
#
# GNIRS IFU Calibration and Science data reduction script
#
# Observation UT dates:  April 9, 2004
#                        January 29, 2005  (twighlights)
#
# This data was taken during the commissioning of the GNIRS IFU mode.
#
# Required packages:
#   gemini
#   gnirs
#

# Location of the raw data
set rawdir=/export/data/gemini_testdata/gnirs_testdata/ifu/

# Load required tasks/packages
gemini
gemtools
gnirs

#--------------setup for gnirs-----------
# set default parameters
# use gnirs defaults except where specified on command lines below
unlearn gemini
unlearn gemtools
unlearn gnirs

# set display
set stdimage=imt1024

# set logfile and database directory
gnirs.logfile="GS-ENG20040409-CAL.log"
gnirs.database="GS-ENG20040409-CAL_database"

# set raw path
nsprepare.rawpath="rawdir$"

# load gnirs header keywords
nsheaders gnirs

#######################################################
###### Process the GCAL flats for science targets #####
#######################################################
#   Note: INDEF values for shiftx and shifty trigger the automatic
#         cross-correlation algorithm.

gemlist range="154,161,172" root="S20040409S" > flat.lis
nsprepare @flat.lis shiftx=INDEF shifty=INDEF
# shift found: 12 pixels
nsreduce n@flat.lis fl_cut+ fl_nsappw- fl_sky- fl_dark- fl_flat-

# The data should be reduced with the closest flat to account
# for any differences due to flexure.  The flats for the standard
# star and the science data are therefore created separately.
#
# Create flat field for telluric standard
nsflat rnS20040409S0154 darks="" flatfile="" darkfile="" fl_save_dark+ \
   process="fit" thr_flo=0.15 thr_fup=1.55

# Create flat field for the science frames
#   Note: The flats were taken before and after the science frames.  
#	  Since in this case the effect of flexure between first science frame
#         and the last is minimal, the two flats are combined.
nsflat rnS20040409S0161,rnS20040409S0172 darks="" flatfile="" \
   darkfile="" fl_save_dark+ process="fit" thr_flo=0.15 thr_fup=1.55

#####################################################
###### Process the GCAL flats for the twilights #####
#####################################################

gemlist range="5-9" root="S20050129S" > flat2.lis
nsprepare @flat2.lis shiftx=INDEF shifty=INDEF
# shift found: 14.9 pixels
nsreduce n@flat2.lis fl_cut+ fl_nsappw- fl_sky- fl_dark- fl_flat-

# Create flat field for the twilights
nsflat rn@flat2.lis darks="" flatfile="" darkfile="" fl_save_dark+ \
    process="fit" thr_flo=0.15 thr_fup=1.55

#################################
##### Process the twilights #####
#################################

# shiftimage points to a GCAL flat for the twilights to ensure that
# the same MDF shift is used (the 14.9-pixel shift found above).

gemlist range="1-3" root="S20050129S" > twi.lis
nsprepare @twi.lis shiftx=INDEF shiftimage="rnS20050129S0005" \
    bpm="rnS20050129S0005_flat_bpm.pl"
nsreduce n@twi.lis fl_cut+ fl_nsappw- fl_sky- fl_dark- fl_flat-

##########################################################
##### Derive and apply slit functions to flat fields #####
##########################################################

# For the telluric standard flat
nsslitfunction rn@twi.lis "rnS20040409S0154_sflat" \
    flat="rnS20040409S0154_flat" flexflat="rnS20050129S0005_flat" \
    dark="" order=1 verbose+

# For the science flat
nsslitfunction rn@twi.lis "rnS20040409S0161_sflat" \
    flat="rnS20040409S0161_flat" flexflat="rnS20050129S0005_flat" \
    dark="" order=1 verbose+

###########################################################
##### Trace the spatial curvature in the flat fields  #####
###########################################################

# For the science flat
nfflt2pin nS20040409S0161 outimage=nS20040409S0161_pin
nscut nS20040409S0161_pin
nssdist snS20040409S0161_pin fwidth=20 thresh=20 #fl_inter+

# Do the same for the telluric standard flat

################################################################
##### Reduce the arc and determine the wavelength solution #####
################################################################

gemlist range="155" root="S20040409S" > arc.lis
nsprepare @arc.lis shiftimage="nS20040409S0154" \
    bpm="rnS20040409S0154_flat_bpm.pl"
nsreduce n@arc.lis fl_cut+ fl_nsappw+ fl_sky- fl_dark- fl_flat-
nswavelength rn@arc.lis nsum=1 nfound=10 match=-1 thresh=100 fl_inter+

####################################################
##### Basic reduction of the science exposures #####
####################################################

# shiftimage is set to one of the GCAL flat images that
# were used to create the flat field.  The MDF shift to
# apply will
# be taken from that image
gemlist range="162-171" root="S20040409S" > sci.lis
nsprepare @sci.lis shiftimage="calibdir$nS20040409S0161" \
    bpm="calibdir$rnS20040409S0161_flat_bpm.pl"
nsreduce n@sci.lis fl_cut+ fl_nsappwave+ fl_dark- fl_sky+ fl_flat+ \
    skyimages="" skyrange=INDEF flatimage="calibdir$rnS20040409S0161_sflat"

###################################################
##### Combine spectra taken at same pointings #####
###################################################

#  Note: 5 images, 3 nod positions
#        Nod position #1 : S20040409S0162       (single => no need to stack)
#        Nod position #2 : S20040409S0165 to 166
#        Nod position #3 : S20040409S0169 to 170

imdelete ("rnS20040409S0165_stack,rnS20040409S0166_stack",
    verify-, >& "dev$null")
delete ("stack.lis", verify-, >& "dev$null")
gemlist range="165-166,169-170" root="rnS20040409S" > stack.lis
nsstack @stack.lis tolerance=0.1 rejtype="none" lthreshold=-100 hthreshold=150

###################################################
##### Resample the spectra onto a linear grid #####
###################################################

# Need local copy of frames used for calibration
imdelete wrnS20040409S0155,snS20040409S0161_pin ver- >& dev$null
copy calibdir$wrnS20040409S0155.fits .
copy calibdir$snS20040409S0161_pin.fits .

# If calibration database in another directory (e.g. without write permission),
# then need to copy database
#
#if (access("GS-ENG20040409-CAL_database"))
#   delete GS-ENG20040409-CAL_database/* ver- >& dev$null
#else
#   !mkdir GS-ENG20040409-CAL_database
#copy calibdir$GS-ENG20040409-CAL_database/* GS-ENG20040409-CAL_database/

imdelete ("frnS20040409S0162,frnS20040409S0165_stack,frnS20040409S0166_stack",
    verify-, >& "dev$null")
nsfitcoords rnS20040409S0162,rnS20040409S0165_stack,rnS20040409S0166_stack \
    lamp=wrnS20040409S0155.fits sdist=snS20040409S0161_pin \
    database=GS-ENG20040409-CAL_database

imdelete ("tfrnS20040409S0162,tfrnS20040409S0165_stack,  \
    tfrnS20040409S0166_stack", verify-, >& "dev$null")
nstransform frnS20040409S0162,frnS20040409S0165_stack,frnS20040409S0166_stack \
    database=GS-ENG20040409-CAL_database

##################################################
##### Stack the 2D spectra into 3D datacubes #####
##################################################

nfcube tfrnS20040409S0162,tfrnS20040409S0165_stack,tfrnS20040409S0166_stack \
    process="rotate"
