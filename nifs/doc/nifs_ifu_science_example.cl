# Copyright(c) 2006-2012 Association of Universities for Research in Astronomy, Inc.
#
# Gemini NIFS data reduction script
# Reduction for:  SCIENCE DATA
# 
# Processed data taken on : 2006Feb10
# Data processor          : tbeck
# Data reduction date     : 2006mar23
#
# File rootname: N20060210S*
# Files for Science Reduction
#   178-190 Science frames
#   195     Shift File (from Basecalib Reduction)
#   195     Flat  (from Basecalib Reduction)
#   389     Ronchi Flats  (from Basecalib Reduction)
#   195     BPM  (from Basecalib Reduction)
#   191     Arc  (from Basecalib Reduction)

###########################################################################
# DATA REDUCTION HINT-   		                                  #
#                                                                         #
# For a general science reduction:                                        #
# Users need to change a few lines in "STEP 2" for their data.  Change    #
# the file names for the flat, arc, Ronchi calibration image, telluric    #
# standard and the first on-source science exposure reference frame.      #
# Also update the directory links for the raw data and calibration data,  #
# and change the listing information for making the science list.         #
#                                                                         #
# Run the routine by cutting and pasting the reduction steps into your    # 
# IRAF session. Highlighting and copying the whole reduction at once      #
# seems to work fine.                                                     #
# (Defining and running the script as a reduction "task" doesn't work     #
# well because the loading of the gemini package in "STEP 1" overwrites   #
# the defined task!)                                                      #
#                                                                         #
# The sky subtraction list of files will need edited if you do not have   #
# the same number of sky files as science frames.  (See below notes).     #
#                                                                         #
# Note that this reduction does the extraction of the science spectra and #
# telluric correction in the "nftelluric" step interactively, and requires# 
# user input at this stage.                                               #
#                                                                         #
###########################################################################

###############################################################
# STEP 1:  PREPARE IRAF                                       #
###############################################################

gemini
gemtools
gnirs
nifs

unlearn gemini gemtools gnirs nifs

set stdimage=imt2048
nsheaders("nifs")

###############################################################
# STEP 2: SET REDUCTION FILE NAMES AND PATHS                  #
###############################################################

string dark, calflat,  adata, skyflat, arc, arcdark, ronchiflat,
       sciencedark, scienceframes, flatdark, raw_data, reference,
       skyfile, scifile, cal_data, telluric
int    junk

user_clobber  = envget("clobber")
reset clobber = yes
log_file      = "nifs.log"

# Change the "raw_data" and "cal_data" to correspond to the
# directory paths for your raw science frames and processed 
# calibrations,respectively.  If all of your data is in your 
# working directory, change both values to "" (null strings).

     raw_data="/net/archie/staging/perm/"
     cal_data=""

     calflat    = "N20060210S0195"    # Only the first frame if many
     arc        = "N20060210S0191"    # Only the first frame if many
     ronchiflat = "N20060210S0389"    # Only the first frame if many
     reference  = "N20060210S0178"    # First on-source science file for
                                      # sky subtraction reference
     telluric   = "N20060210S0142"     # filename for the telluric
                                      # calibrator star
#If your telluric data were reduced in a different way (i.e., because of 
#different offset patterns, for example), change the below prefix to
#correspond to your data.

telluric="gxtfbrsn"//telluric

del sciencelist ver-
gemlist N20060210S 178-185,187-190 > sciencelist

###########################################################################
# STEP 3:  Get the Calibrations for the Reduction                         #
# NOTE: Do this step only if your calibration directory is different      #
# from your current working directory.                                    #
###########################################################################

copy(cal_data//"rgn"//ronchiflat//".fits",output="./")
copy(cal_data//"wrgn"//arc//".fits",output="./")
mkdir("database")
copy(cal_data//"database/*",output="./database/")

###########################################################################
# STEP 4:  Reduce the Science Data                                        #
###########################################################################

nfprepare("@sciencelist", rawpath=raw_data,
   shiftimage=cal_data//"s"//calflat,fl_vardq+,
   bpm=cal_data//"rn"//calflat//"_sflat_bpm.pl",
   logfile="nifs.log")

############################################################
#  DATA REDUCTION HINT -                                   #
#                                                          #
# In the below call to "gemoffsetlist", to construct files #
# with the sky offset information, change the "distance"   #
# parameter to correspond to the offset distance (in       #
# arcseconds) for your observations.                       #
#                                                          #
############################################################

gemoffsetlist("n@sciencelist",targetli="scilist",
   offsetli="skylist", reffile="n"//reference//".fits", 
   distance=10.0, age=INDEF,logfile="nifs.log")
 
#############################################################
#  DATA REDUCTION HINT -                                    #
#                                                           #
# At the present time, we found that there are problems     #
# with the WCS coordinates.  The automatic sky              #
# frame ID and subtraction does not work very well in       #
# "nsreduce" for NIFS data reductions if the number of sky  #
# frames does not equal the number if science frames.  As   #
# a result, the sky subtraction in this script was set up   #
# to work outside of the "nsreduce" call.  This should work #
# for most modes of science acquisition.  However you do    #
# have to edit the "skylist" generated by the above         #
# "gemoffsetlist" call to make sure that the number of      #
# sky images and science images is the same in both lists.  #
#                                                           #
#############################################################

#If you have more than 1 science position per sky field, make
#multiple input science lists and from here down re-run each 
#step for each science position.

# Do sky subtraction:
# Make sure to edit the "skylist" to have the same number
# of images in the list as the "scilist"
string *scilis
string *skylis
 scilis = "scilist"
 skylis = "skylist"
while (fscan (scilis, scifile) != EOF) {
    junk = fscan (skylis, skyfile)
    gemarith (scifile, "-", skyfile, "g"//scifile, fl_vardq+, 
    logfile="nifs.log")
}

#flat field and cut the data
nsreduce("g@scilist", fl_cut+, fl_nsappw+, fl_dark-, fl_sky-, 
   fl_flat+, flatimage=cal_data//"rn"//calflat//"_flat",
   fl_vardq+,logfile="nifs.log")

#interpolate over bad pixels flagged in the DQ plane
nffixbad("rg@scilist",logfile="nifs.log")

#derive the 2D to 3D spatial/spectral transformation
nsfitcoords("brg@scilist",lamptransf="wrgn"//arc, 
   sdisttransf="rgn"//ronchiflat,logfile="nifs.log")

#apply the transformation determined in the nsfitcoords step
nstransform("fbrg@scilist", logfile="nifs.log")

#correct the data for telluric absorption features
#NOTE:  the default for nftelluric is to run the task interactivly
# to determine the shift/scale values.  See the help file
# for nftelluric for additional information.
nftelluric("tfbrg@scilist", telluric, logfile="nifs.log")

#reformat the data into a 3-D datacube
nifcube ("atfbrg@scilist", logfile="nifs.log")

###########################################################################
# Reset to user defaults                                                  #
###########################################################################
if (user_clobber == "no")
   set clobber = no
;

###########################################################################
#          End of the Science Data Reduction                              #
#                                                                         #
# The output of this reduction is a set of 3-D data cubes that have been  #
# sky subtracted, flat fielded, cleaned for bad pixels, telluric          #
# corrected and rectified into a cohesive datacube format.  In the case   #
# of this reduction, the final output files are called: catfbrgn+science, #
# or:   catfbrgnN20060210S0178                                            #
#       catfbrgnN20060210S0180                                            #
#       catfbrgnN20060210S0181                                            #
#       catfbrgnN20060210S0183                                            #
#       catfbrgnN20060210S0184                                            #
#       catfbrgnN20060210S0187                                            #
#       catfbrgnN20060210S0188                                            #
#       catfbrgnN20060210S0190                                            #
#                                                                         #
# The meaning of the output prefixes are described below:                 #
#                                                                         #
# g = gemcombined   n=nfprepared  s=skysubtracted   r=nsreduced           #
# b = bad pixel corrected  f= run through nsfitcoords                     # 
# t = nstransformed   a = corrected for telluric absorption features      #
# c = rectified to a 3D datacube                                          #
#                                                                         #
# This script is meant to be a guideline as a method of a typical data    #
# reduction for NIFS frames.  Of course, NIFS PIs can add or skip steps   #
# in this reduction as they deem fit in order to reduce their particular  #
# datasets.                                                               #
#                                                                         #
###########################################################################

