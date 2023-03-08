# Copyright(c) 2006-2012 Association of Universities for Research in Astronomy, Inc.
#
# Gemini NIFS data reduction script
# Reduction for:  TELLURIC STANDARD CALIBRATIONS
#
# Processed data taken on : 2006Feb10
# Data processor          : tbeck
# Data reduction date     : 2006Mar23
#
# File rootname: N20060210S 
# Files 142-146 : H-band data taken with a Box5 offset pattern.
#
###########################################################################
# DATA REDUCTION HINT-   		                                  #
#                                                                         #
# For a general telluric calibration reduction:                           #
# Users need to change a few lines in "STEP 2" for their data.  Change    #
# the flat, arc and Ronchi calibration file names.  Change the directory  #
# links for the raw data and calibration data, and update the listing     #
# information for the telluric list.                                      #
#                                                                         #
# Run the routine by cutting and pasting the reduction steps into your    # 
# IRAF session. Highlighting and copying the whole reduction at once      #
# seems to work fine.                                                     #
# (Defining and running the script as a reduction "task" doesn't work     #
# well because the loading of the gemini package in "STEP 1" overwrites   #
# the defined task!)                                                      #
#                                                                         #
# In this script,               [2]       [5]                             #
# standard star                                        q                  #
# measurements are                   [1]               ^                  #
# assumed to be taken                                  |                  #
# with a box5 pattern.          [3]       [4]     p <--+                  #
#                                                                         #
# A sky frame is constructed by median combining frames [2]-[5].  Editing #
# the way the sky subtraction is done should be easy if telluric data     #
# were obtained by offsetting to the sky (In this case, just look at the  #
# way the NIFS Science reduction is constructed).                         #
#                                                                         #
# Note that this reduction does the extraction of the telluric spectra    #
# in the "nfextract" step interactively, and requires user input at this  #
# stage.                                                                  #
#                                                                         #
###########################################################################

###########################################################################
# STEP 1:  Prepare IRAF  		                                  #
###########################################################################

gemini
gemtools
gnirs
nifs

unlearn gemini gemtools gnirs nifs

set stdimage=imt2048
nsheaders("nifs",logfile="nifs.log")

###########################################################################
# STEP 2:  Define Variables and Reduction Lists                           #
###########################################################################

string calflat, arc, ronchiflat, cal_data, telluric,
       raw_data, user_clobber, band, log_file

user_clobber  = envget("clobber")
reset clobber = yes
log_file      = "nifs.log"

# Change the below 3 lines for your own data!

# Change the "raw_data" and "cal_data" to correspond to the
# directory paths for your raw science frames and processed 
# calibrations,respectively.  If all of your data is in your 
# working directory, change both values to "" (null strings).

raw_data = "/net/archie/staging/perm/"
cal_data   = ""

     calflat    = "N20060210S0195"    # Only the first frame if many
     arc        = "N20060210S0191"    # Only the first frame if many
     ronchiflat = "N20060210S0389"    # Only the first frame if many

del telluriclist ver-
del skylist ver-
gemlist N20060210S 142-146 > telluriclist
system.tail "telluriclist" nlines=4  > skylist
type("telluriclist")       | scan(telluric)

###########################################################################
# STEP 3:  Get the Calibrations for the Reduction                         #
###########################################################################

copy(cal_data//"rgn"//ronchiflat//".fits",output="./")
copy(cal_data//"wrgn"//arc//".fits",output="./")
mkdir("database")
copy(cal_data//"database/*",output="./database/")

###########################################################################
# STEP 4:  Reduce the Telluric Standard                                   #
###########################################################################

#prepare the data
nfprepare("@telluriclist",rawpath=raw_data,shiftim=cal_data//"s"//calflat,
  bpm=cal_data//"rn"//calflat//"_sflat_bpm.pl",fl_vardq+,
  fl_int+,fl_corr-,fl_nonl-)

#make a median combined sky from the offset frames
gemcombine("n//@skylist",output="n"//telluric//"_sky",
  fl_dqpr+,fl_vardq+,masktype="none",combine="median",logfile="nifs.log")

#do the sky subtraction on all the frames
string *telllis
string tellfile
 telllis = "telluriclist"
while (fscan (telllis, tellfile) != EOF) {
gemarith("n"//tellfile,"-","n"//telluric//"_sky","sn"//tellfile,
  fl_vardq+,logfile=log_file)
}

#reduce and flat field the data
nsreduce("sn@telluriclist",outpref="r", 
  flatim=cal_data//"rn"//calflat//"_flat",
  fl_cut+,fl_nsappw-,fl_vardq+,fl_sky-,fl_dark-,fl_flat+,logfile=log_file)

#fix bad pixels from the DQ plane
nffixbad("rsn@telluriclist",outpref="b",logfile=log_file)

#derive the 2D to 3D spatial/spectral transformation
nsfitcoords("brsn@telluriclist",outpref="f",
    lamptr="wrgn"//arc, sdisttr="rgn"//ronchiflat,logfile=log_file)

#apply the transformation determined in the nsfitcoords step
nstransform("fbrsn@telluriclist",outpref="t",logfile=log_file)

#extract 1D spectra from the 2D data
#NOTE: In order to run nfextract interactively you need an image 
#window open. To run interactively, hit any to mark the aperture
#position and display the extracted spectra.  Hit "q" to continue
#to the next spectrum.
nfextract("tfbrsn@telluriclist",outpref="x",diameter=0.5, fl_int+,
   fl_zval+,z1=0,z2=10000.0, logfile=log_file)

#combine all the 1D spectra to one final output file
gemcombine("xtfbrsn//@telluriclist",output="gxtfbrsn"//telluric,
   statsec="[*]", combine="median",logfile=log_file,masktype="none",
   fl_vardq+)


###########################################################################
# Reset to user defaults                                                  #
###########################################################################
if (user_clobber == "no")
   set clobber = no
;
###########################################################################
#          End of the Telluric Calibration Data Reduction                 #
#                                                                         #
#  The output of this reduction script is a 1-D spectrum used for         #
# telluric calibration of NIFS science data.  For this particular         #
# reduction the output file name is "gxtfbrsn"+telluric, or:              #
# gxtbrsnN20060210S0142.  The file prefixes are described below.          #
#                                                                         #
# g = gemcombined/gemarithed   n=nfprepared  s=skysubtracted              #
# r=nsreduced  b = bad pixel corrected  f= run through nsfitcoords        # 
# t = nstransformed   x = extracted to a 1D spectrum                      #
#                                                                         #
# This script is meant to be a guideline as a method of a typical data    #
# reduction for NIFS frames.  Of course, NIFS PIs can add or skip steps   #
# in this reduction as they deem fit in order to reduce their particular  #
# datasets.                                                               #
#                                                                         #
###########################################################################
