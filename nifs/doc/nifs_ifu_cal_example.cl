# Copyright(c) 2006-2012 Association of Universities for Research in Astronomy, Inc.
#
# Gemini NIFS data reduction script
# Reduction for:  GENERAL BASELINE CALIBRATIONS
#
#
# Processed data taken on : 2006Feb10
# Data processor          : tbeck
# Data reduction date     : 2006Mar23
#
# File rootname: N20060210S 
# Files 191-197,205,389-390 for Calibrations
#   195-197     Lamps on Flats
#   389-390     Ronchi Flats
#   192-194     Lamps off Flats
#       191     Arcs
#       205     Darks for arc

###########################################################################
# DATA REDUCTION HINT-   		                                  #
#                                                                         #
# For a general, automatic baseline calibration reduction:                #
# In STEP 2, edit the data directory keyword "raw_dir" to be your         #
# calibration directory. Also in STEP 2, in the gemlist calls, change the #
# rootname (N20060210S) and the file numbers for your flat fields, lamps  #
# off flats, arcs, arcdarks and Ronchi files.                             #
#                                                                         #
# Run the routine by cutting and pasting the reduction steps into your    # 
# IRAF session. Highlighting and copying the whole reduction at once      #
# seems to work fine.                                                     #
# (Defining and running the script as a reduction "task" doesn't work     #
# well because the loading of the gemini package in "STEP 1" overwrites   #
# the defined task!)                                                      #
#                                                                         #
# On our fast linux machines at Gemini, the calibration reduction for     # 
# one full wavelength setting takes about 15 minutes to run.              #
#                                                                         #
###########################################################################
# Current limitations:  The NIFS Baseline calibration reductions have     #
# not been well tested on data that was obtained with non-standard        #
# wavelength configurations.  (i.e., a different central wavelength       #
# setting than the default Z, J, H and K-band values of 1.05, 1.25, 1.65  #
# and 2.20 microns).                                                      #
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

string calflat, flatdark, arc, arcdark, ronchiflat,
       raw_data, clist, user_clobber, band

user_clobber=envget("clobber")
reset clobber=yes

# For a general reduction of baseline calibrations, change the below 
# 6 lines and then run the full script.

raw_data = "/net/archie/staging/permi/"

gemlist N20060210S   195-197  > flatlist
gemlist N20060210S   192-194  > flatdarklist
gemlist N20060210S       191  > arclist
gemlist N20060210S       205  > arcdarklist
gemlist N20060210S   389-390  > ronchilist

cat flatlist     | scan(calflat)    # catch only the first file
cat flatdarklist | scan(flatdark)   # catch only the first file
cat arclist      | scan(arc)        # catch only the first file
cat arcdarklist  | scan(arcdark)    # catch only the first file
cat ronchilist   | scan(ronchiflat) # catch only the first file

###########################################################################
# STEP 3:  Determine the shift to the MDF file                            #
###########################################################################

nfprepare(calflat,rawpath=raw_data,outpref="s", shiftx=INDEF,  
          shifty=INDEF,fl_vardq-,fl_corr-,fl_nonl-)

###########################################################################
# STEP 4:  Make the Flat field and BPM                                    #
###########################################################################

nfprepare("@flatlist",    
   rawpath=raw_data,shiftim="s"//calflat,fl_vardq+,fl_int+,fl_corr-,fl_nonl-)
nfprepare("@flatdarklist",
   rawpath=raw_data,shiftim="s"//calflat,fl_vardq+,fl_int+,fl_corr-,fl_nonl-)

gemcombine("n//@flatlist",output="gn"//calflat,
   fl_dqpr+,fl_vardq+,masktype="none",logfile="nifs.log")
gemcombine("n//@flatdarklist",output="gn"//flatdark,
   fl_dqpr+,fl_vardq+,masktype="none",logfile="nifs.log")

nsreduce ("gn"//calflat, 
   fl_cut+,fl_nsappw+,fl_vardq+,fl_sky-,fl_dark-,fl_flat-,logfile="nifs.log")
nsreduce ("gn"//flatdark, 
   fl_cut+,fl_nsappw+,fl_vardq+,fl_sky-,fl_dark-,fl_flat-,logfile="nifs.log")

# creating flat image, final name = rnN....._sflat.fits
nsflat("rgn"//calflat,darks="rgn"//flatdark,flatfile="rn"//calflat//"_sflat",
   darkfile="rn"//flatdark//"_dark",fl_save_dark+,process="fit",
   thr_flo=0.15,thr_fup=1.55,fl_vardq+,logfile="nifs.log") 

#rectify the flat for slit function differences - make the final flat.
nsslitfunction("rgn"//calflat,"rn"//calflat//"_flat",
   flat="rn"//calflat//"_sflat",dark="rn"//flatdark//"_dark",combine="median",
   order=3,fl_vary-,logfile="nifs.log")

###########################################################################
# STEP 5:  Reduce the Arc and determine the wavelength solution           #
###########################################################################

nfprepare("@arclist",    
   rawpath=raw_data, shiftimage="s"//calflat,
   bpm="rn"//calflat//"_sflat_bpm.pl",fl_vardq+,fl_corr-,fl_nonl-)
nfprepare("@arcdarklist",
   rawpath=raw_data, shiftimage="s"//calflat,
   bpm="rn"//calflat//"_sflat_bpm.pl",fl_vardq+,fl_corr-,fl_nonl-)

# Determine the number of input arcs and arc darks so that the
# routine runs automatically for single or multiple files.
int nfiles
count("arclist") | scanf("%d",nfiles)
if (nfiles > 1) 
   gemcombine("n//@arclist",output="gn"//arc,
      fl_dqpr+,fl_vardq+,masktype="none",logfile="nifs.log")
else
   copy("n"//arc//".fits","gn"//arc//".fits")

count("arcdarklist") | scanf("%d",nfiles)
if (nfiles > 1)
   gemcombine("n//@arcdarklist",output="gn"//arcdark,
      fl_dqpr+,fl_vardq+,masktype="none",logfile="nifs.log")
else
   copy("n"//arcdark//".fits","gn"//arcdark//".fits")

nsreduce("gn"//arc,outpr="r",darki="gn"//arcdark,flati="rn"//calflat//"_flat",
   fl_vardq-, fl_cut+, fl_nsappw+, fl_sky-, fl_dark+,fl_flat+, 
   logfile="nifs.log")

###########################################################################
#  DATA REDUCTION HINT -                                                  # 
# For the nswavelength call, the different wavelength settings            #
# use different vaues for some of the parameters. For optimal auto        #
# results, use:                                                           #
#                                                                         #
# K-band: thresho=50.0, cradius=8.0   -->  (gives rms of 0.1 to 0.3)      # 
# H-band: thresho=100.0, cradius=8.0  -->  (gives rms of 0.05 to 0.15)    #
# J-band: thresho=100.0               -->  (gives rms of 0.03 to 0.09)    # 
# Z-band: Currently not working very well for non-interactive mode        #
#                                                                         # 
# Note that better RMS fits can be obtained by running the wavelength     #
# calibration interactively and identifying all of the lines              #
# manually.  Tedious, but will give more accurate results than the        #
# automatic mode (i.e., fl_inter-).  Use fl_iner+ for manual mode.        #
#                                                                         # 
###########################################################################

# Determine the wavelength of the observation and set the arc coordinate
# file.  If the user wishes to change the coordinate file to a different
# one, they need only to change the "clist" variable to their line list
# in the coordli= parameter in the nswavelength call.

real my_thresh
imgets("rgn"//arc//"[0]","GRATING")
band=substr(imgets.value,1,1)
if (band == "Z") {
     clist="nifs$data/ArXe_Z.dat"
     my_thresh=100.0
} else if (band == "K") {
     clist="nifs$data/ArXe_K.dat"
     my_thresh=50.0
} else {
     clist="gnirs$data/argon.dat"
     my_thresh=100.0
}

nswavelength("rgn"//arc, coordli=clist, nsum=10, thresho=my_thresh, 
   trace=yes,fwidth=2.0,match=-6,cradius=8.0,fl_inter-,nfound=10,nlost=10,
   logfile="nifs.log")

###########################################################################
# STEP 6:                                                                 #
#  Trace the spatial curvature and spectral distortion in the Ronchi flat #
###########################################################################
nfprepare("@ronchilist",rawpath=raw_data, shiftimage="s"//calflat,
   bpm="rn"//calflat//"_sflat_bpm.pl", fl_vardq+,fl_corr-,fl_nonl-)

# Determine the number of input Ronchi calibration mask files so that
# the routine runs automatically for single or multiple files.
count("ronchilist") | scan(nfiles)
if (nfiles > 1)
   gemcombine("n//@ronchilist",output="gn"//ronchiflat,fl_dqpr+,
           masktype="none",fl_vardq+,logfile="nifs.log")
else
   copy("n"//ronchiflat//".fits","gn"//ronchiflat//".fits")

nsreduce("gn"//ronchiflat, outpref="r", 
   dark="rn"//flatdark//"_dark",
   flatimage="rn"//calflat//"_flat",
   fl_cut+, fl_nsappw+, fl_flat+, fl_sky-, fl_dark+, fl_vardq-,
   logfile="nifs.log")

nfsdist("rgn"//ronchiflat, 
   fwidth=6.0, cradius=8.0, glshift=2.8, 
   minsep=6.5, thresh=2000.0, nlost=3, 
   fl_inter-,logfile="nifs.log")

###########################################################################
# Reset to user defaults                                                  #
###########################################################################
if (user_clobber == "no")
   set clobber = no
;

###########################################################################
#		End of the Baseline Calibration reduction                 #
###########################################################################
#	                                                                  #
#  The final output files created from this script for later science      #
#  reduction have prefixes and file names of:                             #
#     1. Shift reference file:  "s"+calflat                               #
#     2. Flat field:  "rn"+calflat+"_flat"                                #
#     3. Flat BPM (for DQ plane generation):  "rn"+calflat+"_flat_bpm.pl" #
#     4. Wavelength referenced Arc:  "wrn"+arc                            #
#     5. Spatially referenced Ronchi Flat:  "rn"+ronchiflat               #
#     For this reduction,                                                 #
#        Shift ref. file =   sN20060210S0324.fits                         #
#        Flat field      =  rnN20060210S0324_flat.fits                    #
#        Flat BPM        =  rnN20060210S0324_sflat_bpm.pl                 #
#        Arc frame       =  wrnN20060210S0329.fits                        #
#        Ronchi flat     =  rnN20060210S0331.fits                         #
#	                                                                  #
#  NOTE:  Other important information for reducing the science data is    #
#    included in the "database" directory that is created and edited      #
#    within the above "nswavelength" and "nfsdist" IRAF calls. For a      #
#    proper science reduction to work (particularly the "nsfitcoords"     #
#    step), the science data must either be reduced in the same directory #
#    as the calibrations, or the whole "database" directory created by    #
#    this script must be copied into the working science reduction        #
#    directory.                                                           #
#                                                                         #
###########################################################################
