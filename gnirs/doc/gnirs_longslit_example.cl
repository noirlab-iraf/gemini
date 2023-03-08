# Copyright(c) 2004-2015 Association of Universities for Research in Astronomy, Inc.

###############################################################################
# Gemini GNIRS example data reduction script                                  #
# Typical reduction for: Longslit Science and Calibration Data                #
#                                                                             #
# This script is provided to guide the user through the GNIRS data reduction  #
# process and may not be optimised to give the best results. It shows the     #
# reduction steps using GNIRS data and provides explanatory comments at each  #
# step.                                                                       #
#                                                                             #
# It is strongly recommended that the user read the embedded comments and     #
# understand the processing steps, since the optimum steps for a given        #
# dataset may differ, e.g., improved cleaning of cosmic rays and bad pixels   #
# and improved signal-to-noise will most likely be possible. The user can     #
# then edit this script to match their dataset.                               #
#                                                                             #
# Once this script has been properly edited, it can be run by copying and     #
# pasting each command into an IRAF or PyRAF session, or by defining this     #
# script as a task by typing:                                                 #
#                                                                             #
#   ecl> task $thisscript="thisscript.cl"                                     #
#   ecl> thisscript                                                           #
#                                                                             #
# in the IRAF or PyRAF session. It is NOT recommended to run this script      #
# using redirection i.e., cl < thisscript.cl.                                 #
#                                                                             #
# Note that this script is designed to be re-run as needed, so each step is   #
# preceded by a command to delete the files created in that step.             #
###############################################################################

# The data files have been separated by filter and exposure time, where
# appropriate. This information can be found in the primary header unit (PHU)
# of each data file. The imhead or fitsutil.fxhead tasks can be used to view
# the header information in the PHU (or any other extension) of the data file. 
# The hselect task can be used to obtain specific keyword values from the
# headers. Read the help files for these tasks for more information.
#
# Observation UT date : 2004 Jun 15
# Data filename prefix: S20040615S 
#
# File numbers:
#       Science target : 156-183  (GS-2004A-SV-8)
#       Telluric       : 145-152
#       IR lamp flats  : 184-194
#       Argon arcs     : 195-196
#
# This is an example of a basic calibration and science longslit reduction for
# moderate resolution data (R~1700) of an extended source with central
# wavelength=2.2microns. The calibration data consist of GCAL flats, an arc
# and a telluric standard. The science data were taken with an obj-sky dither
# pattern, i.e., only 50% of data are taken on-source (as opposed to dithering
# along the slit).
#
###############################################################################
#
# This script can also be used as a guide for reduction of longslit/111
# (R~6000) data. The reduction steps are identical; what differences 
# there are are in parameter settings and these are identified in comments
# starting with "LS111:".
#
# LS111 Dataset:
#
# Observation UT date : 2004 Jun 13
# Data filename prefix: S20040613S
#
# File numbers:
#       Science target : 124-155  (GS-2004A-SV-15)
#       Telluric       : 161-168
#       IR lamp flats  : 171-180
#       Argon arcs     : 182
#       Pinhole        : 183
#
# To run the LS111 case, modify the gemlist commands below for these datasets.

###############################################################################
# STEP 1: Initialize the required packages                                    #
###############################################################################

# Load the required packages
gemini
gnirs

# If copying and pasting these commands into an interactive PyRAF session, use
# the following lines to import the required packages
#from pyraf.iraf import gemini
#from pyraf.iraf import gnirs

# Use the default parameters except where specified on command lines below
print ("EXAMPLE: Unlearning tasks")
unlearn ("gemini")
unlearn ("gemtools")
unlearn ("gnirs")

###############################################################################
# STEP 2: Define any variables and the logfile                                #
###############################################################################

# Define any variables (not required if copying and pasting into an interactive
# PyRAF session)
string rawdir, image
struct *scanfile

# Define the logfile
gnirs.logfile = "gnirs_longslit_example.log"

# To start from scratch, delete the existing logfile
printf ("EXAMPLE: Deleting %s\n", gnirs.logfile)
delete (gnirs.logfile, verify=no)

# Define the database directory
gnirs.database = "gnirs_longslit_database/"

# To start from scratch, delete the existing database files

# If copying and pasting these commands into an interactive PyRAF session, use
# something similar to the following example instead of using the uncommented
# lines below.
#
#     >>> if (iraf.access(gnirs.database)):
#     ...    print "EXAMPLE: Deleting contents of %s" % (gnirs.database)
#     ...    iraf.delete (gnirs.database + "*", verify=no)

if (access(gnirs.database)) {
    printf ("EXAMPLE: Deleting contents of %s\n", gnirs.database)
    delete (gnirs.database//"*", verify=no)
}
;

# Define the directory where the raw data is located
rawdir = "./"
printf ("EXAMPLE: Raw data is located in %s\n", rawdir)

# Load the header keywords for GNIRS
nsheaders ("gnirs")

# Set the display
set stdimage=imt1024

###############################################################################
# STEP 3: Create the reduction lists                                          #
###############################################################################

delete ("obj.lis,sky.lis,telluric.lis,flats.lis,arc.lis,all.lis", verify=no)

# The user should edit the parameter values in the gemlist calls below to match
# their own dataset.
print ("EXAMPLE: Creating the reduction lists")
gemlist "S20040615S" "156,159,160,163,164,167,168,171,172,175,176,179,180,183"\
    > "obj.lis"
gemlist "S20040615S" "157,158,161,162,165,166,169,170,173,174,177,178,181,182"\
    > "sky.lis"
gemlist "S20040615S" "145-152" > "telluric.lis"
gemlist "S20040615S" "184-194" > "flats.lis"
gemlist "S20040615S" "195" > "arc.lis"
concat ("obj.lis,sky.lis,telluric.lis,flats.lis,arc.lis", "all.lis")

###############################################################################
# STEP 4: Visually inspect the data                                           #
###############################################################################

# Visually inspect all the data. In addition, all data should be visually
# inspected after every processing step. Once the data has been prepared, it is
# recommended to use the syntax [EXTNAME,EXTVER] e.g., [SCI,1], when defining
# the extension.

# Please make sure a display tool (e.g., ds9, ximtool) is already open.

# If copying and pasting these commands into an interactive PyRAF session, use
# something similar to the following example instead of using the uncommented
# lines below.
#
#     >>> file = open("all.lis", "r")
#     >>> for line in file:
#     ...    image = line.strip() + "[1]"
#     ...    iraf.display(image, 1)
#     ...    iraf.sleep(5)
#     >>> file.close()

scanfile = "all.lis"
while (fscan(scanfile, image) != EOF) {
    display (image//"[1]", 1)
    sleep 5
}
scanfile = ""

###############################################################################
# STEP 5: Prepare the data                                                    #
###############################################################################

# If the data show a strong vertical striping pattern (due to variable bias
# problem on readout), run nvnoise on the raw data to correct this. This task
# subtracts a constant value from each pixel and so does not add any noise.
#
# If nvnoise is run, all subsequent files will have a 'v' preceding the file
# prefix (by default).
#
#imdelete ("v@all.lis", verify=no)
#nvnoise ("@all.lis")

# Run nsprepare to correct for non-linearity, determine MDF shifts, fix headers
# and create the variance and data quality planes. Use fl_forcewcs=yes to
# correct the WCS keywords in the headers. 

# At the time of writing, it is not recommend to use the nonlinearity
# correction provided by nsprepare on data taken with GNIRS at Gemini North. In
# this case, specify fl_correct=no in the call to nsprepare below.

# It is best to nsprepare the science data with the calibration data so that a
# flat can be used to determine the MDF shift (shiftx=INDEF and shifty=INDEF). 
# If the science data is prepared separately, use the "shiftimage" parameter to
# match the shift to the flat data. It is important that the calibration and
# science data have the same shift.

# Specify the correct BPM
#     For GNIRS at Gemini South, bpm = "gnirs$data/gnirs_2005sep24_bpm.fits"
#     For GNIRS at Gemini North, bpm = "gnirs$data/gnirsn_2010oct12_bpm.fits"
#                             OR bpm = "gnirs$data/gnirsn_2011apr07_bpm.fits"
#                             OR bpm = "gnirs$data/gnirsn_2012dec05_bpm.fits"
#     (the default value is gnirs$data/gnirsn_2012dec05_bpm.fits)

imdelete ("n@all.lis", verify=no)
nsprepare ("@all.lis", rawpath=rawdir//"$", shiftx=INDEF, shifty=INDEF, \
    fl_forcewcs=yes, bpm="gnirs$data/gnirs_2005sep24_bpm.fits")

###############################################################################
# STEP 6: Generate the normalised flat                                        #
###############################################################################

# Print (to the screen and to the logfile) the statistics of the flat to check
# for anomalous mean or standard deviation values. Often the first flat in a
# set should be discarded. Note: these steps create and delete a file called
# "tmpflat"

printlog ("-------------------------------------------- ", gnirs.logfile, \
    verbose=yes)
delete ("tmpflat", verify=no)

gemextn "n@flats.lis" proc="expand" extname="SCI" extver="1" > "tmpflat"
printlog ("Flat statistics: ", gnirs.logfile, verbose=yes)
imstatistic "@tmpflat" | tee (gnirs.logfile, append=yes)
delete ("tmpflat", verify=no)

# The mean of the first two flats deviate from the mean of the rest of the
# flats, so exclude these flats when creating the normalised flat

delete ("flats.lis", verify=no)
gemlist "S20040615S" "186-194" > "flats.lis"

# It is not required to cut longslit data, but it avoids a large number of
# "bad" pixels from the unilluminated edges in the DQ array.

imdelete ("rn@flats.lis", verify=no)
nsreduce ("n@flats.lis", fl_sky=no, fl_cut=yes, fl_flat=no, fl_dark=no, \
    fl_nsappwave=no)

# No darks are needed for such short exposure times (only add noise)
imdelete ("final_flat.fits,final_flat_bpm.pl", verify=no)
nsflat ("rn@flats.lis", flatfile="final_flat.fits")

display ("final_flat.fits[sci,1]", 1)

###############################################################################
# STEP 7: Reduce the arcs                                                     #
###############################################################################

imdelete ("rn@arc.lis", verify=no)
nsreduce ("n@arc.lis", fl_sky=no, fl_flat=no)

###############################################################################
# STEP 8: Obtain the wavelength solution                                      #
###############################################################################

# It is recommended to run nswavelength interactively, even though for lower
# dispersion, this task often works well non-interactively.

# Interactive Instructions:
#     If there are very few lines in your arc, reduce order (default=4). 
#     Default threshold=100, may want to reduce for very faint lines (but run 
#     interactively to be sure not to identify on noise!)

# Specify the correct coord list
#     For GNIRS at Gemini North, coordlist="gnirs$data/Ar_Xe.dat"
#     For GNIRS at Gemini South, coordlist="gnirs$data/lowresargon.dat" 
#         (this is the default value for the coordlist parameter)

imdelete ("wrn@arc.lis", verify=no)
nswavelength ("rn@arc.lis", coordlist="gnirs$data/lowresargon.dat", \
    fl_inter=yes)

# For LS111: 
#     Need to increase cradius for 111 resolution (cradius=20).
#nswavelength ("rn@arc.lis", cradius=20., threshold=50., order=2, fl_inter=yes)

###############################################################################
# STEP 9: Sky subtract the telluric and science data                          #
###############################################################################

# There are 3 ways to determine the correct frames for sky subtraction in
# nsreduce. If the frames are evenly spaced in time (i.e., ABBA dither
# pattern), the easiest is "skyrange=INDEF" (the default), which will determine
# the best sky frame for each image from the input list. One can also set the
# skyrange (in seconds) by hand, or provide a list of sky images directly. See
# nsreduce help for more information. If two or more images meet the sky
# selection criteria (skyrange and nodsize), they will be combined to create
# the sky frame that is subtracted.

# The following steps can be used to determine the spacing in time between
# images, if needed to set skyrange.

delete ("tmpsky", verify=no)
printlog ("--------------------------------------------", gnirs.logfile, \
    verbose=yes)
gemextn "n@telluric.lis" proc="expand" index="0" > "tmpsky"
printlog ("Telluric Exposure Times: ", gnirs.logfile, verbose=yes)
hselect "@tmpsky" "$I,UT" yes | tee (gnirs.logfile, append=yes)
printlog ("--------------------------------------------", gnirs.logfile, \
    verbose=yes)
delete ("tmpsky", verify=no)

# ~20 seconds will include neighbours for sky; use this value for the skyrange
# parameter.

imdelete ("rn@telluric.lis", verify=no)
nsreduce ("n@telluric.lis", fl_nsappwave=no, fl_sky=yes, skyrange=20, \
    fl_flat=yes, flatimage="final_flat.fits")

# For LS111 data:
#nsreduce ("n@telluric.lis", fl_nsappwave=no, fl_sky=yes, skyrange=25, \
#    fl_flat=yes, flatimage="final_flat.fits")

delete ("tmpsky", verify=no)
printlog ("--------------------------------------------", gnirs.logfile, \
    verbose=yes)
gemextn "n@obj.lis" proc="expand" index="0" > "tmpsky"
printlog ("Science Exposure Times: ", gnirs.logfile, verbose=yes)
hselect "@tmpsky" "$I,UT" yes | tee (gnirs.logfile, append=yes)
printlog ("--------------------------------------------", gnirs.logfile, \
    verbose=yes)
delete ("tmpsky", verify=no)

imdelete ("rn@obj.lis", verify=no)
nsreduce ("n@obj.lis", fl_nsappwave=no, fl_sky=yes, skyimages="n@sky.lis", \
    skyrange=180, fl_flat=yes, flatimage="final_flat.fits")

###############################################################################
# STEP 10: Apply the wavelength solution to the telluric and science data     #
###############################################################################

# The nsfitcoords task is used to determine the final solution (consisting of
# the wavelength solution) to be applied to the data. The nstransform task is
# used to apply this final solution. nsfitcoords is best run interactively.

# IMPORTANT: be sure to apply the same solution for the telluric and the
#            science data. 

# Spatial rectification (s-distortion correction) is not usually needed with 
# GNIRS long slit data. If it is desired, first call nssdist (see XD example).

imdelete ("frn@telluric.lis", verify=no)
nsfitcoords ("rn@telluric.lis", lamp="wrn@arc.lis")

imdelete ("tfrn@telluric.lis", verify=no)
nstransform ("frn@telluric.lis")

imdelete ("frn@obj.lis", verify=no)
nsfitcoords ("rn@obj.lis", lamp="wrn@arc.lis")

imdelete ("tfrn@obj.lis", verify=no)
nstransform ("frn@obj.lis")

###############################################################################
# STEP 11: Combine the telluric and science data                              #
###############################################################################

# FYI:
#
# NSSSTACK combines things at the same offset position (no shifting) and gives
#     1 output file per position. 
# NSCOMBINE shifts all data spatially to a common offset position and gives 1
#     combined output file.

imdelete ("tell_comb.fits", verify=no)
nscombine ("tfrn@telluric.lis", output="tell_comb")

display ("tell_comb.fits", 1)

imdelete ("obj_comb.fits", verify=no)
nscombine ("tfrn@obj.lis", output="obj_comb")

display ("obj_comb.fits", 1)

###############################################################################
# STEP 12: Extract the telluric and science data                              #
###############################################################################

imdelete ("xtell_comb.fits", verify=no)
nsextract ("tell_comb.fits")

# Final s/n ~50
splot ("xtell_comb.fits[sci,1]")

imdelete ("xobj_comb.fits", verify=no)
nsextract ("obj_comb.fits")

# Final s/n ~6
splot ("xobj_comb.fits[sci,1]")

###############################################################################
# STEP 13: Apply the telluric correction to the science data                  #
###############################################################################

# Note that this telluric has not been corrected to remove intrinsic stellar
# features; this will leave a false emission feature in the final spectrum at
# the position of the telluric Br-gamma line.

# Interactive execution of nstelluric may produce better results.

imdelete ("axobj_comb.fits", verify=no)
nstelluric ("xobj_comb.fits", "xtell_comb", thresh=0.01, high=5, fitord=15) 

splot ("axobj_comb.fits[sci,1]")
specplot ("xobj_comb.fits[sci,1],axobj_comb.fits[sci,1],\
    xtell_comb.fits[sci,1]")

###############################################################################
# FINISHED!                                                                   #
###############################################################################
