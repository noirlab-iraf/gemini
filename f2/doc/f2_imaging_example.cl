# Copyright(c) 2010-2013 Association of Universities for Research in Astronomy, Inc.

###############################################################################
# Gemini F2 example data reduction script                                     #
# Typical reduction for: Imaging Science and Calibration Data                 #
#                                                                             #
# This script is provided to guide the user through the F2 data reduction     #
# process and may not be optimised to give the best results. It shows the     #
# reduction steps using F2 data and provides explanatory comments at each     #
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
# Science Dataset:
#
# Observation UT date : 2011 Dec 11
# Data filename prefix: S20111211S
# File numbers:
#       Object (NGC 2442)          : 297-306 (J-band, 20s)
#       Sky                        : 307-316 (J-band, 20s)
#       Darks for the object       : 642-652 (J-band, 20s)
#       Short darks for the BPM    : 576-586 (J-band, 2s)
#
# Observation UT date : 2011 Dec 12
# Data filename prefix: S20111212S
# File numbers:
#       Twilight flats             : 009-012 (J-band, 12s)
#       Darks for the flats        : 379-389 (J-band, 12s)

###############################################################################
# STEP 1: Initialize the required packages                                    #
###############################################################################

# Load the required packages
gemini
f2

# If copying and pasting these commands into an interactive PyRAF session, use
# the following lines to import the required packages
#from pyraf.iraf import gemini
#from pyraf.iraf import f2

# Use the default parameters except where specified on command lines below
print ("\nEXAMPLE: Unlearning tasks")
unlearn ("gemini")
unlearn ("f2")
unlearn ("niri")
unlearn ("gnirs")
unlearn ("gemtools")

###############################################################################
# STEP 2: Define any variables, the database and the logfile                  #
###############################################################################

# Define any variables (not required if copying and pasting into an interactive
# PyRAF session)
string rawdir, image
struct *scanfile

# Define the logfile
f2.logfile = "f2_imaging_example.log"

# To start from scratch, delete the existing logfile
printf ("EXAMPLE: Deleting %s\n", f2.logfile)
delete (f2.logfile, verify=no)

# Define the directory where the raw data is located
# Don't forget the trailing slash!
rawdir = "./"
printf ("EXAMPLE: Raw data is located in %s\n", rawdir)

# Load the header keywords for F2
nsheaders ("f2", logfile=f2.logfile)

# Set the display
set stdimage=imt2048

###############################################################################
# STEP 3: Create the reduction lists                                          #
###############################################################################

delete ("obj.lis,sky.lis,flats.lis,darks.lis,flatdarks.lis,shortdarks.lis,\
    calib.lis,all.lis", verify=no)

# The user should edit the parameter values in the gemlist calls below to match
# their own dataset.
print ("EXAMPLE: Creating the reduction lists")
gemlist "S20111211S" "297-306" > "obj.lis"
gemlist "S20111211S" "307-316" > "sky.lis"
gemlist "S20111211S" "642-652" > "darks.lis"
gemlist "S20111212S" "009-012" > "flats.lis"
gemlist "S20111212S" "379-389" > "flatdarks.lis"
gemlist "S20111211S" "576-586" > "shortdarks.lis"
concat ("flats.lis,darks.lis,flatdarks.lis,shortdarks.lis", "calib.lis")
concat ("calib.lis,obj.lis,sky.lis", "all.lis")

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
#     ...    iraf.display(rawdir + image, 1)
#     ...    iraf.sleep(5)
#     ...
#     >>> file.close()

scanfile = "all.lis"
while (fscan(scanfile, image) != EOF) {
    display (rawdir//image//"[1]", 1)
    sleep 5
}
scanfile = ""

###############################################################################
# STEP 5: f2prepare the calibration data                                      #
###############################################################################

# Run F2PREPARE first on the calibration data to update the headers, derive 
# variance and data quality (DQ) planes, correct for non-linearity (not yet 
# implemented) and flag saturated and non-linear pixels in the DQ plane.

imdelete ("f@calib.lis", verify=no)
f2prepare ("@calib.lis", rawpath=rawdir, fl_vardq=yes, fl_correct=yes, \
    fl_saturated=yes, fl_nonlinear=yes)

###############################################################################
# STEP 6: Create the normalised flat field and BPM                            #
###############################################################################

# Construct the normalised flat field using flat field images and short dark 
# images to identify bad pixels. The imaging flats are derived from darks 
# subtracted from equal exposure twilight flats.

# nsheaders sets the following default values for niflat:
#
#     niflat.thresh_flo = 0.70
#     niflat.thresh_fup = 1.20
#     niflat.thresh_dlo = -50.
#     niflat.thresh_dup = 600.
#     niflat.statsec = "[300:1748,300:1748]"
#
# The thresh_flo and thresh_fup values may need adjusting depending on the
# filter used. For this reason, it is recommended that niflat be run
# interactively. For the J-band flat data in this example, these default values
# produce a reasonable result.

imdelete ("flat.fits,f2_bpm.pl", verify=no)
delete ("fflats.lis,fflatdarks.lis,fshortdarks.lis", verify=no)
sections "f@flats.lis" > "fflats.lis"
sections "f@flatdarks.lis" > "fflatdarks.lis"
sections "f@shortdarks.lis" > "fshortdarks.lis"
niflat ("@fflats.lis", flatfile="flat.fits", lampsoff="@fflatdarks.lis", \
    darks="@fshortdarks.lis", bpmfile="f2_bpm.pl", fl_inter=yes)

# For twilight flats, niflat.fl_rmstars should be set to yes and niflat.scale
# should be set appropriately. Please see the help file for niflat for more
# details.

#niflat ("@fflats.lis", flatfile="flat.fits", lampsoff="@fflatdarks.lis", \
#    darks="@fshortdarks.lis", bpmfile="f2_bpm.pl", scale="median", \
#    fl_rmstars=yes, fl_inter=yes)

# Visually inspect the BPM. If inappropriate thresholds were used, the BPM will
# be poor (e.g., too many good pixels may be flagged), which can cause problems
# later on in the processing.

display ("f2_bpm.pl", 1)

###############################################################################
# STEP 7: f2prepare the science data                                          #
###############################################################################

# Run F2PREPARE on the science data to update the headers, derive variance and 
# data quality (DQ) planes, correct for non-linearity (not yet implemented) 
# and flag saturated, non-linear and bad pixels in the DQ plane.

imdelete ("f@obj.lis", verify=no)
f2prepare ("@obj.lis", rawpath=rawdir, bpm="f2_bpm.pl", fl_vardq=yes, \
    fl_correct=yes, fl_saturated=yes, fl_nonlinear=yes)

imdelete ("f@sky.lis", verify=no)
f2prepare ("@sky.lis", rawpath=rawdir, bpm="f2_bpm.pl", fl_vardq=yes, \
    fl_correct=yes, fl_saturated=yes, fl_nonlinear=yes)

###############################################################################
# STEP 8: Create the sky frame                                                #
###############################################################################

# Construct a sky frame by identifying objects in each image, removing them, 
# and averaging the remaining good pixels. The object masks created by nisky 
# are saved so that they can be checked to ensure that the task is masking 
# objects in the images appropriately.

# Create the dark frame for the sky images

imdelete ("dark.fits", verify=no)
delete ("fdarks.lis", verify=no)
sections "f@darks.lis" > "fdarks.lis"
gemcombine ("@fdarks.lis", "dark.fits", combine="average", fl_vardq=yes, \
    logfile=f2.logfile, fl_dqprop=yes)

# Dark subtract the prepared sky images
# nsheaders sets nireduce.statsec = "[300:1748,300:1748]"

delete ("fsky.lis", verify=no)
sections "f@sky.lis" > "fsky.lis"
imdelete ("df@sky.lis", verify=no)
nireduce ("@fsky.lis", outprefix="d", fl_sky=no, fl_autosky=no, \
    fl_scalesky=no, fl_dark=yes, darkimage="dark.fits", fl_flat=no)

# Create the sky frame using the dark subtracted sky images
# nsheaders sets nisky.statsec = "[300:1748,300:1748]"

imdelete ("sky.fits", verify=no)
delete ("dfskymsk.lis", verify=no)
sections "df@sky.lis//msk.pl" > "dfskymsk.lis"
imdelete ("@dfskymsk.lis", verify=no)
delete ("dfsky.lis", verify=no)
sections "df@sky.lis" > "dfsky.lis"
nisky ("@dfsky.lis", outimage="sky.fits", combtype="median", fl_keepmasks=yes)

# Flat divide the sky frame

imdelete ("fsky.fits", verify=no)
nireduce ("sky.fits", outprefix="f", fl_sky=no, fl_autosky=no, \
    fl_scalesky=no, fl_dark=no, fl_flat=yes, flatimage="flat.fits")

###############################################################################
# STEP 9: Reduce the science data                                             #
###############################################################################

# Reduce the raw science images by subtracting the dark frame, dividing by the 
# normalised flat field image and subtracting the sky frame.

# Since the sky images and science images have the same exposure time, the same
# dark frame can be used for both.

# Dark subtract the prepared science images
# nsheaders sets nireduce.statsec = "[300:1748,300:1748]"

delete ("fobj.lis", verify=no)
sections "f@obj.lis" > "fobj.lis"
imdelete ("df@obj.lis", verify=no)
nireduce ("@fobj.lis", outprefix="d", fl_sky=no, fl_autosky=no, \
    fl_scalesky=no, fl_dark=yes, darkimage="dark.fits", fl_flat=no)

# Flat divide the dark subtracted science images

delete ("dfobj.lis", verify=no)
sections "df@obj.lis" > "dfobj.lis"
imdelete ("fdf@obj.lis", verify=no)
nireduce ("@dfobj.lis", outprefix="f", fl_sky=no, fl_autosky=no, \
    fl_scalesky=no, fl_dark=no, fl_flat=yes, flatimage="flat.fits")

# Sky subtract the flat divided, dark subtracted science images

delete ("fdfobj.lis", verify=no)
sections "fdf@obj.lis" > "fdfobj.lis"
imdelete ("rfdf@obj.lis", verify=no)
nireduce ("@fdfobj.lis", outprefix="r", fl_sky=yes, skyimage="fsky.fits", \
    fl_autosky=yes, fl_scalesky=yes, fl_dark=no, fl_flat=no)

###############################################################################
# STEP 10: Combine the science data                                           #
###############################################################################

imdelete ("obj_comb.fits", verify=no)
delete ("rfdfobj.lis", verify=no)
sections "rfdf@obj.lis" > "rfdfobj.lis"
imcoadd ("@rfdfobj.lis", outimage="obj_comb.fits", rotate=no, \
    geofitgeom="shift", niter=1, statsec="[300:1748,300:1748]", \
    badpix="f2_bpm.pl", fl_fixpix=yes, fl_find=yes, fl_map=yes, fl_trn=yes,
    fl_med=yes, fl_add=yes, fl_avg=yes, fl_scale=yes, fl_overwrite=yes,
    logfile=f2.logfile)

###############################################################################
# STEP 11: Tidy up                                                            #
###############################################################################

delete ("obj.lis,sky.lis,flats.lis,darks.lis,flatdarks.lis,shortdarks.lis,\
    calib.lis,all.lis,fflats.lis,fflatdarks.lis,fshortdarks.lis,fdarks.lis,\
    fsky.lis,dfskymsk.lis,dfsky.lis,fobj.lis,dfobj.lis,fdfobj.lis,\
    rfdfobj.lis", verify=no)

###############################################################################
# Finished!                                                                   #
###############################################################################
