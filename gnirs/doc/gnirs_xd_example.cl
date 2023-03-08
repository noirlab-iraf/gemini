# Copyright(c) 2004-2015 Association of Universities for Research in Astronomy, Inc.

###############################################################################
# Gemini GNIRS example data reduction script                                  #
# Typical reduction for: Cross-Dispersed Science and Calibration Data         #
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
# Science Dataset:
#
# Observation UT date : 2006 May 07
# Data filename prefix: S20060507S
# File numbers:
#       Object        : 062-069
#       Telluric      : 075-082
#       IR lamp flats : 128-136
#       QH lamp flats : 138-146
#       Argon arcs    : 070-071
#       Pinhole       : 125

# Further information on reducing GNIRS Cross-Dispersed data can be found at
# http://www.gemini.edu/sciops/instruments/gnirs/data-format-and-reduction/reducing-xd-spectra
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
unlearn ("gnirs")

###############################################################################
# STEP 2: Define any variables and the logfile                                #
###############################################################################

# Define any variables (not required if copying and pasting into an interactive
# PyRAF session)
string rawdir, image
struct *scanfile

# Define the logfile
gnirs.logfile = "gnirs_xd_example.log"

# To start from scratch, delete the existing logfile
printf ("EXAMPLE: Deleting %s\n", gnirs.logfile)
delete (gnirs.logfile, verify=no)

# Define the database directory
gnirs.database = "gnirs_xd_database/"

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

delete ("obj.lis,telluric.lis,flat1.lis,flat2.lis,arc.lis,pinhole.lis,\
    allflats.lis,all.lis", verify=no)

# This example dataset uses the GCAL NIR balance filter for the flats. 
# Therefore, there are only two sets of flats. Flats taken before ~April 2006
# or any flats using the GCAL ND filters will have three sets of flats. For
# example, for older flats taken with GNIRS at Gemini South, there will be two
# sets of QH lamp flats, where the short exposure QH lamps are for orders 4-5
# and the longer exposure QH lamps are for orders 6-8.

# The user should edit the parameter values in the gemlist calls below to match
# their own dataset.
print ("EXAMPLE: Creating the reduction lists")
gemlist "S20060507S" "062-069" > "obj.lis"
gemlist "S20060507S" "075-082" > "telluric.lis"
# IR lamp: for order 3 
gemlist "S20060507S" "128-136" > "flat1.lis"
# QH lamp: for orders 4-8
gemlist "S20060507S" "138-146" > "flat2.lis"
gemlist "S20060507S" "070-071" > "arc.lis"
gemlist "S20060507S" "125" > "pinhole.lis"
concat ("flat1.lis,flat2.lis", "allflats.lis")
concat ("allflats.lis,obj.lis,telluric.lis,arc.lis,pinhole.lis", "all.lis")

###############################################################################
# STEP 4: Visually inspect the data                                           #
###############################################################################

# Visually inspect all the data. In addition, all data should be visually
# inspected after every processing step. Once the data have been prepared, it
# is recommended to use the syntax [EXTNAME,EXTVER] e.g., [SCI,1], when
# defining the extension.

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

# Note that nsprepare looks for a well exposed order 3 to obtain the shift, so
# shiftx=INDEF and shifty=INDEF will work best with IR lamp flats (nominally
# the first flat set).

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

# For old flats obtained before ~April 2006, or any flats using GCAL ND filters
# flat1.lis contains the IR lamp flats, good for order 3 and possibly order 4
# flat2.lis contains the short exposure QH lamp flats, good for orders 4 and 5
# flat3.lis contains the longer exposure QH lamp flats, orders 3-5 saturated,
# good for orders 6-8 

# For newer flats using GCAL NIR Balance filter there are two sets of flats
# flat1.lis contains the IR lamp flats, good for order 3 (the QH lamp has a
# K-band absorption feature) 
# flat2.lis contains the QH lamp flats, good for orders 4-8

imdelete ("rn@allflats.lis", verify=no)
nsreduce ("n@allflats.lis", fl_sky=no, fl_cut=yes, fl_flat=no, fl_dark=no, \
    fl_nsappwave=no, fl_corner=yes)

# Check the results of the cutting procedure by displaying a flat using
# nxdisplay

scanfile = "allflats.lis"
while (fscan(scanfile, image) != EOF) {
    nxdisplay ("rn"//image, 1)
    sleep 5
}
scanfile = ""

# Print (to the screen and to the logfile) the statistics on a selection of the
# flat orders to check for anomalous mean or standard deviation values. Often
# the first flat in a set should be discarded. Note: these steps create and
# delete a file called "tmpflat"

printlog ("-------------------------------------------- ", gnirs.logfile, \
    verbose=yes)
delete ("tmpflat", verify=no)

gemextn "rn@flat1.lis" proc="expand" extname="SCI" extver="1" > "tmpflat"
printlog ("Order 3 Flats: ", gnirs.logfile, verbose=yes)
imstatistic "@tmpflat" | tee (gnirs.logfile, append=yes)
delete ("tmpflat", verify=no)

gemextn "rn@flat2.lis" proc="expand" extname="SCI" extver="3" > "tmpflat"
printlog ("Order 5 Flats: ", gnirs.logfile, verbose=yes)
imstatistic "@tmpflat" | tee (gnirs.logfile, append=yes)
delete ("tmpflat", verify=no)


gemextn "rn@flat2.lis" proc="expand" extname="SCI" extver="4" > "tmpflat"
# Old flats (3 sets):
#gemextn "rn@flat3.lis" proc="expand" extname="SCI" extver="4" > "tmpflat"
printlog ("Order 6 Flats: ", gnirs.logfile, verbose=yes)
imstatistic "@tmpflat" | tee (gnirs.logfile, append=yes)
delete ("tmpflat", verify=no)

# Intermediate flats are called flat1, flat2, flat3. The final flat is called
# final_flat.

# IMPORTANT: It is recommended to use process="fit" for XD data ("auto" also
#            defaults to "fit" as of GEMINI v1.9). process="fit" gives a much
#            better final flat ("trace" does not work adequately).

# To fit the normalization manually, set fl_inter=yes; to define the lower and
# upper thresholds for good pixels interactively, set fl_range=yes
#
# The S/N given in the far right column by nsflat is incorrect for XD flats. 
# The S/N given in the middle column (under "Lamps") is more accurate.
# 
# A warning that "DQ for flat is poor" is commonly seen with XD data. This is
# because all of the un-illuminated pixels are set to "bad". You can ignore
# this warning (but check that the final flat looks good and there is real
# signal within the illuminated regions of each order).

# flat1 (order 3): the following parameters should work well non-interactively.
# An order=10 fit is needed to fit the spectral signature of the lamp in the
# K-band. If the fitting is done interactively, all 8 orders will be displayed.
# However, only order 3 is used to create the final flat. Therefore, the other
# orders can be ignored (for orders > 4, there is almost no signal).

delete ("flat1.fits,flat1_bpm.pl", verify=no)
nsflat ("rn@flat1.lis", flatfile="flat1.fits", fl_inter=no, fl_corner=yes, \
    process="fit", fitsec="MDF", order=10, lthresh=100., thr_flo=0.35, \
    thr_fup=1.5)
display ("flat1.fits[sci,1]", 1, zr=yes, zs=yes)

# flat2 (order 4-5 for old flats, order 4-8 for new flats): for the higher
# orders, it is better to fit each curve interactively (use fl_inter=yes). A
# high value for the upper threshold limit (thr_fup=4.0) is needed to prevent
# many good pixels from exceeding the limit (the calculated mean is too
# low). Alternatively, set fl_range=yes to set the good pixel limits by hand.

delete ("flat2.fits,flat2_bpm.pl", verify=no)
nsflat ("rn@flat2.lis", flatfile="flat2.fits", fl_inter=yes, fl_corner=yes, \
    process="fit", fitsec="MDF", order=5, lthresh=50., thr_flo=0.35, \
    thr_fup=4.0)

# flat3 (orders 6-8 for old flats only)
#delete ("flat3.fits,flat3_bpm.pl", verify=no)
#nsflat ("rn@flat3.lis", flatfile="flat3.fits", fl_inter=yes, fl_corner=yes, \
#    process="fit", fitsec="MDF", order=5, lthresh=50., thr_flo=0.35, \
#    thr_fup=4.0)

# New flats
# Create the final flat: order: 3   -> flat1 -> position 1-3, (sci: 1)
#                               4-8 -> flat2 -> position 4-18 (sci: 2-6)

# Old flats 
# Create the final flat: order: 3   -> flat1 -> position 1-3,  (sci: 1)
#                               4-5 -> flat2 -> position 4-9,  (sci: 2-3)
#                               6-8 -> flat3 -> position 10-18 (sci: 4-6)

# The fxinsert commands below uses extension numbers, not names, and so assumes
# a certain MEF configuration. The assumed configuration is the one created by
# this example: 6 extracted orders with variance and data quality planes, and
# no attached MDF.  It is a good idea to use gemextn to verify the extensions
# in your generated flat files before using the fxinsert command to create the
# final flat. The output from gemextn should match the following. 
#
#gemextn ("flat1.fits")
#flat1[0]
#flat1[1][SCI,1]
#flat1[2][VAR,1]
#flat1[3][DQ,1]
#flat1[4][SCI,2]
#flat1[5][VAR,2]
#flat1[6][DQ,2]
#flat1[7][SCI,3]
#flat1[8][VAR,3]
#flat1[9][DQ,3]
#flat1[10][SCI,4]
#flat1[11][VAR,4]
#flat1[12][DQ,4]
#flat1[13][SCI,5]
#flat1[14][VAR,5]
#flat1[15][DQ,5]
#flat1[16][SCI,6]
#flat1[17][VAR,6]
#flat1[18][DQ,6]

# NOTE: the MDF files for XD data extract 6 orders (3-8). Order 9 is not
# extracted; to add order 9, modify the MDF file in gnirs$data or the MDF file
# attached to each data file in nsprepare.

imdelete ("final_flat.fits", verify=no)
fxcopy ("flat1.fits", "final_flat.fits", groups="0-3", new_file=yes)
fxinsert ("flat2.fits", "final_flat.fits[3]", groups="4-18")
# for old flats
#fxinsert ("flat2.fits", "final_flat.fits[3]", groups="4-9")
#fxinsert ("flat3.fits", "final_flat.fits[9]", groups="10-18")

nxdisplay ("final_flat.fits", 1)

# Use the lines below to set the pixel values to 1 in any order, if required.
#imreplace ("final_flat.fits[sci,5]", 1.)
#imreplace ("final_flat.fits[var,5]", 0.)
#imreplace ("final_flat.fits[dq,5]", 0.)
#
#imreplace ("final_flat.fits[sci,6]", 1.)
#imreplace ("final_flat.fits[var,6]", 0.)
#imreplace ("final_flat.fits[dq,6]", 0.)

###############################################################################
# STEP 7: Reduce the arcs and pinholes                                        #
###############################################################################

imdelete ("rn@arc.lis", verify=no)
nsreduce ("n@arc.lis", fl_cut=yes, fl_nsappwave=no, fl_dark=no, fl_sky=no, \
    fl_flat=no, fl_corner=yes)

imdelete ("rn@pinhole.lis", verify=no)
nsreduce ("n@pinhole.lis", fl_cut=yes, fl_nsappwave=no, fl_dark=no, \
    fl_sky=no, fl_flat=no, fl_corner=yes)

# Rename the pinhole to a generic name for use below (assumes only 1 pinhole
# image)
imdelete ("pinhole.fits", verify=no)
imrename ("rn@pinhole.lis", "pinhole")

###############################################################################
# STEP 8: Obtain the s-distortion solution                                    #
###############################################################################

# The s-distortion solution is used by nstransform later for spectral
# rectification. It can also be used by nswavelength if "fl_median=no". This
# command is best run interactively. 

# When presented with the cross section, hit "m" on the 5 strongest peaks and
# accept the pixel values given. Do not use any faint peaks to the left (i.e.,
# for data taken with GNIRS at Gemini South; it's partially vignetted by the
# decker). For each order, check if the offset between the peak and the listed
# value is consistent for all features, and if the automated finding does not
# associate the same value for two features. The offset will be different for
# different orders, and the way to override a value is to give a "better
# guess". nlost=0. is used to stop the task where the pinhole ends; this is
# needed to avoid "tracing into the noise" in the higher orders.

# Check the gnirs$data/README file and specify the correct pinhole file for the
# coordlist parameter
#page gnirs$data/README

nssdist ("pinhole", coordlist="gnirs$data/pinholes-short-dense.lis", \
    fl_inter=yes, function="legendre", order=5, minsep=5, thresh=1000, \
    nlost=0.)

###############################################################################
# STEP 9: Obtain the wavelength solution                                      #
###############################################################################

# Only one arc image is needed; if multiple arcs are taken, they can be
# combined to improve signal to noise. Using nscombine (even though there is no
# shift), because nsstack does not have an option to name output image.

imdelete ("arc_comb", verify=no)
nscombine ("rn@arc.lis", output="arc_comb")

# With this data set, the arcs were taken with a wide slit (0.675 arcsec),
# therefore the RMS are fairly large (~1-3angstroms). More accurate results can
# be obtained with a 2-pixel (0.3 arcsec) slit.

# Accurate line identification is usually obtained automatically for orders 
# 3-7; order 8 often needs to be corrected. All orders should be checked 
# interactively.

# If fl_median=yes, a reference arc for XD data is created by median filtering
# the arc in each order. This option is possible because the arc lines are
# approximately parallel to the rows. It also avoids calling nsfitcoords and
# nstransform in nswavelength to spatially rectify the arc. An alternative
# method is to spatially rectify the arc by using fl_median=no and setting
# sdist=pinhole, but be sure to apply the same spatial transformation here as
# is used later to rectify the telluric and science data. Setting fl_median=yes
# provides as good a wavelength solution as setting fl_median=no in the
# cross-dispersed data tested to date.

imdelete ("warc_comb", verify=no)
nswavelength ("arc_comb", coordlist="gnirs$data/lowresargon.dat", \
    fl_median=yes, fl_inter=yes, threshold=300., nlost=10, fwidth=5.)

###############################################################################
# STEP 10: Sky subtract the telluric and science data                         #
###############################################################################

# Data taken before July 2005 are littered with radiation events. For these
# data, cosmic ray rejection is needed earlier in the reduction process. Data
# taken after July 2005 will look much better as they do not suffer from
# excessive radiation events.

# There are 3 ways to determine the correct frames for sky subtraction in
# nsreduce. If the frames are evenly spaced in time (i.e., ABBA dither
# pattern), the easiest is "skyrange=INDEF" (the default), which will determine
# the best sky frame for each image from the input list. One can also set the
# skyrange (in seconds) by hand, or provide a list of sky images directly. See
# nsreduce help for more information. If two or more images meet the sky
# selection criteria (skyrange and nodsize), they will be combined to create
# the sky frame that is subtracted.

imdelete ("rn@telluric.lis", verify=no)
nsreduce ("n@telluric.lis", fl_corner=yes, fl_nsappwave=no, fl_sky=yes, \
    skyrange=INDEF, fl_flat=yes, flatimage="final_flat.fits")

imdelete ("rn@obj.lis", verify=no)
nsreduce ("n@obj.lis", fl_corner=yes, fl_nsappwave=no, fl_sky=yes, \
    skyrange=INDEF, fl_flat=yes, flatimage="final_flat.fits", nodsize=3.0

###############################################################################
# STEP 11: Combine the telluric and science data                              #
###############################################################################

imdelete ("tell_comb.fits", verify=no)
nscombine ("rn@telluric.lis", output="tell_comb")

nxdisplay ("tell_comb.fits", 1)

imdelete ("obj_comb.fits", verify=no)
nscombine ("rn@obj.lis", output="obj_comb")

nxdisplay ("obj_comb.fits", 1)

###############################################################################
# STEP 12: Apply s-distortion and wavelength solutions to the telluric and    #
#          science data                                                       #
###############################################################################

# The nsfitcoords task is used to determine the final solution (consisting of
# either the s-distortion solution, the wavelength solution, or both) to be
# applied to the data. The nstransform task is used to apply this final
# solution. nsfitcoords is best run interactively. If only the wavelength
# solution is applied, the spectra will not be "straightened", but they will be
# flipped in y. To also straighten, add the pinhole (sdisttrans="pinhole").

# The interactive part of nsfitcoords uses fitcoords, which plots the fit and
# the residuals for the 2-dimensional xy fit for each transformation
# (wavelength solution = lamp and s-distortion solution = sdist). See help
# pages for nsfitcoords, nstransform, fitcoords and transform for more
# information. For the wavelength (lamp) solution, xorder=2 and yorder=3 works
# best, and outlying points should be deleted. For the s-distortion (sdist;
# pinhole) solution, use xorder=4 and yorder=4 for XD32 data (lower orders
# e.g., 2-3 are fine for XD111 data). Usually the first row of y points and the
# last few rows of y points in each order should be deleted. Hit 'f' to redo
# the fit. If the pinhole trace is lost, which may happen for the higher orders
# (especially pre-July 2005 data), more upper y data points may need to be
# deleted.

# IMPORTANT: be sure to apply the same solution for the telluric and the
#            science data. 

imdelete ("ftell_comb.fits", verify=no)
nsfitcoords ("tell_comb.fits", lamptrans="warc_comb", sdisttrans="pinhole", \
    fl_inter=yes, lxorder=2, lyorder=3, sxorder=4, syorder=4)

imdelete ("tftell_comb.fits", verify=no)
nstransform ("ftell_comb.fits")

nxdisplay ("tftell_comb.fits", 1)

imdelete ("fobj_comb.fits", verify=no)
nsfitcoords ("obj_comb.fits", lamptrans="warc_comb", sdisttrans="pinhole", \
    fl_inter=yes, lxorder=2, lyorder=3, sxorder=4, syorder=4)

imdelete ("tfobj_comb.fits", verify=no)
nstransform ("fobj_comb.fits")

nxdisplay ("tfobj_comb.fits", 1)

###############################################################################
# STEP 13: Extract the telluric and science data                              #
###############################################################################

imdelete ("xtftell_comb.fits", verify=no)
nsextract ("tftell_comb.fits", line=750, nsum=20, upper=6, low=-6, \
    fl_inter=yes, fl_apall=yes, fl_trace=yes)

imdelete ("xtfobj_comb.fits", verify=no)
nsextract ("tfobj_comb.fits", line=750, nsum=20, upper=6, low=-6, \
    fl_inter=yes, fl_trace=yes, tr_nsum=5, tr_step=2)

###############################################################################
# FINISHED!                                                                   #
###############################################################################
