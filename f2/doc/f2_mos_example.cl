# Copyright(c) 2012 Association of Universities for Research in Astronomy, Inc.

###############################################################################
# Gemini F2 example data reduction script                                     #
# Typical reduction for: MOS Science and Calibration Data                     #
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
# This example shows how difficult it can be to reduce FLAMINGOS-2 MOS data if
# the slits are too close together. To prevent these difficulties, the spectra
# on the array should be separated by at least 2 pixels.
#
# Science Dataset:
#
# Observation UT date : 2012 Jan 04
# Data filename prefix: S20120104S
# File numbers:
#   Darks for the arc            : 394-397 (45s)
#
# Observation UT date : 2012 Jan 05
# Data filename prefix: S20120105S
# File numbers:
#   GCAL flats shutter open      : 143     (Ks-band, R3K, 2pix-slit, 120s)
#   Arc                          : 140     (Ks-band, R3K, 2pix-slit, 45s)
#   Telluric (HIP 31481)         : 136-139 (Ks-band, R3K, 2pix-slit, 120s)
#
# Observation UT date : 2012 Jan 06
# Data filename prefix: S20120106S
# File numbers:
#   GCAL flats shutter open      :  26     (Ks-band, R3K, GS2011SQ600-02, 120s)
#   Darks for the telluric/flats : 202-205 (120s)
#   Arc                          :  23     (Ks-band, R3K, GS2011SQ600-02, 45s)
#   Object (47 Tuc)              :  15-22  (Ks-band, R3K, GS2011SQ600-02, 300s)
#   Darks for the object         : 232-235 (300s)

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
unlearn ("gnirs")
unlearn ("gemtools")

###############################################################################
# STEP 2: Define any variables, the database and the logfile                  #
###############################################################################

# Define any variables (not required if copying and pasting into an interactive
# PyRAF session)
string rawdir, image
int num
struct *scanfile

# Define the logfile
f2.logfile = "f2_mos_example.log"

# To start from scratch, delete the existing logfile
printf ("EXAMPLE: Deleting %s\n", f2.logfile)
delete (f2.logfile, verify=no)

# Define the database directory
f2.database = "f2_mos_database/"

# To start from scratch, delete the existing database files

# If copying and pasting these commands into an interactive PyRAF session, use
# something similar to the following example instead of using the uncommented
# lines below.
#
#     >>> if (iraf.access(f2.database)):
#     ...    print "EXAMPLE: Deleting contents of %s" % (f2.database)
#     ...    iraf.delete (f2.database + "*", verify=no)

if (access(f2.database)) {
    printf ("EXAMPLE: Deleting contents of %s\n", f2.database)
    delete (f2.database//"*", verify=no)
}
;

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

delete ("flatobj.lis,flattel.lis,flatdark.lis,arcobj.lis,arctel.lis,\
arcdark.lis,obj.lis,objdark.lis,tel.lis,teldark.lis,all.lis", verify=no)

# The user should edit the parameter values in the gemlist calls below to match
# their own dataset.
print ("EXAMPLE: Creating the reduction lists")
gemlist "S20120106S" "26"      > "flatobj.lis"
gemlist "S20120105S" "143"     > "flattel.lis"
gemlist "S20120106S" "202-205" > "flatdark.lis"
gemlist "S20120106S" "23"      > "arcobj.lis"
gemlist "S20120105S" "140"     > "arctel.lis"
gemlist "S20120104S" "394-397" > "arcdark.lis"
gemlist "S20120106S" "15-22"   > "obj.lis"
gemlist "S20120106S" "232-235" > "objdark.lis"
gemlist "S20120105S" "136-139" > "tel.lis"
gemlist "S20120104S" "400-403" > "teldark.lis"

concat ("flatobj.lis,flattel.lis,flatdark.lis,arcobj.lis,arctel.lis,\
arcdark.lis,obj.lis,objdark.lis,tel.lis,teldark.lis", "all.lis")

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
# STEP 5: f2prepare all the data                                              #
###############################################################################

# Run F2PREPARE on all the data to update the headers, derive variance and data
# quality (DQ) planes, correct for non-linearity (not yet implemented) and flag
# saturated and non-linear pixels in the DQ plane.

imdelete ("f@all.lis", verify=no)
f2prepare ("@all.lis", rawpath=rawdir, fl_vardq=yes, fl_correct=yes, \
    fl_saturated=yes, fl_nonlinear=yes)

###############################################################################
# STEP 6: Create the necessary dark images                                    #
###############################################################################

delete ("fflatdark.lis", verify=no)
imdelete ("flatdark.fits", verify=no)
sections "f@flatdark.lis" > "fflatdark.lis"
gemcombine ("@fflatdark.lis", "flatdark.fits", combine="average", \
    fl_vardq=yes, logfile=f2.logfile)

delete ("farcdark.lis", verify=no)
imdelete ("arcdark.fits", verify=no)
sections "f@arcdark.lis" > "farcdark.lis"
gemcombine ("@farcdark.lis", "arcdark.fits", combine="average", fl_vardq=yes, \
    logfile=f2.logfile)

delete ("fobjdark.lis", verify=no)
imdelete ("objdark.fits", verify=no)
sections "f@objdark.lis" > "fobjdark.lis"
gemcombine ("@fobjdark.lis", "objdark.fits", combine="average", fl_vardq=yes, \
    logfile=f2.logfile)

delete ("fteldark.lis", verify=no)
imdelete ("teldark.fits", verify=no)
sections "f@teldark.lis" > "fteldark.lis"
gemcombine ("@fteldark.lis", "teldark.fits", combine="average", fl_vardq=yes, \
    logfile=f2.logfile)

###############################################################################
# STEP 7: Trace the slit edges using the flat field                           #
###############################################################################

# Subtract the dark from the flat images prior to cutting.

imdelete ("df@flatobj.lis", verify=no)

# If copying and pasting these commands into an interactive PyRAF session, use
# something similar to the following example instead of using the uncommented
# lines below.
#
#     >>> file = open("flatobj.lis", "r")
#     >>> for line in file:
#     ...    image = line.strip()
#     ...    iraf.gemarith ("f" + image, "-", "flatdark.fits", "df" + image, \
#                fl_vardq=yes, logfile=f2.logfile)
#     ...
#     >>> file.close()

scanfile = "flatobj.lis"
while (fscan(scanfile, image) != EOF) {
    gemarith ("f"//image, "-", "flatdark.fits", "df"//image, fl_vardq=yes, \
        logfile=f2.logfile)
}
scanfile = ""

# Combine the flat files (if there is more than one flat file). This combined
# flat will be used to trace the edges of the slits in f2cut.

imdelete ("combflat.fits", verify=no)
delete ("dfflatobj.lis", verify=no)
sections "df@flatobj.lis//.fits" > "dfflatobj.lis"

# If copying and pasting these commands into an interactive PyRAF session, use
# something similar to the following example instead of using the uncommented
# lines below.
#
#     >>> count = 0
#     >>> file = open("flatobj.lis", "r")
#     >>> for line in file:
#     ...     count += 1
#     ...
#     >>> if count == 1:
#     ...    iraf.copy ("@dfflatobj.lis", "combflat.fits")
#     ... else:
#     ...    iraf.gemcombine ("@dfflatobj.lis", "combflat.fits", fl_vardq=yes)
#     ...
#     >>> file.close()

count ("flatobj.lis") | scan (num)
if (num == 1) {
    copy ("@dfflatobj.lis", "combflat.fits")
} else {
    gemcombine ("@dfflatobj.lis", "combflat.fits", fl_vardq=yes)
}

# First cut the combined flat so that it can be used as a reference image.
# Since some of the slits in the flat are very close together, 2 of the slits
# are completed ignored, since no slit edges are found in the gradient image
# that can be associated to those slits in the MDF. In addition, only one slit
# edge is associated to 9 of the slits in the MDF (EXTVER = 14, 16, 17, 19, 28,
# 29, 30, 34 and 35). These 9 slits will have poor / no s-distortion
# corrections and so the final spectra will most likely be unusable.

imdelete ("cutcombflat.fits,slits.fits", verify=no)
f2cut ("combflat.fits", outimages="cutcombflat.fits", \
    gradimage="combflat.fits", edgediff=3.5)

imdelete ("cdf@flatobj.lis", verify=no)
f2cut ("df@flatobj.lis", refimage="cutcombflat.fits")

f2display ("cutcombflat.fits", 1)

###############################################################################
# STEP 8: Create the normalised flat field and BPM for the telluric data      #
###############################################################################

# Subtract the dark from the flat images prior to cutting.

imdelete ("df@flattel.lis", verify=no)

# If copying and pasting these commands into an interactive PyRAF session, use
# something similar to the following example instead of using the uncommented
# lines below.
#
#     >>> file = open("flattel.lis", "r")
#     >>> for line in file:
#     ...    image = line.strip()
#     ...    iraf.gemarith ("f" + image, "-", "flatdark.fits", "df" + image, \
#                fl_vardq=yes, logfile=f2.logfile)
#     ...
#     >>> file.close()

scanfile = "flattel.lis"
while (fscan(scanfile, image) != EOF) {
    gemarith ("f"//image, "-", "flatdark.fits", "df"//image, fl_vardq=yes, \
        logfile=f2.logfile)
}
scanfile = ""

imdelete ("cdf@flattel.lis", verify=no)
f2cut ("df@flattel.lis")

# Construct the normalised flat field. The flats are derived from images taken 
# with the calibration unit (GCAL) shutter open ("lamps-on"). It is recommended
# to run nsflat interactively.

imdelete ("flattel.fits,f2_ls_bpm.pl", verify=no)
nsflat ("cdf@flattel.lis", flatfile="flattel.fits", bpmfile="f2_ls_bpm.pl", \
    thr_flo=0.35, thr_fup=3.0, fl_inter=yes, order=18)

###############################################################################
# STEP 9: Create the normalised flat field and BPM for the science data       #
###############################################################################

# Construct the normalised flat field. The flats are derived from images taken 
# with the calibration unit (GCAL) shutter open ("lamps-on"). It is recommended
# to run nsflat interactively. The response should look similar for all slits
# (and similar to the telluric flat, since they are taken with the same filter
# / grism combination). However, for this data set, the response (and fit)
# looks particularly bad for EXTVER = 28, 29 and 31. However, this is to be
# expected, since f2cut had difficulty cutting these slits. 

imdelete ("flatobj.fits,f2_mos_bpm.pl", verify=no)
nsflat ("cdf@flatobj.lis", flatfile="flatobj.fits", bpmfile="f2_mos_bpm.pl", \
    thr_flo=0.35, thr_fup=3.0, fl_inter=yes, order=18)

f2display ("flatobj.fits", 1)

###############################################################################
# STEP 10: Reduce the arc for the telluric data and determine the wavelength  #
#          solution                                                           #
###############################################################################

# The quality of the fit of the wavelength solution is often improved when the
# arcs are flat fielded. However, for this dataset, flat fielding is not
# required; when the arc is flat fielded, no lines are rejected from the fit
# and the rms = 0.2784 Angstroms, but when the arc is not flat fielded, no
# lines are rejected and the rms = 0.2728 Angstroms.

# Subtract the dark from the arc images prior to cutting and flat dividing.

imdelete ("df@arctel.lis", verify=no)
nsreduce ("f@arctel.lis", outprefix="d", fl_cut=no, fl_process_cut=no, \
    fl_dark=yes, darkimage="arcdark.fits", fl_sky=no, fl_flat=no)

# Cut the arc images and divide by the normalised flat field image.

imdelete ("rdf@arctel.lis", verify=no)
nsreduce ("df@arctel.lis", fl_cut=yes, fl_dark=no, fl_sky=no, fl_flat=no)

# Combine the arc files (if there is more than one arc file)

imdelete ("arctel.fits", verify=no)
delete ("rdfarctel.lis", verify=no)
sections "rdf@arctel.lis//.fits" > "rdfarctel.lis"

# If copying and pasting these commands into an interactive PyRAF session, use
# something similar to the following example instead of using the uncommented
# lines below.
#
#     >>> count = 0
#     >>> file = open("arctel.lis", "r")
#     >>> for line in file:
#     ...     count += 1
#     ...
#     >>> if count == 1:
#     ...    iraf.copy ("@rdfarctel.lis", "arctel.fits")
#     ... else:
#     ...    iraf.gemcombine ("@rdfarctel.lis", "arctel.fits", fl_vardq=yes)
#     ...
#     >>> file.close()

count ("arctel.lis") | scan (num)
if (num == 1) {
    copy ("@rdfarctel.lis", "arctel.fits")
} else {
    gemcombine ("@rdfarctel.lis", "arctel.fits", fl_vardq=yes)
}

# Now determine the wavelength solution. It is recommended to run nswavelength
# interactively. The default settings work well for most filter / grism
# combinations. However, for Y band data, the following additional parameters
# should be set: threshold=50, nfound=3, nsum=1.

imdelete ("warctel.fits", verify=no)
nswavelength ("arctel.fits", fl_inter=yes)

###############################################################################
# STEP 11: Reduce the arc for the science data and determine the wavelength   #
#          solution                                                           #
###############################################################################

# The quality of the fit of the wavelength solution is improved when the arcs 
# are flat fielded. 

# Subtract the dark from the arc images prior to cutting and flat dividing.

imdelete ("df@arcobj.lis", verify=no)
nsreduce ("f@arcobj.lis", outprefix="d", fl_cut=no, fl_process_cut=no, \
    fl_dark=yes, darkimage="arcdark.fits", fl_sky=no, fl_flat=no)

# Cut the arc images and divide by the normalised flat field image.

imdelete ("rdf@arcobj.lis", verify=no)
nsreduce ("df@arcobj.lis", fl_cut=yes, refimage="cutcombflat.fits", \
    fl_dark=no, fl_sky=no, fl_flat=yes, flatimage="flatobj.fits")

# Combine the arc files (if there is more than one arc file)

imdelete ("arcobj.fits", verify=no)
delete ("rdfarcobj.lis", verify=no)
sections "rdf@arcobj.lis//.fits" > "rdfarcobj.lis"

# If copying and pasting these commands into an interactive PyRAF session, use
# something similar to the following example instead of using the uncommented
# lines below.
#
#     >>> count = 0
#     >>> file = open("arcobj.lis", "r")
#     >>> for line in file:
#     ...     count += 1
#     ...
#     >>> if count == 1:
#     ...    iraf.copy ("@rdfarcobj.lis", "arcobj.fits")
#     ... else:
#     ...    iraf.gemcombine ("@rdfarcobj.lis", "arcobj.fits", fl_vardq=yes)
#     ...
#     >>> file.close()

count ("arcobj.lis") | scan (num)
if (num == 1) {
    copy ("@rdfarcobj.lis", "arcobj.fits")
} else {
    gemcombine ("@rdfarcobj.lis", "arcobj.fits", fl_vardq=yes)
}

# Determine the s-distortion correction. When running nssdist interactively, a
# window will pop up showing two peaks that correspond to slit edges found in
# the gradient image by f2cut. The expected positions of the peaks (slit edges)
# are printed to screen (under "Coordinate list:"). Press m on one of the two
# visible peaks. The first two of the three numbers that appear at the bottom
# of the plot window indicate the pixel position of peak that was just marked.
# The third of the three numbers (the one in parentheses) will match with one
# of the expected positions of the slit edges. Check that these numbers match
# well (they should be within ~1 pixel) and press return. Repeat for the second
# peak. Press q to move to the next extension.

# For this data set, for extensions with EXTVER = 14, 16, 17, 19, 28, 29, 30,
# 34 and 35, only one peak should be visible in the plot. The f2cut task prints
# to screen whether the left slit edge or the right slit edge was used.

# For data with well separated slits, it is not necessary to run nssdist
# interactively and the following call works well
#nssdist ("slits.fits", nsum=5, fwidth=4.0, cradius=5.0, minsep=2.0)

# However, since f2cut had difficulty detecting slit edges in this troublesome
# data set, the following call works better. Running interactively and pressing
# e on the plot to automatically mark the peaks, only one peak is marked for
# slits with EXTVER = 9, 11, 14 (use m to mark right peak), 16 (right peak
# correctly marked), 17 (delete left peak), 19 (delete left peak), 20 is
# skipped, 22 is incorrectly marked, 24, 25, 28, 29, 30, 31, 32 (use m to mark
# first right peak), 33 is skipped, 34, 35, 36 is skipped.

nssdist ("slits.fits", fl_inter=yes, nsum=2, fwidth=4.0, cradius=5.0, \
    minsep=2.0)

# Apply the s-distortion correction to the arc

imdelete ("farcobj.fits", verify=no)
nsfitcoords ("arcobj.fits", sdisttransf="slits.fits")

imdelete ("tfarcobj.fits", verify=no)
nstransform ("farcobj.fits")

f2display ("tfarcobj.fits", 1)

# Now determine the wavelength solution. It is recommended to run nswavelength
# interactively. The default settings work well for most filter / grism
# combinations. However, for Y band data, the following additional parameters
# should be set: threshold=50, nfound=3, nsum=1.

# For this data set, the lines were not automatically identified for only two
# extensions (EXTVER = 22, 28), which is expected since the s-distortion
# solutions determined for these extensions are poor.

imdelete ("wtfarcobj.fits", verify=no)
nswavelength ("tfarcobj.fits", sdist="slits.fits", fl_inter=yes, step=2, \
    threshold=350)

###############################################################################
# STEP 12: Reduce the telluric data                                           #
###############################################################################

# Subtract the dark from the telluric images prior to cutting and flat
# dividing.

imdelete ("df@tel.lis", verify=no)
nsreduce ("f@tel.lis", outprefix="d", fl_cut=no, fl_process_cut=no, \
    fl_dark=yes, darkimage="teldark.fits", fl_sky=no, fl_flat=no)

imdelete ("rdf@tel.lis", verify=no)
nsreduce ("df@tel.lis", fl_cut=yes, fl_dark=no, fl_sky=yes, fl_flat=yes, \
    flatimage="flattel.fits")

###############################################################################
# STEP 13: Combine the telluric data                                          #
###############################################################################

imdelete ("tel_comb.fits", verify=no)
nscombine ("rdf@tel.lis", output="tel_comb.fits", fl_shiftint=no, fl_cross=yes)

display ("tel_comb.fits[SCI,1]", 1)

###############################################################################
# STEP 14: Wavelength calibrate the telluric data                             #
###############################################################################

# The nsfitcoords task is used to determine the final solution (consisting of
# the wavelength solution) to be applied to the data. The nstransform task is
# used to apply this final solution. nsfitcoords is best run interactively.

# IMPORTANT: be sure to apply the same solution for the telluric and the
#            science data. 

imdelete ("ftel_comb.fits", verify=no)
nsfitcoords ("tel_comb.fits", lamptransf="warctel.fits")

imdelete ("tftel_comb.fits", verify=no)
nstransform ("ftel_comb.fits")

display ("tftel_comb.fits[SCI,1]", 2)

###############################################################################
# STEP 15: Extract the telluric spectrum                                      #
###############################################################################

imdelete ("xtftel_comb.fits", verify=no)
nsextract ("tftel_comb.fits", fl_apall=yes, fl_findneg=no, fl_inter=no, \
    fl_trace=yes)

splot ("xtftel_comb.fits[SCI,1]")

###############################################################################
# STEP 16: Reduce the science data                                            #
###############################################################################

# Subtract the dark from the science images prior to cutting and flat dividing.

imdelete ("df@obj.lis", verify=no)
nsreduce ("f@obj.lis", outprefix="d", fl_cut=no, fl_process_cut=no, \
    fl_dark=yes, darkimage="objdark.fits", fl_sky=no, fl_flat=no)

imdelete ("rdf@obj.lis", verify=no)
nsreduce ("df@obj.lis", fl_cut=yes, refimage="cutcombflat.fits", fl_dark=no, \
    fl_sky=yes, fl_flat=yes, flatimage="flatobj.fits")

###############################################################################
# STEP 17: Combine the science data                                           #
###############################################################################

imdelete ("obj_comb.fits", verify=no)
nscombine ("rdf@obj.lis", output="obj_comb.fits", fl_shiftint=no, \
    fl_cross=yes, rejtype="minmax")

f2display ("obj_comb.fits", 1)

###############################################################################
# STEP 18: Wavelength calibrate the science data                              #
###############################################################################

# The nsfitcoords task is used to determine the final solution (consisting of
# the wavelength solution and the spatial distortion solution) to be applied to
# the data. The nstransform task is used to apply this final solution.
# nsfitcoords is best run interactively.

# IMPORTANT: be sure to apply the same solution for the telluric and the
#            science data. 

imdelete ("fobj_comb.fits", verify=no)
nsfitcoords ("obj_comb.fits", lamptransf="wtfarcobj.fits", \
    sdisttransf="slits.fits")

imdelete ("tfobj_comb.fits", verify=no)
nstransform ("fobj_comb.fits")

f2display ("tfobj_comb.fits", 2)

###############################################################################
# STEP 19: Extract the science spectrum                                       #
###############################################################################

# When extracting MOS data, nsextract is best run interactively.
imdelete ("xtfobj_comb.fits", verify=no)
nsextract ("tfobj_comb.fits", nsum=2, fl_apall=yes, fl_findneg=no, \
    fl_inter=yes, fl_trace=yes, tr_nsum=2, tr_step=2)

splot ("xtfobj_comb.fits[SCI,1]")

###############################################################################
# STEP 20: Apply the telluric correction to the science spectrum              #
###############################################################################

# Note that this telluric has not been corrected to remove intrinsic stellar
# features; this will leave false emission features in the final spectrum.

imdelete ("axtfobj_comb.fits", verify=no)
nstelluric ("xtfobj_comb.fits", "xtftel_comb", fitorder=15, threshold=0.01, \
    fl_inter=yes)

splot ("axtfobj_comb.fits[SCI,1]", ymin=-50, ymax=100)
specplot ("xtfobj_comb.fits[sci,1],axtfobj_comb.fits[sci,1],\
    xtftel_comb.fits[sci,1]", fraction=0.015, yscale=yes, ymin=-100, ymax=600)

###############################################################################
# STEP 21: Tidy up                                                            #
###############################################################################

delete ("flatobj.lis,flattel.lis,flatdark.lis,fflatdark.lis,dfflatobj.lis,\
arcobj.lis,arctel.lis,arcdark.lis,farcdark.lis,rdfarcobj.lis,rdfarctel.lis,\
obj.lis,objdark.lis,fobjdark.lis,tel.lis,teldark.lis,fteldark.lis,all.lis", \
    verify=no)

###############################################################################
# Finished!                                                                   #
###############################################################################
