# Copyright(c) 2010-2012 Association of Universities for Research in Astronomy, Inc.

###############################################################################
# Gemini F2 example data reduction script                                     #
# Typical reduction for: Longslit Science and Calibration Data                #
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
# Observation UT date : 2012 Jan 04
# Data filename prefix: S20120104S
# File numbers:
#       GCAL flats shutter open   :  50     (J-band, R3K, 4pix-slit, 10s)
#       Darks for the flats       : 381-386 (20s)
#
# Observation UT date : 2012 Jan 05
# Data filename prefix: S20120105S
# File numbers:
#       Arc                       : 187     (J-band, R3K, 4pix-slit, 60s)
#       Darks for the arc         : 246-251 (60s)
#       Object (oj-131)           : 183-186 (J-band, R3K, 4pix-slit, 180s)
#       Darks for the object      : 222-227 (180s)
#
# Observation UT date : 2012 Jan 14
# Data filename prefix: S20120114S
# File numbers:
#       Telluric (HD21999)        :  45-48  (J-band, R3K, 4pix-slit, 5s)
#
# Observation UT date : 2012 Jan 15
# Data filename prefix: S20120115S
# File numbers:
#       Darks for the telluric    :   2-4   (5s)

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
f2.logfile = "f2_longslit_example.log"

# To start from scratch, delete the existing logfile
printf ("EXAMPLE: Deleting %s\n", f2.logfile)
delete (f2.logfile, verify=no)

# Define the database directory
f2.database = "f2_longslit_database/"

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

delete ("flat.lis,flatdark.lis,arc.lis,arcdark.lis,obj.lis,objdark.lis,\
tel.lis,teldark.lis,all.lis", verify=no)

# The user should edit the parameter values in the gemlist calls below to match
# their own dataset.
print ("EXAMPLE: Creating the reduction lists")
gemlist "S20120104S" "50"      > "flat.lis"
gemlist "S20120104S" "381-386" > "flatdark.lis"
gemlist "S20120105S" "187"     > "arc.lis"
gemlist "S20120105S" "246-251" > "arcdark.lis"
gemlist "S20120105S" "183-186" > "obj.lis"
gemlist "S20120105S" "222-227" > "objdark.lis"
gemlist "S20120114S" "45-48"   > "tel.lis"
gemlist "S20120115S" "2-4"     > "teldark.lis"

concat ("flat.lis,flatdark.lis,arc.lis,arcdark.lis,obj.lis,objdark.lis,\
tel.lis,teldark.lis", "all.lis")

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
# STEP 7: Create the normalised flat field and BPM                            #
###############################################################################

# Subtract the dark from the flat images prior to cutting.

imdelete ("df@flat.lis", verify=no)

# If copying and pasting these commands into an interactive PyRAF session, use
# something similar to the following example instead of using the uncommented
# lines below.
#
#     >>> file = open("flat.lis", "r")
#     >>> for line in file:
#     ...    image = line.strip()
#     ...    iraf.gemarith ("f" + image, "-", "flatdark.fits", "df" + image, \
#                fl_vardq=yes, logfile=f2.logfile)
#     ...
#     >>> file.close()

scanfile = "flat.lis"
while (fscan(scanfile, image) != EOF) {
    gemarith ("f"//image, "-", "flatdark.fits", "df"//image, fl_vardq=yes, \
        logfile=f2.logfile)
}
scanfile = ""

imdelete ("cdf@flat.lis", verify=no)
f2cut ("df@flat.lis")

# Construct the normalised flat field. The flats are derived from images taken 
# with the calibration unit (GCAL) shutter open ("lamps-on"). It is recommended
# to run nsflat interactively.

imdelete ("flat.fits,f2_ls_bpm.pl", verify=no)
nsflat ("cdf@flat.lis", flatfile="flat.fits", bpmfile="f2_ls_bpm.pl", \
    thr_flo=0.35, thr_fup=3.0, fl_inter=yes, order=18)

###############################################################################
# STEP 8: Reduce the arc and determine the wavelength solution                #
###############################################################################

# The quality of the fit of the wavelength solution is improved when the arcs 
# are flat fielded. For example, for this dataset, when the arc is flat
# fielded, 3 lines are rejected from the fit and the rms = 0.05576 Angstroms,
# but when the arc is not flat fielded, no lines are rejected, but the rms =
# 0.1855 Angstroms. 

# Subtract the dark from the arc images prior to cutting and flat dividing.

imdelete ("df@arc.lis", verify=no)
nsreduce ("f@arc.lis", outprefix="d", fl_cut=no, fl_process_cut=no, \
    fl_dark=yes, darkimage="arcdark.fits", fl_sky=no, fl_flat=no)

# Cut the arc images and divide by the normalised flat field image.

imdelete ("rdf@arc.lis", verify=no)
nsreduce ("df@arc.lis", fl_cut=yes, fl_dark=no, fl_sky=no, fl_flat=yes, \
    flatimage="flat.fits")

# Combine the arc files (if there is more than one arc file)

imdelete ("arc.fits", verify=no)
delete ("rdfarc.lis", verify=no)
sections "rdf@arc.lis//.fits" > "rdfarc.lis"

# If copying and pasting these commands into an interactive PyRAF session, use
# something similar to the following example instead of using the uncommented
# lines below.
#
#     >>> count = 0
#     >>> file = open("arc.lis", "r")
#     >>> for line in file:
#     ...     count += 1
#     ...
#     >>> if count == 1:
#     ...    iraf.copy ("@rdfarc.lis", "arc.fits")
#     ... else:
#     ...    iraf.gemcombine ("@rdfarc.lis", "arc.fits", fl_vardq=yes)
#     ...
#     >>> file.close()

count ("arc.lis") | scan (num)
if (num == 1) {
    copy ("@rdfarc.lis", "arc.fits")
} else {
    gemcombine ("@rdfarc.lis", "arc.fits", fl_vardq=yes)
}

# Now determine the wavelength solution. It is recommended to run nswavelength
# interactively. The default settings work well for most filter / grism
# combinations. However, for Y band data, the following additional parameters
# should be set: threshold=50, nfound=3, nsum=1.

imdelete ("warc.fits", verify=no)
nswavelength ("arc.fits", fl_inter=yes)

###############################################################################
# STEP 9: Reduce the telluric data                                            #
###############################################################################

# Subtract the dark from the telluric images prior to cutting and flat
# dividing.

imdelete ("df@tel.lis", verify=no)
nsreduce ("f@tel.lis", outprefix="d", fl_cut=no, fl_process_cut=no, \
    fl_dark=yes, darkimage="teldark.fits", fl_sky=no, fl_flat=no)

imdelete ("rdf@tel.lis", verify=no)
nsreduce ("df@tel.lis", fl_cut=yes, fl_dark=no, fl_sky=yes, fl_flat=yes, \
    flatimage="flat.fits")

###############################################################################
# STEP 10: Combine the telluric data                                          #
###############################################################################

imdelete ("tel_comb.fits", verify=no)
nscombine ("rdf@tel.lis", output="tel_comb.fits", fl_shiftint=no, fl_cross=yes)

display ("tel_comb.fits[SCI,1]", 1)

###############################################################################
# STEP 11: Wavelength calibrate the telluric data                             #
###############################################################################

# The nsfitcoords task is used to determine the final solution (consisting of
# the wavelength solution) to be applied to the data. The nstransform task is
# used to apply this final solution. nsfitcoords is best run interactively.

# IMPORTANT: be sure to apply the same solution for the telluric and the
#            science data. 

# Spatial rectification (s-distortion correction) is not usually needed with 
# F2 longslit data. If it is desired, first call nssdist.

imdelete ("ftel_comb.fits", verify=no)
nsfitcoords ("tel_comb.fits", lamptransf="warc.fits")

imdelete ("tftel_comb.fits", verify=no)
nstransform ("ftel_comb.fits")

###############################################################################
# STEP 12: Extract the telluric spectrum                                      #
###############################################################################

imdelete ("xtftel_comb.fits", verify=no)
nsextract ("tftel_comb.fits", fl_apall=yes, fl_findneg=no, fl_inter=no, \
    fl_trace=yes)

splot ("xtftel_comb.fits[SCI,1]")

###############################################################################
# STEP 13: Reduce the science data                                            #
###############################################################################

# Subtract the dark from the science images prior to cutting and flat dividing.

imdelete ("df@obj.lis", verify=no)
nsreduce ("f@obj.lis", outprefix="d", fl_cut=no, fl_process_cut=no, \
    fl_dark=yes, darkimage="objdark.fits", fl_sky=no, fl_flat=no)

imdelete ("rdf@obj.lis", verify=no)
nsreduce ("df@obj.lis", fl_cut=yes, fl_dark=no, fl_sky=yes, fl_flat=yes, \
    flatimage="flat.fits")

###############################################################################
# STEP 14: Combine the science data                                           #
###############################################################################

imdelete ("obj_comb.fits", verify=no)
nscombine ("rdf@obj.lis", output="obj_comb.fits", fl_shiftint=no, \
    fl_cross=yes, rejtype="minmax")

display ("obj_comb.fits[SCI,1]", 1)

###############################################################################
# STEP 15: Wavelength calibrate the science data                              #
###############################################################################

# The nsfitcoords task is used to determine the final solution (consisting of
# the wavelength solution) to be applied to the data. The nstransform task is
# used to apply this final solution. nsfitcoords is best run interactively.

# IMPORTANT: be sure to apply the same solution for the telluric and the
#            science data. 

# Spatial rectification (s-distortion correction) is not usually needed with 
# F2 longslit data. If it is desired, first call nssdist.

imdelete ("fobj_comb.fits", verify=no)
nsfitcoords ("obj_comb.fits", lamptransf="warc.fits")

imdelete ("tfobj_comb.fits", verify=no)
nstransform ("fobj_comb.fits")

###############################################################################
# STEP 16: Extract the science spectrum                                       #
###############################################################################

# For faint science spectra, the telluric can be used as a reference when
# extracting the science spectra; set trace=tftel_comb.fits.

imdelete ("xtfobj_comb.fits", verify=no)
nsextract ("tfobj_comb.fits", fl_apall=yes, fl_findneg=no, fl_inter=no, \
    fl_trace=yes)

splot ("xtfobj_comb.fits[SCI,1]")

###############################################################################
# STEP 17: Apply the telluric correction to the science spectrum              #
###############################################################################

# Note that this telluric has not been corrected to remove intrinsic stellar
# features; this will leave false emission features in the final spectrum.

imdelete ("axtfobj_comb.fits", verify=no)
nstelluric ("xtfobj_comb.fits", "xtftel_comb", fitorder=12, threshold=0.01, \
    fl_inter=yes)

splot ("axtfobj_comb.fits[SCI,1]", ymin=0, ymax=1500)
specplot ("xtfobj_comb.fits[sci,1],axtfobj_comb.fits[sci,1],\
    xtftel_comb.fits[sci,1]", fraction=0.05, yscale=yes, ymin=-100, ymax=4000)

###############################################################################
# STEP 18: Tidy up                                                            #
###############################################################################

delete ("flat.lis,flatdark.lis,fflatdark.lis,arc.lis,arcdark.lis,farcdark.lis,\
rdfarc.lis,obj.lis,objdark.lis,fobjdark.lis,tel.lis,teldark.lis,fteldark.lis,\
all.lis", verify=no)

###############################################################################
# Finished!                                                                   #
###############################################################################
