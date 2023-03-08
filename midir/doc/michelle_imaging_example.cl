# Copyright(c) 2003-2006 Association of Universities for Research in Astronomy, Inc.
#
# Michelle example reduction script: Typical reduction of image data
#
# This data processing was done to aid the data quality assessment performed
# by the Gemini staff. The data processing is not designed to give the best
# possible result. Better signal-to-noise and better cleaning for cosmic-ray
# hits and bad pixels will most likely be possible. The user of these data
# is encouraged to use the provided co-added images only as guide lines and
# to re-reduce the data to obtain the best possible reduction.
#
# Gemini Michelle data reduction script
# Observation UT date: 2003nov09 
#
# December, 10, 2003  TLB
#
# Brief data description:  Engineering: Silicate filter set of Vega images
#                          GN-CAL20031109
#
#	Frame 0058 = 7.9 micron image
#	Frame 0059 = 8.8 micron image
#	Frame 0060 = 9.7 micron image
#	Frame 0061 = 10.3 micron image
#	Frame 0062 = 11.6 micron image
#	Frame 0063 = 12.5 micron image
#

# Define directory containing raw data and calibration data
string name,rawdir
name = "N20031109S"
rawdir = "/net/tyl/export/data/gemini_testdata/michelle_testdata"

# Set the logfile name
midir.logfile="GN-CAL20031109.log"


delete inlist ver-
for(i=58; i<=63; i+=1) {
  print(name//"0000"+i, >> "inlist")
}

# To view the Michelle data prior to its preparation for further processing, 
# use the task MVIEW.  This script displays each extension of a Michelle image, 
# if the "fl_inter" parameter is set it will display an extension and wait for 
# interactive input from the user.  After the data is processed using MPREPARE, 
# the task MIVIEW should be used to view the extensions.

mview("@inlist", rawpath=rawdir)

# Prepare the data for further reduction
mprepare("@inlist", rawpath=rawdir)

# View the images to test their quality (use fl_inter+ for interactive
# mode to identify bad frames).

sections m@inlist > minlist
miview("@minlist")

# Stack the nod sets for each individual image.  If MIVIEW was run
# interactively, a "v" or appropriate file prefix must be added to
# the input list for MISTACK.

mistack("@minlist")

# Alternatively, the data can be coadded using registration before
# the data are stacked:

miregister("@minlist")

# Main imaging reduction is complete - bad pixel mask and flat fielding
# is not presently available for Michelle data.

# To do all of the above reduction in one single step, you can use
# the below call to MIREDUCE which will prepare and stack the frames,
# and provides the option to view them if the keyword fl_view is
# set.

mireduce("@inlist", rawpath=rawdir)
