# Copyright(c) 2003-2006 Association of Universities for Research in Astronomy, Inc.
#
# T-ReCS example reduction script: Typical reduction of standard star observations
#                                 in imaging mode
#
# This data processing was done to aid the data quality assessment performed
# by the Gemini staff. The data processing is not designed to give the best
# possible result. The user is encouraged to carry out a more careful reduction
# to improve the resulting signal-to-noise ratio of the final images
#
# Gemini T-ReCS data reduction script
# Observation UT date: 2003sep16
#
# Brief data description:  Engineering data. Standard star observations from
#                          GS-2003B-SV-101-eng
#
# In this script it is assumed that the raw data files are in directory
#
# 	/net/tyl/export/data/gemini_testdata/trecs_testdata
#
# and that the script is run in a different directory (to avoid changing the 
# raw data files).  Change the directory name to whatever value is appropriate
# for your setup.
#
#
# One needs to load the gemini and midir packages before carrying 
# out these processing steps.  In IRAF it looks something like this (enter 
# "gemini", then "midir" at the IRAF prompt):
#
#  cl> gemini
#        flamingos.   gmos.        niri.        quirc.       
#        gemtools.    midir.       oscir.       
#  ge> midir
#        midirinfo          miregister         mistack            tbackground
#        mbackground        midobpm            miview             tprepare    
#        mibackground       midoflat           mbackground        tview
#        mibpm              miflat             mprepare                   
#        midirexamples      mireduce           mview                    
#  mi> 
#
#

# Set the logfile name
midir.logfile="20030916_GS-2003B-SV-101-eng.log"

# The following does the step by step reduction for one file out of the 
# test data file set.
#
# Step 1: check for bad frames automaticly (it will not find any bad frames
#         in this data file)
#
# This produces file bS20030916S0067.fits.
#

tbackground S20030916S0067 rawpath=/net/tyl/export/data/gemini_testdata/trecs_testdata

#
# Step 2: look at the frames interactively to detect bad frames
#         (again, not needed in this case since all the frames are good)
#
#         Type "x" in the display window to exit.
#
# The "fl_delete-" flag means that there will be an output file whether
# or not any frames are marked as "bad" (so that the subsequent steps 
# can use this file).
#
# It produces file vbS20030916S0067.fits.

tview bS20030916S0067 zrange=yes zscale=no ztrans="linear" fl_delete-

#
# Step 3: "prepare" the file (stack up frames for each NOD position)
# 
# This produces file tvbS20030916S0067.fits

tprepare vbS20030916S0067

#
# Step 4: stack the frames for the different NODs
#
# This produces file stvbS20030916S0067.fits

mistack tvbS20030916S0067

#
# Step 5: mask bad pixels
#         (using the package bad pixel mask in this case)
#
# This produces file fstvbS20030916S0067.fits.

mireduce stvbS20030916S0067 outpref="f" fl_mask+ fl_flat- fl_background- fl_view- stackoption="stack" fl_display=no

#
# for convenience: rename the file
#

imrename fstvbS20030916S0067 frame67

#
# look at the resulting image in the stdimage display, linear scale
#

display frame67[1] 1 zrange=yes zscale=no ztrans="linear"

#
# Optional step: make a flat field frame (uses two of the test data files)
#
# Actually the flat field file that results is not useful since the 
# attempted sky flat observations used here are not of good quality.
#
# The output file is named flat10p4.fits.

miflat rS20030714S0141 rS20030714S0142 flat10p4 rawpath=/net/tyl/export/data/gemini_testdata/trecs_testdata

#
# Optional step: apply the flat field using "mireduce"
#
# Since the "flat10p4" file is not of good quality, this does not
# actually improve the resulting image.
#
# The output file name is specified as flat_frame67.fits by the command

mireduce frame67 outimage="flat_frame67" flatfieldfile="flat10p4" fl_flat+ fl_background- fl_view- fl_mask- stackoption="stack" fl_display=no logfile="flat.log"

#
# Minimum standard procssing of all the files: it "prepares" and then "stacks"
# all the raw data files in the test data directory.
#
# The output files are named "reducedS20030916S0067.fits", etc.

mireduce S20030916* rawpath=/net/tyl/export/data/gemini_testdata/trecs_testdata outpref="reduced" fl_background+ fl_view- fl_mask- fl_flat- stackoption="stack" fl_display=no logfile="process.log"

#
# If one wishes to register the frames rather than stacking them, one can set 
# the "stackoption" to "register"
#
# This produces files ttS20030916S0074.fits and regS20030916S0074.fits, the 
# first one being the registered "prepared" image, the second one being 
# the final registered image.

mireduce S20030916S0074 rawpath=/net/tyl/export/data/gemini_testdata/trecs_testdata outpref="reg" fl_background- fl_view- fl_mask- fl_flat- stackoption="register" fl_display=no logfile="process.log"

# If it all ran successfully, you should have the following files in the 
# directory:
#
# bS20030916S0067.fits
# files.list
# flat10p4.fits
# flat_frame67.fits
# flat.log
# frame67.fits
# midir.log
# process.log
# reducedS20030916S0067.fits
# reducedS20030916S0068.fits
# reducedS20030916S0069.fits
# reducedS20030916S0070.fits
# reducedS20030916S0071.fits
# reducedS20030916S0072.fits
# reducedS20030916S0073.fits
# reducedS20030916S0074.fits
# regS20030916S0074.fits
# stvbS20030916S0067.fits
# tS20030916S0067.fits
# tS20030916S0068.fits
# tS20030916S0069.fits
# tS20030916S0070.fits
# tS20030916S0071.fits
# tS20030916S0072.fits
# tS20030916S0073.fits
# tS20030916S0074.fits
# ttS20030916S0074.fits
# tvbS20030916S0067.fits
# vbS20030916S0067.fits
#
