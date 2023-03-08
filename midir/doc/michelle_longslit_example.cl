# Copyright(c) 2006-2009 Association of Universities for Research in Astronomy, Inc.

# Reducing a spectrum with MSREDUCE
#
# Note that quite often, the line identification must be done interactively.
# The same applies to the telluric correction.
#
# For more details, see the "Data Reduction" section of the MIchelle webpage:
#   http://www.gemini.edu/sciops/instruments/michelle/MichIndex.html


# Set directory to raw data
set rawdir = "/data/michdata/raw/"

# Reset configuration
unlearn midir
unlearn gnirs

midir.logfile = "spectro.log"

# Run NSHEADERS to configure the GNIRS package for, in this case, Michelle

nsheaders michelle

# Run MSREDUCE
#   MSREDUCE will call MPREPARE, and whichever tasks necessary to do what
#   the several options say it should be doing.  The flat field and the
#   bias frames are prepared and reduced by MSREDUCE.
#
#   When MSABSFLUX calls TELLURIC, one should reset the shift to 0.0, then
#   iterate:
#
#       :shift 0.0
#
#   The 'scale' often needs to be modified too (e.g. :scale 1.).  
#   The goal is to get as smooth a spectrum as possible around 
#   9.5 microns where the ozone feature is seen in the original spectrum.
#
#   In this example, de-fringing is turn off simply because it was not
#   needed for that spectrum.  The default behaviour is to de-fringe.
#
#   WARNING: At this time, MSREDUCE does not allow the user to set output
#            prefixes, or output file names.  Please refer to the MSREDUCE
#            help page.

msreduce N20051218S0083 outtype="fnu" rawpath="rawdir$" fl_std+ std="HD57423" \
    fl_flat+ flat=N20051218S0062 bias=N20051218S0063 fl_telluric+ \
    fl_wavelength+ fl_extract+ fl_defringe-


print ("Please see the Gemini midir webpage for a full example reduction examples")

print ("http://www.gemini.edu/sciops/instruments/midir-resources/data-reduction/spectroscopy-reduction")

