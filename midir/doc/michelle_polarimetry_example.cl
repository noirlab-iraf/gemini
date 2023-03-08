# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.
#
# Polarimetry example:
#
# This carries out a standard reduction of a polarimetry observation.
# Two approaches are presented.
#
# For more details, see the "Data Reduction" section of the Michelle webpage:
#   http://www.gemini.edu/sciops/instruments/michelle/MichIndex.html



# Set directory to raw data
set rawdir = "/data/michpol/raw/"

# Reset configuration
unlearn midir

midir.logfile = "michpol.log"

# Preparing the data.
#
# The only unusual thing is that the initial step needs to have the
# fl_rescue flat active since the observation was stopped part way
# through.  This would not normally be needed.

mprepare N20060105S0211 rawpath="rawdir$" fl_rescue+

#----------------------------------
# First Approach

# Stack each of the four positions of the waveplate in Michelle
# A registration source has been identified in region [92:122,132:152]

mipstack mN20060105S0211 frametype="dif" combine="average" fl_register+ \
    regions="[92:122,132:152]" fl_stair+

# Calculate the Stokes parameters from the stacked Michelle polarimetry file
# A registration source has been identified in region [92:122,132:152]

miptrans smN20060105S0211 fl_register+ regions="[92:122,132:152]"

# Output file: psmN20060105S0211.fits
#----------------------------------

#----------------------------------
# Second Approach

# Calculate the Stokes parameters for each AB NOD.
# A registration source has been identified in region [92:122,132:152]

mipstokes mN20060105S0211 frametype="dif" combine="average" fl_register+ \
    regions="[92:122,132:152]" fl_stair+ fl_mask-

# Stack each set of Stokes parameters into the final Stokes images.
# A registration source has been identified in region [92:122,132:152]

mipsstk zmN20060105S0211 fl_register+ regions="[92:122,132:152]" \
    fl_variance+ fl_stair+

# Output file: azmN20060105S0211.fits
#----------------------------------
