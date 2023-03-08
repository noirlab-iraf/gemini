# Copyright(c) 2010-2017 Association of Universities for Research in Astronomy, Inc.
#
# Package script for the f2 package
#
# Version: Oct  11, 2013 EH     Release v1.12
#          Jan  15, 2015 KL     Release v1.13
#	   Dec   7, 2015 KL	Release v1.13.1
#          Jul  20, 2017 KL     Release v1.14
#
# Load necessary packages. The niri package loads the gemtools and gnirs
# packages.
niri

? niri
print ""
print "Loading the f2 package:"

package f2

# Generic preparations
task    f2prepare=f2$f2prepare.cl

# Spectroscopic tasks
task    f2cut=f2$f2cut.cl
task    f2display=f2$f2display.cl

# Cookbooks
task    f2info=f2$f2info.cl
task    f2infoimaging=f2$f2infoimaging.cl
task    f2infols=f2$f2infols.cl
task    f2infomos=f2$f2infomos.cl

# Examples
task    f2examples=f2$f2examples.cl

clbye()
