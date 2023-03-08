# Copyright(c) 2000-2017 Association of Universities for Research in Astronomy, Inc.
#
# Package script for the niri package
#
# Version: Sept 14, 2002 BR,IJ  Release v1.4
#          Jan   9, 2004 KL     Release v1.5
#          Apr  19, 2004 KL     Release v1.6
#          Oct  25, 2004 KL     Release v1.7
#          May   6, 2005 KL     Release v1.8
#          Jul  28, 2006 KL     Release v1.9
#          Jul  28, 2009 JH     Release v1.10
#          Jan  13, 2011 EH,KL  Beta release v1.11beta
#          Dec  30, 2011 EH     Release v1.11
#          Mar  28, 2012 EH     Release v1.11.1
#          Dec  13, 2012 EH     Beta release v1.12beta
#          May  14, 2013 EH     Beta release v1.12beta2
#          Oct  11, 2013 EH     Release v1.12
#          Jan  30, 2015 KL     Release v1.13
#	   Dec   7, 2015 KL	Release v1.13.1
#          Jul  20, 2017 KL     Release v1.14
#
# Modifications : Package niri created
#          Oct 29, 2003 KL  ns* tasks moved to gnirs package
#          Sep 16, 2004 JJ  added nresidual, changed nisky,niflat
#          Feb. 4, 2005 JJ  added nirotate
#          May 27, 2005 JJ  added nisupersky, changed nisky
#
# Load necessary packages - if anything else than loaded
# by gemini should be necessary
gemtools
gnirs

# spectroscopy packages for longslit - may be moved to gemini later
twodspec
longslit

print ""
print "Loading the gnirs package:"
? gnirs
print ""
print "Loading the niri package:"

package niri

# Generic preparations
task      nprepare=niri$nprepare.cl
task      nresidual=niri$nresidual.cl

# Image reductions
task      nifastsky=niri$nifastsky.cl
task      niflat=niri$niflat.cl
task      nisky=niri$nisky.cl
task      nireduce=niri$nireduce.cl
task      nirotate=niri$nirotate.cl
task      nisupersky=niri$nisupersky.cl

# Spectroscopy reductions through GNIRS package

# Cook book task
task      niriinfo=niri$niriinfo.cl
task      niriexamples=niri$niriexamples.cl

clbye()
