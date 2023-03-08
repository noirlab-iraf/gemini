# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.
#
# Package script for the flamingos package
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
# Load necessary packages - if anything else than loaded
# by gemini should be necessary
gemtools
niri

? niri
print ""
print "Loading the flamingos package:"

package flamingos

task      fprepare=flamingos$fprepare.cl
task      flamingosinfo=flamingos$flamingosinfo.cl

clbye()
