# Copyright(c) 2006-2017 Association of Universities for Research in Astronomy, Inc.
#
# Package script for the nifs package
#
# Version: Mar  24, 2006 TB,KL,IS,CA  Beta release v1.9beta
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
# Load necessary pacakges - if anything else than loaded by gemini
# should be necessary

gemtools
gnirs

# spectroscopy packages twodspec & longslit loaded by gnirs

print ""
print "Loading the gnirs package:"
? gnirs
print ""
print "Loading the nifs package:"

package nifs

# Task definitions
task    nifcube=nifs$nifcube.cl
task    nfwcs=nifs$nfwcs.cl
task    nfprepare=nifs$nfprepare.cl
task    nfimage=nifs$nfimage.cl
task    nfacquire=nifs$nfacquire.cl
task    nfdispc=nifs$nfdispc.cl
task    nfmap=nifs$nfmap.cl
task    nfpad=nifs$nfpad.cl
task    nffixbad=nifs$nffixbad.cl
task    nfsdist=nifs$nfsdist.cl
task    nfextract=nifs$nfextract.cl
task    nftelluric=nifs$nftelluric.cl
hidetask nfwcs nfpad

# Cook book tasks
task    nifsinfo=nifs$nifsinfo.cl
task    nifsexamples=nifs$nifsexamples.cl

clbye()
