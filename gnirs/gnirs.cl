# Copyright(c) 2003-2017 Association of Universities for Research in Astronomy, Inc.
#
# Package script for the gnirs package
#
# Version: Oct  25, 2004 KL     Release v1.7
#          May   6, 2005 KL     Release v1.8
#          Jul  28, 2006 KL     Release v1.9
#          Jul  28, 2009 JH     Release v1.10
#          Jan  13, 2011 EH,KL  Beta release v1.11beta
#          Aug  19, 2011 EH     Beta release v1.11beta2
#          Dec  30, 2011 EH     Release v1.11
#          Mar  28, 2012 EH     Release v1.11.1
#          Dec  13, 2012 EH     Beta release v1.12beta
#          May  14, 2013 EH     Beta release v1.12beta2
#          Oct  11, 2013 EH     Release v1.12
#          Jan  30, 2015 KL     Release v1.13
#	   Dec   7, 2015 KL	Release v1.13.1
#          Jul  20, 2017 KL     Release v1.14
#
# load necessary packages - if anything else than loaded
# by gemini should be necessary
gemtools
twodspec
longslit
imred # needed for irred
crutil # needed for nsedge
irred # needed for irlincor in nsreduce
st4gem.analysis.fourier # needed for nsoffset

package gnirs

# Spectroscopy tasks
task nsheaders=gnirs$nsheaders.cl
task nsmdfhelper=gnirs$nsmdfhelper.cl
task nsprepare=gnirs$nsprepare.cl
task nsflat=gnirs$nsflat.cl
task nsedge=gnirs$nsedge.cl
task nsslitfunction=gnirs$nsslitfunction.cl
task nsreduce=gnirs$nsreduce.cl
task nscut=gnirs$nscut.cl
task nssky=gnirs$nssky.cl
task nsressky=gnirs$nsressky.cl
task nssdist=gnirs$nssdist.cl
task nswavelength=gnirs$nswavelength.cl
task nswhelper=gnirs$nswhelper.cl
task nsappwave=gnirs$nsappwave.cl
task nswedit=gnirs$nswedit.cl
task nsfitcoords=gnirs$nsfitcoords.cl
task nstransform=gnirs$nstransform.cl
task nsextract=gnirs$nsextract.cl
task nscombine=gnirs$nscombine.cl
task nschelper=gnirs$nschelper.cl
task nsstack=gnirs$nsstack.cl
task nsoffset=gnirs$nsoffset.cl
task nfquick=gnirs$nfquick.cl
task nfflt2pin=gnirs$nfflt2pin.cl
task peakhelper=gnirs$x_gnirs.e
task nfcube=gnirs$nfcube.cl
task nxdisplay=gnirs$nxdisplay.cl
task nvnoise=gnirs$nvnoise.cl
task nstelluric=gnirs$nstelluric.cl

# Cookbooks
task gnirsinfo=gnirs$gnirsinfo.cl
task gnirsinfols=gnirs$gnirsinfols.cl
task gnirsinfoxd=gnirs$gnirsinfoxd.cl
task gnirsinfoifu=gnirs$gnirsinfoifu.cl

# Examples
task gnirsexamples=gnirs$gnirsexamples.cl

hidetask nsmdfhelper
hidetask nswhelper
hidetask nschelper
hidetask nssky
hidetask peakhelper

clbye()
