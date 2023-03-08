# Copyright(c) 2003-2017 Association of Universities for Research in Astronomy, Inc.
#
# Package script for the midir package
#
# Version: Jan   9, 2004 KL     Release v1.5
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
onedspec
gnirs

print ""
print "Loading the gnirs package:"
? gnirs
print ""
print "Loading the midir package:"

package midir


# Michelle specific tasks
task      mprepare=midir$mprepare.cl
task      mview=midir$mview.cl
task      mcheckheader=midir$mcheckheader.cl

# TReCS specific tasks
task      tprepare=midir$tprepare.cl
task      tbackground=midir$tbackground.cl
task      tview=midir$tview.cl
task      tcheckstructure=midir$tcheckstructure.cl

# Imaging tasks
task      miview=midir$miview.cl
task      miflat=midir$miflat.cl
task      mistack=midir$mistack.cl
task      miregister=midir$miregister.cl
task      mireduce=midir$mireduce.cl
task      miclean=midir$miclean.cl
task      mipsf=midir$mipsf.cl
task      mistdflux=midir$mistdflux.cl

# Michelle polarimetry reductions
task      mipsplit=midir$mipsplit.cl
task      mipstack=midir$mipstack.cl
task      mipstokes=midir$mipstokes.cl
task      miptrans=midir$miptrans.cl
task      mipsstk=midir$mipsstk.cl
task      mipql=midir$mipql.cl

# Spectroscopy tasks
task      msreduce=midir$msreduce.cl
task      mstelluric=midir$mstelluric.cl
task      msabsflux=midir$msabsflux.cl
task      msdefringe=midir$msdefringe.cl
task      msflatcor=midir$msflatcor.cl
task      msslice=midir$msslice.cl

# Cook book tasks
task      midirinfo=midir$midirinfo.cl
task      midirexamples=midir$midirexamples.cl

# Local tasks

clbye()
