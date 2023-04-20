# Copyright(c) 2000-2017 Association of Universities for Research in Astronomy, Inc.
#
# Package script for the gemtools package
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
#          Apr  21, 2017 KL     Commissioning release v1.14comm
#          Jul  20, 2017 KL     Release v1.14
#
# load necessary packages - if anything else than loaded
# by gemini should be necessary
specred

# gemseeing requires noao.obsutil.psfmeasure (noao loaded in gemini.cl)
obsutil

package gemtools

# Co-addition of direct images, same task in all packages
task    imcoadd=gemtools$imcoadd.cl

# Co-addition of spectra aperture by aperture
task    gemscombine=gemtools$gemscombine.cl

# MONGO like graphs
task    mgograph=gemtools$mgograph.cl

# Seeing determination
task    gemseeing=gemtools$gemseeing.cl

# Generic arithmetic for MEF images/spectra
# MEF combine
task    gemcombine=gemtools$gemcombine.cl

# DQ decomposition
task    gemdqexpand=gemtools$gemdqexpand.cl

# Headers
task    gemhead=gemtools$gemhead.cl

# MEF handling 
task    wmef=gemtools$wmef.cl
task    gemwcscopy=gemtools$gemwcscopy.cl

task    date=gemtools$date.cl
hidetask('date')

hidetask("gemwcscopy")

task    gemvsample=gemtools$gemvsample.cl
task    gemhedit=gemtools$gemhedit.cl
task    gemqa=gemtools$gemqa.cl
task    gemlogname=gemtools$gemlogname.cl
task    gemoffsetlist=gemtools$gemoffsetlist.cl
task    gemlist=gemtools$gemlist.cl
task    ckinput=gemtools$ckinput.cl
task    ckcal=gemtools$ckcal.cl
task    addbpm=gemtools$addbpm.cl
task    gemcrspec=gemtools$gemcrspec.cl
task    gemfix=gemtools$gemfix.cl
task    growdq=gemtools$growdq.cl

# Hidden tasks
task    gemdate=gemtools$gemdate.cl
task    getfakeUT=gemtools$getfakeUT.cl    # UT date string for data set names
task    gextverify=gemtools$gextverify.cl
task    gimverify=gemtools$gimverify.cl
task    gemsecchk=gemtools$gemsecchk.cl
task    gsetsec=gemtools$gsetsec.cl
task    printlog=gemtools$printlog.cl

hidetask ("gemdate","getfakeUT","gextverify","gimverify","gemsecchk", \
    "gsetsec","printlog")

# Compiled tasks 
task    cnvtsec     = "gemtools$x_gemtools.e"
reset   gemcube     = "gemtools$pkg/gemcube/"
task    gemcube     = "gemtools$x_gemtools.e"
task    gemextn     = "gemtools$x_gemtools.e"
task    gemisnumber = "gemtools$x_gemtools.e"
task    gfwcs       = "gemcube$gflib/gfwcs.par"

hidetask("gfwcs")

task    gloginit    = "gemtools$x_gemtools.e"
task    glogprint   = "gemtools$x_gemtools.e"
task    glogclose   = "gemtools$x_gemtools.e"
task    glogextract = "gemtools$x_gemtools.e"
task    glogfix     = "gemtools$x_gemtools.e"
task    glogpars    = "gemtools$glogpars.par"
task    ldisplay    = "gemtools$x_gemtools.e"

hidetask("gloginit","glogprint","glogclose")
hidetask("cnvtsec","gemisnumber")

# GEMEXPR compiled into gemexpr
task    gemarith    = "gemtools$x_gemexpr.e"
task    gemexpr     = "gemtools$x_gemexpr.e"
task    gemexprpars = "gemtools$gemexprpars.par"
task    mimexprpars = "gemtools$mimexprpars.par"

hidetask("gemexprpars","mimexprpars")

clbye()
