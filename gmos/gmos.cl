# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.
#
# Package script for the gmos package
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
#          Jul  22, 2014 KL     Commissioningn release v1.13 GMOS Ham
#          Jan  30, 2015 KL     Release v1.13
#	   Dec   7, 2015 KL	Release v1.13.1
#          Apr  21, 2017 KL     Commisioning release v1.14comm
#          Jul  20, 2017 KL     Release v1.14
#
# load necessary packages - gemini loads most of the packages
gemtools

# spectroscopy packages for longslit - may be moved to gemini later
twodspec
longslit

package gmos

# mostools
set  mostools      = "gmos$mostools/"
task mostools.pkg  = mostools$mostools.cl

# Needed for IFU headers
set min_lenuserarea = 400000

# Setup tasks and hartmann task

# Generic tasks
task gdisplay=gmos$gdisplay.cl
task gprepare=gmos$gprepare.cl
task gbias=gmos$gbias.cl
task gbpm=gmos$gbpm.cl
task gmosaic=gmos$gmosaic.cl
task ggain=gmos$ggain.cl
task ggdbhelper=gmos$ggdbhelper.cl
task gmultiamp=gmos$gmultiamp.cl
task goversub=gmos$goversub.cl
task gqecorr=gmos$gqecorr.cl
task gretroi=gmos$gretroi.cl
task gsat=gmos$gsat.cl
task gtile=gmos$gtile.cl
hidetask ggain ggdbhelper gsat gtile gretroi
hidetask gmultiamp goversub

# Imaging tasks
task giflat=gmos$giflat.cl
task gifringe=gmos$gifringe.cl
task girmfringe=gmos$girmfringe.cl
task gireduce=gmos$gireduce.cl

# Spectroscopy tasks
task gsappwave=gmos$gsappwave.cl
task gscalibrate=gmos$gscalibrate.cl
task gscrrej=gmos$gscrrej.cl
task gscrmask=gmos$gscrmask.cl
task gscut=gmos$gscut.cl
task gsdrawslits=gmos$gsdrawslits.cl
task gsextract=gmos$gsextract.cl
task gsflat=gmos$gsflat.cl
task gsreduce=gmos$gsreduce.cl
task gsscatsub=gmos$gsscatsub.cl
task gsskysub=gmos$gsskysub.cl
task gsstandard=gmos$gsstandard.cl
task gstransform=gmos$gstransform.cl
task gswavelength=gmos$gswavelength.cl

# N&S tasks
task gnsskysub=gmos$gnsskysub.cl
task gnscombine=gmos$gnscombine.cl
task gnsdark=gmos$gnsdark.cl

# Special IFU tasks
task gfapsum=gmos$gfapsum.cl
task gfscatsub=gmos$gfscatsub.cl
task gfcube=gmos$x_gmos.e
task gfdisplay=gmos$gfdisplay.cl
task gfextract=gmos$gfextract.cl
task gffindblocks=gmos$gffindblocks.cl
task gfquick=gmos$gfquick.cl
task gfreduce=gmos$gfreduce.cl
task gfresponse=gmos$gfresponse.cl
task gfskysub=gmos$gfskysub.cl
task gftransform=gmos$gftransform.cl
task gfunexwl=gmos$x_gmos.e
hidetask gfunexwl

# Cookbooks
task gmosinfo=gmos$gmosinfo.cl
task gmosinfoimag=gmos$gmosinfoimag.cl
task gmosinfospec=gmos$gmosinfospec.cl
task gmosinfoifu=gmos$gmosinfoifu.cl

# Examples
task gmosexamples=gmos$gmosexamples.cl

clbye()
