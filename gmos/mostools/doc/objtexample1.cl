# Copyright(c) 2001-2005 Association of Universities for Research in Astronomy, Inc.
#
# ----------------------------------------------------------------------
# File: objtexample1.cl
# Author: Inger Jorgensen, Gemini Observatory
# Date: November 27, 2001
#       February  1, 2002  version for release 1.3
# ----------------------------------------------------------------------
#
# Example IRAF script which shows how to make a valid Object Table
# from an output file from apphot, allstar or nstar.
# The output file from daofind may also be converted to a valid Object
# Table using app2objt.
# The data used in this example are GMOS SV data of the galaxy  
# cluster RXJ0142+2131. 
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Example on how to turn output from apphot.phot into a valid
# Object Table 
# ----------------------------------------------------------------------

# Find the objects using daofind
daofind mrgN20011021S104_add[1] output=mrgN20011021S104.coo \
fwhm=12 threshold=100 verify- ccdread="RDNOISE" gain="GAIN" sigma=14.

# Aperture photometry using apphot.phot
photpars.apertures="10"
phot mrgN20011021S104_add[1] coords=mrgN20011021S104.coo \
output=mrgN20011021S104.mag  ccdread="RDNOISE" gain="GAIN" \
sigma=14. verify- verbose- inter-

# Turn the output from apphot.phot into a valid Object Table using app2objt
# app2objt will remove any objects with mag=INDEF since the GMOS Mask Making
# Software cannot handle these values
delete mrgN20011021S104_OT.fits ver-
app2objt mrgN20011021S104.mag verbose+ image=mrgN20011021S104_add.fits priority="2" 

# The resulting Object Table mrgN20011021S104_OT.fits can be loaded into the
# GMOS Mask Making Software

flpr

# ----------------------------------------------------------------------
# Example on how to turn output from daophot.allstar into a valid
# Object Table.  
# Obviously allstar is not the optimal photometry program to use on this 
# galaxy cluster field. However, this example is just intended to show the 
# different steps required.
# The example requires that daofind and phot from the previous example 
# have been run.
# ----------------------------------------------------------------------

# Load the daophot package
daophot

# Make the PSF - interactive use
# Because of the way the daophot package handles MEF files the output names 
# have the root name mrgN20011021S104_add1
psf mrgN20011021S104_add[1] photfile=mrgN20011021S104.mag pstfile="" \
psfimage=default opstfile=default groupfile=default ccdread="RDNOISE" gain="GAIN" \
verify-

# Run allstar on the image
allstar mrgN20011021S104_add[1] photfile=mrgN20011021S104.mag \
psfimage="default" rejfile="default" subimag="default" \
allstarfile=mrgN20011021S104.all ccdread="RDNOISE" gain="GAIN" sigma=14. \
verify- verbose+ maxiter=5

# Turn the output from daophot.allstar into a valid Object Table using app2objt
delete mrgN20011021S104_OT.fits ver-
app2objt mrgN20011021S104.all verbose+ image=mrgN20011021S104_add.fits priority="2" 

# The resulting Object Table mrgN20011021S104_OT.fits can be loaded into the
# GMOS Mask Making Software

flpr

# ----------------------------------------------------------------------
# Example on how to turn output from daophot.nstar into a valid
# Object Table.  
# Obviously nstar is not the optimal photometry program to use on this 
# galaxy cluster field. However, this example is just intended to show the 
# different steps required.
# The example requires that daofind, phot and psf from the previous 
# examples have been run.
# ----------------------------------------------------------------------

# Run group on the image
group mrgN20011021S104_add[1] photfile=mrgN20011021S104.mag \
psfimage="default" groupfile="default" \
ccdread="RDNOISE" gain="GAIN" sigma=14. verify- verbose+ maxiter=5

# Run nstar on the image
nstar mrgN20011021S104_add[1] groupfile="default" \
psfimage="default" nstarfile=mrgN20011021S104.nst rejfile="default" \
ccdread="RDNOISE" gain="GAIN" sigma=14. verify- verbose+ maxiter=5

# Turn the output from daophot.nstar into a valid Object Table using app2objt
delete mrgN20011021S104_OT.fits ver-
app2objt mrgN20011021S104.nst verbose+ image=mrgN20011021S104_add.fits priority="2" 
