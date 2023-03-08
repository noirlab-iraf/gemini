# Copyright(c) 2001-2012 Association of Universities for Research in Astronomy, Inc.
#
# ----------------------------------------------------------------------
# File: objtexample2.cl
# Author: Inger Jorgensen, Gemini Observatory
# Date: November 26, 2001
#       February  1, 2002  version for release 1.3
# ----------------------------------------------------------------------
#
# Example IRAF script which shows how to make a valid Object Table
# from one or more output tables from SExtractor
# The data used in this example are GMOS SV data of the galaxy  
# cluster RXJ0142+2131. 
# ----------------------------------------------------------------------

# Before running any commands, please make sure the 'tbltype' parameter for
# the tasks: tcreate and tmerge are set to "default" otherwise there may be
# issues with the ".fits" extensions.

# Run SExtractor with output to simple FITS table
# Do not include any arrays in the selected output parameters
# In this example the output tables are
# mrgN20011021S104tab.fits   Photometry in the r-filter
# mrgN20011021S115tab.fits   Photometry in the g-filter
# mrgN20011021S112tab.fits   Photometry in the i-filter

# Example on how to get rid of some of the non-real objects 
# at the edges of the image
tselect mrgN20011021S104tab.fits RXJ_r.fits \
"x_image>975 && x_image<5325 && y_image>75 && y_image<4450"
tselect mrgN20011021S112tab.fits RXJ_i.fits \
"x_image>905 && x_image<5255 && y_image>75 && y_image<4450"
tselect mrgN20011021S115tab.fits RXJ_g.fits \
"x_image>905 && x_image<5255 && y_image>145 && y_image<4520"

# RA in hours, DEC in deg  -- this step is required
tcalc RXJ_r.fits "RA_J2000" "alpha_J2000/15." colfmt="%12.2h" colunit="H"
tcalc RXJ_r.fits "DEC_J2000" "delta_J2000" colfmt="%12.2h" colunit="deg"

# Turn one of these tables into a valid Object Table using stsdas2objt
# The task stsdas2objt needs to be declared first
task stsdas2objt=stsdas2objt.cl
delete RXJ_r_OT.fits ver-
stsdas2objt RXJ_r id_col=number x_col=x_image y_col=y_image \
mag_col=mag_best ra_col=RA_J2000 dec_col=DEC_J2000 \
verbose+ image=mrgN20011021S104_add.fits priority="2" \
other_col="class_star,theta_image,ellipticity" 

# The resulting Object Table RXJ_r_OT.fits can be loaded into the
# GMOS Mask Making Software

flpr

# ----------------------------------------------------------------------
# The rest of the example deals with how to merge the photometry
# in the 3 filters, set priorites and set slit length and width
# ----------------------------------------------------------------------

# Sort the tables to make the matching of the tables a bit faster
tsort  RXJ_r.fits alpha_J2000
tsort  RXJ_i.fits alpha_J2000
tsort  RXJ_g.fits alpha_J2000

# Change column names to get unique names. 
# The filter name is appended to each of the column names from the table
# matching the filter
tlcol  RXJ_i.fits > tmpdat
list="tmpdat"
while(fscan(list,s1)!=EOF) {
  tchcol("RXJ_r.fits",s1,s1//"r","","")
  tchcol("RXJ_g.fits",s1,s1//"g","","")
  tchcol("RXJ_i.fits",s1,s1//"i","","")
}
list="" ; delete tmpdat ver-

# Match r-filter and i-filter photometry with a 1.0" tolerance on object positions. 
# tmpri.diag contains info on double matches and unmatched objects.
tmatch RXJ_r.fits RXJ_i.fits tmpri.fits alpha_J2000r,delta_J2000r \
alpha_J2000i,delta_J2000i maxnorm=0.00027777778 sphere+ diag=tmpri.diag \
nmcol1="numberr,alpha_J2000r,delta_J2000r,mag_bestr,class_starr" \
nmcol2="numberi,alpha_J2000i,delta_J2000i,mag_besti,class_stari"

# Match the common objects with the g-filter photometry with a 1.0" tolerance  
# on object positions
# tmpgri.diag contains info on double matches and unmatched objects.
tmatch tmpri.fits RXJ_g.fits RXJ0142p2131.fits alpha_J2000r,delta_J2000r \
alpha_J2000g,delta_J2000g maxnorm=0.00027777778 sphere+ diag=tmpgri.diag \
nmcol1="numberr,alpha_J2000r,delta_J2000r,mag_bestr,class_starr" \
nmcol2="numberg,alpha_J2000g,delta_J2000g,mag_bestg,class_starg"

# The table RXJ0142p2131.fits now contains all objects detected in all 3 bands

# Separate stars and galaxies
tselect RXJ0142p2131.fits RXJ0142p2131gal.fits "class_starr<0.95"
tselect RXJ0142p2131.fits RXJ0142p2131star.fits "class_starr>=0.95"

# Select the sample: Require an r-magnitude brighter than 21.6 mag
tselect  RXJ0142p2131gal.fits  RXJ0142p2131sample.fits "mag_bestr<=21.6"

# Select a few stars as acquisition objects and merge those with the sample objects
tselect RXJ0142p2131star.fits tmpstar.fits "numberr==578 || numberr==949 || numberr==2824"
tmerge  RXJ0142p2131sample.fits,tmpstar.fits RXJ0142p2131out.fits append

# RA in hours, DEC in deg -- this step is required
tcalc RXJ0142p2131out.fits "RA_J2000" "alpha_J2000r/15." colfmt="%12.2h" colunit="H"
tcalc RXJ0142p2131out.fits "DEC_J2000" "delta_J2000r" colfmt="%12.2h" colunit="deg"

# Calculate colors from the aperture magnitudes
tcalc RXJ0142p2131out.fits g_r "mag_aperg-mag_aperr" colfmt="f6.3"
tcalc RXJ0142p2131out.fits g_i "mag_aperg-mag_aperi" colfmt="f6.3"
tcalc RXJ0142p2131out.fits r_i "mag_aperr-mag_aperi" colfmt="f6.3"

# PA with the correct conventions, limit the slittilt to 30 deg
tcalc RXJ0142p2131out.fits "slittilt" "theta_imager-90." colfmt="f6.1"
tcalc RXJ0142p2131out.fits "slittilt" "if slittilt<-90 then slittilt+180. else slittilt"
tcalc RXJ0142p2131out.fits "slittilt" "if class_starr>=0.95 then 0 else slittilt"
tcalc RXJ0142p2131out.fits "slittilt" "if abs(slittilt)>30 then 0 else slittilt"

# Derive approximative object size and set slit length (slitsize_y)
tcalc RXJ0142p2131out.fits "xsizer" "xmax_imager-xmin_imager" colfmt="f4.0"
tcalc RXJ0142p2131out.fits "ysizer" "ymax_imager-ymin_imager" colfmt="f4.0"
tcalc RXJ0142p2131out.fits "sizer" "sqrt(xsizer**2+ysizer**2)" colfmt="f6.1"
tcalc RXJ0142p2131out.fits "slitsize_y" "min(max(0.8*sizer*0.0727,5.),10.)" \
colfmt="f6.1" colunit="arcsec"
# Adjust slitsize_y for the brightest and largest objects
tcalc RXJ0142p2131out.fits "slitsize_y" "if mag_bestr<17.6 then 15. else slitsize_y"
# Set slit width (slitsize_x)
tcalc RXJ0142p2131out.fits "slitsize_x" "0.75" colfmt="f6.2" colunit="arcsec"

# Add a priority column based on the magnitude in the r-filter and the position on
# the detector array
tcalc RXJ0142p2131out.fits "userprior" "2" colfmt="i3"
tcalc RXJ0142p2131out.fits "userprior" \
"if userprior<3 && mag_bestr<=19.5 && x_imager>2110 && x_imager<4110 then 1 else userprior"
tcalc RXJ0142p2131out.fits "userprior" "if x_imager<1510 || x_imager>4710 then 3 else userprior"
tcalc RXJ0142p2131out.fits "userprior" "if class_starr>0.95 then 0 else userprior"
tcalc RXJ0142p2131out.fits "userprior" "if userprior>0 && mag_bestr>21.3 then 3 else userprior"

# Convert userprior column to a format that can be used by the GMOS Mask Making Software
tprint RXJ0142p2131out.fits col="userprior" showh- showr- > tmpdat
# Any objects with priority 4 are set to "omitted" - the current example has none of these
translit tmpdat "." " " | translit STDIN "4" "X" > tmpdat2   
print("priority ch*1 %d\n") | \
tcreate("tmppri.fits","STDIN","tmpdat2",hist-,tbltype="default")
tmerge RXJ0142p2131out.fits,tmppri.fits tmpout.fits merge
!\mv tmpout.fits RXJ0142p2131out.fits
delete tmp* ver-

flpr
# Example of adjustments to the table
# Use slitpos_y to offset the two central galaxies a bit in the slit to get them both
tcalc RXJ0142p2131out.fits slitpos_y "if numberr==801 then -1 else 0" colfmt="f6.2" colunit="arcsec"
tcalc RXJ0142p2131out.fits slitpos_y "if numberr==802 then 1 else slitpos_y" 

# Convert this table to an Object Table using stsdas2objt.cl
# The task stsdas2objt needs to be declared first, if this has not already been done
task stsdas2objt=stsdas2objt.cl
delete RXJ0142p2131out_OT.fits ver-
stsdas2objt RXJ0142p2131out id_col=numberr x_col=x_imager y_col=y_imager \
mag_col=mag_bestr ra_col=RA_J2000 dec_col=DEC_J2000 \
verbose+ image=mrgN20011021S104_add.fits priority="none" \
other_col="slitsize_y,slitsize_x,slittilt,slitpos_y,class_starr,mag_bestg,mag_besti,g_i,r_i,g_r,\
theta_imager,ellipticityr,priority" 

# Sort the file according to the object ID - makes it easier to navigate
# in the GMOS Mask Making Software
tsort RXJ0142p2131out_OT.fits id

# List the column names in the Object Table
tlcol RXJ0142p2131out_OT.fits

# Print the Object Table to the screen, this includes printing the header parameters
tprint RXJ0142p2131out_OT.fits prpar+

# The Object Table RXJ0142p2131out_OT.fits is now ready to be loaded into the
# GMOS Mask Making Software

