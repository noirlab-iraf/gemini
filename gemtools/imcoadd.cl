# Copyright(c) 1994-2005 Inger Jorgensen
# Copyright(c) 2000-2017 Association of Universities for Research in Astronomy, Inc.

procedure imcoadd(images)

# Combine science images, includes cleaning for cosmic ray events
# This script was originally written for ground based optical CCD data
# It should work on any digital linear imaging data.
#
# First image in the list is used as reference image 'ref'
#
# if fl_find==yes daofind is used for identifying objects in 'ref'
#   else coordinate file 'ref'_pos
#
# if fl_map==yes the transformation is determined:
#    One or two common object(s) is/are marked on each image
#    Each transformation is determined interactively using geomap
# else name convention used for input files for transformation
#
# if fl_trn ==yes the transformation is made ==> 'image'_trn
# if fl_med ==yes  'ref' + SUM 'image'_trn = 'ref'_med
# if fl_add ==yes  'ref'_add = mean of masked images
# if fl_avg ==yes  'ref'_avg = mean of images, no masking
# 
# NAME CONVENTION
#
# 'ref'_pos     : daofind file for 'ref'
# 'ref'_cen     : centered positions from apphot.center
# 'ref'_mag.tab : STSDAS table with magnitudes used for scaling
# 'image'_pos   : approximate positions using x_offset,y_offset
# 'image'_cen   : centered positions from apphot.center
# 'image'_trn   : position input file for geomap, entry in database
# 'image'_trn   : transformed image
# 'ref'_med     : median filtered combined image of all input images
# 'ref'_add     : mean image cleaned for CR, scaled and weighted
#                 with relint if fl_scale=yes
# 'ref'_avg     : mean image, NOT cleaned for CR, scaled and weighted
#                 with relint if fl_scale=yes
#
# If you encounter unexplained problems, unlearn critical tasks
# geomap, geotran, all tasks related to apphot
#
# Version  Oct 11, 2000  IJ
#          May 15, 2001  IJ, MEF implementation
#          May 17, 2001  IJ, Fixpix all input files if DQ,immasks or
#                        badpix file available
#                        Use DQ as immasks if available (default yes)
#          May 21, 2001  IJ, really made it use DQ
#          May 31, 2001  IJ, changed gimverify call to use new version of gimverify
#          Jun  7, 2001  IJ, trapped failure to center objects
#                            trapped no objects in one or more of *_trn
#                            fixed rejection of all objects if rms=0.
#                            delete left over tmp images badpix.{fits,imh,hhh}
#          Jun  8, 2001  IJ, exit if fl_obj=no alignmethod!=header
#          Aug 30, 2001  IJ, changed badpix growth to prevent bus errors
#                            for large images
#                            fixed use of individual *.pl for ref image
#          Nov 1, 2001   BM, generalized for images with different gain/ron
#          Nov 6, 2001   BM, bug fix
#          Feb 28, 2002  IJ  v1.3 release
#          May 15, 2002  IJ  fix deletion of input masks for MEF not using DQ
#          Aug 21, 2002  IJ  alignmethod=header generalized + automatic for
#                            H/Q, NIRI and GMOS data, badpixfile=default
#                            status, cleaned code, parameter encoding
#                            explicit output image name for ref_add
#                            datamax from key_sat if present, stable for few obj
#          Aug 22, 2002  IJ  alignmethod=wcs implemented, revised align=header
#                            to use real instrument offsets only, cleaned code
#          Aug 26, 2002  IJ  wcs for rotate/scale more reliable, imtype!=fits for MEF caught
#          Sept 25, 2002 IJ  wcs for no objects should be allowed
#                            handle limit correctly for sigma limits (bug was introduced
#                            Nov 1, 2001 with the handling of different ron/gain)
#          Sept 25, 2002 IJ  v1.4 release
#          Dec 18, 2002  IJ  make sure nmisc.fixpix is picked up
#          Feb 17, 2003  BM  generalize for GMOS-S
#          Apr 10, 2003  BM  change default sign convention for GMOS-S to get header option to work, movefiles->rename for OLDP
#          Jun 3, 2003   IJ  GMOS=GMOS-N to have support for old data
#          Aug 13, 2003  KL  IRAF2.12 - list of tasks w/ new parameters
#                              hedit: addonly
#                              imstat: nclip,lsigma,usigma,cache
#                              imcombine: rejmask->rejmasks, plfile->nrejmasks
#                                         header,bpmask,expmask,outlimits
#                              txdump: replaced by pdump
#                              phot & center: wcsin,wscout,cache
#                              daofind: wscin,cache
#                              geomap: maxiter
#  
# Reference: Jorgensen I., 2003, AJ, submitted
#
# You may modify this script, in which case the comments in the
# final script should contain the two copyright lines in the start
# of the script.
#

char images {prompt="Images to be combined, first image as reference"}     # OLDP-1-input-primary-combine-suffix=_add
char outimage {"",prompt="Output image name if not derived from reference"} # OLDP-1-output
char sci_ext {"SCI",prompt="Science extension name if MEF input"}          # OLDP-3
char var_ext {"VAR",prompt="Variance extension name if MEF input"}         # OLDP-3
char dq_ext {"DQ",prompt="Data quality extension name if MEF input"}       # OLDP-3
char immasks  {"DQ",prompt="Masks of bad pixels for each image"}           # OLDP-3-xinput
char database {"imcoadd.dat",prompt="Database file for transformations"}   # OLDP-3
real threshold {20.,prompt="Threshold for daofind, sigma above sky"}       # OLDP-2
real fwhm {7.,min=1,prompt="FWHMPSF for daofind, in pixels (ref)"}         # OLDP-2
real box {20.,min=5,prompt="Box size for centering (ref)"}                 # OLDP-2
char alignmethod{"wcs",enum="wcs|user|twodx|header",prompt="Method for rough alignment of images (wcs|user|twodx|header)"} # OLDP-2
char asection{"default",prompt="Image sections (2) (alignmethod=twodx)"}   # OLDP-3
int  xwindow {181,min=11,prompt="Window size to search if using twodx"}    # OLDP-3
char key_inst {"INSTRUME",prompt="Header keyword for instrument"}          # OLDP-3
char key_camera {"CAMERA",prompt="Header keyword for camera (NIRI only)"}  # OLDP-3
char key_inport {"INPORT",prompt="Header keyword for ISS port"}            # OLDP-3
char key_xoff {"default",prompt="Header keyword for instrument X-offset (alignmethod=header)"}  # OLDP-3
char key_yoff {"default",prompt="Header keyword for instrument Y-offset (alignmethod=header)"}  # OLDP-3
real instscale {1.,min=0.,prompt="Scale for offsets arcsec/unit (alignmethod=header)"}          # OLDP-3
char xsign {"default",enum="default|negative|positive",prompt="XOFFSET sign relative to detector (alignmethod=header)"} # OLDP-3
char ysign {"default",enum="default|negative|positive",prompt="YOFFSET sign relative to detector (alignmethod=header)"} # OLDP-3
char key_pixscale {"PIXSCALE",prompt="Header keyword for pixel scale in arcsec/pixel"} # OLDP-3
real pixscale {1.,min=0.000001,prompt="Pixel scale arcsec/pixel"}          # OLDP-3
real dispmag {1.0,prompt="Magnify for display"}                            # OLDP-3
bool rotate {no,prompt="Assume rotation > 0.5 deg"}                        # OLDP-3
bool scale  {no,prompt="Allow significant scale differences"}              # OLDP-3
char geofitgeom{"rscale",enum="shift|xyscale|rotate|rscale|rxyscale|general",prompt="Fitting geometry for geomap"} # OLDP-3
int  order{3,min=2,prompt="Fitting order for geomap (geofitgeom=general)"} # OLDP-3
real sigfit{2.5,min=0,prompt="Sigma rejection for geomap"}                 # OLDP-3
int  niter{5,min=0,prompt="Iterations for geomap"}                         # OLDP-3
real coolimit{0.3,min=0,prompt="Sigma in pixels for converged geomap"}     # OLDP-3
char geointer {"linear",enum="nearest|linear|poly3|poly5|spline3",prompt="Interpolation for geotran"} # OLDP-2
int  geonxblock {2048,prompt="geotran:X dim of working block size in pixels"} # OLDP-4
int  geonyblock {2048,prompt="geotran:Y dim of working block size in pixels"} # OLDP-4
char key_ron {"RDNOISE",prompt="Header keyword for read noise in electrons"}  # OLDP-3
char key_gain {"GAIN",prompt="Header keyqword for gain in electrons/ADU"}  # OLDP-3
real ron {1.,prompt="Read out noise [electrons]"}                          # OLDP-3
real gain {1.,prompt="Gain [electrons/ADU]"}                               # OLDP-3
real datamin {-1000.,prompt="Minimum data value below sky for good pixels"} # OLDP-3
char key_sat {"SATURATI",prompt="Header keyword for maximum data value"}   # OLDP-3
real datamax {50000.,prompt="Maximum data value for good pixels"}          # OLDP-3
real aperture {30.,prompt="Aperture radius for magnitudes for scaling"}    # OLDP-3
real limit {15.,prompt="Intensity limit for CR cleaning"}                  # OLDP-3
bool key_limit {yes,prompt="Is limit given in sigmas above sky?"}          # OLDP-3
real lowsigma {7.,prompt="Sigma rejection limit for CRs"}                  # OLDP-3
real lowlimit {500.,prompt="Absolute rejection limit for CRs"}             # OLDP-3
real scalenoise {0.,prompt="Noise term proportional to signal above sky"}  # OLDP-3
int  growthrad {1,min=0,max=1,prompt="Growth radius for CR, 0 or 1"}       # OLDP-3
char statsec {"default",prompt="Statistics section for sky determination"} # OLDP-3
char badpixfile {"default",prompt="Badpixel file or image"}                # OLDP-2-xinput
char dq_maskexpr {"default", prompt="Expression to use to create mask from input DQ planes (MEFs only)."}
bool fl_inter {no,prompt="Inspect fits interactively"}                     # OLDP-3
bool fl_refmark {no,prompt="Mark objects on display for ref image"}        # OLDP-3
bool fl_mark {no,prompt="Mark objects on display for all images"}          # OLDP-3
bool fl_fixpix {yes,prompt="Fixpix input images using DQ (MEF only)"}      # OLDP-3
bool fl_find {yes,prompt="Find objects with daofind"}                      # OLDP-3
bool fl_map {yes,prompt="Determine the transformation"}                    # OLDP-3
bool fl_trn {yes,prompt="Make transformation and subtract sky"}            # OLDP-3
bool fl_med {yes,prompt="Calculate median image"}                          # OLDP-3
bool fl_add {yes,prompt="Calculate mean of masked images"}                 # OLDP-3
bool fl_avg {no,prompt="Calculate mean of images, no masking"}             # OLDP-3
bool fl_scale {yes,prompt="Combine mean images scaled with RELINT"}        # OLDP-3
bool fl_overwrite {no,prompt="Overwrite previous results"}                 # OLDP-3
char logfile {"imcoadd.log",prompt="Logfile"}                              # OLDP-1
bool verbose {yes,prompt="Verbose extended output"}                        # OLDP-4
int  status {0,prompt="Exit status (0=good)"}                              # OLDP-4
struct *inimag {prompt="List-directed struct, do not change"}              # OLDP-4
struct *inmask {prompt="List-directed struct, do not change"}              # OLDP-4
struct *input  {prompt="List-directed struct, do not change"}              # OLDP-4

begin

char n_ref, n_images, n_coords, n_database, n_realinim, n_inim, n_outimage
char l_sci_ext, l_var_ext, l_dq_ext, l_zero_ext, n_sci_ext, wcs_header
char n_statsec, n_sigim, n_mask, n_badpixfile, n_immasks, n_geointer
char n_geofitgeom, n_alignmethod, n_alignsec
int  n_xwindow
char n_key_xoff, n_key_yoff, n_key_pixscale
real n_xsign, n_ysign
char n_key_inst, n_key_camera, n_key_inport, n_inst, n_camera 
int  n_inport
real n_pixscale, n_instscale
real x0, y0, x11, y11, x21, y21, x12, y12, x22, y22, xc1, yc1, xc2, yc2
real x11_ref, y11_ref
real xsize1, ysize1, xsize2, ysize2
real x_offset, y_offset, n_dispmag
real n_thres, n_box, n_fwhm, l_ron, l_gain
real n_ron[500], n_gain[500], roneff, gaineff
char n_key_gain, n_key_ron, n_key_sat
real n_lim, n_datamin, n_datamax, n_aperture, n_sigfit, n_coolimit
real n_lsig, n_llim, n_usig, n_ulim, n_scnoise
int  Xmax, Ymax, n_order, n_files, n_i, n_growthrad, n
int n_xsamp, n_ysamp, l_nxblock, l_nyblock, xbin
bool n_fl_find, n_fl_trn, n_fl_med, n_fl_map, n_fl_add, n_fl_scale
bool n_fl_msk, n_fl_avg, n_rotate, n_scale, n_fl_overwrite, n_fl_inter
bool l_fl_fixpix
bool l_verbose
char l_logfile

real n_medsky, n_insky, n_refint, n_inint, n_totint
real n_sigx, n_sigy, n_boxfac
int  n_lines, n_count, n_iter, n_ordoff, n_sec, n_ii
bool n_interact, n_break, n_badim, fl_obj, l_fl_mef, n_fl_del
char n_imtype, n_str, n_test
char n_correlation

char tmpcen, tmpcoo, tmpref, tmpnew, tmpim, tmplis, tmpout, tmpbad 
char tmpmli, tmpmag, n_m0, n_m1, n_m2, n_m3, n_m4, tmpba2, tmpdq
char tmpres, tmpin, n_tmpbase, tmplas
char tmpshift, tmpsh1, tmpsh2, tmpimlas, l_dq_maskexpr, dqexpr, dq_inst

string dettype, bpmversion
int  namps

struct n_struct

status=0

cache("tstat","imgets","tinfo","gimverify","keypar", "gemdate")
l_logfile=logfile
if(l_logfile=="" || l_logfile==" ") {
  l_logfile="imcoadd.log"
  printlog("WARNING - IMCOADD: Using default logfile imcoadd.log",l_logfile,yes)
}
l_verbose = verbose 

printlog("-----------------------------------------------------------------------------",l_logfile,l_verbose)
date | scan(n_struct)
printlog("IMCOADD -- "//n_struct,l_logfile,l_verbose)
printlog("",l_logfile,l_verbose)

# make sure xray.ximages is NOT loaded as it contains another imcalc
if(defpac("ximages")) {
 printlog("ERROR - IMCOADD: Package ximages should be unloaded before running",l_logfile,yes)
 goto crash
}

# get output imagetype
show("imtype") | scan(n_imtype)
n_imtype="."//n_imtype

n_images = images ; n_immasks = immasks ; n_outimage=outimage
l_sci_ext = "["//sci_ext//"]" ; l_var_ext = "["//var_ext//"]" 
l_dq_ext = "["//dq_ext//"]" ; l_zero_ext = "[0]"
n_sci_ext = sci_ext
n_database = database
n_dispmag = dispmag
n_thres = threshold ; n_box = box ; n_fwhm = fwhm ; n_rotate=rotate
n_scale = scale ; n_geofitgeom=geofitgeom ; n_alignmethod=alignmethod
n_alignsec = asection ; n_xwindow=xwindow
n_order = order ; n_sigfit=sigfit ; n_geointer=geointer ; n_iter=niter
n_coolimit=coolimit
l_fl_fixpix=fl_fixpix ; n_fl_find = fl_find
n_fl_map = fl_map ; n_fl_trn = fl_trn
n_fl_med = fl_med ; n_fl_add = fl_add ; n_fl_avg = fl_avg
n_fl_scale=fl_scale
n_fl_overwrite = fl_overwrite ; n_fl_inter=fl_inter
n_statsec = statsec ; n_badpixfile = badpixfile
n_key_inst = key_inst ; n_key_camera=key_camera ; n_key_inport=key_inport
n_key_gain = key_gain ; n_key_ron = key_ron ; n_key_sat=key_sat
l_gain = gain ; l_ron = ron ; roneff=0.0 ; gaineff=0.0
n_lim = limit ; n_growthrad = growthrad ; n_scnoise=scalenoise
n_lsig = lowsigma ; n_llim=lowlimit 
n_usig = 500000. ; n_ulim=100000000.
n_datamin = datamin ; n_datamax = datamax ; n_aperture = aperture
l_dq_maskexpr = dq_maskexpr

# temporary list of images
tmplis = mktemp("tmplis")   # list of images to combine
tmpmli = mktemp("tmpmli")   # list of masks of images
tmpout = mktemp("tmpout")   # 
tmpbad = mktemp("tmpbad")   # badpixel file
tmpba2 = mktemp("tmpba2")   # badpixel file  no 2, used for input
tmpdq  = mktemp("tmpdq")    # Mask files from DQ (MEF only)
n_sigim = mktemp("tmpsig")  # sigma image for masking
n_mask = mktemp("tmpmsk")   # mask image
n_m0 = mktemp("tmpm0")      # mask image
n_m1 = mktemp("tmpm1")      # mask image
n_m2 = mktemp("tmpm2")      # mask image
n_m3 = mktemp("tmpm3")      # mask image
n_m4 = mktemp("tmpm4")      # mask image
tmpres = mktemp("tmpres")   # results from geomap, iterative run
tmpin  = mktemp("tmpin")    # filtered input for geomap, interative run

tmpcen = mktemp("tmpcen")
tmpcoo = mktemp("tmpcoo")
tmpref = mktemp("tmpref")
tmplas = mktemp("tmplas")   # used in mapping loop, only
tmpimlas = mktemp("tmpiml") # twodx
tmpshift = mktemp("tmpshi") # twodx
tmpsh1 = mktemp("tmpsh1")   # twodx
tmpsh2 = mktemp("tmpsh2")   # twodx
tmpnew = mktemp("tmpnew")
tmpim = mktemp("tmpim")
tmpmag = mktemp("tmpmag")

if(n_alignmethod!="user") {
  fl_mark=no ; fl_refmark=no
}

# n_images is either a list of images separated with commas or
# a list in a file in which case the name starts with @

# make the list of images into a list in a file
if(substr(n_images,1,1)!="@") 
  files(n_images,sort-, > tmplis)
else if (access(substr(n_images,2,strlen(n_images))))
 fields( (substr(n_images,2,strlen(n_images))),"1",lines="1-", quit-, print-) | words((substr(n_images,2,strlen(n_images))), > tmplis)
else {
 printlog("ERROR - IMCOADD: Input file "//substr(n_images,2,strlen(n_images))//" not found",l_logfile,yes)
 goto crash
}

# cut off get any extensions (.fits, .imh) on input images 
inimag=tmplis
while(fscan(inimag,n_inim)!=EOF) {
 if( substr(n_inim,strlen(n_inim)-4,strlen(n_inim))==".fits")
  n_inim=substr(n_inim,1,strlen(n_inim)-5)
 if( substr(n_inim,strlen(n_inim)-3,strlen(n_inim))==".imh")
  n_inim=substr(n_inim,1,strlen(n_inim)-4)
 print(n_inim, >> tmpmli)
}
inimag=""
delete(tmplis,verify-)
rename(tmpmli,tmplis,field="all")

# n_immasks is either a list of masks separated with commas,
# a list in a file in which case the name starts with @, or "DQ"

# make the list of immasks into a list in a file
if(n_immasks!="DQ") {
  if(substr(n_immasks,1,1)!="@")
    files(n_immasks,sort-, > tmpmli)
  else if (access(substr(n_immasks,2,strlen(n_immasks))))
   fields( (substr(n_immasks,2,strlen(n_immasks))),"1",lines="1-", quit-, print-) | words((substr(n_immasks,2,strlen(n_immasks))), > tmpmli)
  else {
   printlog("ERROR - IMCOADD: Input file "//substr(n_immasks,2,strlen(n_immasks))//" not found",l_logfile,yes)
   goto crash
  }
}

count(tmplis) | scan(n_files)
if(n_files < 2) {
  printlog("ERROR - IMCOADD: Only one image in list",l_logfile,yes)
  goto crash
} else if (n_files > 500) {
  printlog("ERROR - IMCOADD: The number of files exceeds the maximum allowable number of 500",l_logfile,yes)
  goto crash
}

# Check mask list
n_i=0 ; n_fl_msk = yes
if(access(tmpmli))
 count(tmpmli) | scan(n_i)
if(n_i==0)
 n_fl_msk=no
if(n_i>0 && n_i!=n_files) {
 printlog("ERROR - IMCOADD: Same number of images and masks required. Use none for no mask",l_logfile,yes)
 goto crash
}

printlog("Images (and masks) in list",l_logfile,l_verbose)
if(n_fl_msk)
  joinlines(tmplis//","//tmpmli, list2="", output="STDOUT", delim=" ",
    missing="Missing", maxchars=161, shortest+, verbose+, >> l_logfile)
else
  type(tmplis, map_cc+, device="terminal", >> l_logfile)

if(l_verbose) {
if(n_fl_msk)
  joinlines(tmplis//","//tmpmli, list2="", output="STDOUT", delim=" ",
    missing="Missing", maxchars=161, shortest+, verbose+)
else
  type(tmplis, map_cc+, device="terminal")
}
printlog("Alignment method: "//n_alignmethod,l_logfile,l_verbose)

# Checking images and files
l_fl_mef = no  # Default is not MEF
inimag = tmplis
n_files = 0
while(fscan(inimag,n_inim) != EOF) {
    n_files += 1
  if(n_files == 1) 
    n_ref = n_inim

  gimverify(n_inim)
  if (gimverify.status==1) {
    printlog("ERROR - IMCOADD: Input image missing: "//n_inim,l_logfile,yes)
    goto crash
  }
  if ((gimverify.status==0 && n_files>1 && !l_fl_mef) || \
      (gimverify.status>1 && n_files>1 && l_fl_mef) ) {
    printlog("ERROR - IMCOADD: Mix of MEF and other image formats",
      logfile,yes)
    goto crash
  }
  if (gimverify.status==0) {
    l_fl_mef=yes
    if(no==imaccess(n_inim//l_sci_ext) ) {
      printlog("ERROR - IMCOADD: Cannot access image : "//n_inim//"["//l_sci_ext//"]",
        l_logfile,yes)
      goto crash
    }
  } else {
    l_sci_ext="" ; l_var_ext="" ; l_dq_ext="" ; l_zero_ext=""
  }

# Exit if imtype!=fits and input are MEF
  if(l_fl_mef && n_imtype!=".fits") {
    printlog("ERROR - IMCOADD: Input is MEF, requires imtype=fits",l_logfile,yes)
    goto crash
  }

  # ---- MEF: If DQ available: fixpix, use DQ as immasks
   if(l_fl_mef && imaccess(n_inim//l_dq_ext)) {
     if(l_fl_fixpix || n_fl_med) {
       if(n_files==1)
         printlog("Making individual masks from DQ",l_logfile,l_verbose)

         if (l_dq_maskexpr == "default") {
             dq_inst = "INDEF"
             keypar(n_inim//l_zero_ext, n_key_inst, silent+)
             if (keypar.found) {
                 dq_inst = keypar.value
             }

             if (strstr("GMOS",dq_inst) == 1) {
                 dqexpr = "(((int(a) & 2) || (int(a) & 4)) > 0 ? 0 : a) > 0"
             } else {
                 dqexpr = "(int(a/2.)*2!=a)"
             }
         } else {
             dqexpr = l_dq_maskexpr
         }

       imexpr(dqexpr//" ? 1 : 0",tmpdq//str(n_files)//".pl",
        n_inim//l_dq_ext, dims="auto", intype="auto", outtype="auto",
        bwidth=0, btype="nearest", bpixval=0., rangecheck+, verbose-,
        exprdb="none")
     }
     if(l_fl_fixpix) {
       if(n_files==1)
         printlog("Fixpix input images using DQ",l_logfile,l_verbose)
       proto.fixpix(n_inim//l_sci_ext,tmpdq//str(n_files)//".pl",
           linterp="INDEF", cinterp="INDEF", verbose-, pixels-)
       gemhedit(n_inim//l_zero_ext,"IMCOFIX","done","Imcoadd: fixpix using DQ")
     }

     if(n_immasks=="DQ") {
        print(tmpdq//str(n_files), >> tmpmli)
        n_fl_msk=yes
     }
   }

  # ----

  if(n_files==1 && n_fl_find && access(n_inim//"_pos") && !n_fl_overwrite ) {
    printlog("ERROR - IMCOADD: Output file "//n_inim//"_pos exists",l_logfile,yes)
   goto crash
  } else if (n_files==1 && n_fl_find && access(n_inim//"_pos") && n_fl_overwrite ) 
    delete(n_inim//"_pos",verify-)

  if(n_files>1 && n_fl_map && access(n_inim//"_pos") && !n_fl_overwrite ) {
    printlog("ERROR - IMCOADD: Output file "//n_inim//"_pos exists",l_logfile,yes)
   goto crash
  } else if (n_files>1 && n_fl_map && access(n_inim//"_pos") && n_fl_overwrite ) 
    delete(n_inim//"_pos",verify-)

  if(n_fl_map && (access(n_inim//"_cen") || access(n_inim//"_trn")) && !n_fl_overwrite ) {
    printlog("ERROR - IMCOADD: Output file(s) "//n_inim//"_cen and/or "//n_inim//"_trn exist",l_logfile,yes)
   goto crash
  } else if (n_fl_map && (access(n_inim//"_cen") || access(n_inim//"_trn") ) && n_fl_overwrite) 
    delete(n_inim//"_cen,"//n_inim//"_trn",verify-, >>& "dev$null")

  if(n_fl_trn && imaccess(n_inim//"_trn") && !n_fl_overwrite ) {
    printlog("ERROR - IMCOADD: Output image "//n_inim//"_trn exists",l_logfile,yes)
    goto crash
  } else if (n_fl_trn && imaccess(n_inim//"_trn") && n_fl_overwrite )
    imdelete(n_inim//"_trn",verify-)

  # Check ron and gain in image
  imgets(n_inim//l_zero_ext,n_key_ron, >>& "dev$null")
  if (imgets.value == "0") {
   printlog("WARNING - IMCOADD: Adopting default read noise of "//l_ron//" e-",
   l_logfile,yes)
   n_ron[n_files]=l_ron
  } else {
   n_ron[n_files]=real(imgets.value)
  }
  roneff=roneff+n_ron[n_files]**2

  imgets(n_inim//l_zero_ext,n_key_gain, >>& "dev$null")
  if (imgets.value == "0") {
   printlog("WARNING - IMCOADD: Adopting default gain of "//l_gain//" e-/ADU",
   l_logfile,yes)
   n_gain[n_files]=l_gain
  } else {
   n_gain[n_files]=real(imgets.value)
  }
  gaineff=gaineff+n_gain[n_files]
  #print(n_files,n_gain[n_files],n_ron[n_files],gaineff,roneff)

}
inimag = ""

# check for existing badpix files and existing magnitude output files
if(n_fl_med) {
  if(access(n_ref//"_mag.tab") && !n_fl_overwrite) {
     printlog("ERROR - IMCOADD: Output file "//n_ref//"_mag.tab exists",l_logfile,yes)
     goto crash
  } else if (access(n_ref//"_mag.tab") && n_fl_overwrite)
     delete(n_ref//"_mag.tab",verify-)
     
  if(access(n_ref//"_mag.fits") && !n_fl_overwrite) {
     printlog("ERROR - IMCOADD: Output file "//n_ref//"_mag.fits exists",l_logfile,yes)
     goto crash
  } else if (access(n_ref//"_mag.fits") && n_fl_overwrite)
     delete(n_ref//"_mag.fits",verify-) 

  inimag=tmplis
  while(fscan(inimag,n_inim)!=EOF) {
     if( imaccess(n_inim//"badpix.pl") && !n_fl_overwrite) {
       printlog("ERROR - IMCOADD: Bad pixel image "//n_inim//"badpix exists",l_logfile,yes)
       goto crash
     } else if ( imaccess(n_inim//"badpix.pl") && n_fl_overwrite)
       imdelete(n_inim//"badpix.pl",verify-)
     if ( imaccess(n_inim//"badpix") ) {
       printlog("WARNING - IMCOADD: Deleting left over temporary image "//n_inim//"badpix",
         l_logfile,yes)
       imdelete(n_inim//"badpix",verify-)
     }

     if (access(n_inim//"_trn_mag") && !n_fl_overwrite ) {
       printlog("ERROR - IMCOADD: Output file "//n_inim//"_trn_mag exists",l_logfile,yes)
       goto crash
     } else if ( access(n_inim//"_trn_mag") && n_fl_overwrite )
       delete(n_inim//"_trn_mag",verify-)

  }
  inimag=""
}

# Define and check the access to the badpix file - first get the instrument
n_inst=""
n_camera=""

keypar(n_ref//l_zero_ext, n_key_inst, silent+)
if (keypar.found) {
    n_inst = keypar.value
    if (n_inst == "NIRI") {
        keypar(n_ref//l_zero_ext, n_key_camera, silent+)
        if (keypar.found)
            n_camera = keypar.value
        else {
            n_camera = "f6"
            printlog("WARNING - IMCOADD: No camera definition in header, \
                assuming f6", l_logfile, verbose+)
        }
    } else if ((n_inst == "F2") || (n_inst =="Flam")) {
        n_inst = "F2"
    } else if ((n_inst!="GMOS") && (n_inst!="Hokupaa+QUIRC") && \
        (n_inst!="GMOS-N") && (n_inst!="GMOS-S") && (n_inst!="F2")) {
        n_inst=""
        n_camera=""
        printlog("WARNING - IMCOADD: Instrument definition in header not \
            recognized", l_logfile, verbose+)
    }
} else {
    printlog("WARNING - IMCOADD: No instrument definition in header", \
        l_logfile, verbose+)
}

n_inport=0
if(n_inst!="") {
  keypar(n_ref//l_zero_ext,n_key_inport,silent+)
  if(keypar.found)
    n_inport=int(keypar.value)
  else
    printlog("WARNING - IMCOADD: No ISS port definition in header",l_logfile,yes)
}

printlog("Instrument + camera      : "//n_inst//"  "//n_camera,l_logfile,l_verbose)
printlog("Telescope ISS port number: "//n_inport,l_logfile,l_verbose)

# Define where the WCS is located
if (n_inst == "F2") {
    wcs_header = l_sci_ext
} else {
    wcs_header = l_zero_ext
}

# Get datamax if from header
if(n_key_sat!="" && n_key_sat!=" ") {
  keypar(n_ref//l_zero_ext,n_key_sat,silent+)
  if(keypar.found)
    n_datamax=real(keypar.value)
}

n_badim=no
if(n_badpixfile=="default" && n_inst=="") 
  n_badpixfile=""

if(n_badpixfile!="" && n_badpixfile!=" ") {
  if(n_badpixfile=="default") {
     if(n_inst=="Hokupaa+QUIRC")
       n_badpixfile="quirc$data/quirc_bpm.pl"
     if(n_inst=="NIRI")
       n_badpixfile="niri$data/niri_bpm.pl"
     if(n_inst=="GMOS" || n_inst=="GMOS-N" || n_inst=="GMOS-S") {
       xbin=1
       keypar(n_ref//l_sci_ext,"CCDSUM",silent+)
       if(keypar.found)
          print(keypar.value) | scan(xbin)

       n_badpixfile = "gmos$data/"//strlwr(n_inst)

       if (n_inst == "GMOS") {
           n_badpixfile = n_badpixfile//"-n"
       }

       keypar (n_ref//l_zero_ext, "DETTYPE", silent-)
       if (keypar.value == "SDSU II CCD") {
           dettype = "EEV"
       } else if (keypar.value == "SDSU II e2v DD CCD42-90") {
           dettype = "e2v"
       } else if (keypar.value == "S10892" || keypar.value == "S10892-N") {
           dettype = "HAM"
       }

       keypar (n_ref//l_zero_ext, "NAMPS", silent-)
       namps = int(keypar.value) * 3

       keypar (n_ref//l_zero_ext, "DETECTOR", silent-)
       if (str(keypar.value) == "GMOS + Blue1 + new CCD1") {
           # GMOS-S; newer EEV CCDs
           bpmversion = "v2"
       } else {
           bpmversion = "v1"
       }

       n_badpixfile = n_badpixfile//"_bpm_"//dettype//"_"//\
           str(xbin)//str(xbin)//"_"//str(namps)//"amp_"//bpmversion//"_"//\
           "mosaic.pl"
     }
  }
  if(no==access(n_badpixfile) && no==imaccess(n_badpixfile)) {
   printlog("WARNING - IMCOADD: Bad pixel file "//n_badpixfile//" not found",l_logfile,yes)
   n_badpixfile="gemtools$badpix.none"
  } else if(imaccess(n_badpixfile)) {
    n_badim=yes    # It's an image
  } else
    n_badim=no     # It's a file
} else
   n_badpixfile="gemtools$badpix.none"
printlog("Bad pixel file "//n_badpixfile,l_logfile,l_verbose)

printf("Read noise: %6.2f e-   Gain: %6.2f e-/ADU\n",n_ron[1],n_gain[1]) | scan(n_struct)
printlog(n_struct,l_logfile,l_verbose)

# Check access to median image
if(n_fl_med && imaccess(n_ref//"_med") && !n_fl_overwrite ) {
  printlog("ERROR - IMCOADD: Output image "//n_ref//"_med exists",l_logfile,yes)
  goto crash
} else if (n_fl_med && imaccess(n_ref//"_med") && n_fl_overwrite)
  imdelete(n_ref//"_med",verify-)

# Check access to cleaned average image
if(n_fl_add && imaccess(n_ref//"_add") && !n_fl_overwrite ) {
  printlog("ERROR - IMCOADD: Output image "//n_ref//"_add exists",l_logfile,yes)
  goto crash
} else if(n_fl_add && imaccess(n_ref//"_add") && n_fl_overwrite )
  imdelete(n_ref//"_add",verify-)

if(n_outimage!="") {
  gimverify(n_outimage)
  n_outimage=gimverify.outname
  if(n_fl_add && imaccess(n_outimage) && !n_fl_overwrite) {
    printlog("ERROR - IMCOADD: Output image "//n_outimage//" exists",l_logfile,yes)
    goto crash
  } else if(n_fl_add && imaccess(n_outimage) && n_fl_overwrite)
    imdelete(n_outimage,verify-)
}

if(n_fl_add) {
  if(n_outimage!="")
  printlog("Co-added cleaned output image: "//n_outimage,l_logfile,l_verbose)
  else
  printlog("Co-added cleaned output image: "//n_ref//"_add",l_logfile,l_verbose)
}

# Check access to uncleaned average image
if(n_fl_avg && imaccess(n_ref//"_avg") && !n_fl_overwrite ) {
  printlog("ERROR - IMCOADD: Output image "//n_ref//"_avg exists",l_logfile,yes)
  goto crash
} else if (n_fl_avg && imaccess(n_ref//"_avg") && n_fl_overwrite )
  imdelete(n_ref//"_avg",verify-)

# Check access to database if required
if(no==n_fl_map && no==access(n_database) && (n_fl_trn || n_fl_med) ) {
  printlog("ERROR - IMCOADD: Database file missing: "//n_database,l_logfile,yes)
  goto crash
}

# Check access to median image if required
if(no==n_fl_med && n_fl_add && no==imaccess(n_ref//"_med")) {
  printlog("ERROR - IMCOADD: Median image "//n_ref//"_med does not exist",l_logfile,yes)
  printlog("                 Use fl_med=yes to create the median image",l_logfile,yes)
  goto crash
}

# check alignsec if needed
if(n_alignsec!="default" && n_fl_map && n_alignmethod=="twodx") {
  print(n_alignsec) | tokens("STDIN",ignore+,begin_com="#",
    end_comm="eol",newlines-) | \
    match(":","STDIN",stop-,print_file-,meta+) | count("STDIN") | scan(n_sec)
  if(n_sec!=4) {
    printlog("ERROR - IMCOADD: alignsec must be given as ",l_logfile,yes)
    printlog("       [x11:x12,y11:y12] [x21:x22,y21:y22] or set to default",l_logfile,yes)
    goto crash
  }
}

# check if xoffset,yoffset related stuff is needed and set the parameters
if(n_alignmethod=="header" && n_fl_map) {
  n_key_pixscale=key_pixscale
  if(n_key_pixscale!="" && n_key_pixscale!=" ") {
    keypar(n_ref//l_zero_ext,n_key_pixscale,silent+)
    if(keypar.found)
      n_pixscale=real(keypar.value)
    else
      n_pixscale=pixscale
  } else
      n_pixscale=pixscale
printlog("Pixel scale: "//n_pixscale,l_logfile,l_verbose)

if (key_xoff == "default" || key_yoff == "default") {

    n_key_xoff = "XOFFSET"
    n_key_yoff = "YOFFSET"

    if (n_inst == "GMOS" || n_inst == "GMOS-N")  {

        n_instscale = 1.

        if(n_inport != 1) {
            n_xsign = 1
            n_ysign = 1
        } else {
            n_xsign = -1
            n_ysign = 1
        }

    } else if (n_inst == "GMOS-S") {

        n_instscale = 1.

        if (n_inport != 1) {
            # Good for port 3 after changing sign convention in seqexec
            n_xsign = -1
            n_ysign = -1
        } else {
            # Not tested
            n_xsign = -1
            n_ysign = 1
        }

    } else if (n_inst == "Hokupaa+QUIRC") {

        n_xsign = -1
        n_ysign = -1
        n_instscale = 1.611444

    } else if (n_inst == "NIRI") {

        n_key_xoff = "YOFFSET"
        n_key_yoff = "XOFFSET"
        n_instscale = 1.

        if (n_inport == 1 && n_camera == "f6") {
            n_xsign = -1
            n_ysign = 1
        }
        if (n_inport != 1 && n_camera == "f6") {
            n_xsign = 1
            n_ysign = 1 
        }
        if (n_inport == 1 && n_camera != "f6") {
            n_xsign = 1
            n_ysign = -1
        }
        if (n_inport != 1 && n_camera != "f6") {
            n_xsign = -1
            n_ysign = -1
        }

    } else if (n_inst == "F2") {

        n_instscale = 1.

        if (n_inport == 5) {
            n_xsign = -1
            n_ysign = -1
        } else {
            # Not tested
            n_xsign = -1
            n_ysign = -1
        }

    } else {
        # Take the input parameters
        n_instscale = instscale

        if (xsign == "default" || xsign == "positive")
            n_xsign = 1
        else
            n_xsign = -1

        if (ysign == "default" || ysign == "positive")
            n_ysign = 1
        else
            n_ysign = -1
    }
# ---------
} else {
  n_key_xoff=key_xoff ; n_key_yoff=key_yoff
    if(xsign=="default" || xsign=="positive")
      n_xsign=1
    else
      n_xsign=-1
    if(ysign=="default" || ysign=="positive")
      n_ysign=1
    else
      n_ysign=-1
    n_instscale=instscale
}

if(n_rotate) {
  printlog("WARNING - IMCOADD: Rotation >0.5 deg not handled based on ",l_logfile,verbose+)
  printlog("                   xoffset and yoffset from header. Setting rotate=no",l_logfile,verbose+)
  n_rotate=no
}

printlog("key_xoff ="//n_key_xoff//" key_yoff="//n_key_yoff,l_logfile,l_verbose)
printlog("xsign ="//n_xsign//"  ysign="//n_ysign//"  instscale="//str(n_instscale),
  l_logfile,l_verbose)

} # ---- end of offset parameter settings

if(n_alignmethod=="wcs" && n_fl_map) {
   keypar(n_ref//wcs_header,"CD1_1",silent+)
   if(no==keypar.found) {
      printlog("ERROR - IMCOADD: WCS missing",l_logfile,verbose+)
      goto crash
   }
   n_pixscale=real(keypar.value)
   keypar(n_ref//wcs_header,"CD1_2",silent+)
   if(no==keypar.found) {
      printlog("ERROR - IMCOADD: WCS missing",l_logfile,verbose+)
      goto crash
   }
   n_pixscale=sqrt(n_pixscale**2+real(keypar.value)*real(keypar.value))*3600.
   printf("%8.4f\n",n_pixscale) | scan(n_pixscale)
   printlog("Pixel scale: "//n_pixscale,l_logfile,l_verbose)
}


# cut one line off tmplis - this is the reference image
fields(tmplis,"1",lines="2-9999", quit_if_miss-, print_file_n-, > tmpout)
delete(tmplis,verify=no)
rename(tmpout,tmplis,field="root")

# size of image
imgets( n_ref//l_sci_ext, "i_naxis1" )
Xmax= int(imgets.value)
imgets( n_ref//l_sci_ext, "i_naxis2" )
Ymax= int(imgets.value)

# set statsec if default
if(n_statsec=="default") 
  n_statsec="[100:"//str(Xmax-100)//",100:"//str(Ymax-100)//"]"
printlog("Statistics section: "//n_statsec,l_logfile,l_verbose)

# Hardwire n_coords
n_coords=n_ref//"_pos"

# set fl_obj for a start
fl_obj=yes

# Find objects using daofind or check coordinate file existing
if(n_fl_find) {

# Get sky in order to estimate sigma of the sky - 
#  accurate value not required
imstat(n_ref//l_sci_ext//n_statsec,fields="midpt",
  lower=n_datamin,upper=n_datamax,nclip=0,lsigma=INDEF,usigma=INDEF,
  binwidth=0.1,format=no,cache=no) | scan(n_medsky)

printlog("Finding objects in "//n_ref,l_logfile,l_verbose)
     daofind(n_ref//l_sci_ext, output=n_ref//"_pos", starmap="", skymap="",
        boundary="nearest", constant=0., interactive=no, icommands="",
        gcommands="", wcsout="logical", cache=no, verify=no, update=")_.update",
        verbose=")_.verbose", graphics=")_.graphics", display=")_.display",
        datapars.scale=1., datapars.fwhmpsf=n_fwhm, datapars.emission=yes,
        datapars.sigma=(sqrt((n_ron[1]/n_gain[1])**2+n_medsky/n_gain[1])),
        datapars.datamin=n_datamin, datapars.datamax=n_datamax,
        datapars.noise="poisson", datapars.ccdread="", datapars.gain="",
        datapars.readnoise=0., datapars.epadu=1., datapars.exposure="",
        datapars.airmass="", datapars.filter="", datapars.obstime="",
        datapars.itime=1., datapars.xairmass=INDEF, datapars.ifilter="INDEF",
        datapars.otime="INDEF",
        findpars.threshold=n_thres, findpars.nsigma=2., findpars.ratio=1.,
        findpars.theta=0., findpars.sharplo=0.2, findpars.sharphi=1.,
        findpars.roundlo=-1., findpars.roundhi=1., findpars.mkdetections-)
     n_coords = n_ref//"_pos"
}
else {
  if(n_fl_map) {
    if(no==access(n_coords)) {
      printlog("ERROR - IMCOADD: Coordinate file missing: "//n_coords,l_logfile,yes)
      goto crash
    }
  }
}

# Only some values work for sampling
n_xsamp  = 1 ; n_ysamp =1
if(Xmax>500)
  n_xsamp = 10
if(Ymax>500)
  n_ysamp = 10
if(Xmax>1000)
  n_xsamp = 20
if(Ymax>1000)
  n_ysamp = 20
l_nxblock=geonxblock ; l_nyblock=geonyblock
if(n_fl_trn) {
  printlog("Sampling for geotran  : "//n_xsamp//"  "//n_ysamp,l_logfile,l_verbose)
  printlog("Block size for geotran: "//l_nxblock//"  "//l_nyblock,l_logfile,l_verbose)
}

# start of mapping IF

if(n_fl_map) {
printlog("Centering objects in "//n_ref,l_logfile,l_verbose)
apphot.center( n_ref//l_sci_ext, coords=n_coords, output= n_ref//"_cen" , 
    plotfile="", interactive=no, radplots=no, icommands="", gcommands="",
    wcsin="logical", wcsout="logical", cache=no, verify=no, update=")_.update",
    verbose=")_.verbose", graphics=")_.graphics", display=")_.display",
    datapars.scale=1., datapars.fwhmpsf=2.5, datapars.emission=yes,
    datapars.sigma=INDEF, datapars.datamin=n_datamin,
    datapars.datamax=n_datamax, datapars.noise="poisson", datapars.ccdread="",
    datapars.gain="", datapars.readnoise=0., datapars.epadu=1.,
    datapars.exposure="", datapars.airmass="", datapars.filter="",
    datapars.obstime="", datapars.itime=1., datapars.xairmass=INDEF,
    datapars.ifilter="INDEF", datapars.otime="INDEF",
    centerpars.calgorithm="centroid", centerpars.cbox=n_box,
    centerpars.cthreshold=0., centerpars.minsnratio=1., centerpars.cmaxiter=10,
    centerpars.maxshift=n_box, centerpars.clean=no, centerpars.rclean=1.,
    centerpars.rclip=2., centerpars.kclean=3., centerpars.mkcenter=no)

if(no==access(n_ref//"_cen")) {
 n_test="0"
} else {
# Select only objects where first centering is OK
  pdump(n_ref//"_cen","xcenter,ycenter,cerror,id","cerror=='NoError'",
    headers=no,parameters=no) | match("No","STDIN",stop=no, > tmpref)

# --------------- changing settings to catch few objects -------------
# Settings to catch too few objects and/or n_iter=0
# Do not try to iterate if nobj<2, use n_iter=1 and no rejection
  n_test="10"
  count(tmpref) | scan(n_test)
}

# If no objects, the only working setting is
# n_iter=0, n_fl_scale=no, n_alignmethod="header"
fl_obj=yes  # there are objects
if(int(n_test)==0)  {
   printlog("WARNING - IMCOADD: No objects sucessfully centered",
    l_logfile,verbose+)
   if(n_alignmethod!="header" && n_alignmethod!="wcs") {
     printlog("ERROR - IMCOADD: Only alignmethod=header or wcs may work",
       l_logfile,verbose+)
     if(n_immasks=="DQ" && n_fl_msk) 
        imdelete("@"//tmpmli,verify=no, >>& "dev$null")
     goto crash
   } else {
     printlog("WARNING - IMCOADD: Setting niter=0 and fl_scale=no",
      l_logfile,verbose+)
     n_iter=0 ; n_fl_scale=no ; fl_obj=no
   }
}

if(int(n_test)<3 && n_iter>1) {
  printlog("WARNING - IMCOADD: Too few objects. Setting niter=1, sigfit=100",
     l_logfile,verbose+)
  n_iter=1
  n_sigfit=100
}
if(int(n_test)<6 && n_iter>2) {
  printlog("WARNING - IMCOADD: Too few objects. Setting niter=2",
     l_logfile,verbose+)
  n_iter=2
}

# Check the fit geometry depending on the number of objects, be safe
# shift|xyscale|rotate|rscale|rxyscale|general
if(int(n_test)==1) {
  if(n_geofitgeom!="shift")
    printlog("WARNING - IMCOADD: Too few objects. Setting geofitgeom=shift",
      l_logfile,verbose+)
  n_geofitgeom="shift"
  if(n_rotate) {
    printlog("WARNING - IMCOADD: Cannot determine rotation, setting rotate=no",
      l_logfile,verbose+)
    n_rotate=no
  }
}

if(int(n_test)==2 && n_geofitgeom!="shift" && n_geofitgeom!="rotate") {
    printlog("WARNING - IMCOADD: Too few objects. Setting geofitgeom=rotate",
      l_logfile,verbose+)
  n_geofitgeom="rotate"
}


if(int(n_test)<6 && n_geofitgeom!="shift" && n_geofitgeom!="rotate" \
 && n_geofitgeom!="xyscale") {
  if(n_geofitgeom!="rscale")
    printlog("WARNING - IMCOADD: Too few objects. Setting geofitgeom=rscale",
      l_logfile,verbose+)
  n_geofitgeom="rscale"
}

if(n_iter==0)
  tmpres=n_database

printlog("Number of objects successfully centered in reference image "//n_test,
   l_logfile,l_verbose)

# --------------- end of changing settings ------------------

# set tmplas for the first round
tmplas=tmpref
# - - - - - - - 
printlog("Alignment method: "//n_alignmethod,l_logfile,l_verbose)

if(n_alignmethod=="user") {  # ---- start of user alignment for refim
print(" ==> Reference image: ", n_ref )
display( n_ref//l_sci_ext, 1, xmag=n_dispmag, ymag=n_dispmag )
if(fl_refmark) {
  print("Marking found objects in reference image")
fields(tmpref,fields="1,2,4",lines="1-",quit_if_miss=no,
  print_file_n=no) | tvmark(1,"STDIN",logfile="",autolog-,outimage="",
  deletions="",commands="",mark="point",color=202,pointsize=1,
  txsize=1,nxoffset=5,nyoffset=5,interactive=no,tolerance=1.5,
  label+,number-,font="raster")
}
if(fl_refmark)
  print("Objects found and centered are marked")

if (no == n_rotate && no == n_scale) {
  print("Point to one common object in reference image" )
  print("    strike any key")
  = fscan( imcur, x11, y11 )
} else {
  print("Point to first common object in reference image" )
  print("    strike any key")
  = fscan( imcur, x11, y11 )
  print("Point to second common object in reference image" )
  print("    strike any key")
  = fscan( imcur, x12, y12 )
}

} # ---- end of user alignment for reference image

if(n_alignmethod=="header") {  # ----- start of header alignment refim
   keypar(n_ref//l_zero_ext,n_key_xoff,silent+)
   if(keypar.found) {
     x11_ref=real(keypar.value)
     x11=Xmax/2
  } else {
    printlog("ERROR - IMCOADD: XOFFSET not defined for image "//n_ref,l_logfile,verbose+)
    if(n_immasks=="DQ" && n_fl_msk) 
       imdelete("@"//tmpmli,verify=no, >>& "dev$null")
    goto crash
  }
  keypar(n_ref//l_zero_ext,n_key_yoff,silent+)
  if(keypar.found) {
     y11_ref=real(keypar.value)
     y11=Ymax/2
  } else {
    printlog("ERROR - IMCOADD: YOFFSET not defined for image "//n_ref,l_logfile,verbose+)
    if(n_immasks=="DQ" && n_fl_msk) 
       imdelete("@"//tmpmli,verify=no, >>& "dev$null")
    goto crash
  }
} # ----- end of header alignment refim

if(n_alignmethod=="twodx") {  # ----- start of twodx alignment refim
 tmpimlas = n_ref
 # setup sections for xregister
 if(n_alignsec=="default") {
   n_alignsec="["//int(0.25*Xmax)+1//":"//int(0.375*Xmax)//","//int(0.25*Ymax)+1//":"//int(0.375*Ymax)//"]"
   n_alignsec=n_alignsec//" ["//int(0.625*Xmax)+1//":"//int(0.75*Xmax)//","//int(0.625*Ymax)+1//":"//int(0.75*Ymax)//"]"
   xc1 = (int(0.25*Xmax)+1+int(0.375*Xmax))/2.
   yc1 = (int(0.25*Ymax)+1+int(0.375*Ymax))/2.
   xc2 = (int(0.625*Xmax)+1+int(0.75*Xmax))/2.
   yc2 = (int(0.625*Ymax)+1+int(0.75*Ymax))/2.
   xsize1=int(0.375*Xmax)-(int(0.25*Xmax)+1)
   ysize1=int(0.375*Ymax)-(int(0.25*Ymax)+1)
   xsize2=int(0.75*Xmax)-(int(0.625*Xmax)+1)
   ysize2=int(0.75*Ymax)-(int(0.625*Ymax)+1)
 } else {
  for(n_ii=1;n_ii<=4;n_ii+=1) {
  print(n_alignsec) | tokens("STDIN",ignore+,begin_com="#",
    end_comm="eol",newlines-) | \
    match(":","STDIN",stop-,print_file-,meta+) | \
    fields("STDIN","1",lines=str(n_ii),qui-,print-) | scan(n_str)
    x11=int( substr(n_str,1,stridx(":",n_str)-1) )
    x12=int( substr(n_str,stridx(":",n_str)+1,strlen(n_str)) )
    if(n_ii==1) {
       xc1=(x11+x12)/2. ; xsize1=x12-x11 }
    if(n_ii==2) {
       yc1=(x11+x12)/2. ; ysize1=x12-x11 }
    if(n_ii==3) {
       xc2=(x11+x12)/2. ; xsize2=x12-x11 }
    if(n_ii==4) {
       yc2=(x11+x12)/2. ; ysize2=x12-x11 }
  }
 } # end of if-else
 if( n_xwindow > min(xsize1,ysize1,xsize2,ysize2) ) {
   n_xwindow=min(xsize1,ysize1,xsize2,ysize2)
   printlog("WARNING - IMCOADD: Sections for X-correlation too small",l_logfile,yes)
   printlog("                   Resetting xwindow = "//str(n_xwindow),l_logfile,yes)
 } 
} #  ---- end twodx alignment refim


# SCAN trough the other images
n=1
tmpimlas=n_ref  # this image is reference for the first image
inimag = tmplis
while(fscan(inimag,n_inim) != EOF) {
n=n+1

if(n_alignmethod=="user") {  # ---- start of user alignment

print(" ==> Image to be transformed: ", n_inim )
display( n_inim//l_sci_ext, 1, xmag=n_dispmag, ymag=n_dispmag )

if(no==n_rotate && no==n_scale) {
  print("Point to one common object in image to be transformed" )
  print ("    coordinates for last image: ", x11, y11)
  print ("    strike any key")
  = fscan( imcur, x21, y21 )
  printf("%7.1f %7.1f %7.1f %7.1f\n",x11, y11,x21, y21, > tmpin)
  geomap(tmpin,tmpres,xmin=1,xmax=Xmax,ymin=1,ymax=Ymax,
   transforms=n_inim//"_trn",results="",fitgeometry="shift",
   function="polynomial", xxorder=2, xyorder=2, xxterms="half", yxorder=2,
   yyorder=2, yxterms="half", maxiter=0, reject=3., calctype="real",
   verbose=l_verbose, interactive=no, graphics="stdgraph", cursor="")
# shift the variables
if(n_iter>0) {
  x11 = x21 ; y11 = y21
}
} else {
  print("Point to first common object in image to be transformed" )
  print ("    coordinates for last image: ", x11, y11)
  print ("    strike any key")
  = fscan( imcur, x21, y21 )
  print("Point to second common object in image to be transformed" )
  print ("    coordinates for last image: ", x12, y12)
  print ("    strike any key")
  = fscan( imcur, x22, y22 )
# transform ref coordinates to next image's coordinates
# allow rotation and/or scale change
  printf("%7.1f %7.1f %7.1f %7.1f\n",x11, y11,x21, y21, > tmpin)
  printf("%7.1f %7.1f %7.1f %7.1f\n",x12, y12,x22, y22, >> tmpin)
  if(n_rotate && !n_scale) {
    geomap(tmpin,tmpres,xmin=1,xmax=Xmax,ymin=1,ymax=Ymax,
     transforms=n_inim//"_trn",results="",fitgeometry="rotate",
     function="polynomial", xxorder=2, xyorder=2, xxterms="half", yxorder=2,
     yyorder=2, yxterms="half", maxiter=0, reject=3., calctype="real",
     verbose=l_verbose, interactive=no, graphics="stdgraph", cursor="")
  } else if (n_scale) {
    geomap(tmpin,tmpres,xmin=1,xmax=Xmax,ymin=1,ymax=Ymax,
     transforms=n_inim//"_trn",results="",fitgeometry="rscale",
     function="polynomial", xxorder=2, xyorder=2, xxterms="half", yxorder=2,
     yyorder=2, yxterms="half", maxiter=0, reject=3., calctype="real",
     verbose=l_verbose, interactive=no, graphics="stdgraph", cursor="")
  }
# shift the variables
if(n_iter>0) {
  x11 = x21 ; y11 = y21
  x12 = x22 ; y12 = y22
}

}

}  # ---- end of user alignment

if(n_alignmethod=="header") {  # ----- start of header alignment 
  keypar(n_inim//l_zero_ext,n_key_xoff,silent+)
  if(keypar.found) {
       x21=n_xsign*(real(keypar.value)-x11_ref)*n_instscale/n_pixscale+Xmax/2.
  } else {
    printlog("ERROR - IMCOADD: XOFFSET not defined for image "//n_inim,
     l_logfile,verbose+)
    goto crash
  }
  keypar(n_inim//l_zero_ext,n_key_yoff,silent+)
  if(keypar.found) {
       y21=n_ysign*(real(keypar.value)-y11_ref)*n_instscale/n_pixscale+Ymax/2.
  } else {
    printlog("ERROR - IMCOADD: YOFFSET not defined for image "//n_inim,
     l_logfile,verbose+)
    goto crash
  }
  printf("%7.1f %7.1f %7.1f %7.1f\n",x11, y11,x21, y21, > tmpin)
  geomap(tmpin,tmpres,xmin=1,xmax=Xmax,ymin=1,ymax=Ymax,
   transforms=n_inim//"_trn",results="",fitgeometry="shift",
   function="polynomial", xxorder=2, xyorder=2, xxterms="half", yxorder=2,
   yyorder=2, yxterms="half", maxiter=0, reject=3., calctype="real",
   verbose=l_verbose, interactive=no, graphics="stdgraph", cursor="")
# shift the variables if niter>0
if(n_iter>0) {
  x11 = x21 ; y11 = y21
}

} # ----- end of header alignment

if(n_alignmethod=="wcs") { # --------- start of wcs alignment
   wcsmap(n_inim//l_sci_ext, tmpimlas//l_sci_ext, tmpres,
    transforms=n_inim//"_trn", results="", xmin=1, xmax=Xmax, ymin=1,
    ymax=Ymax, nx=10, ny=10, wcs="world", transpose=no, xformat="%10.3f",
    yformat="%10.3f", wxformat="", wyformat="", fitgeometry="rscale",
    function="polynomial", xxorder=2, xyorder=2, xxterms="half", yxorder=2,
    yyorder=2, yxterms="half", reject=INDEF, calctype="real",
     verbose=l_verbose, interactive=no, graphics="stdgraph", gcommands="")
} # ----- end of wcs alignment

if(n_alignmethod=="twodx") { # --- start of twodx alignment
 printlog("xregister on the sections "//n_alignsec,l_logfile,l_verbose)
 printlog("Window size for search "//n_xwindow//" pixels",
   l_logfile,l_verbose)

 if(n_xwindow>21)
   n_correlation="fourier"
 else
   n_correlation="discrete"

 images.immatch.xregister(n_inim//l_sci_ext,tmpimlas//l_sci_ext,
   n_alignsec,tmpshift, output="", databasefmt+, append-, 
   records="",coords="",xlag=0,ylag=0,dxlag=0,
   dylag=0,background="none",border=INDEF,loreject=INDEF,hireject=INDEF,
   apodize=0.,filter="none",correlation=n_correlation,
   xwindow=n_xwindow,ywindow=n_xwindow,function="centroid",
   xcbox=n_box,ycbox=n_box,
   interp_type = "linear",boundary_typ = "nearest",constant = 0.,
   interactive=no,verbose=l_verbose,graphics = "stdgraph",display = "stdimage",
   gcommands = "",icommands = "") | \
    tee(l_logfile,out_type="text",append+)

 printlog("",l_logfile,no)
 printlog("Output from xregister:",l_logfile,no)
 type(tmpshift, >> l_logfile)
 printlog("",l_logfile,no)

 fields(tmpshift,"3,4",lines="7", quit-, print-, > tmpsh1)
 tcalc(tmpsh1,"c1","c1+"//str(xc1)) ; tcalc(tmpsh1,"c2","c2+"//str(yc1))
 tcalc(tmpsh1,"c3",str(xc1)) ; tcalc(tmpsh1,"c4",str(yc1))
 fields(tmpshift,"3,4",lines="8", quit-, print-, > tmpsh2)
 tcalc(tmpsh2,"c1","c1+"//str(xc2)) ; tcalc(tmpsh2,"c2","c2+"//str(yc2))
 tcalc(tmpsh2,"c3",str(xc2)) ; tcalc(tmpsh2,"c4",str(yc2))
 fields(tmpsh1//","//tmpsh2,"1-4",lines="1-999", quit-, print-, > tmpin)
 if(no==n_rotate && no==n_scale) {
  geomap(tmpin,tmpres,xmin=1,xmax=Xmax,ymin=1,ymax=Ymax,
   transforms=n_inim//"_trn",results="",fitgeometry="shift",
   function="polynomial", xxorder=2, xyorder=2, xxterm="half", yxorder=2,
   yyorder=2, yxterms="half", maxiter=0, reject=3., calctype="real",
   verbose=l_verbose,interactive=no, graphics="stdgraph", cursor="")
 }
 if(n_rotate && !n_scale) {
    geomap(tmpin,tmpres,xmin=1,xmax=Xmax,ymin=1,ymax=Ymax,
     transforms=n_inim//"_trn",results="",fitgeometry="rotate",
     function="polynomial",xxorder=2, xyorder=2, xxterm="half", yxorder=2,
   yyorder=2, yxterms="half", maxiter=0,reject=3., calctype="real",
     verbose=l_verbose,interactive=no, graphics="stdgraph", cursor="")
 } 
 if (n_scale) {
    geomap(tmpin,tmpres,xmin=1,xmax=Xmax,ymin=1,ymax=Ymax,
     transforms=n_inim//"_trn",results="",fitgeometry="rscale",
     function="polynomial",xxorder=2, xyorder=2, xxterm="half", yxorder=2,
   yyorder=2, yxterms="half", maxiter=0,reject=3., calctype="real",
     verbose=l_verbose,interactive=no, graphics="stdgraph", cursor="")
 }
 delete(tmpshift//","//tmpsh1//","//tmpsh2,verify-)
 if(n_iter>0) {
   tmpimlas=n_inim  # this image is reference for the next
 }
}  # ---- end of twodx aligment

# - - - - - - - - - 
delete(tmpin,verify=no, >>& "dev$null")
if(fl_obj) {
  geoxytran(tmplas,n_inim//"_pos",tmpres,n_inim//"_trn",
    geometry="geometric", direction="forward",
    xref=INDEF,yref=INDEF,xmag=INDEF,ymag=INDEF, xrotation=INDEF,
    yrotation=INDEF,
    xout=INDEF,yout=INDEF,xshift=INDEF,yshift=INDEF,xcolumn=1,ycolumn=2,
    calctype="real",xformat="",yformat="",min_sigdigit=7)
  if(tmpres!=n_database)
    delete(tmpres,verify=no, >>& "dev$null")
# Use this image as the next reference
  tmplas=n_inim//"_pos"
  if(n_iter>0) {
    tmpimlas=n_inim  # this image is reference for the next, used by wcs
  }

  printlog("Centering objects in "//n_inim,l_logfile,l_verbose)
  printlog("",l_logfile,l_verbose)
  apphot.center( n_inim//l_sci_ext, coords= n_inim//"_pos",
    output=tmpcen, plotfile="", interactive=no, radplots=no, icommands="",
    gcommands="", wcsin="logical", wcsout="logical", cache=no, verify=no,
    update=")_.update", verbose=")_.verbose", graphics=")_.graphics",
    display=")_.display",
    datapars.scale=1., datapars.fwhmpsf=2.5, datapars.emission=yes,
    datapars.sigma=INDEF, datapars.datamin=n_datamin,
    datapars.datamax=n_datamax, datapars.noise="poisson",
    datapars.ccdread="", datapars.gain="", datapars.readnoise=0.,
    datapars.epadu=1., datapars.exposure="", datapars.airmass="",
    datapars.filter="", datapars.obstime="", datapars.itime=1.,
    datapars.xairmass=INDEF, datapars.ifilter="INDEF", datapars.otime="INDEF",
    centerpars.calgorithm="centroid", centerpars.cbox=2.*n_box,
    centerpars.cthreshold=0., centerpars.minsnratio=1.,
    centerpars.cmaxiter=10, centerpars.maxshift=n_box, centerpars.clean=no,
    centerpars.rclean=1., centerpars.rclip=2., centerpars.kclean=3.,
    centerpars.mkcenter=no)
  pdump(tmpcen,"xcenter,ycenter,id",yes,headers=no,param=no, > tmpcoo)
  apphot.center( n_inim//l_sci_ext, coords=tmpcoo, output= n_inim//"_cen",
    plotfile="", interactive=no, radplots=no, icommands="", gcommands="",
    wcsin="logical", wcsout="logical", cache=no, verify=no,
    update=")_.update", verbose=")_.verbose", graphics=")_.graphics",
    display=")_.display",
    datapars.scale=1., datapars.fwhmpsf=2.5, datapars.emission=yes,
    datapars.sigma=INDEF, datapars.datamin=n_datamin,
    datapars.datamax=n_datamax, datapars.noise="poisson",
    datapars.ccdread="", datapars.gain="", datapars.readnoise=0.,
    datapars.epadu=1., datapars.exposure="", datapars.airmass="",
    datapars.filter="", datapars.obstime="", datapars.itime=1.,
    datapars.xairmass=INDEF, datapars.ifilter="INDEF", datapars.otime="INDEF",
    centerpars.calgorithm="centroid", centerpars.cbox=n_box,
    centerpars.cthreshold=0.,centerpars.minsnratio=1.,
    centerpars.cmaxiter=10, centerpars.maxshift=n_box/1.5, centerpars.clean=no,
    centerpars.rclean=1., centerpars.rclip=2., centerpars.kclean=3.,
    centerpars.mkcenter=no)
  delete(tmpcen//","//tmpcoo,verify=no)
  
  pdump(n_inim//"_cen","xcenter,ycenter,cerror",yes,headers=no,
    param=no, > tmpnew)
  # print("Ignore this warning")
  joinlines   (tmpref, tmpnew, output="STDOUT", delim=" ", missing="Missing",
    maxchars=161, shortest+, verbose+) | \
    match  ("NoError?*NoError", "STDIN", stop-, print+, meta+,> tmpcoo)
  fields (tmpcoo, "1,2,5,6", lines="1-", quit-, print-) | \
    unique("STDIN", > n_inim//"_trn")

# Have to check that centering went ok
tinfo(n_inim//"_trn",ttout-)
if(tinfo.nrows==0 && n_iter>0) {
  printlog("ERROR - IMCOADD: Failed to center any objects",l_logfile,verbose+)
  printlog("                 Adjust datamax and box as needed. May have to use niter=0.",
     l_logfile,verbose+)
  printlog("                 Check "//n_inim//"_cen and "//n_inim//"_trn",
     l_logfile,verbose+)
  delete(tmpcoo//","//tmpnew,verify-, >>& "dev$null")
  goto crash

} else {
  if(tmpres!=n_database)
    delete(tmpres,verify=no, >>& "dev$null")
}

if(fl_mark) {
print("Marking objects used for transformation")
fields(tmpcoo,fields="5,6,4",lines="1-",quit_if_miss=no,
  print_file_n=no) | tvmark(1,"STDIN",logfile="",autolog-,outimage="",
  deletions="",commands="",mark="point",color=202,pointsize=1,
  txsize=1,nxoffset=5,nyoffset=5,interactive=no,tolerance=1.5,
  label+,number-,font="raster")
}

delete(tmpcoo//","//tmpnew,verify=no, >>& "dev$null")

}

endcenter:  # jump point for failed centering  (is in fact not used anymore)

printf("Xmax = %5d   Ymax = %d\n",Xmax,Ymax) | scan(n_struct)
printlog(n_struct,l_logfile,l_verbose)
if(n_iter>0) {
  printlog("Entering geomap",l_logfile,l_verbose)
  printlog("Fitting geometry "//n_geofitgeom,l_logfile,l_verbose)
  if(n_geofitgeom=="general")
    printlog("Fitting order "//n_order,l_logfile,l_verbose)
  printlog("Iterating a maximum of "//str(n_iter)//" times",l_logfile,l_verbose)
} else 
  printlog("No iterations done with geomap",l_logfile,l_verbose)

# Iterate for geomap - geomap cannot do this on its own...
n_interact=no ; n_break=no ; n_tmpbase=tmpin
if(n_order>2 && n_iter>=3) 
  n_ordoff=n_order-2   # first 2 iterations with 2nd order fits
else
  n_ordoff=0

for(n_count=1;n_count<=n_iter;n_count+=1) {
  if(n_count>=3 || n_break)
    n_ordoff=0   # return to correct order
  if(n_count==n_iter || n_break) {
    n_interact=n_fl_inter ; n_tmpbase=n_database ; n_break=yes
    if(n_fl_inter) {
    print(" ") ; print("Last iteration, number "//n_count)
    }
  }
if((n_break || n_count==n_iter) && n_interact) {
print("---------------------------------------------------")
print(" x   graph x residuals ")
print(" y   graph y residuals ")
print(" d   delete nearest point ")
print(" u   undelete nearest point ")
print(" f   make new fit ")
print(" g   map data points ")
print(" q   exit curve fitting ")
print("---------------------------------------------------")
}

geomap(n_inim//"_trn",n_tmpbase,xmin=1,xmax=Xmax,ymin=1,ymax= Ymax,
  transforms=n_inim//"_trn",results=tmpres,fitgeometry=n_geofitgeom,
  function="legendre",xxorder=(n_order-n_ordoff),
  xyorder=(n_order-n_ordoff),xxterms="half",yxorder=(n_order-n_ordoff),
  yyorder=(n_order-n_ordoff),yxterms="half", maxiter=1,reject=n_sigfit, 
  calctype="real",verbose=l_verbose,interactive=n_interact, graphics="stdgraph",
  cursor="")
if(n_break)
  goto endgeoiter
# use these 4 lines if debugging with interactive=yes for all fits
# head(tmpres,nlines=5) | tokens(ignore-) | \
#   fields("STDIN","1",lines="29") | scan(n_sigx)
# head(tmpres,nlines=5) | tokens(ignore-) | \
#   fields("STDIN","1",lines="30") | scan(n_sigy)
#
# Use match instead of head, since head command fails if a path is used as 
# part of the input file name; tokens splits at "\" so sigma ends up not on 
# lines 37 and 38
match("rms:",tmpres,stop-,print_file-,meta+) | \
    tokens(ignore-, begin_commen="#", end_comment="eol", newlines=yes) | \
    fields("STDIN","1",lines="8", quit-, print-) | scan(n_sigx) 
match("rms:",tmpres,stop-,print_file-,meta+) | \
    tokens(ignore-, begin_commen="#", end_comment="eol", newlines=yes) | \
    fields("STDIN","1",lines="9", quit-, print-) | scan(n_sigy) 
#head(tmpres,nlines=6) | tokens(ignore-, begin_commen="#", end_comment="eol",
#    newlines=yes) | \
#    fields("STDIN","1",lines="37", quit-, print-) | scan(n_sigx)
#head(tmpres,nlines=6) | tokens(ignore-, begin_commen="#", end_comment="eol",
#    newlines=yes) | \
#    fields("STDIN","1",lines="38", quit-, print-) | scan(n_sigy)
delete(tmpin,verify=no, >>& "dev$null")
# stop when required accuracy is reached
if(n_sigx<=n_coolimit && n_sigy<=n_coolimit)
  n_break=yes

tselect(tmpres,tmpin//".fits", \
  "abs(c7) <= "//str(n_sigx*n_sigfit)//" && abs(c8) <= "//str(n_sigy*n_sigfit))

tdump (tmpin//".fits", cdfile="", pfile="", datafile="STDOUT", columns="", \
    rows="-", >> tmpin)

delete (tmpin//".fits", verify-, >& "dev$null")

count(tmpin) | scan(n_lines)
if(n_lines>=25) {
  delete(n_inim//"_trn",verify=no, >>& "dev$null")
  fields(tmpin,"1-4",lines="1-"//str(n_lines), quit-, print-, > n_inim//"_trn")
} else
 n_break=yes

delete(tmpres//","//tmpin,verify-, >>& "dev$null")

}

# jump point for converged geomap
endgeoiter:
# delete last results file if there is such a thing
if(tmpres!=n_database)
  delete(tmpres,verify-, >>& "dev$null")

}  # End of SCAN through the other images
inimag = ""  # close file
#n_inim = ""  # close file
delete(tmpref,verify=no, >>& "dev$null")
}  # End of mapping IF

# start transformation IF
if(n_fl_trn) {

# Get sky for first image
imstat(n_ref//l_sci_ext//n_statsec,fields="midpt",
  lower=n_datamin,upper=n_datamax,nclip=0,lsigma=INDEF,usigma=INDEF,
  binwidth=0.1,format=no,cache=no) | scan(n_medsky)
# changed to get better sky estimate  27.11.95
imstat(n_ref//l_sci_ext//n_statsec,fields="midpt",
  lower=max( n_datamin,n_medsky-5.*sqrt((n_ron[1]/n_gain[1])**2+n_medsky/n_gain[1]) ),
  upper=max( 50.,n_medsky+5.*sqrt((n_ron[1]/n_gain[1])**2+n_medsky/n_gain[1]) ),
  nclip=0,lsigma=INDEF,usigma=INDEF,binwidth=0.1,
  format=no,cache=no ) | scan(n_medsky)
printlog("",l_logfile,l_verbose)
printf("Median sky level for %s: %8.1f\n",n_ref,n_medsky) | scan(n_struct)
printlog(n_struct,l_logfile,l_verbose)
imarith(n_ref//l_sci_ext,"-",str(n_medsky),n_ref//"_trn", title="",
    divzero=0., hparams="", pixtype="real", calctype="real", verbose-, noact-)
gemdate ()
gemhedit(n_ref//"_trn","GEM-TLM",gemdate.outdate,
   "UT Last modification with GEMINI")
gemhedit(n_ref//"_trn","IMCOADD",gemdate.outdate,"UT Time stamp for imcoadd")
gemhedit(n_ref//"_trn","MED_SKY",str(n_medsky),"Median sky level")

# SCAN trough the other images  - transform and subtract off sky
n=1
inimag = tmplis
while(fscan(inimag,n_inim) != EOF) {
n=n+1
printlog("Transforming "//n_inim//" to "//n_inim//"_trn",l_logfile,l_verbose)
# Check if input image is short !! 
  n_realinim = n_inim
  imgets(n_inim//l_sci_ext,"i_pixtype")
  if(imgets.value == "3") {
    chpixtype(n_inim//l_sci_ext,tmpim,"real",oldpixtype="all",verbose=no)
    n_realinim = tmpim
  }
  geotran( n_realinim//l_sci_ext, n_inim//"_trn", n_database, n_inim//"_trn", 
    geometry="geometric", xin=INDEF,yin=INDEF,xshift=INDEF,yshift=INDEF,
    xout=INDEF,yout=INDEF,xmag=INDEF,ymag=INDEF,xrotation=INDEF,
    yrotation=INDEF, xmin= 1, xmax= Xmax, ymin= 1, ymax= Ymax,
    xscale=1.,yscale=1., ncols=INDEF,nlines=INDEF,
    xsample=n_xsamp, ysample=n_ysamp, interpolant=n_geointer,
    boundary="nearest",constant=0.,fluxconserve=yes,
    nxblock=l_nxblock,nyblock=l_nyblock, verbose=l_verbose )
  imdelete(tmpim,verify=no, >>& "dev$null")

# sky level
  imstat(n_realinim//l_sci_ext//n_statsec,fields="midpt",
    lower=n_datamin,upper=n_datamax,nclip=0,lsigma=INDEF,usigma=INDEF,
    binwidth=0.1,format=no,cache=no) | scan(n_insky)
# changed to get better sky estimate. 27.11.95
  imstat(n_realinim//l_sci_ext//n_statsec,fields="midpt",
    lower=max( n_datamin,n_insky-5.*sqrt((n_ron[n]/n_gain[n])**2+n_insky/n_gain[n]) ),
    upper=max( 50.,n_insky+5.*sqrt((n_ron[n]/n_gain[n])**2+n_insky/n_gain[n]) ),
    nclip=0,lsigma=INDEF,usigma=INDEF,binwidth=0.1,
    format=no,cache=no ) | scan(n_insky)
#  imstat(n_inim//n_statsec,fields="midpt",format=no,
#    lower=n_datamin,upper=(max(50.,3.*n_insky)) ) | scan(n_insky)
  printf("Median sky level for %s: %8.1f\n",n_inim,n_insky) | scan(n_struct)
  printlog(n_struct,l_logfile,l_verbose)

  imarith(n_inim//"_trn"//n_imtype,"-",str(n_insky),n_inim//"_trn"//n_imtype,
      title="", divzero=0., hparams="", pixtype="real", calctype="real",
      verbose-, noact-)
  gemdate ()
  gemhedit(n_inim//"_trn","GEM-TLM",gemdate.outdate,
    "UT Last modification with GEMINI")
  gemhedit(n_inim//"_trn","IMCOADD",gemdate.outdate,"UT Time stamp for imcoadd")
  gemhedit(n_inim//"_trn","MED_SKY",str(n_insky),"Median sky level")

}    # end of SCAN trough the images
inimag = ""
}    # end of transformation IF

# combine to get median image IF
# first subtract backgrounds as midpt for each image
# the mean of the backgrounds is added to the median filtered image

if(n_fl_med) {

# make the badpixel mask for the first image, or copy it if already
# and image
if(no==n_badim) {
# Is the badpixfile for trimmed or untrimmed image
type(n_badpixfile) | scan(tmpcen,tmpcoo)
printlog("",l_logfile,l_verbose)
printlog("Bad pixel file: "//n_badpixfile//" is for "//tmpcoo//" images",l_logfile,l_verbose)
if(tmpcoo == 'untrimmed') {
 imgets(n_ref,"TRIMSEC")
 x0 = real( substr(imgets.value,2,(stridx(":",imgets.value)-1) ) )
 imgets.value=substr(imgets.value,(stridx(",",imgets.value)+1),strlen(imgets.value) )
 y0 = real( substr(imgets.value,1,(stridx(":",imgets.value)-1) ) )
 fields(n_badpixfile,"1-4",lines="1-9999", > tmpbad)
 tcalc(tmpbad,"c1","c1+1-"//str(x0) )
 tcalc(tmpbad,"c2","c2+1-"//str(x0) )
 tcalc(tmpbad,"c3","c3+1-"//str(y0) )
 tcalc(tmpbad,"c4","c4+1-"//str(y0) )
 n_badpixfile=tmpbad
}
ccdred.badpiximage(n_badpixfile,n_ref//l_sci_ext,n_ref//"badpix",goodvalue=0,
  badvalue=1)
delete(tmpbad,verify=no, >>& "dev$null")

} else 
imcopy(n_badpixfile,n_ref//"badpix.pl",verbose-)

# let tmpcen and tmpcoo have some resonable values again
tmpcen = mktemp("tmpcen")
tmpcoo = mktemp("tmpcoo")

# Get sky for first image from the header, 
#  write image to list of images to combine
imgets(n_ref//"_trn","med_sky") ; n_medsky=real(imgets.value)
printf("Median sky level for %s: %8.1f\n",n_ref,n_medsky) | scan(n_struct)
printlog(n_struct,l_logfile,l_verbose)
n_files = 1
print(n_ref//"_trn", > tmpout)

# put badpixel mask in first image
gemhedit(n_ref//"_trn","BPM",n_ref//"badpix.pl","Bad pixel mask")
# change pixel type to be able to transform this image 
chpixtype(n_ref//"badpix",n_ref//"badpix","real",
   oldpixtype="all",verbose=no)

# sky level for other images, transform their badpixel mask
inimag = tmplis
while(fscan(inimag,n_inim) != EOF) {

  imgets(n_inim//"_trn","med_sky") ; n_insky=real(imgets.value)
  printf("Median sky level for %s: %8.1f\n",n_inim,n_insky) | scan(n_struct)
  printlog(n_struct,l_logfile,l_verbose)
  n_medsky += n_insky 
  n_files += 1
  print(n_inim//"_trn", >> tmpout)   # list of images to combine

# make individual mask if required
  if(n_fl_msk) {
   n_fl_del=no
   fields(tmpmli,"1",lines=str(n_files), quit-, print-) | scan(tmpba2)
   if(tmpba2!="none" && imaccess(tmpba2) ) {
#     n_mask=tmpba2
     imarith(tmpba2,"+",n_ref//"badpix"//n_imtype,n_mask, title="",
        divzero=0., hparams="", pixtype="real", calctype="real",
        verbose-, noact-)
   } else if(tmpba2!="none" && access(tmpba2) ) { 
     print("Making individual mask for "//n_inim)
     type(tmpba2) | scan(tmpcen,tmpcoo)
     printlog("Mask file is for "//tmpcoo//" image",l_logfile,l_verbose)
     fields(tmpba2,"1-4",lines="1-9999", quit-, print-, > tmpbad)
     if(tmpcoo == 'untrimmed') {
       imgets(n_inim,"TRIMSEC")
       x0 = real( substr(imgets.value,2,(stridx(":",imgets.value)-1) ) )
       imgets.value=substr(imgets.value,(stridx(",",imgets.value)+1),strlen(imgets.value) )
       y0 = real( substr(imgets.value,1,(stridx(":",imgets.value)-1) ) )
       tcalc(tmpbad,"c1","c1+1-"//str(x0) )
       tcalc(tmpbad,"c2","c2+1-"//str(x0) )
       tcalc(tmpbad,"c3","c3+1-"//str(y0) )
       tcalc(tmpbad,"c4","c4+1-"//str(y0) )
     }
     badpiximage(tmpbad,n_inim//l_sci_ext,n_mask,goodvalue=0,badvalue=1)
     delete(tmpbad,verify-, >>& "dev$null")
     imarith(n_mask,"+",n_ref//"badpix"//n_imtype,n_mask,title="",
        divzero=0., hparams="", pixtype="real", calctype="real",
        verbose-, noact-)
   } else
     n_mask = n_ref//"badpix"

  } else
  n_mask = n_ref//"badpix"

# let tmpcen, tmpcoo, and tmpbad have some resonable values again
  tmpcen = mktemp("tmpcen")
  tmpcoo = mktemp("tmpcoo")
  tmpbad = mktemp("tmpbad")

# transform badpixelmask
if(n_mask!=n_ref//"badpix")  
  printlog("Transforming "//n_ref//"badpix +mask to "//n_inim//"badpix",l_logfile,l_verbose)
else
  printlog("Transforming "//n_ref//"badpix to "//n_inim//"badpix",l_logfile,l_verbose)
  geotran(n_mask, n_inim//"badpix", n_database, 
    n_inim//"_trn", geometry="geometric", xin=INDEF, yin=INDEF, xshift=INDEF,
    yshift=INDEF, xout=INDEF,yout=INDEF,xmag=INDEF,ymag=INDEF,
    xrotation=INDEF,yrotation=INDEF,xmin= 1, xmax= Xmax, ymin= 1, ymax= Ymax, 
    xscale=1.,yscale=1., ncols=INDEF,nlines=INDEF, xsample=n_xsamp,
    ysample=n_ysamp, interpolant="linear",
    boundary="constant", constant=1.,fluxconserve=yes,
    nxblock=l_nxblock,nyblock=l_nyblock, verbose=l_verbose )
#  imarith(n_inim//"badpix","+",0.75,n_inim//"badpix",pixtype="real")
#  chpixtype(n_inim//"badpix",n_inim//"badpix","int",
#    oldpixtype="all",verbose=no)
#  imcopy(n_inim//"badpix"//n_imtype,n_inim//"badpix.pl",verbose=no)
  imcalc(n_inim//"badpix"//n_imtype,n_inim//"badpix.pl","im1+0.75",
    pixtype="int", nullval=0., verbose=no)
  imdelete(n_inim//"badpix"//n_imtype,verify=no)

  gemhedit(n_inim//"_trn"//n_imtype,"BPM",n_inim//"badpix.pl",
    "Bad pixel mask")

if(n_mask!=n_ref//"badpix")  # delete mask if not n_refbadpix
  imdelete(n_mask,verify=no)

# let n_mask have tmp value again
n_mask = mktemp("tmpmsk")

}             # end of scan of inimag
inimag = ""
n_medsky = n_medsky/n_files

# make individual mask file for 'ref'
if(n_fl_msk) {
 fields(tmpmli,"1",lines="1", quit-, print-) | scan(tmpba2)
   if(tmpba2!="none" && imaccess(tmpba2) ) {
     n_mask=tmpba2
    imarith(n_ref//"badpix","+",n_mask,n_ref//"badpix", title="",
        divzero=0., hparams="", pixtype="", calctype="",
        verbose-, noact-)
} else if(tmpba2!="none" && access(tmpba2) ) {
  print("Making individual mask for "//n_ref)
  type(tmpba2) | scan(tmpcen,tmpcoo)
  printlog("Mask pixel file is for "//tmpcoo//" image",l_logfile,l_verbose)
  fields(tmpba2,"1-4",lines="1-9999", quit-, print-, > tmpbad)
  if(tmpcoo == 'untrimmed') {
    imgets(n_inim,"TRIMSEC")
    x0 = real( substr(imgets.value,2,(stridx(":",imgets.value)-1) ) )
    imgets.value=substr(imgets.value,(stridx(",",imgets.value)+1),strlen(imgets.value) )
    y0 = real( substr(imgets.value,1,(stridx(":",imgets.value)-1) ) )
    tcalc(tmpbad,"c1","c1+1-"//str(x0) )
    tcalc(tmpbad,"c2","c2+1-"//str(x0) )
    tcalc(tmpbad,"c3","c3+1-"//str(y0) )
    tcalc(tmpbad,"c4","c4+1-"//str(y0) )
  }
  badpiximage(tmpbad,n_inim,n_mask,goodvalue=0,badvalue=1)
  imarith(n_ref//"badpix","+",n_mask,n_ref//"badpix", title="", divzero=0.,
    hparams="", pixtype="", calctype="", verbose-, noact-)
  imdelete(n_mask,verify=no)
  delete(tmpbad,verify=no)
 }
}
# let tmpcen, tmpcoo and tmpbad have some resonable values again
tmpcen = mktemp("tmpcen")
tmpcoo = mktemp("tmpcoo")//".fits"
tmpbad = mktemp("tmpbad")
n_mask = mktemp("tmpmsk")
imcopy(n_ref//"badpix"//n_imtype,n_ref//"badpix.pl",verbose=no)
imdelete(n_ref//"badpix"//n_imtype,verify=no)

#------------------------------------------------------
# check that there are objects in all *_trn
inimag=tmplis   # no suffix, ref image not included
while(fscan(inimag,n_inim)!=EOF) {
  tinfo(n_inim//"_trn",ttout-)
  if(tinfo.nrows==0)
    fl_obj=no
}
inimag=""

if(fl_obj) {  # ---- start of getting intensities of images from objects
inimag=tmplis   # no suffix, ref image not included
n_files = 1
while(fscan(inimag,n_inim) != EOF) {

n_files += 1

# make tables of _trn files and join using y-coordinate for
# objects in ref image  ==> one table

print("x_n1 r f6.2", > tmpbad)
print("y_n1 r f6.2", >> tmpbad)
print("x_n"//str(n_files)//" r f6.2", >> tmpbad)
print("y_n"//str(n_files)//" r f6.2", >> tmpbad)

if(n_files==2)
tcreate(tmpcoo,tmpbad,n_inim//"_trn",uparfile="",nskip=0, nlines=0, nrows=0,
    hist=no, extrapar=5, tbltype="default", extracol=0)
else {
  tcreate(tmpcen//".fits",tmpbad,n_inim//"_trn",uparfile="",nskip=0, nlines=0,
    nrows=0, hist=no, extrapar=5, tbltype="default", extracol=0)
#
  tjoin(tmpcoo,tmpcen//".fits",tmpnew//".fits","y_n1","y_n1", extrarows="neither",
    tolerance="0.0", casesens=yes)
#
#  tmatch(tmpcoo,tmpcen,tmpnew,"x_n1,y_n1","x_n1,y_n1",
#    maxnorm=0.00000001)

  tdelete(tmpcoo//","//tmpcen//".fits", go=yes, verify=no, default=yes)
  tchcol(tmpnew//".fits","x_n1_1","x_n1","","",verbose=no)
  tchcol(tmpnew//".fits","x_n1_2","junk","","",verbose=no)
# this is the fix!
#  tproject(tmpnew,tmpcoo,"x*,y*",uniq=yes)
  tproject(tmpnew//".fits",tmpcoo,"x_n1,y_n1",uniq=yes)
  tdelete(tmpnew//".fits", go=yes, verify=no, default=yes)
}
delete(tmpbad,verify=no)

}
inimag=""   # close file

# get magnitudes and make one table with output
# coordinates are always as for ref image because we use
# the transformed sky subtracted images
tdump(tmpcoo,cdfile="STDOUT",pfile="STDOUT",datafile=tmpcen,
  columns="x_n1,y_n1",rows="-", pwidth=-1, >& "dev$null")

for(n_i=1 ; n_i<=n_files ; n_i+=1) {
# tmpout is list with suffix  _trn
fields(tmpout,"1",lines=str(n_i), quit-, print-) | scan(n_inim)
apphot.phot(n_inim//""//n_imtype, skyfile="", coords=tmpcen,
    output=n_inim//"_mag", plotfile="", interactive=no, radplots=no,
    icommands="", gcommands="", wcsin="logical",wcsout="logical",
    cache=no,verify=no, update=")_.update", verbose=")_.verbose",
    graphics=")_.graphics", display=")_.display",
    datapars.scale=1., datapars.fwhmpsf=2.5, datapars.emission=yes,
    datapars.sigma=INDEF, datapars.datamin=n_datamin, 
    datapars.datamax=n_datamax, datapars.noise="poisson", datapars.ccdread="",
    datapars.gain="", datapars.readnoise=n_ron[n_i], datapars.epadu=n_gain[n_i],
    datapars.exposure="", datapars.airmass="", datapars.filter="",
    datapars.obstime="", datapars.itime=1., datapars.xairmass=INDEF,
    datapars.ifilter="INDEF", datapars.otime="INDEF",
    centerpars.calgorithm="none", centerpars.cbox=5., centerpars.cthreshold=0.,
    centerpars.minsnratio=1., centerpars.cmaxiter=10, centerpars.maxshift=1.,
    centerpars.clean=no, centerpars.rclean=1., centerpars.rclip=2.,
    centerpars.kclean=3., centerpars.mkcenter=no,
    fitskypars.salgorithm="constant", fitskypars.annulus=10.,
    fitskypars.dannulus=10., fitskypars.skyvalue=0., fitskypars.smaxiter=10,
    fitskypars.sloclip=0., fitskypars.shiclip=0., fitskypars.snreject=50,
    fitskypars.sloreject=3., fitskypars.shireject=3., fitskypars.khist=3.,
    fitskypars.binsize=0.1, fitskypars.smooth=no, fitskypars.rgrow=0.,
    fitskypars.mksky=no,
    photpars.weighting="constant", photpars.apertures=str(n_aperture),
    photpars.zmag=22., photpars.mkapert=no)

print("x_n"//str(n_i)//" r f6.2", > tmpbad)
print("y_n"//str(n_i)//" r f6.2", >> tmpbad)
print("mag_n"//str(n_i)//" r f6.2", >> tmpbad)
print("perr_n"//str(n_i)//" ch*15 a", >> tmpbad)

if(n_i==1) {
  pdump(n_inim//"_mag","xcenter,ycenter,mag,perror",yes,headers-,param-) | \
   tcreate(tmpmag//".fits",tmpbad,"STDIN",uparfile="",nskip=0, nlines=0, nrows=0,
    hist=no, extrapar=5, tbltype="default", extracol=0)
  delete(tmpbad,verify=no)
} else {
  pdump(n_inim//"_mag","xcenter,ycenter,mag,perror",yes,headers-,param-) | \
   tcreate(tmpnew//".fits",tmpbad,"STDIN",uparfile="",nskip=0,nlines=0, nrows=0,
    hist=no, extrapar=5, tbltype="default", extracol=0)
  delete(tmpbad,verify=no)
  tmerge(tmpmag//".fits,"//tmpnew//".fits",tmpbad//".fits","merge", \
      allcols=yes, tbltype="default", \
      allrows=100, extracol=0)

  tdelete(tmpmag//".fits,"//tmpnew//".fits", go=yes, verify=no, default=yes)
  rename(tmpbad//".fits",tmpmag//".fits", field="all")
}

}
delete(tmpcen,verify=no)
tdelete(tmpcoo,verify=no)

# select stars for which everything is ok
for(n_i=1 ; n_i<=n_files ; n_i+=1) {
  tselect(tmpmag//".fits",tmpnew//".fits",
   "mag_n"//str(n_i)//"<9999. && perr_n"//str(n_i)//"=='NoError'")
  tdelete(tmpmag//".fits", go=yes, verify=no, default=yes)
  rename(tmpnew//".fits",tmpmag//".fits", field="all")
}

# check output table one last time in case all objects are bad
tinfo(tmpmag//".fits",ttout-, >>& "dev$null")
if(tinfo.nrows==0) {
  printlog("WARNING - IMCOADD: No good objects in common. Setting fl_scale=no",l_logfile,verbose+)
  delete(tmpmag//".fits",verify-, >>& "dev$null")
  n_fl_scale=no
  goto noscaling
}

printlog("",l_logfile,l_verbose)
# Calculate mean intensity
n_totint=0.
for(n_i=1 ; n_i<=n_files ; n_i+=1) {
  fields(tmpout,"1",lines=str(n_i), quit-, print-) | scan(n_inim)  # image name
  tcalc(tmpmag//".fits","tot_"//str(n_i),"10**(0.4*(22.-mag_n"//str(n_i)//"))" )
  tstat(tmpmag//".fits","tot_"//str(n_i),outtable="STDOUT", lowlim=INDEF, highlim=INDEF,
    rows="-", >& "dev$null")
  if(n_i==1)
    n_refint = tstat.mean
  n_inint = tstat.mean
  printf("Mean intensity for   %s:  (N,absolute,relative)= %4d %8.0f %6.3f\n",
     n_inim,tstat.nrows,n_inint,(n_inint/n_refint) ) | scan(n_struct)
  printlog(n_struct,l_logfile,l_verbose)
# put the relative intensity in the image header
# the images are *multiplied* by this factor, so it has to be the
# inverse relative intensity
  gemhedit(n_inim//""//n_imtype,"RELINT",(n_refint/n_inint),
    "Inverse of the relative intensity")
}
# save the magnitudes for control
tcopy (tmpmag//".fits",n_ref//"_mag.fits", verbose=no)
gemhedit (n_ref//"_mag.fits[0]", "FILENAME", n_ref//"_mag.fits", "", \
    delete-, upfile="")
delete (tmpmag//".fits", verify-, >& "dev$null")

} else { # ---- end of getting intensities of images from objects
# if no objects, turn off scaling
  if(n_fl_scale) {
     n_fl_scale=no
     printlog("WARNING - IMCOADD: No objects. Setting fl_scale=no",l_logfile,verbose+)
  }
}
#------------------------------------------------------
noscaling:

# threshold not used for combining!
# the median image is on the scale of the reference image
if(n_fl_scale==yes) {
  imcombine("@"//tmpout,n_ref//"_med",headers="",bpmasks="",
    rejmasks="",nrejmasks="",expmasks="",sigmas="",logfile=l_logfile,
    combine="median",reject="none",project=no,outtype="real",outlimits="",
    offsets="none",masktype="goodvalue",maskvalue=0.,blank=0.,
    scale="!RELINT",zero="none",weight="none",statsec="",expname="",
    lthreshold=INDEF,hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1, mclip=yes,
    lsigma=3., hsigma=3., rdnoise="0.", gain="1.", snoise="0.", sigscale=0.1,
    pclip=-0.5, grow=0.)
  gemdate ()
  gemhedit(n_ref//"_med","GEM-TLM",gemdate.outdate,
    "UT Last modification with GEMINI")
  gemhedit(n_ref//"_med","IMCOADD",gemdate.outdate,"UT Time stamp for imcoadd")
  hedit(n_ref//"_med","COMBSC","imcoadd combine scaled with RELINT",
    add=yes,addonly=no,delete=no,verify=no,show=no,update=yes)
} else {
  imcombine("@"//tmpout,n_ref//"_med",headers="",bpmasks="",
    rejmasks="",nrejmasks="",expmasks="",sigmas="",logfile=l_logfile,
    combine="median",reject="none",project=no,outtype="real",outlimits="",
    offsets="none",masktype="goodvalue",maskvalue=0.,blank=0.,
    scale="none",zero="none",weight="none",statsec="",expname="",
    lthreshold=INDEF,hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1, mclip=yes,
    lsigma=3., hsigma=3., rdnoise="0.", gain="1.", snoise="0.",
    sigscale=0.1, pclip=-0.5, grow=0.)
  gemdate ()
  gemhedit(n_ref//"_med","GEM-TLM",gemdate.outdate,
    "UT Last modification with GEMINI")
  gemhedit(n_ref//"_med","IMCOADD",gemdate.outdate,"UTTime stamp for imcoadd")
  hedit(n_ref//"_med","COMBSC","imcoadd combine not scaled with RELINT",
    add=yes,addonly=no,delete=no,verify=no,show=no,update=yes)
}

gemhedit(n_ref//"_med","med_sky",str(n_medsky),"Median sky level")
# remove reference to BPM for reference image
hedit(n_ref//"_med","BPM",add=no,addonly=no,delete=yes,verify=no,
  show=no,update=yes)

# add sky to median image if cleaned mean images not to be made
if(no==n_fl_add) {
imarith(n_ref//"_med","+",str(n_medsky),n_ref//"_med", title="", divzero=0.,
    hparams="", pixtype="real", calctype="real", verbose-, noact-)
}

# Update gain and readnoise keywords
imgets(n_ref//"_med","GAINORIG", >& "dev$null")
if (imgets.value=="0") {
   gemhedit(n_ref//"_med","GAINORIG",n_gain[1],"Input gain")
}
imgets(n_ref//"_med","RONORIG", >& "dev$null")
if (imgets.value=="0") {
   gemhedit(n_ref//"_med","RONORIG",n_ron[1],"Input read-noise")
}
# Effective ron and gain for combined image
# ron in e-/pix  as expected by other iraf tasks
if(n_key_gain=="" || n_key_gain==" ")
  n_key_gain="GAIN"
gemhedit(n_ref//"_med",n_key_gain,(2.*gaineff/3.),"Effective gain [e-/ADU]")

if(n_key_ron=="" || n_key_ron==" ")
  n_key_ron="RDNOISE"
gemhedit(n_ref//"_med",n_key_ron,(sqrt(2.*roneff/3.)),"Effective read-noise [e-]")

}  # end of fl_med IF

# check if n_fl_add=yes and n_fl_med=no
# if so assume that sky has to be subtracted off median image
if(n_fl_add && !n_fl_med) {
 print("Subtracting sky level off median image")
 imgets(n_ref//"_med","med_sky")
 imarith(n_ref//"_med","-",real(imgets.value),n_ref//"_med", title="",
    divzero=0., hparams="", pixtype="real", calctype="real", verbose-, noact-)

} 

# make the list of images if not made during fl_med
# get the median sky -- not safe to get it from median image
# as it may not exist if fl_med+ was never run
if(no==n_fl_med && (n_fl_add || n_fl_avg)) {
  print(n_ref//"_trn", > tmpout)
  imgets(n_ref//"_trn","med_sky")
  n_medsky=real(imgets.value) ; n_files=1
  inimag = tmplis
    while(fscan(inimag,n_inim) != EOF) {
      print(n_inim//"_trn", >> tmpout)   # list of images to combine
      imgets(n_inim//"_trn","med_sky")
      n_medsky=n_medsky+real(imgets.value) ; n_files+=1
    }
  inimag=""
n_medsky=n_medsky/n_files
print("Median sky: "//str(n_medsky))
}


# make masks for individual images by comparison to median image
# then calculate mean of all images
if(n_fl_add) {

printlog("",l_logfile,l_verbose)

n=0
inimag = tmpout   # suffix _trn
while(fscan(inimag,n_inim) != EOF) {
  n=n+1

  n_ron[n] = (n_ron[n]/n_gain[n])**2          # RON^2 in ADU
  n_scnoise = n_scnoise**2           # Noise scaling ^2  (fraction)
  # make n_lim absolute if key_limit=yes
  if(key_limit)
    n_lim=limit*sqrt(n_ron[n]+n_medsky/n_gain[n])

  printlog("Masking cosmic ray events in "//n_inim,l_logfile,l_verbose)
# difference image, get the scale, remember the header contains
# the inverse scale
  if(n_fl_scale) {
     imgets(n_inim,"RELINT") ; n_inint = 1./real(imgets.value)
  } else 
     n_inint = 1.

  if (imaccess(n_inim//"_dif") && n_fl_overwrite) {
        imdelete(n_inim//"_dif", verify-)
  }
  imcalc(n_inim//","//n_ref//"_med",n_inim//"_dif",
   "im1/"//str(n_inint)//"-im2",pixtype="real",
   nullval=0.,verbose=no)
# sigma image, get sky value
  imgets(n_inim,"MED_SKY") ; n_insky = real(imgets.value)

print("sqrt(( 1./("//str(n_inint)//"**2)+\
1.5625/"//str(n_files)//")*"//str(n_ron[n])//" + \
( (im1+"//str(n_insky)//")/("//str(n_inint)//"**2) + \
(im2+"//str(n_medsky)//")*1.5625/"//str(n_files)//" ) /"//str(n_gain[n])//\
" + im2**2*"//n_scnoise//" )", 
> tmpbad)

  imcalc(n_inim//","//n_ref//"_med",n_sigim,"@"//tmpbad,pixtype="real",
    nullval=0., verbose=no)
  delete(tmpbad,verify=no)

# mask image  one: deviation, positive only
#            zero: no deviation
print("(im3<="//n_lim//" && (im1 >"//n_lsig//"*im2 || im1 >"//n_llim//"))",
 > tmpbad)
# || \
# (im3>"//n_lim//" && (im1 >"//n_usig//"*im2 || im1 >"//n_ulim//"))",
#  > tmpbad )

  imcalc(n_inim//"_dif,"//n_sigim//","//n_ref//"_med",n_mask//".pl",
    "@"//tmpbad,pixtype="int", nullval=0., verbose=no)

delete(tmpbad,verify=no)
# if growthrad=1 shift n_mask around and add
if(n_growthrad>0) {
  printlog("Applying growth radius "//n_growthrad,l_logfile,l_verbose)

# Code changed Aug 30, 2001 to solve bus errors for large images
# Make mask 1 and 0 only
  imcalc(n_mask//".pl",n_m4//".pl","if im1>0 then 1 else 0",pixtype="old",
    nullval=0., verbose-)
# Make a mask with huge values in badpix
  imarith(n_m4//".pl","*","1000.",n_m0//".pl", title="", divzero=0.,
    hparams="", pixtype="", calctype="", verbose-, noact-)
# Combine and use output plfile as the growth
  imcombine(n_mask//".pl,"//n_m0//".pl",n_m1//".pl",headers="",
    bpmasks="",rejmasks="",nrejmasks=n_m2//".pl",expmasks="",sigmas="",
    logfile="",combine="average",reject="minmax",project=no,outtype="real",
    outlimits="",offsets="none",masktype="none",maskvalue=0.,blank=1.,
    scale="none",zero="none",weight="none",statsec="",expname="", 
    lthreshold=INDEF,hthreshold=INDEF,nlow=0,nhigh=1,nkeep=1,mclip=yes,
    lsigma=3., hsigma=3., rdnoise="0.", gain="1.", snoise="0.", sigscale=0.1,
    pclip=-0.5, grow=1.)           
   
# Add the original mask and the growth   
  imcalc(n_mask//".pl,"//n_m2//".pl",n_m3//".pl",
      "if im1>0 || im2>1 then 1 else 0", pixtype="old", nullval=0.0, verbose-)
# Clean up
  imdelete(n_mask//".pl,"//n_m0//".pl,"//n_m1//".pl,"//n_m2//".pl,"//n_m4//".pl",verify-)
  imrename(n_m3//".pl",n_mask//".pl",verbose-)
}

# add the mask image to the badpixmask for the image
n_inim = substr( n_inim,1,(strlen(n_inim)-4) )
   imarith(n_inim//"badpix.pl","+",n_mask//".pl",n_inim//"badpix.pl",
     title="", divzero=0., hparam="", pixtype="", calctype="", verbose-, noact-)
# delete  mask and sigma and difference
   imdelete(n_mask//".pl,"//n_sigim//","//n_inim//"_trn_dif",verify=no)

}
inimag = ""  # close file, end of masking

# combine the images with masking

# threshold not used for combining, scale if fl_scale=yes
if(n_fl_scale) {
  imcombine("@"//tmpout,n_ref//"_add",headers="",bpmasks="",rejmasks="",
    nrejmasks="",expmasks="",sigmas="",logfile=l_logfile,combine="average",
    reject="none",project=no,outtype="real",outlimits="",offsets="none",
    masktype="goodvalue",maskvalue=0.,blank=0.,scale="!RELINT",
    zero="none",weight="none",statsec=n_statsec,expname="",
    lthreshold=INDEF,hthreshold=INDEF)
} else {
  imcombine("@"//tmpout,n_ref//"_add",headers="",bpmasks="",rejmasks="",
    nrejmasks="",expmasks="",sigmas="",logfile=l_logfile,combine="average",
    reject="none",project=no,outtype="real",outlimits="",offsets="none",
    masktype="goodvalue",maskvalue=0.,blank=0.,scale="none",
    zero="none",weight="none",statsec=n_statsec,expname="",
    lthreshold=INDEF,hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1, mclip=yes,
    lsigma=3., hsigma=3., rdnoise="0.", gain="1.", snoise="0.", sigscale=0.1,
    pclip=-0.5, grow=0.)
}

# Add median sky to mean image
imarith(n_ref//"_add","+",str(n_medsky),n_ref//"_add", title="", divzero=0.,
    hparams="", pixtype="real", calctype="real", verbose-, noact-)

# Pack up as MEF if input is MEF
if(l_fl_mef) {
 tmpimlas=mktemp("tmpimlas")
 imrename(n_ref//"_add",tmpimlas,verbose-)
 wmef(tmpimlas,n_ref//"_add",extnames=n_sci_ext,phu=n_ref,verbose-)
 imdelete(tmpimlas,verify-)
}

gemdate ()
gemhedit(n_ref//"_add"//l_zero_ext,"GEM-TLM",gemdate.outdate,
  "UT Last modification with GEMINI")
gemhedit(n_ref//"_add"//l_zero_ext,"IMCOADD",gemdate.outdate,"UT Time stamp for imcoadd")
gemhedit(n_ref//"_add"//l_zero_ext,"IMCOBPM",n_badpixfile,"Detector bad pixel file")

if(n_fl_scale) {
hedit(n_ref//"_add"//l_zero_ext,"COMBSC","imcoadd combine scaled with RELINT",
  add=yes,addonly=no,delete=no,verify=no,show=no,update=yes)
} else {
hedit(n_ref//"_add"//l_zero_ext,"COMBSC","imcoadd combine not scaled with RELINT",
  add=yes,addonly=no,delete=no,verify=no,show=no,update=yes)

}

# Add median sky  to median image 
imarith(n_ref//"_med","+",str(n_medsky),n_ref//"_med", title="", divzero=0.,
    hparams="", pixtype="real", calctype="real", verbose-, noact-)

gemhedit(n_ref//"_add"//l_zero_ext,"med_sky",str(n_medsky),"Median sky level")

# various header info in n_ref//"_add"

# Input ron and gain
# Update gain and readnoise keywords
imgets(n_ref//"_add"//l_zero_ext,"GAINORIG", >& "dev$null")
if (imgets.value=="0") {
   gemhedit(n_ref//"_add"//l_zero_ext,"GAINORIG",n_gain[1],"Input gain")
   gemhedit(n_ref//"_add"//l_sci_ext,"GAINORIG",n_gain[1],"Input gain")
}
imgets(n_ref//"_add"//l_zero_ext,"RONORIG", >& "dev$null")
if (imgets.value=="0") {
   gemhedit(n_ref//"_add"//l_zero_ext,"RONORIG",(sqrt(n_ron[1])*n_gain[1]),"Input read-noise")
   gemhedit(n_ref//"_add"//l_sci_ext,"RONORIG",(sqrt(n_ron[1])*n_gain[1]),"Input read-noise")
}

# Effective ron and gain for combined image
# ron in e-/pix  as expected by other iraf tasks
if(n_key_gain=="" || n_key_gain==" ")
  n_key_gain="GAIN"
gemhedit(n_ref//"_add"//l_zero_ext,n_key_gain,gaineff,"Effective gain [e-/ADU]")
gemhedit(n_ref//"_add"//l_sci_ext,n_key_gain,gaineff,"Effective gain [e-/ADU]")

if(n_key_ron=="" || n_key_ron==" ")
  n_key_ron="RDNOISE"
gemhedit(n_ref//"_add"//l_zero_ext,n_key_ron,(sqrt(roneff)),"Effective read-noise [e-]")
gemhedit(n_ref//"_add"//l_sci_ext,n_key_ron,(sqrt(roneff)),"Effective read-noise [e-]")

gemhedit(n_ref//"_add"//l_zero_ext,"SCNOISE",str(sqrt(n_scnoise)),
  "imcoadd scnoise")
gemhedit(n_ref//"_add"//l_zero_ext,"LIMIT",str(n_lim),"imcoadd limit for cleaning")
gemhedit(n_ref//"_add"//l_zero_ext,"LOWSIG",str(n_lsig),"imcoadd sigma rejection")
gemhedit(n_ref//"_add"//l_zero_ext,"LOWLIM",str(n_llim),
  "imcoadd absolute rejection limit")
gemhedit(n_ref//"_add"//l_zero_ext,"GROWTHR",str(n_growthrad),
  "imcoadd growth radius")
gemhedit(n_ref//"_add"//l_zero_ext,"GEOGEOM",n_geofitgeom,
  "imcoadd geomap geometry")
gemhedit(n_ref//"_add"//l_zero_ext,"GEOINTER",n_geointer,
  "imcoadd geotrans interpolation")
# remove reference to BPM for reference image
hedit(n_ref//"_add"//l_sci_ext,"BPM",add=no,addonly=no,delete=yes,
  verify=no,show=no,update=yes)

}  # end of fl_add

# uncleaned mean image - fl_avg=yes
if(n_fl_avg) {
if(n_fl_scale) {
imcombine("@"//tmpout,n_ref//"_avg",headers="",bpmasks="",rejmasks="",
  nrejmasks="",expmasks="",sigmas="",logfile=l_logfile,combine="average",
  reject="none",project=no,outtype="real",outlimits="",offsets="none",
  masktype="none",maskvalue=0.,blank=0.,scale="!RELINT",zero="none",
  weight="none",statsec=n_statsec,expname="",
  lthreshold=INDEF,hthreshold=INDEF, nlow=1, nhigh=1, mclip=yes, lsigma=3.,
  hsigma=3., rdnoise="0.", gain="1.", snoise="0.", sigscale=0.1, pclip=-0.5,
  grow=0.)
gemdate ()
gemhedit(n_ref//"_avg","GEM-TLM",gemdate.outdate,
   "UT Last modification with GEMINI")
gemhedit(n_ref//"_avg","IMCOADD",gemdate.outdate,"UT Time stamp for imcoadd")
hedit(n_ref//"_avg","COMBSC","imcoadd combine scaled with RELINT",add=yes,
  addonly=no,delete=no,verify=no,show=no,update=yes)
} else {
imcombine("@"//tmpout,n_ref//"_avg",headers="",bpmasks="",rejmasks="",
  nrejmasks="",expmasks="",sigmas="",logfile=l_logfile,combine="average",
  reject="none",project=no,outtype="real",outlimits="",offsets="none",
  masktype="none",maskvalue=0.,blank=0.,scale="none",zero="none",
  weight="none",statsec=n_statsec,expname="",
  lthreshold=INDEF,hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1, mclip=yes,
  lsigma=3., hsigma=3., rdnoise="0.", gain="1.", snoise="0.", sigscale=0.1,
  pclip=-0.5, grow=0.)
gemdate ()
gemhedit(n_ref//"_avg","GEM-TLM",gemdate.outdate,
  "Last modification with GEMINI")
gemhedit(n_ref//"_avg","IMCOADD",gemdate.outdate,"Time stamp for imcoadd")
hedit(n_ref//"_avg","COMBSC","imcoadd combine not scaled with RELINT",add=yes,
  addonly=no,delete=no,verify=no,show=no,update=yes)
}
imarith(n_ref//"_avg","+",str(n_medsky),n_ref//"_avg", title="", divzero=0.,
    hparams="", pixtype="real", calctype="real", verbose-, noact-)
gemhedit(n_ref//"_avg","MED_SKY",str(n_medsky),"Median sky level")
# remove reference to BPM for reference image
hedit(n_ref//"_avg","BPM",add=no,addonly=no,delete=yes,verify=no,show=no,
  update=yes)

# Input ron and gain
# Update gain and readnoise keywords
imgets(n_ref//"_avg","GAINORIG", >& "dev$null")
if (imgets.value=="0") {
   gemhedit(n_ref//"_avg","GAINORIG",n_gain[1],"Input gain")
}
imgets(n_ref//"_avg","RONORIG", >& "dev$null")
if (imgets.value=="0") {
   gemhedit(n_ref//"_avg","RONORIG",(sqrt(n_ron[1])*n_gain[1]),"Input read-noise")
}

# Effective ron and gain for combined image
# ron in e-/pix  as expected by other iraf tasks
# n_ron before this statement is RON^2 in ADU/pix
if(n_key_gain=="" || n_key_gain==" ")
  n_key_gain="GAIN"
gemhedit(n_ref//"_avg",n_key_gain,gaineff,"Effective gain [e-/ADU]")

if(n_key_ron=="" || n_key_ron==" ")
  n_key_ron="RDNOISE"
gemhedit(n_ref//"_avg",n_key_ron,(sqrt(roneff)),"Effective read-noise [e-]")

}  # end of fl_avg


# Delete temporary masks from DQ
if(n_immasks=="DQ" && n_fl_msk) 
  imdelete("@"//tmpmli,verify=no, >>& "dev$null")

delete(tmplis//","//tmpmli,verify=no, >>& "dev$null")

# update output headers to contain input images
# tmpout contains  image_trn of all images
if(n_fl_med || n_fl_add || n_fl_avg) {
  n_i=1
  inimag=tmpout
  while(fscan(inimag,n_inim)!=EOF) {
   if(n_i<10) {
     if(n_fl_med)
       gemhedit(n_ref//"_med","IMAGE0"//str(n_i),
         substr(n_inim,1,strlen(n_inim)-4),"Input image for imcoadd")
     if(n_fl_add)
       gemhedit(n_ref//"_add"//l_zero_ext,"IMAGE0"//str(n_i),
         substr(n_inim,1,strlen(n_inim)-4),"Input image for imcoadd")
     if(n_fl_avg)
       gemhedit(n_ref//"_avg","IMAGE0"//str(n_i),
         substr(n_inim,1,strlen(n_inim)-4),"Input image for imcoadd")
   } else {
     if(n_fl_med)
       gemhedit(n_ref//"_med","IMAGE"//str(n_i),
         substr(n_inim,1,strlen(n_inim)-4),"Input image for imcoadd")
     if(n_fl_add)
       gemhedit(n_ref//"_add"//l_zero_ext,"IMAGE"//str(n_i),
         substr(n_inim,1,strlen(n_inim)-4),"Input image for imcoadd")
     if(n_fl_avg)
       gemhedit(n_ref//"_avg","IMAGE"//str(n_i),
         substr(n_inim,1,strlen(n_inim)-4),"Input image for imcoadd")
   }
   n_i+=1
  }
  inimag=""
  
  delete(tmpout,verify=no, >>& "dev$null")

  if(n_outimage!="" && n_imtype==".fits" && n_ref != n_outimage) {
#    movefiles(n_ref//"_add.fits",n_outimage//".fits")
      rename(n_ref//"_add.fits",n_outimage//".fits", field="all")
  } else if (n_outimage!="" && n_imtype==".imh") {
      imrename(n_ref//"_add",n_outimage)
  }


}
goto clean

crash:
  status=1
  goto clean

clean:
inimag="" ; inmask="" ; input=""
delete(tmplis//","//tmpmli//","//tmpref,verify-, >>& "dev$null")
if(status==0)
  printlog("IMCOADD exit status: good",l_logfile,l_verbose)
else
  printlog("IMCOADD exit status: error",l_logfile,l_verbose)
printlog("-----------------------------------------------------------------------------",l_logfile,l_verbose)

end
