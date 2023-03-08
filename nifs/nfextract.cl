# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

procedure nfextract (inimages)

# Routine to extract a 1D spectrum from a transformed NIFS MEF file.
#
# Version December, 2005  I. Song 
#
# History: 
#
# 25-MAR-05  PJM 		- Creation based on nsextract.cl.
#
# DEC-05     I. Song	- modify to work with the Gemini NIFS package
#                       - create SCI, VAR, DQ MEF file for the output.
#
# MAR-06     I. Song    - eliminate a few unused variables.
#                       - clean up the code.
#                       - add an option to extract spectrum automatically
#                         at user input position (xc,yc)
#
# (xc, yc) can be determined from an output from 'nfimage' 
#          xc = nfimage_x / 5.0, yc = nfimage_y / 2.0
#
# MAR-06     C.Aspin    - added option to plot image with z1,z2 specified
#
# Still to do:
#   - add more detailed reduction history into headers?
#

char inimages    {           prompt="Input NIFS images or list"}
char outspectra  {"",        prompt="Output spectra"}
char outprefix   {"x",       prompt="Prefix for output spectra"}
real diameter    {0.5,       prompt="Aperture diameter (arcsec)\n"}

bool fl_inter    {yes,       prompt="Mark source position interactively"}
real xc          {INDEF,     prompt="X-position (pixels) for automatic extraction"}
real yc          {INDEF,     prompt="Y-position (pixels) for automatic extraction\n"}

bool fl_zval     {no,        prompt="Input z1,z2 values for image display scaling"}
real z1          {0.0,       prompt="z1 display scaling value"}
real z2          {10000.0,   prompt="z2 display scaling value"}
int  dispaxis    {1,         prompt="Dispersion axis",min=1,max=2}

char logfile     {"",        prompt="Logfile name"}
bool verbose     {yes,       prompt="Verbose output"}
int  status      {0,         prompt="Exit status (0=good)"}
struct *scanfile1{"",        prompt="For internal use only"}
struct *scanfile2{"",        prompt="For internal use only"}

begin
# Define local variables.
int  l_dispaxis
real l_xc, l_yc
real l_diameter, l_z1, l_z2
char l_inimages, l_outspectra, l_outprefix, l_logfile
char l_sci_ext, l_var_ext, l_dq_ext
bool l_verbose, l_fl_inter, l_fl_zval

# Define run time variables.
int  nbad, nin, nout, junk, i, j, numext, nx, ny
int  nxpix, nypix, nzpix, n1, n2
real rad, dx, dy, x, y, rad2, dist
char s_empty, img, phu,  extn, text, imgin, imgout
char tmpfile, tmpextn, tmpinlist, tmprootlist, tmpoutlist, tmplist, tmpcube
char scispe, ver, sec, secvar, secphu, tmpver, tmpsec, tmpmsk, tmp, tmpimg
char varspec, dqspec, tmpvarver, tmpdqver
bool used, first

struct sdate

# Set local variable values.
l_inimages     = inimages
l_outspectra   = outspectra
l_outprefix    = outprefix
l_diameter     = diameter
l_xc           = xc
l_yc           = yc
l_fl_inter     = fl_inter
l_fl_zval      = fl_zval
l_z1           = z1
l_z2           = z2
l_dispaxis     = dispaxis
l_logfile      = logfile
l_verbose      = verbose
l_sci_ext      = nsheaders.sci_ext
l_var_ext      = nsheaders.var_ext
l_dq_ext       = nsheaders.dq_ext

if (!l_fl_inter && (INDEF == l_xc || INDEF == l_yc)) {
   printlog("-----------------------------------------------------------------"//
     "-----------",l_logfile,l_verbose)
   printlog("NFEXTRACT -- "//sdate,l_logfile,l_verbose)
   printlog("xc and yc should be defined if not in the interactive mode!",l_logfile,l_verbose)
   goto crash
}

# Make temporary file names.
tmp            = mktemp("tmp$nfs")
tmpfile        = mktemp("tmp$nfs")
tmpextn        = mktemp("tmp$nfs")
tmpinlist      = mktemp("tmp$nfs")
tmprootlist    = mktemp("tmp$nfs")
tmpoutlist     = mktemp("tmp$nfs")
tmplist        = mktemp("tmp$nfs")
tmpcube        = mktemp("tmp$nfs")
tmpimg         = mktemp("tmp$nfs")
tmpmsk         = mktemp("tmp$nfs")
scispe         = mktemp("tmp$nfs")
varspec        = mktemp("tmp$nfs")
dqspec         = mktemp("tmp$nfs")

# Initialize variables.
status = 0

# Keep task parameters from changing from the outside.
cache ("gemextn", "keypar", "gemdate") 

# Test the logfile.
s_empty=""; print(l_logfile) | scan(s_empty); l_logfile=s_empty
if (l_logfile == "" || stridx(" ",l_logfile) > 0) {
  l_logfile = nifs.logfile
  if (l_logfile == "" || stridx(" ",l_logfile) > 0) {
    l_logfile = "nifs.log"
    printlog("WARNING - NFEXTRACT: Both nfextract.logfile and "//
      "nifs.logfile are empty.",l_logfile,l_verbose)
    printlog("                     Using default file nifs.log.",l_logfile,
      l_verbose)
  }
}

# Start logging.
date | scan(sdate)
printlog("-----------------------------------------------------------------"//
  "-----------",l_logfile,l_verbose)
printlog("NFEXTRACT -- "//sdate,l_logfile,l_verbose)
printlog("",l_logfile,l_verbose)

# Logs the relevant parameters:
printlog("Input images or list = "//l_inimages,l_logfile,l_verbose)
printlog("Output spectra       = "//l_outspectra,l_logfile,l_verbose)
printlog("Output prefix        = "//l_outprefix,l_logfile,l_verbose)
printlog("Dispersion axis      = "//l_dispaxis,l_logfile,l_verbose)
if (!l_fl_inter) {
  printlog("Auto extract at X    = "//l_xc,l_logfile,l_verbose)
  printlog("Auto extract at Y    = "//l_yc,l_logfile,l_verbose)
}
printlog("",l_logfile,l_verbose)

# Make lists of input images, putting all input images in a temporary file.
gemextn(l_inimages,proc="none",check="mef,exists",index="",extname="",
  extver="",ikparam="",replace="",omit="extension",outfile=tmpfile,logfile="",
  glogpars="",verbose-)
if (gemextn.fail_count != 0 || gemextn.count == 0) {
  printlog("ERROR - NFEXTRACT: Input images are missing or are not MEF.",
    l_logfile,yes)
  goto crash
}

# Filter/expand to tmpinlist, if suitably processed.
nbad = 0
nin = 0
scanfile1 = tmpfile
while (fscan(scanfile1,img) != EOF) {
  phu = img//"[0]"
  keypar(phu,"PREPARE", silent+)
  if (!keypar.found) {
    printlog("ERROR - NFEXTRACT: Image "//img//" has not been NFPREPARE'd.",
      l_logfile,yes)
    nbad += 1
  } else {
    used = yes
    delete(tmpextn,ver-, >& "dev$null")
    gemextn(img,proc="expand",check="exists",index="",extname=l_sci_ext,
      extver="1-",omit="",replace="",ikparams="",outfile=tmpextn,
      logfile="dev$null",glogpars="",verbose-)
    numext = gemextn.count
    scanfile2 = tmpextn
    while (fscan(scanfile2,extn) != EOF) {
      text = substr(extn,stridx(",",extn)+1,stridx("]",extn)-1)
      version = int(text)
      keypar(img//"["//l_sci_ext//","//version//",inherit]","NSTRANSF",silent+)
      if (!keypar.found) {
        printlog("WARNING - NFEXTRACT: Extension \
          "//img//"["//l_sci_ext//","//version//",inherit] has not been "//
          "run through NSTRANSFORM before.",l_logfile,yes)
        nbad += 1
        used = no
      }
    }
    if (used) {
      print(img, >> tmprootlist)
      print(img//" "//numext, >> tmpinlist)
      nin += 1
    }
  }
}
if (nbad > 0) {
  printlog("ERROR - NFEXTRACT: "//nbad//" image(s) has not been "//
    "NFPREPARE'd and NSTRANSFORM'd.",l_logfile,yes)
  goto crash
}
if (nin == 0) {
  printlog("ERROR - NFEXTRACT: No input images to process.",l_logfile,yes)
  goto crash
}
delete(tmpfile,ver-, >& "dev$null") 

# Check output images.
gemextn(l_outspectra,check="absent",process="none",index="",extname="",
  extversion="",ikparams="",omit="kernel,exten",replace="",outfile=tmpoutlist,
  logfile="",glogpars="",verbose=l_verbose)
if (gemextn.fail_count != 0) {
  printlog("ERROR - NFEXTRACT: Existing or incorrectly formatted output "//
    "files.",l_logfile,yes)
  goto crash
}
# If tmpoutlist is empty, the output files names should be created with the 
# prefix parameter. (We could have a separate message for outprefix="", but 
# that will trigger an error in gemextn anyway.)
if (gemextn.count == 0) {
  gemextn("%^%"//l_outprefix//"%"//"@"//tmprootlist,check="absent",process="none",
    index="",extname="",extversion="",ikparams="",omit="kernel,exten",
    replace="",outfile=tmpoutlist,logfile="",glogpars="",verbose=l_verbose)
  if (gemextn.fail_count != 0 || gemextn.count == 0) {
    printlog("ERROR - NFEXTRACT: No or incorrectly formatted output files",
      l_logfile,yes)
    goto crash
  }
}
nout = gemextn.count
delete(tmprootlist,ver-, >& "dev$null")

# Check number of input and output images.
if (nin != nout) {
  printlog("ERROR - NFEXTRACT: Different number of input and output images.", 
    l_logfile,yes)
  goto crash
}

# Add output files to temporary file list.
scanfile1 = tmpinlist
scanfile2 = tmpoutlist
while (fscan(scanfile1,imgin,numext) != EOF) {
  junk = fscan(scanfile2,imgout)
  print(imgin//" "//numext//" "//imgout, >> tmplist)
}

# Process input files.
rad = l_diameter/2.0
scanfile1 = tmplist
while (fscan(scanfile1,imgin,numext,imgout) != EOF) {

  printlog("",l_logfile,yes)

  ver    = ","//version
  sec    = "["//l_sci_ext//ver//"]"
  secvar = "["//l_var_ext//ver//"]"
  secphu = "["//l_sci_ext//ver//",inherit]"

  keypar(imgin//secphu,nsheaders.key_dispaxis,silent+)
  if (keypar.found) {
    l_dispaxis = int(keypar.value)
  } else {
    l_dispaxis = dispaxis
  }

  # Copy PHU to the output file.
  imcopy(imgin//"[0]",imgout,verbose-)

  # Form a temporary 3D data cube.
  for (i=1; i<=numext; i+=1) {
    if (l_dispaxis == 1) {
      imgets(imgin//"["//l_sci_ext//","//i//"]","i_naxis1", >& "dev$null")
    } else {
      imgets(imgin//"["//l_sci_ext//","//i//"]","i_naxis2", >& "dev$null")
    }
    nzpix=int(imgets.value)
    n1=1
    n2=min(2048,nzpix)
    print(imgin//"["//l_sci_ext//","//i//"]["//n1//":"//n2//",*]", >> tmpfile)
  }
  imstack("@"//tmpfile,tmpcube,title="*",pixtype="double")
  delete(tmpfile,verify-, >& "dev$null")
  if (l_dispaxis == 1) {
    imgets(tmpcube,"i_naxis3", >& "dev$null")
    nxpix=int(imgets.value)
    dx=0.103
    imgets(tmpcube,"i_naxis2", >& "dev$null")
    nypix=int(imgets.value)
    imgets(tmpcube,"CD2_2", >& "dev$null")
    dy=real(imgets.value)
  } else {
    imgets(tmpcube,"i_naxis3", >& "dev$null")
    nxpix=int(imgets.value)
    dx=0.103
    imgets(tmpcube,"i_naxis1", >& "dev$null")
    nypix=int(imgets.value)
    imgets(tmpcube,"CD1_1", >& "dev$null")
    dy=real(imgets.value)
  }

  # Project 3D cube to 2D image.
  improject(tmpcube,tmpimg,projaxis=l_dispaxis,average-,
    highcut=100000000.0,lowcut=-100000000.0,pixtype="double",verbose-)
  # Clunky bug-fix.
  #hedit(tmpimg,"WAT2_001","",del+,update+,show-,verify-, >& "dev$null")
  rotate(tmpimg,tmpimg,90.0,xin=INDEF,yin=INDEF,xout=INDEF,yout=INDEF,
    ncols=0,nlines=0,interpolant="linear",boundary="nearest",constant=0,
    nxblock=512,nyblock=512, >& "dev$null")
  gemhedit (tmpimg, "CTYPE2", "LINEAR", "", delete-, >& "dev$null")
  imcopy(tmpimg,tmpmsk,ver-, >& "dev$null")
  imreplace(tmpmsk,value=1.0,imaginary=0.0,lower=INDEF,upper=INDEF,radius=0.0)

  xc = l_xc   # set xc to the user input position for automatic source extraction.
  yc = l_yc   # ..for interactive extraction, (xc,yc) will be overwritten later.

  if (l_fl_inter) {
     if (l_fl_zval) {
        # Display 2D image with z1,z2 scaling
        printlog("Displaying image with z1="//l_z1//" and z2="//l_z2,l_logfile,l_verbose)
        display(tmpimg,1,bpmask="BPM",bpdisplay="overlay",bpcolor="red",overlay="",
          erase+,border-,select+,repeat-,fill-,zscale-,zrange-,xmag=5.0,ymag=2.0,
          ztrans="linear",z1=l_z1,z2=l_z2, >& "dev$null")
     } else {
        # Display 2D image with zscale+
        printlog("Displaying image with zscale+",l_logfile,l_verbose)
        display(tmpimg,1,bpmask="BPM",bpdisplay="overlay",bpcolor="red",overlay="",
          erase+,border-,select+,repeat-,fill-,zscale+,zrange-,xmag=5.0,ymag=2.0,
          ztrans="linear", >& "dev$null")
     }
     # Define aperture center coordinates.
     xc = real(nxpix)/2.0
     yc = real(nypix)/2.0
     printlog("Press any key to mark aperture center with cursor",
       l_logfile,l_verbose)
     junk = fscan(imcur,xc,yc)
     delete(tmp,ver-, >& "dev$null")
     print(xc,yc,"100","v", >> tmp)
     x=xc+rad/dx
     print(x,yc,"100","v", >> tmp)
     tvmark(1,coords="",logfile="",autolog-,outimage="",deletions="",
        commands=tmp,mark="circle",radii=rad/dx,lengths=0,font="raster",
        color=208,label-,number-,nxoff=0,nyoff=0,point=3,txsize=1,
        tol=1.5,inter-)
  }
  printlog("Extracting aperture at position xc="//xc//", yc="//yc,l_logfile,
    l_verbose)

  # Extract spectrum from a circular aperture.
  first=yes
  rad2=rad**2
  xc=real(nxpix+1)-xc
  for (i=1; i<=nxpix; i+=1) {
    for (j=1; j<=nypix; j+=1) {
      dist=( (real(i)-xc)*dx )**2+( (real(j)-yc)*dy )**2
      if (dist < rad2) {
        tmpver    = "["//l_sci_ext//","//i//"]"
        tmpvarver = "["//l_var_ext//","//i//"]"
        tmpdqver  = "["//l_dq_ext//","//i//"]"
        if (l_dispaxis == 1) {
          tmpsec = "[*,"//j//":"//j//"]"
        } else {
          tmpsec = "["//j//":"//j//",*]"
        }
        if (first) {
          imcopy(imgin//tmpver//tmpsec,   scispe,  >& "dev$null")
          imcopy(imgin//tmpvarver//tmpsec,varspec, >& "dev$null")
          imcopy(imgin//tmpdqver//tmpsec, dqspec,  >& "dev$null")
          first = no
        } else {
           imarith(scispe, "+",imgin//tmpver//tmpsec,   scispe)
           imarith(varspec,"+",imgin//tmpvarver//tmpsec,varspec)
           imarith(dqspec, "+",imgin//tmpdqver//tmpsec, dqspec)
        }
        imreplace(tmpmsk//"["//i//":"//i+1//","//j//":"//j+1//"]",value=0.0,
          imaginary=0.0,lower=INDEF,upper=INDEF,radius=0.0)
      }
    }
  }

  # Create the output spectrum.
  imcopy(scispe, imgout//"["//l_sci_ext//",1,append]",ver-, >& "dev$null")
  imcopy(varspec,imgout//"["//l_var_ext//",1,append]",ver-, >& "dev$null")
  imcopy(dqspec, imgout//"["//l_dq_ext //",1,append]",ver-, >& "dev$null")

  # Update header of output image (needs much more!).
  gemdate ()
  gemhedit (imgout//"[0]", "NFEXTRAC", gemdate.outdate,
    "UT Time stamp for NFEXTRACT", delete-)
  gemhedit (imgout//"[0]", "GEM-TLM", gemdate.outdate,
    "UT Last modification with GEMINI", delete-)

  # Clean up.
  delete(tmp,ver-, >& "dev$null")
  imdelete(tmpcube,ver-, >& "dev$null")
  imdelete(tmpimg,ver-, >& "dev$null")
  imdelete(tmpmsk,ver-, >& "dev$null")
  imdelete(scispe,ver-, >& "dev$null")
  imdelete(dqspec,ver-, >& "dev$null")
  imdelete(varspec,ver-, >& "dev$null")

  # Plot the spectrum.
  if (l_fl_inter) {
     printlog("Plotting the spectrum.  Press q to quit and continue.",l_logfile,l_verbose)
  splot(imgout//"["//l_sci_ext//"]",units="",options="auto,zero,xydraw,histogram",
        xmin=INDEF, xmax=INDEF, ymin=INDEF, ymax=INDEF)
  }

}

goto clean

###############################################################################
# Exit with error.
crash:
  status=1

# Clean up and exit.
clean:
  scanfile1 = ""
  scanfile2 = ""

  delete(tmp//","//tmpinlist//","//tmpoutlist//","//tmpfile,ver-, >& "dev$null")
  delete(tmpextn//","//tmprootlist,ver-, >& "dev$null")

  printlog("",l_logfile,l_verbose)
  if (status == 0) {
    printlog("NFEXTRACT exit status:  good.",l_logfile,l_verbose)
  } else {
    printlog("NFEXTRACT exit status:  error.",l_logfile,l_verbose)
  }
  printlog("-------------------------------------------------------------"//
    "---------------",l_logfile,l_verbose)
  
end
