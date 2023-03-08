# Copyright(c) 2000-2011 Association of Universities for Research in Astronomy, Inc.

procedure qfastsky(inimages,outimage)

# Make a quick sky image for QUIRC images.
# Version  Jan 22, 2001  JJ
#          Feb 15, 2001  JJ, added min for nlow and nhigh
#          Aug 20, 2003  KL, IRAF2.12 - new/modified parameters
#                              hedit: addonly
#                              imcombine: headers,bpmasks,expmasks,outlimits
#                                 rejmask->rejmasks, plfile->nrejmasks

char    inimages  {prompt="Raw QUIRC images to combine"}
char    outimage  {prompt="Output sky image"}
char    outtitle  {"default",prompt="Title for output image"}
char    key_exptime {"EXPTIME",prompt="Header keyword for exposure time"}
char    combtype  {"default",prompt="Type of combine operation",
                    enum="default|median|average"}
char    rejtype   {"minmax",prompt="Type of rejection", enum="none|minmax"}
char    logfile   {"",prompt="Name of log file"}
int     nlow      {0,min=0,prompt="minmax: Number of low pixels to reject"}
int     nhigh     {1,min=0,prompt="minmax: Number of high pixels to reject"}
bool    verbose   {no,prompt="Verbose actions"}
int     status    {0,prompt="Exit status (0=good)"}
struct* scanfile  {prompt="Internal use only"}

begin

char    l_inimages, l_outimage, l_combtype, l_rejtype, l_logfile
char    l_key_exptime, l_outtitle
int     l_nlow, l_nhigh
bool    l_verbose
int     i, nimages
char    temp, l_extension, img, tmpfile1, tmpfile2
struct  l_struct
real    l_expzero, l_expone

status=0

# Set local variables
l_inimages = inimages ; l_outimage = outimage ; l_key_exptime = key_exptime
l_combtype = combtype ; l_rejtype = rejtype ; l_outtitle = outtitle
l_logfile = logfile ; l_verbose=verbose
l_nlow = nlow ; l_nhigh = nhigh
cache("quirc","imgets", "gemdate")
if((l_logfile=="") || (l_logfile==" ")) {
   l_logfile=quirc.logfile
   if((l_logfile=="") || (l_logfile==" ")) {
      l_logfile="quirc.log"
      printlog("WARNING - QFASTSKY:  Both qfastsky.logfile and quirc.logfile are empty.",logfile=l_logfile, verbose+)
      printlog("                     Using default file quirc.log.",
         logfile=l_logfile, verbose+)
   }
}

# Open log file
date | scan(l_struct)
printlog("----------------------------------------------------------------------------",logfile=l_logfile, verbose=l_verbose)
printlog("QFASTSKY -- "//l_struct, logfile=l_logfile, verbose=l_verbose)
printlog(" ",logfile=l_logfile, verbose=l_verbose)

# Make temporary files
tmpfile1 = mktemp("tmpfl")
tmpfile2 = mktemp("tmpfl")

if(l_inimages=="" || l_inimages==" ") {
  printlog("ERROR - QFASTSKY: No input images defined",logfile=l_logfile,verbose+)
  status=1
  goto clean
}
if(l_outimage=="" || l_outimage==" ") {
  printlog("ERROR - QFASTSKY: No output image defined",logfile=l_logfile,verbose+)
  status=1
  goto clean
}

show imtype | scan(l_extension)
#check to see if output image already exists; strip .fits if present
if(substr(l_outimage,strlen(l_outimage)-4,strlen(l_outimage)) == ".fits") {
  l_outimage=substr(l_outimage,1,(strlen(l_outimage)-5))
}

if(imaccess(l_outimage) ) {
  printlog("ERROR - QFASTSKY: Output image "//l_outimage//" already exists",
     logfile=l_logfile, verbose+)
  status=1
  goto clean
  }

# Put all images in a temporary file: tmpfile1
if(substr(l_inimages,1,1)=="@") 
  type(substr(l_inimages,2,strlen(l_inimages)), > tmpfile1)
else if (stridx("*", l_inimages)>0)
  files(l_inimages,sort-, > tmpfile1)
else 
  files(l_inimages,sort-, > tmpfile1)

# Verify that input images actually exist. (Redundant for * though)
# at the same time check the exposure times 
l_expzero = -999.
scanfile = tmpfile1
while(fscan(scanfile,img) != EOF) {
  if(imaccess(img)) {
    print(img, >> tmpfile2)
    imgets(img,l_key_exptime, >& "dev$null")
      if(imgets.value=="0") {
        printlog("WARNING - QFASTSKY: No exposure time defined for "//img,
           logfile=l_logfile,verbose+)
      } else {
    l_expone=real(imgets.value)
      if((abs(l_expzero-l_expone)>0.1) && (l_expzero>-1))
        printlog("WARNING - QFASTSKY: Input images have different exposure times",
           logfile=l_logfile,verbose+)
      l_expzero=l_expone
      }
    }
  else
    printlog("WARNING - QFASTSKY: image "//img//" does not exist",
       logfile=l_logfile,verbose+)
}
scanfile=""

printlog("Using input files:",logfile=l_logfile, verbose=l_verbose)
if(l_verbose) type(tmpfile2)
type(tmpfile2, >> l_logfile)
printlog("Output image: "//l_outimage,logfile=l_logfile, verbose=l_verbose)

# Get the number of images
nimages=0
if(access(tmpfile2))
  count(tmpfile2) | scan(nimages) 
  
if (nimages == 1) {
  printlog("ERROR - QFASTSKY: Cannot combine a single image",
     logfile=l_logfile,verbose+)
  status=1
  goto clean
} else if (nimages == 0) {
  printlog("ERROR - QFASTSKY: No images to combine!",logfile=l_logfile,verbose+)
  status=1
  goto clean
}

# save the user's parameters for imcombine
delete("uparm$imhimcome.par.org",ver-, >& "dev$null")
if(access("uparm$imhimcome.par"))
  copy("uparm$imhimcome.par","uparm$imhimcome.par.org",verbose=yes)

cache("imcombine")
imcombine.headers=""
imcombine.bpmasks=""
imcombine.rejmasks=""
imcombine.nrejmasks=""
imcombine.expmasks=""
imcombine.sigmas=""
imcombine.logfile = "STDOUT"   # so that it can be tee'd later
# imcombine.combine defined later
# imcombine.reject defined later
imcombine.project = no
imcombine.outtype = "real"
imcombine.outlimits = ""
imcombine.offsets = "none"
imcombine.masktype = "none"
imcombine.maskvalue=0.
imcombine.blank = 0.
imcombine.scale = "none"
imcombine.zero = "median"
imcombine.weight = "none"
imcombine.statsec = ""
imcombine.expname = ""
imcombine.lthreshold = INDEF
imcombine.hthreshold = INDEF
# imcombine.nlow defined later
# imcombine.nhigh defined later
imcombine.nkeep=1
imcombine.grow = 0

if(l_combtype=="default") {
  if (nimages < 5) {
     imcombine.combine = "average"
     imcombine.reject = "minmax"
     imcombine.nlow=0
     imcombine.nhigh=1
     printlog("WARNING - QFASTSKY:  Averaging "//nimages//" images with 1 high pixel rejected",logfile=l_logfile,verbose+)
  } else if(nimages < 8) {
     imcombine.combine = "median"
     imcombine.reject = "minmax"
     imcombine.nlow=1
     imcombine.nhigh=1
  } else {
     imcombine.combine = "median"
     imcombine.reject = "minmax"
     imcombine.nlow=1
     imcombine.nhigh=2
  } 
} else {
  imcombine.combine = l_combtype
  imcombine.reject= l_rejtype
  imcombine.nlow = l_nlow
  imcombine.nhigh = l_nhigh
  if(nimages < 5) {
     printlog("WARNING - QFASTSKY:  Combining 4 or fewer images using "//l_combtype,logfile=l_logfile,verbose+)
     if(l_rejtype != "none") {
        printlog("                      with "//l_nlow//" low and "//l_nhigh//" high pixels rejected.", logfile=l_logfile,verbose+)
     } else {
        printlog("                      with no pixels rejected.",
           logfile=l_logfile,verbose+)
     }
  } # end if(nimages < 5)
  if((nimages <= (l_nlow+l_nhigh)) && (l_rejtype=="minmax")) {
     printlog("ERROR - QFASTSKY: Cannot reject more pixels than the number of images.",logfile=l_logfile,verbose+)
     status=1
     goto clean
  }
} # end not-default section

# Do the combine, quietly
imcombine("@"//tmpfile2,l_outimage, >& "dev$null")

printlog("Combining "//str(nimages)//" images, using "//imcombine.combine,
   logfile=l_logfile, verbose=l_verbose)
printlog("Rejection is "//imcombine.reject//", with "//imcombine.nlow//" low and "//imcombine.nhigh//" high values rejected",
   logfile=l_logfile,verbose=l_verbose)

# time stamp output
gemdate ()
gemhedit(l_outimage,"GEM-TLM",gemdate.outdate,"UT Last modification with GEMINI")
gemhedit(l_outimage,"QFASTSKY",gemdate.outdate,"UT Time stamp for qfastsky")

# fix the title
if(l_outtitle=="default")
  gemhedit (l_outimage//".fits", "i_title", "SKY IMAGE from gemini.quirc.qfastsky",
     "", delete-)
else
  gemhedit (l_outimage//".fits", "i_title", l_outtitle, "", delete-)

# Clean up
clean:
{
  if(status==0) {
    printlog("QFASTSKY exit status:  good",logfile=l_logfile, verbose=l_verbose)
  }
   printlog("----------------------------------------------------------------------------", logfile=l_logfile, verbose=l_verbose)

   delete(tmpfile1//","//tmpfile2,ver-, >& "dev$null")

# return to default parameters for imcombine
   unlearn("imcombine")
# restore the user's parameters for imcombine
   if(access("uparm$imhimcome.par.org"))
     rename("uparm$imhimcome.par.org","uparm$imhimcome.par",field="all")
}

end
