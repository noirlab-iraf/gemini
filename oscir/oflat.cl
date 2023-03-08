# Copyright(c) 2001-2011 Association of Universities for Research in Astronomy, Inc.

procedure oflat(skyflat,polyflat,flatimage)

# Derive flat fields for OSCIR images.
#
# Version: Sept 14, 2002 IJ  Release v1.4
#          Aug 20, 2003  KL  IRAF2.12 - new/modified parameters
#                              hedit: addonly
#                              imcombine: headers,bpmasks,expmasks,outlimits
#                                 rejmask->rejmasks, plfile->nrejmasks
#                              imstat: nclip,lsigma,usigma,cache

char  skyflat     {prompt="Input OSCIR sky flatfield"}
char  polyflat    {prompt="Input OSCIR polystyrene flatfield"}
char  flatimage   {prompt="Output flat field image"}
char  flattitle   {"",prompt="Title for output flat image"}
char  logfile     {"",prompt="Logfile"}
bool  verbose     {no,prompt="Verbose"}
int   status      {0,prompt="Exit status (0=good)"}
struct* scanfile

begin

char l_skyflat, l_polyflat, l_logfile, l_flatimage, l_flattitle
bool l_verbose
struct l_struct

int  n_saveset, n_i
char tmpin, tmpimage, l_test, l_junk
real adc_sat, n_norm

# Detector parameters
#adc_dark = 1.315E4  # This must be for a certain frame time ?
adc_sat = 5.793E4

status=0
tmpin = mktemp("tmpin")
tmpimage = mktemp("tmpim")

# cache imgets - used throughout the script
cache("imgets","tstat", "gemdate")

# set the local variables
l_skyflat=skyflat ; l_polyflat=polyflat ; l_flatimage=flatimage
l_flattitle=flattitle
l_verbose=verbose ; l_logfile=logfile 

# Check for package log file or user-defined log file
cache("oscir")
if((l_logfile=="") || (l_logfile==" ")) {
   l_logfile=oscir.logfile
   if((l_logfile=="") || (l_logfile==" ")) {
      l_logfile="oscir.log"
      printlog("WARNING - OFLAT:  Both oreduce.logfile and oscir.logfile are empty.",logfile=l_logfile, verbose+)
      printlog("                  Using default file oscir.log.",
         logfile=l_logfile, verbose+)
   }
}
# Open log file
date | scan(l_struct)
printlog("----------------------------------------------------------------------------",
   logfile=l_logfile, verbose=l_verbose)
printlog("OFLAT -- "//l_struct, logfile=l_logfile, verbose=l_verbose)
printlog(" ",logfile=l_logfile, verbose=l_verbose)

#-----------------------------------------------------------------------
# check input and output names

if(l_skyflat=="" || l_skyflat==" ") {
   printlog("ERROR - OFLAT: Skyflat defined as empty string.",logfile=l_logfile,verbose+)
   status=1
   goto clean
}
if(l_polyflat=="" || l_polyflat==" ") {
   printlog("ERROR - OFLAT: Polystyrene image defined as empty string.",logfile=l_logfile,verbose+)
   status=1
   goto clean
}
if(l_flatimage=="" || l_flatimage==" ") {
   printlog("ERROR - OFLAT: Output flat image defined as empty string.",logfile=l_logfile,verbose+)
   status=1
   goto clean
}


if(!imaccess(l_skyflat)) {
   printlog("ERROR - OFLAT: Input image "//l_skyflat//" not found.",logfile=l_logfile,verbose+)
   status=1
   goto clean
}
if(!imaccess(l_polyflat)) {
   printlog("ERROR - OFLAT: Input image "//l_polyflat//" not found.",logfile=l_logfile,verbose+)
   status=1
   goto clean
}
if(imaccess(l_flatimage)) {
   printlog("ERROR - OFLAT: Output image "//l_flatimage//" exists.",logfile=l_logfile,verbose+)
   status=1
   goto clean
}

# Use this to figure out if the input files are valid FITS files
l_junk="" ; l_test=""
imhead(l_skyflat,imlist = "*.imh,*.fits,*.pl,*.qp,*.hhh",long-,
   userfields+) |& scan(l_junk,l_test)
if(l_test=="Negative") {
  printlog("ERROR - OFLAT: Image "//l_skyflat//" not a valid FITS file",
  l_logfile,verbose+)
  status=1
  goto clean
}
l_junk="" ; l_test=""
imhead(l_polyflat,imlist = "*.imh,*.fits,*.pl,*.qp,*.hhh",long-,
   userfields+) |& scan(l_junk,l_test)
if(l_test=="Negative") {
  printlog("ERROR - OFLAT: Image "//l_polyflat//" not a valid FITS file",
  l_logfile,verbose+)
  status=1
  goto clean
}


printlog("Sky flat field    : "//l_skyflat,logfile=l_logfile,verbose=l_verbose)
printlog("Polystyrene image : "//l_polyflat,logfile=l_logfile,verbose=l_verbose)
printlog("Output flat field : "//l_flatimage,logfile=l_logfile,verbose=l_verbose)

#--------------------------------------------------------------------------
# The math and bookkeeping

# Check/Get the dimensions
imgets(l_skyflat,"i_naxis3")
  if(imgets.value!="1") {
   printlog("ERROR - OFLAT: Number of chop positions != 1 for image"//l_skyflat,l_logfile,verbose+)
   status=1 
   goto clean
  }
imgets(l_skyflat,"i_naxis5")
  if(imgets.value!="1") {
   printlog("ERROR - OFLAT: Number of nod positions != 1 for image"//l_skyflat,l_logfile,verbose+)
   status=1 
   goto clean
  }
imgets(l_skyflat,"i_naxis6")
  if(imgets.value!="1") {
   printlog("ERROR - OFLAT: Number of nod sets != 1 for image"//l_skyflat,l_logfile,verbose+)
   status=1 
   goto clean
  }

imgets(l_polyflat,"i_naxis3")
  if(imgets.value!="1") {
   printlog("ERROR - OFLAT: Number of chop positions != 1 for image"//l_polyflat,l_logfile,verbose+)
   status=1 
   goto clean
  }
imgets(l_polyflat,"i_naxis5")
  if(imgets.value!="1") {
   printlog("ERROR - OFLAT: Number of nod positions != 1 for image"//l_polyflat,l_logfile,verbose+)
   status=1 
   goto clean
  }
imgets(l_polyflat,"i_naxis6")
  if(imgets.value!="1") {
   printlog("ERROR - OFLAT: Number of nod sets != 1 for image"//l_polyflat,l_logfile,verbose+)
   status=1 
   goto clean
  }

imgets(l_skyflat,"i_naxis4")
  if(imgets.value=="0") {
   printlog("ERROR - OFLAT: No savesets",l_logfile,verbose+)
   status=1 
   goto clean
  } else
   n_saveset=int(imgets.value)

# Combine all frames for skyflat
for(n_i=1;n_i<=n_saveset;n_i+=1) {
  print(l_skyflat//"[*,*,1,"//str(n_i)//",1,1]", >> tmpin)
}
imcombine("@"//tmpin,l_flatimage,headers="",bpmasks="",rejmasks="",
   nrejmasks="",expmasks="",sigmas="",logfile="",combine="average",
   reject="none",project=no,outtype="real",outlimits="",offsets="none",
   masktype="none",maskvalue=0.,blank=0.,scale="none",zero="none",weight="none",
   statsec="",expname="",lthreshold=INDEF,hthreshold=INDEF)
delete(tmpin,ver-)
imstat(l_flatimage,fields="mean",lower=INDEF,upper=INDEF,nclip=0,
   lsigma=INDEF,usigma=INDEF,binwidth=0.1,format-,cache-) | scan(n_norm)
imgets(l_skyflat,"FRMCOADD")
if(n_norm>adc_sat*real(imgets.value)) {
  printlog("ERROR - OFLAT: Skyflat is saturated",l_logfile,verbose+)
  print(n_norm)
  status=1
  goto clean
}

# Combine all frames for polyflat
imgets(l_polyflat,"i_naxis4")
  if(imgets.value=="0") {
   printlog("ERROR - OFLAT: No savesets",l_logfile,verbose+)
   status=1 
   goto clean
  } else
   n_saveset=int(imgets.value)

# Combine all frames for skyflat
for(n_i=1;n_i<=n_saveset;n_i+=1) {
  print(l_polyflat//"[*,*,1,"//str(n_i)//",1,1]", >> tmpin)
}
imcombine("@"//tmpin,tmpimage,headers="",bpmasks="",rejmasks="",nrejmasks="",
   expmasks="",sigmas="",logfile="",combine="average",reject="none",project=no,
   outtype="real",outlimits="",offsets="none",masktype="none",maskvalue=0.,
   blank=0.,scale="none",zero="none",weight="none",statsec="",expname="",
   lthreshold=INDEF,hthreshold=INDEF)
imstat(tmpimage,fields="mean",lower=INDEF,upper=INDEF,nclip=0,lsigma=INDEF,
   usigma=INDEF,binwidth=0.1,format-,cache-) | scan(n_norm)
imgets(l_polyflat,"FRMCOADD")
if(n_norm>adc_sat*real(imgets.value)) {
  printlog("ERROR - OFLAT: Polystyrene image is saturated",l_logfile,verbose+)
  status=1
  goto clean
}
imarith(l_flatimage,"-",tmpimage,l_flatimage,title="",divzero=0,hparams="",
  pixtype="real",calctype="real",verbose-,noact-)
imdelete(tmpimage,verify-)
delete(tmpin,verify-)

#-------------
# Normalization
imstat(l_flatimage,fields="mean",lower=INDEF,upper=INDEF,nclip=0,lsigma=INDEF,
   usigma=INDEF,binwidth=0.1,format-,cache-) | scan(n_norm)
imarith(l_flatimage,"/",n_norm,l_flatimage,title="",divzero=0,hparams="",pixtype="",
  calctype="",verbose-,noact-)

printf("Normalization value = %10.0f\n",abs(n_norm)) | scan(l_struct)
printlog(l_struct,l_logfile,l_verbose)

# title
if(l_flattitle=="" || l_flattitle==" " || l_flattitle=="default")
  gemhedit (l_flatimage, "i_title", "Flatfield from OFLAT", "", delete-)
else
  gemhedit (l_flatimage, "i_title", l_flattitle, "", delete-)

# Update the header
gemdate ()
# date stamp the modification
gemhedit(l_flatimage,"GEM-TLM",gemdate.outdate,"UT Last modification with GEMINI")
gemhedit(l_flatimage,"OFLAT",gemdate.outdate,"UT Time stamp for oflat")
gemhedit(l_flatimage,"OFLTNORM",abs(n_norm),"Normalization value")

#---------------------------------------------------------------------------
# Clean up
clean:
{
  if(status==0) {
    printlog("OFLAT exit status:  good.",logfile=l_logfile, verbose=l_verbose)
  }
   printlog("----------------------------------------------------------------------------", logfile=l_logfile, verbose=l_verbose)
  scanfile=""
  delete(tmpin,ver-, >& "dev$null")
}
end
