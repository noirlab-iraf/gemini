# Copyright(c) 2000-2006 Association of Universities for Research in Astronomy, Inc.

procedure qreduce(inimages)

# Basic reductions of QUIRC images. Requires a flat and sky image.
# A constant sky value is added back on (fl_autosky+ or skylevel!=0).
#
# Flatfield is selected on the basis key_filter
# 6 different flat fields can be used.
#
# Version  Jan 26, 2001  JJ
#          Aug 20, 2003  KL IRAF2.12 - new parameters
#                              imstat: nclip,lsigma,usigma,cache

char  inimages     {prompt="Input QUIRC image(s)"}
char  outimages    {"",prompt="Output image(s)"}
char  outprefix    {"r",prompt="Prefix for output image(s)"}
char  skyimage     {"",prompt="Sky image to subtract"}
real  skylevel     {0.0,prompt="Constant sky level to add"}
bool  fl_autosky   {yes,prompt="Add median of the sky frame as a constant?"} 
char  logfile      {"",prompt="Logfile"}
bool  fl_sky       {yes,prompt="Do sky subtraction"}
bool  fl_flat      {yes,prompt="Do flat fielding"}
char  key_filter   {"FILTER",prompt="Keyword for filter id"}
char  flatimage1   {"flatJ",prompt="Flat field image no.1"}
char  filter1      {"J",prompt="Filter for flat field no.1"}
char  flatimage2   {"flatH",prompt="Flat field image no.2"}
char  filter2      {"H",prompt="Filter for flat field no.2"}
char  flatimage3   {"flatK",prompt="Flat field image no.3"}
char  filter3      {"K",prompt="Filter for flat field no.3"}
char  flatimage4   {"",prompt="Flat field image no.4"}
char  filter4      {"",prompt="Filter for flat field no.4"}
char  flatimage5   {"",prompt="Flat field image no.5"}
char  filter5      {"",prompt="Filter for flat field no.5"}
char  flatimage6   {"",prompt="Flat field image no.6"}
char  filter6      {"",prompt="Filter for flat field no.6"}
bool  verbose      {no,prompt="Verbose"}
int   status       {0,prompt="Exit status (0=good)"}
struct* scanfile

begin

char l_inimages, l_outimages, l_skyimage, l_flatimage, l_filter
char l_expression, l_prefix, l_logfile, l_temp, tmpin, tmpfile
char in[1000], out[1000]
char l_flatimage1, l_flatimage2, l_flatimage3, l_flatimage4, l_flatimage5 
char l_flatimage6, l_filter1, l_filter2, l_filter3, l_filter4, l_filter5
char l_filter6, l_keyfilter
int  i, nimages, noutimages, maxfiles
bool l_fl_sky, l_fl_flat, l_verbose, l_fl_skytemp, l_fl_flattemp
bool flatok, l_fl_first, l_fl_autosky, l_fl_dollar
real l_skylevel, l_mid, l_sig
struct l_struct

status=0
maxfiles=1000
tmpfile = mktemp("tmpin")
tmpin = mktemp("tmpin")

# cache imgets - used throughout the script
cache("imgets","gemdate")

# set the local variables
l_inimages=inimages ; l_outimages=outimages ; l_skyimage=skyimage
l_skylevel=skylevel ; l_fl_sky=fl_sky 
l_fl_flat=fl_flat ; l_verbose=verbose ; l_prefix=outprefix
l_logfile=logfile ; l_fl_autosky=fl_autosky
l_filter1=filter1 ; l_filter2=filter2 ; l_filter3=filter3 ; l_filter4=filter4
l_filter5=filter5 ; l_filter6=filter6 ; l_keyfilter=key_filter
l_flatimage1=flatimage1 ; l_flatimage2=flatimage2 ; l_flatimage3=flatimage3 
l_flatimage4=flatimage4 ; l_flatimage5=flatimage5 ; l_flatimage6=flatimage6

# Check for package log file or user-defined log file
cache("quirc")
if((l_logfile=="") || (l_logfile==" ")) {
   l_logfile=quirc.logfile
   if((l_logfile=="") || (l_logfile==" ")) {
      l_logfile="quirc.log"
      printlog("WARNING - QREDUCE:  Both qreduce.logfile and quirc.logfile are empty.",logfile=l_logfile, verbose+)
      printlog("                 Using default file quirc.log.",
         logfile=l_logfile, verbose+)
   }
}
# Open log file
date | scan(l_struct)
printlog("----------------------------------------------------------------------------",logfile=l_logfile, verbose=l_verbose)
printlog("QREDUCE -- "//l_struct, logfile=l_logfile, verbose=l_verbose)
printlog(" ",logfile=l_logfile, verbose=l_verbose)

#-----------------------------------------------------------------------
# Check for consistent sky/sky level logic
if(l_fl_autosky) {
  if(l_skylevel != 0.0) {
     printlog("ERROR - QREDUCE:  You have specified both a sky constant to add to the",logfile=l_logfile,verbose+)
     printlog("                  final image AND set the fl_autosky flag to determine the",logfile=l_logfile,verbose+)
     printlog("                  sky constant from the sky image.  These options are",logfile=l_logfile,verbose+)
     printlog("                  incompatible.",logfile=l_logfile,verbose+)
     status=1
     goto clean
   } else if (!l_fl_sky) {
     printlog("ERROR - QREDUCE:  You have set the fl_autosky flag to determine the sky",logfile=l_logfile,verbose+)
     printlog("                  constant from the sky image, but have the fl_sky sky",logfile=l_logfile,verbose+)
     printlog("                  subtraction flag off.  These options are incompatible.",logfile=l_logfile,verbose+)
     status=1
     goto clean
   } 
} else {
  if((l_skylevel == 0.0) && (l_fl_sky)) {
     printlog("WARNING - QREDUCE:  You have set the sky constant to 0.0 and have the",logfile=l_logfile,verbose+)
     printlog("                    fl_autosky flag off.  No constant will be added",logfile=l_logfile,verbose+)
     printlog("                    after the sky image is subtracted.",logfile=l_logfile,verbose+)
  } else if ((l_skylevel != 0.0) && (!l_fl_sky)) {
     printlog("WARNING - QREDUCE:  You have specified a sky constant to add on, but",logfile=l_logfile,verbose+)
     printlog("                    the sky image subtraction flag fl_sky is off.  The",logfile=l_logfile,verbose+)
     printlog("                    sky constant "//l_skyconst//" will be added anyway!",logfile=l_logfile,verbose+)
  }
}

#-----------------------------------------------------------------------
# Load up arrays of input name lists

# Make list if * in inimages
if(stridx("*",l_inimages)>0) {
  files(l_inimages, > tmpin)
  l_inimages="@"//tmpin
}

#
nimages=0
if(substr(l_inimages,1,1)=="@") 
  scanfile=substr(l_inimages,2,strlen(l_inimages))
else {
 files(l_inimages,sort-, > tmpfile)
 scanfile=tmpfile
}

l_fl_dollar=no
while(fscan(scanfile,l_temp) != EOF) {
  if(!imaccess(l_temp))
     printlog("WARNING - QREDUCE: Input image "//l_temp//" not found.",logfile=l_logfile,verbose+)
  else {
  nimages=nimages+1
  if(nimages > maxfiles) {
     printlog("ERROR - QREDUCE: Maximum number of input images exceeded ("//str(maxfiles)//")",logfile=l_logfile,verbose+)
     status=1
     goto clean
     }
  in[nimages]=l_temp 
# Catch $  and / in input
    if(stridx("$",in[nimages])!=0)
      l_fl_dollar=yes
    if(stridx("/",in[nimages])!=0)
      l_fl_dollar=yes
  }
}
printlog("Processing "//nimages//" file(s).",logfile=l_logfile,verbose=l_verbose)
scanfile="" ; delete(tmpfile//","//tmpin,ver-, >& "dev$null")

# Now for the output images
# outimages could contain legal * if it is of a form like %st%stX%*.imh

noutimages=0
if(l_outimages!="" ) {
  if(substr(l_outimages,1,1)=="@") 
    scanfile=substr(l_outimages,2,strlen(l_outimages))
  else if (stridx("*",l_outimages)>0)  {
    files(l_outimages,sort-) | 
       match(".hhd",stop+,print-,metach-, > tmpfile)
    scanfile=tmpfile
  } else {
    files(l_outimages,sort-, > tmpfile)
    scanfile=tmpfile
  }

  while(fscan(scanfile,l_temp) != EOF) {
    noutimages=noutimages+1
    if(noutimages > maxfiles) {
       printlog("ERROR - QREDUCE: Maximum number of input images exceeded ("//str(maxfiles)//").",logfile=l_logfile,verbose+)
       status=1
       goto clean
     }
    out[noutimages]=l_temp 
    if(imaccess(out[noutimages])) {
      printlog("ERROR - QREDUCE: Output image "//out[noutimages]//" exists.",logfile=l_logfile,verbose+)
      status=1
      goto clean
    }
  }
}
scanfile="" ; delete(tmpfile,ver-, >& "dev$null")

# if there are too many or too few output images, and any defined
# at all at this stage - exit with error
if(nimages!=noutimages && l_outimages!="") {
  printlog("ERROR - QREDUCE: Number of input and output images mismatch.",logfile=l_logfile,verbose+)
  status=1
  goto clean
}

# If prefix is to be used instead
if(l_outimages=="" ) {
  if((l_prefix=="") || (l_prefix==" ")) {
    printlog("ERROR - QREDUCE: Neither output image name or output prefix is defined.",logfile=l_logfile,verbose+)
    status=1
    goto clean
  }
  if(l_fl_dollar) {
    printlog("ERROR - QREDUCE: Cannot use outprefix with path as part of input image names",
      l_logfile,verbose+)
    printlog("                 Set outimages to avoid this error",
      l_logfile,verbose+)
    status=1
    goto clean
  }
i=1
  while(i<=nimages) {
    out[i]=l_prefix//in[i]
    if(imaccess(out[i])) {
      printlog("ERROR - QREDUCE: Output image "//out[i]//" exists.",logfile=l_logfile,verbose+)
      status=1
      goto clean
    }
    i=i+1
  }
}

#-------------------------------------------------------------------------
# Check for existence of sky image, if needed

if(l_fl_sky) {
  if(!imaccess(l_skyimage) && l_skyimage!="" && stridx(" ",l_skyimage)<=0) {
    printlog("ERROR - QREDUCE: Sky image "//l_skyimage//" not found.",logfile=l_logfile,verbose+)
    status=1
    goto clean
  }
  else if (l_skyimage=="" || stridx(" ",l_skyimage)>0 ) {
   printlog("ERROR - QREDUCE: Sky image defined either as an empty string or contains spaces.",logfile=l_logfile,verbose+)
   status=1
   goto clean
 }
}

# Check if key_filter is set
if((l_keyfilter=="" || l_keyfilter==" ") && l_fl_flat) {
  printlog("WARNING - QREDUCE: key_filter not set. \
     Flat-fielding not performed.",l_logfile,verbose+)
  l_fl_flat=no
}

if(!l_fl_flat && !l_fl_sky) {
  printlog("ERROR - QREDUCE: No reduction steps selected.",l_logfile,verbose+)
  status=1
  goto clean
}

#--------------------------------------------------------------------------
# The math and bookkeeping:  (MAIN LOOP)

printlog(" ",logfile=l_logfile,verbose=l_verbose)
printlog("n  input --> output (sky,flat,flat filter,sky constant)",
   logfile=l_logfile,verbose=l_verbose)

l_fl_first=yes  # flag for some warnings
i=1
while(i<=nimages) {
  l_expression="im1"

#----------------
# check for previous sky subtraction and turn it off if necessary
  if(l_fl_sky) {
    imgets(in[i],"SKYIMAGE",>& "dev$null")
    if(imgets.value != "0") {
      l_fl_skytemp=no
      l_skyimage="none"
      printlog("WARNING - QREDUCE: Image "//in[i]//" has already been sky-subtracted.",logfile=l_logfile,verbose+)
      printlog("                   by qreduce.  Sky-subtraction not performed.",logfile=l_logfile,verbose+)
    }
    else {
      l_fl_skytemp=yes
      l_expression="("//l_expression//"-im2)"
# find the sky background level in the sky image, ignoring bad pixels
      if(l_fl_autosky) {
         imstat(in[i],fields="midpt,stddev",lower=INDEF,upper=INDEF,nclip=0,
	    lsigma=INDEF,usigma=INDEF,binwidth=0.01,format-,cache-) | \
	    scan(l_mid,l_sig)
         imstat(in[i],fields="midpt",lower=(l_mid-4*l_sig),
	    upper=(l_mid+4*l_sig),nclip=0,lsigma=INDEF,usigma=INDEF,
            binwidth=0.01,format-,cache-) | scan(l_mid)
         l_skylevel=l_mid
      }
    } # end else
  } # end if (l_fl_sky)
  else {
    l_skyimage="none"
    l_fl_skytemp=no
  }

#----------------
# check for previous flat fielding, and turn it off if necessary
l_filter="none" ; l_flatimage="none"  # for logging
  if(l_fl_flat) {
    imgets(in[i],"FLATIMAG",>& "dev$null")
    if(imgets.value != "0") {
      l_fl_flattemp=no
      l_flatimage="none"
      printlog("WARNING - QREDUCE: Image "//in[i]//" has already been flat-fielded.",logfile=l_logfile,verbose+)
      printlog("                   by qreduce.  Flat-fielding not performed.",logfile=l_logfile,verbose+)
      l_filter="none"  # set l_filter to something, for logging
    }
    else {
      l_fl_flattemp=yes
      imgets(in[i],l_keyfilter,>& "dev$null")
      l_filter=imgets.value
      if(i==1)
        l_temp=imgets.value

# Find the right flat field by comparing the filter name from the header
      flatok=no
      if(l_filter == filter1 && flatimage1!="") {
         l_flatimage=flatimage1 ; flatok=yes
         goto flatset
      }
      if(l_filter == filter2 && flatimage2!="") {
         l_flatimage=flatimage2 ; flatok=yes
         goto flatset
      }
      if(l_filter == filter3 && flatimage3!="") {
         l_flatimage=flatimage3 ; flatok=yes
         goto flatset
      }
      if(l_filter == filter4 && flatimage4!="") {
         l_flatimage=flatimage4 ; flatok=yes
         goto flatset
      }
      if(l_filter == filter5 && flatimage5!="") {
         l_flatimage=flatimage5 ; flatok=yes
         goto flatset
      }
      if(l_filter == filter6 && flatimage6!="") {
         l_flatimage=flatimage6 ; flatok=yes
         goto flatset
      }

flatset:
      if(!flatok) {
         printlog("WARNING - QREDUCE: Image "//in[i]//" is taken in filter "//l_filter,logfile=l_logfile,verbose+)
         printlog("                   Flat field not defined for this filter.",logfile=l_logfile,verbose+)
         printlog("                   Flat-fielding not performed.",logfile=l_logfile,verbose+)
         l_fl_flattemp = no ; l_flatimage="none"
      } else {  # START of FLATOK

        if(l_temp != l_filter && l_fl_sky && l_fl_first) {
          printlog("WARNING - QREDUCE: The list of input images contains data taken with ",logfile=l_logfile,verbose+)
          printlog("                   different filters!  If sky-subtracting, the same",logfile=l_logfile,verbose+)
          printlog("                   sky frame will be subtracted from ALL the input",logfile=l_logfile,verbose+)
          printlog("                   data.  Probably a bad idea!",logfile=l_logfile,verbose+)       
          l_fl_first=no
        }
  
        if((i==1) || (l_temp != l_filter)) {
          if(!imaccess(l_flatimage)) {
            printlog("ERROR - QREDUCE: Flat field "//l_flatimage//" not found.",logfile=l_logfile,verbose+)
            status=1
            goto clean
          }
        }
      l_temp=l_filter
      l_expression=l_expression//"/im3"
      }
    }
  }
  else {   # NO flat fielding
    l_fl_flattemp=no ; l_flatimage="none"
  }

#-------------
# add sky constant to expression
l_expression = l_expression//"+"//str(l_skylevel)

#-------------

printlog(i//"  "//in[i]//" --> "//out[i]//" ("//l_skyimage//", "//l_flatimage//", "//l_filter//", "//str(l_skylevel)//")",logfile=l_logfile,verbose=l_verbose)

# Do the math, finally!
    imcalc(in[i]//","//l_skyimage//","//l_flatimage,out[i],l_expression,
      pixtype="real",verbose=no)

# and update the header
    gemdate ()
# date stamp the modification
    gemhedit(out[i],"GEM-TLM",gemdate.outdate,"UT Last modification with GEMINI")
    gemhedit(out[i],"QREDUCE",gemdate.outdate,"UT Time stamp for qreduce")
    if(l_fl_skytemp)
      gemhedit(out[i],"SKYIMAGE",l_skyimage,"Sky image used by qreduce")
    if(l_skylevel!=0.0)
      gemhedit(out[i],"SKYCONST",l_skylevel,"Constant sky level from qreduce")
    if(l_fl_flattemp)
      gemhedit(out[i],"FLATIMAG",l_flatimage,"Flat image used by qreduce")


  i=i+1
}
# end the main loop

#---------------------------------------------------------------------------
# Clean up
clean:
{
  if(status==0) {
    printlog("QREDUCE exit status:  good.",logfile=l_logfile, verbose=l_verbose)
  }
   printlog("----------------------------------------------------------------------------", logfile=l_logfile, verbose=l_verbose)
  scanfile=""
  delete(tmpin,ver-, >& "dev$null")
}
end
