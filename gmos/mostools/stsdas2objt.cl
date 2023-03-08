# Copyright(c) 2001-2012 Association of Universities for Research in Astronomy , Inc.
procedure stsdas2objt(intable)

# Make Object Table from STSDAS or FITS table
# Optional to get WCS info from image 
#
# Version Nov 7, 2001  IJ  First draft, minimal checking of input
#         Nov 26, 2001 IJ  stand-alone release (logfile handling modified in released version)
#
#         Feb 28, 2002 IJ  v1.3 release
#         Jul 23, 2002 ML - some fixes for Flamingos and improved handled of existing priorities
#         Aug 14, 2002 IJ bug fix for handling of priority column
#         Sept 20, 2002 IJ v1.4 release
#         Dec 4, 2002  IJ  fixed priority bug

char      intable   {prompt="STSDAS or FITS table to be converted"}
char      image     {"default",prompt="Image to be used for mask design"}
bool      fl_wcs    {no,prompt="Calculate RA and DEC from image WCS"}
char      outtable  {"default",prompt="Output FITS Object Table"}
char      instrument {"gmos",enum="gmos|flamingos",prompt="Instrument (gmos|flamingos)"} 
char      id_col    {"id",prompt="Input id number column (contains integer number)"}
char      mag_col   {"mag",prompt="Input magnitude column"}
char      x_col     {"x_ccd",prompt="Input column for X coordinate of position (pixels)"}
char      y_col     {"y_ccd",prompt="Input column for Y coordinate of position (pixels)"}
char      ra_col    {"RA",prompt="Input column for RA (hours)"}
char      dec_col   {"DEC",prompt="Input column for DEC (degrees)"}
char      pri_col   {"priority",prompt="Input column for priorities"} 
char      other_col {"",prompt="Comma separated list of other columns to include"}
char      priority  {"2",prompt="Default priority",enum="1|2|3|X|none"}
char      logfile   {"",prompt="Logfile"}
bool      verbose   {no,prompt="Verbose?"}
int       status    {0,prompt="Exit status (0=good)"}
char      *scanfile {"",prompt="For internal use only"}

begin

char l_intable, l_image, l_outtable, l_logfile, l_priority
char l_idcol, l_magcol, l_xcol, l_ycol, l_racol, l_deccol, l_othercol, l_allcol
bool l_verbose, l_fl_wcs
char l_instrument, l_pri_col

char l_rootname, s_empty, tmpdat, tmppri, tmpout, l_task, l_suf
struct l_struct
int l_ii, l_check
bool l_tryadd
real crpix1, crpix2, crval1, crval2, cd1_1, cd1_2, cd2_1, cd2_2

l_intable=intable ; l_image=image ; l_outtable=outtable ; l_instrument=instrument
l_logfile=logfile ; l_verbose=verbose ; l_priority=priority ; l_fl_wcs=fl_wcs
l_idcol=id_col ; l_magcol=mag_col ; l_xcol=x_col ; l_ycol=y_col 
l_racol=ra_col ; l_deccol=dec_col ; l_othercol=other_col ; l_pri_col=pri_col
l_allcol=l_idcol//","//l_racol//","//l_deccol//","//l_magcol//","//l_xcol//","//l_ycol//","//l_othercol//","//l_pri_col
status=0

tmpdat=mktemp("tmpdat")
tmppri=mktemp("tmppri")
tmpout=mktemp("tmpout")

cache("gimverify","imgets","parkey","tinfo")

# Define the name of the logfile
s_empty=""; print(l_logfile) | scan(s_empty); l_logfile=s_empty
if (l_logfile == "") {
  l_logfile = mostools.logfile
  if (l_logfile == "") {
     l_logfile = "gmos.log"
     printlog("WARNING - STSDAS2OBJT: Both stsdas2objt.logfile and \
      mostools.logfile are empty.",l_logfile,l_verbose)
     printlog("                Using default file gmos.log",l_logfile,l_verbose)
  }
}

# Write to the logfile
date | scan(l_struct)
printlog("-----------------------------------------------------------------\
-----------",l_logfile,l_verbose)
printlog("STSDAS2OBJT -- "//l_struct,l_logfile,l_verbose)
printlog("",l_logfile,l_verbose)

# Figure out if the table is STSDAS or FITS
if(substr(l_intable,strlen(l_intable)-4,strlen(l_intable))==".fits") {
  l_suf=".fits"
  l_intable=substr(l_intable,1,strlen(l_intable)-5)
} else if(substr(l_intable,strlen(l_intable)-3,strlen(l_intable))==".tab") {
  l_suf=".tab"
  l_intable=substr(l_intable,1,strlen(l_intable)-4)
} else if (access(l_intable//".tab"))
  l_suf=".tab"
else if (access(l_intable//".fits"))
  l_suf=".fits"
else {
  printlog("ERROR - STSDAS2OBJT: Input table "//l_intable//" not found",l_logfile,yes)
  goto crash
}

l_rootname=l_intable

if(l_image=="" || l_image=="default") {
 l_image=l_rootname
}
l_tryadd=no
if(l_image!=l_intable) {
  gimverify(l_image)
  if(gimverify.status==1) {
      printlog("WARNING - STSDAS2OBJT: Image "//l_image//" not found.",l_logfile,yes)
      l_tryadd=yes
  } else if(gimverify.status>1) {
      printlog("WARNING - STSDAS2OBJT: Image "//l_image//" not a MEF FITS image.",l_logfile,yes)
      l_tryadd=yes
  } 
} else
l_tryadd=yes

if(l_tryadd) {
   printlog("WARNING - STSDAS2OBJT: Trying "//l_image//"_add",l_logfile,yes)
   gimverify(l_image//"_add")
   if(gimverify.status==1) {
       printlog("WARNING - STSDAS2OBJT: Image "//l_image//"_add not found.",l_logfile,yes)
       goto crash
   } else if(gimverify.status>1) {
        printlog("WARNING - STSDAS2OBJT: Image "//l_image//"_add not a MEF FITS image.",l_logfile,yes)
       goto crash
   } 
}
l_image=gimverify.outname

if(l_outtable=="" || l_outtable=="default")
  l_outtable=l_rootname//"_OT"
if(access(l_outtable//".fits")) {
    printlog("ERROR - STSDAS2OBJT: Output table "//l_outtable//".fits exits",
      l_logfile,yes)
    goto crash
}


#################################################################
# Convert the table
printlog("Input table : "//l_intable,l_logfile,l_verbose)
printlog("Image       : "//l_image,l_logfile,l_verbose)
printlog("Output table: "//l_outtable,l_logfile,l_verbose)

# tproject columns - write to STSDAS since the next column
# manipulation does not work in FITS
tproject(l_intable//l_suf,l_outtable//".fits",col=l_allcol)

# Make copies of any non-compliant column names
if(l_xcol!="x_ccd")
  tchcol(l_outtable//".fits",l_xcol,"x_ccd","","pixels",ver-)
#  tcalc(l_outtable//".fits","x_ccd",l_xcol,colfmt="f6.2",colunit="pixel")
if(l_ycol!="y_ccd")
  tchcol(l_outtable//".fits",l_ycol,"y_ccd","","pixels",ver-)
#  tcalc(l_outtable//".fits","y_ccd",l_ycol,colfmt="f6.2",colunit="pixel")
if(l_idcol!="ID")
  tchcol(l_outtable//".fits",l_idcol,"ID","%9d","##",ver-)
#  tcalc(l_outtable//".fits","ID","int("//l_idcol//")",colfmt="i5")
if(l_magcol!="MAG")
  tchcol(l_outtable//".fits",l_magcol,"MAG","","magnitudes",ver-)
#  tcalc(l_outtable//".fits","MAG",l_magcol,colfmt="f6.2",colunit="mag")
if(l_racol!="RA")
  tchcol(l_outtable//".fits",l_racol,"RA","","H",ver-)
#  tcalc(l_outtable//".fits","RA",l_racol,colfmt="%16.6f",colunit="H")
if(l_deccol!="DEC")
  tchcol(l_outtable//".fits",l_deccol,"DEC","","deg",ver-)
#  tcalc(l_outtable//".fits","DEC",l_deccol,colfmt="%16.6f",colunit="deg")
if(l_pri_col!="priority" && l_pri_col!="")
  tchcol(l_outtable//".fits",l_pri_col,"priority","","",ver-)

if(l_fl_wcs) {
# Get the WCS info
imgets(l_image//"[0]","CRPIX1") ; crpix1=real(imgets.value)
imgets(l_image//"[0]","CRPIX2") ; crpix2=real(imgets.value)
imgets(l_image//"[0]","CRVAL1") ; crval1=real(imgets.value)
imgets(l_image//"[0]","CRVAL2") ; crval2=real(imgets.value)
imgets(l_image//"[0]","CD1_1") ; cd1_1=real(imgets.value)
imgets(l_image//"[0]","CD1_2") ; cd1_2=real(imgets.value)
imgets(l_image//"[0]","CD2_1") ; cd2_1=real(imgets.value)
imgets(l_image//"[0]","CD2_2") ; cd2_2=real(imgets.value)


tcalc(l_outtable//".fits","DEC",
"(x_ccd-"//str(crpix1)//")*"//str(cd2_1)//"+(y_ccd-"//str(crpix2)//")*"//str(cd2_2)//"+"//str(crval2),
colfmt="%12.2h",colunit="deg")

tcalc(l_outtable//".fits","RA",
"((x_ccd-"//str(crpix1)//")*"//str(cd1_1)//"+(y_ccd-"//str(crpix2)//")*"//str(cd1_2)//")/cos(DEC/57.2956)+"//str(crval1),
colfmt="%12.2h",colunit="H")

tcalc(l_outtable//".fits","RA","RA/15.")

} # end of WCS from image

# Header info for the table - get most of this from the input image
imgets(l_image//"[0]","GEMPRGID")
parkey(imgets.value,l_outtable//".fits","PID_IMAG",add+)
parkey(imgets.value,l_outtable//".fits","PID_SPEC",add+)
imgets(l_image//"[0]","DATE-OBS")
parkey(imgets.value,l_outtable//".fits","DATEIMAG",add+)
imgets(l_image//"[0]","TIME-OBS")
parkey("x",l_outtable//".fits","TIMEIMAG",add+)
parkey(imgets.value,l_outtable//".fits","TIMEIMAG",add+)
imgets(l_image//"[0]","RA")
parkey(real(imgets.value)/15.,l_outtable//".fits","RA_IMAG",add+)
imgets(l_image//"[0]","DEC")
parkey(imgets.value,l_outtable//".fits","DEC_IMAG",add+)
imgets(l_image//"[0]","FILTER2")
parkey(imgets.value,l_outtable//".fits","FILTER",add+)
imgets(l_image//"[0]","ADCUSED")
parkey(imgets.value,l_outtable//".fits","ADCMODE",add+)
time | scan(l_struct)
parkey(l_struct,l_outtable//".fits","DATE_OT",add+)
parkey(l_struct,l_outtable//".fits","TIME_OT",add+)
parkey("FITS/STSDAS table",l_outtable//".fits","EXTSOFT",add+)
parkey("Gemini IRAF package gmos "//gemini.verno,l_outtable//".fits","GEMSOFT",add+)
imgets(l_image//"[0]","INSTRUME")
parkey(imgets.value,l_outtable//".fits","INSTRUME",add+)
imgets(l_image//"[1]","CCDSUM")
print(imgets.value) | scan(l_ii)
if(l_instrument == "gmos") { 
   parkey(str(0.0727*l_ii),l_outtable//".fits","PIXSCALE",add+)
} else if (l_instrument == "flamingos") { 
  parkey(0.078,l_outtable//".fits","PIXSCALE",add+)
}   

# Set the default priority if needed, but make sure the priority column does
# not already exist - could be in the user's other_col as in the old version
# of this script
if(l_priority!="none" && l_pri_col=="") {
  l_check=0
  tlcol(l_outtable//".fits") | match("priority","STDIN",stop-) | \
        count("STDIN") | scan(l_check)
  if(l_check==0) {
    tinfo(l_outtable//".fits",ttout-)
    for(l_ii=1;l_ii<=tinfo.nrows;l_ii+=1) {
      print(l_priority, >> tmpdat)
    }
    print("priority ch*1 %1d") | \
    tcreate(tmppri//".fits","STDIN",tmpdat,hist-,tbltype="default")
    tmerge(l_outtable//".fits,"//tmppri//".fits",tmpout//".fits","merge")
    delete(l_outtable//".fits",verify-)
    rename(tmpout//".fits",l_outtable//".fits",field="all")
    delete(tmppri//".fits,"//tmpdat,verify-)
  }
}
  
goto clean

crash:
status=1
goto theend

clean:
status=0

theend:
printlog("-----------------------------------------------------------------\
-----------",l_logfile,l_verbose)

end
