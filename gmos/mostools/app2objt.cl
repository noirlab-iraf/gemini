# Copyright(c) 2001-2017 Association of Universities for Research in Astronomy , Inc.

procedure app2objt(intable)

# Make Object Table from apphot or daophot output (ascii)
# Needs image for WCS info
#
# Version Nov 26, 2001 IJ  stand-alone release (logfile handling modified in released version)
#         Feb 28, 2002 IJ  v1.3 release
#         Sept 20, 2002 IJ v1.4 release
#         Sept 18, 2003 KL Bug fix: missing closing bracket at pselect and 
#                              tproject calls.
#   2017-04-21 (mischa): Commented out the part with CCDSUM and thus
#                        the necessity for a MEF

char      intable   {prompt="APPHOT/DAOPHOT output to be converted"}
char      image     {"default",prompt="Image containing the WCS info"}
char      outtable  {"default",prompt="Output FITS Object Table"}
char      priority  {"2",prompt="Default priority",enum="1|2|3|X"}
char      logfile   {"",prompt="Logfile"}
bool      verbose   {no,prompt="Verbose?"}
int       status    {0,prompt="Exit status (0=good)"}
char      *scanfile {"",prompt="For internal use only"}

begin

char l_intable, l_image, l_outtable, l_logfile, l_priority
bool l_verbose

char l_rootname, s_empty, tmpdat, tmppri, tmpout, l_task
struct l_struct
int l_ii
bool l_tryadd
real crpix1, crpix2, crval1, crval2, cd1_1, cd1_2, cd2_1, cd2_2

l_intable=intable ; l_image=image ; l_outtable=outtable
l_logfile=logfile ; l_verbose=verbose ; l_priority=priority
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
     printlog("WARNING - APP2OBJT: both app2objt.logfile and \
      mostools.logfile are empty.",l_logfile,l_verbose)
     printlog("                Using default file gmos.log.",l_logfile,l_verbose)
  }
}

# Write to the logfile
date | scan(l_struct)
printlog("-----------------------------------------------------------------\
-----------",l_logfile,l_verbose)
printlog("APP2OBJT -- "//l_struct,l_logfile,l_verbose)
printlog("",l_logfile,l_verbose)


if(!access(l_intable)) {
  printlog("ERROR - APP2OBJT: Input table "//l_intable//" not found",l_logfile,yes)
  goto crash
}

l_ii=stridx(".",l_intable)
if(l_ii>0)
  l_rootname=substr(l_intable,1,l_ii-1)
else
  l_rootname=l_intable

if(l_image=="" || l_image=="default") {
 l_image=l_rootname
}
l_tryadd=no
gimverify(l_image)
if(gimverify.status==1) {
    printlog("WARNING - APP2OBJT: Image "//l_image//" not found.",l_logfile,yes)
    l_tryadd=yes
  }
#	else if(gimverify.status>1) {
#    printlog("WARNING - APP2OBJT: Image "//l_image//" not a MEF FITS image.",l_logfile,yes)
#    l_tryadd=yes
#} 
if(l_tryadd) {
   printlog("                    Trying "//l_image//"_add instead",l_logfile,yes)
   gimverify(l_image//"_add")
   if(gimverify.status==1) {
       printlog("WARNING - APP2OBJT: Image "//l_image//"_add not found.",l_logfile,yes)
       goto crash
   }
	#	else if(gimverify.status>1) {
    #    printlog("WARNING - APP2OBJT: Image "//l_image//"_add not a MEF FITS image.",l_logfile,yes)
    #   goto crash
    # } 
}
l_image=gimverify.outname

if(l_outtable=="" || l_outtable=="default")
  l_outtable=l_rootname//"_OT"
if(access(l_outtable//".fits")) {
    printlog("ERROR - APP2OBJT: Output table "//l_outtable//".fits exits",
      l_logfile,yes)
    goto crash
}

# Check format of input tabl
l_task="DUMMY"
match("TASK",l_intable,stop-,print_file_n-) | scan(l_task,l_task,l_task,l_task)
if(l_task!="daofind" && l_task!="phot" && l_task!="nstar" && l_task!="allstar") {
   printlog("ERROR - APP2OBJT: Input table "//l_intable//" not from a supported task",
     l_logfile,yes)
   if(l_task!="DUMMY")
     printlog("Input table "//l_intable//" is from task "//l_task,l_logfile,yes)
   printlog("Supported tasks are daofind, phot, allstar and nstar",
     l_logfile,yes)
   goto crash
}


#################################################################
# Convert the table
printlog("Input table : "//l_intable,l_logfile,l_verbose)
printlog("Image       : "//l_image,l_logfile,l_verbose)
printlog("Output table: "//l_outtable,l_logfile,l_verbose)

# Get rid of INDEF magnitudes
pselect(l_intable,tmpout,"mag<10000.")
# convert directly to FITS
pconvert(tmpout,l_outtable//".fits","*",expr="yes",append-)
delete(tmpout,verify-)
# Fix non-compliant column names, and column names that create conflicts
tchcol(l_outtable//".fits","xcenter","x_ccd","","",verbose-)
tchcol(l_outtable//".fits","ycenter","y_ccd","","",verbose-)
tchcol(l_outtable//".fits","RAPERT","aperture","","",verbose-, >>& "dev$null")

# Get the WCS info
imgets(l_image//"[0]","CRPIX1") ; crpix1=real(imgets.value)
imgets(l_image//"[0]","CRPIX2") ; crpix2=real(imgets.value)
imgets(l_image//"[0]","CRVAL1") ; crval1=real(imgets.value)
imgets(l_image//"[0]","CRVAL2") ; crval2=real(imgets.value)
imgets(l_image//"[0]","CD1_1") ; cd1_1=real(imgets.value)
imgets(l_image//"[0]","CD1_2") ; cd1_2=real(imgets.value)
imgets(l_image//"[0]","CD2_1") ; cd2_1=real(imgets.value)
imgets(l_image//"[0]","CD2_2") ; cd2_2=real(imgets.value)

###CRPIX1  =      1501.5550450756 / Ref pix of axis 1
###CRPIX2  =     1145.01924489263 / Ref pix of axis 2
###CRVAL1  =     348.702907594812 / RA at Ref pix in decimal degrees
###CRVAL2  =     2.03074277696353 / DEC at Ref pix in decimal degrees
###CD1_1   = -3.09827815448039E-05 / WCS matrix element 1 1
###CD1_2   = -2.59680639792843E-05 / WCS matrix element 1 2
###CD2_1   = -2.59350194118771E-05 / WCS matrix element 2 1
###CD2_2   = 3.10042969173884E-05 / WCS matrix element 2 2
##
tcalc(l_outtable//".fits","DEC",
"(x_ccd-"//str(crpix1)//")*"//str(cd2_1)//"+(y_ccd-"//str(crpix2)//")*"//str(cd2_2)//"+"//str(crval2),
colfmt="%12.2h",colunit="deg")

tcalc(l_outtable//".fits","RA",
"((x_ccd-"//str(crpix1)//")*"//str(cd1_1)//"+(y_ccd-"//str(crpix2)//")*"//str(cd1_2)//")/cos(DEC/57.2956)+"//str(crval1),
colfmt="%12.2h",colunit="H")

tcalc(l_outtable//".fits","RA","RA/15.")

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
parkey("apphot/daophot",l_outtable//".fits","EXTSOFT",add+)
parkey("Gemini IRAF package gmos "//gemini.verno,l_outtable//".fits","GEMSOFT",add+)
imgets(l_image//"[0]","INSTRUME")
parkey(imgets.value,l_outtable//".fits","INSTRUME",add+)
# imgets(l_image//"[1]","CCDSUM")
# print(imgets.value) | scan(l_ii)
# parkey(str(0.0727*l_ii),l_outtable//".fits","PIXSCALE",add+)

# Set the default priority
tinfo(l_outtable//".fits",ttout-)
for(l_ii=1;l_ii<=tinfo.nrows;l_ii+=1) {
  print(l_priority, >> tmpdat)
}
print("priority ch*1 %1d") | \
tcreate(tmppri//".fits","STDIN",tmpdat,hist-,tbltype="default")
tmerge(l_outtable//".fits,"//tmppri//".fits",tmpout//".fits","merge",\
    tbltype="default")
delete(l_outtable//".fits",verify-)
# Project essential columns - get rid of unneeded string columns,
# and columns containing []. The [] makes the mask making software crash
tproject(tmpout//".fits",l_outtable//".fits", \
    col="x_ccd,y_ccd,mag,id,ra,dec,priority,merr,msky,stdev,aperture,"//\
    "sharpness,chi")
delete(tmpdat//","//tmppri//".fits,"//tmpout//".fits",verify-)

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
