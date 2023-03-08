# Copyright(c) 2001-2012 Association of Universities for Research in Astronomy, Inc.

procedure obackground(image)

# Compute median statistics of reference frames in OSCIR
#   image; print statistics and plot median values as
#   % of full well.  Identify frames outside of
#   user-specified (default=3) sigma of average median.
# Output stats to logfile, plot to .ps file, and
#   list of "Bad" frames to text file
# NOTES: works on single image only
#
# Version: Sept 14, 2002 BR  Release v1.4
#          Aug 20, 2003  KL  IRAF2.12 - new parameters
#                              imstat: nclip, lsigma, usigma, cache

char image 	{prompt="Input OSCIR Image"}
real sigma	{3.0,prompt="Sigma tolerance for bad frames"}
char logfile 	{"",prompt="Logfile"}
char bsetfile	{"",prompt="Bad Frame list file"}
bool verbose	{yes,prompt="Verbose?"}
bool fl_writeps {yes,prompt="Write .ps file?"}
#int  status	{0,prompt="Exit status (0=good)"}

begin

char l_image,l_logfile,l_bsetfile
real l_sigma
bool l_verbose, l_fl_writeps
char tmpstata,tmpstatb, tmpmeana, tmpmeanb, tmpmeans, tmpps, tmpbsets
int naxis3, nsaveset, naxis5, nnodset
int i, j
int frmcoadd, chpcoadd, ncoadd, nref
real ADC_DARK,ADC_SAT
real pymin, pymax
struct l_struct

cache("imgets")

l_image=image
l_sigma=sigma
l_verbose=verbose
l_fl_writeps=fl_writeps
l_logfile=logfile
l_bsetfile=bsetfile

# Check for package log file or user-defined log file
cache("oscir")
if((l_logfile=="") || (l_logfile==" ")) {
   l_logfile=oscir.logfile
   if((l_logfile=="") || (l_logfile==" ")) {
#if logfile not defined, default to oscir.log
      l_logfile="oscir.log"
}
}
printf("OBACKGROUND:  Using logfile %s\n",l_logfile)

# Open log file
date | scan(l_struct)
printlog("----------------------------------------------------------------------------",l_logfile,l_verbose)
printlog("OBACKGROUND -- "//l_struct,l_logfile,l_verbose)
printlog(" ",l_logfile,l_verbose)

#check if image exists
if(!imaccess(l_image)) {
  printlog("ERROR - OBACKGROUND: Image "//l_image//" not found",l_logfile,verbose+)
  bye
}

#get axes
imgets(l_image,"i_naxis3") ; naxis3=int(imgets.value)
imgets(l_image,"i_naxis4") ; nsaveset=int(imgets.value)
imgets(l_image,"i_naxis5") ; naxis5=int(imgets.value)
imgets(l_image,"i_naxis6") ; nnodset=int(imgets.value)
#check if image contains chop-nod data
if ((naxis3!=2)||(naxis5!=2)) {
   printlog("ERROR - OBACKGROUND: Image "//l_image//" is not chop-nod data;",l_logfile,verbose+)
   printlog("                     n_choppos= "//naxis3//", n_nodpos= "//naxis5,l_logfile,l_verbose)
   bye
}

imgets(l_image,"frmcoadd") ; frmcoadd=int(imgets.value)
imgets(l_image,"chpcoadd") ; chpcoadd=int(imgets.value)
ncoadd=frmcoadd*chpcoadd
if (ncoadd <= 0) {
	printlog("ERROR - OBACKGROUND: Header info missing for FRMCOADD or CHPCOADD",l_logfile,verbose+)
	bye
}

#create tmpfiles
tmpstata=mktemp("tmpstat"); tmpstatb=mktemp("tmpstat")
tmpmeana=mktemp("tmpmean"); tmpmeanb=mktemp("tmpmean")
tmpmeans =mktemp("tmpmeans")//".fits"
tmpbsets =mktemp("tmpbsets")//".fits"

#define well parameters for OSCIR (from f6bstat)
ADC_DARK=1.315e4
ADC_SAT=5.793e4

printf("%-25s %-3s %-3s %-3s %-3s\n","Image","n_choppos","n_savesets","n_nodpos","n_nodsets") | scan(l_struct)
printlog(l_struct,l_logfile,l_verbose)
printf("%-25s   %-9d %-10d %-8d %-9d\n",l_image,naxis3,nsaveset,naxis5,nnodset) | scan(l_struct)
printlog(l_struct,l_logfile,l_verbose)

printf("\n     Please wait while I compute frame statistics...\n\n")

# print to file all ref frames for each save & nod set: 
# 	ref1=[chop2,nod1]; ref2=[chop1,nod2]
# Save nod1 and nod2 to seperate lists
#NOTE: f6bstat lists all ref1s then all ref2s per nod set (affects plot)
for(i=1;i<=nnodset;i+=1) {
for(j=1;j<=nsaveset;j+=1) {
  print(l_image//"[*,*,2,"//str(j)//",1,"//str(i)//"]", >> tmpstata)
  print(l_image//"[*,*,1,"//str(j)//",2,"//str(i)//"]", >> tmpstatb)
}
}
# stats do not exactly match f6bstat because f6bstat also demeans(?) result (subtracting 0th moment I think)
#get means for each nod seperately for plotting
imstat("@"//tmpstata,fields="midpt",lower=INDEF,upper=INDEF,nclip=0,
   lsigma=INDEF,usigma=INDEF,binwidth=0.1,format-,cache-, > tmpmeana)
imstat("@"//tmpstatb,fields="midpt",lower=INDEF,upper=INDEF,nclip=0,
   lsigma=INDEF,usigma=INDEF,binwidth=0.1,format-,cache-, > tmpmeanb)

#change to %well: 
tcalc(tmpmeana,"Well",
"(c1/"//str(ncoadd)//"-"//str(ADC_DARK)//")/"//str(ADC_SAT-ADC_DARK),colfmt="f8.4")
tcalc(tmpmeana,"Well","Well*100.",colfmt="f6.2")
tcalc(tmpmeana,"Row","ROWNUM",colfmt="i5")
tcalc(tmpmeanb,"Well",
"(c1/"//str(ncoadd)//"-"//str(ADC_DARK)//")/"//str(ADC_SAT-ADC_DARK),colfmt="f8.4")
tcalc(tmpmeanb,"Well","Well*100.",colfmt="f6.2")
tcalc(tmpmeanb,"Row","ROWNUM",colfmt="i5")
#join tables for statistics: final table = ref1(nodsets A) + ref2(nodsets B)
tmerge(tmpmeana//","//tmpmeanb,tmpmeans,"append", tbltype="default")
tstat(tmpmeans,"Well",outtable="STDOUT", >>& "dev$null")
printf("Signal [percent full well] in reference frames:\n") | scan(l_struct)
printlog(l_struct,l_logfile,l_verbose)
printf("                        Average  = %6.2f\n",tstat.mean) | scan(l_struct)
printlog(l_struct,l_logfile,l_verbose)
printf("                        Stddev   = %8.4f\n",tstat.stddev) | scan(l_struct)
printlog(l_struct,l_logfile,l_verbose)
printf("                        Minimum  = %6.2f\n",tstat.vmin) | scan(l_struct)
printlog(l_struct,l_logfile,l_verbose)
printf("                        Maximum  = %6.2f\n",tstat.vmax) | scan(l_struct)
printlog(l_struct,l_logfile,l_verbose)

#set x axis max value
nref=nsaveset*nnodset
#set y axis values for plot to min/max or +/- 2 percentage points
if (tstat.mean-2.0 < tstat.vmin) {
	pymin=tstat.mean - 2.0
} else {
	pymin=tstat.vmin - tstat.stddev
}
if (tstat.mean+2.0 > tstat.vmax) {
	pymax=tstat.mean + 2.0
} else {
	pymax=tstat.vmax + tstat.stddev
}

#plot results

#get rid of "_" in image title
printf(" %s Background\n",l_image) | translit("STDIN","_"," ") | scan(l_struct)

if (l_fl_writeps){
#plot to file first (takes awhile)
  mgograph(tmpmeana,3,2,
   rows="-",wx1=0,wx2=nref,wy1=pymin,wy2=pymax,excol="",eycol="",logx=no,logy=no,
   labelexp=1.5,boxexp=1.,xlabel="Reference Frame Number",ylabel="% of Full Well",
   title=l_struct, postitle="topcenter",append=no,
   pointmode=no,pattern="solid",crvstyle="straight",lweight=1,color=1,
   mkzero=no,device="psi_land",gkifile="mgo.gki")
   mgograph(tmpmeanb,3,2,
    append=yes,title="",pointmode=no,pattern="solid",crvstyle="straight",lweight=2,
    color=4,mkzero=no,device="psi_land",gkifile="mgo.gki")
  #add legend (position units are for normalized coords: 0-1 in x,y
  igi(initcmd="DRELOCATE 0.35 0.85 LWEIGHT 1 COLOR 1 DDRAW 0.42 0.85 PUTLABEL 6 NodA ;DRELOCATE 0.6 0.85 LWEIGHT 2 COLOR 4 DDRAW 0.67 0.85 PUTLABEL 6 NodB; END",append=yes,device="psi_land", >>&"dev$null")
  gflush
}

mgograph(tmpmeana,3,2,
  rows="-",wx1=0,wx2=nref,wy1=pymin,wy2=pymax,excol="",eycol="",logx=no,logy=no,
  labelexp=1.5,boxexp=1.,xlabel="Reference Frame Number",ylabel="% of Full Well",
  title=l_struct, postitle="topcenter",append=no,
  pointmode=no,pattern="solid",crvstyle="straight",lweight=1,color=1,
  mkzero=no,device="stdgraph",gkifile="mgo.gki")
mgograph(tmpmeanb,3,2,
  append=yes,title="",pointmode=no,pattern="solid",crvstyle="straight",lweight=2,color=4,
  mkzero=no,device="stdgraph",gkifile="mgo.gki")
#add legend
#igi(initcmd="DRELOCATE 0.35 0.85 LWEIGHT 1 COLOR 1 DDRAW 0.42 0.85 PUTLABEL 6 NodA ;DRELOCATE 0.6 0.85 LWEIGHT 2 COLOR 4 DDRAW 0.67 0.85 PUTLABEL 6 NodB; END",append=yes,device="stdgraph", >>&"dev$null")

#Finding bad savesets: median values l_sigma(default 4) below/above average
if ((tstat.vmin < (tstat.mean-tstat.stddev*l_sigma)) || (tstat.vmax > (tstat.mean+tstat.stddev*l_sigma))) {
   printlog("WARNING - OBACKGROUND: Bad frames (median outside "//l_sigma//" sigma) exist.\n",l_logfile,verbose+)
   printf("\n    Looking for bad frames...\n")
   if((l_bsetfile=="") || (l_bsetfile==" ")) {
#if bad frame file not defined, default to image_name.bsets
      l_bsetfile=l_image//".bsets"
   }
#check if file exists
if (access(l_bsetfile)) {
	printlog("WARNING - OBACKGROUND: Appending to existing Bad frames file: "//l_bsetfile,l_logfile,verbose+)
#	delete(l_bsetfile,verify-, >>&"dev$null")
  }
   printf("\n File: "//l_image//"  -  Save and Nod sets outside "//l_sigma//" sigma of mean\n", >>l_bsetfile)
   tselect(tmpmeans,tmpbsets,"Well > ("//str(tstat.mean)//"+"//str(tstat.stddev)//"*"//str(l_sigma)//") || Well < ("//str(tstat.mean)//"-"//str(tstat.stddev)//"*"//str(l_sigma)//")")
#NOTE: calculation of nodset saveset assumes all ref1's followed by ref2's:
   tcalc(tmpbsets,"Nodset","int((Row-1)/"//str(nsaveset)//")+1",datatype="int",colfmt="i5")
   tcalc(tmpbsets,"Saveset","Row-((Nodset-1)*"//str(nsaveset)//")",datatype="int",colfmt="i5")   
   tprint(tmpbsets,prdata+,showrow-,showhdr+,align+,>>l_bsetfile)
   printlog(" Savesets outside "//l_sigma//" sigma of mean written to "//l_bsetfile, l_logfile,l_verbose) 
} else {
   printlog("  No frames with median outside "//l_sigma//" sigma of mean.\n",l_logfile,l_verbose)
}

#delete temporary files
delete(tmpstata//","//tmpstatb,verify-, >>&"dev$null")
delete(tmpmeana//","//tmpmeanb//","//tmpmeans,verify-, >>&"dev$null")
if (access(tmpbsets)) {
	delete(tmpbsets,verify-, >>&"dev$null")
}
if (l_fl_writeps) {
  # Move the PS image from temporary file
  tmpps=""
  files("ps*.eps",sort-) | scan(tmpps)
  if(tmpps!="") {
    if (access(l_image//"_ref.ps")) {
	printlog("WARNING - OBACKGROUND: Overwriting previous .ps output image",l_logfile,verbose+)
	delete(l_image//"_ref.ps",verify-, >>&"dev$null")
    }
    rename(tmpps,l_image//"_ref.ps", field="all")
    printlog("Postscript file of median sky level in ref frames: "//l_image//"_ref.ps",
      l_logfile,l_verbose)
  } else {
    printlog("ERROR - OBACKGROUND: Cannot find temporary .ps output image",
      l_logfile,verbose+)
  }
}

end








