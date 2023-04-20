# Copyright(c) 1997-2006 Inger Jorgensen
# Copyright(c) 2000-2012 Association of Universities for Research in Astronomy, Inc.

procedure mgograph(intab,xcolumn,ycolumn)

# MONGO like graph directly from STSDAS tables or ASCII files
#    using st4gem.stplot.igi
#
# Version  Sept 1, 2000  IJ
#          Sept 20,2002  IJ v1.4 release
#
# You may modify this script, in which case the comments in the
# final script should contain the lines
# This script is based on original script written by
#    Inger Jorgensen, Gemini Observatory, Hilo, Hawaii, USA
#    e-mail: ijorgensen@gemini.edu

char  intab {prompt="Input table"}
char  xcolumn {prompt="X column"}
char  ycolumn {prompt="Y column"}
char  rows {"-",prompt="Rows to be plotted"}
real  wx1 {0.,prompt="X1, lower limit for X-axis"}
real  wx2 {0.,prompt="X2, upper limit for X-axis"}
real  wy1 {0.,prompt="Y1, lower limit for Y-axis"}
real  wy2 {0.,prompt="Y2, upper limit for Y-axis"}
char  excolumn {"",prompt="Error X column"}
char  eycolumn {"",prompt="Error Y column"}
bool  logx {no,prompt="Take log10 of X before plotting"}
bool  logy {no,prompt="Take log10 of Y before plotting"}
real  labelexp {2.,prompt="Expand factor for labels"}
real  boxexp {2.,prompt="Expand factor for box"}
char  xlabel {"",prompt="X label"}
char  ylabel {"",prompt="Y label"}
char  title{"",prompt="Title"}
char  postitle{"topleft",prompt="Position of title (topleft,topcenter,tl,tr,bl,br)",
      min="topleft|topcenter|tl|tr|bl|br"}
bool  append {no,prompt="Append"}
real  awx1 {0.,prompt="Saved X1, do not change"}
real  awx2 {0.,prompt="Saved X2, do not change"}
real  awy1 {0.,prompt="Saved Y1, do not change"}
real  awy2 {0.,prompt="Saved Y2, do not change"}
bool  pointmode {yes,prompt="Plot points instead of lines"}
char  marker {"box",prompt="Marker ([fill,skel]box|cross|plus|star|[fill,skel]triangle|[fill]circle|none)",
      min="box|fillbox|cross|plus|star|triangle|filltriangle|circle|fillcircle|skelbox|skeltriangle|none"}
real  pointexp {0.5,prompt="Expand factor for points"}
char  pcolumn {"",prompt="Point size column for variable point sizes"}
char  pattern {"solid",prompt="Line type (solid|dashed|dotted|dot-dash)",
      min="solid|dashed|dotted|dot-dash"}
char  crvstyle {"straight",prompt="curve style (straight|fullhist|pseudohist)",
      min="straight|fullhist|pseudohist"}
int   lweight {1,prompt="Line weight"}
int   color {1,prompt="Color, 1=black, 2=red, 3=green, 4=blue"}
bool  mkzeroline {no,prompt="Make zero line in y"}
char  device {"stdgraph",prompt="Device"}
char  gkifile {"mgo.gki",prompt="Gki-file if device=gkifile"}

begin

char n_intab, n_xcol, n_ycol, n_excol, n_eycol, n_xlabel, n_ylabel
bool n_point, n_append, n_mkz, n_logx, n_logy
char n_marker, n_pattern, n_device, n_gkifile, n_rows, n_crvstyle
char n_title, n_postitle, n_pcolumn
real n_wx1, n_wx2, n_wy1, n_wy2
real n_labelexp, n_pointexp, n_boxexp, n_pctype
int  n_lweight, n_color

char tmpigi, n_ptype, n_angle, n_ltype
char tmpd, tmpc, tmpdat
char s_fields, n_suf
bool n_ascii
int  n_test

n_intab = intab

# Check for [mdf]
n_test=0 ; n_test=stridx("[",n_intab)
if(n_test!=0) {
  n_suf=substr(n_intab,n_test,strlen(n_intab))
  n_intab=substr(n_intab,1,n_test-1)
} else
  n_suf=""

if((no==access (n_intab)) && (no==access (n_intab//".tab")) && \
  (no==access (n_intab//".fits"))) {
  print("ERROR - MGOGRAPH: Input table "//n_intab//" not found")
  bye
}

n_xcol = xcolumn
n_ycol = ycolumn
n_append = append
n_wx1 = wx1 ; n_wx2 = wx2 ; n_wy1 = wy1 ; n_wy2 = wy2
if(n_append) {
  n_wx1 = awx1 ; n_wx2 = awx2 ; n_wy1 = awy1 ; n_wy2 = awy2
}
n_rows = rows
n_excol = excolumn
n_eycol = eycolumn
n_logx = logx ; n_logy = logy
n_xlabel = xlabel
n_ylabel = ylabel
n_title = title
n_postitle = postitle
n_labelexp = labelexp
n_boxexp = boxexp
n_pointexp = pointexp
n_pcolumn = pcolumn
n_point = pointmode
n_mkz = mkzeroline
n_marker = marker
n_pattern = pattern
n_crvstyle = crvstyle
n_lweight = lweight
n_color = color
n_device = device
n_gkifile = gkifile

n_ptype = "4 0"
n_angle = "0"
n_pctype = 40
if(n_marker == "fillbox") {
  n_ptype = "4 3" ; n_pctype=43
}
if(n_marker == "plus") {
  n_ptype ="4 1" ; n_pctype=41
  n_angle ="45"
}
if(n_marker == "cross" || n_marker == "skelbox") {
  n_ptype = "4 1" ; n_pctype=41
}
if(n_marker == "triangle" ) {
  n_ptype = "3 0" ; n_pctype=30
}
if(n_marker == "skeltriangle" ) {
  n_ptype = "3 1" ; n_pctype=31
}
if(n_marker == "filltriangle") {
  n_ptype = "3 3" ; n_pctype=33
}
if(n_marker == "circle") {
  n_ptype = "25 0" ; n_pctype=250
}
if(n_marker == "fillcircle") {
  n_ptype = "25 3" ; n_pctype=253
}
if(n_marker == "star") {
  n_ptype = "6 1" ; n_pctype=61
}
if(n_marker == "none")
  n_point=no

n_ltype = "0"
if(no==n_point) {
  if(n_pattern == "dotted")
  n_ltype = "1"
  if(n_pattern == "dashed")
  n_ltype = "3"
  if(n_pattern == "dot-dash")
  n_ltype = "4"
}

tmpigi = mktemp("tmpigi")
tmpdat = mktemp("tmpdat")//".fits"
tmpd = mktemp("tmpd")
tmpc = mktemp("tmpc")

n_ascii = no
s_fields=n_xcol//","//n_ycol
if(n_excol!="")
  s_fields=s_fields//","//n_excol
if(n_eycol!="")
  s_fields=s_fields//","//n_eycol
if(n_pcolumn!="")
  s_fields=s_fields//","//n_pcolumn

# Add extension if FITS and no STSDAS table exists with same name
if(access (n_intab//".fits") && !access(n_intab//".tab") )
  n_intab=n_intab//".fits"

# check if table is empty
if(n_intab!="STDIN") {
  cache("tinfo")
  tinfo(n_intab//n_suf,ttout=no)
  if(tinfo.nrows<1) {
    print("ERROR - MGOGRAPH: Input table "//n_intab//" is empty")
    bye
  }
}

# if not-STSDAS or FITS make tmp-STSDAS
if(no==access (n_intab//".tab") && !access (n_intab//".fits") && \
  substr(n_intab,strlen(n_intab)-3,strlen(n_intab))!=".tab" \
  && substr(n_intab,strlen(n_intab)-3,strlen(n_intab))!="fits" ) {
    if(n_rows=="-")
      n_rows="1-999999"
    fields(n_intab,s_fields,lines=n_rows, >> tmpd)
    print("x r f8.3", > tmpc)
    print("y r f8.3", >> tmpc)
    if(n_excol!="" && n_excol!=" ") {
      print("e_x r f8.3", >> tmpc)
      n_excol="e_x"
    }
    if(n_eycol!="" && n_eycol!=" ") {
      print("e_y r f8.3", >> tmpc)
      n_eycol="e_y"
    }
    if(n_pcolumn!="" && n_pcolumn!=" ") {
      print("pcol r f8.3", >> tmpc)
      n_pcolumn="pcol"
    }
    tcreate(tmpdat,cdfile=tmpc,datafile=tmpd,upar="",tbltype="default")
    n_xcol="x"
    n_ycol="y"
    n_intab = tmpdat
    n_ascii = yes
} else {
  if(access (n_intab//".tab") && access (n_intab//".fits") ) {
    print("WARNING - MGOGRAPH: STSDAS table and FITS table with identical \
      names exist")
    print("                    plotting from STSDAS table")
    n_intab=n_intab//".tab"
  }
  tdump(n_intab//n_suf,cdfile=tmpc,datafile=tmpd,pfile="",
    columns=s_fields,rows=n_rows,pwidth=158, >> "dev$null")
  tcreate(tmpdat,cdfile=tmpc,datafile=tmpd,uparfile="",nskip=0,nlines=0,
    nrows=0,hist=no,extrapar=5,tbltype="default",extracol=0)
  n_intab = tmpdat
}
# add ptype to pcolumn to get final pcolumn
if(n_pcolumn!="") {
  tstat(n_intab, n_pcolumn, outtable="", >>& "dev$null")
  tcalc(n_intab, n_pcolumn, "("//n_pcolumn//"-"//str(tstat.vmin)//")*0.949/"\
    //str(tstat.vmax-tstat.vmin)//"+0.05+"//str(n_pctype))
}

if(n_logx) {
  tcalc(tmpdat,"lx","log10("//n_xcol//")",datatype="real",
    colunits="",colfmt="")
  n_xcol="lx"
  if(n_excol!="" && n_excol!=" ") {
    n_excol=""
    print("X-axis: Errorbars not supported for log-plots")
  }
}
if(n_logy) {
  tcalc(tmpdat,"ly","log10("//n_ycol//")",datatype="real",
    colunits="",colfmt="")
  n_ycol="ly"
  if(n_eycol!="" && n_eycol!=" ") {
    n_eycol=""
    print("Y-axis: Errorbars not supported for log-plots")
  }
}

# Auto scaling 

if(n_wx1==0 && n_wx2==0) {
  cache("tstat")
  tstat(n_intab,n_xcol, outtable="", >> "dev$null")
if(tstat.nrows>1 && tstat.stddev>0) {
  n_wx1 = real(tstat.vmin)-abs(real(tstat.stddev))
  n_wx2 = real(tstat.vmax)+abs(real(tstat.stddev))
} else {
  n_wx1 = real(tstat.vmin)-0.1*abs(real(tstat.vmin))
  n_wx2 = real(tstat.vmax)+0.1*abs(real(tstat.vmax))
}
}
if(n_wy1==0 && n_wy2==0) {
  tstat(n_intab,n_ycol, outtable="", >> "dev$null")
if(tstat.nrows>1 && tstat.stddev>0) {
  n_wy1 = real(tstat.vmin)-abs(real(tstat.stddev))
  n_wy2 = real(tstat.vmax)+abs(real(tstat.stddev))
} else {
  n_wy1 = real(tstat.vmin)-0.1*abs(real(tstat.vmin))
  n_wy2 = real(tstat.vmax)+0.1*abs(real(tstat.vmax))
}
}

# Make input file for igi
print("location 0.2 0.9 0.2 0.9", > tmpigi)
printf("%s %9.6e %9.6e %9.6e %9.6e \n","limits ",n_wx1,n_wx2,n_wy1,n_wy2, \
  >> tmpigi)
printf("color %s\n",n_color, >> tmpigi)
print("data "//n_intab, >> tmpigi)
print("xcolumn "//n_xcol, >> tmpigi)
print("ycolumn "//n_ycol, >> tmpigi)
print("lweight "//n_lweight, >> tmpigi)
if(no==n_append) {
  print("expand "//n_boxexp, >> tmpigi)
  print("box", >> tmpigi)
  print("expand "//n_labelexp, >> tmpigi)
  if(n_xlabel != "")
    print("xlabel "//n_xlabel, >> tmpigi)
  if(n_ylabel != "")
    print("ylabel "//n_ylabel, >> tmpigi)
  }
print("expand "//n_pointexp, >> tmpigi)
if(n_point) {
  if(n_pcolumn!="") 
    print("pcolumn "//n_pcolumn, >> tmpigi)
  else
    print("ptype "//n_ptype, >> tmpigi)
  print("angle "//n_angle, >> tmpigi)
  print("points ", >> tmpigi)
  print("angle 0", >> tmpigi)
}
else {
  print("ltype "//n_ltype, >> tmpigi)
  if(n_crvstyle=="straight")
    print("connect", >> tmpigi)
  if(n_crvstyle=="fullhist")
    print("histogram", >> tmpigi)
  if(n_crvstyle=="pseudohist")
    print("step", >> tmpigi)
}
if(n_excol != "") {
  print("ecolumn "//n_excol, >> tmpigi)
  print("errorbar -1", >> tmpigi)
  print("errorbar 1", >> tmpigi)
}
if(n_eycol != "") {
  print("ecolumn "//n_eycol, >> tmpigi)
  print("errorbar -2", >> tmpigi)
  print("errorbar 2", >> tmpigi)
}
if(n_mkz) {
  print("ltype 4", >> tmpigi)
  print("relocate "//n_wx1//" 0", >> tmpigi)
  print("draw "//n_wx2//" 0", >> tmpigi)
  print("ltype "//n_ltype, >> tmpigi)
}
if(n_title!="" && n_title!=" ") {
  print("expand "//n_labelexp, >> tmpigi)
  print("limits 0 1 0 1", >> tmpigi)
  if(n_postitle == "topleft") 
    print("relocate 0 1.07", >> tmpigi)
  if(n_postitle == "topcenter") 
    print("relocate 0.5 1.07", >> tmpigi)
  if(n_postitle == "tl")
    print("relocate 0.1 0.9", >> tmpigi)
  if(n_postitle == "bl")
    print("relocate 0.1 0.1", >> tmpigi)
  if(n_postitle=="topleft" || n_postitle=="tl" || n_postitle=="bl")
    print("putlabel 6 "//n_title, >> tmpigi)
  if(n_postitle == "tr")
    print("relocate 0.9 0.9", >> tmpigi)
  if(n_postitle == "br")
    print("relocate 0.9 0.1", >> tmpigi)
  if(n_postitle=="tr" || n_postitle=="br")
    print("putlabel 4 "//n_title, >> tmpigi)
  if(n_postitle=="topcenter")
    print("putlabel 5 "//n_title, >> tmpigi)
  
}
print("expand 1", >> tmpigi)
print("end", >> tmpigi)


if(n_device == "gkifile") 
  type(tmpigi) | igi(device="stdgraph",append=n_append, >>G n_gkifile)
else
  type(tmpigi) | igi(device=n_device,append=n_append)

delete(tmpigi,verify=no)
delete(tmpdat//","//tmpc//","//tmpd,verify=no, >>& "dev$null")

# Save cuts for appending

awx1=n_wx1
awx2=n_wx2
awy1=n_wy1
awy2=n_wy2

end
