# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gsdrawslits(image,mdf)

# Draw slits
# Version Dec 8, 2003 IJ Preparing task for release
# Version May 15, 2015 BM Many generalizations, added option to print slit ID labels

char image    {prompt="GMOS image"}
char mdf      {prompt="MDF for mask design"}
bool fl_disp  {yes,prompt="Display image"}
bool fl_slits {yes,prompt="Draw slits"}
bool fl_label {no,prompt="Label slits with ID"}
bool fl_spec  {no,prompt="Draw spectra"}
int  lwidth   {4,prompt="Line width on slits in pixels"}
int  frame    {1,prompt="Display frame"}
char sci_ext  {"1",prompt="Name of science extension"}
char drawprior {"all",prompt="Priorities to draw"}
struct *scanfile {"",prompt="Internal use only"}

begin

char n_image, n_mdf
char inst,dettype
bool n_fl_disp, n_fl_slits, n_fl_spec, n_fl_label
int  n_frame, n_lwidth

real x_ccd,y_ccd,slitsize_x,slitsize_y,ratio,slittilt, slitpos_y
real pixscale,pixscmdf,rscale
int  nrows, n_i, priority, x_max, n_j, id, xbin, ybin, mdfbin
bool n_prior[4]  # drawing flags for priority 0,1,2,3
char tmpin
char l_sci_ext

print("")
print("GSDRAWSLITS: GMOS slit drawing task -- prototype")
print("     this task has no logfile; tilted slits not drawn tilted")
print("")
cache("tinfo","keypar","fparse","gimverify")

n_image=image ; n_mdf=mdf ; n_fl_disp=fl_disp; n_frame=frame
n_fl_slits=fl_slits; n_fl_spec=fl_spec ; n_lwidth=lwidth
n_fl_label=fl_label
l_sci_ext=sci_ext

tmpin=mktemp("tmpin")

# Make sure table has .fits, don't loose the directory 
fparse(n_mdf) ; n_mdf=fparse.directory//fparse.root//".fits"

# Access checking
gimverify(n_image)
if(gimverify.status!=0) {
  print("GSDRAWSLITS -- ERROR: Input image "//n_image//" does not exist")
  bye
}
if (!access(n_mdf)) {
  print("GSDRAWSLITS -- ERROR: Input MDF "//n_mdf//" does not exist")
  bye
}

# Set drawing flags
if(drawprior=="all") {
   for(n_i=1;n_i<=4;n_i+=1)
     n_prior[n_i]=yes
} else {
   for(n_i=1;n_i<=4;n_i+=1)
     n_prior[n_i]=no 
  files(drawprior, > tmpin)
  scanfile=tmpin
  while(fscan(scanfile,n_i)!=EOF) {
     n_prior[n_i+1]=yes
  }
  scanfile="" ; delete(tmpin,verify-)
}

if(n_fl_disp) 
  display(n_image//"["//l_sci_ext//"]",n_frame)

keypar(n_image//"[0]","INSTRUME")
inst = keypar.value

keypar(n_image//"[0]","DETTYPE")
dettype = keypar.value

ybin = 1
xbin = 1
keypar(n_image//"["//l_sci_ext//"]","CCDSUM")
print(keypar.value) | scan(xbin,ybin) 

keypar(n_image//"[0]","PIXSCALE")
if(keypar.found) {
    pixscale=real(keypar.value)
} else {
    if (dettype == "SDSU II CCD") { # EEV CCDs
        if (inst == "GMOS-S") {
            pixscale=0.073*xbin
        } else {
            pixscale=0.0727*xbin
        }
    } else if (dettype == "SDSU II e2v DD CCD42-90") { # e2vDD CCDs
        pixscale=0.07288*xbin
    } else if (dettype == "S10892") { # Hamamatsu CCDs  SOUTH
        pixscale=0.0800*xbin  ##M PIXEL_SCALE
    } else if (dettype == "S10892-N")  {  #Hamamatsu NORTH
        pixscale=0.0807*xbin  
    }
}

keypar(n_image//"["//l_sci_ext//"]","i_naxis1")
x_max=int(keypar.value)

# pixel scale, MDF
keypar(n_mdf,"PIXSCALE")
if(keypar.found) {
  pixscmdf=real(keypar.value)
} else {
  # Pre-image most likely is 2x2 binned
  pixscmdf=0.0727*2.0
}
# Try to homogenize the pixel scales
# Have to use pixel scales assumed by GMMPS
if (dettype == "S10892") {
    mdfbin = nint(pixscmdf/0.080)
} else if (dettype == "S10892-N") {
    mdfbin = nint(pixscmdf/0.0807)  # KL get pixscale from GMMPS.
} else {
    mdfbin = nint(pixscmdf/0.0727)
}

rscale=real(xbin)/real(mdfbin)

tinfo(n_mdf,ttout-) ; nrows=tinfo.nrows
print("red    prior 1")
print("green  prior 2")
print("blue   prior 3")
print("yellow prior 0")
for(n_i=1;n_i<=nrows;n_i+=1) {
slittilt=0
slitpos_y=0
  tprint(n_mdf,col="id,x_ccd,y_ccd,slitsize_x,slitsize_y,priority",
   showr-,showh-,row=str(n_i),option="plain",align-) | \
   scan(id,x_ccd,y_ccd,slitsize_x,slitsize_y,priority)
  tprint(n_mdf,col="slitpos_y",
   showr-,showh-,row=str(n_i),option="plain",align-) | \
   scan(slitpos_y)
 ratio=slitsize_y/slitsize_x
 slitsize_x=int(slitsize_x/pixscale+0.5)
 slitsize_y=int(slitsize_y/pixscale+0.5)
 x_ccd = x_ccd/rscale
 y_ccd = y_ccd/rscale + slitpos_y/pixscale # Correct offset slits
# 203 = white
# 204 = red     prior 1
# 205 = green   prior 2
# 206 = blue    prior 3
# 207 = yellow  acq objects
if(n_fl_slits) {
 if(priority==0 && n_prior[1]) {
  for(n_j=1;n_j<=n_lwidth;n_j+=1) {
  ratio = real(slitsize_y-(n_j-1))/real(slitsize_x-(n_j-1))
   printf("%8.2f %8.2f %d\n",x_ccd,y_ccd,id) | tvmark(n_frame,"STDIN",
    mark="rectangle",label=n_fl_label,nxoffset=(10/xbin),txsize=3,
    lengths=str(slitsize_x-(n_j-1))//" "//str(ratio),
    color=207)
  }
 } else if (priority==1 && n_prior[2]) {
  for(n_j=1;n_j<=n_lwidth;n_j+=1) {
  ratio = real(slitsize_y-(n_j-1))/real(slitsize_x-(n_j-1))
   printf("%8.2f %8.2f %d\n",x_ccd,y_ccd,id) | tvmark(n_frame,"STDIN",
    mark="rectangle",label=n_fl_label,nxoffset=(10/xbin),txsize=3,
    lengths=str(slitsize_x-(n_j-1))//" "//str(ratio),
    color=204)
  }
 } else if (priority==2 && n_prior[3]) {
  for(n_j=1;n_j<=n_lwidth;n_j+=1) {
  ratio = real(slitsize_y-(n_j-1))/real(slitsize_x-(n_j-1))
   printf("%8.2f %8.2f %d\n",x_ccd,y_ccd,id) | tvmark(n_frame,"STDIN",
    mark="rectangle",label=n_fl_label,nxoffset=(10/xbin),txsize=3,
    lengths=str(slitsize_x-(n_j-1))//" "//str(ratio),
    color=205)
  }
 } else if (n_prior[4]) {
  for(n_j=1;n_j<=n_lwidth;n_j+=1) {
  ratio = real(slitsize_y-(n_j-1))/real(slitsize_x-(n_j-1))
   printf("%8.2f %8.2f %d\n",x_ccd,y_ccd,id) | tvmark(n_frame,"STDIN",
    mark="rectangle",label=n_fl_label,nxoffset=(10/xbin),txsize=3,
    lengths=str(slitsize_x-(n_j-1))//" "//str(ratio),
    color=206)
  }
 }
} # end of slits
if(n_fl_spec) {
  for(n_j=1;n_j<=n_lwidth;n_j+=1) {
   printf("%8.2f %8.2f 101 s\n",1,int(y_ccd-slitsize_y/2.+0.5+(n_j-1)), >> tmpin) 
   printf("%8.2f %8.2f 101 s\n",x_max,int(y_ccd-slitsize_y/2.+0.5+(n_j-1)), >> tmpin) 
   printf("%8.2f %8.2f 101 o\n",1,int(y_ccd-slitsize_y/2.+0.5-(n_j-1)), >> tmpin) 
   printf("%8.2f %8.2f 101 s\n",1,int(y_ccd+slitsize_y/2.+0.5-(n_j-1)), >> tmpin) 
   printf("%8.2f %8.2f 101 s\n",x_max,int(y_ccd+slitsize_y/2.+0.5-(n_j-1)), >> tmpin) 
   printf("%8.2f %8.2f 101 o\n",1,int(y_ccd-slitsize_y/2.+0.5+(n_j-1)), >> tmpin) 
  }
   printf("q\n", >> tmpin)
 if(priority==0 && n_prior[1]) {
    print("0 0") | tvmark(n_frame,"STDIN",commands=tmpin,color=207,label-)
 } else if (priority==1 && n_prior[2]) {
    print("0 0") | tvmark(n_frame,"STDIN",commands=tmpin,color=204,label-)
 } else if (priority==2 && n_prior[3]) {
    print("0 0") | tvmark(n_frame,"STDIN",commands=tmpin,color=205,label-)
 } else if (n_prior[4]) {
    print("0 0") | tvmark(n_frame,"STDIN",commands=tmpin,color=206,label-)
 }
delete(tmpin,verify-)
} # end of spectra

}

end
