# Copyright(c) 2015-2017 Association of Universities for Research in Astronomy, Inc.

procedure gfscatsub(image,mask)

# Subtracts background/dark from IFU frames
# Bryan Miller
# 2005apr15 - created

string image	{prompt="Image for background subtraction"}
string mask     {prompt="Mask file with inter-gap regions (fndblocks)"}
string outimage {"",prompt="Output image with background subtracted"}
string prefix   {"b",prompt="Output image prefix"}
#int    xorder   {1,prompt="X order for surface fit"}
#int    yorder   {3,prompt="Y order for surface fit"}
string  xorder  {"1",prompt="X order for surface fit"}
string  yorder  {"3",prompt="Y order for surface fit"}
bool   cross    {no,prompt="Fix cross terms"}
bool    fl_inter {no, prompt="display surface fit for visual inspection?"}
int    status   {0,prompt="Exit status (0=good)"}
struct* scanfile {"",prompt="Internal use only"}

begin

string l_image,l_mask,l_outimage,l_prefix
string bkgimg,tmporder
int nsci,ii,jj,nc,dum
int l_xorder[12], l_yorder[12]
bool l_cross
bool l_fl_inter

# Define local input variables
l_image=image ; l_mask=mask ; l_outimage=outimage ; l_prefix=prefix
l_cross=cross
l_fl_inter=fl_inter

status = 0

cache("imgets","gimverify")

tmporder=mktemp("tmporder")

gimverify(l_image)
if (gimverify.status != 0) {
    print("ERROR - input image not a valid MEF")
    goto error
}
l_image=gimverify.outname//".fits"

imgets(l_image//"[0]","NSCIEXT")
nsci=int(imgets.value)

# Check whether already run
imgets(l_image//"[0]","GFSCATSUB", >& "dev$null")
if (imgets.value != "0") {
    print("ERROR - GFSCATSUB: "//l_image//" already processed with GFSCATSUB. Exiting.")
    goto clean
}

if (l_outimage == "") {
   l_outimage = l_prefix//l_image
}
gimverify(l_outimage)
if (gimverify.status != 1) {
    print("ERROR - GFSCATSUB: "//l_outimage//" already exists.")
     goto error
}
l_outimage = gimverify.outname//".fits"

files(xorder,sort-, > tmporder)
count(tmporder) | scan(nc)
if (nc == 1) {
    for (ii=1;ii<=nsci;ii+=1) {
	l_xorder[ii]=int(xorder)
    }
} else if (nc == nsci) {
    ii=0
    scanfile=tmporder
    while(fscan(scanfile,dum) != EOF) {
	ii=ii+1
	l_xorder[ii]=dum
    }
    scanfile=""
} else {
    print("ERROR - number of xorder entries must be one or NSCI")
    goto error
}
delete(tmporder,verify-)

files(yorder,sort-, > tmporder)
count(tmporder) | scan(nc)
if (nc == 1) {
    for (ii=1;ii<=nsci;ii+=1) {
	l_yorder[ii]=int(yorder)
    }
} else if (nc == nsci) {
    ii=0
    scanfile=tmporder
    while(fscan(scanfile,dum) != EOF) {
	ii=ii+1
	l_yorder[ii]=dum
    }
    scanfile=""
} else {
    print("ERROR - number of yorder entries must be one or NSCI")
    goto error
}
delete(tmporder,verify-)

copy(l_image,l_outimage,ver-)

for (ii=1;ii<=nsci;ii+=1) {
    bkgimg=mktemp("tmpbkgimg")
    imsurfit(l_outimage//"[SCI,"//ii//"]",bkgimg,l_xorder[ii],l_yorder[ii],
	cross=l_cross,regions="sections",type_output="fit",function="leg",
	sections=l_mask,upper=3.,lower=3.,niter=3,xmedian=1,ymedian=1,
        ngrow=0,rows="*",columns="*",border=50)
    if (l_fl_inter) {
        imexamine(bkgimg,1)
    }
#    display(bkgimg,1)
    imarith(l_outimage//"[SCI,"//ii//"]","-",bkgimg,
	l_outimage//"[SCI,"//ii//",overwrite+]")
    # If a VAR plane should add the square of bkgimg
    imdelete(bkgimg,ver-)
}
gemdate()
gemhedit (l_outimage//"[0]", "GFSCATSUB", gemdate.outdate,
        "UT Time stamp for GFSCATSUB")

goto clean

error:
  # here on error
  status = 1

clean:
   delete(tmporder,verify-,>>& "dev$null")

end
