# Copyright(c) 2001-2006 Association of Universities for Research in Astronomy, Inc.

procedure ohead(inimages)

# Print frame information for OSCIR files
# Information is printed to screen only, not logged.
# NOTES: 16nov01- need to print more header info for checking by DP
#
# Version: Sept 14, 2002 BR  Release v1.4

char inimages	 {prompt="Input OSCIR image(s)"}
bool fl_showqa	 {yes,prompt="Show QA header values"}
#bool fl_imstat	 {no,prompt="Show imstat info"}
struct* scanfile

begin

char l_inimages 
bool l_fl_showqa,l_fl_imstat
char l_temp, tmpstat, tmpin, tmpfile, tmpfile2
char in[100], n_image
int i, nimages 
int naxis3, naxis4, naxis5, naxis6

l_inimages=inimages
l_fl_showqa=fl_showqa 
#l_fl_imstat=fl_imstat

tmpin = mktemp("tmpin")
tmpfile = mktemp("tmpfile")
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

while(fscan(scanfile,l_temp) != EOF) {
  if(substr(l_temp,strlen(l_temp)-4,strlen(l_temp))==".fits" )
    l_temp=substr(l_temp,1,strlen(l_temp)-5)

  if(!imaccess(l_temp))
     print("WARNING - OHEAD: Input image "//l_temp//" not found.")
  else {
     nimages=nimages+1
  
     in[nimages]=l_temp 

# Catch $  and / in input
    if((stridx("$",in[nimages])!=0)||(stridx("/",in[nimages])!=0)) {
      printf("ERROR - OHEAD: Error in input file name")
      goto clean
    }
 }
}
printf("Processing "//nimages//" file(s).\n")
scanfile="" ; 
if(nimages==0) {
  print("ERROR - OHEAD: No existing input images defined")
  goto clean
}

cache("imgets")

#--------------------------------------------------------------------------
#MAIN LOOP

i=1
while(i<=nimages) {
  n_image=in[i]

  imgets(n_image,"i_naxis3") ; naxis3=int(imgets.value)
  imgets(n_image,"i_naxis4") ; naxis4=int(imgets.value)
  imgets(n_image,"i_naxis5") ; naxis5=int(imgets.value)
  imgets(n_image,"i_naxis6") ; naxis6=int(imgets.value)

  if (i==1) printf("%-25s %-3s %-3s %-3s %-3s %-3s\n","Image","n_choppos","n_savesets","n_nodpos","n_nodsets","tot_frames")
  printf("%-25s   %-9d %-10d %-8d %-9d %-9d\n",n_image,naxis3,naxis4,naxis5,naxis6,naxis3*naxis4*naxis5*naxis6)

  #if (l_fl_imstat) {
  #   tmpstat =mktemp("tmpstat")
  #   printf("\nFirst and last frame statistics:\n")
#
#     print(n_image//"[*,*,1,1,1,1]",>> tmpstat)
#     print(n_image//"[*,*,"//naxis3//","//naxis4//","//naxis5//","//naxis6//"]",>> tmpstat)
#
#     imstat ("@"//tmpstat)
#     delete (tmpstat,verify-)
#  }

  i=i+1

}

  if (l_fl_showqa) {
	printf("%-25s %-14s %-15s %-15s %-15s %-15s\n",
	"Filename","OBSID","RAWBG","RAWCC","RAWIQ","RAWWV",yes)
	hselect(l_inimages,"$I,OBSID,RAWBG,RAWCC,RAWIQ,RAWWV",yes)
  }

clean:
{
 delete(tmpfile//","//tmpin,ver-, >>& "dev$null")
}

end
