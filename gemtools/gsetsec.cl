# Copyright(c) 2001-2006 Association of Universities for Research in Astronomy, Inc.

procedure gsetsec(image)

# Reset datasec on MEF images to current image size
#
# Version Oct  12, 2001  IJ  v1.2 release
#         Sept 20, 2002  IJ  v1.4 release

char image      {"",prompt="MEF image to update"}
char key_datsec  {"DATASEC",prompt="Header keyword for data section"}
struct* scanfile {"",prompt="Internal use only."}

begin

char l_image,  l_key_datsec

char tmpin, l_datsec
int n_i, Xmin, Xmax, Ymin, Ymax
struct l_struct

l_image=image ; l_key_datsec=key_datsec 
tmpin=mktemp("tmpin")

cache("imgets","gimverify")
date | scan(l_struct)

gimverify(l_image)
if(gimverify.status==1) {
  print("ERROR - GSETSEC: Cannot access image "//l_image)
  bye
}
if(gimverify.status!=0) {
  print("ERROR - GSETSEC: Image "//l_image//" is not a MEF file")
  bye
}

fxhead(l_image) | match("IMAGE","STDIN",stop-, > tmpin)

scanfile=tmpin
while(fscan(scanfile,n_i)!=EOF) {

# Get the size of the axis
 imgets(l_image//"["//str(n_i)//"]","i_naxis1")
 Xmax=int(imgets.value)
 imgets(l_image//"["//str(n_i)//"]","i_naxis2")
 Ymax=int(imgets.value)

 gemhedit(l_image//"["//str(n_i)//"]",l_key_datsec,
   "[1:"//str(Xmax)//",1:"//str(Ymax)//"]","Data section")

}
scanfile=""

delete(tmpin,verify-, >>& "dev$null")

end
