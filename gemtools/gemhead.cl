# Copyright(c) 2000-2011 Association of Universities for Research in Astronomy, Inc.

procedure gemhead(images)

# List headers of MEF images
#
# Version Feb 28, 2002  IJ v1.3 release
#         Sept 20,2002  IJ v1.4 release

char images {prompt="Images"}
bool longheader  {yes,prompt="List long headers"}
struct *scanfile {prompt="Internal use only"}

begin

char l_images, l_input
bool l_long
int  n_i, n_ext
char tmpin

l_images=images
l_long=longheader
tmpin=mktemp("tmpin")

cache("gimverify")
files(l_images, > "uparm$"//tmpin)
scanfile="uparm$"//tmpin

while(fscan(scanfile,l_input)!=EOF) {

n_ext=1
gimverify(l_input)
if(gimverify.status!=0 && gimverify.status!=4)
   print("GEMHEAD - WARNING: Image "//l_input//" not found or is not FITS")
else {
fxhead(l_input, > "uparm$"//tmpin//"2")  
system.tail ("uparm$"//tmpin//"2", nlines=1) | scan(n_ext)
delete("uparm$"//tmpin//"2",verify-, >>& "dev$null")

  for(n_i=0;n_i<=n_ext;n_i+=1) 
    imhead(l_input//"["//str(n_i)//"]",long=l_long)

}

}
scanfile=""
delete("uparm$"//tmpin,verify-, >>& "dev$null")

end

