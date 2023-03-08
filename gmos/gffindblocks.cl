# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.

procedure gffindblocks(image,extspec,mask)

# Finds gaps between IFU blocks
# Bryan Miller
# 2005apr15 - created
# 2005dec07 - make strips narrower
# 2013mar15 - JT: rewrote to use measurements from both slits
# 2014sep15 - JT: adjust expected blocks / increment to allow for N+S flats

string image	   {"",prompt="Unextracted image"}
string extspec     {"",prompt="Extracted spectra"}
string mask	   {"",prompt="Out mask 'regions' file"}
int    status      {0,prompt="Exit status (0=good)"}
struct *scanfile   {"",prompt="Internal use only"}

begin

string l_image,l_extspec,l_mask
string apnum
int  i, i1, j, k, nap1, nap2, ninc, ngap, ngaps, nx, ny, num1[750], num2[750]
int  bef1, bef2, aft1, aft2, s2off, lastap1, lastap2, glo, ghi
real x, y, ylo1[750], yhi1[750], ylo2[750], yhi2[750], glor, ghir, gcen, gsiz


# define variables
l_image=image
l_extspec=extspec
l_mask=mask

cache("gimverify")

status=0

gimverify(l_image)
if (gimverify.status != 0) {
  print("ERROR - "//l_image//" doesn't exist or isn't a MEF file")
  status=1
  bye
}
gimverify(l_extspec)
if (gimverify.status != 0) {
  print("ERROR - "//l_extspec//" doesn't exist or isn't a MEF file")
  status=1
  bye
}

if (access(l_mask)) {
  print("Mask file "//l_mask//" already exists")
  status=1
  bye
}

# Header info
imgets(l_image//"[sci,1]","i_naxis1")
nx=int(imgets.value)
imgets(l_image//"[sci,1]","i_naxis2")
ny=int(imgets.value)

# Read the aperture (fibre) limits for the first slit:
apnum=mktemp("tmpapnum")
imhead(l_extspec//"[sci,1]",l+) | match("APNUM","STDIN") \
 | translit("STDIN", "APNUM'=", "", del+, > apnum)

scanfile=apnum
while(fscan(scanfile,i,j,k,x,y) != EOF) {
    num1[i]=j
    ylo1[i]=x
    yhi1[i]=y
}
scanfile=""
nap1=i

# If a second slit was extracted, read its aperture limits, otherwise
# we'll use the first set again:
s2off=0
if (imaccess(l_extspec//"[sci,2]")) {
    delete(apnum,verify-)
    imhead(l_extspec//"[sci,2]",l+) | match("APNUM","STDIN") \
     | translit("STDIN", "APNUM'=", "", del+, > apnum)
    s2off=750
}

scanfile=apnum
while(fscan(scanfile,i,j,k,x,y) != EOF) {
    num2[i]=j
    ylo2[i]=x
    yhi2[i]=y
}
scanfile=""
nap2=i

delete(apnum,verify-)

# Distinguish the number of blocks in N+S modes from the usual 15 per slit
# based on the number of apertures found in the header:
if (nap1 > 350) {
  ngaps = 16
  ninc = 50
} else {
  ngaps = 8
  ninc = 100
}

# Use a fixed region 6 pixels below the first fibre as the first "gap":
ghi = int(max(1, min(ylo1[1], ylo2[1]) - 6))
glo = max(1, ghi - 4)

print("1 "//nx//" "//glo//" "//ghi)
print("1 "//nx//" "//glo//" "//ghi, > l_mask)

# Loop over a known number of gaps -- this isn't really much less general
# than Bryan's method, as we were already assuming 50-fibre blocks anyway,
# and this makes it easier to work with 2 slits at once.
lastap1=0; lastap2=s2off
for (ngap=2; ngap < ngaps; ngap+=1) {

    # Identify the nominal fibre right before this gap:
    i1 = (ngap-1)*50
    lastap1 += ninc
    lastap2 += ninc
    bef1=i1
    bef2=i1

    # Adjust for missing fibres:
    # bef1/bef2 are the fibre array indices (not numbers) before the gap
    for (i=i1; num1[i] > lastap1; i=i-1);
    bef1=i
    for (i=i1; num2[i] > lastap2; i=i-1);
    bef2=i

    aft1=bef1+1
    aft2=bef2+1

    # print("*", num1[bef1], num1[aft1], num2[bef2], num2[aft2])

    # Determine the gap in common between the 2 slits:
    glor = max(yhi1[bef1], yhi2[bef2])
    ghir = min(ylo1[aft1], ylo2[aft2])
    # print ("limits1 ", glor, ghir)
    gcen = 0.5 * (glor+ghir)
    gsiz = max(1, 0.5 * (ghir-glor) - 5)   # 5 pix buffer from bounding fibs
    # Here we deliberately round down because the GMOS-S distortions otherwise
    # put the regions right at the high end of the gap for slit 2:
    glo = int(gcen - gsiz)
    ghi = int(gcen + gsiz)

    # Record the region to use in the output:
    print("1 "//nx//" "//glo//" "//ghi)
    print("1 "//nx//" "//glo//" "//ghi, >> l_mask)

}

# For the last "gap" at the top, use a fixed region 6 pixels above the
# last fibre:
glo = int(min(ny, max(yhi1[nap1], yhi2[nap2]) + 6))
ghi = min(ny, glo + 4)

# Only use this last region if there's really a gap at the top (for GMOS-S
# there often isn't one):
if (glo < ny) {
    print("1 "//nx//" "//glo//" "//ghi)
    print("1 "//nx//" "//glo//" "//ghi, >> l_mask)
    print(ngap//" gaps found")
} else {
    print((ngap-1)//" gaps found")
}

end
