# Copyright(c) 2002-2006 Association of Universities for Research in Astronomy, Inc.

procedure gemvsample ()

# Check the validity of a string specifying a sample of pixels in 1D
# (eg. y1:y2,y3:y4 or "*") and ensure that the pixels are within bounds.
#
# Inside apextract, this is done in the subroutine ap_skyeval (apskyeval.x).
# However, we will adopt a more stringent set of criteria, since apall
# accepts nonsensical strings with little indication of an error.
#
# Numbers & ranges can be out of ascending order or repeated here (except
# that an asterisk is only allowed on its own), since apall seems to cope
# with that. Only fixed point notation is allowed.
#
# The parameter zeropt allows co-ordinates to be specified either as absolute
# pixel indices, or relative to the image centre (in y). It is not possible
# to specify co-ordinates relative to an aperture centre, since the MDF does
# not contain an actual object offset for long-slit mode. Hence when checking
# a sample string for apall, the boundary check is approximate and assumes
# that the object is centred along the slit length. For this reason, error
# code 1 should only generate a warning in the calling task.
#
# Status: 0 Good
#         1 Sample out of bounds (possibly, if relative to an aperture)
#         2 Invalid sample string
#         3 Can't open image extension, MDF or header keyword
#
# Error status 2 takes precedence over 3, so that the string can be checked
# without requiring an input image extension.
#
# Version  Jun 06, 2002  JT  v1.4 release

string sample {"", prompt="Sample string to check"}
string image  {"", prompt="Corresponding image (extension)"}
string zeropt {"center", enum="center|firstpix",
               prompt="Co-ordinate zero point"}
int status {0, prompt="Exit status"}

begin

# Local variables:
string l_sample, l_image, l_zeropt, thischar, thisnum
int n, sslen, nsec, ybin, ydim, tint
real minval, maxval, mingood, maxgood
bool firstnum, newnum, newrange, gotdigit, gotpoint, gotcolon, gotast,
     gotcomma

# Avoid parameters changing externally:
cache("imgets","gimverify")

# Read & save parameters:
l_sample = sample
l_image = image
l_zeropt = zeropt
sample = l_sample
image = l_image
zeropt = l_zeropt

# Default status:
status = 2

# Get sample string length:
sslen=strlen(l_sample)

# Initialize logic:
newnum = yes
newrange = yes
gotdigit = no
gotpoint = no
gotcolon = no
gotast = no
gotcomma = yes
firstnum = yes
thisnum = ""
minval = 0
maxval = 0

# Loop through the characters of the sample string; verify validity and
# store min/max numbers:
for (n=1; n <= sslen; n+=1) {

  thischar = substr(l_sample,n,n)

  # Copy next digit or compare last num with min/max; exit if error:
  switch(thischar) {

    case '0','1','2','3','4','5','6','7','8','9': {
      if (gotast) break
      thisnum = thisnum//thischar
      gotdigit = yes
      gotcomma = no
      newnum = no
      newrange = no
    }

    case '-': {
      if (no == newnum) break
      thisnum = thisnum//thischar
      gotcomma = no
      newnum = no
      newrange = no
    }

    case '.': {
      if (gotpoint || gotast) break
      thisnum = thisnum//thischar
      gotcomma = no
      gotpoint = yes
      newnum = no
      newrange = no
    }

    case ':': {
      if (no == gotdigit || gotcolon || gotast) break
      if (firstnum || real(thisnum) < minval) {
        minval = real(thisnum)
      }
      if (firstnum || real(thisnum) > maxval) {
        maxval = real(thisnum)
        firstnum = no
      }
      newnum = yes
      gotdigit = no
      gotpoint = no
      gotcolon = yes
      thisnum=""
    }

    case ',',' ','\t': {
      if (no == (gotdigit || gotast) && no == newrange) break
      if (thischar == ',') {
        if (gotcomma) break
        else gotcomma = yes
      }
      if (gotdigit) {
        if (firstnum || real(thisnum) < minval)
          minval = real(thisnum)
        if (firstnum || real(thisnum) > maxval) {
          maxval = real(thisnum)
          firstnum = no
        }
      }
      newnum = yes
      newrange = yes
      gotdigit = no
      gotpoint = no
      gotcolon = no
      thisnum=""
    }

    case '*': {
      if (no == firstnum || no == newnum) break
      gotast = yes
      gotcomma = no
      newnum = no
      newrange = no
      firstnum = no
    }

    default:
      break

  } # end switch (thischar)

} # end for (n <= sslen)


# Continue checks / set status if loop finished without err:
if (n > sslen) {

  # Error if unfinished number at end:
  if (((no == newnum || gotcolon) && no == (gotdigit || gotast)) || gotcomma)
    goto error

  # Compare final number against min/max:
  if (gotdigit) {
    # Determine final min/max and compare with image:
    if (firstnum || real(thisnum) < minval)
      minval = real(thisnum)
    if (firstnum || real(thisnum) > maxval)
      maxval = real(thisnum)
  }

  # New default status:
  status = 3

  # Check input image/extension is readable:
  if (no == imaccess(l_image)) goto error

  # Check whether file is MEF:
  gimverify(l_image)

  # If MEF, check that an extension was specified:
  # (must be an extension if '[' present and imaccess==yes)
  if (gimverify.status == 0 && stridx("[", l_image) == 0) goto error

  # Check that the extension was not 0:
  imgets(l_image, "EXTEND", >& "dev$null")
  if(imgets.value =="T" ) goto error

  # Get dimensions of input extension:
  imgets(l_image, "naxis2", >& "dev$null")
  ydim = int(imgets.value)

  # Set co-ordinate bounds according to zero-point:
  if (l_zeropt == "center") {
	mingood = 0.5-0.5*real(ydim)
	maxgood = -mingood
  }
  else if (l_zeropt == "firstpix") {
    mingood = 1.0
	maxgood = real(ydim)
  }
		
  # Set status depending on range (centre +min/max):
  if (no == gotast && (minval < mingood || maxval > maxgood))
    status = 1
  else
    status = 0

}

# Jump here in event of error:
error: ;

end
