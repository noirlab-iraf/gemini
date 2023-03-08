# Copyright(c) 2000-2012 Association of Universities for Research in Astronomy, Inc.

procedure gimverify(image)

# status:
# 0 - exists and is a MEF file
# 1 - does not exist
# 2 - exists and is an imh file
# 3 - exists and is an hhh file
# 4 - exists and is a simple fits file
#
# Version: Oct 12, 2001 MT,IJ  v1.2 release
#          Feb 28, 2002 MT,IJ  v1.3 release
#          Sept 20,2002 MT,IJ  v1.4 release

# This task will check for common input errors and adjust the outname
# accordingly. If using this script to verify the presence of an image one
# should use the outname value beyond the validity check, as the input name and
# the output name may differ beyond the normal difference in the file
# extension. - MS
#
# NOTE: gemextn may be a better choice to validate a file.

string  image        {prompt="Input image"}
string  outname      {"", prompt="Output name w/o suffix"}
int     status       {0, prompt="Status"}

begin

string l_image, rear, envvar, tmpstring, tmpin
int    dlocation, slocation, sspos, l_test
bool   mef, debug, dslash

# Initiliaze param
status = 0
debug = no
dslash = no

# Initaiate variables
dlocation = 0
slocation = 0
tmpstring = ""

# Query parameters
l_image = image

if (debug) {
    print ("GIMVERIFY: l_image (1) is: "//l_image)
}

if(stridx("[",l_image)>0)
     l_image = substr(l_image,1,stridx("[",l_image)-1)

####

# Check for // in filename
sspos = strstr("//", l_image)
if (debug) {
    print ("GIMVERIFY: sspos is: "//sspos)
}

l_test = sspos
while (l_test > 0) {

    if (sspos == 1) {
        l_image = substr(l_image,2,strlen(l_image))
    } else if (sspos > 1) {
        l_image = substr(l_image,1,sspos)//\
            substr(l_image,sspos+2,strlen(l_image))
    }
    l_test = strstr("//", l_image)
    if (debug) print ("GIMVERIFY: sspos (2) is: "//l_image)
}

if (debug) {
    print ("GIMVERIFY: l_image (4) is: "//l_image)
}

# Check for environmental variables in the image name
#

# Locate any (last) "$" (closest to filename where extra "/" is likely to be
# placed)
dlocation = strldx("$",l_image)

# Locate last (if any) "/"
slocation = strldx("/",l_image)

if (debug) {
    print ("GIMVERIFY: l_image (2) is: "//l_image)
    print ("GIMVERIFY: dlocation is: "//dlocation)
    print ("GIMVERIFY: slocation is: "//slocation)
}

# This doesn't work if there is more than one environmental variable set in the
# name - MS
if (dlocation > 0) {

    # Assume the environmental varibale is either preceeded by / or is at the
    # begining of the string
    if (slocation > dlocation) {
        # Check for "$/"
        if (slocation == dlocation + 1) {
            dslash = yes
        }
        slocation = 0
    }

    tmpstring = substr(l_image, slocation+1, dlocation-1)
    tmpin = tmpstring//"$"

    if (debug) {
        print ("GIMVERIFY: tmpstring (1) is: "//tmpstring)
        print ("GIMVERIFY: tmpin (1) is: "//tmpin)
    }

    if (tmpstring != "") {

        # Read the environmental variable
        envvar = ""

        if (!defvar(tmpstring)) {

            # Envrionmental varibale is either a null string or doesn't exist
            # Skip the following tests and allow the imaccess calls etc. later
            # take care of that fact, i.e., to check if the file exits - MS
            goto skip_expansion

        } else {
            show (tmpstring) | scan (envvar)
        }

        if (envvar == "") {
            goto skip_expansion
        }

        # Form the expanded version
        if (slocation > 0) {
            tmpstring = substr(l_image, 1, slocation)//envvar
        } else {
            tmpstring = envvar
        }

        if (debug) {
            print ("GIMVERIFY: tmpstring (2) is: "//tmpstring)
        }

        # Check access - will work if directory and no "/" at the end
        # If no access assume the envronmental variable is just a prefix
        if (access(tmpstring)) {
            if (debug) {
                print ("GIMVERIFY: accessed "//tmpstring//" dslash is: "//\
                    dslash)
            }
            # Check for a slash
            if (substr(tmpstring,strlen(tmpstring),\
                strlen(tmpstring)) != "/") {

                if (!dslash) {
                    tmpin = tmpin//"/"
                }
            } else if (dslash) {
                    dlocation += 1
            }
        }
        if (debug) {
            print ("GIMVERIFY: tmpstring (4) is: "//tmpstring)
            print ("GIMVERIFY: tmpin (4) is: "//tmpin)
        }

        l_image = tmpin//substr(l_image, dlocation+1, strlen(l_image))
    }
}
####

skip_expansion:

if (debug) {
    print ("GIMVERIFY: l_image (3) is: "//l_image)
}

# Check existance status of l_image
rear = substr(l_image,strlen(l_image)-3,strlen(l_image))
if(no == imaccess(l_image)) {
     status=1
     if (rear=="fits") l_image = substr(l_image,1,strlen(l_image)-5)
     if (rear==".imh" || rear==".hhh") l_image = substr(l_image,1,strlen(l_image)-4)
} else if (rear==".imh" || imaccess(l_image//".imh")) {
     status=2
} else if (rear==".hhh" || imaccess(l_image//".hhh")) {
     status=3
} else {

     if (rear!="fits") rear = "fits"
     else if (rear=="fits") l_image = substr(l_image,1,strlen(l_image)-5)

     mef=no
     cache imgets
     imgets(l_image//"[0]","EXTEND", >& "dev$null")
     if(imgets.value =="T" ) mef=yes
     if (no == mef || (no == imaccess(l_image//"[1]") && no == imaccess(l_image//"[2]"))  ) {
          status = 4
     }
}
outname = l_image

end
