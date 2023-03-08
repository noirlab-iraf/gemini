# Copyright(c) 2005-2011 Association of Universities for Research in Astronomy, Inc.

procedure miptrans (inimages)

# This script takes an output file from "mipstack" with the four 
# individual waveplate polarimetry frames and converts these to 
# Stokes I, Q, and U images.
#
# Version:  September 13, 2005  KV write original script
#
#           The input file from "mipstack" has four extensions.  The 
#           output file has three extenstions, each 320 by 240 pixel images.
#           Extension 1 is the I image, extension 2 is the U image, and 
#           extension 3 is the Q image.  The I, Q, and U Stokes parameters 
#           describe the linear polarization.  See texts on polarimetry 
#           for a discussion of this.
#
#           Jan  10, 2006   KV added the WCS parameters to the image extension from the 
#                              primary header.
#

char    inimages    {prompt="Input stacked Michelle polarimetry image(s)"}   # OLDP-1-input-primary-single-prefix=s
char    outimages   {"",prompt="Output image(s)"}                       # OLDP-1-output
char    outpref     {"p",prompt="Prefix for output image(s)"}           # OLDP-4
char    rawpath     {"",prompt="Path for input raw images"}             # OLDP-4
bool    fl_register {no,prompt="Register images"}                       # OLDP-4
char    regions     {"[*,*]",prompt="Reference image regions used for registration [xregister]"}    # OLDP-4
char    logfile     {"",prompt="Logfile"}                               # OLDP-1
bool    verbose     {yes,prompt="Verbose?"}                             # OLDP-4
int     status      {0,prompt="Exit status: (0=good, >0=bad)"}          # OLDP-4
struct  *scanfile   {"",prompt="Internal use only"}                     # OLDP-4

begin

    char    l_inputimages, l_outputimages, l_filename, l_prefix, l_logfile
    char    l_rawpath, header
    bool    l_verbose, l_fl_register
    char    in[100], out[100], paramstr
    char    tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    char    tmpfile, tmpshifts
    char    l_regions, keyfound
    int     i, j, nimages, maximages, noutimages, nbad
    real    sigma

    l_verbose = verbose 
    l_inputimages = inimages
    l_outputimages = outimages
    l_logfile = logfile
    l_prefix = outpref
    l_fl_register = fl_register
    l_regions = regions
    l_rawpath = rawpath

    cache ("gloginit")

    nimages = 0
    maximages = 100
    status = 0

    tmpfile = mktemp ("tmpfile")
    tmpshifts = mktemp ("tmpshifts")

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "fl_register    = "//fl_register.p_value//"\n"
    paramstr += "regions        = "//regions.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "miptrans", "midir", paramstr, fl_append+,
        verbose=l_verbose )
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    #----------------------------------------------------------------------
    # Check the rawpath name for a final /
    if (substr (l_rawpath, (strlen(l_rawpath)), (strlen(l_rawpath))) != "/")
        l_rawpath = l_rawpath//"/"
    if ((l_rawpath == "/") || (l_rawpath == " "))
        l_rawpath = ""

    #-----------------------------------------------------------------------
    # Load up arrays of input name lists
    # This version handles both *s, / and commas in l_inputimages

    # check that list file exists
    if (substr(l_inputimages,1,1) == "@") {
        l_temp = substr (l_inputimages, 2, strlen(l_inputimages))
        if (!access(l_temp) && !access(l_rawpath//l_temp)) {
            glogprint (l_logfile, "miptrans", "status", type="error",
                errno=101, str="Input file "//l_temp//" not found.", verbose+)
            status = 101
            goto clean
        }
    }

    # Count the number of in images
    # First, generate the file list if needed

    if (stridx("*",l_inputimages) > 0) {
        files (l_inputimages, > tmpfile)
        l_inputimages = "@"//tmpfile
    }

    if (substr(l_inputimages,1,1) == "@")
        scanfile = substr (l_inputimages, 2, strlen(l_inputimages))
    else {
        files (l_inputimages, sort-, > tmpfile)
        scanfile = tmpfile
    }

    i = 0
    while ((fscan(scanfile,l_filename) != EOF) && (i <= 10)) {

        i = i+1

        if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) == ".fits")
            l_filename = substr (l_filename, 1, strlen(l_filename)-5)

        if (!imaccess(l_filename) && !imaccess(l_rawpath//l_filename)) {
            glogprint (l_logfile, "miptrans", "status", type="error",
                errno=101, str="Input image"//l_filename//" was not found.",
                verbose+)
            status = 1
            goto clean
        } else {
            nimages = nimages + 1
            if (nimages > maximages) {
                glogprint (l_logfile, "miptrans", "status", type="error",
                    errno=121, str="Maximum number of input images \
                    exceeded:"//maximages, verbose+)
                status = 1 
                goto clean
            }

            if ((l_rawpath == "") || (l_rawpath == " "))
                in[nimages] = l_filename
            else
                in[nimages] = l_rawpath//l_filename
        }
    }

    scanfile = ""
    delete (tmpfile, verify-, >& "dev$null")

    if (nimages == 0) {
        glogprint (l_logfile, "miptrans", "status", type="error", errno=121,
            str="No input images defined.", verbose+)
        status = 121
        goto clean
    }

    # Now, do the same counting for the out file

    nbad = 0
    noutimages = 0
    if ((l_outputimages != "") && (l_outputimages != " ")) {
        if (substr(l_outputimages,1,1) == "@")
            scanfile = substr (l_outputimages, 2, strlen(l_outputimages))
        else {
            if (stridx("*",l_outputimages) > 0) {
                files (l_outputimages, sort-) | \
                    match (".hhd", stop+, print-, metach-, > tmpfile)
                scanfile=tmpfile
            } else {
                files (l_outputimages, sort-, > tmpfile)
                scanfile = tmpfile
            }
        }

        while (fscan(scanfile,l_filename) != EOF) {
            noutimages = noutimages + 1
            if (noutimages > maximages) {
                glogprint (l_logfile, "miptrans", "status", type="error",
                    errno=121, str="Maximum number of output images \
                    exceeded:"//maximages, verbose+)
                status=121
                goto clean
            }
            if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) != ".fits") 
                l_filename = l_filename//".fits"
                
            out[noutimages] = l_filename
            if (imaccess(out[noutimages])) {
                glogprint (l_logfile, "miptrans", "status", type="error",
                    errno=102, str="Output image "//l_filename//" already \
                    exists.", verbose+)
                nbad += 1
            }
        }
        if (noutimages != nimages) {
            glogprint (l_logfile, "miptrans", "status", type="error",
                errno=121, str="Different number of in images ("//nimages//") \
                and out images ("//noutimages//")", verbose+)
            status = 121
            goto clean
        }

        scanfile = ""
        delete (tmpfile, verify-, >& "dev$null")

    } else {    # If prefix is to be used instead of filename

        print (l_prefix) | scan (l_prefix)
        if ((l_prefix == "") || (l_prefix == " ")) {
            glogprint (l_logfile, "miptrans", "status", type="error",
                errno=121, str="Neither output image name nor output prefix \
                is defined.", verbose+)
            status = 121
            goto clean
        }

        i = 1
        while (i <= nimages) {
            fparse (in[i])
            out[i] = l_prefix//fparse.root//".fits"

            if (imaccess(out[i])) {
                glogprint (l_logfile, "miptrans", "status", type="error", 
                    errno=102, str="Output image "//out[i]//" already exists.",
                    verbose+)
                nbad += 1
            }
            i = i+1
        }
    }

    if (nbad > 0) {
        glogprint (l_logfile, "miptrans", "status", type="error", errno=102,
            str=nbad//" image(s) already exist.", verbose+)
        status = 1
        goto clean
    }

    nbad = 0
    i = 1
    while (i <= nimages) {
        j = 1
        header = in[i]//"[0]"
        imgets (in[i]//"[0]", "MIPSTACK", >& "dev$null")
        if (imgets.value == "0") {
            glogprint (l_logfile, "miptrans", "status", type="warning", 
                errno=123, str="File "//in[i]//" has NOT been stacked using \
                mipstack.", verbose=l_verbose)
            goto nextimage
        }
        imgets (in[i]//"[0]", "MIPTRANS", >& "dev$null")
        if (imgets.value != "0") {
            glogprint (l_logfile, "miptrans", "status", type="warning", 
                errno=123, str="File "//in[i]//" has already been transformed \
                by MIPTRANS.", verbose=l_verbose)
            goto nextimage
        }
        tmp1 = mktemp ("tmp1")
        tmp2 = mktemp ("tmp2")
        tmp3 = mktemp ("tmp3")
        if (!l_fl_register) {
            imarith (in[i]//"[2]", "-", in[i]//"[4]", tmp1, title="",
                divzero=0., hparams="", pixtype="", calctype="", verbose-,
                noact-)
            # tmp1 is the Q stokes parameter, goes in the 3rd extension
            imarith (in[i]//"[1]", "-", in[i]//"[3]", tmp2, title="",
                divzero=0., hparams="", pixtype="", calctype="", verbose-,
                noact-)
            # tmp2 is the U stokes parameter, goes in the 2nd extension
            imcombine(in[i]//"[1],"//in[i]//"[2]"//in[i]//"[3],"//in[i]//"[4]",
                tmp3, combine="sum", headers="", bpmasks="", rejmasks="",
                expmask="", sigmas="", logfile="STDOUT", reject="none",
                project-, outtype="double", outlimits="", offsets="none",
                masktype="none", maskvalue=0., blank=0., scale="none",
                zero="none", weight="none", statsec="", expname="",
                lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1, mclip+,
                lsigma=3., hsigma=3., rdnoise="0.", gain="1.", snoise="0.",
                sigscale=0.1, pclip= -0.5, grow=0., >& "dev$null")
            imarith (tmp3, "/", 2.0, tmp3, title="", divzero=0., hparams="",
                pixtype="", calctype="", verbose-,noact-)
            # tmp3 is the I stokes parameter, goes in the first extension
            
        } else {
            tmp4 = mktemp ("tmp4")
            tmp5 = mktemp ("tmp5")
            tmp6 = mktemp ("tmp6")
            images.immatch.xregister (in[i]//"[2]", in[i]//"[1]", 
                regions=l_regions, shifts=tmpshifts, output=tmp4, 
                background="none", loreject=INDEF, hireject=INDEF, apodize=0.,
                filter="none", append+, records="", correlation="discrete",
                xwindow=11, ywindow=11, xcbox=11, ycbox=11,
                function="centroid", interp_type="poly5", interact-, xlag=0,
                ylag=0, dxlag=0, dylag=0) 
            images.immatch.xregister (in[i]//"[3]", in[i]//"[1]", 
                regions=l_regions, shifts=tmpshifts, output=tmp5,
                background="none", loreject=INDEF, hireject=INDEF, apodize=0.,
                filter="none", append+, records="", correlation="discrete",
                xwindow=11, ywindow=11, xcbox=11, ycbox=11,
                function="centroid", interp_type="poly5", interact-, xlag=0,
                ylag=0, dxlag=0, dylag=0) 
            images.immatch.xregister (in[i]//"[4]", in[i]//"[1]", 
                regions=l_regions, shifts=tmpshifts, output=tmp6,
                background="none", loreject=INDEF, hireject=INDEF, apodize=0.,
                filter="none", append+, records="", correlation="discrete",
                xwindow=11, ywindow=11, xcbox=11, ycbox=11,
                function="centroid", interp_type="poly5", interact-, xlag=0,
                ylag=0, dxlag=0, dylag=0) 
            imarith (tmp4, "-", tmp6, tmp1, title="", divzero=0., hparams="",
                pixtype="", calctype="", verbose-, noact-)
            # tmp1 is the Q stokes parameter, goes in the 3rd extension
            imarith (in[i]//"[1]", "-", tmp5, tmp2, title="", divzero=0.,
                hparams="", pixtype="", calctype="", verbose-, noact-)
            # tmp2 is the U stokes parameter, goes in the 2nd extension
            imcombine (in[i]//"[1],"//tmp4//","//tmp5//","//tmp6,
                tmp3, combine="sum", headers="", bpmasks="", rejmasks="",
                expmask="", sigmas="", logfile="STDOUT", reject="none",
                project-, outtype="double", outlimits="", offsets="none",
                masktype="none", maskvalue=0., blank=0., scale="none",
                zero="none", weight="none", statsec="", expname="",
                lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1, mclip+,
                lsigma=3., hsigma=3., rdnoise="0.", gain="1.", snoise="0.",
                sigscale=0.1, pclip= -0.5, grow=0., >& "dev$null")
            imarith (tmp3, "/", 2.0, tmp3, title="", divzero=0., hparams="",
                pixtype="", calctype="", verbose-, noact-)
            # tmp3 is the I stokes parameter, goes in the first extension

            # Clean-up of the tmp files used
            imdelete (tmp4, verify-, >& "dev$null")
            imdelete (tmp5, verify-, >& "dev$null")
            imdelete (tmp6, verify-, >& "dev$null")
        }
        wmef (tmp3, out[i], extname="SCI", phu=header, verbose-, >& "dev$null")
        fxinsert (tmp2, out[i]//"["//j//"]", "", verbose=l_verbose,
            >& "dev$null")
        j = j+1
        fxinsert (tmp1, out[i]//"["//j//"]", "", verbose=l_verbose,
            >& "dev$null")
        imdelete (tmp1//","//tmp2//","//tmp3, verify-, >& "dev$null")

        #-----------------------------------
        ###### Update header keywords ######
        #-----------------------------------

        # Copy WCS information to the output image
        
        keyfound = ""
        hselect (header, "CTYPE", yes) | scan (keyfound)
        if (keyfound != "") {
            for (k=1; k<=3; k=k+1) {
                gemhedit (out[i]//"["//str(k)//"]", "WCSAXES", 2,
                    comment="Number of WCS axes in the image", delete-)
                gemhedit (out[i]//"["//str(k)//"]", "CTYPE1", keyfound,
                    comment="R.A. in tangent plane projection", delete-)
            }
        }
        keyfound = ""
        hselect (header, "CRPIX1", yes) | scan (keyfound)
        if (keyfound != "") {
            for (k=1; k<=3; k=k+1) {
                gemhedit (out[i]//"["//str(k)//"]", "CRPIX1", keyfound,
                    comment="Ref pix of axis 1", delete-)
            }
        }
        keyfound = ""
        hselect (header, "CRVAL1", yes) | scan (keyfound)
        if (keyfound != "") {
            for (k=1; k<=3; k=k+1) {
                gemhedit (out[i]//"["//str(k)//"]", "CRVAL1", keyfound,
                    comment="RA at Ref pix in decimal degrees", delete-)
            }
        }
        keyfound = ""
        hselect (header, "CTYPE2", yes) | scan (keyfound)
        if (keyfound != "") {
            for (k=1; k<=3; k=k+1) {
                gemhedit (out[i]//"["//str(k)//"]", "CTYPE2", keyfound,
                    comment="DEC. in tangent plane projection", delete-)
            }
        }
        keyfound = ""
        hselect (header, "CRPIX2", yes) | scan (keyfound)
        if (keyfound != "") {
            for (k=1; k<=3; k=k+1) {
                gemhedit (out[i]//"["//str(k)//"]", "CRPIX2", keyfound,
                    comment="Ref pix of axis 2", delete-)
            }
        }
        keyfound = ""
        hselect (header, "CRVAL2", yes) | scan (keyfound)
        if (keyfound != "") {
            for (k=1; k<=3; k=k+1) {
                gemhedit (out[i]//"["//str(k)//"]", "CRVAL2", keyfound,
                    comment="DEC at Ref pix in decimal degrees", delete-)
            }
        }
        keyfound = ""
        hselect (header, "CD1_1", yes) | scan (keyfound)
        if (keyfound != "") {
            for (k=1; k<=3; k=k+1) {
                gemhedit (out[i]//"["//str(k)//"]", "CD1_1", keyfound,
                    comment="WCS matrix element 1 1", delete-)
            }
        }
        keyfound = ""
        hselect (header, "CD1_2", yes) | scan (keyfound)
        if (keyfound != "") {
            for (k=1; k<=3; k=k+1) {
                gemhedit (out[i]//"["//str(k)//"]", "CD1_2", keyfound,
                    comment="WCS matrix element 1 2", delete-)
            }
        }
        keyfound = ""
        hselect (header, "CD2_1", yes) | scan (keyfound)
        if (keyfound != "") {
            for (k=1; k<=3; k=k+1) {
                gemhedit (out[i]//"["//str(k)//"]", "CD2_1", keyfound,
                    comment="WCS matrix element 2 1", delete-)
            }
        }
        keyfound = ""
        hselect (header, "CD2_2", yes) | scan (keyfound)
        if (keyfound != "") {
            for (k=1; k<=3; k=k+1) {
                gemhedit (out[i]//"["//str(k)//"]", "CD2_2", keyfound,
                    comment="WCS matrix element 2 2", delete-)
            }
        }
        keyfound = ""
        hselect (header, "RADECSYS", yes) | scan (keyfound)
        if (keyfound != "") {
            for (k=1; k<=3; k=k+1) {
                gemhedit (out[i]//"["//str(k)//"]", "RADECSYS", keyfound,
                    comment="R.A./DEC. coordinate system reference",
                    delete-)
            }
        }

        for (k=1; k<=3; k=k+1) {
            gemhedit (out[i]//"["//str(k)//"]", "EXTNAME", "SCI",
                comment="Extension name", delete-)
            gemhedit (out[i]//"["//str(k)//"]", "EXTVER", k,
                comment="Extension version", delete-)
        }

        # Time stamps
        gemdate ()
        gemhedit (out[i]//"[0]", "MIPTRANS", gemdate.outdate,
            "UT Time stamp for MIPTRANS", delete-)
        gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)

nextimage:
        i = i+1

    }

clean:
    scanfile = ""
    imdelete (tmpfile, verify-, >& "dev$null")
    delete ("tmpshifts*", verify-, >& "dev$null")
    
    if (status == 0)
        glogclose (l_logfile, "miptrans", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "miptrans", fl_success-, verbose=l_verbose)

exitnow:
    ;

end
