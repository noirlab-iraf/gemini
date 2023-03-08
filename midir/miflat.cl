# Copyright(c) 2003-2011 Association of Universities for Research in Astronomy, Inc.

procedure miflat(inimage1,inimage2,outimage)

# This routine makes a flat field file from two images, which must be 
# the output of either "tprepare" or "mprepare".
#
# Version:  Sept 29, 2003  KV wrote the first version.
#           Oct   3, 2003  KV made small syntax changes in the parameters
#	    Oct 29, 2003  KL IRAF2.12 - new parameters
#			       imstat: nclip, lsigma, usigma, cache

char    inimage1    {"",prompt="First input T-ReCS or Michelle flat-field image"}
char    inimage2    {"",prompt="Second input T-ReCS or Michelle flat-field image"}
char    rawpath     {"",prompt="Path for input raw images"}
char    outimage    {"",prompt="Output image name"}
char    logfile     {"",prompt="Logfile name"}
bool    verbose     {yes,prompt="Verbose logging yes/no?"}
int     status      {0,prompt="Exit error status: (0=good, >0=bad)"}
struct* scanfile    {"",prompt="Internal use only"}

begin

    char    l_inputimage1 = ""
    char    l_inputimage2 = ""
    char    l_outputimage = ""
    char    l_rawpath = ""
    char    l_logfile = ""
    bool    l_verbose

    char    tmpimage, header1, header2, tmpstring, tmpsci
    char    paramstr
    int     naxis[4], iaxis, inst1, inst2, i, j, junk, idx
    real    mean1, mean2

    junk = fscan ( inimage1, l_inputimage1 )
    junk = fscan ( inimage2, l_inputimage2 )
    junk = fscan ( outimage, l_outputimage )
    junk = fscan ( rawpath, l_rawpath )
    junk = fscan ( logfile, l_logfile )
    l_verbose = verbose

    cache ("imgets", "gloginit", "gemdate")

    # Initialize
    status = 0
    
    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimage1       = "//inimage1.p_value//"\n"
    paramstr += "inimage2       = "//inimage2.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "outimage       = "//outimage.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value
    
    # Assign a logfile name, if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "miflat", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Check the rawpath name for a final /
    if (substr(l_rawpath,(strlen(l_rawpath)),(strlen(l_rawpath))) != "/")
        l_rawpath=l_rawpath//"/"
    if (l_rawpath=="/" || l_rawpath==" ")
        l_rawpath=""

    # Check inputs
    if (!imaccess(l_rawpath//l_inputimage1//"[0]")) {
        glogprint (l_logfile, "miflat", "status", type="error", errno=101,
            str="Input image "//l_rawpath//l_inputimage1//" was not found.",
            verbose+)
        status = 1
        goto exit
    } else {
        idx = strldx (".", l_inputimage1)
        if ( substr(l_inputimage1, idx+1, strlen(l_inputimage1)) == "fits" ) {
            l_inputimage1 = substr(l_inputimage1, 1, idx-1)
        }
    }
    if (!imaccess(l_rawpath//l_inputimage2//"[0]")) {
        glogprint (l_logfile, "miflat", "status", type="error", errno=101,
            str="Input image "//l_rawpath//l_inputimage2//" was not found.",
            verbose+)
        status=1
        goto exit
    } else {
        idx = strldx (".", l_inputimage2)
        if ( substr(l_inputimage2, idx+1, strlen(l_inputimage2)) == "fits" ) {
            l_inputimage2 = substr(l_inputimage2, 1, idx-1)
        }
    }

    if (imaccess(l_outputimage)) {
        glogprint (l_logfile, "miflat", "status", type="error", errno=102,
            str="Output image "//l_outputimage//" already exists.",verbose+)
        status = 1
        goto exit
    } else {
        idx = strldx (".", l_outputimage)
        if ( substr(l_outputimage, idx+1, strlen(l_outputimage)) == "fits" ) {
            l_outputimage = substr(l_outputimage, 1, idx-1)
        }
    }


    header1 = l_rawpath//l_inputimage1//"[0]"
    header2 = l_rawpath//l_inputimage2//"[0]"

    imgets (header1,"INSTRUME", >& "dev$null")

    inst1 = 0
    if (imgets.value == "TReCS")
        inst1 = 1
    else if (imgets.value == "Michelle")
        inst1 = 2

    imgets (header2,"INSTRUME", >& "dev$null")
    inst2 = 0
    if (imgets.value == "TReCS")
        inst2=1
    else if (imgets.value == "Michelle")
        inst2=2

    if (inst1 == 0) {
        glogprint (l_logfile, "miflat", "status", type="error", errno=123,
            str="Input image "//l_inputimage1//" does not seem to be a T-ReCS \
            or Michelle file.",verbose+)
        status = 1
        goto exit
    }

    if (inst2 == 0) {
        glogprint (l_logfile, "miflat", "status", type="error", errno=123,
            str="Input image "//l_inputimage2//" does not seem to be a T-ReCS \
            or Michelle file.",verbose+)
        status = 1
        goto exit
    }

    if (inst1 != inst2) {
        glogprint (l_logfile, "miflat", "status", type="error", errno=121,
            str="Input images "//l_inputimage1//" and "//l_inputimage2//\
            " do not come from the same instrument!.",verbose+)
        status = 1
        goto exit
    }

    if (inst1 == 1) {
        imgets (header1, "TPREPARE", >& "dev$null")
        if (imgets.value == "0") {
            glogprint (l_logfile, "miflat", "status", type="error", errno=121,
                str="Input image "//l_inputimage1//" has not been processed \
                with TPREPARE.",verbose+)
            status = 1
            goto exit
        }
        imgets (header2, "TPREPARE", >& "dev$null")
        if (imgets.value == "0") {
            glogprint (l_logfile, "miflat", "status", type="error", errno=121,
                str="Input image "//l_inputimage2//" has not been processed \
                with TPREPARE.",verbose+)
            status = 1
            goto exit
        }
    }

    if (inst1 == 2) {
        imgets (header1, "MPREPARE", >& "dev$null")
        if (imgets.value == "0") {
            glogprint (l_logfile, "miflat", "status", type="error", errno=121,
                str="Input image "//l_inputimage1//" has not been processed \
                with MPREPARE.",verbose+)
            status = 1
            goto exit
        }
        imgets (header2, "MPREPARE", >& "dev$null")
        if (imgets.value == "0") {
            glogprint (l_logfile, "miflat", "status", type="error", errno=121,
                str="Input image "//l_inputimage2//" has not been processed \
                with MPREPARE.",verbose+)
            status = 1
            goto exit
        }
    }

    header1 = l_rawpath//l_inputimage1//"[1]"
    if (imaccess(header1)) {
        imgets(header1,"i_naxis")
        iaxis = int(imgets.value)
        if (iaxis < 2 || iaxis > 4) {
            glogprint (l_logfile, "miflat", "status", type="error", errno=123,
                str="Input image "//l_inputimage1//" has a bad number of \
                dimensions ("//str(iaxis)//") for extension 1.",verbose+)
            status = 1
            goto exit
        }
        naxis[1] = 0
        naxis[2] = 0
        naxis[3] = 0
        naxis[4] = 0

        for (j=1; j <= iaxis; j+=1) {
            tmpstring = "naxis"//str(j)
            imgets(header1,tmpstring)
            naxis[j] = int(imgets.value)
            if (naxis[j] < 0) {
                glogprint (l_logfile, "miflat", "status", type="error",
                    errno=123, str="Image "//l_inputimage1//" has negative \
                    image dimension(s) in extension 1.",verbose+)
                status = 1
                goto exit
            }
        }
        if (iaxis > 2 || naxis[1] != 320 || naxis[2] != 240) {
            glogprint (l_logfile, "miflat", "status", type="error", errno=123,
                str="Image "//l_inputimage1//" has bad image dimension(s) \
                in extension 1.",verbose+)
            status = 1
            goto exit
        }
    } else {
        glogprint (l_logfile, "miflat", "status", type="error", errno=101,
            str="A problem occured when reading extension 1 of input image "//\
            l_inputimage1//".",verbose+)
        status = 1
        goto exit
    }

    imstat(header1,field="midpt",lower=INDEF,upper=INDEF,nclip=0,lsigma=INDEF,
        usigma=INDEF,binwidth=0.1,format-,cache-) | scanf("%f",mean1)

    header2 = l_rawpath//l_inputimage2//"[1]"
    if (imaccess(header2)) {
        imgets(header2,"i_naxis")
        iaxis = int(imgets.value)
        if (iaxis < 2 || iaxis > 4) {
            glogprint (l_logfile, "miflat", "status", type="error", errno=123,
                str="Input image "//l_inputimage2//" has a bad number of \
                dimensions ("//str(iaxis)//") for extension 1.",verbose+)
            status = 1
            goto exit
        }
        naxis[1] = 0
        naxis[2] = 0
        naxis[3] = 0
        naxis[4] = 0

        for (j=1; j <= iaxis; j+=1) {
            tmpstring = "naxis"//str(j)
            imgets(header2,tmpstring)
            naxis[j] = int(imgets.value)
            if (naxis[j] < 0) {
                glogprint (l_logfile, "miflat", "status", type="error",
                    errno=123, str="Image "//l_inputimage2//" has negative \
                    image dimension(s) in extension 1.",verbose+)
                status = 1
                goto exit
            }
        }
        if (iaxis > 2 || naxis[1] != 320 || naxis[2] != 240) {
            glogprint (l_logfile, "miflat", "status", type="error", errno=123,
                str="Image "//l_inputimage2//" has bad image dimension(s) \
                in extension 1.",verbose+)
            status = 1
            goto exit
        }
    } else {
        glogprint (l_logfile, "miflat", "status", type="error", errno=101,
            str="A problem occured when reading extension 1 of input image "//\
            l_inputimage2//".",verbose+)
        status = 1
        goto exit
    }

    imstat (header2, field="midpt", lower=INDEF, upper=INDEF, nclip=0,
        lsigma=INDEF, usigma=INDEF, binwidth=0.1, format-, cache-) | \
        scanf ("%f",mean2)
#    printf("mean value 1 = %f\n",mean2)

    tmpimage = mktemp("tmpimage")
	tmpsci = mktemp("tmpsci")
    if (mean1 > mean2)
        imarith(header1//"[*,*]","-",header2//"[*,*]",tmpimage,verbose-)
    else
        imarith(header2//"[*,*]","-",header1//"[*,*]",tmpimage,verbose-)

    imstat(tmpimage,field="midpt",lower=INDEF,upper=INDEF,nclip=0,lsigma=INDEF,
        usigma=INDEF,binwidth=0.1,format-,cache-) | scanf("%f",mean1)

    if (mean1 == 0.) {
        glogprint (l_logfile, "miflat", "status", type="error", errno=1,
            str="A problem occured when making the flat field image (median \
            value is zero).",verbose+)
        status = 1
        goto exit
    }

    imarith (tmpimage, "/", mean1, tmpsci, verbose-)

	# Add PHU from first image to create MEF, name/number extensions
    fxcopy (l_rawpath//l_inputimage1//".fits", l_outputimage//".fits", 
        groups="0", new+, ver-)
    fxinsert (tmpsci//".fits", l_outputimage//".fits[0]", group="", verbose-)
    gemhedit (l_outputimage//"[0]", "EXTVER", 0, "Extension version")
	gemhedit (l_outputimage//"[1]", "EXTVER", 1, "Extension version")
	gemhedit (l_outputimage//"[1]", "EXTNAME", "SCI", "Extension name")

    gemdate()
    gemhedit (l_outputimage//"[0]", "MIFLAT", gemdate.outdate,
        "UT Time stamp for MIFLAT", delete-)
    gemhedit (l_outputimage//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)
    
    
exit:
    delete ("tmpimage*.fits", ver-, >& "dev$null")
	delete ("tmpsci*.fits", ver-, >& "dev$null")

    if (status == 0)
        glogclose (l_logfile, "miflat", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "miflat", fl_success-, verbose=l_verbose)

exitnow:
    ;

end
