# Copyright(c) 2002-2012 Association of Universities for Research in Astronomy, Inc.

procedure gscrmask (inimage, outimage)

# Create a bad pixel/CR image by thresholding an input image
# (typically called by gscrrej.cl). Doesn't handle MEF directly.
#
# If the optional noise image already exists, it is used in place of
# calculated noise estimates, otherwise it is created by this task. 
# Hence it is possible to call gscrmask repeatedly without recalculating
# the noise, which is rather slow.
#
# Date    Sept 20, 2002 JT v1.4 release
#         Aug 26, 2003  KL IRAF2.12 - new parameters
#                             imstat: nclip, lsigma, usigma. cache

string  inimage     {prompt="Input image with real features removed"}
string  outimage    {prompt="Output pixel mask image, flagging CRs etc."}
real    nsigma      {"10", min=1.0, prompt="Sigma clipping threshold"}
int     boxsize     {"11", min=3, prompt="Box size (pix) for local noise est."}
string  varimage    {"", prompt="Optional noise image to use or create\n"}

bool    verbose     {no, prompt="Verbose output?"}
string  logfile     {"", prompt="Log file"}
int     status      {0, prompt="Exit status (0=good)"}

begin

    string  l_inimage, l_outimage, l_varimage, l_logfile, tstr
    real    l_nsigma
    int     l_boxsize
    bool    l_verbose

    string  timage, timage2, tvarimage
    real    imav
    real    varlim
    int     n
    bool    usevarim, mkvarim

    # Keep imgets parameters from changing by outside world:
    cache ("gimverify", "gemlogname")

    # Read the parameters:
    l_inimage = inimage
    l_outimage = outimage
    l_nsigma = nsigma
    l_boxsize = boxsize
    l_varimage = varimage
    l_verbose = verbose
    l_logfile = logfile

    # Default exit status:
    status = 1

    # Request temp image names:
    timage = mktemp("tmpim")
    timage2 = mktemp("tmpim")
    tvarimage = mktemp("tmpim")

    # Determine the logfile name:
    gemlogname (logpar=l_logfile, package="gmos")
    l_logfile = gemlogname.logname

    # Warn if logfile name had to revert to default:
    if (gemlogname.status == 1) {
        printlog ("WARNING - GSCRMASK: both gscrmask.logfile and gmos.logfile \
            are empty;", l_logfile, verbose+)
        printlog ("                    using "//l_logfile, l_logfile, verbose+)
    } else if (gemlogname.status == 2)
        printlog ("WARNING - GSCRMASK: bad logfile name, "//\
            gemlogname.logpar//" - using "//l_logfile, l_logfile, verbose+)

    # Check that input file exists and is not MEF:
    gimverify (l_inimage)
    if (gimverify.status == 1) {
        printlog ("ERROR - GSCRMASK: "//l_inimage//" does not exist",
            l_logfile, verbose+)
        goto error
    } else if (gimverify.status == 0) {
        printlog ("ERROR - GSCRMASK: task does not handle MEF files directly",
            l_logfile, verbose+)
        goto error
    }

    # Check that the output filename is valid:
    if (l_outimage == "" || (stridx(" ", l_outimage) > 0)) {
        printlog ("ERROR - GSCRMASK: invalid output filename", l_logfile, 
            verbose+)
        goto error
    }

    # Check that output file does not exist, before wasting time processing:
    gimverify(l_outimage)
    if (gimverify.status != 1) {
        printlog ("ERROR - GSCRMASK: "//l_outimage//" already exists",
            l_logfile, verbose+)
        goto error
    }

    # Check whether noise image specified & whether it exists:
    print (l_varimage) | scanf ("%s ", tstr)
    if (tstr != "") {
        usevarim = yes
        gimverify (l_varimage)
        if (gimverify.status == 0) {
            printlog ("ERROR - GSCRMASK: "//l_varimage//" is a MEF file",
                l_logfile, verbose+)
            goto error
        } else if (gimverify.status == 1)
            mkvarim = yes
        else
            mkvarim = no
    } else {
        usevarim = no
        mkvarim = no
    }

    # Compute image background as mode (quicker than iteratively estimating an
    # accurate mean).
    imstat (l_inimage, fields="mode", lower=INDEF, upper=INDEF, nclip=0,
        lsigma=INDEF, usigma=INDEF, binwidth=0.1, format-, cache-) |\
        scanf ("%f", imav)

    # Subtract image "background" level for generality (but expect ~0):
    imarith (l_inimage, "-", imav, timage, title="", divzero=0, hparams="",
        pixtype="", calctype="", verbose-, noact-)

    # Square the pixel values to get square residuals:
    imexpr ("a**2", timage2, timage, dims="auto", intype="auto", 
        outtype="auto", refim="auto", bwidth=0, btype="constant", bpixval=0.0,
        rangecheck+, verbose-, exprdb="none")
    imdelete (timage, go_ahead+, verify-)

    # Estimate variance or use supplied estimate:
    if (!usevarim || mkvarim) {
        #printlog("Estimating local noise in ~"//l_boxsize//"x"//l_boxsize//"\
        #    box", l_logfile, l_verbose)

        if (!mkvarim) l_varimage = tvarimage

        # Median filter the local residuals to get robust noise estimates:
        median (timage2, l_varimage, l_boxsize, l_boxsize, zloreject=INDEF,
            zhireject=INDEF, boundary="constant", constant=0.0, verbose-)
        # NB. tried replacing this call by fmedian, but only saved ~1 minute
        # out of 17 and the results were worse.

        #if (mkvarim)
        #    printlog("Saved variance estimate for later ("//l_varimage//")",
        #        l_logfile, l_verbose)
    }
    #else {
    #  printlog("Using previous variance estimate ("//l_varimage//")",
    #           l_logfile, l_verbose)
    #}

    #printlog("Creating mask with "//l_nsigma//" sigma clipping",
    #    l_logfile, l_verbose)

    # Make mask by clipping residuals with ref. to noise image:
    imexpr ("(a < -b*c || a > b*c)", l_outimage, timage2, l_nsigma*l_nsigma,
        l_varimage, dims="auto", intype="auto", outtype="auto", refim="auto",
        bwidth=0, btype="constant", bpixval=0.0, rangecheck+,verbose-,
        exprdb="none")

    # Set exit status good:
    status = 0

error:
    # Jump here directly in the event of an error:

    # Clear up
    if (imaccess(timage)) imdelete (timage, go_ahead+, verify-)
    if (imaccess(timage2)) imdelete (timage2, go_ahead+, verify-)
    if (imaccess(tvarimage)) imdelete (tvarimage, go_ahead+, verify-)

end
