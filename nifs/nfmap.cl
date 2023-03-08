# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.
#
# nfmap -- Extract a channel map from a NIFS data cube
#
# Original author: Kevin Volk  [15Mar2006]

procedure nfmap (inimage, wavelength, npix)

char    inimage     {prompt="Input NIFS datacube file"}   # OLDP-1-input-primary
real    wavelength  {0.0, prompt="Wavelength of line center"}
int     npix        {9, prompt="Number of pixels to extract (normally odd)"}
real    velocity    {0.0, prompt="Optional velocity offset for region to extract (km/s)"}
char    outimage    {"", prompt="Output image name"}
char    outpref     {"m", prompt="Output prefix (if no output image name is given)"}
bool    fl_cont     {yes, prompt="Subtract continuum from image"}
char    cont1       {"-5:-8", prompt="Continuum range 1 (pixels, from line CENTER)"}
char    cont2       {"5:8", prompt="Continuum range 2 (pixels, from line CENTER)"}
bool    fl_interp   {no, prompt="Interpolate continuum rather than average"}
char    logfile     {"", prompt="Log file name"}
bool    verbose     {yes, prompt="verbose logging?"}
int     status      {0, prompt="Exit error status: (0=good, >0=bad)"}
struct  *scanfile   {"", prompt="Internal use only"}

begin

    char    l_inimage = ""
    char    l_outimage = ""
    char    l_outpref = ""
    char    l_cont1 = ""
    char    l_cont2 = ""
    char    l_logfile = ""
    int     l_npix
    real    l_wavelength, l_velocity
    bool    l_fl_interp, l_fl_cont, l_verbose
    
    char    tmpwork, tmpstring, tmpcont1, tmpcont2, tmplist
    char    tmpfinal, tmpspectrum, tmpcursor
    char    tmpfile, tmpinimg, tmpoutimg
    char    paramstr, errmsg, outputstr, keyfound
    int     junk, ninimages, noutimages
    int     i, j, k, l, iaxis, naxis[3], imin, imax, icont[4], icenter
    real    pixel0, wavel0, delwl, pixel1, xlow, xhigh, scale, ymin, ymax
    real    xval1, xval2, yval1, yval2
    struct  l_struct
    
    junk = fscan (inimage, l_inimage)
    l_wavelength    = wavelength
    l_npix          = npix
    l_velocity      = velocity
    junk = fscan (outimage, l_outimage)
    junk = fscan (outpref, l_outpref)
    l_fl_cont       = fl_cont
    junk = fscan (cont1, l_cont1)
    junk = fscan (cont2, l_cont2)
    l_fl_interp     = fl_interp
    junk = fscan (logfile, l_logfile)
    l_verbose       = verbose
    
    cache ("gloginit", "gemdate")
    
    # Initialize
    status = 0
    
    tmpfile = mktemp ("tmpfile")
    tmpoutimg = mktemp ("tmpoutimg")
    tmpinimg = mktemp ("tmpinimg")

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimage        = "//inimage.p_value//"\n"
    paramstr += "wavelength     = "//wavelength.p_value//"\n"
    paramstr += "npix           = "//npix.p_value//"\n"
    paramstr += "velocity       = "//velocity.p_value//"\n"
    paramstr += "outimage       = "//outimage.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "fl_cont        = "//fl_cont.p_value//"\n"
    paramstr += "cont1          = "//cont1.p_value//"\n"
    paramstr += "cont2          = "//cont2.p_value//"\n"
    paramstr += "fl_interp      = "//fl_interp.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "nfmap", "nifs", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Check input
    gemextn (l_inimage, check="exists,mef", process="none", index="",
        extname="", extver="", ikparam="", omit="extension", replace="",
        outfile=tmpinimg, logfile=l_logfile, verbose=l_verbose)
    ninimages = gemextn.count
    
    if ((gemextn.fail_count > 0) || (ninimages == 0) || \
        (ninimages > 1)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" input cube(s) not found."
            status = 101
        } else if (ninimages == 0) {
            errmsg = "No input cube defined."
            status = 121
        } else if (ninimages > 1) {
            errmsg = "Task works on only one input cube at a time."
            status = 121
        }
        
        glogprint (l_logfile, "nfmap", "status", type="error", errno=status,
            str=errmsg, verbose+)
        goto clean
    } else {
        scanfile = tmpinimg
        junk = fscan (scanfile, l_inimage)
    }

    # Check output
    if (l_outimage != "")
        outputstr = l_outimage
    else if (l_outpref != "") {
        gemextn (l_inimage, check="", process="none", index="", extname="",
            extversion="", ikparams="", omit="path", replace="",
            outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
        outputstr = l_outpref//"@"//tmpoutimg
        delete (tmpfile, verify-, >& "dev$null")
    } else {
        status = 121
        glogprint (l_logfile, "nfmap", "status", type="error", errno=status,
            str="Neither output image name nor output prefix defined.",
            verbose+)
        goto clean
    }

    gemextn (outputstr, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension", replace="",
        outfile=tmpfile, logfile=l_logfile, verbose=l_verbose)
    delete (tmpoutimg, verify-, >& "dev$null")
    gemextn ("@"//tmpfile, check="absent", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
    noutimages = gemextn.count
    delete (tmpfile, verify-, >& "dev$null")
    
    if ((gemextn.fail_count > 0) || (noutimages == 0) || \
        (noutimages != ninimages)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" image(s) already exist(s)."
            status = 102
        } else if (noutimages == 0) {
            errmsg = "No output images defined."
            status = 121
        } else if (noutimages != ninimages) {
            errmsg = "Different number of input images ("//ninimages//") and \
                output images ("//noutimages//")."
            status = 121
        }
        
        glogprint (l_logfile, "nfmap", "status", type="error", errno=status,
            str=errmsg, verbose+)
        goto clean
    } else {
        scanfile = tmpoutimg
        junk = fscan (scanfile, l_outimage)
    }
    delete (tmpinimg, verify-, >& "dev$null")
    delete (tmpoutimg, verify-, >& "dev$null")
    
    # Check input for required processing and keywords
    keyfound = ""
    hselect (l_inimage//"[0]", "INSTRUME", yes) | scan (keyfound)
    if (keyfound != "NIFS") {
        glogprint (l_logfile, "nfmap", "status", type="error", errno=121,
            str="The input cube is not NIFS data.", verbose+)
        goto clean
    }

    keyfound = ""
    hselect (l_inimage//"[0]", "NSTRANSF", yes) | scan (keyfound)
    if (keyfound == "") {
        glogprint (l_logfile, "nfmap", "status", type="error", errno=121,
            str="The input cube has not been processed with NSTRANSFORM.",
            verbose+)
        goto clean
    }

    keyfound = ""
    hselect (l_inimage//"[1]", "CTYPE3", yes) | scan (keyfound)
    if (keyfound != "LINEAR") {
        glogprint (l_logfile, "nfmap", "status", type="error", errno=121,
            str="The spectra are not linearized.", verbose+)
        goto clean
    }

    keyfound = ""
    hselect (l_inimage//"[1]", "CRVAL3", yes) | scan (keyfound)
    if (keyfound == "")  status = 131
    else                wavel0 = real(keyfound)
    
    keyfound = ""
    hselect (l_inimage//"[1]", "CD3_3", yes) | scan (keyfound)
    if (keyfound == "")  status = 131
    else                delwl = real(keyfound)
    
    keyfound = ""
    hselect (l_inimage//"[1]", "CRPIX3", yes) | scan (keyfound)
    if (keyfound == "")  status = 131
    else                pixel0 = real(keyfound)
    
    if (status != 0) {
        errmsg = "Could not obtain the wavelength parameters from the header."
        glogprint (l_logfile, "nfmap", "status", type="error", errno=status,
            str=errmsg, verbose+)
        goto clean
    }

    scale = 299792.458 * delwl / l_wavelength
    printf ("The pixel scale of %6.3f Angstroms corresponds to %8.3f km/s \
        at %10.3f Angstroms\n", delwl, scale, l_wavelength) | scan (l_struct)
    glogprint (l_logfile, "nfmap", "status", type="string", str=l_struct,
        verbose=l_verbose)

    if (l_fl_cont) {

        k = stridx (":", l_cont1)

        if ((k < 2) || (k == strlen(l_cont1))) {
            glogprint (l_logfile, "nfmap", "status", type="error", errno=121,
                str="Bad cont1 parameter: "//l_cont1, verbose+)
            goto clean
        } else {
            l_cont1 = substr (l_cont1, 1, k-1)//" "//\
                substr (l_cont1, k+1, strlen(l_cont1))
        }
        k = stridx (":", l_cont2)

        if ((k < 2) || (k == strlen(l_cont2))) {
            glogprint (l_logfile, "nfmap", "status", type="error", errno=121,
                str="Bad cont2 parameter: "//l_cont2, verbose+)
            goto clean
        } else {
            l_cont2 = substr (l_cont2, 1, k-1)//" "//\
                substr (l_cont2, k+1, strlen(l_cont2))
        }
    }

    keyfound = "0"
    hselect (l_inimage//"[1]", "i_naxis", yes) | scan (keyfound)
    iaxis = int(keyfound)

    if (iaxis < 2 || iaxis > 3) {
        tmpstring = "Image extension 1 has NAXIS = "//str(iaxis)//\
            ".\nThe expected number is 3."
        glogprint (l_logfile, "nfmap", "status", type="error", errno=121,
            str=tmpstring, verbose+)
        goto clean
    } else {
        naxis[1] = 0
        naxis[2] = 0
        naxis[3] = 0

        for (k=1; k <= iaxis; k+=1) {
            tmpstring = "naxis"//str(k)
            
            keyfound = ""
            hselect (l_inimage//"[1]", tmpstring, yes) | scan (keyfound)
            if (keyfound == "") {
                errmsg = "Axis "//k//" length not found."
                status = 132
                glogprint (l_logfile, "nfmap", "status", type="error",
                    errno=status, str=errmsg, verbose+)
                goto clean
            } else
                naxis[k] = int(keyfound)
        }
    }

    pixel1 = ((l_wavelength-wavel0) / delwl) + pixel0
    icenter = int(pixel1+0.5)

    if (l_velocity != 0.) {
        scale = sqrt((1.+l_velocity/299792.458) / (1.-l_velocity/299792.458))
        l_wavelength = l_wavelength * scale
        scale = ((l_wavelength-wavel0) / delwl) + pixel0
        scale = scale - pixel1
        if (abs(scale) > 100.) {
            tmpstring = "Specified velocity shift is too large (> 100 pixels)"
            glogprint (l_logfile, "nfmap", "status", type="error", errno=121,
                str=tmpstring, verbose+)
            goto clean
        }
        pixel1 = ((l_wavelength-wavel0) / delwl) + pixel0
        i = int(pixel1+0.5)
        printf ("The shifted wavelength specified, %10.3f Angstroms, \
            corresponds to pixel %7.2f.\n", l_wavelength,pixel1) | \
            scan(l_struct)
        glogprint (l_logfile, "nfmap", "status", type="string", str=l_struct,
            verbose+)
    } else {
        printf ("The wavelength specified, %10.3f Angstroms, corresponds to \
            pixel %7.2f.\n",l_wavelength,pixel1) | scan(l_struct)
        glogprint (l_logfile, "nfmap", "status", type="string", str=l_struct,
            verbose+)
    }
    i = int(pixel1+0.5)

    if ((pixel1 < 4) || (pixel1 > naxis[3]-3)) {
        glogprint (l_logfile, "nsmap", "status", type="error", errno=121,
            str="The wavelength requested is out of range.", verbose+)
        goto clean
    }

    glogprint (l_logfile, "nsmap", "status", type="string",
        str="Pixel nearest specified wavelength: "//str(i), verbose=l_verbose)

    tmpwork = mktemp ("tmpspectrum")
    tmpspectrum = mktemp ("tmpspectrum")
    blkavg (l_inimage//"[1]", tmpwork, 100, 100, 1, >& "dev$null")
    imcopy (tmpwork//"[1,1,*]", tmpspectrum, verbose-, >& "dev$null")
    imdelete (tmpwork, verify-, >& "dev$null")

    j = l_npix / 2
    if (l_npix == (j*2)) {
        imin = i - j+1
        imax = i + j
    } else {
        imin = i - j
        imax = i + j
    }

    if (l_fl_cont) {
        l = fscanf (l_cont1, "%d %d", j, k)
        icont[1] = j
        icont[2] = k
        l = fscanf (l_cont2, "%d %d", j, k)
        icont[3] = j
        icont[4] = k
        xlow = (real(icont[1])+real(icont[2])) / 2.
        xhigh = (real(icont[3])+real(icont[4])) / 2.
        scale = xlow * xhigh
        if (scale >= 0.) {
            glogprint (l_logfile, "nsmap", "status", type="error", errno=121,
                str="The continuum regions do not bracket the line.", verbose+)
            goto clean
        }

        icont[1] = icont[1] + icenter
        icont[2] = icont[2] + icenter
        icont[3] = icont[3] + icenter
        icont[4] = icont[4] + icenter

        if (icont[1] > icont[2]) {
            k = icont[1]
            icont[1] = icont[2]
            icont[2] = k
        }
        if (icont[3] > icont[4]) {
            k = icont[3]
            icont[3] = icont[4]
            icont[4] = k
        }

        if ((icont[2] >= imin) || (icont[3] <= imax)) {
            glogprint (l_logfile, "nsmap", "status", type="error", errno=121,
                str="The continuum regions overlap the line region.", verbose+)
            goto clean
        }

        xval2 = wavel0 + (delwl * (real(icont[1])-pixel0)) - delwl/2.
        xval1 = wavel0 + (delwl * (real(icont[2])-pixel0)) + delwl/2.
        imstat (tmpspectrum//"["//str(icont[1]-5)//":"//str(icont[4]+5)//"]",
            format-, fields="min", lower=INDEF, upper=INDEF) | \
            scanf ("%f",ymin)
        imstat (tmpspectrum//"["//str(icont[1]-5)//":"//str(icont[4]+5)//"]",
            format-, fields="max", lower=INDEF, upper=INDEF) | \
            scanf ("%f",ymax)
        tmpcursor = mktemp ("tmpcursor")
        yval1 = ymin + (ymax-ymin)/10.
        yval2 = ymin + (ymax-ymin)*9./10.
        printf (" 0 0 1 :xydraw yes\n %f %f 1 x\n %f %f 1 x\n", xval1, yval1,
            xval2, yval1, >> tmpcursor)
        xval2 = wavel0 + (delwl * (real(icont[3])-pixel0)) - delwl/2.
        xval1 = wavel0 + (delwl * (real(icont[4])-pixel0)) + delwl/2.
        printf (" 0 0 1 :xydraw yes\n %f %f 1 x\n %f %f 1 x\n", xval1, yval1,
            xval2, yval1, >> tmpcursor)
        xval2 = wavel0 + (delwl * (real(imin)-pixel0)) - delwl/2.
        xval1 = wavel0 + (delwl * (real(imax)-pixel0)) + delwl/2.
        printf (" 0 0 1 :xydraw yes\n %f %f 1 x\n %f %f 1 x\n", xval1, yval2,
            xval2, yval2, >> tmpcursor)
        printf (" 0 0 1 q", >> tmpcursor)
        bplot (tmpspectrum//"["//str(icont[1]-5)//":"//str(icont[4]+5)//"]",
            cursor=tmpcursor)
        imdelete (tmpspectrum, verify-, >& "dev$null")
        delete (tmpcursor, verify-, >& "dev$null")

        tmplist = " "

        for (k=icont[1]; k <= icont[2]; k=k+1) {
            tmpwork = mktemp ("tmpwork")
            imcopy (l_inimage//"[1][*,*,"//str(k)//"]", tmpwork, verbose-,
                >& "dev$null")
            if (tmplist == " ") {
                tmplist = tmpwork
            } else {
                tmplist = tmplist//","//tmpwork
            }
        }

        tmpcont1 = mktemp ("tmpcont")
        tmpcont2 = mktemp ("tmpcont")

        imcombine (tmplist, tmpcont1, combine="average", weight="none",
            lthreshold=INDEF, hthreshold=INDEF, reject="none", zero="none",
            scale="none", outtype="real", outlimits="", offsets="none",
            masktype="none", headers="", bpmasks="", rejmasks="", nrejmasks="",
            expmasks="", sigmas="", logfile="STDOUT", project-, blank=0.,
            >& "dev$null")

        tmplist = " "

        for (k=icont[3]; k <= icont[4]; k=k+1) {
            tmpwork = mktemp ("tmpwork")
            imcopy (l_inimage//"[1][*,*,"//str(k)//"]", tmpwork, verbose-,
                >& "dev$null")
            if (tmplist == " ") {
                tmplist = tmpwork
            } else {
                tmplist = tmplist//","//tmpwork
            }
        }

        imcombine (tmplist, tmpcont2, combine="average", weight="none",
            lthreshold=INDEF, hthreshold=INDEF, reject="none", zero="none",
            scale="none", outtype="real", outlimits="", offsets="none",
            masktype="none", headers="", bpmasks="", rejmasks="", nrejmasks="",
            expmasks="", sigmas="", logfile="STDOUT", project-, blank=0.,
            >& "dev$null")

        tmpfinal = mktemp ("tmpcont")

        if (!l_fl_interp) {
            imcombine (tmpcont1//","//tmpcont2, tmpfinal, combine="average",
                weight="none", lthreshold=INDEF, hthreshold=INDEF,
                reject="none", zero="none", scale="none", outtype="real",
                outlimits="", offsets="none", masktype="none", headers="",
                bpmasks="", rejmasks="", nrejmasks="", expmasks="", sigmas="",
                logfile="STDOUT", project-, blank=0., >& "dev$null")
        } else {

            # The following only works when xlow < 0. and xhigh > 0. as is 
            # assumed here

            scale = abs(xlow) / (abs(xlow) + abs(xhigh))
            imarith (tmpcont2, "*", scale, tmpcont2, verbose-, >& "dev$null")
            scale = 1. - scale
            imarith (tmpcont1, "*", scale, tmpcont1, verbose-, >& "dev$null")

            imcombine (tmpcont1//","//tmpcont2, tmpfinal, combine="sum",
                weight="none", lthreshold=INDEF, hthreshold=INDEF,
                reject="none", zero="none", scale="none", outtype="real",
                outlimits="", offsets="none", masktype="none", headers="",
                bpmasks="", rejmasks="", nrejmasks="", expmasks="", sigmas="",
                logfile="STDOUT", project-, blank=0., >& "dev$null")
        }
    } else {
        xval2 = wavel0 + (delwl * (real(imin)-pixel0)) - delwl/2.
        xval1 = wavel0 + (delwl * (real(imax)-pixel0)) + delwl/2.
        imstat (tmpspectrum//"["//str(imin-5)//":"//str(imax+5)//"]",
            format-, fields="min", lower=INDEF, upper=INDEF) | \
            scanf ("%f",ymin)
        imstat (tmpspectrum//"["//str(imin-5)//":"//str(imax+5)//"]",
            format-, fields="max", lower=INDEF, upper=INDEF) | \
            scanf ("%f",ymax)
        tmpcursor = mktemp ("tmpcursor")
        yval1 = ymin + (ymax-ymin)/10.
        yval2 = ymin + (ymax-ymin)*9./10.
        printf (" 0 0 1 :xydraw yes\n %f %f 1 x\n %f %f 1 x\n", xval1, yval2,
            xval2, yval2, >> tmpcursor)
        printf (" 0 0 1 q", >> tmpcursor)
        bplot (tmpspectrum//"["//str(imin-10)//":"//str(imax+10)//"]",
            cursor=tmpcursor)
        imdelete (tmpspectrum, verify-, >& "dev$null")
        delete (tmpcursor, verify-, >& "dev$null")
    }

    tmplist = " "

    tmpstring = "Line extraction region: pixels "//str(imin)//" to "//str(imax)
    glogprint (l_logfile, "nsmap", "status", type="string", str=tmpstring,
        verbose=l_verbose)

    if (l_fl_cont) {
        tmpstring = "Continuum region 1: pixels "//str(icont[1])//" to "//\
            str(icont[2])
        glogprint (l_logfile, "nsmap", "status", type="string", str=tmpstring,
            verbose=l_verbose)
        tmpstring = "Continuum region 2: pixels "//str(icont[3])//" to "//\
            str(icont[4])
        glogprint (l_logfile, "nsmap", "status", type="string", str=tmpstring,
            verbose=l_verbose)
    }

    for (k=imin; k <= imax; k=k+1) {
        tmpwork = mktemp ("tmpwork")
        imcopy (l_inimage//"[1][*,*,"//str(k)//"]", tmpwork, verbose-,
            >& "dev$null")
        if (tmplist == " ") {
            tmplist = tmpwork
        } else {
            tmplist = tmplist//","//tmpwork
        }
    }

    tmpcont1 = mktemp ("tmpline")

    imcombine (tmplist, tmpcont1, combine="average", weight="none",
        lthreshold=INDEF, hthreshold=INDEF, reject="none", zero="none",
        scale="none", outtype="real", outlimits="", offsets="none",
        masktype="none", headers="", bpmasks="", rejmasks="", nrejmasks="",
        expmasks="", sigmas="", logfile="STDOUT", project-, blank=0.,
        >& "dev$null")

    if (l_fl_cont) {
        imarith (tmpcont1, "-", tmpfinal, tmpfinal, >& "dev$null")
        wmef (tmpfinal, l_outimage, extname="SCI", phu=l_inimage//"[0]", \
            verbose-, >& "dev$null")
    } else {
        wmef (tmpcont1, l_outimage, extname="SCI", phu=l_inimage//"[0]", \
            verbose-, >& "dev$null")
    }


    # Time stamps
    gemdate ()
    gemhedit (l_outimage//"[0]", "NFMAP", gemdate.outdate,
        "UT Time stamp for NFMAP", delete-)
    gemhedit (l_outimage//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)

    # clean-up of temp files

    delete ("tmpcont*", verify-, >& "dev$null")
    delete ("tmpwork*", verify-, >& "dev$null")
    delete ("tmpline*", verify-, >& "dev$null")

clean:
    delete (tmpfile, verify-, >& "dev$null")
    delete (tmpinimg, verify-, >& "dev$null")
    delete (tmpoutimg, verify-, >& "dev$null")

    if (status == 0)
        glogclose (l_logfile, "nfmap", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "nfmap", fl_success-, verbose=l_verbose)

exitnow:
    ;

end
