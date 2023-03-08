# Copyright(c) 2001-2013 Association of Universities for Research in Astronomy, Inc.

procedure nstransform (inimages) 

# Transform GNIRS/NIRI/NIFS  to rectified and wavelength
# calibrated format.  This requires the 2D dispersion and distortion functions
# to be previously fit by NSFITCOORDS.
#
# Version: Sept 20, 2002 IJ,MT,JJ v1.4 release
#          Aug 19, 2003  KL  IRAF2.12 - new parameter, addonly, in hedit
#          Oct 29, 2003  KL  moved from niri to gnirs package
#          Dec 12, 2007  NZ  Separate the previous version into two tasks
#                            nsfitcoords.cl and nstransform.cl
#          Jan 30, 2007  NZ  Incorporate ckinput (former nfcheck) to replace
#                            gemextn calls.

char    inimages    {prompt = "Input GNIRS/NIRI spectra"}               # OLDP-1-primary-single-prefix=t
char    outspectra  {"", prompt = "Output spectra"}                     # OLDP-1-output
char    outprefix   {"t", prompt = "Prefix for output spectra"}         # OLDP-4
int     dispaxis    {1, min = 1, max = 2, prompt = "Dispersion axis if not defined in the header"}  # OLDP-3
char    database    {"", prompt = "Directory for calibration files"}    # OLDP-3
bool    fl_stripe   {no, prompt = "Stripe in spatial direction (for XD wavelength reference)?"}     # OLDP-3
char    interptype  {"poly3", enum = "nearest|linear|poly3|poly5|spline3", prompt = "Interpolation type for transform"} # OLDP-2
bool    xlog        {no, prompt = "Logarithmic x coordinate for transform"}                         # OLDP-2
bool    ylog        {no, prompt = "Logarithmic y coordinate for transform"}                         # OLDP-2
real    pixscale    {1., prompt = "Pixel scale in arcsec/pixel if not defined in the header"}       # OLDP-3
char    logfile     {"", prompt = "Logfile"}                            # OLDP-1
bool    verbose     {yes, prompt = "Verbose"}                           # OLDP-2
bool    debug       {no, prompt = "Very verbose"}                       # OLDP-2
int     status      {0, prompt = "Exit status (0=Good)"}                # OLDP-4
struct* scanfile1   {"", prompt = "Internal use"}                       # OLDP-4
struct* scanfile2   {"", prompt = "Internal use"}                       # OLDP-4

begin

    char    l_inimages = ""
    char    l_outspectra = ""
    char    l_outprefix = ""
    int     l_dispaxis
    char    l_database = ""
    bool    l_fl_stripe
    char    l_interptype = ""
    bool    l_xlog
    bool    l_ylog
    real    l_pixscale
    char    l_logfile = ""
    bool    l_verbose
    bool    l_debug

    char    l_key_dispaxis = ""
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_pixscale = ""
    char    l_key_mode = ""

    char    badhdr, modetxt, tmpoutlist, tmpcoords, tmpdqint, tmpwave
    char    l_fitname, l_fitname_lamp, l_fitname_sdist, line, img, imgin
    char    imgout, sec, phu, key, mdf, extn[3], tmpmedian, outimg, instrume
    int     nver, version, lcount, junk, iextn, nx, ny, nxx, nyy, nspatial
    int     ndisp, ifupass, dqscale
    real    dx, dy, x1, x2, y1, y2, pixsc, wcsxlo, wcsylo, wcsxhi, wcsyhi
    real    margin, width, correcn
    bool    intdbg, first, ifuproc, prevproc, nifs
    struct  l_date, sline

    # Define temporary files
    tmpoutlist = mktemp ("tmpout") 
    tmpcoords = mktemp ("tmpcoords")
    tmpwave = mktemp ("tmpwave")

    # Initialize exit status
    status = 1
    intdbg = no
    ifuproc = no

    dqscale = 1000

    # margin is used in the IFU calculations below.
    # if a value of 0 is used then the code will generate a cube that
    # contains the exact aperture (plus padding to get an integer
    # number of pixels at the correct pixel scale).
    # but there's no way you can get data in the edge pixel on either
    # side in this case, unless you have an aperture that lines up
    # *exactly* with the detector. having a margin of -0.9 discards
    # those "guaranteed empty" pixels.
    margin = -0.9

    # cache parameter files
    cache ("niri", "gemhedit", "gemextn", "nsmdfhelper", "gemdate") 

    # Set the local variables

    junk = fscan (inimages, l_inimages)
    junk = fscan (outspectra, l_outspectra)
    junk = fscan (outprefix, l_outprefix)
    l_dispaxis = dispaxis
    junk = fscan (database, l_database)
    l_fl_stripe = fl_stripe
    junk = fscan (interptype, l_interptype)
    l_xlog = xlog
    l_ylog = ylog
    l_pixscale = pixscale
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose
    l_debug = debug

    badhdr = ""
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"
    junk = fscan (nsheaders.key_mode, l_key_mode)
    if ("" == l_key_mode) badhdr = badhdr + " key_mode"

    extn[1] = l_sci_ext
    extn[2] = l_var_ext
    extn[3] = l_dq_ext

    if (l_logfile == "") {
        junk = fscan (gnirs.logfile, l_logfile)
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSTRANSFORM: Both nstransform.logfile \
                and gnirs.logfile", l_logfile, verbose+) 
            printlog ("                        are empty. \
                Using " // l_logfile, l_logfile, verbose+) 
        }
    }

    printlog ("----------------------------------------------------------\
        ---------------------", l_logfile, l_verbose) 
    date | scan (l_date) 
    printlog ("NSTRANSFORM -- " // l_date, l_logfile, l_verbose) 
    printlog ("", l_logfile, l_verbose) 

    if (l_database == "") {
        junk = fscan (gnirs.database, l_database)
        if (l_database == "") {
            l_database = "database"
            printlog ("WARNING - NSTRANSFORM: Both nstransform.database \
                and gnirs.database are", l_logfile, verbose+)
            printlog ("                       undefined. Using " \
                // l_database, l_logfile, verbose+)
        }
    }

    # Remove any trailing slashes from the database name
    if (strlen (l_database) > 1 && \
        substr (l_database, strlen (l_database), strlen (l_database)) \
        == "/") {
        l_database = substr (l_database, 1, strlen (l_database) - 1)
    }

    if ("" != badhdr) {
        printlog ("ERROR - NSTRANSFORM: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }

    # Expand and verify input
        
    if (intdbg) print ("expansion and verification of input")

    # It is possible that some GNIRS XD data will have already been run through
    # NSTRANSFORM (in this case it is ok to run NSTRANSFORM twice), so don't
    # set procflag="NSTRANSF"
    ckinput (input=l_inimages, output=l_outspectra, prefix=l_outprefix,
        outlist=tmpoutlist, name="NSTRANSFORM", dependflag="NSFITCOO",
        procflag="", sci_ext=l_sci_ext, logfile=l_logfile, verbose=l_verbose,
        vverbose=intdbg) 
    if (ckinput.status == 1)
        goto clean

    nver = ckinput.nver

    # Main loop for transformation

    lcount = 0
    scanfile1 = tmpoutlist
    while (fscan (scanfile1, imgin, imgout) != EOF) {

        if (intdbg) print ("processing " // imgin)
        phu = imgin // "[0]"

        keypar (phu, "NSTRANSF", silent+)
        if (keypar.found) {
            printlog ("WARNING - NSTRANSFORM: Image " // imgin // " has been \
                been run though NSTRANSFORM before", l_logfile, verbose+) 
        }

        l_dispaxis = dispaxis
        hselect (phu, l_key_dispaxis, yes) | scan (l_dispaxis)

        pixsc = l_pixscale
        hselect (phu, l_key_pixscale, yes) | scan (pixsc)

        modetxt = "LS"
        prevproc = ifuproc
        hselect (phu, l_key_mode, yes) | scan (modetxt)
        ifuproc = ("IFU" == strupr (modetxt)) && no

        if (no == (ifuproc == prevproc) && ifuproc) {
            printlog ("NSTRANSFORM: No trace, so doing sub-pixel \
                alignment for IFU.", l_logfile, verbose+) 
        }

        x1 = INDEF
        hselect (phu, "FCX1", yes) | scan (x1)
        x2 = INDEF
        hselect (phu, "FCX2", yes) | scan (x2)
        dx = INDEF
        hselect (phu, "FCDX", yes) | scan (dx)
        nx = INDEF
        hselect (phu, "FCNX", yes) | scan (nx)
        y1 = INDEF
        hselect (phu, "FCY1", yes) | scan (y1)
        y2 = INDEF
        hselect (phu, "FCY2", yes) | scan (y2)
        dy = INDEF
        hselect (phu, "FCDY", yes) | scan (dy)
        ny = INDEF
        hselect (phu, "FCNY", yes) | scan (ny)

        hselect (phu, "INSTRUME", yes) | scan (instrume)
        nifs = (instrume == "NIFS")

        for (version = 1; version <= nver; version = version + 1) {

            if (intdbg) print ("version " // version)

            img = imgin // "[" // l_sci_ext // "," // version // "]"

            first = (1 == version)

            if (intdbg) print ("checking for files")
            if (first) {
                if (intdbg) print ("copying phu")
                imcopy (imgin // "[0]", imgout, verbose-)
                gemhedit (imgout, "FCNLAMPS", "", "", delete+)
                gemhedit (imgout, "FCNSDIST", "", "", delete+)
                gemhedit (imgout, "FCS_ENG", "", "", delete+)
                gemhedit (imgout, "FCX1", "", "", delete+)
                gemhedit (imgout, "FCX2", "", "", delete+)
                gemhedit (imgout, "FCY1", "", "", delete+)
                gemhedit (imgout, "FCY2", "", "", delete+)
            }

            if (intdbg) print ("get fitcoords information")

            # DB name is already checked above. The header kw FCDB
            # is not used. The value needs to be passed using the 
            # database par.

            #l_database = ""
            #hselect (img, "FCDB", yes) | scan (l_database)
            #if (l_database == "") {
            #    printlog ("ERROR - NSTRANSFORM: No wavelength or \
            #        S-distortion transformations defined", l_logfile, \
            #        verbose+)
            #    goto clean
            #}
            l_fitname_lamp = ""
            l_fitname_sdist = ""
            hselect (img, "FCFIT1", yes) | scan (l_fitname_lamp)
            hselect (img, "FCFIT2", yes) | scan (l_fitname_sdist)
            if (l_fitname_lamp == "" && l_fitname_sdist != "")
                l_fitname = l_fitname_sdist
            else if (l_fitname_sdist == "" && l_fitname_lamp != "")
                l_fitname = l_fitname_lamp
            else if (l_fitname_sdist != "" && l_fitname_lamp != "")
                l_fitname = l_fitname_lamp // "," // l_fitname_sdist
            else
                l_fitname = ""

            if (l_fitname == "") {
                printlog ("WARNING - NSTRANSFORM: No wavelength or \
                    s-distortion information for " // img // ".", l_logfile,
                    verbose+)
                printlog ("                       The WCS of the image will \
                    be used in the transform.", l_logfile, verbose+)
            } else {
                printlog ("NSTRANSFORM: Fitnames for " // imgin \
                    // ": " // l_fitname, l_logfile, l_verbose)
            }

            # Image size
            img = imgin // "[" // l_sci_ext // "," // version // "]"
            hselect (img, "i_naxis" // l_dispaxis, yes) \
                | scan (ndisp)
            hselect (img, "i_naxis" // (3 - l_dispaxis), yes) \
                | scan (nspatial)

            if (intdbg) print ("using files")
            if (intdbg) print (imgin)
            if (intdbg) print ("dispaxis: " // l_dispaxis) 

            longslit.dispaxis = l_dispaxis
            dx = INDEF
            dy = INDEF

            if (l_fitname_lamp != "" &&
                access (l_database// "/fc" //l_fitname_lamp) == no) {
                printlog ("ERROR - NSTRANSFORM: Wavelength \
                    transformations not found", l_logfile, verbose+)
                goto clean
            }
            if (l_fitname_sdist != "" &&
                access (l_database// "/fc" //l_fitname_sdist) == no) {
                printlog ("ERROR - NSTRANSFORM: S-distortion \
                    transformations not found", l_logfile, verbose+)
                goto clean
            }

            # Image size

            img = imgin // "[" // l_sci_ext // "," // version // "]"
            hselect (img, "i_naxis" // l_dispaxis, yes) \
                | scan (ndisp)
            hselect (img, "i_naxis" // (3 - l_dispaxis), yes) \
                | scan (nspatial)

            if (!nifs) {
                if (no == ifuproc || first) {
                    if (1 == l_dispaxis) {
                        nx = INDEF
                        nxx = 1
                        ny = nspatial
                        nyy = ny
                        x1 = INDEF
                        x2 = INDEF
                        y1 = 1
                        y2 = ny
                    } else {
                        nx = nspatial
                        nxx = nx
                        ny = INDEF
                        nyy = 1
                        x1 = 1
                        x2 = nx
                        y1 = INDEF
                        y2 = INDEF
                    }
                }
            } 
            if (ifuproc) {

                # Spatial direction

                gemextn (imgin, check="exists,table", \
                    process="expand", \
                    index="", extname="MDF", extver="", ikparams="", \
                    omit="", replace="%\[%.fits[%", outfile="STDOUT", \
                    logfile="dev$null", glogpars="", verbose-) \
                    | scan (mdf)
                if (no == (1 == gemextn.count) \
                    || 0 < gemextn.fail_count) {
                    
                    printlog ("ERROR - NSTRANSFORM: Cannot access \
                        MDF for " // img, l_logfile, verbose+) 
                    goto clean
                }

                nsmdfhelper (mdf, version, img, area="spectrum", \
                    logfile=l_logfile, logname="NSTRANFORM", \
                    verbose-)
                if (no == (0 == nsmdfhelper.status)) {
                    printlog ("ERROR - NSTRANSFORM: Cannot read \
                        MDF.", l_logfile, verbose+) 
                    goto clean
                }

                # bleagh. here we need to spatially align slices.
                # nsmdfhelper tells us which pixels were selected
                # from the original, and also what the full extent
                # of the aperture was. we want to reassemble that
                # full aperture, plus some margin.

                # so if the aperture started at ixlo, but the correct
                # start was xlo, then we we have data starting at ixlo
                # rather than xlo. so we need to shift the origin
                # to compensate. note an extra "1" because iraf
                # pixels are 1-based.

                # we can do the same to the right - we want xhi, 
                # corrected by the same origin shift.

                # but then we end up with a width that isn't an
                # integer number of pixels. so we need an extra 
                # correction if we want to stay with the same
                # pixel scale.

                if (1 == l_dispaxis) {
                    y1 = (nsmdfhelper.ylo - nsmdfhelper.iylo) + \
                        1 - margin
                    y2 = (nsmdfhelper.yhi - nsmdfhelper.iylo) + \
                        1 + margin
                    width = y2 - y1 + 1
                    ny = int (width + 1)
                    nyy = ny
                    correcn = (ny - width) / 2.0
                    y1 = y1 - correcn
                    y2 = y2 + correcn
                } else {
                    x1 = (nsmdfhelper.xlo - nsmdfhelper.ixlo) + \
                        1 - margin
                    x2 = (nsmdfhelper.xhi - nsmdfhelper.ixlo) + \
                        1 + margin
                    width = x2 - x1 + 1
                    nx = int (width + 1)
                    nxx = nx
                    correcn = (nx - width) / 2.0
                    x1 = x1 - correcn
                    x2 = x2 + correcn
                }
            }
            # Loop to here for IFU second pass
            ifupass = 1
            while (ifupass <= 2 || (first && ifuproc)) {
                # Transform the data using the fit information

                # Copy to local header - don't want to use "inherit"
                # because we're messing with WCS and the PHU has 
                # its own WCS that might confuse things.

                img = imgin // "[" // l_sci_ext // "," // version // "]"
                gemhedit (img, l_key_dispaxis, l_dispaxis, "Dispersion axis", \
                    delete-)

                if (intdbg) print ("transform sci") 

                if (ifuproc && first) {
                    outimg = mktemp ("tmpifu")
                } else {
                    if (ifuproc && 1 == version) {
                        imdelete (outimg, verify-, >& "dev$null")
                    }
                    outimg = imgout // "[" // l_sci_ext // "," // \
                        version // ",append]"
                }

                ## Reduce this to avoid mis-matches and errors
                ## By this point in the process all images will be 2D,
                ## although some may have started larger
                #gemhedit (img, "WCSDIM", 2, "", delete+)

                # Work around for WCS problem.
                junk = 2
                hselect (img, "WCSDIM", yes) | scan (junk)
                if (junk != 2)
                    gemhedit (img, "WCSDIM", "", "", delete+)

                #if (intdbg && no == isindef (x1) && no == isindef (y1))
                #    print (x1 // "," // y1)
                #if (intdbg && no == isindef (x2) && no == isindef (y2))
                #    print (x2 // "," // y2)
                #if (intdbg && no == isindef (nx) && no == isindef (ny))
                #    print (nx // "," // ny)
                #if (intdbg) print (nxx // "," // nyy)

                transform (img, outimg, fitnames = l_fitname, \
                    database = l_database, interptype = l_interptype, \
                    xlog = l_xlog, ylog = l_ylog, flux = yes, \
                    logfiles = l_logfile, dx=dx, dy=dy, \
                    x1=x1, x2=x2, y1=y1, y2=y2, nx=nx, ny=ny, blank=0) 
                 
                gemhedit (outimg, "FCDB", "", "", delete+)
                gemhedit (outimg, "FCFIT1", "", "", delete+)
                gemhedit (outimg, "FCFIT2", "", "", delete+)
                gemhedit (outimg, "FCNX", "", "", delete+)
                gemhedit (outimg, "FCNY", "", "", delete+)
                gemhedit (outimg, "FCX1", "", "", delete+)
                gemhedit (outimg, "FCX2", "", "", delete+)
                gemhedit (outimg, "FCY1", "", "", delete+)
                gemhedit (outimg, "FCY2", "", "", delete+)

                if (junk != 2)
                    gemhedit (img, "WCSDIM", junk, "", delete-)

                if (nifs) {
                    # Convert WCS to 3D for use with NIFCUBE.
                    nfwcs (outimg, mdf="!mdf_file",xmdfname="DEC",
                        ymdfname="RA",
                        dymdfname="slitsize_x", pixscale="!pixscale",
                        xoffset="0.", yoffset="0.", keepwcs+)

                }
                # Measure and repeat for IFU
                if (ifuproc && first) {

                    if (1 == l_dispaxis) {
                        print ("1 " // (nspatial/2), >> tmpcoords)
                        print (ndisp // " " // (nspatial/2), >> tmpcoords)
                    } else {
                        print ((nspatial/2) // " 1", >> tmpcoords)
                        print ((nspatial/2) // " " // ndisp, >> tmpcoords)
                    }
                    if (intdbg) type (tmpcoords)
                    wcsctran (tmpcoords, tmpwave, outimg, "logical", \
                        "world", columns="1 2 3 4 5 6 7", units="native", \
                        formats="", min_sigdigits=8, verbose-)
                    if (intdbg) type (tmpwave)
                    scanfile2 = tmpwave
                    junk = fscan (scanfile2, wcsxlo, wcsylo)
                    junk = fscan (scanfile2, wcsxhi, wcsyhi)
                    if (1 == l_dispaxis) {
                        x1 = wcsxlo - 0.1 * (wcsxhi - wcsxlo)
                        x2 = wcsxhi + 0.1 * (wcsxhi - wcsxlo)
                        nx = int (1.2 * real (ndisp))
                        nxx = 1
                    } else {
                        y1 = wcsylo - 0.1 * (wcsyhi - wcsylo)
                        y2 = wcsyhi + 0.1 * (wcsyhi - wcsylo)
                        ny = int (1.2 * real (ndisp))
                        nyy = 1
                    }

                    first = no
                } else {
                    ifupass = 2
                }
                ifupass = ifupass + 1
            }

            if (ifuproc) {
                printlog ("NSTRANSFORM: IFU resampled to [" // x1 \
                    // "," // y1 // ":" // x2 // "," // y2 // "]", \
                    l_logfile, l_verbose)
            }

            if (fl_stripe) {
                tmpmedian = mktemp ("tmpmedian") 
                imcopy (imgout // "[" // l_sci_ext // \
                    "," // version // "]", tmpmedian, >& "dev$null")
                if (intdbg) print ("striping " // nxx // "," // nyy)
                median (tmpmedian, tmpmedian, \
                    nxx, nyy, boundary = "wrap", zlo = 0.001, \
                    zhi = INDEF, constant = 0, >& "dev$null")
                imcopy (tmpmedian, imgout // "[" \
                    // l_sci_ext // "," // version // ",overwrite]", \
                    >& "dev$null")
                imdelete (tmpmedian, verify-, >& "dev$null")
            }

            img = imgin // "[" // l_var_ext // "," // version // "]"
            if (imaccess (img)) {
                if (intdbg) print ("transform var") 
                gemhedit (img, l_key_dispaxis, l_dispaxis, "Dispersion axis", \
                    delete-)
                gemhedit (img, "WCSDIM", 2, "", delete+)
                transform (img, imgout // "[" // l_var_ext // "," \
                    // version // ",append]", fitnames = l_fitname, \
                    database = l_database, interptype = l_interptype, \
                    xlog = l_xlog, ylog = l_ylog, flux = yes, \
                    logfiles = l_logfile, dx=dx, dy=dy, \
                    x1=x1, x2=x2, y1=y1, y2=y2, nx=nx, ny=ny,
                    blank=0) 
            }

            img = imgin // "[" // l_dq_ext // "," // version // "]"
            if (imaccess (img)) {
                if (intdbg) print ("transform dq") 
                tmpdqint = mktemp ("tmpdqintimg") 

                # need to scale dq because otherwise interpolation
                # is in integers (or convert to real)
                # and back again!
                imarith (img, "*", dqscale, tmpdqint, title="", \
                    divzero=0, hparams="", pixtype="", \
                    calctype="", verbose-, noact-)

                # transform (must be linear to avoid bleeding)
                gemhedit (tmpdqint, l_key_dispaxis, l_dispaxis, \
                    "Dispersion axis", delete-)
                gemhedit (tmpdqint, "WCSDIM", 2, "", delete+)
                transform (tmpdqint, tmpdqint, fitnames = l_fitname, \
                    database = l_database, interptype = "linear", \
                    xlog = l_xlog, ylog = l_ylog, flux = yes, \
                    logfiles = l_logfile, dx=dx, dy=dy, \
                    x1=x1, x2=x2, y1=y1, y2=y2, nx=nx, ny=ny,
                    blank=1) 

                # anything non-zero has been contaminated
                # so set to 1 (all we can do, really)
                imexpr ("a > 0 ? 1 : 0", \
                    imgout // "[" // l_dq_ext // "," \
                    // version // ",append]", tmpdqint, \
                    dims = "auto", intype = "auto", outtype = "ushort", \
                    refim = "auto", bwidth = 0, btype = "nearest", \
                    bpixval = 0, rangecheck+, verbose-, \
                    exprdb = "none")

                imdelete (tmpdqint, verify-) 
            }

            # Copy across nsappwave info if no wavelength transform

            if (l_fitname_lamp == "") {
                for (iextn = 1; iextn <= 3; iextn = iextn+1) {

                    sec = "[" // extn[iextn] // "," // version // "]"

                    if (imaccess (imgin // sec) && \
                        imaccess (imgout // sec)) {

                        if (intdbg) print ("nsappwave info for " // sec)
                        if (intdbg) print ("from " // imgin)
                        if (intdbg) print ("to " // imgout)

                        key = "CTYPE" // l_dispaxis
                        hselect (imgin // sec, key, yes) \
                            | scan (sline)
                        if (fscan (sline, line) > 0) {
                            gemhedit (imgout // sec, key, \
                                sline, "", delete-)
                        }

                        key = "CRPIX" // l_dispaxis
                        hselect (imgin // sec, key, yes) \
                            | scan (sline)
                        if (fscan (sline, line) > 0) {
                            gemhedit (imgout // sec, key, \
                                sline, "", delete-)
                        }

                        key = "CRVAL" // l_dispaxis
                        hselect (imgin // sec, key, yes) \
                            | scan (sline)
                        if (fscan (sline, line) > 0) {
                            gemhedit (imgout // sec, key, \
                                sline, "", delete-)
                        }

                        key = "CD" // l_dispaxis // "_" // l_dispaxis
                        hselect (imgin // sec, key, yes) \
                            | scan (sline)
                        if (fscan (sline, line) > 0) {
                            gemhedit (imgout // sec, key, \
                                sline, "", delete-)
                        }

                        key = "WAT" // l_dispaxis // "_000"
                        hselect (imgin // sec, key, yes) \
                            | scan (sline)
                        if (fscan (sline, line) > 0) {
                            gemhedit (imgout // sec, key, \
                                sline, "", delete-)
                        }

                        key = "WAT" // l_dispaxis // "_001"
                        hselect (imgin // sec, key, yes) \
                            | scan (sline)
                        if (fscan (sline, line) > 0) {
                            gemhedit (imgout // sec, key, \
                                sline, "", delete-)
                        }

                        key = "WAT" // l_dispaxis // "_002"
                        hselect (imgin // sec, key, yes) \
                            | scan (sline)
                        if (fscan (sline, line) > 0) {
                            gemhedit (imgout // sec, key, \
                                sline, "", delete-)
                        }
                    }
                }
            }
        }

        # Update the header

        flpr
        gemdate () 
        gemhedit (imgout // "[0]", "NSTRANSF", gemdate.outdate, \
            "UT Time stamp for NSTRANSFORM", delete-)
        gemhedit (imgout // "[0]", "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI", delete-)

        if (intdbg) print ("and copy MDF")
        gemextn (imgin, check="exists,table", process="expand", \
            index="1-", extname="MDF", extversion="", ikparams="", \
            omit="name,version", replace="", outfile="STDOUT", \
            logfile="", glogpars="", verbose-) | scan (sline)
        if (0 == gemextn.count || no == (0 == gemextn.fail_count)) {
            printlog ("WARNING - NSTRANSFORM: No MDF in " \
                // imgin, l_logfile, verbose+)
        } else {
            fxinsert (sline, imgout // ".fits[0]", groups="", \
                verbose=intdbg)
        }
        delete (tmpwave, verify-, >& "dev$null")
    }
        
        
    # Completed successfully
    status = 0

clean:
    if (status == 0) {
        printlog ("", l_logfile, verbose-) 
        printlog ("NSTRANSFORM: Exit status good", l_logfile, l_verbose) 
        printlog ("-----------------------------------------------------------\
            --------------------", l_logfile, l_verbose) 
    }
    scanfile1 = ""
    scanfile2 = ""
    delete (tmpoutlist, verify-, >& "dev$null") 
    delete (tmpcoords, verify-, >& "dev$null") 
    delete (tmpwave, verify-, >& "dev$null") 

end
