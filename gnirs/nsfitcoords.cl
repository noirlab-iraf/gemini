# Copyright(c) 2009-2013 Association of Universities for Research in Astronomy, Inc.

# NSFITCOORDS -- This task takes arc and s-distortion calibration images,
# which have previously been traced, and fits dispersion and s-distortion
# surface functions.  The fits are stored in the database directory in a
# file for each input 2D spectrum.  Note that the fitting is redundant
# because the fits for many spectra will use the same calibration data
# but the overhead is low and allows unique surface functions for each
# image.  The resulting functions are recorded in the image headers for
# later use; e.g. by NSTRANSFORM and NIFCUBE.

procedure nsfitcoords (inimages) 

# Determine 2D wavelength and S-distortions functions
#
# Version:
#          Dec 01, 2007  NZ  Use nstransform, remove the 'transform' calls.
#          Jan 31, 2008  NZ  Merge step. Add ckcal (old nfckcal) step
#                            that replaces the initial verification steps
#                            of the old nstransform.
#          Mar 02, 2009  EH  Use ckinput instead of nfcheck.

char    inimages    {prompt = "Input GNIRS/NIRI spectra"}               # OLDP-1-primary-single-prefix=t
char    outspectra  {"", prompt = "Output spectra"}                     # OLDP-1-output
char    outprefix   {"f", prompt = "Prefix for output spectra"}         # OLDP-4
char    lamptransf  {"", prompt = "Name of arc lamp spectrum used to determine wavelength calibration"} # OLDP-1-input
char    sdisttransf {"", prompt = "Name of image used to determine S-distortion calibration"}       # OLDP-1-input
int     dispaxis    {1, min = 1, max = 2, prompt = "Dispersion axis if not defined in the header"}  # OLDP-3
char    database    {"", prompt = "Directory for calibration files"}    # OLDP-3
bool    fl_inter    {no, prompt = "Examine fitted transformations interactively"}                   # OLDP-4
bool    fl_align    {no, prompt = "Align images in the spatial direction?"}                         # OLDP-3
char    function    {"chebyshev", min = "legendre|chebyshev|spline1|spline3", prompt = "Coordinate fitting function"}   # OLDP-2
int     lxorder     {2, min = 1, prompt = "X order of lamp fitting function"}                 # OLDP-2
int     lyorder     {2, min = 1, prompt = "Y order of lamp fitting function"}                 # OLDP-2
int     sxorder     {2, min = 1, prompt = "X order of S-distortion fitting function"}                   # OLDP-2
int     syorder     {2, min = 1, prompt = "Y order of S-distortion fitting function"}                   # OLDP-2
real    pixscale    {1., prompt = "Pixel scale in arcsec/pixel if not defined in the header"}       # OLDP-3
char    logfile     {"", prompt = "Logfile"}                            # OLDP-1
bool    verbose     {yes, prompt = "Verbose"}                           # OLDP-2
bool    debug       {no, prompt = "Very verbose"}                       # OLDP-2
bool    force       {no, prompt = "Force use with earlier IRAF versions?"}                          # OLDP-3
int     status      {0, prompt = "Exit status (0=Good)"}                # OLDP-4
struct *scanfile1   {"", prompt = "Internal use"}                       # OLDP-4
struct *scanfile2   {"", prompt = "Internal use"}                       # OLDP-4


begin

    char    l_inimages = ""
    char    l_outspectra = ""
    char    l_outprefix = ""
    char    l_lamptransf = ""
    char    l_sdisttransf = ""
    int     l_dispaxis
    char    l_database = ""
    bool    l_fl_inter
    bool    l_fl_align
    char    l_function = ""
    int     l_lxorder
    int     l_lyorder
    int     l_sxorder
    int     l_syorder
    real    l_pixscale
    char    l_logfile = ""
    bool    l_verbose
    bool    l_debug
    bool    l_force

    char    l_key_dispaxis = ""
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_xoff = ""
    char    l_key_yoff = ""
    char    l_key_pixscale = ""
    char    l_key_mode = ""

    char    badhdr, modetxt, tmpinlist, tmpoutlist, l_fitname, l_fitname_lamp
    char    l_fitname_sdist, lampsec, sdistsec, lampsec2, sdistsec2, line, mdf
    char    img, imgin, imgout, lamp, sdist, sec, sec2, sec4, phu, key, refimg
    char    extn[3], tmplog, instrume, tmpfit, dbfitfile
    int     nlamps, nsdist, nver, version, xo, yo, nfiles, lcount, junk, iextn
    int     nx, ny, nxx, nyy, nspatial, ndisp, count, nfail
    real    dx, dy, x1, x2, y1, y2, pixsc, shift, shiftx, shifty, small
    real    wmin, wmax, smin, smax, margin, coeff, width, correcn
    bool    l_fl_lamptransf, l_fl_sdisttransf, intdbg, ifuproc, prevproc
    bool    beforesurface, nifs 
    struct  l_date, sline

    # Define temporary files
    tmpinlist = mktemp ("tmpinl") 
    tmpoutlist = mktemp ("tmpout") 
    tmplog = mktemp ("tmplog")

    # Initialize exit status
    status = 1
    intdbg = no
    ifuproc = no

    small = 0.001

    # margin is used in the IFU calculations below.
    # if a value of 0 is used then the code will generate a cube that
    # contains the exact aperture (plus padding to get an integer
    # number of pixels at the correct pixel scale).
    # but there's no way you can get data in the edge pixel on either
    # side in this case, unless you have an aperture that lines up
    # *exactly* with the detector. Having a margin of -0.9 discards
    # those "guaranteed empty" pixels.
    margin = -0.9

    # cache parameter files
    cache ("niri", "gemextn", "gemoffsetlist", "gemwcscopy", \
        "nsmdfhelper", "gemdate") 

    # Set the local variables

    junk = fscan (inimages, l_inimages)
    junk = fscan (outspectra, l_outspectra)
    junk = fscan (outprefix, l_outprefix)
    junk = fscan (lamptransf, l_lamptransf)
    junk = fscan (sdisttransf, l_sdisttransf)
    l_dispaxis = dispaxis
    junk = fscan (database, l_database)
    l_fl_inter = fl_inter
    l_fl_align = fl_align
    junk = fscan (function, l_function)
    l_lxorder = lxorder
    l_lyorder = lyorder
    l_sxorder = sxorder
    l_syorder = syorder
    l_pixscale = pixscale
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose
    l_debug = debug
    l_force = force

    badhdr = ""
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_xoff, l_key_xoff)
    if ("" == l_key_xoff) badhdr = badhdr + " key_xoff"
    junk = fscan (nsheaders.key_yoff, l_key_yoff)
    if ("" == l_key_yoff) badhdr = badhdr + " key_yoff"
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
            printlog ("WARNING - NSFITCOORDS: Both nsfitcoords.logfile \
                and gnirs.logfile are", l_logfile, verbose+) 
            printlog ("                        undefined. Using " \
                // l_logfile, l_logfile, verbose+) 
        }
    }

    printlog ("----------------------------------------------------------\
        ---------------------", l_logfile, l_verbose) 
    date | scan (l_date) 
    printlog ("NSFITCOORDS -- " // l_date, l_logfile, l_verbose) 
    printlog ("", l_logfile, l_verbose) 

    if (l_database == "") {
        junk = fscan (gnirs.database, l_database)
        if (l_database == "") {
            l_database = "database"
            printlog ("WARNING - NSFITCOORDS: Both nsfitcoords.database \
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
        printlog ("ERROR - NSFITCOORDS: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }

    # Expand and verify input
    if (intdbg) print ("expansion and verification of input")

    # It is possible that some GNIRS XD data will have already been run through
    # NSFITCOORDS (in this case it is ok to run NSFITCOORDS twice), so don't
    # set procflag="NSFITCOO"
    ckinput (input=l_inimages, output=l_outspectra, prefix=l_outprefix,
        outlist=tmpoutlist, name="NSFITCOORDS", dependflag="*PREPAR*",
        procflag="", sci_ext=l_sci_ext, logfile=l_logfile, verbose=l_verbose,
        vverbose=intdbg)
    if (ckinput.status == 1)
        goto clean

    if (l_fl_align) {
        # tmpoutlist (created by ckinput) contains two columns: the input file
        # name and the output file name. Only the input file name is required
        # for use in gemoffsetlist below.
        fields (tmpoutlist, "1", lines="1-", quit_if_missing=no,
            print_file_names=no > tmpinlist)
        head (tmpinlist, nlines=1) | scan (refimg)
    }

    nfiles = ckinput.nver
    nver= nfiles

    # Check number of lamps and sdist calibration files.
    ckcal (lampfiles=l_lamptransf, sdistfiles=l_sdisttransf,
        outlist=tmpoutlist, nfiles=nfiles, nver=nver, sci_ext=l_sci_ext,
        database=l_database, name="NSFITCOORDS", logfile=l_logfile,
        verbose=l_verbose, vverbose=intdbg)
    if (ckcal.status == 1)
        goto clean

    nlamps = ckcal.nlamps
    nsdist = ckcal.nsdist

    if ((nsdist == 0) && l_fl_align) {
        printlog ("WARNING - NSFITCOORD: no spatial transform, so \
             setting fl_align=no", l_logfile, verbose+)
        l_fl_align = no
    }

    if (l_fl_align) {
        gemoffsetlist ("@"//tmpinlist, refimg, distance=0.0, age=INDEF, \
            targetlist = "dev$null", offsetlist = "dev$null", \
            fl_nearer+, direction = 3, fl_younger+, fl_noref-, \
            wcs_source="phu", key_xoff=l_key_xoff, key_yoff=l_key_yoff, \
            key_date="", key_time="", logfile=l_logfile, verbose=l_debug, \
            force = l_force)

        if (gemoffsetlist.status != 0) {
            printlog ("ERROR - NSFITCOORDS: Error in GEMOFFSETLIST. \
                Can't determine offsets.", l_logfile, verbose+) 
            goto clean
        }
    }

    wmin = INDEF
    wmax = INDEF
    smin = INDEF
    smax = INDEF

    # Main loop for transformation

    lcount = 0
    scanfile1 = tmpoutlist
    l_fl_sdisttransf = (nsdist > 0)
    l_fl_lamptransf = (nlamps > 0)
    while (fscan (scanfile1, imgin, imgout, lamp, sdist) != EOF) {

        if (intdbg) print ("processing " // imgin //" "// imgout)

        # If no output image requested then need to fix scan.
        # Otherwise we need to make a copy of the input.
        if (l_outspectra == "" && l_outprefix == "") {
            sdist = lamp
            lamp = imgout
            imgout = imgin
        }

        phu = imgin // "[0]"

        keypar (phu, "NSFITCOO", silent+)
        if (keypar.found) {
            printlog ("WARNING - NSFITCOORDS: Image " // imgin // " has been \
                been run though NSFITCOORDS before", l_logfile, verbose+) 
        }

        l_dispaxis = dispaxis
        hselect (phu, l_key_dispaxis, yes) | scan (l_dispaxis)

        hselect (phu, "INSTRUME", yes) | scan(instrume)
        nifs = (instrume == "NIFS")

        nfail = 0

        pixsc = l_pixscale
        hselect (phu, l_key_pixscale, yes) | scan (pixsc)

        modetxt = "LS"
        prevproc = ifuproc
        hselect (phu, l_key_mode, yes) | scan (modetxt)
        ifuproc = ("IFU" == strupr (modetxt)) && (no == l_fl_sdisttransf)

        if (no == (ifuproc == prevproc) && ifuproc) {
            printlog ("NSFITCOORDS: No trace, so doing sub-pixel \
                alignment for IFU.", l_logfile, verbose+) 
        }

        if (small > pixsc && l_fl_align) {

            printlog ("WARNING - NSFITCOORDS: Zero pixscale, no \
                alignment for " // imgin, l_logfile, verbose+) 

        } else if (l_fl_align) {

            hselect (phu, l_key_xoff, yes) | scan (shiftx)
            if (0 == nscan ()) {
                printlog ("WARNING - NSFITCOORDS: No X offset for " \
                    // phu, l_logfile, verbose+) 
                shiftx = 0.0
            }

            hselect (phu, l_key_yoff, yes) | scan (shifty)
            if (0 == nscan ()) {
                printlog ("WARNING - NSFITCOORDS: No X offset for " \
                    // phu, l_logfile, verbose+) 
                shifty = 0.0
            }

            if ((1 == l_dispaxis && small < abs (shiftx)) ||
                (2 == l_dispaxis && small < abs (shifty))) {

                printlog ("WARNING - NSFITCOORDS: Non-zero shift in \
                    dispersion direction for", l_logfile, verbose+) 
                printlog ("                       " // imgin \
                    // " (ignoring)", l_logfile, verbose+) 
            }

            if (1 == l_dispaxis) shift = shifty / pixsc
            else shift = shiftx / pixsc
            if (intdbg) print ("shifty " // shifty)
            if (intdbg) print ("shiftx " // shiftx)
            if (intdbg) print ("shift " // shift // ", " // pixsc)

            printlog ("NSFITCOORDS: Aligning " // imgin // " by " \
                // shift // " pixels", l_logfile, l_verbose)

        }

        for (version = 1; version <= nver; version = version + 1) {

            if (intdbg) print ("version " // version)

            sec = "_" // l_sci_ext // "_" // version // "_"
            sec2 = "[" // l_sci_ext // "," // version // "]"
            sec4 = "[" // l_sci_ext // "," // version // ",inherit]"

            lampsec = lamp // sec
            sdistsec = sdist // sec
            lampsec2 = lamp // sec2
            sdistsec2 = sdist // sec2

            if (intdbg) print ("checking for files")
            if ((l_fl_lamptransf && \
                (! access (l_database // "/id" // lampsec))) || \
                (l_fl_sdisttransf && \
                (! access (l_database // "/id" // sdistsec)))) {

                printlog ("WARNING - NSFITCOORDS: Incomplete \
                    transform information for " // imgin // sec, \
                    l_logfile, verbose+) 
                nfail = nfail + 1

            } else {

                if (intdbg) print ("using files")

                if (intdbg) print (imgin // sec4)
                if (intdbg) print ("dispaxis: " // l_dispaxis) 

                longslit.dispaxis = l_dispaxis
                dx = INDEF
                dy = INDEF


                # Fit wavelength solution to input image

                l_fitname_lamp = ""
                if (l_fl_lamptransf) {

                    l_fitname_lamp = imgout // sec // "lamp"
                    if (intdbg)
                        print ("fitting wavelength to " \
                            // lampsec2 // ", " // l_fitname_lamp) 

                    delete (l_database // "/" // l_fitname_lamp, \
                        verify-, >& "dev$null")

                    # juggle file names around here
                    delete (l_database // "/id" // lamp, verify-, \
                        >& "dev$null")
                    copy (l_database // "/id" // lampsec, \
                        l_database // "/id" // lamp)

                    if (nifs) {
                        # Check for maximum order allowed.
                        match ("features", l_database//"/id"//lamp, stop-) |
                            fields ("STDIN", 2, lines="1-") |
                            average (opstring="new_sample") |
                            scan (x1, x2, nx)
                        if (l_dispaxis == 1) {
                            xo = max (2, min (l_lxorder, int (x1)-1))
                            yo = max (2, min (l_lyorder, nx)-1)
                        } else {
                            xo = max (2, min (l_lxorder, nx)-1)
                            yo = max (2, min (l_lyorder, nint(x1)-1))
                        }
                    } else {
                        xo = l_lxorder
                        yo = l_lyorder
                    }

                    if (intdbg) print ("fitcoords...")
                    # Work around for WCS problem.
                    junk = 2
                    hselect (lampsec2, "WCSDIM", yes) | scan (junk)
                    if (junk != 2)
                        gemhedit (lampsec2, "WCSDIM", "", "", delete=yes)

                    fitcoords (lampsec2, fitname = l_fitname_lamp, \
                        interactive = l_fl_inter, combine = yes, \
                        database = l_database, deletions = "", \
                        function = l_function, xorder = xo, \
                        yorder = yo, logfiles = tmplog, \
                        plotfile = "")

                    if (junk != 2)
                        gemhedit (lampsec2, "WCSDIM", junk, "", delete=no)

                    concatenate (tmplog, l_logfile, append+)
                    delete (l_database // "/id" // lamp, verify-, \
                        >& "dev$null")

                    # Extract limits.
                    scanfile2 = tmplog
                    while (fscan (scanfile2, sline) != EOF) {
                        if (fscanf (sline, "    (%d, %d) = %g  (%d, %d) = %g",
                                    junk, junk, x1, junk, junk, x2) != 6) {
                            next
                        }
                        if (wmin == INDEF) {
                            wmin = min (x1, x2)
                            wmax = max (x1, x2)
                        } else {
                            wmin = min (wmin, x1, x2)
                            wmax = max (wmax, x1, x2)
                        }
                    }
                    scanfile2 = ""
                    delete (tmplog, verify-)
                    if (intdbg) print ("...done")
                }


                # Fit the s-distortion solution to input image

                l_fitname_sdist = ""
                if (l_fl_sdisttransf) {

                    l_fitname_sdist = imgout // sec // "sdist"
                    if (intdbg)
                        print ("fitting spatial to " \
                            // sdistsec2 // ", " // l_fitname_sdist) 

                    delete (l_database // l_fitname_sdist, \
                        verify-, >& "dev$null")

                    # juggle file names around here
                    delete (l_database // "/id" // sdist, verify-, \
                        >& "dev$null")
                    copy (l_database // "/id" // sdistsec,
                        l_database // "/id" // sdist)

                    if (nifs) {
                        # Check for maximum order allowed.
                        match ("features", l_database//"/id"//sdist, stop-) |
                            fields ("STDIN", 2, lines="1-") |
                            average (opstring="new_sample") |
                            scan (x1, x2, nx)
                        if (l_dispaxis == 1) {
                            xo = max (2, min (l_sxorder, nx))
                            yo = max (2, min (l_syorder, nint(x1)))
                        } else {
                            xo = max (2, min (l_sxorder, nint (x1)))
                            yo = max (2, min (l_syorder, nx))
                        }
                    } else {
                        xo = l_sxorder
                        yo = l_syorder
                    }

                    if (intdbg) print ("fitcoords...")

                    # Work around for WCS problem.
                    junk = 2
                    hselect (sdistsec2, "WCSDIM", yes) | scan (junk)
                    if (junk != 2)
                        gemhedit (sdistsec2, "WCSDIM", "", "", delete=yes)

                    fitcoords (sdistsec2, fitname = l_fitname_sdist, \
                        interactive = l_fl_inter, combine = yes, \
                        database = l_database, deletions = "", \
                        function = l_function, xorder = xo, \
                        yorder = yo, logfiles = tmplog, plotfile = "") 

                    if (junk != 2)
                        gemhedit (sdistsec2, "WCSDIM", junk, "", delete=no)

                    concatenate (tmplog, l_logfile, append+)
                    delete (l_database // "/id" // sdist, verify-, \
                        >& "dev$null")

                    # Extract limits.
                    scanfile2 = tmplog
                    while (fscan (scanfile2, sline) != EOF) {
                        if (fscanf (sline, "    (%d, %d) = %g  (%d, %d) = %g",
                            junk, junk, x1, junk, junk, x2) != 6) {

                            next
                        }
                        if (smin == INDEF) {
                            smin = min (x1, x2)
                            smax = max (x1, x2)
                        } else {
                            smin = min (smin, x1, x2)
                            smax = max (smax, x1, x2)
                        }
                    }
                    scanfile2 = ""
                    delete (tmplog, verify-)
                    if (intdbg) print ("...done")

                    if (small < pixscale && l_fl_align) {

                        tmpfit = mktemp ("tmpfit")
                        dbfitfile = l_database // "/fc" // imgout \
                            // sec // "sdist"
                        #dbfitfile = l_database // "/" // imgout // "sdist"

                        rename (dbfitfile, tmpfit, field="all")

                        scanfile2 = tmpfit
                        beforesurface = yes
                        count = 0
                        while (fscan (scanfile2, sline) != EOF) {
                            if (beforesurface) {
                                beforesurface = \
                                strstr ("surface", sline) == 0
                            } else {
                                count = count + 1
                                if (count == 9) {
                                    coeff = real (sline) + shift
                                    printf ("%f\n", coeff) | scan (sline)
                                }
                            }
                            print (sline, >> dbfitfile)
                        }

                        delete (tmpfit, verify-, >& "dev$null")
                    }
                }

                # If order along spatial direction = 1 then we
                # need to fix the scale (only correcting shape, 
                # not scale)
                # But only on final fit! If we do this for spatial
                # alone then the wavelength cal cannot be combined
                # because the x size is different.

                if (l_fl_lamptransf && l_fl_sdisttransf) {
                    if (1 == l_dispaxis) {
                        if (1 == l_syorder) {
                            dy = 1
                            printlog ("WARNING - NSFITCOORDS: Keeping \
                                spatial (y) scale constant for " \
                                // imgin // sec, l_logfile, l_verbose) 
                        }
                    } else if (2 == l_dispaxis) {
                        if (1 == l_sxorder) {
                            dx = 1
                            printlog ("WARNING - NSFITCOORDS: Keeping \
                                spatial (x) scale constant for " \
                                // imgin // sec, l_logfile, l_verbose) 
                        }
                    }
                }

                if (l_fitname_lamp == "" && l_fitname_sdist != "")
                    l_fitname = l_fitname_sdist
                if (l_fitname_sdist == "" && l_fitname_lamp != "")
                    l_fitname = l_fitname_lamp
                if (l_fitname_sdist != "" && l_fitname_lamp != "")
                    l_fitname = l_fitname_lamp // "," // l_fitname_sdist

                printlog ("NSFITCOORDS: Fitnames for " // imgin \
                    // ": " // l_fitname, l_logfile, l_verbose) 

                # Image size

                img = imgin // "[" // l_sci_ext // "," // version // "]"
                hselect (img, "i_naxis" // l_dispaxis, yes) \
                    | scan (ndisp)
                hselect (img, "i_naxis" // (3 - l_dispaxis), yes) \
                    | scan (nspatial)

                if (version == 1) {
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
                # Add the WCS and FITCOORDS information to the header.
                # This is only done to the science header. The other
                # headers will get this information their science headers.
                # Copy the input to the output if needed.

                img = imgout // "[" // l_sci_ext // "," // version // "]"
                if (!imaccess (img))
                    copy (imgin//".fits", imgout//".fits", verbose-)

                if (intdbg) print ("fitcoords sci")
                gemhedit (img, l_key_dispaxis, l_dispaxis, "", delete=no)

                if (nifs) {
                    nfwcs (img, mdf="!mdf_file", xmdfname="DEC", ymdfname="RA",
                        dymdfname="slitsize_x", pixscale="!pixscale",
                        xoffset="0.", yoffset="0.", keepwcs-)
                }

                gemhedit (img, "FCDB", l_database, "", delete=no)
                if (l_fitname_lamp != "")
                    gemhedit (img, "FCFIT1", l_fitname_lamp, "", delete=no)
                if (l_fitname_sdist != "")
                    gemhedit (img, "FCFIT2", l_fitname_sdist, "", delete=no)
                if (!isindef(x1))
                    gemhedit (img, "FCX1", x1, "", delete=no)
                if (!isindef(y1))
                    gemhedit (img, "FCY1", y1, "", delete=no)
                if (!isindef(x2))
                    gemhedit (img, "FCX2", x2, "", delete=no)
                if (!isindef(y2))
                    gemhedit (img, "FCY2", y2, "", delete=no)
                if (!isindef(dx))
                    gemhedit (img, "FCDX", dx, "", delete=no)
                if (!isindef(dy))
                    gemhedit (img, "FCDY", dy, "", delete=no)
                if (!isindef(nx))
                    gemhedit (img, "FCNX", nx, "", delete=no)
                if (!isindef(ny))
                    gemhedit (img, "FCNY", ny, "", delete=no)

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
                        
                        printlog ("ERROR - NSFITCOORDS: Cannot access \
                            MDF for " // img, l_logfile, verbose+) 
                        goto clean
                    }

                    nsmdfhelper (mdf, version, img, area="spectrum", \
                        logfile=l_logfile, logname="NSFITCOORDS", \
                        verbose-)
                    if (no == (0 == nsmdfhelper.status)) {
                        printlog ("ERROR - NSFITCOORDS: Cannot read \
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

                # Copy across nsappwave info if no wavelength transform

                if (no == l_fl_lamptransf) {
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
                                gemhedit (imgout // sec, key, sline, "",
                                    delete=no)
                            }

                            key = "CRPIX" // l_dispaxis
                            hselect (imgin // sec, key, yes) \
                                | scan (sline)
                            if (fscan (sline, line) > 0) {
                                gemhedit (imgout // sec, key, sline, "",
                                    delete=no)
                            }

                            key = "CRVAL" // l_dispaxis
                            hselect (imgin // sec, key, yes) \
                                | scan (sline)
                            if (fscan (sline, line) > 0) {
                                gemhedit (imgout // sec, key, sline, "",
                                    delete=no)
                            }

                            key = "CD" // l_dispaxis // "_" // l_dispaxis
                            hselect (imgin // sec, key, yes) \
                                | scan (sline)
                            if (fscan (sline, line) > 0) {
                                gemhedit (imgout // sec, key, sline, "",
                                    delete=no)
                            }

                            key = "WAT" // l_dispaxis // "_000"
                            hselect (imgin // sec, key, yes) \
                                | scan (sline)
                            if (fscan (sline, line) > 0) {
                                gemhedit (imgout // sec, key, sline, "",
                                    delete=no)
                            }

                            key = "WAT" // l_dispaxis // "_001"
                            hselect (imgin // sec, key, yes) \
                                | scan (sline)
                            if (fscan (sline, line) > 0) {
                                gemhedit (imgout // sec, key, sline, "",
                                    delete=no)
                            }

                            key = "WAT" // l_dispaxis // "_002"
                            hselect (imgin // sec, key, yes) \
                                | scan (sline)
                            if (fscan (sline, line) > 0) {
                                gemhedit (imgout // sec, key, sline, "",
                                    delete=no)
                            }
                        }
                    }
                }
            }
        }

        # Update the header

        # Set PHU keywords.
        phu = imgout // "[0]"
        if (!imaccess (phu))
            copy (imgin, imgout, verbose-)
        gemhedit (phu, "FCNLAMPS", nlamps, "", delete=no)
        gemhedit (phu, "FCNSDIST", nsdist, "", delete=no)
        if (l_dispaxis == 1) {
            if (!isindef(wmin))
                gemhedit (phu, "FCX1", wmin, "", delete=no)
            if (!isindef(wmax))
                gemhedit (phu, "FCX2", wmax, "", delete=no)
            if (!isindef(smin))
                gemhedit (phu, "FCY1", smin, "", delete=no)
            if (!isindef(smax))
                gemhedit (phu, "FCY2", smax, "", delete=no)
        } else {
            if (!isindef(wmin))
                gemhedit (phu, "FCY1", wmin, "", delete=no)
            if (!isindef(wmax))
                gemhedit (phu, "FCY2", wmax, "", delete=no)
            if (!isindef(smin))
                gemhedit (phu, "FCX1", smin, "", delete=no)
            if (!isindef(smax))
                gemhedit (phu, "FCX2", smax, "", delete=no)
        }

        flpr
        gemdate () 
        gemhedit (imgout // "[0]", "NSFITCOO", gemdate.outdate, \
            "UT Time stamp for NSFITCOORDS", delete=no)
        gemhedit (imgout // "[0]", "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI", delete=no) 

        # If aligned, move reference WCS to shifted PHU

        if (l_fl_align) {
            if (intdbg)
                print ("wcs from " // refimg // " to " // imgout)

            # This used to require a patch to work correctly...
            gemwcscopy (imgout // "[0]", refimg // "[0]", verbose-, \
                logfile=l_logfile)
            # ...so check with gemoffsetlist
            gemoffsetlist (infiles=imgout, reffile=refimg, distance=1, \
                age=INDEF, targetlist="dev$null", offsetlist="dev$null", \
                fl_nearer+, direction=3, fl_younger+, fl_noref+, \
                wcs_source="phu", key_xoff=l_key_xoff, \
                key_yoff=l_key_yoff, key_date="", key_time="",
                logfile=l_logfile, verbose=l_debug, force=l_force)

            if (gemoffsetlist.status != 0 || gemoffsetlist.count != 1) {
                printlog ("WARNING - NSFITCOORDS: The WCS information \
                    in the aligned file has not been", l_logfile, verbose+)
                printlog ("                       corrected. Stack the \
                    results with gemcombine, not ", l_logfile, verbose+)
                printlog ("                       nscombine or nsstack \
                    (see `help nstransform`).", l_logfile, verbose+)
            }
        }

        if (nfail > 0) {
            hselect (imgout // "[0]", "NSCIEXT", yes) | scan (junk)
            gemhedit (imgout // "[0]", "NSCIEXT", (junk - nfail), "",
                delete=no)
        }
        #nz avoid the extra mdf extension
        #if (intdbg) print ("and copy MDF")
        #gemextn (imgin, check="exists,table", process="expand", \
        #    index="1-", extname="MDF", extversion="", ikparams="", \
        #    omit="name,version", replace="", outfile="STDOUT", \
        #    logfile="", glogpars="", verbose-) | scan (sline)
        #if (0 == gemextn.count || no == (0 == gemextn.fail_count)) {
        #    printlog ("WARNING - NSFITCOORDS: No MDF in " \
        #        // imgin, l_logfile, verbose+)
        #} else {
        #    fxinsert (sline, imgout // ".fits[0]", groups="", \
        #        verbose=intdbg)
        #}
    }

    # Completed successfully
    status = 0

clean:
    if (status == 0) {
        printlog ("", l_logfile, verbose-) 
        printlog ("NSFITCOORDS: Exit status good", l_logfile, l_verbose) 
        printlog ("-----------------------------------------------------------\
            --------------------", l_logfile, l_verbose) 
    }
    scanfile1 = ""
    scanfile2 = ""
    delete (tmpinlist, verify-, >& "dev$null") 
    delete (tmpoutlist, verify-, >& "dev$null") 

end
