# Copyright(c) 2001-2012 Association of Universities for Research in Astronomy, Inc.

procedure nsextract (inimages) 

# Extract GNIRS/NIRI spectra to 1D
# 
# Version  Sept 20, 2002 MT,JJ v1.4 release
#          Aug 19,  2003 KL  IRAF2.12 new/modified parameters
#                            hedit: addonly
#                            imcombine: headers,bpmasks,expmasks,outlimits
#                                rejmask->rejmasks, plfile->nrejmasks
#          Oct 29,  2003 KL  moved from niri to gnirs

char    inimages    {prompt = "Input images"}                       # OLDP-1-primary-single-prefix=x
char    outspectra  {"", prompt = "Output images"}                      # OLDP-1-output
char    outprefix   {"x", prompt = "Prefix for output spectra"}         # OLDP-4
int     dispaxis    {1, min = 1, max = 2, prompt = "Dispersion axis (if not read from header)"} # OLDP-3
char    database    {"", prompt = "Directory for calibration files"}    # OLDP-3
int     line        {500, prompt = "Line to identify apertures"}        # OLDP-2
int     nsum        {10, prompt = "Number of dispersion lines to sum or median"} # OLDP-2
real    ylevel      {0.1, prompt = "Aperture width selection (see apall)"} # OLDP-2
real    upper       {10, prompt = "Upper aperture limit relative to center"} # OLDP-2
real    lower       {-10, prompt = "Lower aperture limit relative to center"} # OLDP-2
char    background  {"average", prompt = "Background to subtract", enum = "none|average|median|minimum|fit"} # OLDP-2
bool    fl_vardq    {yes, prompt = "Propagate VAR/DQ planes?"} 
bool    fl_addvar   {yes, prompt = "Add average variance to data before extraction?"} # OLDP-2
bool    fl_skylines {yes, prompt = "Allow variance to vary in the dispersion direction?"} # OLDP-2
bool    fl_inter    {no, prompt = "Fit interactively?"}                 # OLDP-2
bool    fl_apall    {yes, prompt = "Use apall for extraction?"}         # OLDP-2
bool    fl_trace    {no, prompt = "Retrace apall reference?"}           # OLDP-2
char    aptable     {"gnirs$data/apertures.fits", prompt = "Table of aperture data"} # OLDP-3
bool    fl_usetabap {no, prompt = "Use upper and lower apertures specified in aptable?"}
bool    fl_flipped  {yes, prompt = "Flipped in dispersion direction?"}  # OLDP-3
bool    fl_project  {no, prompt = "Project spectra to 1D if not using apall?"} # OLDP-2
bool    fl_findneg  {no, prompt = "Find and extract a negative spectrum?"} # OLDP-2
char    bgsample    {"*", prompt = "Default sample for background for apall"} # OLDP-3
char    trace       {"", prompt = "Reference trace for apall"}          # OLDP-2
int     tr_nsum     {10, prompt = "Number of dispersion lines to sum in trace"} # OLDP-3
int     tr_step     {10, prompt = "Tracing step"}                       # OLDP-3
int     tr_nlost    {3, prompt = "Number of consecutive times profile can be lost"} # OLDP-3
char    tr_function {"legendre", prompt = "S-distortion fitting function"} # OLDP-3
int     tr_order    {5, prompt = "S-distortion fitting order"}          # OLDP-3
char    tr_sample   {"*", prompt = "Trace sample"}                      # OLDP-3
int     tr_naver    {1, prompt = "Trace average or median"}             # OLDP-3
int     tr_niter    {0, prompt = "Trace rejection iterations"}          # OLDP-3
real    tr_lowrej   {3., prompt = "Trace lower rejection threshold"}    # OLDP-3
real    tr_highrej  {3., prompt = "Trace upper rejection threshold"}    # OLDP-3
real    tr_grow     {0., prompt = "Trace rejection growing radius"}     # OLDP-3
char    weights     {"variance", enum="none|variance", prompt = "Extraction weights (none|variance)"}
char    logfile     {"", prompt = "Logfile name"}                       # OLDP-1
bool    verbose     {yes, prompt = "Verbose output?"}                   # OLDP-4
int     status      {0, prompt = "Exit status (0=good)"}                # OLDP-4
struct  *scanfile1  {"", prompt = "For internal use only"}              # OLDP-4
struct  *scanfile2  {"", prompt = "For internal use only"}              # OLDP-4

begin
        
    char    l_inimages = ""
    char    l_outspectra = ""
    char    l_prefix = ""
    char    l_key_dispaxis = ""
    int     l_dispaxis
    char    l_database = ""
    int     l_line, l_nsum
    real    l_ylevel, l_upper, l_lower
    char    l_background = ""
    bool    l_fl_vardq, l_fl_addvar, l_fl_skylines, l_fl_inter, l_fl_apall
    bool    l_fl_trace
    char    l_aptable = ""
    bool    l_fl_usetabap
    bool    l_fl_flipped
    bool    l_fl_project, l_fl_findneg
    char    l_bgsample = ""
    char    l_trace = ""
    int     l_tr_nsum, l_tr_step, l_tr_nlost
    char    l_tr_function = ""
    int     l_tr_order
    char    l_tr_sample = ""
    int     l_tr_naver, l_tr_niter
    real    l_tr_lowrej, l_tr_highrej
    real    l_tr_grow
    char    l_weights = ""
    char    l_logfile = ""
    bool    l_verbose

    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_arrayid = ""
    char    l_key_prism = ""
    char    l_key_decker = ""
    char    l_key_fpmask = ""

    int     junk, ndisp, apslo, apshi
    bool    debug, first, used, usetable, dotrace, firstext
    char    badhdr, img, phu, sci, secphu, imgsec, tmpimg, secvar, secdq
    int     nin, nout, nbad, version
    char    tmpextn, extn, prev, ver, text, tmpsn
    char    prism, decker, fpmask, aprow, search
    char    arrayid = ""
    int     ap_line, ap_order, ap_upper, ap_lower
    char    ap_sample, tracename, dbtracename, bg_sample
    char    tmpvar, tmpdq, tmpscivar
    real    medn, noise, backup_lower, backup_upper

    char    speca, specb, specvara, specvarb, specdqa, specdqb
    char    imlog, imgin, imgout
    char    tsci, tvar, tdq, tmplog, sdum
    char    scispec, varspec, dqspec, varb, dqb, unused
    char    l_temp, l_logfilesave, l_databasesave
    bool    l_verbosesave, useinputvar, useinputdq, stored
    int     l_dispaxissave, numext
    char    tmpfile1, tmpinlist, tmprootlist, tmpoutlist
    char    sec, imcombtdq, keyfound
    int     varpl, nmax, appix
    int     count, nfiles, oldcount, lcount
    real    ycp, xcp, xcn, ycn, y1, y2, npixr
    real    xmin, xmax, ymin, ymax
    bool    fl_runbefore
    struct  sdate, expr, sline

    debug = no
    stored = no
    status = 1
    unused = "(unused)"

    tmpfile1    = mktemp ("tmpfile1") 
    tmpinlist   = mktemp ("tmpinl") 
    tmprootlist = mktemp ("tmproot") 
    tmpoutlist  = mktemp ("tmpout") 
    tmplog      = mktemp ("tmplog") 
    tdq         = mktemp ("tmpdq") 
    imcombtdq   = mktemp ("tmpimcombdq") 
    tmpextn     = mktemp ("tmpextn")
    tmpsn       = ""

    speca    = mktemp ("tmpspeca")
    specb    = mktemp ("tmpspecb")
    varb     = mktemp ("tmpvarb")
    dqb      = mktemp ("tmpdqb")
    tsci     = mktemp ("tmpsci")
    scispec  = mktemp ("tmpscispec") 
    varspec  = mktemp ("tmpvarspec") 
    dqspec   = mktemp ("tmpdqspec") 

    junk = fscan (  inimages, l_inimages)
    junk = fscan (  outspectra, l_outspectra)
    junk = fscan (  outprefix, l_prefix)
    # dispaxis read later if required
    junk = fscan (  database, l_database)
    l_line =        line
    l_nsum =        nsum
    l_ylevel =      ylevel
    l_lower =       lower
    l_upper =       upper
    junk = fscan (  background, l_background)
    l_fl_vardq =    fl_vardq
    l_fl_addvar =   fl_addvar
    l_fl_skylines = fl_skylines
    l_fl_inter =    fl_inter
    l_fl_apall =    fl_apall
    l_fl_trace =    fl_trace
    junk = fscan (  aptable, l_aptable)
    l_fl_usetabap = fl_usetabap
    l_fl_flipped =  fl_flipped
    l_fl_project =  fl_project
    l_fl_findneg =  fl_findneg
    l_bgsample =    bgsample # Can contain spaces
    junk = fscan (  trace, l_trace)
    l_tr_nsum =     tr_nsum
    l_tr_step =     tr_step
    l_tr_nlost =    tr_nlost
    junk = fscan (  tr_function, l_tr_function)
    l_tr_order =    tr_order
    l_tr_sample =   tr_sample # Can contain spaces
    l_tr_naver =    tr_naver
    l_tr_niter =    tr_niter
    l_tr_lowrej =   tr_lowrej
    l_tr_highrej =  tr_highrej
    l_tr_grow =     tr_grow
    junk = fscan (  weights, l_weights)
    junk = fscan (  logfile, l_logfile)
    l_verbose =     verbose

    badhdr = ""
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_arrayid, l_key_arrayid)
    if ("" == l_key_arrayid) badhdr = badhdr + " key_arrayid"
    junk = fscan (nsheaders.key_prism, l_key_prism)
    if ("" == l_key_prism) badhdr = badhdr + " key_prism"
    junk = fscan (nsheaders.key_decker, l_key_decker)
    if ("" == l_key_decker) badhdr = badhdr + " key_decker"
    junk = fscan (nsheaders.key_fpmask, l_key_fpmask)
    if ("" == l_key_fpmask) badhdr = badhdr + " key_fpmask"

    if (l_weights == "variance") {
        if (l_background == "none") varpl = 3
        else varpl = 4
    }

    # Keep some parameters from changing by outside world
    cache ("gemhedit", "specred", "gemextn", "keypar", "gemdate") 

    # Check logfile
    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile)
        if ("" == l_logfile) {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSEXTRACT: Both nsextract.logfile and \
                gnirs.logfile are", l_logfile, verbose+) 
            printlog ("                     undefined.  Using " \
                // l_logfile, l_logfile, verbose+)
        }
    }

    date | scan (sdate) 
    printlog ("----------------------------------------------------------\
        ---------------------", l_logfile, l_verbose) 
    printlog ("NSEXTRACT -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NSEXTRACT: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }

    if (l_database == "") {
        junk = fscan (gnirs.database, l_database)
        if (l_database == "") {
            l_database = "database"
            printlog ("WARNING - NSEXTRACT: Both nsextract.database \
                and gnirs.database are", l_logfile, verbose+)
            printlog ("                     undefined.  Using " \
                // l_database, l_logfile, verbose+)
        }
    }

    # Remove any trailing slashes from the database name
    if (strlen (l_database) > 1 && \
        substr (l_database, strlen (l_database), strlen (l_database)) \
        == "/") {
        l_database = substr (l_database, 1, strlen (l_database) - 1)
    }

    # Store default specred values
    l_logfilesave = specred.logfile
    l_databasesave = specred.database
    l_dispaxissave = specred.dispaxis
    l_verbosesave = specred.verbose
    stored = yes

    specred.logfile = l_logfile
    specred.database = l_database
    specred.verbose = l_verbose

    # Check background logic
    if ("none" == l_background && l_fl_addvar && l_fl_apall) {
        printlog ("ERROR - NSEXTRACT: No background subtraction, \
            but variance added to data", l_logfile, verbose+)
        printlog ("                   (change background or fl_addvar).", \
            l_logfile, verbose+)
        goto clean
    }

    # Make lists of input spectra
    # Put all input images in a temporary file: tmpfile1

    gemextn (l_inimages, proc="none", check="mef, exists", \
        index="", extname="", extver="", ikparam="", replace="", \
        omit="extension", outfile=tmpfile1, logfile="", glogpars="",
        verbose-)

    if (0 != gemextn.fail_count || 0 == gemextn.count) {
        printlog ("ERROR - NSEXTRACT: Problem with input images \
            (missing or not MEF).", l_logfile, verbose+) 
        goto clean
    }

    # Filter/expand to tmpinlist if not already processed
    nbad = 0
    nin = 0

    scanfile1 = tmpfile1
    while (fscan (scanfile1, img) != EOF) {

        phu = img // "[0]"

        keyfound = ""
        hselect(phu, "*PREPAR*", yes) | scan(keyfound)
        if (keyfound == "") {

            printlog ("ERROR - NSEXTRACT: Image " // img // " not \
                PREPAREd.", l_logfile, verbose+) 
            nbad += 1

        } else {

            used = no
            delete (tmpextn, verify-, >& "dev$null")
            gemextn (img, proc="expand", check="exists", index="", \
                extname=l_sci_ext, extver="1-", omit="", replace="", \
                ikparams="", outfile=tmpextn, logfile="dev$null", 
                glogpars="", verbose-)
            numext = gemextn.count

            scanfile2 = tmpextn
            firstext = yes
            while (fscan (scanfile2, extn) != EOF) {
                text = substr (extn, stridx (",", extn) + 1, \
                stridx ("]", extn) - 1)
                if (debug) print ("extension: " // text // " (" // extn // ")")
                version = int (text)
                sci = img // "[" // l_sci_ext // "," // version \
                    // ",inherit]"
                keypar (sci, "NSEXTRAC", silent+)
                if (keypar.found) {
                    printlog ("WARNING - NSEXTRACT: Extension " // sci \
                        // " has been run through nsextract before.", \
                        l_logfile, verbose+) 
                } else {
                    used = yes
                    print (img // " " // version // " " // (numext > 1), \
                        >> tmpinlist)
                }
                # Check if variance and dq planes exist, if required
                useinputvar = no
                useinputdq = no
                if (l_fl_vardq) {
                    if (l_weights == "variance") {
                        if (l_fl_apall) {
                            if (firstext)
                                printlog ("WARNING - NSEXTRACT: Output \
                                    variance planes derived from apall rather \
                                    than input variance plane.", l_logfile, \
                                    verbose+)
                        } else {
                            # weights = "variance" is only used in apall
                            if (firstext) {
                                printlog ("WARNING - NSEXTRACT: fl_apall=no, \
                                    so cannot use extraction weighted by ", 
                                    l_logfile, verbose+)
                                printlog ("                     variance, \
                                    setting weights=none.", l_logfile, \
                                    verbose+)
                            }
                            l_weights = "none"
                        }
                    }
                    if (l_weights == "none") {
                        gemextn (img, check="exists", process="expand", 
                            index="", extname=l_var_ext, extversion="1-", 
                            ikparams="", omit="", replace="", 
                            outfile="dev$null", logfile="dev$null", 
                            glogpars="", verbose-)
                        if (0 != gemextn.fail_count || 0 == gemextn.count) {
                            printlog ("ERROR - NSEXTRACT: Variance plane \
                                does not exist in input image. Please set \
                                fl_vardq=no or if fl_apall=yes, please set \
                                weights=variance.", l_logfile, verbose+) 
                            goto clean
                        } else {
                            if (firstext) {
                                printlog ("WARNING - NSEXTRACT: Weights \
                                    parameter set to none.", l_logfile, \
                                    verbose+) 
                                printlog ("                     Creating \
                                    output variance planes from input \
                                    image.", l_logfile, verbose+)
                            }
                            useinputvar = yes
                        }
                    }
                    gemextn (img, check="exists,mef", process="expand", 
                        index="", extname=l_dq_ext, extversion="1-", 
                        ikparams="", omit="", replace="", outfile="dev$null",
                        logfile="dev$null", glogpars="", verbose-)
                    if (0 != gemextn.fail_count || 0 == gemextn.count) {
                        if (firstext) 
                            printlog ("WARNING - NSEXTRACT: DQ planes do not \
                                exist in input image. Setting output DQ \
                                planes to zero.", l_logfile, verbose+)
                    } else 
                        useinputdq = yes
                }
                firstext = no
            }

            if (used) {
                print (img, >> tmprootlist)
                nin = nin + 1
            }

            if (l_fl_apall && l_fl_addvar) {
                gemextn (img, check="exists,mef", process="expand", 
                    index="", extname=l_var_ext, extversion="1-", 
                    ikparams="", omit="", replace="", outfile="dev$null",
                    logfile="dev$null", glogpars="", verbose-)
                if (0 == gemextn.count) {
                    printlog ("WARNING - NSEXTRACT: No " // l_var_ext \
                        // " extension in " // img, l_logfile, verbose+) 
                    printlog ("                     Setting \
                        fl_addvar=no", l_logfile, verbose+) 
                    l_fl_addvar = no
                }
            }
        }
    }

    if (nbad > 0) {
        printlog ("ERROR - NSEXTRACT: " // nbad // " image(s) have not \
            been prepared.", l_logfile, verbose+) 
        goto clean
    }
    if (nin == 0) {
        printlog ("ERROR - NSEXTRACT: No input images to process.", 
            l_logfile, verbose+) 
        goto clean
    }

    printlog ("Using input files:", l_logfile, l_verbose) 
    if (l_verbose) type (tmpinlist) 
    delete (tmpfile1, verify-, >& "dev$null") 

    # Check output images

    if (debug) print ("output check")

    gemextn (l_outspectra, check="absent", process="none", index="", \
        extname="", extversion="", ikparams="", omit="kernel, exten", 
        replace="", outfile=tmpoutlist, logfile="", glogpars="",
        verbose=l_verbose)
    if (gemextn.fail_count != 0) {
        printlog ("ERROR - NSEXTRACT: Existing or incorrectly formatted \
            output files", l_logfile, verbose+) 
        goto clean
    }

    # If tmpoutlist is empty, the output files names should be 
    # created with the prefix parameter
    # (we could have a separate msg for outprefix="", but that will
    # trigger an error in gemextn anyway)
    if (gemextn.count == 0) {

        if (debug) print ("output from substitution")

        gemextn ("%^%" // l_prefix // "%" // "@" // tmprootlist, \
            check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", omit="kernel,exten", \
            replace="", outfile=tmpoutlist, logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - NSEXTRACT: No or incorrectly formatted \
                output files", l_logfile, verbose+) 
            goto clean
        }
    }
    nout = gemextn.count

    printlog ("Using output files:", l_logfile, l_verbose) 
    if (l_verbose) type (tmpoutlist) 
    type (tmpoutlist, >> l_logfile) 
    delete (tmpfile1, verify-, >& "dev$null") 

    # Check number of input output images
    if (nin != nout) {
        printlog ("ERROR - NSEXTRACT: Different number of input and \
            output spectra", l_logfile, verbose+) 
        goto clean
    }

    # Paste input and output files and put it in tmpfile1

    scanfile1 = tmpinlist
    scanfile2 = tmpoutlist
    tmpimg = ""
    while (fscan (scanfile1, imgin, version, usetable) != EOF) {
        if (imgin != tmpimg) {
            tmpimg = imgin
            junk = fscan (scanfile2, imgout)
        }
        print (imgin // " " // version // " " // usetable // " " \
            // imgout, >> tmpfile1)
    }

    # Read header info for tracing, if a separate file is used
    if (l_fl_apall && "" != l_trace) {

        if (debug) print ("separate trace, so read header")
        phu = l_trace // "[0]"

        keypar (phu, l_key_prism, silent+)
        print (unused == l_key_prism)
        if (keypar.found) {
            prism = keypar.value
        } else {
            if (unused != l_key_prism) {
                printlog ("WARNING - NSEXTRACT: No " // l_key_prism \
                    // " in " // phu // ".", l_logfile, verbose+)
            }
            prism = ""
        }
        if (debug) print ("prism: " // prism)

        keypar (phu, l_key_decker, silent+)
        if (keypar.found) {
            decker = keypar.value
        } else {
            if (unused != l_key_decker) {
                printlog ("WARNING - NSEXTRACT: No " // l_key_decker \
                    // " in " // phu // ".", l_logfile, verbose+)
            }
            decker = ""
        }
        if (debug) print ("decker: " // decker)

        keypar (phu, l_key_fpmask, silent+)
        if (keypar.found) {
            fpmask = keypar.value
        } else {
            if (unused != l_key_fpmask) {
                printlog ("WARNING - NSEXTRACT: No " // l_key_fpmask \
                    // " in " // phu // ".", l_logfile, verbose+)
            }
            fpmask = ""
        }
        if (debug) print ("fpmask: " // fpmask)
        
        keypar (phu, l_key_arrayid, silent+)
        if (keypar.found) {
            arrayid = keypar.value
        } else {
            if (unused != l_key_arrayid) {
                printlog ("WARNING - NSEXTRACT: No " // l_key_arrayid \
                    // " in " // phu // ".", l_logfile, verbose+)
            }
            arrayid = ""
        }
        if (debug) print ("arrayid: " // arrayid)
    }

    # Run apall or extract a line

    prev = ""
    scanfile1 = tmpfile1
    while (fscan (scanfile1, imgin, version, usetable, imgout) != EOF) {

        first = prev != imgin
        prev = imgin

        ver = "," // version
        sec = "[" // l_sci_ext // ver // "]"
        secvar = "[" // l_var_ext // ver // "]"
        secdq = "[" // l_dq_ext // ver // "]"
        secphu = "[" // l_sci_ext // ver // ",inherit]"

        keypar (imgin // secphu, l_key_dispaxis, silent+)
        if (keypar.found) l_dispaxis = int (keypar.value)
        else l_dispaxis = dispaxis
        if (debug) print ("dispaxis: " // l_dispaxis)

        specred.dispaxis = l_dispaxis

        if (first) imcopy (imgin // "[0]", imgout, verbose-) 

        # tmp files used within this loop
        speca = mktemp ("tmpspeca")
        specb = mktemp ("tmpspecb")
        if (l_fl_vardq) {
            specvara = mktemp ("tmpspecvara")
            specvarb = mktemp ("tmpspecvarb")
            specdqa  = mktemp ("tmpspecdqa")
            specdqb  = mktemp ("tmpspecdqb")
        }
        varb = mktemp ("tmpvarb")
        dqb = mktemp ("tmpdqb")
        tsci = mktemp ("tmpsci")
        scispec = mktemp ("tmpscispec") 
        varspec = mktemp ("tmpvarspec") 
        dqspec = mktemp ("tmpdqspec") 

        printlog ("", l_logfile, l_verbose) 
 
        ap_line = l_line
        ap_order = l_tr_order
        ap_sample = l_tr_sample
        ap_upper = l_upper
        ap_lower = l_lower
        bg_sample = l_bgsample

        if (l_fl_apall) {

            # If we are tracing using the same file we are extracting
            # from, then we need the header parameters
            if (first && "" == l_trace) {

                if (debug) print ("self-trace, so read header")
                phu = imgin // "[0]"

                keypar (phu, l_key_prism, silent+)
                if (keypar.found) {
                    prism = keypar.value
                } else {
                    if (unused != l_key_prism) {
                        printlog ("WARNING - NSEXTRACT: No " \
                            // l_key_prism // " in " // phu // ".", \
                            l_logfile, verbose+)
                    }
                    prism = ""
                }
                if (debug) print ("prism: " // prism)

                keypar (phu, l_key_decker, silent+)
                if (keypar.found) {
                    decker = keypar.value
                } else {
                    if (unused != l_key_decker) {
                        printlog ("WARNING - NSEXTRACT: No " \
                            // l_key_decker // " in " // phu // ".", \
                            l_logfile, verbose+)
                    }
                    decker = ""
                }
                if (debug) print ("decker: " // decker)

                keypar (phu, l_key_fpmask, silent+)
                if (keypar.found) {
                    fpmask = keypar.value
                } else {
                    if (unused != l_key_fpmask) {
                        printlog ("WARNING - NSEXTRACT: No " \
                            // l_key_fpmask // " in " // phu // ".", \
                            l_logfile, verbose+)
                    }
                    fpmask = ""
                }
                if (debug) print ("fpmask: " // fpmask)
                
                keypar (phu, l_key_arrayid, silent+)
                if (keypar.found) {
                    arrayid = keypar.value
                } else {
                    if (unused != l_key_arrayid) {
                        printlog ("WARNING - NSEXTRACT: No " // l_key_arrayid \
                            // " in " // phu // ".", l_logfile, verbose+)
                    }
                    arrayid = ""
                }
                if (debug) print ("arrayid: " // arrayid)
            }

            # Read params from aptable for multiple versions

            if (usetable) {

                aprow = mktemp ("aprow")//".fits"
                search = "prism == '" // prism // "' && " \
                // "slit == '" // fpmask // "' && " \
                // "decker == '" // decker // "' && " \
                // "arrayid == '" // arrayid // "' && " \
                // "order == " // version
                if (debug) print (search)
                tselect (l_aptable, aprow, search)

                if (debug) tprint (aprow)
                tinfo (table=aprow, >& "dev$null")
                usetable = (1 == tinfo.nrows)
                if (debug) print (tinfo.nrows // ", " // usetable)

                if (no == usetable) {
                    printlog ("WARNING - NSEXTRACT: No aperture \
                        parameters in " // l_aptable, l_logfile, verbose+)
                    printlog ("                     for " // search, \
                        l_logfile, verbose+)
                    printlog ("                     so using default \
                        params.", l_logfile, verbose+)
                } else {
                    if (l_fl_usetabap) {
                        tprint (aprow, prparam-, prdata+, \
                            showrow-, showhdr-, showunits-, \
                            col="apline,aporder,apsample,apupper,aplower,\
                            bgsample", rows=1, pwidth=160) \
                            | scan (ap_line, ap_order, ap_sample, \
                            ap_upper, ap_lower, bg_sample)
                    } else {
                        tprint (aprow, prparam-, prdata+, \
                            showrow-, showhdr-, showunits-, \
                            col="apline,aporder,apsample,bgsample", rows=1, \
                            pwidth=160) | scan (ap_line, ap_order, ap_sample, \
                            bg_sample)
                    }
                }
                delete (aprow, verify-, >& "dev$null")

                if (l_fl_flipped) {
                    keypar (imgin // secphu, "i_naxis" // l_dispaxis, \
                    silent+)
                    if (keypar.found) {
                        if (debug)
                            print ("before flip: "//ap_line//", "//ap_sample)
                        ndisp = int (keypar.value)
                        ap_line = ndisp - ap_line
                        if (2 == \
                            fscanf (ap_sample, "%d:%d", apslo, apshi)) {
                            
                            apslo = ndisp - apslo
                            apshi = ndisp - apshi
                            ap_sample = apshi // ":" // apslo
                        }
                        if (debug)
                            print ("after flip: "//ap_line//", "//ap_sample)
                    } else {
                        printlog ("WARNING - NSEXTRACT: Couldn't flip \
                            apline, apsample.", l_logfile, verbose+)
                    }
                }

            }

            if ("" != l_trace) {
                tracename = l_trace//"["//l_sci_ext//","//version//"]"
            } else {
                tracename = ""
            }
            if (debug) print ("tracename: " // tracename)

            # Check whether the reference has already been traced (if
            # it has, then we don't trace again unless fl_trace is set)

            if ("" != l_trace) {
                dbtracename = \
                    "ap"//l_trace//"_"//l_sci_ext//"_"//version//"_"
            } else {
                dbtracename = \
                    "ap"//imgin//"_"//l_sci_ext//"_"//version//"_"
            }
            if (debug) print ("dbtracename: " // dbtracename)
            dbtracename = l_database // "/" // dbtracename

            dotrace = yes
            if (no == l_fl_trace) {
                dotrace = ! access (dbtracename)
                if (no == dotrace) {
                    printlog ("NSEXTRACT: Using previous trace from " \
                        // dbtracename, l_logfile, l_verbose)
                    printlog ("           (use fl_trace+ to force \
                        retrace)", l_logfile, l_verbose)
                }
            }

            # Need to trace separately if not using same file (don't
            # know why)
            if ("" != l_trace && dotrace) {

                if (debug) print ("calling apall")
                apall (input=tracename, nfind=1, output="", apertures="", \
                    format="multispec", references="", profiles="", \
                    interactive=l_fl_inter, find+, recenter+, resize+, \
                    edit=l_fl_inter, trace=dotrace, fittrace=dotrace && \
                    l_fl_inter, extract-, extras-, review=l_fl_inter, \
                    line=ap_line, nsum=l_nsum, lower=ap_lower, \
                    upper=ap_upper, apidtable="", b_function="chebyshev", \
                    b_order=1, b_sample=bg_sample, b_naverage=-3, \
                    b_niterate=0, b_low_reject=3, b_high_reject=3, b_grow=0, \
                    width=5, radius=10, threshold=0, minsep=5, maxsep=1000, \
                    order="increasing", aprecenter="", npeaks=1, shift+, \
                    llimit=ap_lower, ulimit=ap_upper, ylevel=l_ylevel, peak+, \
                    bkg+, r_grow=0, avglimits-, t_nsum=l_tr_nsum, \
                    t_step=l_tr_step, t_nlost=l_tr_nlost, \
                    t_funct=l_tr_function, t_order=ap_order, \
                    t_sample=ap_sample, t_naver=l_tr_naver, \
                    t_niter=l_tr_niter, t_low_r=l_tr_lowrej, \
                    t_high_r=l_tr_highrej, t_grow=l_tr_grow, \
                    background=l_background, skybox=1, weights=l_weights, \
                    pfit="fit1d", clean-, saturation=INDEF, readnoise=0, \
                    gain=1, lsigma=4, usigma=4, nsubaps=1)
                if (debug) print ("apall done")
                dotrace = no
            }

            # Handle variance propagation

            if (l_fl_addvar) {
                tmpscivar = mktemp ("tmpscivar")
                if (l_fl_skylines) {
                    tmpvar = mktemp ("tmpvar")
                    if (debug) print ("skyline variance: " // tmpvar)
                    improject (imgin // secvar, tmpvar, \
                        projaxis=(3-l_dispaxis), average+, highcut=0.0, \
                        lowcut=0.0, pixtype="real", verbose-)
                    # clunky bug-fix
                    gemhedit (tmpvar, "WAT2_001", "", "", delete+)
                    if (debug) print ("dispaxis: " // l_dispaxis)
                    if (1 == l_dispaxis) {
                        keypar (imgin // secphu, "i_naxis2", silent+)
                        if (debug) print ("naxis: " // keypar.value)
                        imstack (tmpvar, tmpscivar, title="*", pixtype="*")
                        blkrep (tmpscivar, tmpscivar, 1, int (keypar.value))
                    } else {
                        keypar (imgin // secphu, "i_naxis1", silent+)
                        if (debug) print ("naxis: " // keypar.value)
                        imstack (tmpvar, tmpscivar, title="*", pixtype="*")
                        blkrep (tmpscivar, tmpscivar, 1, int (keypar.value))
                        imtranspose (tmpscivar // "[*,-*]", tmpscivar)
                    }
                    imdelete (tmpvar, verify-, >& "dev$null")
                } else {
                    imstat (imgin // secvar, fields="midpt", lower=INDEF, \
                        upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, \
                        binwidth=0.1, format-, cache-) | scan (medn)
                    if (debug) print ("constant variance: " // medn)
                    imcalc (imgin // secvar, tmpscivar, str (medn), \
                        pixtype="old", nullval=0.0, verbose-)
                }
                imarith (imgin // sec, "+", tmpscivar, tmpscivar, \
                    divzero=0.0, hparams="", pixtype="", calctype="", \
                    verbose-, noact-)
            } else {
                tmpscivar = imgin // sec
            }

            # All direct values match defaults except:
            # npeaks
            # recenter
            # Should t_nsum be a separate parameter?
            if (debug) print ("calling apall")
            if (debug) print ("dotrace: " // dotrace)
            apall (input=tmpscivar, nfind=1, output=speca, apertures="", \
                format="multispec", references=tracename, profiles="", \
                interactive=l_fl_inter, find+, recenter-, resize+, \
                edit=l_fl_inter, trace=dotrace, fittrace=dotrace && \
                l_fl_inter, extract+, extras+, review=l_fl_inter, \
                line=ap_line, nsum=l_nsum, lower=ap_lower, upper=ap_upper, \
                apidtable="", b_function="chebyshev", b_order=1, \
                b_sample=bg_sample, b_naverage=-3, b_niterate=0, \
                b_low_reject=3, b_high_reject=3, b_grow=0, width=5, \
                radius=10, threshold=0, minsep=5, maxsep=1000, \
                order="increasing", aprecenter="", npeaks=1, shift=INDEF, \
                llimit=ap_lower, ulimit=ap_upper, ylevel=l_ylevel, peak+, \
                bkg+, r_grow=0, avglimits-, t_nsum=l_tr_nsum, \
                t_step=l_tr_step, t_nlost=l_tr_nlost, t_funct=l_tr_function, \
                t_order=ap_order, t_sample=ap_sample, t_naver=l_tr_naver, \
                t_niter=l_tr_niter, t_low_r=l_tr_lowrej, \
                t_high_r=l_tr_highrej, t_grow=l_tr_grow, \
                background=l_background, skybox=1, weights=l_weights, \
                pfit="fit1d", clean-, saturation=INDEF, readnoise=0, gain=1, \
                lsigma=4, usigma=4, nsubaps=1)
            if (debug) print ("apall done")

            if (useinputvar) {
                # Run apall without weights for variance plane
                # Otherwise they are the same call as above
                if (debug) print ("calling apall for var")
                apall (input=imgin // secvar, nfind=1, output=specvara, \
                    apertures="", format="multispec", references=tmpscivar, \
                    profiles="", interactive-, find-, recenter-, resize-, \
                    edit-, trace-, fittrace-, extract+, extras-, review-, \
                    line=INDEF, nsum=10, lower=-5., upper=5., apidtable="", \
                    b_function="chebyshev", b_order=1, b_sample=bg_sample, \
                    b_naverage=-3, b_niterate=0, b_low_reject=3, \
                    b_high_reject=3, b_grow=0, width=5, radius=10, \
                    threshold=0, minsep=5, maxsep=1000, order="increasing", \
                    aprecenter="", npeaks=1, shift=INDEF, llimit=ap_lower, \
                    ulimit=ap_upper, ylevel=l_ylevel, peak+, bkg+, r_grow=0, \
                    avglimits-, t_nsum=l_tr_nsum, t_step=l_tr_step, \
                    t_nlost=l_tr_nlost, t_funct=l_tr_function, \
                    t_order=ap_order, t_sample=ap_sample, t_naver=l_tr_naver, \
                    t_niter=l_tr_niter, t_low_r=l_tr_lowrej, \
                    t_high_r=l_tr_highrej, t_grow=l_tr_grow, \
                    background="none", skybox=1, weights="none", \
                    pfit="fit1d", clean-, saturation=INDEF, readnoise=0, \
                    gain=1, lsigma=4, usigma=4, nsubaps=1)
                if (debug) print ("apall done for var")
            }
 
            if (useinputdq) {
                # Run apall without weights for DQ plane
                # Otherwise they are the same call as above
                if (debug) print ("calling apall for dq")
                apall (input=imgin // secdq, nfind=1, output=specdqa, \
                    apertures="", format="multispec", references=tmpscivar, \
                    profiles="", interactive-, find-, recenter-, resize-, \
                    edit-, trace-, fittrace-, extract+, extras-, review-, \
                    line=INDEF, nsum=10, lower=-5., upper=5., apidtable="", \
                    b_function="chebyshev", b_order=1, b_sample=bg_sample, \
                    b_naverage=-3, b_niterate=0, b_low_reject=3, \
                    b_high_reject=3, b_grow=0, width=5, radius=10, \
                    threshold=0, minsep=5, maxsep=1000, order="increasing", \
                    aprecenter="", npeaks=1, shift=INDEF, llimit=ap_lower, \
                    ulimit=ap_upper, ylevel=l_ylevel, peak+, bkg+, r_grow=0, \
                    avglimits-, t_nsum=l_tr_nsum, t_step=l_tr_step, \
                    t_nlost=l_tr_nlost, t_funct=l_tr_function, \
                    t_order=ap_order, t_sample=ap_sample, t_naver=l_tr_naver, \
                    t_niter=l_tr_niter, t_low_r=l_tr_lowrej, \
                    t_high_r=l_tr_highrej, t_grow=l_tr_grow, \
                    background="none", skybox=1, weights="none", \
                    pfit="fit1d", clean-, saturation=INDEF, readnoise=0, \
                    gain=1, lsigma=4, usigma=4, nsubaps=1)
                if (debug) print ("apall done for dq")
            }

            if (l_fl_addvar) {
                imdelete (tmpscivar, verify-, >& "dev$null")
                # rename self-trace
                if ("" == l_trace && dotrace) {
                    rename (l_database // "/ap" // tmpscivar, \
                        dbtracename, field="all")
                }
            }

            if (access (dbtracename) ) {
                if (debug) print ("scanning for centre")
                match ("center", dbtracename, stop-) | scan (sdum, xcp, ycp) 
                if (debug) print ("center: " // xcp // ", " // ycp)
                printlog ("NSEXTRACT: Aperture centre at " // xcp \
                    // "," // ycp, l_logfile, l_verbose)
            } else {
                printlog ("ERROR - NSEXTRACT: Cannot find database file " \
                    // dbtracename, l_logfile, verbose+) 
                goto clean
            }

            if (l_fl_findneg) {

                # Find negative spectra

                if (debug) print ("finding negative spectra")

                printlog ("", l_logfile, l_verbose) 
                imdelete ("neg_" // imgin, verify-, >& "dev$null")
                imarith (imgin // sec, "*", -1.0, "neg_" // imgin, verbose-) 

                if (l_fl_addvar) {
                    tmpscivar = mktemp ("tmpscivar")
                    if (l_fl_skylines) {
                        tmpvar = mktemp ("tmpvar")
                        if (debug) print ("skyline variance: " // tmpvar)
                        improject ("neg_" // imgin // secvar, tmpvar, \
                            projaxis=(3-l_dispaxis), average+, highcut=0.0, \
                            lowcut=0.0, pixtype="real", verbose-)
                        # clunky bug-fix
                        gemhedit (tmpvar, "WAT2_001", "", "", delete+)
                        if (debug) print ("dispaxis: " // l_dispaxis)
                        if (1 == l_dispaxis) {
                            keypar (imgin // secphu, "i_naxis2", silent+)
                            if (debug) print ("naxis: " // keypar.value)
                            imstack (tmpvar, tmpscivar, title="*", pixtype="*")
                            blkrep (tmpscivar, tmpscivar, 1, \
                                int (keypar.value))
                        } else {
                            keypar (imgin // secphu, "i_naxis1", silent+)
                            if (debug) print ("naxis: " // keypar.value)
                            imstack (tmpvar, tmpscivar, title="*", pixtype="*")
                            blkrep (tmpscivar, tmpscivar, 1, \
                                int (keypar.value))
                            imtranspose (tmpscivar // "[*,-*]", tmpscivar)
                        }
                        imdelete (tmpvar, verify-, >& "dev$null")
                    } else {
                        imstat (imgin // secvar, fields="midpt", lower=INDEF, \
                            upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, \
                            binwidth=0.1, format-, cache-) | scan (medn)
                        if (debug) print ("constant variance: " // medn)
                        imcalc (imgin // secvar, tmpscivar, str (medn), \
                            pixtype="old", nullval=0.0, verbose-)
                    }
                    imarith ("neg_" // imgin // sec, "+", tmpscivar, \
                        tmpscivar, divzero=0.0, hparams="", pixtype="", \
                        calctype="", verbose-, noact-)
                } else {
                    # Should this be "neg_" // imgin ???
                    tmpscivar = imgin // sec
                }

                if (debug) print ("calling apall")
                apall (input=tmpscivar, nfind=1, output=specb, apertures="", \
                    format="multispec", references=tracename, profiles="", \
                    interactive=l_fl_inter, find+, recenter-, resize+, \
                    edit=l_fl_inter, trace=dotrace, fittrace=dotrace && \
                    l_fl_inter, extract+, extras+, review=l_fl_inter, \
                    line=ap_line, nsum=l_nsum, lower=ap_lower, \
                    upper=ap_upper, apidtable="", b_function="chebyshev", \
                    b_order=1, b_sample=bg_sample, b_naverage=-3, \
                    b_niterate=0, b_low_reject=3, b_high_reject=3, b_grow=0, \
                    width=5, radius=10, threshold=0, minsep=5, maxsep=1000, \
                    order="increasing", aprecenter="", npeaks=1, shift+, \
                    llimit=INDEF, ulimit=INDEF, ylevel=l_ylevel, peak+, bkg+, \
                    r_grow=0, avglimits-, t_nsum=l_tr_nsum, t_step=l_tr_step, \
                    t_nlost=l_tr_nlost, t_funct=l_tr_function, \
                    t_order=ap_order, t_sample=ap_sample, t_naver=l_tr_naver, \
                    t_niter=l_tr_niter, t_low_r=l_tr_lowrej, \
                    t_high_r=l_tr_highrej, t_grow=l_tr_grow, \
                    background=l_background, skybox=1, weights=l_weights, \
                    pfit="fit1d", clean-, saturation=INDEF, readnoise=0, \
                    gain=1, lsigma=4, usigma=4, nsubaps=1)
                if (debug) print ("apall done")
                imdelete ("neg_" // imgin, verify-, >& "dev$null")

                if (useinputvar) {
                    if (debug) print ("calling apall for var for negative \
                        spectrum")
                    apall (input=imgin // secvar, nfind=1, output=specvarb, \
                        apertures="", format="multispec", \
                        references=tmpscivar, profiles="", interactive-, \
                        find-, recenter-, resize-, edit-, trace-, fittrace-, \
                        extract+, extras-, review-, line=INDEF, nsum=10, \
                        lower=-5., upper=5., apidtable="", \
                        b_function="chebyshev", b_order=1, \
                        b_sample=bg_sample, b_naverage=-3, b_niterate=0, \
                        b_low_reject=3, b_high_reject=3, b_grow=0, width=5, \
                        radius=10, threshold=0, minsep=5, maxsep=1000, \
                        order="increasing", aprecenter="", npeaks=1, shift+, \
                        llimit=INDEF, ulimit=INDEF, ylevel=l_ylevel, peak+, \
                        bkg+, r_grow=0, avglimits-, t_nsum=l_tr_nsum, \
                        t_step=l_tr_step, t_nlost=l_tr_nlost, \
                        t_funct=l_tr_function, t_order=ap_order, \
                        t_sample=ap_sample, t_naver=l_tr_naver, \
                        t_niter=l_tr_niter, t_low_r=l_tr_lowrej, \
                        t_high_r=l_tr_highrej, t_grow=l_tr_grow, \
                        background="none", skybox=1, weights="none", \
                        pfit="fit1d", clean-, saturation=INDEF, readnoise=0, \
                        gain=1, lsigma=4, usigma=4, nsubaps=1)
                }

                if (useinputdq) {
                    if (debug) print ("calling apall for dq for negative \
                        spectrum")
                    apall (input=imgin // secdq, nfind=1, output=specdqb, \
                        apertures="", format="multispec", \
                        references=tmpscivar, profiles="", interactive-, \
                        find-, recenter-, resize-, edit-, trace-, fittrace-, \
                        extract+, extras-, review-, line=INDEF, nsum=10, \
                        lower=-5., upper=5., apidtable="", \
                        b_function="chebyshev", b_order=1, \
                        b_sample=bg_sample, b_naverage=-3, b_niterate=0, \
                        b_low_reject=3, b_high_reject=3, b_grow=0, width=5, \
                        radius=10, threshold=0, minsep=5, maxsep=1000, \
                        order="increasing", aprecenter="", npeaks=1, shift+, \
                        llimit=INDEF, ulimit=INDEF, ylevel=l_ylevel, peak+, \
                        bkg+, r_grow=0, avglimits-, t_nsum=l_tr_nsum, \
                        t_step=l_tr_step, t_nlost=l_tr_nlost, \
                        t_funct=l_tr_function, t_order=ap_order, \
                        t_sample=ap_sample, t_naver=l_tr_naver, \
                        t_niter=l_tr_niter, t_low_r=l_tr_lowrej, \
                        t_high_r=l_tr_highrej, t_grow=l_tr_grow, \
                        background="none", skybox=1, weights="none", \
                        pfit="fit1d", clean-, saturation=INDEF, readnoise=0, \
                        gain=1, lsigma=4, usigma=4, nsubaps=1)
                    if (debug) print ("apall done for vardq for negative \
                        spectrum")
                }

                if (l_fl_addvar) {
                    imdelete (tmpscivar, verify-, >& "dev$null")
                    # delete self-trace (what else?)
                    if ("" == l_trace && dotrace) {
                        delete (l_database // "/ap" // tmpscivar, \
                            verify-, >& "dev$null")
                    }
                }

                sarith (speca, "+", specb, scispec, apertures=1, bands=1, \
                    reverse-, format="multispec", verbose-) 

                # Now get the variance plane

                if (useinputvar) {
                    # In this case, weights=none, so no variance plane is 
                    # created by apall. Use the input variance plane to 
                    # create the output variance plane

                    imarith (specvara, "+", specvarb, varspec, verbose-, \
                        noact-)
                } else {
                    # In this case, the output variance plane is taken 
                    # from the output from apall
                    # varpl is a sigma spectrum, need to square

                    imarith (speca // "[*,1," // varpl // "]", "*", 
                        speca // "[*,1," // varpl // "]", speca, verbose-) 
                    imarith (specb // "[*,1," // varpl // "]", "*", 
                        specb // "[*,1," // varpl // "]", specb, verbose-) 
                    imarith (speca, "+", specb, varspec, verbose-, noact-)
                } 

                if (useinputdq) {
                    # Now get the DQ plane. If more than 50% of the pixels 
                    # in the aperture are bad, mark the output pixel in 
                    # the 1D DQ plane as bad

                    # First, get the aperture size in pixels

                    match ("low", l_database //"/ap"// imgin // "_" // \
                        l_dq_ext // "_" // version // "_", stop-) | \
                        fields ("STDIN", fields="2,3", lines="1", \
                        quit_if_miss=yes, print_file_n=no) | scan (xmin, ymin)

                    match ("high", l_database //"/ap"// imgin // "_" // \
                        l_dq_ext // "_" // version // "_", stop-) | \
                        fields ("STDIN", fields="2,3", lines="1", \
                        quit_if_miss=yes, print_file_n=no) | scan (xmax, ymax)

                    # If dispaxis = 1, the spectrum was extracted along rows,
                    # so the extraction aperture in pixels is equal to 
                    # ymax - ymin. If dispaxis = 2, the spectrum was extracted
                    # along columns, so the extraction aperture in pixels is
                    # equal to xmax - xmin
                    if (l_dispaxis == 1)
                        appix = int(0.5 * (ymax - ymin))
                    if (l_dispaxis == 2)
                        appix = int(0.5 * (xmax - xmin))

                    # Then replace pixels

                    imreplace (specdqa, 0, imaginary=0., lower=INDEF, 
                        upper=appix, radius=0.)
                    imreplace (specdqa, 1, imaginary=0., lower=appix, 
                        upper=INDEF, radius=0.)
                    imreplace (specdqb, 0, imaginary=0., lower=INDEF, 
                        upper=appix, radius=0.)
                    imreplace (specdqb, 1, imaginary=0., lower=appix, 
                        upper=INDEF, radius=0.)
                    imarith (specdqa, "+", specdqb, tdq, noact-, verbose-)
                    imreplace (tdq, 1, imaginary=0., lower=1, 
                        upper=INDEF, radius=0.)

                    imarith (tdq, "*", "1", dqspec, noact-, verbose-, \
                        pixtype="short")

                    imdelete (specvara //","// specvarb //","// specdqa \
                        //","// specdqb //","// tdq, verify-, >& "dev$null") 
                } else {
                    # Create an empty DQ plane if no input DQ plane exists
                    imarith (varspec, "*", 0, dqspec, pixtype="short", \
                        verbose-, noact-)
                }
                imdelete (speca //","// specb, verify-, >& "dev$null")

            } else { # if (l_fl_findneg)
                imcopy (speca // "[*,1,1]", scispec, verbose-) 

                if (useinputvar) {
                    # In this case, weights=none, so no variance plane is 
                    # created by apall. Use the input variance plane to 
                    # create the output variance plane

                    imcopy (specvara, varspec, verbose-)
                } else {
                    # In this case, the output variance plane is taken 
                    # from the output from apall
                    # varpl is a sigma spectrum, need to square

                    imarith (speca // "[*,1," // varpl // "]", "*", 
                        speca // "[*,1," // varpl // "]", varspec, \
                        verbose-, noact-) 
                }

                if (useinputdq) {
                    # Now get the DQ plane. If more than 50% of the pixels 
                    # in the aperture are bad, mark the output pixel in 
                    # the 1D DQ plane as bad

                    # First, get the aperture size in pixels

                    match ("low", l_database //"/ap"// imgin // "_" // \
                        l_dq_ext // "_" // version // "_", stop-) | \
                        fields ("STDIN", fields="2,3", lines="1", \
                        quit_if_miss=yes, print_file_n=no) | scan (xmin, ymin)

                    match ("high", l_database //"/ap"// imgin // "_" // \
                        l_dq_ext // "_" // version // "_", stop-) | \
                        fields ("STDIN", fields="2,3", lines="1", \
                        quit_if_miss=yes, print_file_n=no) | scan (xmax, ymax)

                    # If dispaxis = 1, the spectrum was extracted along rows,
                    # so the extraction aperture in pixels is equal to 
                    # ymax - ymin. If dispaxis = 2, the spectrum was extracted
                    # along columns, so the extraction aperture in pixels is
                    # equal to xmax - xmin
                    if (l_dispaxis == 1)
                        appix = int(0.5 * (ymax - ymin))
                    if (l_dispaxis == 2)
                        appix = int(0.5 * (xmax - xmin))

                    # Then replace pixels

                    imreplace (specdqa, 0, imaginary=0., lower=INDEF, 
                        upper=appix, radius=0.)
                    imreplace (specdqa, 1, imaginary=0., lower=appix, 
                        upper=INDEF, radius=0.)

                    imarith (specdqa, "*", "1", dqspec, noact-, verbose-, \
                        pixtype="short")

                    imdelete (specvara //","// specdqa, verify-, >& "dev$null")
                } else {
                    # Create an empty DQ plane if no input DQ plane exists
                    imarith (varspec, "*", 0, dqspec, pixtype="short", \
                        verbose-, noact-)
                }
                imdelete (speca, verify-, >& "dev$null")
            }

            # Remove the keywords associated with the multispec fits file
            gemhedit (scispec//".fits", "CD3_3", "", "", delete+)
            gemhedit (scispec//".fits", "CRPIX3", "", "", delete+)
            gemhedit (scispec//".fits", "CTYPE3", "", "", delete+)
            gemhedit (scispec//".fits", "LTM3_3", "", "", delete+)
            gemhedit (scispec//".fits", "LTV3", "", "", delete+)
            gemhedit (scispec//".fits", "WAT3_001", "", "", delete+)
            gemhedit (scispec//".fits", "WAXMAP01", "", "", delete+)
            # Edit WCSDIM so that it equals 1 instead of 3
            gemhedit (scispec//".fits", "WCSDIM", "1", "", delete-)

            if (l_fl_vardq) {
                # Remove the keywords associated with the multispec fits file
                gemhedit (varspec//".fits", "CD3_3", "", "", delete+)
                gemhedit (varspec//".fits", "CRPIX3", "", "", delete+)
                gemhedit (varspec//".fits", "CTYPE3", "", "", delete+)
                gemhedit (varspec//".fits", "LTM3_3", "", "", delete+)
                gemhedit (varspec//".fits", "LTV3", "", "", delete+)
                gemhedit (varspec//".fits", "WAT3_001", "", "", delete+)
                gemhedit (varspec//".fits", "WAXMAP01", "", "", delete+)
                # Edit WCSDIM so that it equals 1 instead of 3
                gemhedit (varspec//".fits", "WCSDIM", "1", "", delete-)
                gemhedit (dqspec//".fits", "CD3_3", "", "", delete+)
                gemhedit (dqspec//".fits", "CRPIX3", "", "", delete+)
                gemhedit (dqspec//".fits", "CTYPE3", "", "", delete+)
                gemhedit (dqspec//".fits", "LTM3_3", "", "", delete+)
                gemhedit (dqspec//".fits", "LTV3", "", "", delete+)
                gemhedit (dqspec//".fits", "WAT3_001", "", "", delete+)
                gemhedit (dqspec//".fits", "WAXMAP01", "", "", delete+)
                # Edit WCSDIM so that it equals 1 instead of 3
                gemhedit (dqspec//".fits", "WCSDIM", "1", "", delete-)
            }

        } else { # if (l_fl_apall)

            # Find positive spectra

            if (debug) print ("apfind...")

            imlog = l_database // "/ap" // imgin // "_" // l_sci_ext \
                // "_" // version // "_"
            if (access (imlog))
                delete (imlog, verify-, >& "dev$null")

            # Store current apdefault settings
            backup_lower = apdefault.lower
            backup_upper = apdefault.upper

            # Set the default upper and lower aperture values (as defined by
            # apdefault) equal to the upper and lower aperture values as
            # defined by the user, since they are used by apfind.
            apdefault.lower = ap_lower
            apdefault.upper = ap_upper

            specred.apfind (imgin // sec, 1, apertures="", \
                references="", interact=l_fl_inter, find+, \
                recenter-, resize-, edit+, line=l_line, \
                nsum=l_nsum, minsep=5., maxsep=1000., \
                order="increasing") 

            # Restore old apdefault settings
            apdefault.lower = backup_lower
            apdefault.upper = backup_upper

            if (access (imlog) ) {
                # Use the upper and lower aperture values from the database
                # file. If nsextract was run interactively, the database will
                # contain the aperture values interactively selected by the
                # user. Otherwise, the values will be equal to the upper and
                # lower aperture values as defined by the user in the call to
                # nsextract.
                if (debug)
                    print ("scanning for centre")
                match ("center", imlog, stop-) | scan (sdum, xcp, ycp)
                if (debug) 
                    print ("center: " // xcp // ", " // ycp)

                printlog ("NSEXTRACT: Aperture centre at " // xcp \
                    // "," // ycp, l_logfile, l_verbose)

                if (debug)
                    print ("retrieving upper and lower aperture limits from \
                        the database file")
                match ("low", imlog, stop-) | fields ("STDIN", fields="2,3",
                    lines="1", quit_if_miss=yes, print_file_n=no) | \
                    scan (xmin, ymin) 

                match ("high", imlog, stop-) | fields ("STDIN", fields="2,3",
                    lines="1", quit_if_miss=yes, print_file_n=no) | \
                    scan (xmax, ymax)
                if (debug)
                    print ("xmin = "//xmin//"  xmax = "//xmax//"  ymin = "//\
                        ymin//"  ymax = "//ymax)
            } else {
                printlog ("ERROR - NSEXTRACT: Cannot find database file " \
                    // imlog, l_logfile, verbose+) 
                goto clean
            }

            # If dispaxis = 1, the dispersion axis lies along rows, so the
            # extraction aperture in pixels should be calculated using ycp,
            # ymin and ymax. If dispaxis = 2, the dispersion axis lies along
            # columns, so the extraction aperture should be calculated using
            # xcp, xmin and xmax
            if (l_dispaxis == 1) {
                printlog ("NSEXTRACT: Using aperture " // nint(ymin) // ":" \
                    // nint(ymax), l_logfile, l_verbose)
                ap_lower = nint(ycp + ymin)
                ap_upper = nint(ycp + ymax)
                imgsec = "[*," // ap_lower // ":" // ap_upper // "]"
                keypar (imgin // sec, "i_naxis1", silent+)

            } else {
                printlog ("NSEXTRACT: Using aperture " // nint(xmin) // ":" \
                    // xmax, l_logfile, l_verbose)
                ap_lower = nint(xcp + xmin)
                ap_upper = nint(xcp + xmax)
                imgsec = "[" // ap_lower // ":" // ap_upper // ",*]"
                keypar (imgin // sec, "i_naxis2", silent+)
            }
            nmax = int (keypar.value)

            if (ap_lower < 1 || ap_upper > nmax) {
                printlog ("ERROR - NSEXTRACT: Aperture too large \
                    (reduce upper / lower)", l_logfile, verbose+) 
                printlog ("                   Currently " // imgsec, \
                    l_logfile, verbose+) 
                goto clean
            }

            printlog ("NSEXTRACT: Using section "// imgsec //" for \
                extraction", l_logfile, l_verbose)

            if (debug) print ("copy " // imgin // sec // imgsec)
            imcopy (imgin // sec // imgsec, scispec, verbose-) 

            if (useinputvar) {
                if (debug) print ("copy " // imgin // secvar // imgsec)
                imcopy (imgin // secvar // imgsec, varspec, verbose-) 
            }
            if (useinputdq) {
                if (debug) print ("copy " // imgin // secdq // imgsec)
                imcopy (imgin // secdq // imgsec, dqspec, verbose-)
            }

            if (l_fl_findneg) {

                # Find negative spectra
                # TODO - add variance again

                if (debug) print ("finding negative spectra")

                printlog ("", l_logfile, l_verbose) 
                imdelete ("neg_" // imgin, verify-, >& "dev$null")
                imarith (imgin // sec, "*", -1.0, "neg_" // imgin, verbose-) 

                imlog = l_database // "/ap" // "neg_" // imgin
                if (access (imlog))
                    delete (imlog, verify-, >& "dev$null")

                # Store current apdefault settings
                backup_lower = apdefault.lower
                backup_upper = apdefault.upper

                # Set the default upper and lower aperture values (as defined 
                # by apdefault) equal to the upper and lower aperture values as
                # defined by the user, since they are used by apfind.
                apdefault.lower = ap_lower
                apdefault.upper = ap_upper

                specred.apfind ("neg_" // imgin, 1, apertures="", \
                    references="", interact=l_fl_inter, find+, recenter-, \
                    resize-, edit+, line=l_line, nsum=l_nsum, minsep=5., \
                    maxsep=1000., order="increasing") 

                # Restore old apdefault settings
                apdefault.lower = backup_lower
                apdefault.upper = backup_upper

                if (access (imlog) ) {
                    # Use the upper and lower aperture values from the database
                    # file. If nsextract was run interactively, the database
                    # will contain the aperture values interactively selected
                    # by the user. Otherwise, the values will be equal to the
                    # upper and lower aperture values as defined by the user in
                    # the call to nsextract.
                    if (debug)
                        print ("scanning for centre")
                    match ("center", imlog, stop-) | scan (sdum, xcn, ycn)
                    if (debug) 
                        print ("center: " // xcn // ", " // ycn)

                    printlog ("NSEXTRACT: Aperture centre for negative \
                        spectrum at " // xcn // "," // ycn, l_logfile, 
                        l_verbose)
                    imdelete ("neg_" // imgin, verify-, >& "dev$null")

                    if (debug)
                        print ("retrieving upper and lower aperture limits \
                            from the database file")
                    match ("low", imlog, stop-) | fields ("STDIN", 
                        fields="2,3", lines="1", quit_if_miss=yes, 
                        print_file_n=no) | scan (xmin, ymin) 

                    match ("high", imlog, stop-) | fields ("STDIN", 
                        fields="2,3", lines="1", quit_if_miss=yes,
                        print_file_n=no) | scan (xmax, ymax) 
                    if (debug)
                        print ("xmin = "//xmin//"  xmax = "//xmax//\
                            "  ymin = "//ymin//"  ymax = "//ymax)
                } else {
                    printlog ("ERROR - NSEXTRACT: Cannot find database file " \
                        // imlog, l_logfile, verbose+) 
                    goto clean
                }

                # If dispaxis = 1, the dispersion axis lies along rows, so the
                # extraction aperture in pixels should be calculated using ycn,
                # ymin and ymax. If dispaxis = 2, the dispersion axis lies
                # along columns, so the extraction aperture should be
                # calculated using xcn, xmin and xmax
                if (l_dispaxis == 1) {
                    printlog ("NSEXTRACT: Using aperture " // nint(ymin) //\
                        ":" // nint(ymax), l_logfile, l_verbose)
                    ap_lower = nint(ycn + ymin)
                    ap_upper = nint(ycn + ymax)
                    imgsec = "[*," // ap_lower // ":" // ap_upper // "]"
                    keypar (imgin // sec, "i_naxis1", silent+)

                } else {
                    printlog ("NSEXTRACT: Using aperture " // nint(xmin) //\
                        ":" // xmax, l_logfile, l_verbose)
                    ap_lower = nint(xcn + xmin)
                    ap_upper = nint(xcn + xmax)
                    imgsec = "[" // ap_lower // ":" // ap_upper // ",*]"
                    keypar (imgin // sec, "i_naxis2", silent+)
                }
                nmax = int (keypar.value)

                if (ap_lower < 1 || ap_upper > nmax) {
                    printlog ("ERROR - NSEXTRACT: Aperture too large \
                        (reduce upper / lower)", l_logfile, verbose+) 
                    printlog ("                   Currently " // imgsec, \
                        l_logfile, verbose+) 
                    goto clean
                }

                printlog ("NSEXTRACT: Using section "// imgsec //" for \
                    extraction", l_logfile, l_verbose)

                if (debug) print ("copy " // imgin // sec // imgsec)
                imcopy (imgin // sec // imgsec, specb, verbose-) 

                y1 = ycp-ycn
                y2 = mod (y1, int (y1) ) 
                imshift (specb, specb, 0.0, y2, inter="linear") 
                imarith (scispec, "-", specb, scispec, verbose-, noact-) 
                imdelete (specb, verify-, >& "dev$null")

                if (useinputvar) {
                    if (debug) print ("copy " // imgin // secvar // imgsec)
                    imcopy (imgin // secvar // imgsec, varb, verbose-)
                    imarith (varspec, "+", varb, varspec, verbose-, noact-) 
                    imdelete (varb, verify-, >& "dev$null")
                }
                if (useinputdq) {
                    if (debug) print ("copy " // imgin // secdq // imgsec)
                    imcopy (imgin // secdq // imgsec, dqb, verbose-) 
                    tmpdq = mktemp ("tmpdq")
                    imrename (dqspec, tmpdq, verbose-)
                    addmasks(tmpdq//","//dqb, dqspec, "im1 || im2", flags=" ")
                    imdelete (dqb//","//tmpdq, verify-, >& "dev$null")
                }
            }

            # Collapse 2D spectra into 1D

            if (l_fl_project) {

                # Projection occurs across the highest dimension of the input
                # images (i.e., the y axis). If DISPAXIS = 2, rotate the
                # science data (rotate clockwise so that wavelength increases
                # from x = 1 along the 1D spectrum, otherwise the iraf task
                # telluric (which is run from nstelluric) complains).

                if (1 != l_dispaxis) {
                    if (debug) print ("rotating scispec: " // scispec)
                    imtranspose (scispec // "[-*,*]", scispec)
                }

                # Determine the number of rows that will be summed in the call
                # to imcombine below
                hselect (scispec, "i_naxis2", yes) | scan (npixr)

                imrename (scispec, tsci, verbose-) 

                # Add the BPM keyword to the header of the science extension so
                # that bad pixels are rejected / excluded in imcombine
                if (useinputdq) {

                    # If DISPAXIS = 2, rotate the DQ data clockwise
                    if (1 != l_dispaxis) {
                        if (debug) print ("rotating dqspec: " // dqspec)
                        imtranspose (dqspec // "[-*,*]", dqspec)
                    }

                    imdelete (tdq // ".pl", verify-, >& "dev$null")
                    imcopy (dqspec, tdq // ".pl", verbose-) 
                    imdelete (dqspec, verify-, >& "dev$null")
                    gemhedit (tsci, "BPM", tdq // ".pl", "", delete-)
                }

                # The imcombine.project parameter is set to yes, so the highest
                # dimension elements of tsci are combined to produce a
                # 1-dimensional output image. A sum is required, so the
                # imcombine.combine parameter is set to sum.
                imcombine (tsci, scispec, headers="", bpmasks="", \
                    rejmasks="", nrejmasks=imcombtdq, expmasks="", sigma="", \
                    imcmb="", logfile=tmplog, combine="sum", reject="none", \
                    project=yes, outtype="real", outlimits="", \
                    offsets="none", masktype="goodvalue", maskvalue=0, \
                    blank=0, scale="none", zero="none", weight="none", nkeep=0)
                imdelete (tsci, verify-, >& "dev$null")

                type (tmplog, >> l_logfile) 
                delete (tmplog, verify-, >& "dev$null")
                gemhedit (scispec, "BPM", "", "", delete+)

                # Edit the WCSDIM keyword so that it equals 1 instead of 2 
                gemhedit (scispec//".fits", "WCSDIM", "1", "", delete-)
                gemhedit (scispec//".fits", "LTM1_1", "1", "", delete-)
                gemhedit (scispec//".fits", "WAXMAP01", "", "", delete+)

                # Now the science data is 1D, remove the DISPAXIS keyword from
                # the header so nstelluric doesn't complain that the axis
                # parameter exceeds the image dimensions (FYI, apall removes
                # the DISPAXIS keyword) 
                gemhedit (scispec//".fits", "DISPAXIS", "1", "", delete+)

                if (useinputvar) {

                    # If DISPAXIS = 2, rotate the variance data clockwise
                    if (1 != l_dispaxis) {
                        if (debug) print ("rotating varspec: " // varspec)
                        imtranspose (varspec // "[-*,*]", varspec)
                    }

                    # Project the variance data to 1D in the same way as the
                    # science data
                    if (useinputdq) {
                        tvar = mktemp ("tmpvar") 
                        imrename (varspec, tvar, verbose-)
                        gemhedit (tvar, "BPM", tdq // ".pl", "", delete-)
                    }

                    imcombine (tvar, varspec, headers="", bpmasks="", \
                        rejmasks="", nrejmasks="", expmasks="", sigma="", \
                        imcmb="", logfile="", combine="sum", reject="none", \
                        project=yes, outtype="real", outlimits="", \
                        offsets="none", masktype="goodvalue", maskvalue=0, \
                        blank=0, scale="none", zero="none", weight="none", \
                        nkeep=0)
                    imdelete (tvar, verify-, >& "dev$null")

                    gemhedit (varspec, "BPM", "", "", delete+)

                    # Edit the WCSDIM keyword so that it equals 1 instead of 2 
                    gemhedit (varspec//".fits", "WCSDIM", "1", "", delete-)
                    gemhedit (varspec//".fits", "LTM1_1", "1", "", delete-)
                    gemhedit (varspec//".fits", "WAXMAP01", "", "", delete+)

                    # Now the variance data is 1D, remove the DISPAXIS keyword
                    gemhedit (varspec//".fits", "DISPAXIS", "1", "", delete+)
                }

                if (useinputdq) {
                    # If a pixel in imcombtdq (which gives the number of input
                    # pixels rejected or excluded from the input images in the
                    # first call to imcombine) has a value equal or greater
                    # than 50% of the rows that were summed in imcombine
                    # (npixr), then that pixel in the output from imcombine is
                    # bad; mark it as bad in the DQ extension. 
                    expr = "(a >= int(0.5 * " // npixr // ")) ? 1 : 0"
                    imexpr (expr, dqspec, imcombtdq, dims="auto", \
                        intype="auto", outtype="short", refim="auto", \
                        bwidth=0, btype="nearest", bpixval=0., rangecheck+, \
                        verbose+, exprdb="none", >& "dev$null")
                    delete (tdq // ".pl", verify-, >& "dev$null")

                    # Edit the WCSDIM keyword so that it equals 1 instead of 2 
                    gemhedit (dqspec//".fits", "WCSDIM", "1", "", delete-)
                    gemhedit (dqspec//".fits", "LTM1_1", "1", "", delete-)
                    gemhedit (dqspec//".fits", "WAXMAP01", "", "", delete+)

                    # Now the DQ data is 1D, remove the DISPAXIS keyword
                    gemhedit (dqspec//".fits", "DISPAXIS", "1", "", delete+)
                }
                delete (imcombtdq // ".pl", verify-, >& "dev$null")
            }
        }

        # Create the output image
        imcopy (scispec, imgout // "[" // l_sci_ext // ver // ",append]", \
            verbose-)

        if (l_fl_vardq) {
            imcopy (varspec, imgout // "[" // l_var_ext // ver // ",append]", \
                verbose-) 
            imcopy (dqspec, imgout // "[" // l_dq_ext // ver // ",append]", \
                verbose-) 
        }

        imdelete (scispec // "," // varspec // "," // dqspec, verify-, \
            >& "dev$null")

        # Update header of output image (needs much more!)
        gemdate ()
        gemhedit (imgout // "[0]", "NSEXTRAC", gemdate.outdate, \
            "UT Time stamp for NSEXTRACT", delete-)
        gemhedit (imgout // "[0]", "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI IRAF", delete-)

        # Calculate S/N
        # This can only be done if a variance plane is created
        if (l_fl_vardq && (debug || l_verbose)) {
            printlog ("NSEXTRACT: S/N for " // imgout // ": ",
                l_logfile, l_verbose) 
            tmpsn = mktemp ("tmpsn")
            flpr
            imexpr ("a / sqrt(abs(b))", tmpsn, \
                imgout // "[" // l_sci_ext // ver // "]",
                imgout // "[" // l_var_ext // ver // "]", verbose-)
            imstat (tmpsn, fields="midpt", nclip=1, format-) | scan (noise)
            printf ("%6.2f\n", noise) | scan (sline)
            printlog ("           from variance array:    " \
                // sline, l_logfile, l_verbose) 

            imdelete (tmpsn, verify-, >& "dev$null")
            tmpsn = mktemp ("tmpsn")
            median (imgout // "[" // l_sci_ext // ver // "]", \
                tmpsn, xwindow=100, ywindow=1, zlo=INDEF, zhi=INDEF, \
                boundary="nearest", constant=0, verbose-)
            imarith (imgout // "[" // l_sci_ext // ver // "]", \
                "/", tmpsn, tmpsn, title="", divzero=1e6, \
                pixtype="", calctype="", verbose-, noact-)
            imstat (tmpsn, fields="stddev", nclip=1, format-) | scan (noise)
            noise = 1.0 / max (noise, 1e-6)
            printf ("%6.2f\n", noise) | scan (sline)
            printlog ("           from smoothed spectrum: " \
                // sline, l_logfile, l_verbose) 
            imdelete (tmpsn, verify-, >& "dev$null")
        }
    }

    status = 0

clean:
    scanfile1 = ""
    scanfile2 = ""
    delete (tmpinlist // "," // tmpoutlist // "," // tmpfile1,
        verify-, >& "dev$null")
    delete (tmpextn // "," // tmprootlist, verify-, >& "dev$null")
    delete (l_database // "/aplast", verify-, >& "dev$null")
    imdelete (tdq // "," // tsci // "," // specb // "," // tmpsn,
        verify-, >& "dev$null")

    # Restore default specred values
    if (stored) {
        specred.logfile = l_logfilesave
        specred.database = l_databasesave
        specred.dispaxis = l_dispaxissave
        specred.verbose = l_verbosesave
    }

    printlog ("", l_logfile, l_verbose) 
    if (status == 0) 
        printlog ("NSEXTRACT: Exit status good", l_logfile, l_verbose) 

    printlog ("---------------------------------------------------------\
        ----------------------", l_logfile, l_verbose) 
        
end
