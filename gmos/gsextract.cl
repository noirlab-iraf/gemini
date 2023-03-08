# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gsextract (inimages)

# Extract GMOS MOS/Longslit spectra to 1D
#
# Version     Feb 28, 2002  ML,JT  v1.3 release
#             Jun 06, 2002  JT sky subtraction
#             Sept 20, 2002    v1.4 release
#             Sept 30, 2002 IJ aperture finding pix at center, linear interpolation
#             Mar 20, 2003  BM pixel scales for both instruments
#             May 9, 2003   IJ change in logic for instrument, support for old GMOS-N data
#             Oct 28, 2005  BM Fix extraction of VAR plane from apall output

string  inimages    {prompt="Input images"}
string  outimages   {"",prompt="Output images"}
string  outprefix   {"e",prompt="Output prefix"}
string  refimages   {"",prompt="Reference images for tracing apertures"}
real    apwidth     {1.,prompt="Extraction aperture in arcsec (diameter)"}
bool    fl_inter    {no,prompt="Run interactively?"}
string  database    {"database",prompt="Directory for calibration files"}
bool    find        {yes,prompt="Define apertures automatically?"}
bool    recenter    {yes,prompt="Recenter apertures?"}
bool    trace       {yes,prompt="Trace apertures?"}
# int     nfind       {1,min=1,prompt="Number of spectra to find & extract (/slit)"}
int     coloffset   {0, prompt="Adjustment to column for finding source peaks"}
string  tfunction   {"chebyshev",enum="chebyshev|legendre|spline1|spline3",prompt="Trace fitting function"}
int     torder      {5,min=1,prompt="Trace fitting function order"}
int     tnsum       {20,min=1,prompt="Number of dispersion lines to sum for trace"}
int     tstep       {50,min=1,prompt="Tracing step"}
string  weights     {"none",enum="none|variance",prompt="Extraction weights (none|variance)"}

bool    clean       {no,prompt="Detect and replace bad pixels?"}
real    lsigma      {3.0,min=1.0,prompt="Lower rejection threshold for cleaning"}
real    usigma      {3.0,min=1.0,prompt="Upper rejection threshold for cleaning"}

string  background  {"none", enum="none|average|median|minimum|fit", prompt="Background subtraction method"}
string  bfunction   {"chebyshev", enum="chebyshev|legendre|spline1|spline3", prompt="Background function"}
int     border      {1, min=0, prompt="Order of background fit"}
string  long_bsample    {"*",prompt="LONGSLIT: backgr sample regions, WRT aperture"}
real    mos_bsample {0.9, min=0.01, max=1, prompt="MOS: fraction of slit length to use (bkg+obj)"}
int     bnaverage   {1, prompt="Number of samples to average over"}
int     bniterate   {2, min=0, prompt="Number of rejection iterations"}
real    blow_reject {2.5, min=1.0, prompt="Background lower rejection sigma"}
real    bhigh_reject    {2.5, min=1.0, prompt="Background upper rejection sigma"}
real    bgrow       {0.0, min=0.0, prompt="Background rejection growing radius (pix)"}

bool    fl_vardq    {no,prompt="Propagate VAR/DQ planes? (if yes, must use variance weighting)"}
string  sci_ext     {"SCI",prompt="Name of science extension"}
string  var_ext     {"VAR",prompt="Name of variance extension"}
string  dq_ext      {"DQ",prompt="Name of data quality extension"}
string  key_ron     {"RDNOISE",prompt="Keyword for readout noise in e-"}
string  key_gain    {"GAIN",prompt="Keyword for gain in electrons/ADU"}
real    ron         {3.5,min=0.,prompt="Default readout noise rms in electrons"}
real    gain        {2.2,min=0.00001,prompt="Default gain in e-/ADU"}

string  logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use"}

begin

    # Define local variables
    string  l_inimages, l_outimages, l_outprefix, l_tfunction, l_weights
    string  l_refimages, l_database
    string  l_sci_ext, l_var_ext, l_dq_ext, l_logfile
    string  l_key_ron, l_key_gain, l_background, l_bfunction, l_long_bsample
    int     l_torder, l_nfind, l_tnsum, l_tstep, l_border
    int     l_bnaverage, l_bniterate, l_coloffset
    bool    l_interactive, l_trace, l_recenter, l_verbose, l_find
    bool    l_fl_vardq, l_clean
    real    l_apwidth, l_ron, l_gain, l_lsigma, l_usigma, l_blow_reject
    real    l_bhigh_reject, l_bgrow, l_mos_bsample


    # Other variables used in task
    file    tmpin, tmpout, tmpref, tmpspec, tmpdq, scispec, varspec
    file    dqspec, scilist, mdffile
    struct  sdate
    string  inlist, outlist, reflist, suf, obsmode[200], img, imgnofits, outimg
    string  refimg, tstr, tstr2, refpar, slittype, bsample, l_databasesave
    string  l_logfilesave, l_interpsave, l_nsumsave
    bool    useprefix, mdfexist, useref, useinlist, useoutlist, usereflist
    bool    delete_keywords, stored, l_verbosesave
    int     i, nsciext[200], nim, nout, objectid, num, ybin, nextend
    int     mdfpos[200], nbad, appix
    int     masktype, l_dispaxissave
    int     nin, nref, tint, pixcw, mdfrow, nx, y1, y2, y3, y4, inst, iccd
    real    l_sgain, l_sron, scale, slitlen, objoff, pixscale[2,3], ap1, ap2

    # Constant values:
    real    asecmm = 1.611444
    # Pixel scales [inst,iccd]
    # inst: 1 - GMOS-N, 2 - GMOS-S
    # iccd: 1 - EEV CCDs, 2 - e2vDD CCDs, 3 - Hamamatsu
    pixscale[1,1] = 0.0727
    pixscale[1,2] = 0.07288
    pixscale[1,3] = 0.0807
    pixscale[2,1] = 0.073
    pixscale[2,3] = 0.0800 ##M PIXEL_SCALE

    # Make local variable assignments
    l_inimages = inimages ; l_outimages = outimages ; l_outprefix = outprefix
    l_apwidth = apwidth ; l_interactive = fl_inter ; l_database = database
    l_find = find ; l_trace = trace  # ; l_nfind = nfind
    l_coloffset = coloffset  
    l_tfunction = tfunction ; l_torder = torder ; l_tnsum = tnsum
    l_tstep = tstep ; l_recenter = recenter ; l_weights = weights
    l_clean = clean; l_lsigma = lsigma; l_usigma = usigma
    l_background = background; l_bfunction = bfunction; l_border = border
    l_long_bsample = long_bsample; l_mos_bsample = mos_bsample
    l_bnaverage = bnaverage; l_bniterate = bniterate
    l_blow_reject = blow_reject; l_bhigh_reject = bhigh_reject; l_bgrow = bgrow
    l_fl_vardq = fl_vardq ; l_sci_ext = sci_ext ; l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_key_ron = key_ron ; l_key_gain = key_gain ; l_ron = ron ; l_gain = gain
    l_logfile = logfile ; l_verbose = verbose ; l_refimages = refimages

    l_nfind = 1  # until 2D apall output can be sorted into separate sci exts

    stored = no

    status = 0

    # Make temp files
    tmpin = mktemp("tmpin")
    tmpout = mktemp("tmpout")
    tmpref = mktemp("tmpref")
    scilist = mktemp("tmpscilist")

    # Assign dummy values to tmp file variables
    # Real values will be assigned later
    tmpspec = "dummy"
    tmpdq = "dummy"
    scispec = "dummy"
    varspec = "dummy"
    dqspec = "dummy"
    mdffile = "dummy"

    # Keep some parameters from changing by outside world
    cache ("imgets", "gemhedit", "apall", "gextverify", "gimverify", "tabpar")
    cache ("gemdate", "specred")

    # Check logfile ...
    if (l_logfile == "STDOUT") {
        l_logfile = ""
        l_verbose = yes
    } else if (l_logfile == "" || stridx(" ", l_logfile) > 0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ", l_logfile) > 0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GSEXTRACT: both gsextract.logfile and \
                gmos.logfile are empty.", l_logfile, l_verbose)
            printlog ("                     Using default file gmos.log.",
                l_logfile, l_verbose)
        }
    }

    # Print out the setup to the logfile
    gemdate()
    printlog ("------------------------------------------------------------\
        --------------------", l_logfile, verbose=l_verbose)
    printlog ("GSEXTRACT -- "//gemdate.outdate, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)
    printlog ("inimages     = "//l_inimages, l_logfile, l_verbose)
    printlog ("outimages    = "//l_outimages, l_logfile, l_verbose)
    printlog ("outprefix    = "//l_outprefix, l_logfile, l_verbose)
    printlog ("refimages    = "//l_refimages, l_logfile, l_verbose)
    printlog ("apwidth      = "//l_apwidth, l_logfile, l_verbose)
    printlog ("fl_inter     = "//l_interactive, l_logfile, l_verbose)
    printlog ("find         = "//l_find, l_logfile, l_verbose)
    printlog ("recenter     = "//l_recenter, l_logfile, l_verbose)
    printlog ("trace        = "//l_trace, l_logfile, l_verbose)
    # printlog ("nfind        = "//l_nfind, l_logfile, l_verbose)
    printlog ("coloffset    = "//l_coloffset, l_logfile, l_verbose)
    printlog ("tfunction    = "//l_tfunction, l_logfile, l_verbose)
    printlog ("torder       = "//l_torder, l_logfile, l_verbose)
    printlog ("tnsum        = "//l_tnsum, l_logfile, l_verbose)
    printlog ("tstep        = "//l_tstep, l_logfile, l_verbose)
    printlog ("weights      = "//l_weights, l_logfile, l_verbose)
    printlog ("clean        = "//l_clean, l_logfile, l_verbose)
    if (l_clean) {
        printlog ("lsigma       = "//l_lsigma, l_logfile, l_verbose)
        printlog ("usigma       = "//l_usigma, l_logfile, l_verbose)
    }
    printlog ("background   = "//l_background, l_logfile, l_verbose)
    if (l_background != "none") {
        printlog ("bfunction    = "//l_bfunction, l_logfile, l_verbose)
        printlog ("border       = "//l_border, l_logfile, l_verbose)
        printlog ("long_bsample = "//l_long_bsample, l_logfile, l_verbose)
        printlog ("mos_bsample  = "//l_mos_bsample, l_logfile, l_verbose)
        printlog ("bnaverage    = "//l_bnaverage, l_logfile, l_verbose)
        printlog ("bniterate    = "//l_bniterate, l_logfile, l_verbose)
        printlog ("blow_reject  = "//l_blow_reject, l_logfile, l_verbose)
        printlog ("bhigh_reject = "//l_bhigh_reject, l_logfile, l_verbose)
        printlog ("bgrow        = "//l_bgrow, l_logfile, l_verbose)
    }
    printlog ("fl_vardq     = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("sci_ext      = "//l_sci_ext, l_logfile, l_verbose)
    if (l_fl_vardq) {
        printlog ("var_ext      = "//l_var_ext, l_logfile, l_verbose)
        printlog ("dq_ext       = "//l_dq_ext, l_logfile, l_verbose)
    }
    printlog ("key_ron      = "//l_key_ron, l_logfile, l_verbose)
    printlog ("key_gain     = "//l_key_gain, l_logfile, l_verbose)
    printlog ("ron          = "//l_ron, l_logfile, l_verbose)
    printlog ("gain         = "//l_gain, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    # We start out good :)
    nbad = 0

    # Do the usual basic checks on input parameters and the existence of
    # files, etc.
    tint = fscanf (l_inimages, "%s", tstr)
    if (tstr == "") {
        printlog ("ERROR - GSEXTRACT: input files not specified",
            l_logfile, l_verbose)
        nbad = nbad+1
    }

    # Check existence of input list files
    inlist = "" # so it's not undefined after an err when checking for others
    if (strstr("@",tstr) != 0) {
        tint = stridx ("@", l_inimages)
        inlist = substr (l_inimages, tint+1, strlen(l_inimages))
        useinlist = yes
        if (!access(inlist)) {
            printlog ("ERROR - GSEXTRACT: Input list "//inlist//" not found",
                l_logfile, l_verbose)
            nbad = nbad+1
            inlist = ""
        } else {
            # Reset inlist back to input value
            inlist = tstr
        }
    } else
        useinlist = no

    tint = fscanf (l_refimages, "%s", tstr)
    if (tstr == "")
        useref = no
    else
        useref = yes

    # Check that the reference solution actually exists
    if (useref) {
        if (!access(l_database)) {
            printlog ("ERROR - GSEXTRACT: no database solutions in " // \
                l_database//" for refimages (extract these first)", \
                l_logfile, l_verbose)
            nbad = nbad+1
        }
    }

    # Check existence of ref list files
    reflist = "" # so it's not undefined after an err when checking for others
    if (useref && substr(tstr, 1, 1) == "@") {
        tint = stridx ("@", l_refimages)
        reflist = substr (l_refimages, tint+1, strlen(l_refimages))
        usereflist = yes
        if (!access(reflist)) {
            printlog ("ERROR - GSEXTRACT: Input list "//reflist//" not found",
                l_logfile, l_verbose)
            nbad = nbad+1
            reflist = ""
        }
    } else
        usereflist = no

    # Check whether output filename specified
    useoutlist = no
    tint = fscanf (l_outimages, "%s", tstr)
    if (tstr == "")
        useprefix = yes
    else
        useprefix = no

    # Error if no output filename or prefix
    if (useprefix) {
        tint = fscanf (l_outprefix, "%s %s", tstr, tstr2)
        if (tstr == "" || tint > 1) {
            printlog ("ERROR - GSEXTRACT: outprefix is not specified or \
                contains spaces.", l_logfile, l_verbose)
            nbad = nbad+1
            l_outprefix = ""
        }
        ## Redundant code commented out by JT - remove in next version:
        ## (was it supposed to be outside the if(useprefix)?)
        #
        #  else if(substr(l_outimages, 1, 1) == "@") {
        #     outlist = substr(l_outimages, 2, strlen(l_outimages))
        #     if (!access(outlist)) {
        #       printlog("ERROR - GSEXTRACT: Output list "//outlist//" not \
        #           found", l_logfile, l_verbose)
        #       nbad = nbad + 1
        #     }
        #  }
    }

    # Check if input images exist
    if (useinlist) {
        if (inlist != "")
            sections (inlist, > tmpin)
    } else
        files (l_inimages, sort-, > tmpin)

    nin = 0 # so it's not undefined after an err when checking for other errs
    if (!useinlist || inlist != "") {
        count(tmpin) | scan(nin)
        scanfile = tmpin
        while (fscan(scanfile, img) != EOF) {
            gimverify (img)
            if (gimverify.status > 0) {
                printlog ("ERROR - GSEXTRACT: Input image "//img//" does \
                    not exist or not MEF.", l_logfile, l_verbose)
                nbad = nbad+1
            }
        }
    }

    # Check correct number of reference images exist & construct list
    if (useref && usereflist) {
        if (reflist != "")
            sections ("@"//reflist, > tmpref)
    } else if (useref)
        files (l_refimages, sort-, > tmpref)

    if (useref && (!usereflist || reflist != "")) {
        count(tmpref) | scan(nref)
        scanfile = tmpref
        if (nref != 1 && nref != nin) {
            printlog ("ERROR - GSEXTRACT: Number of refimages does not \
                match number of inimages.", l_logfile, l_verbose)
            nbad = nbad+1
        }
        while (fscan(scanfile, refimg) != EOF) {
            gimverify (refimg)
            if (gimverify.status > 0) {
                printlog ("ERROR - GSEXTRACT: Input image "//refimg//" does \
                    not exist or not MEF.", l_logfile, l_verbose)
                nbad = nbad+1
            }
        }
        # If nref=1 and nin > 1, replicate reference file name in tmpref:
        for (i = nref; i < nin; i += 1)
            print (refimg, >> tmpref)
    } # end if

    scanfile = ""

    # Check if outimages exist and if the the number of input images are equal
    # to the number of output images
    outlist = ""
    if (!useprefix) {
        # Create temporary output list file and check that given list file
        # exists
        tint = fscanf(l_outimages, "%s", tstr)
        if (substr(tstr, 1, 1) == "@") {
            tint = stridx ("@", l_outimages)
            outlist = substr (l_outimages, tint+1, strlen(l_outimages))
            useoutlist = yes
            if (access(outlist))
                sections ("@"//outlist, > tmpout)
            else {
                printlog ("ERROR - GSEXTRACT: Output list "//outlist//" not \
                    found", l_logfile, l_verbose)
                nbad = nbad+1
                outlist = ""
            }
        } else
            files (l_outimages, sort-, > tmpout)

        if (!useoutlist || outlist != "") {
            count(tmpout) | scan(nout)
            scanfile = tmpout
            if (nin != nout) {
                printlog ("ERROR - GSEXTRACT: Number of outfiles does not \
                    match number of infiles.", l_logfile, l_verbose)
                nbad = nbad+1
            }
            while (fscan(scanfile, outimg) != EOF) {
                if (imaccess(outimg)) {
                    printlog ("ERROR - GSEXTRACT: Output image "//outimg//\
                        " exists.", l_logfile, l_verbose)
                    nbad = nbad+1
                }
            }
        }

    } else if (l_outprefix != "" && nin > 0) {
        sections (l_outprefix//"@"//tmpin, > tmpout)
        scanfile = tmpout
        while (fscan(scanfile, outimg) != EOF) {
            if (imaccess(outimg)) {
                printlog ("ERROR - GSEXTRACT: Output image "//outimg//\
                    " exists.", l_logfile, l_verbose)
                nbad = nbad+1
            }
        }
    }

    scanfile = ""

    # Check extension names
    gextverify (l_sci_ext)
    l_sci_ext = gextverify.outext
    if (gextverify.status == 1) {
        printlog ("ERROR - GSEXTRACT: sci_ext is an empty string.",
            l_logfile, l_verbose)
        nbad = nbad+1
    }
    if (l_fl_vardq) {
        gextverify (l_var_ext)
        l_var_ext = gextverify.outext
        if (gextverify.status == 1) {
            printlog ("ERROR - GSEXTRACT: var_ext is an empty string.",
                l_logfile, l_verbose)
            nbad = nbad+1
        }
        gextverify (l_dq_ext)
        l_dq_ext = gextverify.outext
        if (gextverify.status == 1) {
            printlog ("ERROR - GSEXTRACT: dq_ext is an empty string.",
                l_logfile, l_verbose)
            nbad = nbad+1
        }
    }

    # Check that if we want to output the variance spectra that weights is
    # also set to variance
    if (l_fl_vardq && l_weights == "none") {
        printlog ("ERROR - GSEXTRACT: You must set weights=variance for the \
            extraction in order to", l_logfile, l_verbose)
        printlog ("                   output the variance spectrum \
            (fl_vardq=yes).", l_logfile, l_verbose)
        nbad = nbad+1
    }

    # If nbad > 0, go out, else continue
    if (nbad > 0)
        goto error

    # NOW, assuming we're ok so far: let's look over all the input images
    # and make sure they are the right format
    scanfile = tmpin
    i = 0

    while (fscan(scanfile, img) != EOF) {
        i += 1
        gimverify (img)
        img = gimverify.outname//".fits"

        # Which instrument?
        imgets (img//"[0]", "INSTRUME", >& "dev$null")
        if (imgets.value == "0") {
            printlog ("ERROR - GSEXTRACT: Instrument keyword not found.",
                l_logfile, verbose+)
            goto error
        }
        inst = 1 # default is GMOS-N, support for old data
        if (imgets.value == "GMOS-S")
            inst = 2

        # Type of detector?
        imgets (img//"[0]", "DETTYPE", >& "dev$null")
        if (imgets.value == "0") {
            printlog ("ERROR - GSEXTRACT: DETTYPE keyword not found.",
                l_logfile, verbose+)
            goto error
        }

        if (imgets.value == "SDSU II CCD") { # EEV CCDs
            iccd = 1
        } else if (imgets.value == "SDSU II e2v DD CCD42-90") {# New e2vDD CCDs
            iccd = 2
        } else if (imgets.value == "S10892" || imgets.value == "S10892-N") { # Hamamatsu CCDs
            iccd = 3
        }

        imgets (img//"[0]", "MASKTYP", >& "dev$null")
        masktype = int(imgets.value)
        if (masktype != 1) {
            printlog ("ERROR - GSEXTRACT: "//img//" has MASKTYP other than \
                spectral mode.", l_logfile, l_verbose)
            printlog ("                   This task is not meant for image \
                or IFU mode frames.", l_logfile, l_verbose)
            nbad = nbad+1
        }

        tint = fscanf (l_key_gain, "%s", tstr)
        if (tstr != "") {
            tstr = img//"["//l_sci_ext//",1]"
            if (imaccess(tstr)) {
                imgets (tstr, l_key_gain, >& "dev$null")
                if (imgets.value == "0") {
                    printlog ("WARNING - GSEXTRACT: keyword "//l_key_gain//\
                        " not found in "//img, l_logfile, l_verbose)
                    printlog ("                     Using GAIN = "//l_gain,
                        l_logfile, l_verbose)
                    l_sgain = l_gain
                } else
                    l_sgain = real(imgets.value)
            } # (else error later on)
        } else
            l_sgain = l_gain

        tint = fscanf (l_key_ron, "%s", tstr)
        if (tstr != "") {
            tstr = img//"["//l_sci_ext//",1]"
            if (imaccess(tstr)) {
                imgets (tstr, l_key_ron, >& "dev$null")
                if (imgets.value == "0") {
                    printlog ("WARNING - GSEXTRACT: keyword "//l_key_ron//\
                        " not found in "//img, l_logfile, l_verbose)
                    printlog ("                     Using RON = "//l_ron,
                        l_logfile, l_verbose)
                    l_sron = l_ron
                } else
                    l_sron = real(imgets.value)
            } # (else error later on)
        } else
            l_sron = l_ron

        keypar (img//"[0]", "GSSKYSUB", silent+, >& "dev$null")
        if (!keypar.found) {
            keypar (img//"[0]", "GNSSKYSU", silent+, >& "dev$null")
            if (!keypar.found) {
                printlog ("WARNING - GSEXTRACT: Spectra "//img//" have not \
                    been previously sky subtracted.", l_logfile, l_verbose)
                if (l_background == "none") {
                    printlog ("WARNING - GSEXTRACT: Spectra "//img//" will \
                        not be sky subtracted.", l_logfile, l_verbose)
                }
            }
        }

        imgets (img//"[0]", "GSTRANSF", >& "dev$null")
        if (imgets.value == "" || imgets.value == " " || imgets.value == "0")
            printlog ("WARNING - GSEXTRACT: Spectra "//img//" have not been \
                transformed. Run GSTRANSFORM first.", l_logfile, l_verbose)

        # Check for a MDF.  Should be here in the 1st extension, but don't
        # assume that...
        mdfexist = no
        imgets (img//"[0]","NEXTEND", >& "dev$null")
        nextend = int(imgets.value)
        for (j = 1; j <= nextend; j += 1) {
            keypar (img//"["//str(j)//"]", keyword="EXTNAME", silent=yes)
            if (strupr(keypar.value) == "MDF" && keypar.found) {
                mdfexist = yes
                mdfpos[i] = j
            }
        }

        if (!mdfexist) {
            printlog ("ERROR - GSEXTRACT: Input file "//img//" does not have \
                an attached MDF", l_logfile, l_verbose)
            nbad = nbad+1
            imgets (img//"[0]", "MDFFILE", >& "dev$null")
            if ((imgets.value == "") || (imgets.value == " ") || \
                (imgets.value == "0")) {
                printlog ("ERROR - GSEXTRACT: Could not find the MDFFILE \
                    keyword for image "//img, l_logfile, l_verbose)
                nbad = nbad+1
            }
        }
        # Get obsmode
        imgets (img//"[0]", "OBSMODE", >& "dev$null")
        obsmode[i] = imgets.value
        if (imgets.value == "" || imgets.value == " " || imgets.value == "0") {
            printlog ("ERROR - GSSEXTRACT: Could not find the OBSMODE keyword \
                for image "//img, l_logfile, l_verbose)
            nerror = nerror+1
        }

        # Check background sample string is valid, if applicable:
        if (l_background != "none" && obsmode[i] == "LONGSLIT") {
            gemvsample (sample=l_long_bsample, image="img")
            if (gemvsample.status == 2) {
                printlog ("ERROR - GSEXTRACT: invalid sample string \
                    ["//l_long_bsample//"] for "//img, l_logfile, l_verbose)
                nbad = nbad+1
            }
        }

        # Make sure we know how many SCI extensions to expect
        imgets (img//"[0]", "NSCIEXT", >& "dev$null")
        nsciext[i] = int(imgets.value)

        if (nsciext[i] == 0) {
            printlog ("ERROR - GSEXTRACT: Number of SCI extensions unknown \
                for image: "//img, l_logfile, verbose+)
            nbad = nbad+1
        } else {
            for (j = 1; j <= nsciext[i]; j += 1) {
                # Store extension name in temp string:
                tstr = img//"["//l_sci_ext//","//str(j)//"]"

                if (!imaccess(tstr)) {
                    printlog ("ERROR - GSEXTRACT: Could not access "//tstr,
                        l_logfile, l_verbose)
                    nbad = nbad+1
                } else {
                    # Check extension is 2 dimensional:
                    imgets (tstr, "i_naxis", >& "dev$null")
                    if (imgets.value != "2") {
                        printlog ("ERROR - GSEXTRACT: "//tstr//" not \
                            2-dimensional", l_logfile, l_verbose)
                        nbad = nbad+1
                    }

                    # Check if the MDFROW keyword exist inside each extension
                    # (MOS only)
                    if (obsmode[i] == "MOS") {
                        imgets (tstr, "MDFROW", >& "dev$null")
                        if ((imgets.value == "") || (imgets.value == " ") || \
                            (imgets.value == "0")) {
                            printlog ("ERROR - GSEXTRACT: could not find the \
                                MDFROW "//"keyword for image: "//tstr,
                                l_logfile, l_verbose)
                            nbad = nbad+1
                        }
                    }

                    # Check background sample is valid:
                    if (l_background != "none" && obsmode[i] == "LONGSLIT") {
                        gemvsample (sample=l_long_bsample, image=tstr,
                            zeropt="center")
                        if (gemvsample.status == 1) {
                            printlog ("WARNING - GSEXTRACT: background sample \
                                ["//l_long_bsample//"] may be over",
                                l_logfile, l_verbose)
                            printlog ("                     limits for "//\
                                tstr//"\n", l_logfile, l_verbose)
                        }
                    }

                    if (l_fl_vardq) {
                        if (!imaccess(img//"["//l_var_ext//","//str(j)//"]")) {
                            printlog ("ERROR - GSEXTRACT: Could not access "//\
                                img//"["//l_var_ext//","//str(j)//"]",
                                l_logfile, l_verbose)
                            nbad = nbad+1
                        }
                        if (!imaccess(img//"["//l_dq_ext//","//str(j)//"]")) {
                            printlog ("ERROR - GSEXTRACT: Could not access "//\
                                img//"["//l_dq_ext//","//j//"]",
                                l_logfile, l_verbose)
                            nbad = nbad+1
                        }
                    }
                } # end if (image accessible)
            } # end for (j <= nsciext)
        } # end (if got nsciext)

        # Check the binning
        imgets (img//"[0]", "CCDSUM")
        # ybin=int(substr(imgets.value, 3, 3))
        tint = fscanf(imgets.value, "%*d %d", ybin)
        if (tint == 0 || ybin == 0) {
            printlog ("ERROR - GSEXTRACT: Could not determine on-chip binning \
                for image "//img, l_logfile, l_verbose)
            nbad = nbad+1
        }
    } # end of while loop over images

    # Loop over the reference images & check database solutions exist
    if (useref) {
        scanfile = tmpref
        i = 0
        while (fscan(scanfile, refimg) != EOF) {
            i += 1
            gimverify (refimg)
            refimg = gimverify.outname

            for (j = 1; j <= nsciext[i]; j += 1) {
                tstr = l_database//"/ap"//refimg//"_"//l_sci_ext//"_"//j//"_"
                if (!access(tstr)) {
                    printlog ("ERROR - GSEXTRACT: can't find "//tstr,
                        l_logfile, l_verbose)
                    nbad = nbad+1
                }

                # Other checks need to go here, to make sure the input & ref
                # are compatible (2d, same dimensions)
            } # end for j
        }
    }

    # If no error found, continue, else stop
    if (nbad > 0)
        goto error

    scanfile = ""

    # If we made it here then everything must be OK :)
    # The number of science extensions/image is now in the array nsciext[] and
    # the image names (with .fits if needed) is in infile[].

    # Create the scilist, join with the outlist -> proclist
    if (useref) {
        joinlines (tmpin//","//tmpout//","//tmpref, output=scilist, delim=" ",
            missing="", maxchars=161, shortest-, verbose-)
        delete (tmpin//","//tmpout//","//tmpref, verify-, >& "dev$null")
    } else {
        joinlines (tmpin//","//tmpout, output=scilist, delim=" ",missing="",
            maxchars=161, shortest-, verbose-)
        delete (tmpin//","//tmpout, verify-, >& "dev$null")
    }

    scanfile = scilist
    num = 0

    # Store default specred values
    l_dispaxissave = specred.dispaxis
    l_nsumsave = specred.nsum
    l_interpsave = specred.interp
    l_verbosesave = specred.verbose
    l_databasesave = specred.database
    l_logfilesave = specred.logfile
    stored = yes

    specred.dispaxis = 1
    specred.nsum = 10
    specred.interp = "linear"
    specred.verbose = l_verbose
    specred.database = l_database
    specred.logfile = l_logfile

    refimg = ""

    while (fscan(scanfile, img, outimg, refimg) != EOF) {
        nbad = 0
        num = num+1

        # Create tmp FITS file name used within this loop
        mdffile = mktemp("tmpmdf")

        # Make sure we have the .fits extension on the input and output file if
        # it isn't there (fitsutils need it)
        gimverify (img)
        imgnofits = gimverify.outname
        img = gimverify.outname//".fits"

        gimverify (outimg)
        outimg = gimverify.outname//".fits"

        if (useref) {
            gimverify (refimg)
            refimg = gimverify.outname//".fits"
        } else
            refimg = ""

        # From spatial binning in header determine the aperture size in pixels
        # for the extraction, based on the input desired aperture size in
        # arcsec
        imgets (img//"[0]", "CCDSUM")
        tint = fscanf(imgets.value, "%*d %d", ybin)
        scale = pixscale[inst,iccd]*real(ybin)
        l_apwidth = l_apwidth/(2.0*scale)

        # Create output image (copy PHU from input and the MDF)
        tcopy (img//"["//mdfpos[num]//"]", mdffile//".fits", verbose-)
        wmef (mdffile//".fits", outimg, extnames="MDF", verbose-, phu=img,
            >& "dev$null")

        # Loop over the number of SCI extensions for each image
        for (i = 1; i <= nsciext[num]; i += 1) {

            # Create the tmp FITS file names used within this loop.
            tmpspec = mktemp("tmpspec")
            scispec = mktemp("tmpscispec")
            if (l_fl_vardq){
                varspec = mktemp("tmpvarspec")
                dqspec = mktemp("tmpdqspec")
                tmpdq = mktemp("tmpdq")
            }

            # If we're using a reference image, construct the name
            if (refimg == "")
                refpar = ""
            else
                refpar = refimg//"["//l_sci_ext//","//str(i)//"]"

            # Set wavelength reference pixel for centring
            imgets (img//"["//l_sci_ext//","//str(i)//"]", "i_naxis1",
                >& "dev$null")
            if (imgets.value == "0") {
                pixcw = INDEF   # this should never happen
            }
            else {
                nx = int(imgets.value)
                pixcw = nint(real(nx)/2.) + l_coloffset
                pixcw = max(min(pixcw, nx), 1)
            }

            # Determine background sample region if subtracting and slit
            # rectangular:
            if (l_background != "none") {

                # Get MDFROW from the ext header:
                imgets (img//"["//l_sci_ext//","//str(i)//"]", "MDFROW",
                    >& "dev$null")
                mdfrow = int(imgets.value) # NEED ERR CHECK IF MISSING?

                # Determine sample region according to obsmode:
                if (obsmode[num] == "LONGSLIT") {
                    bsample = l_long_bsample
                } else if (obsmode[num] == "MOS") {
                    # Check slit type:
                    tabpar (mdffile//".fits", "slittype", row=mdfrow,
                        >& "dev$null")
                    slittype = tabpar.value

                    # Determine sample region if slit rectangular, else warn
                    # user:
                    if (slittype == "rectangle") {
                        # Determine slit size in pix:
                        tabpar (mdffile//".fits", "slitsize_my", row=mdfrow,
                            >& "dev$null")
                        # (slit length is 1.05 times that defined in the MDF)
                        slitlen = 1.05*real(tabpar.value)*asecmm/scale

                        # Determine object offset along slit, in pix (inverted
                        # in MDF):
                        tabpar (mdffile//".fits", "slitpos_y", row=mdfrow,
                            >& "dev$null")
                        objoff = 0.5*slitlen - real(tabpar.value)/scale

                        # Calculate the points defining the sample & generate
                        # the string: (similar to gsskysub; keep ints for
                        # consistency & simplicity)
                        y1 = int(max(1., 0.5 * ((1.0-l_mos_bsample) * \
                            slitlen)) - objoff)
                        y4 = int(slitlen - y1 - 2 * objoff)
                        y2 = int(max(y1, -l_apwidth))
                        y3 = int(min(y4, l_apwidth))
                        bsample = y1//":"//y2//","//y3//":"//y4

                    } else {
                        printlog ("Slit #"//str(i)//" ; type = "//slittype//\
                            " --- not sky-subtracted", l_logfile, l_verbose)
                        l_background = "none"
                        bsample = ""
                    }

                } # end if (obsmode)
                # end if (l_background != "none")
            } else
                bsample = "*" # dummy value

            printlog ("", l_logfile, no)
            if (l_background == "none")
                printlog ("Extracting spectrum "//img//\
                    "["//l_sci_ext//","//str(i)//"]", l_logfile, l_verbose)
            else
                printlog ("Extracting spectrum "//img//\
                    "["//l_sci_ext//","//str(i)//"] ; backgr=["//bsample//"]",
                    l_logfile, l_verbose)

            if (l_interactive) {
                apall (img//"["//l_sci_ext//","//str(i)//"]", nfind=l_nfind, \
                    output=tmpspec//".fits", apertures="", \
                    format="multispec", references=refpar, profiles="", \
                    interactive=l_interactive, find=l_find, \
                    recenter=l_recenter, resize-, edit+, trace=l_trace, \
                    fittrace+, extract+, extras+, review+, line=pixcw, \
                    nsum=l_tnsum, lower=(-1. * l_apwidth), upper=l_apwidth, \
                    apidtable="", b_function=l_bfunction, b_order=l_border, \
                    b_sample=bsample, b_naverage=l_bnaverage, \
                    b_niterate=l_bniterate, b_low_reject=l_blow_reject, \
                    b_high_reject=l_bhigh_reject, b_grow=l_bgrow, \
                    width=l_apwidth, radius=10., threshold=0., minsep=5., \
                    maxsep=1000., order="increasing", aprecenter="", \
                    npeaks=INDEF, shift+, llimit=INDEF, ulimit=INDEF, \
                    ylevel=0.1, peak+, bkg-, r_grow=0., avglimits-, \
                    t_nsum=l_tnsum, t_step=l_tstep, t_nlost=3, \
                    t_function=l_tfunction, t_order=l_torder, t_sample="*", \
                    t_naverage=1, t_niterate=3, t_low_reject=3., \
                    t_high_reject=3., t_grow=0., background=l_background, \
                    skybox=1, weights=l_weights, pfit="fit1d", clean=l_clean, \
                    saturation=INDEF, readnoise=l_sron, gain=l_sgain, \
                    lsigma=l_lsigma, usigma=l_usigma, nsubaps=1)
            } else {
                apall (img//"["//l_sci_ext//","//str(i)//"]", nfind=l_nfind, \
                    output=tmpspec//".fits", apertures="", \
                    format="multispec", references=refpar, profiles="", \
                    interactive=l_interactive, find=l_find, \
                    recenter=l_recenter, resize-, edit+, trace=l_trace, \
                    fittrace+, extract+, extras+, review+, line=pixcw, \
                    nsum=l_tnsum, lower=(-1. * l_apwidth), upper=l_apwidth, \
                    apidtable="", b_function=l_bfunction, b_order=l_border, \
                    b_sample=bsample, b_naverage=l_bnaverage, \
                    b_niterate=l_bniterate, b_low_reject=l_blow_reject, \
                    b_high_reject=l_bhigh_reject, b_grow=l_bgrow, \
                    width=l_apwidth, radius=10., threshold=0., minsep=5., \
                    maxsep=1000., order="increasing", aprecenter="", \
                    npeaks=INDEF, shift+, llimit=INDEF, ulimit=INDEF, \
                    ylevel=0.1, peak+, bkg-, r_grow=0., avglimits-, \
                    t_nsum=l_tnsum, t_step=l_tstep, t_nlost=3, \
                    t_function=l_tfunction, t_order=l_torder, t_sample="*", \
                    t_naverage=1, t_niterate=3, t_low_reject=3., \
                    t_high_reject=3., t_grow=0., background=l_background, \
                    skybox=1, weights=l_weights, pfit="fit1d", clean=l_clean, \
                    saturation=INDEF, readnoise=l_sron, gain=l_sgain, \
                    lsigma=l_lsigma, usigma=l_usigma, nsubaps=1, >& "dev$null")
            }

            if (fl_vardq) {
                # Run apall without weights for DQ plane
                # Use aperture obtained from science extension
                # Otherwise the parameters used are the same as above
                apall (img//"["//l_dq_ext//","//str(i)//"]", nfind=l_nfind, \
                    output=tmpdq, apertures="", format="multispec", \
                    references=img//"["//l_sci_ext//","//str(i)//"]", \
                    profiles="", interactive-, find-, recenter-, resize-, \
                    edit-, trace-, fittrace-, extract+, extras-, review-, \
                    line=INDEF, nsum=10, lower=-5., upper=5., apidtable="", \
                    b_function="chebyshev", b_order=1, b_sample=bsample, \
                    b_naverage=-3, b_niterate=0, b_low_reject=3, \
                    b_high_reject=3, b_grow=0, width=5, radius=10, \
                    threshold=0, minsep=5, maxsep=1000, order="increasing", \
                    aprecenter="", npeaks=INDEF, shift+, llimit=INDEF, \
                    ulimit=INDEF, ylevel=0.1, peak+, bkg-, r_grow=0., \
                    avglimits-, t_nsum=l_tnsum, t_step=l_tstep, t_nlost=3, \
                    t_function=l_tfunction, t_order=l_torder, t_sample="*", \
                    t_naverage=1, t_niterate=3, t_low_reject=3., \
                    t_high_reject=3., t_grow=0., background="none", skybox=1, \
                    weights="none", pfit="fit1d", clean-, saturation=INDEF, \
                    readnoise=0, gain=1, lsigma=l_lsigma, usigma=l_usigma, \
                    nsubaps=1)
            }

            if (!imaccess(tmpspec//".fits")) {
                printlog ("WARNING - GSEXTRACT: apall failed for "//\
                    img//"["//l_sci_ext//","//str(i)//"]", l_logfile, verbose+)
                nbad += 1
                flprcache
            } else {
                # Unfortunately, fxinsert isn't friendly enough to allow us to
                # insert a specified dimension of a multispec fits file into a
                # MEF extension.  So, we must separate out the science and
                # variance spectra of the multispec file first
                # 2005oct18 - logic fixed by BM
                delete_keywords = no
                if (l_fl_vardq) {
                    imcopy(tmpspec//".fits[*,*,1]", scispec, verbose-)
                    delete_keywords = yes
                    if (l_background == "none") {
                        imcopy(tmpspec//".fits[*,*,3]", varspec, verbose-)
                    } else {
                        imcopy(tmpspec//".fits[*,*,4]", varspec, verbose-)
                    }
                    # Output is a sigma spectrum, need to square
                    imarith(varspec//".fits", "*", varspec//".fits", \
                        varspec//".fits", verbose-)
                } else {
                    if (l_background == "none") {
                        imcopy(tmpspec//".fits", scispec, verbose-)
                    } else {
                        imcopy(tmpspec//".fits[*,*,1]", scispec, verbose-)
                        delete_keywords = yes
                    }
                }

                # Remove the keywords associated with the multispec fits file
                if (delete_keywords) {
                    gemhedit(scispec//".fits", "CD3_3", "", "", delete+)
                    gemhedit(scispec//".fits", "CRPIX3", "", "", delete+)
                    gemhedit(scispec//".fits", "CTYPE3", "", "", delete+)
                    gemhedit(scispec//".fits", "LTM3_3", "", "", delete+)
                    gemhedit(scispec//".fits", "LTV3", "", "", delete+)
                    gemhedit(scispec//".fits", "WAT3_001", "", "", delete+)
                    gemhedit(scispec//".fits", "WAXMAP01", "", "", delete+)
                    # Edit WCSDIM so that it equals 1 instead of 3
                    gemhedit(scispec//".fits", "WCSDIM", "1", "", delete-)
                }

                # Insert extracted spectra into outfile
                imcopy(scispec//".fits", \
                    outimg//"["//l_sci_ext//","//str(i)//",append]", \
                    verbose-, >& "dev$null")

                if (l_fl_vardq) {
                    # Remove the keywords associated with the multispec fits
                    # file that causes gscalibrate to not update the variance
                    # extension
                    gemhedit(varspec//".fits", "CD3_3", "", "", delete+)
                    gemhedit(varspec//".fits", "CRPIX3", "", "", delete+)
                    gemhedit(varspec//".fits", "CTYPE3", "", "", delete+)
                    gemhedit(varspec//".fits", "LTM3_3", "", "", delete+)
                    gemhedit(varspec//".fits", "LTV3", "", "", delete+)
                    gemhedit(varspec//".fits", "WAT3_001", "", "", delete+)
                    gemhedit(varspec//".fits", "WAXMAP01", "", "", delete+)
                    # Edit WCSDIM so that it equals 1 instead of 3
                    gemhedit(varspec//".fits", "WCSDIM", "1", "", delete-)

                    # Add variance plane
                    imcopy(varspec//".fits",
                        outimg//"["//l_var_ext//","//str(i)//",append]",
                        verbose-, >& "dev$null")

                    # Now get the DQ plane. If more than 50% of the pixels
                    # in the aperture are bad, mark the output pixel in
                    # the 1D DQ plane as bad

                    # First, get the aperture size

                    match ("low", l_database //"/ap"// imgnofits // "_" // \
                        l_dq_ext // "_" // str(i) // "_", stop-) | \
                        fields ("STDIN", fields="3", lines="1", \
                        quit_if_miss=yes, print_file_n=no) | scan (ap1)
                    match ("high", l_database //"/ap"// imgnofits // "_" // \
                        l_dq_ext // "_" // str(i) // "_", stop-) | \
                        fields ("STDIN", fields="3", lines="1", \
                        quit_if_miss=yes, print_file_n=no) | scan (ap2)
                    appix = int(0.5 * (ap2 - ap1))

                    # Then replace pixels

                    imreplace (tmpdq, 0, imaginary=0., lower=INDEF,
                        upper=appix, radius=0.)
                    imreplace (tmpdq, 1, imaginary=0., lower=appix,
                        upper=INDEF, radius=0.)
                    imstack(tmpdq,
                        outimg//"["//l_dq_ext//","//str(i)//",append]", \
                        pixtype="short", >& "dev$null")
                    imdelete (tmpdq, verify-, >& "dev$null")
                }
            }

            # end of for-loop over SCI extensions
            # (check file exists, just to be on the safe side!)
            if (imaccess(tmpspec//".fits"))
                imdelete (tmpspec//".fits", verify-, >& "dev$null")
            if (imaccess(scispec//".fits"))
                imdelete (scispec//".fits", verify-, >& "dev$null")

            if (l_fl_vardq) {
                if (imaccess(varspec//".fits"))
                    imdelete (varspec//".fits", verify-, >& "dev$null")
                if (imaccess(dqspec//".fits"))
                    imdelete (dqspec//".fits", verify-, >& "dev$null")
            }
        }

        # Final header update
        gemdate ()
        gemhedit (outimg//"[0]", "GSEXTRAC", gemdate.outdate,
            "UT Time stamp for GSEXTRACT", delete-)
        gemhedit (outimg//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
        gemhedit (outimg//"[0]", "NSCIEXT", nsciext[num],
            "Number of science extensions", delete-)
        if (imaccess(mdffile))
            imdelete (mdffile//".fits", verify-, >& "dev$null")

        if (nbad != 0)
            imdelete (outimg, verify-, >& "dev$null")

    } # end of while-loop over images

    goto clean

error:
    status = 1
    delete (tmpin, verify-, >& "dev$null")
    delete (tmpout, verify-, >& "dev$null")
    delete (tmpref, verify-, >& "dev$null")
    imdelete (tmpspec, verify-, >& "dev$null")
    imdelete (scispec, verify-, >& "dev$null")
    imdelete (varspec, verify-, >& "dev$null")
    imdelete (dqspec, verify-, >& "dev$null")
    imdelete (mdffile, verify-, >& "dev$null")
    printlog ("ERROR - "//nbad//" fatal error(s) found.", l_logfile, l_verbose)
    goto clean

clean:
    # Restore default specred values
    if (stored) {
        specred.dispaxis = l_dispaxissave
        specred.nsum = l_nsumsave
        specred.interp = l_interpsave
        specred.verbose = l_verbosesave
        specred.database = l_databasesave
        specred.logfile = l_logfilesave
    }

    # Close log file
    if (access(scilist))
        delete (scilist, verify-, >& "dev$null")
    printlog (" ", l_logfile, l_verbose)
    gemdate()
    printlog ("GSEXTRACT: Done -- "//gemdate.outdate, l_logfile, l_verbose)
    if (status == 0) {
        printlog ("GSEXTRACT: Exit status GOOD", l_logfile, l_verbose)
    } else {
        printlog ("GSEXTRACT: Exit status ERROR", l_logfile, l_verbose)
    }
    printlog ("------------------------------------------------------------\
        --------------------", l_logfile, l_verbose )
    scanfile = ""

end
