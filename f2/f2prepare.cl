# Copyright(c) 2010-2017 Association of Universities for Research in Astronomy, Inc.

# F2PREPARE - Prepare raw FLAMINGOS-2 data for reduction
# Full details of this task are given in the help file

procedure f2prepare(inimages)

char    inimages        {prompt = "Input FLAMINGOS-2 image(s)"}
char    rawpath         {"", prompt = "Path for input raw images"}
char    outimages       {"", prompt = "Output image(s)"}
char    outprefix       {"f", prompt = "Prefix for output image(s)"}
char    bpm             {"", prompt = "Bad pixel mask file"}

bool    fl_vardq        {yes, prompt = "Create variance and data quality frames?"}
bool    fl_correct      {no, prompt = "Correct for non-linearity in the data? (Not yet implemented)"}
bool    fl_saturated    {yes, prompt = "Flag saturated pixels in DQ? "}
bool    fl_nonlinear    {yes, prompt = "Flag non-linear pixels in DQ?"}

char    arraytable      {"f2$data/f2array.fits", prompt = "Table of array settings"}

bool    fl_addmdf       {yes, prompt = "Add Mask Definition File (MOS mode)"}
char    key_mdf         {"MASKNAME", prompt = "Header keyword for the Mask Definition File"}
char    mdffile         {"", prompt = "Mask Definition File to use if keyword not found"}
char    mdfdir          {"f2$data/", prompt ="Mask Definition File database directory"}
bool    fl_dark_mdf     {no, prompt = "Force the attachment of MDFs to dark frames?"}

char    logfile         {"", prompt = "Logfile"}
bool    verbose         {yes, prompt = "Verbose"}
int     status          {0, prompt = "Exit status (0=good)"}

struct  *scanfile       {"", prompt = "Internal use only"}

begin

    char    l_inimages = ""
    char    l_rawpath = ""
    char    l_outimages = ""
    char    l_outprefix = ""
    char    l_bpm = ""
    bool    l_fl_vardq
    bool    l_fl_correct
    bool    l_fl_saturated
    bool    l_fl_nonlinear
    char    l_arraytable = ""
    bool    l_fl_addmdf
    char    l_key_mdf = ""
    char    l_mdffile = ""
    char    l_mdfdir = ""
    bool    l_fl_dark_mdf = ""
    char    l_logfile = ""
    bool    l_verbose

    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_ron = ""
    char    l_key_gain = ""
    char    l_key_sat = ""
    char    l_key_nonlinear = ""
    char    l_key_filter = ""
    char    l_key_lnrs = ""
    char    l_key_camera = ""
    char    l_key_fpmask = ""
    char    l_key_dispaxis = ""
    char    l_key_section = ""
    char    l_key_pixscale = ""
    char    l_key_xoff = ""
    char    l_key_yoff = ""
    char    l_key_date = ""

    char    l_filter, l_filter1, l_filter2, l_temp, tmpfile, tmpfile2, tmparray
    char    in[1000], inpath[1000], out[1000], expr, keyfound
    char    tmpsci, tmpvar, tmpdq, tmpinimage, tmpheader, tmpchk1, tmpchk2
    int     i, nimages, noutimages, maxfiles, nbad, junk, lines1, lines2
    char    l_varexpression, l_tempexpression, badhdr, lastchar
    char    l_lyot, l_mospos, l_dateobs, l_detid, suf, mdfname[200]
    real    l_ron, l_gain, l_well
    int     l_lnrs, l_sat, l_linlimit, len, linpix
    real    x, l_coeff1, l_coeff2, l_coeff3
    real    l_cd11, l_cd12, l_cd21, l_cd22, l_pixscale
    real    l_linearlimit, l_nonlinlimit, l_xoff, l_yoff
    bool    alreadyfixed[1000], bad, addmdf[200], longslit, mos, update_keyword
    int     bpmxsize, bpmysize, xsize, ysize
    struct  l_struct

    status = 0
    maxfiles = 1000
    noutimages = 0

    # Cache tasks used throughout the script
    cache ("keypar", "gemdate")

    # Set the local variables
    l_inimages = inimages ; l_rawpath = rawpath ; l_outimages = outimages
    l_outprefix = outprefix ; l_bpm = bpm ; l_fl_vardq = fl_vardq
    l_fl_correct = fl_correct ; l_fl_saturated = fl_saturated
    l_fl_nonlinear = fl_nonlinear ; l_arraytable = arraytable
    l_fl_addmdf = fl_addmdf ; l_key_mdf = key_mdf ; l_mdffile = mdffile
    l_mdfdir = mdfdir ; l_fl_dark_mdf = fl_dark_mdf ; l_logfile = logfile
    l_verbose = verbose

    # Open temp files
    tmpfile = mktemp ("tmp1")
    tmpfile2 = mktemp ("tmp2")
    tmparray = mktemp ("tmparray")//".fits"
    tmpchk1 = mktemp ("tmpchk1")
    tmpchk2 = mktemp ("tmpchk2")
    
    # Assign other tmp files variables dummy names.
    # Just need to give the variable a value. Real names will be created later.
    tmpsci = "dummy"
    tmpvar = "dummy"
    tmpdq = "dummy"
    tmpinimage = "dummy"

    # Shared definitions, define parameters from nsheaders
    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_ron, l_key_ron)
    if ("" == l_key_ron) badhdr = badhdr + " key_ron"
    junk = fscan (nsheaders.key_gain, l_key_gain)
    if ("" == l_key_gain) badhdr = badhdr + " key_gain"
    junk = fscan (nsheaders.key_sat, l_key_sat)
    if ("" == l_key_sat) badhdr = badhdr + " key_sat"
    junk = fscan (nsheaders.key_nonlinear, l_key_nonlinear)
    if ("" == l_key_nonlinear) badhdr = badhdr + " key_nonlinear"
    junk = fscan (nsheaders.key_filter, l_key_filter)
    if ("" == l_key_filter) badhdr = badhdr + " key_filter"
    junk = fscan (nsheaders.key_lnrs, l_key_lnrs)
    if ("" == l_key_lnrs) badhdr = badhdr + " key_lnrs"
    junk = fscan (nsheaders.key_camera, l_key_camera)
    if ("" == l_key_camera) badhdr = badhdr + " key_camera"
    junk = fscan (nsheaders.key_fpmask, l_key_fpmask)
    if ("" == l_key_fpmask) badhdr = badhdr + " key_fpmask"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (nsheaders.key_section, l_key_section)
    if ("" == l_key_section) badhdr = badhdr + " key_section"
    junk = fscan (nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"
    junk = fscan (nsheaders.key_xoff, l_key_xoff)
    if ("" == l_key_xoff) badhdr = badhdr + " key_xoff"
    junk = fscan (nsheaders.key_yoff, l_key_yoff)
    if ("" == l_key_yoff) badhdr = badhdr + " key_yoff"
    junk = fscan (nsheaders.key_date, l_key_date)
    if ("" == l_key_date) badhdr = badhdr + " key_date"

    #----------------------------------------------------------------------
    # Test that name of logfile makes sense
    if (l_logfile == "" || l_logfile == " ") {
        junk = fscan (f2.logfile, l_logfile) 
        if (l_logfile == "" || l_logfile == " ") {
            l_logfile = "f2.log"
            printlog ("WARNING - F2PREPARE: Both f2prepare.logfile and \
                f2.logfile are undefined.", l_logfile, verbose+)
            printlog ("                  Using " // l_logfile, l_logfile, 
                verbose+)
        }
    }
    
    # Open logfile
    date | scan (l_struct)
    printlog ("--------------------------------------------------------------\
        --------------", l_logfile, l_verbose)
    printlog ("F2PREPARE -- "// l_struct, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    if ("" != badhdr) {
        printlog ("ERROR - F2PREPARE: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        status = 1
        goto clean
    }
    if (l_fl_vardq && ! l_fl_saturated) {
        printlog ("WARNING - F2PREPARE: Saturated pixels not flagged.", \
            l_logfile, verbose+) 
    }
    if (l_fl_vardq && ! l_fl_nonlinear) {
        printlog ("WARNING - F2PREPARE: Non-linear pixels not flagged.", \
            l_logfile, verbose+) 
    }
    if (! l_fl_vardq && l_fl_saturated && l_fl_nonlinear) {
        printlog ("WARNING - F2PREPARE: Cannot flag saturated and non-linear \
            pixels if", l_logfile, verbose+) 
        printlog ("                     fl_vardq=no. Setting fl_vardq=yes.", \
            l_logfile, verbose+)
        l_fl_vardq = yes
    }
    else if (! l_fl_vardq && l_fl_saturated) {
        printlog ("WARNING - F2PREPARE: Cannot flag saturated pixels if \
            fl_vardq=no.", l_logfile, verbose+) 
        printlog ("                     Setting fl_vardq=yes.", \
            l_logfile, verbose+)
        l_fl_vardq = yes
    }
    else if (! l_fl_vardq && l_fl_nonlinear) {
        printlog ("WARNING - F2PREPARE: Cannot flag non-linear pixels if \
            fl_vardq=no.", l_logfile, verbose+) 
        printlog ("                     Setting fl_vardq=yes.", \
            l_logfile, verbose+)
        l_fl_vardq = yes
    }
    if (! l_fl_correct) {
        printlog ("WARNING - F2PREPARE: Not correcting for non-linearity.", \
            l_logfile, verbose+)
    }

    printlog ("", l_logfile, l_verbose)

    # Test the MDF related parameters
    if (l_fl_addmdf) {
        if (((l_key_mdf == "") || (stridx (" ", l_key_mdf) > 0)) && \
            ((l_mdffile == "") || (stridx (" ", l_mdffile) > 0))) {
            
            printlog ("ERROR - F2PREPARE: Neither the MDF keyword key_mdf or \
                the", l_logfile, verbose+)
            printlog ("                  MDF filename mdffile are defined",
                l_logfile, verbose+)
            status = 1
            goto clean
        }
    }

    #----------------------------------------------------------------------
    # Ensure that rawpath has a trailing slash. Don't add one if rawpath is an
    # iraf variable (ends with a "$") or if it already ends with a trailing
    # slash. 
    lastchar = substr (l_rawpath, (strlen (l_rawpath)), (strlen (l_rawpath)))
    if (lastchar != "/" && lastchar != "$")
        l_rawpath = l_rawpath //"/"
    if (l_rawpath == "/" || l_rawpath == " ")
        l_rawpath = ""

    #----------------------------------------------------------------------
    # Load up arrays of input name lists
    # This version handles both *s, / and commas in l_inimages

    # Check that list file exists
    if (substr (l_inimages, 1, 1) == "@") {
        l_temp = substr (l_inimages, 2, strlen (l_inimages))
        if (!access (l_temp) && !access (l_rawpath // l_temp)) {
            printlog ("ERROR - F2PREPARE: Input file "// l_temp //" not \
                found.", l_logfile, verbose+)
            status = 1
            goto clean
        }
    }

    # Parse wildcard and comma-separated lists
    if (substr (l_inimages, 1, 1) == "@") 
        scanfile = substr (l_inimages, 2, strlen (l_inimages))
    else {
        print (l_inimages, > tmpfile)
        sed ("-e", 's/\,/\ /g', tmpfile) | words ("STDIN", > tmpfile2)
        scanfile = tmpfile2
        delete (tmpfile, verify-, >& "dev$null")
    }
    while (fscan (scanfile, l_temp) != EOF) {
        files (l_rawpath // l_temp, sort-) | unique ("STDIN", >> tmpfile)
    }
    delete (tmpfile2, verify-, >& "dev$null")
    scanfile = tmpfile

    nimages = 0
    nbad = 0
    while (fscan (scanfile, l_temp) != EOF) {
        bad = no
        # Remove rawpath
        l_temp = substr (l_temp, strlen (l_rawpath) + 1, strlen (l_temp))

        # Remove .fits if present
        if (substr (l_temp, strlen (l_temp) - 4, strlen (l_temp)) == ".fits")
            l_temp = substr (l_temp, 1, (strlen (l_temp) - 5))

        # Check to see if the file exists (or if it's not MEF FITS)
        gimverify (l_rawpath // l_temp)
        if (gimverify.status == 1) {
            printlog ("ERROR - F2PREPARE: Input image "// l_rawpath // \
                l_temp //" not found.", l_logfile, verbose+)
            nbad += 1
            bad = yes
        } else if (gimverify.status > 1) {
            printlog ("ERROR - F2PREPARE: Input image not multiple-extension \
                FITS.", l_logfile, verbose+)
            nbad += 1
            bad = yes
        } else {
            nimages = nimages + 1
            if (nimages > maxfiles) {
                printlog ("ERROR - F2PREPARE: Maximum number of input images \
                    exceeded ("// str (maxfiles) //")", l_logfile, verbose+)
                status = 1
                goto clean
            }
            # Check to see if already fixed with f2prepare
            if (!bad) {
                keyfound = ""
                hselect (l_rawpath // l_temp //"[0]", "*PREPAR*", yes) | \
                    scan (keyfound)
                if (keyfound != "") {
                    alreadyfixed[nimages] = yes
                    printlog ("WARNING - F2PREPARE: Image \
                        "// l_rawpath // l_temp //" already fixed using \
                        *PREPARE.", l_logfile, verbose+) 
                    printlog ("                    Data scaling, correction \
                        for number of reads, and header", l_logfile, verbose+)
                    printlog ("                    updating not performed.",
                        l_logfile, verbose+)
                } else
                    alreadyfixed[nimages] = no
                in[nimages] = l_temp
            } # end if (!bad)

            # Trim off path from in[i], if present
            i = strlen (in[nimages]) - 1
            inpath[nimages] = ""
            while (i > 1) {
                if (substr (in[nimages], i, i) == "/") {
                    inpath[nimages] = substr (in[nimages], 1, i)
                    in[nimages] = substr (in[nimages], i + 1, \
                        strlen (in[nimages]))
                    i = 0
                }
                i = i - 1
            } # end while (i > 1) path-chopping loop
        } # end else
    } # end while (fscan) loop

    # Exit if problems found
    if (nbad > 0) {
        printlog ("ERROR - F2PREPARE: "// nbad //" image(s) either do not \
            exist or are not MEF files.", l_logfile, verbose+)
        status = 1
        goto clean
    }

    printlog ("Processing "// nimages //" files", l_logfile, l_verbose)
    scanfile = ""
    delete (tmpfile, verify-, >& "dev$null")

    # Now for the output images
    # Outimages could contain legal * if it is of a form like %st%stX%*.imh

    nbad = 0
    print (l_outimages) | scan (l_outimages)
    if (l_outimages != "" && l_outimages != " ") {
        if (substr (l_outimages, 1, 1) == "@")
            scanfile = substr (l_outimages, 2, strlen (l_outimages))
        else {
            files (l_outimages, sort-) | unique ("STDIN", > tmpfile)
            scanfile = tmpfile
        }

        while (fscan (scanfile, l_temp) != EOF) {
            # Remove .fits if present
            if (substr (l_temp, strlen (l_temp) - 4, strlen (l_temp)) == \
                ".fits")
                l_temp = substr (l_temp, 1, (strlen (l_temp) - 5))
                noutimages = noutimages + 1
                if (noutimages > maxfiles) {
                    printlog ("ERROR - F2PREPARE: Maximum number of output \
                        images exceeded ("// str (maxfiles) //")", l_logfile, \
                        verbose+)
                    status = 1 
                    goto clean
                }
                out[noutimages] = l_temp 
                if (imaccess (out[noutimages])) {
                    printlog ("ERROR - F2PREPARE: Output image \
                        "// out[noutimages] //" exists", l_logfile, verbose+)
                    nbad += 1
                }
        } # end while
    }
    scanfile = ""
    delete (tmpfile, verify-, >& "dev$null")

    # Exit if problems found
    if (nbad > 0) {
        printlog ("ERROR - F2PREPARE: "// nbad //" image(s) already exist.",
            l_logfile, verbose+)
        status = 1
        goto clean
    }

    # If there are too many or too few output images, and any defined
    # at all at this stage - exit with error
    if (nimages != noutimages && l_outimages != "") {
        printlog ("ERROR - F2PREPARE: Number of input and output images \
            unequal.", l_logfile, verbose+)
        status = 1
        goto clean
    }

    # If prefix is to be used instead
    if (l_outimages == "" || l_outimages == " ") {
        print (l_outprefix) | scan (l_outprefix)
        if (l_outprefix == "" || l_outprefix == " ") {
            printlog ("ERROR - F2PREPARE: Neither output image name nor \
                output prefix is defined.", l_logfile, verbose+)
            status = 1
            goto clean
        }
        i = 1
        nbad = 0
        while (i <= nimages) {
            out[i] = l_outprefix // in[i]
            if (imaccess (out[i])) {
                printlog ("ERROR - F2PREPARE: Output image "// out[i] //" \
                    already exists.", l_logfile, verbose+)
                nbad += 1
            }
            i = i + 1
        }
        if (nbad > 0) {
            printlog ("ERROR - F2PREPARE: "// nbad //" image(s) already \
                exist.", l_logfile, verbose+)
            status = 1
            goto clean
        }
    }

    # Check there are no duplications in the output file
    i = 1
    while (i <= nimages) {
        print (out[i], >> tmpchk1)
        i = i + 1
    }
    unique (tmpchk1, > tmpchk2)
    count (tmpchk1) | scan (lines1)
    count (tmpchk2) | scan (lines2)
    if (lines1 != lines2) {
        printlog ("ERROR - F2PREPARE: One or more output files have the same \
            name.", l_logfile, verbose+)
        status = 1
        goto clean
    }
    delete (tmpchk1 // "," // tmpchk2, verify-, >& "dev$null")

    #----------------------------------------------------------------------
    # Check for existence of bad pixel mask
    # Still need to add check for .pl or .fits, and convert .pl to fits
    # if required by addmasks

    if (l_fl_vardq && !imaccess (l_bpm) && l_bpm != "" && \
        stridx (" ", l_bpm) <= 0) {
        printlog ("WARNING - F2PREPARE: Bad pixel mask "// l_bpm //" not \
            found.", l_logfile, verbose+)
        printlog ("                     Not using BPM to generate data \
            quality plane.", l_logfile, verbose+)
        l_bpm = ""
    } else if (l_fl_vardq && (l_bpm == "" || stridx (" ", l_bpm) > 0)) {
        printlog ("WARNING - F2PREPARE: Bad pixel mask is either an empty \
            string or contains", l_logfile, verbose+)
        printlog ("                     spaces. Not using BPM to generate \
            data quality plane.", l_logfile, verbose+)
        l_bpm = ""
    }
    if (l_bpm == "")
        l_bpm = "none"

    # Check to make sure BPM and input images are same size
    if (l_bpm != "none") {
        # Determine size of BPM
        keypar (l_bpm, "i_naxis1", silent+, >& "dev$null")
        if (keypar.found)
            bpmxsize = int (keypar.value)
        keypar (l_bpm, "i_naxis2", silent+, >& "dev$null")
        if (keypar.found)
            bpmysize = int (keypar.value)

        # Determine size of input image
        keypar (l_rawpath // inpath[1] // in[1] //"[1]", "i_naxis1", silent+, \
            >& "dev$null")
        if (keypar.found)
            xsize = int (keypar.value)
        keypar (l_rawpath // inpath[1] // in[1] //"[1]", "i_naxis2", silent+, \
            >& "dev$null")
        if (keypar.found)
            ysize = int (keypar.value)

        if (bpmxsize != xsize || bpmysize != ysize) {
            printlog ("WARNING - F2PREPARE: Input images and BPM are not the \
                same size.", l_logfile, verbose+)
            printlog ("                    Not using BPM to generate data \
                quality planes.", l_logfile, verbose+)
            l_bpm = "none"
        }
    }

    #----------------------------------------------------------------------
    # Start output
    printlog (" ", l_logfile, l_verbose)
    printlog ("  n      input file -->      output file", l_logfile, l_verbose)
    printlog ("                 filter     focal plane       input BPM   RON  \
        gain     sat", l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    # The loop:  create VAR and DQ if fl_vardq=yes
    #            determine read noise, gain, saturation, etc.

    i = 1
    while (i <= nimages) {
        update_keyword = no

        # Create tmp FITS file names used within this loop
        tmpsci = mktemp ("tmpsci")
        tmpvar = mktemp ("tmpvar")
        tmpdq = mktemp ("tmpdq")
        tmpinimage = mktemp ("tmpinimage")

        in[i] = l_rawpath // inpath[i] // in[i]

        #------------------------------------------------------------------
        # Read array info from table

        keypar (in[i] // "[0]", l_key_lnrs, silent+, >& "dev$null")
        if (keypar.found) {
            l_lnrs = int (keypar.value)
        } else {
            printlog ("ERROR - F2PREPARE: Could not read number of \
                non-destructive reads from header", l_logfile, verbose+)
            printlog ("                  of file "// in[i], l_logfile, \
                verbose+)
            status = 1
            goto clean
        }

        if (no == access (l_arraytable)) {
            printlog ("ERROR - F2PREPARE: Array table "// l_arraytable \
                //" not found.", l_logfile, verbose+) 
            status = 1
            goto clean
        } else {

            keypar (in[i] //"[0]", "DETID", silent+, >& "dev$null")
            if (keypar.found) {
                l_detid = int(keypar.value)
            } else {
                printlog ("ERROR - F2PREPARE: Could not read the detector ID \
                    from the header of file", l_logfile, verbose+)
                printlog ("                   "//in[i], l_logfile, verbose+)
                status = 1
                goto clean
            }

            tdelete (tmparray, go_ahead+, verify-, >& "dev$null")

            expr = "detid = "//l_detid //" && lnrs = "// l_lnrs
            tselect (intable=l_arraytable, outtable=tmparray, expr=expr)

            tinfo (table=tmparray, >& "dev$null")
            if (tinfo.nrows != 1) {
                printlog ("WARNING - F2PREPARE: No array calibration data \
                    for number of reads equal to "// l_lnrs, l_logfile, \
                    verbose+)
                printlog ("                   in "// l_arraytable //".", \
                    l_logfile, verbose+) 
                printlog ("                   Assuming number of reads is \
                    equal to 1", l_logfile, verbose+)
                tdelete (tmparray, go_ahead+, verify-, >& "dev$null")
                expr = "detid = "//l_detid //" && lnrs = 1"
                tselect (intable=l_arraytable, outtable=tmparray, expr=expr)
            }

            tabpar (table=tmparray, column="readnoise", row=1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - F2PREPARE: No readnoise for number of \
                    reads equal to "// l_lnrs, l_logfile, verbose+)
                printlog ("                   in "// l_arraytable //".", \
                    l_logfile, verbose+)
                status = 1
                goto clean
            } else
                l_ron = real (tabpar.value)

            tabpar (table=tmparray, column="gain", row=1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - F2PREPARE: No gain for number of reads \
                    equal to "// l_lnrs, l_logfile, verbose+)
                printlog ("                   in "// l_arraytable //".", \
                    l_logfile, verbose+) 
                status = 1
                goto clean
            } else
                l_gain = real (tabpar.value)

            tabpar (table=tmparray, column="well", row=1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - F2PREPARE: No well defined in \
                    "// l_arraytable //".", l_logfile, verbose+) 
                status = 1
                goto clean
            } else
                l_well = real (tabpar.value)

            tabpar (table=tmparray, column="linearlimit", row=1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - F2PREPARE: No linearlimit defined in \
                    "// l_arraytable //".", l_logfile, verbose+) 
                status = 1
                goto clean
            } else
                l_linearlimit = real (tabpar.value)

            tabpar (table=tmparray, column="nonlinlimit", row=1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - F2PREPARE: No nonlinlimit defined in \
                    "// l_arraytable //".", l_logfile, verbose+) 
                status = 1
                goto clean
            } else
                l_nonlinlimit = real (tabpar.value)

            if (l_fl_correct) {

                tabpar (table=tmparray, column="coeff1", row=1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2PREPARE: No coeff1 defined in \
                        "// l_arraytable //".", l_logfile, verbose+) 
                    status = 1
                    goto clean
                } else
                    l_coeff1 = real (tabpar.value)

                tabpar (table=tmparray, column="coeff2", row=1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2PREPARE: No coeff2 defined in \
                        "// l_arraytable //".", l_logfile, verbose+)
                    status = 1
                    goto clean
                } else
                    l_coeff2 = real (tabpar.value)

                tabpar (table=tmparray, column="coeff3", row=1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2PREPARE: No coeff3 defined in \
                        "// l_arraytable //".", l_logfile, verbose+)
                    status = 1
                    goto clean
                } else
                    l_coeff3 = real (tabpar.value)
            } 
        } # end else

        # The well depth (l_well) is in electrons. The gain (l_gain) is in
        # electrons / ADU. Therefore, the saturation (l_sat) is in ADU. 
        l_sat = int (l_well/l_gain)
        if (l_fl_correct)
            l_linlimit = int (l_sat*l_linearlimit)
        else
            l_linlimit = int (l_sat*l_nonlinlimit)

	if (alreadyfixed[i]) {
            # Copy SCI plane of already-prepared image
	    imcopy (in[i]//"[SCI]", tmpinimage, verbose-, >& "dev$null")
	} else {
            # Fix the data if number of non-destructive reads > 1
	    if (l_lnrs != 1) {
                printlog ("Dividing "// in[i] //" by number of non-destructive \
                    reads: "// l_lnrs, l_logfile, l_verbose) 
                imexpr ("a/"// l_lnrs, tmpinimage, in[i] //"[1][*,*,1]",
                    outtype="real", >& "dev$null")
            } else
                # Copy just the first two dimensions of the input image - ignore 
                # third dimension
                imcopy (in[i] //"[1][*,*,1]", tmpinimage, verbose-, >& "dev$null")
	}

        #------------------------------------------------------------------
        # Determine whether the input image is longslit or mos
        keypar (in[i] //"[0]", l_key_fpmask, silent+, >& "dev$null")
        if (keypar.found)
            l_mospos = str(keypar.value)

        longslit = no
        mos = no
        if (substr (l_mospos, 1, (strlen (l_mospos) - 1)) == "mos") {
            mos = yes
        }
        if (substr (l_mospos, strlen (l_mospos) - 3, \
            strlen (l_mospos)) == "slit") {
            longslit = yes
        }

        #------------------------------------------------------------------
        # If adding the MDF, test if the image is of a valid format to receive 
        # an MDF (only MOS data should have MDFs attached).

        if (mos) {
            if (l_fl_addmdf) {
                addmdf[i] = yes
            } else {
                # Warn if not adding MDF to a MOS observation
                printlog ("WARNING - F2PREPARE: Spectroscopic frame "// \
                    in[i] //" will have no MDF, use fl_addmdf=yes", \
                    l_logfile, verbose+)
                addmdf[i] = no
            }
        }

        if (l_fl_addmdf && addmdf[i]) {
            # Check if the input is a dark
            keypar (in[i] //"[0]", "FILTER2", silent+, >& "dev$null")
            if (keypar.found)
                l_filter2 = str(keypar.value)
            else {
                # Try FILT2POS for backwards compatibility with engineering
                # data 
                keypar (in[i] //"[0]", "FILT2POS", silent+, >& "dev$null")
                if (keypar.found)
                    l_filter2 = str(keypar.value)
                else {
                    printlog ("ERROR - F2PREPARE: Could not read the filter \
                        from the header of file", l_logfile, verbose+)
                    printlog ("                   "// in[i], \
                        l_logfile, verbose+)
                    status = 1
                    goto clean
                }
            }

            if (substr (l_filter2, strlen (l_filter2) - 5, \
                strlen (l_filter2) - 4) == "_G")
                l_filter2 = substr (l_filter2, 1, (strlen (l_filter2) - 6))

            if (l_filter2 == "DK" || l_filter2 == "Dark") {
                if (l_fl_dark_mdf) {
                    # Add an MDF to the dark
                    printlog ("WARNING - F2PREPARE: " // in[i] // " is a \
                        dark, but fl_dark_mdf=yes,", l_logfile, verbose+)
                    printlog ("                     so adding MDF", \
                        l_logfile, verbose+)
                    addmdf[i] = yes
                } else {
                    # Don't add an MDF to the dark
                    printlog ("WARNING - F2PREPARE: " // in[i]// " is a dark, \
                        so no MDF will be added", l_logfile, \
                        verbose+)
                    addmdf[i] = no
                }
            }
        }

        # Test if the MDF keyword exists in the header, check if it is the 
        # same as the input file. 

        if (l_fl_addmdf && addmdf[i]) {

            keypar (in[i] //"[0]", l_key_mdf, silent+, >& "dev$null")

            if (keypar.found) {
                # MDF in the header
                mdfname[i] = str(keypar.value)
                printlog ("F2PREPARE: Using MDF defined in the header "//\
                    mdfname[i], l_logfile, l_verbose)

            } else {
                printlog ("F2PREPARE: Using MDF defined in the parameter \
                    list "// l_mdffile, l_logfile, l_verbose)
                    mdfname[i] = l_mdffile
            }

            len = strlen (mdfname[i])
            suf = substr (mdfname[i], len - 4, len)
            if (suf != ".fits")
                mdfname[i] = mdfname[i] //".fits"
            if (access (mdfname[i]))
                printlog ("F2PREPARE: Taking MDF from the current directory.",
                    l_logfile, l_verbose)
            else {
                printlog ("F2PREPARE: Taking MDF from directory "// mdfdir,
                    l_logfile, l_verbose)
                mdfname[i] = l_mdfdir // mdfname[i]
            }

            if ((mdfname[i] == "") || (stridx (" ", mdfname[i]) > 0)) {
                printlog ("ERROR - F2PREPARE: The MDF filename is not \
                    defined.", l_logfile, verbose+)
                status = 1
                goto clean
            }

            if (access (mdfname[i]) == no) {
                printlog ("ERROR - F2PREPARE: The MDF file does not exist.",
                    l_logfile, verbose+)
                status = 1
                goto clean
            }

        } # end of fl_addmdf 

        #------------------------------------------------------------------
        # Correct for non-linearity
        if (l_fl_correct && !alreadyfixed[i]) {
            x = l_sat / 32767.0
            x = (((l_coeff3 * x) + l_coeff2) * x) + l_coeff1
            l_sat = l_sat * x
            x = l_linlimit / 32767.0
            x = (((l_coeff3 * x) + l_coeff2) * x) + l_coeff1
            l_linlimit = int (l_linlimit * x)
            irlincor (tmpinimage, tmpinimage, coeff1=l_coeff1, \
                coeff2=l_coeff2, coeff3=l_coeff3)
        }

        # Since there is currently no linearity correction for F2, check that
        # no more than 5% of the pixels are non-linear
        imstat (tmpinimage, fields="npix", lower=l_linlimit, upper=INDEF, \
            nclip=0, lsigma=3., usigma=3., binwidth=0.1, format-, cache-) | \
            scan (linpix)
        if ((linpix / 4194304) * 100 > 5) {
            printlog ("WARNING - F2PREPARE: More than 5% of the pixels are \
                non-linear", l_logfile, verbose+)
        }

        #------------------------------------------------------------------
        if (l_fl_vardq) {

            # Create the variance frame, if it doesn't already exist
            # The variance frame is generated as:
            #     var = (readnoise / gain)^2 + max (data, 0.0) / gain

            if (!imaccess (in[i] //"["// l_var_ext //"]")) {

                l_varexpression = \
                    "((max (a, 0.0))/"// l_gain //" + \
                    ("// l_ron //"/"// l_gain //")**2)"
                imexpr (l_varexpression, tmpvar, tmpinimage, outtype="real",
                    verbose-)

            } else {
                printlog ("WARNING - F2PREPARE: Variance frame already exists \
                    for "// in[i] //".", l_logfile, verbose+) 
                printlog ("                    New variance frame not \
                    created.", l_logfile, verbose+)
                imcopy (in[i] //"["// l_var_ext //"]", tmpvar, verbose-)
            }

            #--------------------------------------------------------------
            # Create the DQ frame, if it doesn't already exist
            # The preliminary DQ frame is constructed by using the bad pixel
            # mask to set bad pixels to 1, pixels in the non-linear regime 
            # to 2 and saturated pixels to 4.

            if (!imaccess (in[i] //"["// l_dq_ext //"]")) {
                if (l_fl_saturated == no && l_fl_nonlinear == no) {
                    l_tempexpression = "0 * a"
                } else if (l_fl_saturated && l_fl_nonlinear) {
                    l_tempexpression = \
                        "(a > "// l_sat //") ? 4 : \
                        ((a > "// l_linlimit //") ? 2 : 0)"
                } else if (l_fl_saturated) {
                    l_tempexpression = "(a > "// l_sat //") ? 4 : 0"
                } else if (l_fl_nonlinear) {
                    l_tempexpression = "(a > "// linlimit //") ? 2 : 0"
                }
                imexpr (l_tempexpression, tmpsci, tmpinimage, outtype="short",
                    verbose-)

                # If there's no BPM, then just keep the saturated pixels
                if (l_bpm == "none") {
                    addmasks (tmpsci //".fits", tmpdq //".fits", "im1")
                } else {
                    addmasks (tmpsci //".fits, "// l_bpm, tmpdq //".fits",
                        "im1 || im2")
                    update_keyword = yes
                }
                imdelete (tmpsci, verify-, >& "dev$null")
            } else {
                printlog ("WARNING - F2PREPARE: Data quality frame already \
                    exists for "// in[i] //".", l_logfile, verbose+) 
                printlog ("                    New DQ frame not created.",
                    l_logfile, verbose+)
                imcopy (in[i] //"["// l_dq_ext //"]", tmpdq //".fits", \
                    verbose-)
            }

            #--------------------------------------------------------------
            # Pack up the results and clean up

            wmef (tmpinimage //","// tmpvar //".fits,"// tmpdq //".fits", \
                out[i], extnames=l_sci_ext //","// l_var_ext //","// \
                l_dq_ext, phu=in[i] //".fits[0]", verbose-, >& "dev$null")
            if (wmef.status != 0) {
                printlog ("ERROR - F2PREPARE: Could not write final MEF file \
                    (WMEF).", l_logfile, verbose+)
                status = 1
                goto clean
            }
            if (update_keyword)
                gemhedit (out[i] //"[0]", "BPMFILE", l_bpm, \
                    "Input bad pixel mask file", delete-)

            gemhedit (out[i] //"["// l_sci_ext //"]", "EXTVER", 1, "", delete-)
            gemhedit (out[i] //"["// l_var_ext //"]", "EXTVER", 1, "", delete-)
            gemhedit (out[i] //"["// l_dq_ext //"]", "EXTVER", 1, "", delete-)
            imdelete (tmpvar //","// tmpdq, verify-)

            gemhedit (out[i] //"["// l_sci_ext //"]", "AXISLAB3", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_sci_ext //"]", "WCSDIM", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_sci_ext //"]", "LTM3_3", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_sci_ext //"]", "WAXMAP01", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_sci_ext //"]", "WAT3_001", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_var_ext //"]", "AXISLAB3", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_var_ext //"]", "WCSDIM", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_var_ext //"]", "LTM3_3", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_var_ext //"]", "WAXMAP01", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_var_ext //"]", "WAT3_001", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_dq_ext //"]", "AXISLAB3", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_dq_ext //"]", "WCSDIM", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_dq_ext //"]", "LTM3_3", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_dq_ext //"]", "WAXMAP01", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_dq_ext //"]", "WAT3_001", "", "", \
                delete+)

            # end if (l_fl_vardq) 
        } else {
            wmef (tmpinimage, out[i], extnames=l_sci_ext, \
                phu=in[i] //".fits[0]", verbose-, >& "dev$null")
            if (wmef.status != 0) {
                printlog ("ERROR - F2PREPARE: Could not write final MEF file \
                    (WMEF).", l_logfile, verbose+)
                status = 1
                goto clean
            }
            gemhedit (out[i] //"["// l_sci_ext //"]", "EXTVER", 1, "", delete-)
            gemhedit (out[i] //"["// l_sci_ext //"]", "AXISLAB3", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_sci_ext //"]", "WCSDIM", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_sci_ext //"]", "LTM3_3", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_sci_ext //"]", "WAXMAP01", "", "", \
                delete+)
            gemhedit (out[i] //"["// l_sci_ext //"]", "WAT3_001", "", "", \
                delete+)
        }
        
        if (l_fl_addmdf && addmdf[i]) {
            fxinsert (mdfname[i], out[i] //".fits[0]", groups="1", verbose-, \
                >& "dev$null")
        }

        #------------------------------------------------------------------
        # Fix up the headers
        gemhedit (out[i] //"[0]", "WMEF", "", "", delete+)
        gemhedit (out[i] //"[0]", "GEM-TLM", "", "", delete+)

        # Add comment to NEXTEND (added by WMEF)
        keypar (out[i] //"[0]", "NEXTEND", silent+, >& "dev$null")
        if (keypar.found) {
            gemhedit (out[i] //"[0]", "NEXTEND", int (keypar.value), \
                "Number of extensions", delete-)
        }

        if (!alreadyfixed[i]) {

            # Copy the original WCS in the pixel data extensions to the PHU, if
            # it exists. This allows the original WCS to be stored, since
            # nsappwave updates the WCS in the pixel data extensions,
            # overwriting the original WCS in F2 data.
            tmpheader = mktemp("tmpheader")
            keypar (in[i] //"[1]", "CD1_1", silent+, >& "dev$null")
            if (keypar.found) {
                printf ("%-8s= %20s / %-s\n", "CD1_1", keypar.value, \
                    "WCS matrix element (1,1) in (deg/pix)", >> tmpheader)
            }
            keypar (in[i] //"[1]", "CD1_2", silent+, >& "dev$null")
            if (keypar.found) {
                printf ("%-8s= %20s / %-s\n", "CD1_2", keypar.value, \
                    "WCS matrix element (1,2) in (deg/pix)", >> tmpheader)
            }
            keypar (in[i] //"[1]", "CD2_1", silent+, >& "dev$null")
            if (keypar.found) {
                printf ("%-8s= %20s / %-s\n", "CD2_1", keypar.value, \
                    "WCS matrix element (2,1) in (deg/pix)", >> tmpheader)
            }
            keypar (in[i] //"[1]", "CD2_2", silent+, >& "dev$null")
            if (keypar.found) {
                printf ("%-8s= %20s / %-s\n", "CD2_2", keypar.value, \
                    "WCS matrix element (2,2) in (deg/pix)", >> tmpheader)
            }
            keypar (in[i] //"[1]", "CRPIX1", silent+, >& "dev$null")
            if (keypar.found) {
                printf ("%-8s= %20s / %-s\n", "CRPIX1", keypar.value, \
                    "Pixel at which reference applies", >> tmpheader)
            }
            keypar (in[i] //"[1]", "CRPIX2", silent+, >& "dev$null")
            if (keypar.found) {
                printf ("%-8s= %20s / %-s\n", "CRPIX2", keypar.value, \
                    "Pixel at which reference applies", >> tmpheader)
            }
            keypar (in[i] //"[1]", "CRVAL1", silent+, >& "dev$null")
            if (keypar.found) {
                printf ("%-8s= %20s / %-s\n", "CRVAL1", keypar.value, \
                    "Ref pixel value in degrees", >> tmpheader)
            }
            keypar (in[i] //"[1]", "CRVAL2", silent+, >& "dev$null")
            if (keypar.found) {
                printf ("%-8s= %20s / %-s\n", "CRVAL2", keypar.value, \
                    "Ref pixel value in degrees", >> tmpheader)
            }
            keypar (in[i] //"[1]", "CTYPE1", silent+, >& "dev$null")
            if (keypar.found) {
                printf ("%-8s= %-20s / %-s\n", "CTYPE1", \
                    "'"// str(keypar.value) //"'", "RA---TAN X-axis", \
                    >> tmpheader)
            }
            keypar (in[i] //"[1]", "CTYPE2", silent+, >& "dev$null")
            if (keypar.found) {
                printf ("%-8s= %-20s / %-s\n", "CTYPE2", \
                    "'"// str(keypar.value) //"'", "DEC--TAN Y-axis", \
                    >> tmpheader)
            }
            keypar (in[i] //"[1]", "RADECSYS", silent+, >& "dev$null")
            if (keypar.found) {
                printf ("%-8s= %-20s / %-s\n", "RADECSYS", \
                    "'"// str(keypar.value) //"'", "R.A./DEC. coordinate \
                    system reference", >> tmpheader)
            }
            if (access(tmpheader)) {
                mkheader (out[i] //"[0]", tmpheader, append+, verbose-)
            }
            delete (tmpheader, verify-, >& "dev$null")

            update_keyword = no
            keypar (in[i] //"[0]", l_key_ron, silent+, >& "dev$null")
            if (keypar.found) {
                if ((real (keypar.value) == 0.) || \
                    (real (keypar.value) == 0.0)) {
                    # If the keyword value is equal to 0, update the keyword 
                    # value
                    update_keyword = yes

                } else if ((real (keypar.value) != 0.) && \
                    (real (keypar.value) != 0.0) && \
                    (real (keypar.value) != l_ron)) {
                    # If the keyword value is different to the calculated value
                    # (but not equal to 0), copy the keyword value to RAWRDNOI
                    # and update the keyword with the calculated value 
                    gemhedit (out[i] //"[0]", "RAWRDNOI", \
                        real (keypar.value), "Raw "// l_key_ron //" keyword \
                        value", delete-)
                    update_keyword = yes
                }
            } else {
                # If the keyword doesn't exist, update the keyword with the
                # calculated value
                update_keyword = yes
            }

            if (update_keyword) {
                gemhedit (out[i] //"[0]", l_key_ron, l_ron,
                    "Estimated read noise (electrons)", delete-)
            }

            update_keyword = no
            keypar (in[i] //"[0]", l_key_gain, silent+, >& "dev$null")
            if (keypar.found) {
                if ((real (keypar.value) == 0.) || \
                    (real (keypar.value) == 0.0)) {
                    # If the keyword value is equal to 0, update the keyword 
                    # value
                    update_keyword = yes

                } else if ((real (keypar.value) != 0.) && \
                    (real (keypar.value) != 0.0) && \
                    (real (keypar.value) != l_gain)) {
                    # If the keyword value is different to the calculated value
                    # (but not equal to 0), copy the keyword value to RAWGAIN
                    # and update the keyword with the calculated value 
                    gemhedit (out[i] //"[0]", "RAWGAIN", real (keypar.value),
                        "Raw "// l_key_gain //" keyword value", delete-)
                    update_keyword = yes
                }
            } else {
                # If the keyword doesn't exist, update the keyword with the
                # calculated value
                update_keyword = yes
            }
            if (update_keyword) {
                gemhedit (out[i] //"[0]", l_key_gain, l_gain,
                    "Gain (electrons/ADU)", delete-)
            }

            gemhedit (out[i] //"[0]", l_key_sat, l_sat, \
                "Saturation level (ADU)", delete-)
            gemhedit (out[i] //"[0]", l_key_nonlinear, l_linlimit,
                "Non-linear regime (ADU)", delete-)

            # Parse the individual filter entries into a single header keyword
            # Determine the value for the FILTER1 keyword
            keypar (in[i] //"[0]", "FILTER1", silent+, >& "dev$null")
            if (keypar.found)
                l_filter1 = str(keypar.value)
            else {
                # Try FILT1POS for backwards compatibility with engineering
                # data 
                keypar (in[i] //"[0]", "FILT1POS", silent+, >& "dev$null")
                if (keypar.found)
                    l_filter1 = str(keypar.value)
                else {
                    printlog ("ERROR - F2PREPARE: Could not read the filter \
                        from the header of file", l_logfile, verbose+)
                    printlog ("                   "//in[i], \
                        l_logfile, verbose+)
                    status = 1
                    goto clean
                }
            }

            # Determine the value for the FILTER2 keyword
            keypar (in[i] //"[0]", "FILTER2", silent+, >& "dev$null")
            if (keypar.found)
                l_filter2 = str(keypar.value)
            else {
                # Try FILT2POS for backwards compatibility with engineering
                # data 
                keypar (in[i] //"[0]", "FILT2POS", silent+, >& "dev$null")
                if (keypar.found)
                    l_filter2 = str(keypar.value)
                else {
                    printlog ("ERROR - F2PREPARE: Could not read the filter \
                        from the header of file", l_logfile, verbose+)
                    printlog ("                   "// in[i], \
                        l_logfile, verbose+)
                    status = 1
                    goto clean
                }
            }

            # First rename the HK_G0806_good filter to HK_G0806
            if (l_filter1 == "HK_G0806_good")
                l_filter1 = "HK_G0806"
            if (l_filter2 == "HK_G0806_good")
                l_filter2 = "HK_G0806"

            # Reject "Open"
            if (l_filter1 == "0" || l_filter1 == "OPEN" || \
                l_filter1 == "open" || l_filter1 == "Open" || \
                l_filter1 == "INVALID")
                l_filter1 = ""
            if (l_filter2 == "0" || l_filter2 == "OPEN" || \
                l_filter2 == "open" || l_filter2 == "Open" || \
                l_filter2 == "INVALID")
                l_filter2 = ""

            # Strip the Gemini filter numbers before making the combined 
            # filter keyword
            if (substr (l_filter1, strlen (l_filter1) - 5, \
                strlen (l_filter1) - 4) == "_G")
                l_filter1 = substr (l_filter1, 1, (strlen (l_filter1) - 6))
            if (substr (l_filter2, strlen (l_filter2) - 5, \
                strlen (l_filter2) - 4) == "_G")
                l_filter2 = substr (l_filter2, 1, (strlen (l_filter2) - 6))

            # If FILTER2 is the dark filter, set FILTER equal to "Dark"
            if (l_filter2 == "DK" || l_filter2 == "Dark") {
                l_filter1 = ""
                l_filter2 = "Dark"
            }

            printf ("%s\n%s\n", l_filter1, l_filter2) | \
                sort ("STDIN", col=0, ignore+, num-, rev-, > tmpfile2)

            l_filter = ""
            l_temp = ""
            scanfile = tmpfile2
            while (fscan (scanfile, l_temp) != EOF)
                l_filter = l_filter + l_temp

            if (l_filter == "")
                l_filter = "Open"

            update_keyword = no
            keypar (in[i] //"[0]", l_key_filter, silent+, >& "dev$null")
            if (keypar.found) {
                if (keypar.value != l_filter) {
                    # If the keyword value is different to the calculated 
                    # value, copy the keyword value to RAWFILT and update the
                    # keyword with the calculated value
                    gemhedit (out[i] //"[0]", "RAWFILT", str (keypar.value),
                        "Raw "// l_key_filter //" keyword value", delete-)
                    update_keyword = yes
                }
            } else {
                # If the keyword doesn't exist, update the keyword with the
                # calculated value
                update_keyword = yes
            }
            if (update_keyword) 
                gemhedit (out[i] //"[0]", l_key_filter, l_filter,
                    "Filter name combined from both wheels", delete-)

            scanfile=""
            delete(tmpfile2, verify-, >& "dev$null")

            if (l_fl_correct) {
                printf ("%7.5f %7.5f %7.5f\n", l_coeff1, l_coeff2, l_coeff3) \
                    | scan (l_struct)
                gemhedit (out[i] //"[0]", "NONLINCR", l_struct,
                    "Non-linear correction applied", delete-)
            }

            #--------------------------------------------------------------
            # Compute pixscale and add to header

            keypar (in[i] //"[1]", "CD1_1", silent+, >& "dev$null")
            if (keypar.found) {
                l_cd11 = real (keypar.value)
            } else {
                # Look in the PHU for backwards compatibility with engineering
                # data 
                keypar (in[i] //"[0]", "CD1_1", silent+, >& "dev$null")
                if (keypar.found) {
                    l_cd11 = real (keypar.value)
                } else {
                    printlog ("ERROR - F2PREPARE: Could not read CD1_1 from \
                        the header of file", l_logfile, verbose+)
                    printlog ("                   "//in[i], l_logfile, \
                        verbose+)
                    status = 1
                    goto clean
                }
            }

            keypar (in[i] //"[1]", "CD1_2", silent+, >& "dev$null")
            if (keypar.found) {
                l_cd12 = real (keypar.value)
            } else {
                # Look in the PHU for backwards compatibility with engineering
                # data 
                keypar (in[i] //"[0]", "CD1_2", silent+, >& "dev$null")
                if (keypar.found) {
                    l_cd12 = real (keypar.value)
                } else {
                    printlog ("ERROR - F2PREPARE: Could not read CD1_2 from \
                        the header of file", l_logfile, verbose+)
                    printlog ("                   "//in[i], l_logfile, \
                        verbose+)
                    status = 1
                    goto clean
                }
            }

            keypar (in[i] //"[1]", "CD2_1", silent+, >& "dev$null")
            if (keypar.found) {
                l_cd21 = real (keypar.value)
            } else {
                # Look in the PHU for backwards compatibility with engineering
                # data 
                keypar (in[i] //"[0]", "CD2_1", silent+, >& "dev$null")
                if (keypar.found) {
                    l_cd21 = real (keypar.value)
                } else {
                    printlog ("ERROR - F2PREPARE: Could not read CD2_1 from \
                        the header of file", l_logfile, verbose+)
                    printlog ("                   "//in[i], l_logfile, \
                        verbose+)
                    status = 1
                    goto clean
                }
            }

            keypar (in[i] //"[1]", "CD2_2", silent+, >& "dev$null")
            if (keypar.found) {
                l_cd22 = real (keypar.value)
            } else {
                # Look in the PHU for backwards compatibility with engineering
                # data 
                keypar (in[i] //"[0]", "CD2_2", silent+, >& "dev$null")
                if (keypar.found) {
                    l_cd22 = real (keypar.value)
                } else {
                    printlog ("ERROR - F2PREPARE: Could not read CD2_2 from \
                        the header of file", l_logfile, verbose+)
                    printlog ("                   "//in[i], l_logfile, \
                        verbose+)
                    status = 1
                    goto clean
                }
            }

            # Calculate the pixel scale
            l_pixscale = 3600. * (sqrt (l_cd11*l_cd11 + l_cd12*l_cd12) + \
                sqrt (l_cd21*l_cd21 + l_cd22*l_cd22)) / 2.
            l_pixscale = real (int (l_pixscale*10000.)/10000.)

            if (l_pixscale != 0.) {

                keypar (in[i] //"[0]", l_key_pixscale, silent+, >& "dev$null")
                if (keypar.found) {
                    gemhedit (out[i] //"[0]", "RAWPIXSC", real (keypar.value),
                        "Raw "// l_key_pixscale //" keyword value", delete-)
                }

                gemhedit (out[i] //"[0]", l_key_pixscale, l_pixscale,
                    "Pixel scale in arcsec/pixel", delete-)
            }

            #--------------------------------------------------------------
            # Set DISPAXIS for both longslit and mos data. The DISAXIS keyword
            # is written to the PHU.
            if (longslit || mos) {
                gemhedit (out[i] //"[0]", l_key_dispaxis, 2, \
                    "Dispersion axis", delete-)
            }

            # Fix offsets in header
            keypar (in[i] //"[0]", l_key_xoff, silent+, >& "dev$null")
            if (keypar.found)
                l_xoff = real (keypar.value)
            keypar (in[i] //"[0]", l_key_yoff, silent+, >& "dev$null")
            if (keypar.found)
                l_yoff = real (keypar.value)

            gemhedit (out[i] //"[0]", l_key_xoff, -l_yoff, "", delete-)
            gemhedit (out[i] //"[0]", l_key_yoff, -l_xoff, "", delete-)

            # Fix date keyword in header
            keypar (in[i] //"[0]", l_key_date, silent+, >& "dev$null")
            if (!keypar.found) {
                keypar (in[i] //"[0]", "DATE", silent+, >& "dev$null")
                if (keypar.found) {
                    l_dateobs = keypar.value
                    gemhedit (out[i] //"[0]", l_key_date, l_dateobs, 
                        "UT Date of observation (YYYY-MM-DD)", delete-)
                }
            }

            # Time stamps
            gemdate ()
            gemhedit (out[i] //"[0]", "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI", delete-)
            gemhedit (out[i] //"[0]", "PREPARE", gemdate.outdate,
                "UT Time stamp for F2PREPARE", delete-)

            # Print output to logfile
            printf ("%3.0d %15s --> %16s \n", i, in[i], out[i]) | \
                scan (l_struct)
            printlog (l_struct, l_logfile, l_verbose)
            printf ("      %17s %15s %15s %5.1f %5.1f %7.0d \n", l_filter,
                l_mospos, l_bpm, l_ron, l_gain, l_sat) | scan (l_struct)
            printlog (l_struct, l_logfile, l_verbose)

        } # end if (!alreadyfixed[i])

        # Delete tmp files created in this loop
        imdelete (tmpinimage //","// tmpsci //","// tmpvar //","// tmpdq, \
            verify-, >& "dev$null")
        
        i = i + 1
    } # end loop


clean:
    #----------------------------------------------------------------------
    # Clean up
    if (status == 1) {
        if (noutimages > 0) {
            for (i = 0; i <= noutimages; i += 1) {
                imdelete (out[i], verify-, >& "dev$null")
            }
        }
        printlog (" ", l_logfile, verbose+)
        printlog ("F2PREPARE exit status: fail.", l_logfile, verbose+)
    }

    if (status == 0) {
        printlog (" ", l_logfile, verbose+)
        printlog ("F2PREPARE exit status: good.", l_logfile, verbose+)
    }

    printlog ("--------------------------------------------------------------\
        --------------", l_logfile, l_verbose)

    scanfile = "" 
    delete (tmpfile // "," // tmpfile2 // "," // tmpchk1 // "," // tmpchk2, \
        verify-, >& "dev$null")
    tdelete (tmparray, go_ahead+, verify-, >& "dev$null")
    imdelete (tmpinimage // "," // tmpsci // "," // tmpvar // "," // tmpdq, \
        verify-, >& "dev$null")

end
