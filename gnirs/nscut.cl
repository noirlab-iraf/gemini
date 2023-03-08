# Copyright(c) 2001-2012 Association of Universities for Research in Astronomy, Inc.

procedure nscut(inimages)

# Version Sept 20, 2002 JJ v1.4 release
#         Aug 19, 2003  KL IRAF2.12 - new parameter, addonly, in hedit
#         Oct 29, 2003  KL moved from niri to gnirs
#         Nov 5, 2003   AC initial reorganisation (as earlier work)
#                       - reformatting
#                       - replaced "default" for keysection with default value
#                         (can be reset using unlearn)
#                       - replaced "default" for outprefix with fl_autoprefix
#                       - only prompt for variables required
#                       - move temp images to tmp$
#                       - clarified key_section/section logic
#                         (BUG? - if only l_section defined, was not used)
#         Nov 5, 2003   AC removed calls to gemoffset and nsappwave
#         Nov 5, 2003   AC added gemextn calls (requires logging integration)
#         Nov 6, 2003   AC added (untested) support for sections in MDF extn
#         Nov 7, 2003   AC further support for extensions
#         Nov 8, 2003   AC removed index selection for input - assumes [SCI,1]

char    inimages    {prompt = "Input files"}                        # OLDP-1-primary-single-prefix=s
char    outspectra  {"", prompt = "Output spectra list"}            # OLDP-1-output
char    outprefix   {"s", prompt = "Prefix for output spectra"}     # OLDP-4
char    section     {"", prompt = "Alternative section or keyword (blank for MDF)"} # OLDP-3
bool    fl_corner   {yes, prompt = "Blank corners (if info in MDF)?"}   # OLDP-3
char    logfile     {"", prompt = "Logfile"}                        # OLDP-1
bool    verbose     {yes, prompt = "Verbose output?"}               # OLDP-2
bool    debug       {no, prompt = "Very verbose output?"}           # OLDP-2
int     status      {0, prompt = "Exit status (0=good)"}            # OLDP-4

struct  *scanin1    {"", prompt="Internal use only"}                # OLDP-4
struct  *scanin2    {"", prompt="Internal use only"}                # OLDP-4

begin

    char    l_inimages, l_outspectra, l_outprefix
    char    l_section = ""
    bool    l_fl_corner
    char    l_key_pixscale = ""
    char    l_key_order = ""
    char    l_key_waveorder = ""
    char    l_key_cut_section = ""
    char    l_key_instrument = ""
    char    l_key_dark = ""
    char    l_val_dark = ""
    char    l_sci_ext, l_var_ext, l_dq_ext
    char    l_logfile = ""
    bool    l_verbose
    bool    l_debug

    char    tmproot, tmprootuncut, tmpall, tmpout, tmpsci, tmpvar, tmpdq

    char    secn, inimg, outimg, imdest, imsrc, subsec, outver
    char    phu, badhdr, line, keyfound, maskexpr, unused
    char    instrument, valdark
    int     nx, ny, sx, sy, nsecns, isecn, junk, x1, x2, y1, y2, specorder
    int     nfiles, corner, nbad, dx, xover
    bool    intdbg, loopkeyw, first, notfirst
    struct  sdate, sline

    real    pixscale, slitwidth
    int     slitid, undefined_order

    bool    secn_mdf, secn_text, secn_keyw

    cache ("keypar", "gemextn", "tinfo", "nsmdfhelper")
    status = 1 # Error on early exit
    intdbg = no
    undefined_order = -9999 # -1 and larger values are ok
    unused = "(unused)"

    secn_mdf = no
    secn_text = no
    secn_keyw = no

    l_inimages = inimages; l_outspectra = outspectra; l_outprefix = outprefix
    l_verbose = verbose; l_fl_corner = fl_corner; l_debug = debug

    junk = fscan (section, l_section)
    junk = fscan (logfile, l_logfile) 

    # Shared definitions
    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_instrument, l_key_instrument)
    if ("" == l_key_instrument) badhdr = badhdr + " key_instrument"
    junk = fscan (nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"
    junk = fscan (nsheaders.key_order, l_key_order)
    if ("" == l_key_order) badhdr = badhdr + " key_order"
    junk = fscan (nsheaders.key_waveorder, l_key_waveorder)
    if ("" == l_key_waveorder) badhdr = badhdr + " key_waveorder"
    junk = fscan (nsheaders.key_cut_section, l_key_cut_section)
    if ("" == l_key_cut_section) badhdr = badhdr + " key_cut_section"
    junk = fscan (nsheaders.key_dark, l_key_dark)
    if ("" == l_key_dark) badhdr = badhdr + " key_dark"
    junk = fscan (nsheaders.val_dark, l_val_dark)
    if ("" == l_val_dark) badhdr = badhdr + " val_dark"

    # Make temp files:

    tmproot = mktemp ("tmproot") 
    tmprootuncut = mktemp ("tmprootuncut") 
    tmpout = mktemp ("tmpout")
    tmpall = mktemp ("tmpall")

    # Start logging to file

    if (l_logfile == "") {
        junk = fscan (gnirs.logfile, l_logfile) 
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSCUT: Both nscut.logfile and \
                gnirs.logfile are", l_logfile, verbose+) 
            printlog ("                 undefined. Using " \
                // l_logfile, l_logfile, verbose+) 
        }
    }

    date | scan (sdate) 
    printlog ("----------------------------------------------------\
        --------------------------", l_logfile, l_verbose) 
    printlog ("NSCUT -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NSCUT: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }

    # Choose source for section information
    if (intdbg) print ("choose source: " // l_section)

    # Decide what the source of section information should be
    if (l_section == "") {
        # Use MDF
        secn_mdf = yes
    } else if (stridx ("[", l_section) == 1) {
        # Use explicit section "[...]"
        secn_text = yes
    } else {
        # Use keyword (direct, or with index appended)
        secn_keyw = yes
    }

    # Quick check that all input files exist and have a [SCI,1] extension and,
    # if an MDF will be used for cutting, an [MDF] extension
    if (intdbg) print ("exist check")

    # TODO - integrate gemextn errfile with logging
    gemextn (l_inimages, check="exists,mef,image", process="append", \
        index="", extname=l_sci_ext, extversion="1", ikparams="", \
        omit="", replace="", outfile="dev$null", logfile="", glogpars="",
        verbose=l_verbose)
    if (gemextn.fail_count != 0 || gemextn.count == 0) {
        printlog ("ERROR - NSCUT: Missing science data", l_logfile, 
            verbose+) 
        goto clean
    }
    if (secn_mdf) {
        gemextn (l_inimages, check="exists,mef,table", process="append", \
            index="", extname="MDF", extversion="", ikparams="", omit="", \
            replace="", outfile="dev$null", logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - NSCUT: Missing MDF data", \
                l_logfile, verbose+) 
            goto clean
        }
    }

    # Now expand without kernel + extensions to get root names
    gemextn (l_inimages, check="", process="append", index="", 
        extname=l_sci_ext, extversion="1", ikparams="", omit="kernel,exten",
        replace="", outfile=tmproot, logfile="", glogpars="",
        verbose=l_verbose)

    if (intdbg) print ("nscut check")

    # Filter files already processed and check for PREPARE
    # TODO - Requires utility to check headers
    nbad = 0
    scanin1 = tmproot
    while (fscan (scanin1, inimg) != EOF) {
        keypar (inimg // "[0]", "NSCUT", silent+) 
        if (keypar.found) {
            printlog ("WARNING - NSCUT: File " // inimg \
                // " already run through NSCUT.", l_logfile, verbose+) 
        } else {
            print (inimg, >> tmprootuncut) 
        }

        keyfound = ""
        hselect(inimg // "[0]", "*PREPAR*", yes) | scan(keyfound)
        if (keyfound == "") {
            printlog ("ERROR - NSCUT: Image " // inimg \
                // " not PREPAREd.", l_logfile, verbose+)
            nbad += 1
        }
    }
    if (nbad > 0) {
        printlog ("ERROR - NSCUT: " // nbad // " image(s) \
            have not been run through *PREPARE", l_logfile, verbose+) 
        goto clean
    }
    if (!access (tmprootuncut)) {
        printlog ("ERROR - NSCUT: no input files to process.", \
            l_logfile, verbose+) 
        goto clean
    }

    # Generate the output images

    if (intdbg) print ("output check")

    gemextn (l_outspectra, check="absent", process="none", index="", \
        extname="", extversion="", ikparams="", omit="kernel,exten", 
        replace="", outfile=tmpout, logfile="", glogpars="",
        verbose=l_verbose)
    if (gemextn.fail_count != 0) {
        printlog ("ERROR - NSCUT: Existing or incorrectly formatted \
            output files", l_logfile, verbose+) 
        goto clean
    }

    # If tmpout is empty, the output files names should be 
    # created with the prefix parameter
    # (we could have a separate msg for outprefix="", but that will
    # trigger an error in gemextn anyway)
    if (gemextn.count == 0) {

        if (intdbg) print ("output from substitution")

        gemextn ("%^%" // l_outprefix // "%" // "@" // tmprootuncut, \
            check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", omit="kernel,exten", \
            replace="", outfile=tmpout, logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - NSCUT: No or incorrectly formatted \
                output files", l_logfile, verbose+) 
            goto clean
        }
    }

    # Rather than repeat processing, we calculate the sections here,
    # write info to tmpall, and then process once all done (the
    # idea is to detect errors before processing).

    # For each file in turn
    scanin1 = tmprootuncut
    scanin2 = tmpout
    nfiles = 0
    while (fscan (scanin1, inimg) != EOF) {

        nfiles = nfiles + 1
        if (intdbg) print ("inimg: " // inimg)

        if (fscan (scanin2, outimg) == EOF) {
            printlog ("ERROR - NSCUT: No output image for " \
                // inimg // ".", l_logfile, verbose+) 
            goto clean
        }

        # Read image size
        imsrc = inimg // "[" // l_sci_ext // ",1]"
        keypar (imsrc, "i_naxis1", silent+)
        if (no == keypar.found) {
            printlog ("ERROR - NSCUT: No i_naxis1 in " // imsrc // ".", \
                l_logfile, verbose+)
            goto clean
        } else {
            nx = int (keypar.value)
        }
        keypar (imsrc, "i_naxis2", silent+)
        if (no == keypar.found) {
            printlog ("ERROR - NSCUT: No i_naxis2 in " // imsrc // ".", \
                l_logfile, verbose+)
            goto clean
        } else {
            ny = int (keypar.value)
        }

        # Find the number of sections
        nsecns = 1
        phu = inimg // "[0]"
        if (secn_mdf) {
            if (intdbg) print ("number of sections via mdf")
            tinfo (inimg // ".fits[MDF]", ttout-, >& "dev$null")
            if (tinfo.tbltype != "fits") {
                printlog ("ERROR - NSCUT: Bad MDF in " // inimg // ".",
                    l_logfile, verbose+)
                goto clean
            } else {
                nsecns = tinfo.nrows
            }
        } else if (secn_keyw) {
            if (intdbg) print ("number of sections via " // l_section)
            keypar (phu, l_section, silent+)
            loopkeyw = ! keypar.found
            if (loopkeyw) {
                nsecns = 0
                keypar (phu, l_section // (nsecns+1), silent+)
                while (keypar.found) {
                    nsecns = nsecns + 1
                    keypar (phu, l_section // (nsecns+1), silent+)
                }
                if (0 == nsecns) {
                    # If NIRI dark without a SPECSEC issue explanation
                    # and exit. Otherwise, it is a real error.
                    hselect (phu, l_key_instrument, yes) | scan (instrument)
                    hselect (phu, l_key_dark, yes) | scan (valdark)
                    if ((instrument == "NIRI") && (valdark == l_val_dark)) {
                        printlog ("WARNING - NSCUT: No header value for " \
                            // l_section // " in dark frame " // inimg // ".",
                            l_logfile, verbose+)
                        printlog ("", l_logfile, verbose+)
                        printlog ("                 This is not unusual. It \
                            is recommended *not* to cut the NIRI", 
                            l_logfile, verbose+)
                        printlog ("                 darks before NSFLAT, and \
                            let NSFLAT handle the cutting", l_logfile, \
                            verbose+)
                        printlog ("                 of the darks.",
                            l_logfile, verbose+)
                    } else {
                        printlog ("ERROR - NSCUT: No header value for " \
                            // l_section // " in " // inimg // ".", \
                            l_logfile, verbose+)
                    }
                    goto clean
                }
            }
        }

        # Information needed for MDF calculation
        if (secn_mdf) {
            keypar (phu, l_key_pixscale, silent+)
            if (keypar.found) {
                pixscale = real (keypar.value)
            } else {
                printlog ("ERROR - NSCUT: No header value for " \
                    // l_key_pixscale // " in " // phu // ".", \
                    l_logfile, verbose+)
                goto clean
            }
            if (intdbg) print ("MDF pixscale: " // pixscale)
        }

        # Repeat for each section
        # (possibly more than one if using MDF or header keywords)
        for (isecn = 1; isecn <= nsecns; isecn += 1) {

            if (intdbg) print ("section: " // isecn)
            specorder = isecn
            corner = 0
            xover = 0
            slitwidth = 0
            x2 = 0
            x1 = 0

            # Get the section information
            if (secn_text) {
                secn = l_section
            }

            # For NIRI, this used to be SPECSEC1 by default
            if (secn_keyw) {
                if (loopkeyw) keypar (phu, l_section // isecn, silent+)
                else keypar (phu, l_section, silent+)
                # No need to check - checked above for nsecns
                secn = keypar.value
            }

            if (secn_mdf) {
                
                nsmdfhelper (inimg // ".fits[MDF]", isecn, \
                    inimg // "[" // l_sci_ext // ",1]", \
                    area="spectrum", logfile=l_logfile, \
                    logname="NSCUT", \
                    verbose=((nfiles == 1 && l_verbose) || l_debug))
                if (no == (0 == nsmdfhelper.status)) {
                    printlog ("ERROR - NSCUT: Could not read MDF.", \
                        l_logfile, verbose+)
                    goto clean
                }

                specorder = nsmdfhelper.specorder
                corner = nsmdfhelper.corner
                x1 = nsmdfhelper.ixlo
                y1 = nsmdfhelper.iylo
                x2 = nsmdfhelper.ixhi
                y2 = nsmdfhelper.iyhi
                xover = nsmdfhelper.ixovershoot
                slitwidth = nsmdfhelper.slitwidth

                # Magic value - use header (eg IFU)
                # Note: NIFS does not used key_waveorder
                
                if ((-1 == specorder) && (unused != l_key_waveorder)) {
                    keypar (phu, l_key_waveorder, silent+)
                    if (no == keypar.found) {
                        printlog ("ERROR - NSCUT: No " // l_key_waveorder \
                            // " in " // inimg // ".", \
                            l_logfile, verbose+)
                        goto clean
                    } else {
                        specorder = int (keypar.value)
                    }
                }

                secn = "[" // x1 // ":" // x2 // "," \
                    // y1 // ":" // y2 // "]"

                if (intdbg) print ("scanned from mdf: "//secn)
            }

            # Check that the section lies within the image data

            if (intdbg) print ("size check: " // secn)

            if (stridx (":", secn) > 0 && \
                stridx (":", secn) < stridx (",", secn)) {

                sx = int (substr (secn, stridx (":", secn) +1,
                    stridx (",", secn) -1) )
            } else if (stridx ("*", secn) < stridx (",", secn) \
                && stridx ("*", secn) > 0) {

                sx = nx
            } else {
                printlog ("ERROR - NSCUT: Cannot parse section: " \
                    // secn // ".", l_logfile, verbose+) 
                goto clean
            }

            subsec = substr (secn, stridx (",", secn) +1, strlen (secn) -1)
            if (stridx (":", subsec) > 0) {
                sy = int (substr (subsec, stridx (":", subsec) +1, \
                    strlen (subsec) ) )
            } else if (stridx ("*", subsec) > 0) {
                sy = ny
            } else {
                printlog ("ERROR - NSCUT: Cannot parse section: " \
                    // secn // ".", l_logfile, verbose+) 
                goto clean
            }

            if (nx<sx || ny<sy) {
                printlog ("WARNING - NSCUT: section of " // imsrc \
                    // " is larger than image.", l_logfile, verbose+) 
            }

            # And write to tmpall
            line = inimg // " " // secn // " " // outimg // " " // isecn \
                // " " // specorder // " " // corner // " " \
                // (x2-x1) // " " // slitwidth // " " // xover

            print (line, >> tmpall)
            if (intdbg) print (line)

        } # More sections

    } # More images

    if (no ==access (tmpall) ) {
        printlog ("ERROR - NSCUT: No images will be chopped.", l_logfile,
            verbose+)
        goto clean
    }

    if (intdbg) type (tmprootuncut) 
    if (intdbg) type (tmpout) 

    # Finally, copy the subsection from input to output and update
    # headers etc

    first = yes
    notfirst = no
    scanin1 = tmpall
    printlog ("NSCUT: Processing " // nfiles // " files\n", \
        l_logfile, l_verbose) 
    while (fscan (scanin1, inimg, secn, outimg, isecn, specorder, \
        corner, nx, slitwidth, xover) != EOF) {

        # Make temporary files
        tmpsci = mktemp ("tmpsci")
        tmpvar = mktemp ("tmpvar")
        tmpdq = mktemp ("tmpdq")

        # Won't exist on first section
        if (no == imaccess (outimg)) {

            printf ("%s --> %s\n", inimg, outimg) | scan (sline)
            printlog (sline, l_logfile, l_verbose)

            # allow printing of sections on first pass
            if (notfirst) first = no
            if (first) notfirst = yes

            if (intdbg) print ("copying [0], [MDF] for " // inimg)

            # Copy PHU and MDF
            imcopy (inimg // ".fits[0]", outimg // ".fits", verbose-)
            if (secn_mdf) {
                if (intdbg) print ("copy mdf")
                tcopy (inimg // ".fits[MDF]", outimg // ".fits[MDF]", \
                    verbose-)
            }

            # Set global header info here (extensions below)

            if (intdbg) print ("setting global header info")

            imsrc = inimg // "[" // l_sci_ext // ",1]"
            imdest = outimg // "[0]"
            keypar (imsrc, "i_naxis1", silent+) 
            gemhedit (imdest, "ORIGXSIZ", int (keypar.value), 
                "Original size in X", delete-)
            keypar (imsrc, "i_naxis2", silent+)
            gemhedit (imdest, "ORIGYSIZ", int (keypar.value),
                "Original size in Y", delete-)
            
            gemdate ()
            gemhedit (imdest, "NSCUT", gemdate.outdate, 
                "UT Time stamp for NSCUT", delete-)
            gemhedit (imdest, "GEM-TLM", gemdate.outdate, 
                "UT Last modification with GEMINI IRAF", delete-) 

        }

        if (intdbg) {
            print "gemextn " // specorder
            gemextn (outimg)
        }

        # Destination version
        outver = "," // isecn

        # Copy the section data and update headers (note that
        # nscut used to use NSCUTSEC in [SCI], but Claudia's 
        # script used DETSEC instead, for all extensions. Here we
        # use l_key_cut_section

        if (intdbg) print ("copying data for " // inimg)

        # Use temporary files (tmpsci, tmpvar and tmpdq) instead of the 
        # imdest = outimg // "[" // l_sci_ext // outver // ",append]" syntax to
        # solve a memory stack issue - EH
        imsrc = inimg // "[" // l_sci_ext // ",1]" // secn
        if (first) {
            printlog (" [" // l_sci_ext // ",1]" // secn // " --> " \
                // "[" // l_sci_ext // outver // "]", l_logfile, l_verbose)
        }

        if (l_fl_corner && 0 != corner && (nx + xover - corner) > 0) {
            maskexpr = "(I < " // corner // "*J/" // ny // " || " \
                // "I - " // xover // " > " // slitwidth // "+" // corner // \
                "*J/" // ny // ")"
            if (intdbg) print (maskexpr)
            imexpr (maskexpr // "? 0 : a", tmpsci, imsrc, verbose=intdbg)
        } else {
            if (intdbg) print ("no corner: " // l_fl_corner // ", " \
                // corner // ", " // nx)
            imcopy (imsrc, tmpsci, verbose-)
        }
        gemhedit (tmpsci, l_key_cut_section, secn,
            "Region extracted by NSCUT", delete-)
        gemhedit (tmpsci, l_key_order, str(specorder), "Spectral order",
            delete-)
        imcopy (tmpsci, \
            outimg // "[" // l_sci_ext // outver // ",append]", verbose-)

        imsrc = inimg // "[" // l_var_ext // ",1]" // secn
        if (imaccess (imsrc)) {
            if (intdbg) print ("Working on VAR...")
            if (l_fl_corner && 0 != corner && nx > 0) {
                imexpr (maskexpr // "? 0 : a", tmpvar, imsrc, \
                    verbose=intdbg)
            } else {
                imcopy (imsrc, tmpvar, verbose-)
            }
            gemhedit (tmpvar, l_key_cut_section, secn,
                "Region extracted by NSCUT", delete-)
            gemhedit (tmpvar, l_key_order, str(specorder), "Spectral order",
                delete-)
            imcopy (tmpvar, \
                outimg // "[" // l_var_ext // outver // ",append]", verbose-)
            if (intdbg) print ("Done working on VAR")
        }

        imsrc = inimg // "[" // l_dq_ext // ",1]" // secn
        if (imaccess (imsrc) ) {
            if (intdbg) print ("Working on DQ ...")
            if (l_fl_corner && 0 != corner && nx > 0) {
                imexpr (maskexpr // "? (a | 8) : a", tmpdq, imsrc, \
                    verbose=intdbg)
            } else {
                imcopy (imsrc, tmpdq, verbose-)
            }
            gemhedit (tmpdq, l_key_cut_section, secn,
                "Region extracted by NSCUT", delete-)
            gemhedit (tmpdq, l_key_order, str(specorder), "Spectral order",
                delete-)
            imcopy (tmpdq, \
                outimg // "[" // l_dq_ext // outver // ",append]", verbose-)
            if (intdbg) print ("Done working on DQ.")
        }

        gemhedit (outimg//"[0]", "NSCIEXT", str(isecn), 
            "Number of science extensions", delete-)
        if (intdbg) print ("Done with output: "//outimg)

        imdelete (tmpsci // "," // tmpvar // "," // tmpdq, verify-, \ 
            >& "dev$null") 
    }

    # If we have managed to arrive here, everything has worked correctly
    status = 0

clean:
    scanin1 = ""
    scanin2 = ""
    delete (tmproot // "," // tmprootuncut // "," // tmpall // "," \
        // tmpout, verify-, >& "dev$null")
    imdelete (tmpsci // "," // tmpvar // "," // tmpdq, verify-, >& "dev$null")
    if (status == 0) {
        printlog (" ", l_logfile, l_verbose) 
        printlog ("NSCUT exit status: good.", l_logfile, l_verbose) 
    }
    printlog ("-------------------------------------------------------------\
        -----------------", l_logfile, l_verbose) 

end
