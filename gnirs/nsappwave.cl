# Copyright(c) 2001-2012 Association of Universities for Research in Astronomy, Inc.

procedure nsappwave (inimages) 

# Put an approximate wavelength system in headers of GNIRS/NIRI spectra
#
# Old version log:
#          Sept 20, 2002 JJ v1.4 release
#          Aug 19, 2003  KL IRAF2.12 - new parameter, addonly, in hedit
#          Oct 29, 2003  KL moved from niri to gnirs package

char    inimages    {prompt = "Input GNIRS/NIRI spectra"}           # OLDP-1-primary-single-prefix=w
char    outspectra  {"", prompt = "Output spectra"}                 # OLDP-1-output
char    outprefix   {"w", prompt = "Prefix for output spectra"}     # OLDP-4
char    nsappwavedb {"gnirs$data/nsappwave.fits", prompt = "nsappwave calibration table"}    # OLDP-3

real    crval       {INDEF, prompt = "Central wavelength"}          # OLDP-3
real    cdelt       {INDEF, prompt = "Resolution in wavelength per pixel"}  # OLDP-3
real    cradius     {INDEF, prompt = "Search radius for nswavelength"}      # OLDP-3
bool    fl_phu      {no, prompt = "Force changes to PHU for single science data?"}          # OLDP-3

char    logfile     {"", prompt = "Logfile"}                        # OLDP-1
bool    verbose     {yes, prompt = "Verbose"}                       # OLDP-2
bool    debug       {no, prompt = "Very verbose"}                   # OLDP-3
int     status      {0, prompt = "Exit status (0=good)"}            # OLDP-4

struct  *scanin1    {prompt = "Internal use"}                       # OLDP-4
struct  *scanin2    {prompt = "Internal use"}                       # OLDP-4
struct  *scanin3    {prompt = "Internal use"}                       # OLDP-4

begin
    char    l_inimages = ""
    char    l_outspectra = ""
    char    l_outprefix = ""
    char    l_nsappwavedb = ""
    real    l_crval, l_cdelt, l_cradius
    bool    l_fl_phu
    char    l_key_instrument = ""
    char    l_key_camera = ""
    char    l_key_grating = ""
    char    l_key_filter = ""
    char    l_key_prism = ""
    char    l_key_fpmask = ""
    char    l_key_order = ""
    char    l_key_arrayid = ""
    char    l_key_dispaxis = ""
    char    l_key_wave = ""
    char    l_key_delta = ""
    char    l_key_waveorder = ""
    char    l_key_wavevar = ""
    char    l_key_cradius = ""
    char    l_key_mode = ""
    char    l_key_pixscale = ""
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_logfile = ""
    bool    l_verbose
    bool    l_debug

    int     junk, nbad, nin, istart, iend, version, count, nfiles, idxlastextn
    int     dispaxis, spatial, iorder, exti
    char    tmpin, tmpout, tmpexpand, tmpall, tmpextn
    char    inimg, outimg, extn, line, inphu, outphu, inherit, inextn
    char    camera, filter, prism, grating, fpmask, order, arrayid
    real    wave, wavevar
    int     waveorder
    char    extns[3], badhdr, keyfound, unused
    bool    intdbg, single, first, usetable, dowarn, isifu, f2
    struct  sdate, sline
    real    lambda, delta, center, hdelta, radius
    real    pixscale
    char    filter1, filter2, lyot

    junk = fscan (inimages, l_inimages)
    junk = fscan (outspectra, l_outspectra)
    junk = fscan (outprefix, l_outprefix)
    junk = fscan (nsappwavedb, l_nsappwavedb)
    l_crval = crval
    l_cdelt = cdelt
    l_cradius = cradius
    l_fl_phu = fl_phu
    junk = fscan (logfile, l_logfile)
    l_debug = debug
    l_verbose = verbose || l_debug
    
    status = 1

    badhdr = ""
    junk = fscan (nsheaders.key_instrument, l_key_instrument)
    if ("" == l_key_instrument) badhdr = badhdr + " key_instrument"
    junk = fscan (nsheaders.key_camera, l_key_camera)
    if ("" == l_key_camera) badhdr = badhdr + " key_camera"
    junk = fscan (nsheaders.key_grating, l_key_grating)
    if ("" == l_key_grating) badhdr = badhdr + " key_grating"
    junk = fscan (nsheaders.key_filter, l_key_filter)
    if ("" == l_key_filter) badhdr = badhdr + " key_filter"
    junk = fscan (nsheaders.key_prism, l_key_prism)
    if ("" == l_key_prism) badhdr = badhdr + " key_prism"
    junk = fscan (nsheaders.key_fpmask, l_key_fpmask)
    if ("" == l_key_fpmask) badhdr = badhdr + " key_fpmask"
    junk = fscan (nsheaders.key_order, l_key_order)
    if ("" == l_key_order) badhdr = badhdr + " key_order"
    junk = fscan (nsheaders.key_arrayid, l_key_arrayid)
    if ("" == l_key_arrayid) badhdr = badhdr + " key_arrayid"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (nsheaders.key_wave, l_key_wave)
    if ("" == l_key_wave) badhdr = badhdr + " key_wave"
    junk = fscan (nsheaders.key_delta, l_key_delta)
    if ("" == l_key_delta) badhdr = badhdr + " key_delta"
    junk = fscan (nsheaders.key_wavevar, l_key_wavevar)
    if ("" == l_key_wavevar) badhdr = badhdr + " key_wavevar"
    junk = fscan (nsheaders.key_waveorder, l_key_waveorder)
    if ("" == l_key_waveorder) badhdr = badhdr + " key_waveorder"
    junk = fscan (nsheaders.key_cradius, l_key_cradius)
    if ("" == l_key_cradius) badhdr = badhdr + " key_cradius"
    junk = fscan (nsheaders.key_mode, l_key_mode)
    if ("" == l_key_mode) badhdr = badhdr + " key_mode"
    junk = fscan (nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"    
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"

    intdbg = no
    usetable = no
    extns[1] = l_sci_ext
    extns[2] = l_var_ext
    extns[3] = l_dq_ext
    unused = "(unused)"
    f2 = no

    tmpin = mktemp ("tmpin")
    tmpout = mktemp ("tmpout")
    tmpexpand = mktemp ("tmpexpand")
    tmpall = mktemp ("tmpall")
    
    tmpextn = "" # set later in for-loop

    cache ("keypar", "nswedit", "gemextn", "gemdate")


    # Start logging

    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile) 
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSAPPWAVE: Both nsappwave.logfile and \
                gnirs.logfile are empty.", l_logfile, verbose+) 
            printlog ("                     Using default file " \
                // l_logfile // ".", l_logfile, verbose+) 
        }
    }
    date | scan (sdate) 
    printlog ("--------------------------------------------------------\
        ----------------------", l_logfile, verbose = l_verbose) 
    printlog ("NSAPPWAVE -- " // sdate, l_logfile, verbose = l_verbose) 
    printlog (" ", l_logfile, verbose = l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NSAPPWAVE: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }


    # Check files exist, are MEF, and have been through nsprepare 
    # (and set order if required)

    if (intdbg) print ("running gemextn on file lists")

    gemextn (l_inimages, check="exists,mef", process="none", \
        index="", extname="", extversion="", ikparams="", \
        omit="kernel,exten", replace="", outfile=tmpin, logfile="",
        glogpars="", verbose=l_verbose)
    nbad = gemextn.fail_count
    nin = gemextn.count
    if (0 == nin) {
        printlog ("ERROR - NSAPPWAVE: No input images found for " \
            // l_inimages, l_logfile, verbose+) 
        goto clean
    }

    if (intdbg) print ("check for nsprepare")

    scanin1 = tmpin
    while (fscan (scanin1, inimg) != EOF) {

        if (intdbg) print ("inimg: " // inimg)

        keyfound=""
        hselect (inimg//"[0]", "*PREPAR*", yes) | scan(keyfound)
        if (keyfound == "") {
            printlog ("ERROR - NSAPPWAVE: Image " // inimg \
                // " not PREPAREd.", l_logfile, l_verbose) 
            nbad += 1
        }
    }

    if (nbad > 0) {
        printlog ("ERROR - NSAPPWAVE: " // nbad // " image(s) either do \
            not exist, are not MEF files, or", l_logfile, verbose+) 
        printlog ("                   have not been run through \
            *PREPARE", l_logfile, verbose+) 
        goto clean
    }


    # Generate the output images

    if (intdbg) print ("output check")

    gemextn (l_outspectra, check="absent", process="none", index="", \
        extname="", extversion="", ikparams="", omit="kernel,exten", 
        replace="", outfile=tmpout, logfile="", glogpars="",
        verbose=l_verbose)
    if (gemextn.fail_count != 0) {
        printlog ("ERROR - NSAPPWAVE: Existing or incorrectly formatted \
            output files", l_logfile, verbose+) 
        goto clean
    }

    # If tmpout is empty, the output files names should be 
    # created with the outprefix parameter
    # (we could have a separate msg for outprefix="", but that will
    # trigger an error in gemextn anyway)
    if (gemextn.count == 0) {

        if (intdbg) print ("output from substitution")

        gemextn ("%^%" // l_outprefix // "%" // "@" // tmpin, \
            check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", omit="kernel,exten", \
            replace="", outfile=tmpout, logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - NSAPPWAVE: No or incorrectly formatted \
                output files", l_logfile, verbose+) 
            goto clean
        }
    }

    if (nin != gemextn.count) {
        printlog ("ERROR - NSAPPWAVE: Numberof input and output files \
            do not match", l_logfile, verbose+) 
        goto clean
    }


    # Expand input images into extensions and merge information
    # into a single list before processing.

    if (intdbg) print ("expansion and collation")

    scanin1 = tmpin
    scanin2 = tmpout
    nfiles = 0
    count = 0
    nbad = 0
    while (fscan (scanin1, inimg) != EOF) {
        junk = fscan (scanin2, outimg) # already checked numbers ok
        nfiles = nfiles + 1
        first = yes

        delete (tmpexpand, verify-, >& "dev$null")
        gemextn (inimg, check="exists", process="expand", index="", \
            extname=l_sci_ext, extversion="1-", ikparams="", \
            omit="index,name,params,section", replace="", 
            outfile=tmpexpand, logfile="", glogpars="", verbose=l_verbose)

        if (0 == gemextn.count) {
            printlog ("ERROR - NSAPPWAVE: No science data in " // inimg, \
                l_logfile, l_verbose) 
            nbad += 1
        } else {
            single = (1 == gemextn.count && l_fl_phu)

            scanin3 = tmpexpand
            while (fscan (scanin3, extn) != EOF) {
                count = count + 1

                # ugly - need a parse task or change to gemextn output
                istart = stridx ("=", extn)
                iend = stridx ("]", extn)
                version = int (substr (extn, istart+1, iend-1))

                line = inimg // " " // outimg // " " // version // " " \
                    // single // " " // first
                if (intdbg) print (line)
                print (line, >> tmpall)

                first = no
            }
        }
    }

    if (nbad > 0) {
        printlog ("ERROR - NSAPPWAVE: " // nbad // " image(s) do not \
            contain science data", l_logfile, verbose+) 
        goto clean
    }
    printlog ("NSAPPWAVE: Processing " // count // " extension(s) from " \
        // nfiles // " file(s).", l_logfile, l_verbose) 

    printlog ("  Camera     Grating    Filter     Prism      FPMask     \
        Wave   Delta Order Axis", l_logfile, l_verbose)

    # For each, process

    if (intdbg) print ("process")

    scanin1 = tmpall
    while (fscan (scanin1, inimg, outimg, version, single, first) != EOF) {
        extn = inimg // "[" // l_sci_ext // "," // version // "]"
        inherit = inimg // "[" // l_sci_ext // "," // version \
            // ",inherit]"

        if (first) {
            # Create output file from input's PHU
            fxcopy (inimg//".fits", outimg//".fits", group="0", new_file+,
                verbose-)

            # Propagate the MDF, if present
            gemextn (inimg//".fits[MDF]", check="exists", process="none",
                index="", extname="", extver="", ikparams="", omit="",
                replace="", outfile="dev$null", logfile="dev$null", verbose-)
            if (gemextn.count == 1) {
                # There is an MDF. Copy to output.
                fxhead (inimg//".fits", format="", long=no, count=no) | \
                    match ("MDF", "STDIN", stop-, print+, meta+) | scan(exti)
                fxinsert (inimg//".fits", outimg//".fits[0]", group=exti, \
                    verbose+)
                idxlastextn = 1
            } else {
                # Reset output extensions counter
                idxlastextn = 0
            }
            

            # Read information from PHU
            inphu = inimg//"[0]"

            keypar (inphu, l_key_mode, silent+)
            isifu = no
            if (keypar.found)
                isifu = "IFU" == keypar.value

            keypar (inphu, l_key_instrument, silent+)
            if (keypar.found) {
                if (keypar.value == "F2" || keypar.value == "Flam")
                    f2 = yes
            } else if (unused != l_key_instrument) {
                printlog ("WARNING - NSAPPWAVE: No value found for " \
                    // l_key_instrument // " (key_instrument) in " \
                    // inphu // ".", l_logfile, l_verbose)
            }

            # KL 20220217 - F2 headers have modified. Broken. Fixing.
            camera = ""
            if (f2) {
                keypar (inphu, l_key_pixscale, silent+)
            } else {
                keypar (inphu, l_key_camera, silent+)
            }
            if (keypar.found) {
                if (f2) {
                    pixscale = real( keypar.value )
                    if ( (pixscale > 0.10) && (pixscale < 0.20) ) {
                        camera = "f/16_G5830"
                    } else {
                        keypar (inphu, l_key_camera, silent+)
                        if (keypar.found) {
                            camera = keypar.value
                        }
                    }
                } else {
                    camera = keypar.value
                }
                usetable = yes
                printlog ("NSAPPWAVE: camera found " // camera, \
                    l_logfile, l_debug)
            } else if (unused != l_key_camera) {
                printlog ("WARNING - NSAPPWAVE: No value found for " \
                    // l_key_camera // " (key_camera) in " // inphu // ".", \
                    l_logfile, l_verbose)
            }

            if (camera == "NIFS")
                l_key_fpmask = "(unused)"

            keypar (inphu, l_key_grating, silent+)
            grating = ""
            if (keypar.found) {
                grating = keypar.value
                usetable = yes
                printlog ("NSAPPWAVE: grating found " // grating, \
                    l_logfile, l_debug)
            } else if (unused != l_key_grating) {
                printlog ("WARNING - NSAPPWAVE: No value found for " \
                    // l_key_grating // " (key_grating) in " // inphu \
                    // ".", l_logfile, l_verbose)
            }

            # KL 20220217
            # Now get back that component ID that was removed in f2prepare...
            # Also, check the LYOT keyword, they are now putting filters in there.

            filter = ""
            if (f2) {
                keypar (inphu, l_key_filter, silent+)
                if (keypar.found) {
                    filter = keypar.value
                    usetable = yes
                    
                    keypar (inphu, l_key_filter//"1", silent+)
                    if (keypar.found) {
                        filter1 = keypar.value
                    } else { # early ENG data
                        keypar (inphu, "FILT1POS", silent+)
                        if (keypar.found) {
                            filter1 = keypar.value
                        }
                    }
                    keypar (inphu, l_key_filter//"2", silent+)
                    if (keypar.found) {
                        filter2 = keypar.value
                    } else {  #  early ENG data
                        keypar (inphu, "FILT2POS", silent+)
                        if (keypar.found) {
                            filter2 = keypar.value
                        }
                    }
                    keypar (inphu, "LYOT", silent+)
                    if (keypar.found) {
                        lyot = keypar.value
                    } else { # early ENG data
                        keypar (inphu, "LYOTPOS", silent+)
                        if (keypar.found) {
                            lyot = keypar.value
                        }
                    }

                    if ( (filter == "Open") && \
                         (  (substr(lyot, 1, 1) != "f") && \
                            (substr(lyot, 1, 4) != "GEMS") && \
                            (substr(lyot, 1, 4) != "Hart") \
                          ) ) {
                        filter = lyot
                    } else if ( (substr(filter1, 1, (strlen (filter1) - 6)) == filter) || \
                         (substr(filter1, 1, (strlen (filter1) - 11)) == filter) ) {
                        filter = filter1
                    } else if ( (substr(filter2, 1, (strlen (filter2) - 6)) == filter) || \
                         (substr(filter2, 1, (strlen (filter2) - 11)) == filter) ) {
                        filter = filter2
                    } 
                    printlog ("NSAPPWAVE: filter found " // filter, \
                        l_logfile, l_debug)
                } else if (unused != l_key_filter) {
                    printlog ("WARNING - NSAPPWAVE: No value found for " \
                        // l_key_filter // " (key_filter) in " // inphu // ".", \
                        l_logfile, l_verbose)
                }

            } else {
                keypar (inphu, l_key_filter, silent+)
                if (keypar.found) {
                    filter = keypar.value
                    usetable = yes
                    printlog ("NSAPPWAVE: filter found " // filter, \
                        l_logfile, l_debug)
                } else if (unused != l_key_filter) {
                    printlog ("WARNING - NSAPPWAVE: No value found for " \
                        // l_key_filter // " (key_filter) in " // inphu // ".", \
                        l_logfile, l_verbose)
                }
            }

            keypar (inphu, l_key_prism, silent+)
            prism = ""
            if (keypar.found) {
                prism = keypar.value
                usetable = yes
                printlog ("NSAPPWAVE: prism found " // prism, \
                    l_logfile, l_debug)
            } else if (unused != l_key_prism) {
                printlog ("WARNING - NSAPPWAVE: No value found for " \
                    // l_key_prism // " (key_prism) in " // inphu // ".", \
                    l_logfile, l_verbose)
            }

            keypar (inphu, l_key_fpmask, silent+)
            fpmask = ""
            if (keypar.found) {
                fpmask = keypar.value
                usetable = yes
                printlog ("NSAPPWAVE: fpmask found " // fpmask, \
                    l_logfile, l_debug)
            } else if (unused != l_key_fpmask) {
                printlog ("WARNING - NSAPPWAVE: No value found for " \
                    // l_key_fpmask // " (key_fpmask) in " // inphu // ".", \
                    l_logfile, l_verbose)
            }
            
            keypar (inphu, l_key_arrayid, silent+)
            arrayid = ""
            if (keypar.found) {
                arrayid = keypar.value
                usetable = yes
                printlog ("NSAPPWAVE: arrayid found " // arrayid, \
                    l_logfile, l_debug)
            } else if (unused != l_key_arrayid) {
                printlog ("WARNING - NSAPPWAVE: No value found for " \
                    // l_key_arrayid // " (key_arrayid) in " // inphu // ".", \
                    l_logfile, l_verbose)
            }

            keypar (inphu, l_key_dispaxis, silent+)
            if (intdbg) print ("dispaxis: " // keypar.value)
            dispaxis = 1
            if (keypar.found)
                dispaxis = int (keypar.value)
            else {
                printlog ("WARNING - NSAPPWAVE: No value found for " \
                    // l_key_dispaxis // " (key_dispaxis) in " // inphu \
                    // ".", l_logfile, l_verbose)
            }
            spatial = 3 - dispaxis

            keypar (inphu, l_key_wave, silent+)
            wave = 0.0
            if (keypar.found)
                wave = real (keypar.value)
            else if (unused != l_key_wave) {
                printlog ("WARNING - NSAPPWAVE: No value found for " \
                    // l_key_wave // " (key_wave) in " // inphu // ".", \
                    l_logfile, l_verbose)
            }

            keypar (inphu, l_key_delta, silent+)
            hdelta = 0.0
            if (keypar.found)
                hdelta = real (keypar.value)
            else if (unused != l_key_delta) {
                printlog ("WARNING - NSAPPWAVE: No value found for " \
                    // l_key_delta // " (key_delta) in " // inphu // ".", \
                    l_logfile, l_verbose)
            }

            keypar (inphu, l_key_waveorder, silent+)
            waveorder = INDEF
            if (keypar.found)
                waveorder = int (keypar.value)
            else if (unused != l_key_waveorder) {
                printlog ("WARNING - NSAPPWAVEORDER: No value found for " \
                    // l_key_waveorder // " (key_waveorder) in " // inphu \
                    // ".", l_logfile, l_verbose)
            }

            usetable = usetable && "" != l_nsappwavedb

            printf ("%20s -> %20s (%2d)\n", inimg, outimg, \
                version) | scan (sline)
            printlog (sline, l_logfile, l_verbose) 

            printf ("  %10s %10s %10s %10s %10s %8.2f %5.2f %2d %1d\n", \
                camera, grating, filter, prism, fpmask, \
                wave, hdelta, waveorder, dispaxis) | scan (sline)
            printlog (sline, l_logfile, l_verbose) 

            dowarn = yes # once per image
        }

        # short-cut for IFU data

        if (first || no == isifu) {

            keypar (inherit, l_key_wavevar, silent+)
            wavevar = INDEF
            if (keypar.found)
                wavevar = real (keypar.value)
            else if (unused != l_key_wavevar) {
                printlog ("WARNING - NSAPPWAVE: No value found for " \
                    // l_key_wavevar // " (key_wavevar) in " // inherit \
                    // ".", l_logfile, l_verbose)
            }

            keypar (extn, l_key_order, silent+)
            iorder = -1
            if (keypar.found) {
                iorder = int (keypar.value)
                usetable = yes
            } else if (unused != l_key_wavevar) {
                printlog ("WARNING - NSAPPWAVE: No value found for " \
                    // l_key_order // " (key_order) in " // extn // ".", \
                    l_logfile, l_verbose)
                printlog ("                     Assuming -1.", \
                    l_logfile, l_verbose)
            }
            if (isindef (iorder))
                iorder = waveorder
            if (-1 == iorder)
                iorder = waveorder

            # If the values for the crval, cdelt and cradius parameters are
            # not INDEF, these value will be used for lambda (the central
            # wavelengths, in Angstroms), delta (the dispersion, in Angstroms
            # per pixel) and radius (the search radius, in pixels),
            # respectively, regardless of any values defined in the header or
            # the lookup table.
            lambda = l_crval
            delta = l_cdelt
            radius = l_cradius

            # If lambda, delta and radius are INDEF, check if any values are
            # defined in the lookup table 
            if (usetable) {
                if (intdbg) print ("trying table")
                # (minmatch = 2 to hit order plus something else)
                nswedit (table=l_nsappwavedb, camera=camera, grating=grating, \
                    filter=filter, prism=prism, mask=fpmask, arrayid=arrayid, \
                    order=iorder, \
                    centre=wavevar, range=INDEF, minmatch=2, append-, \
                    overwrite-, create-, logfile=l_logfile, verbose=l_debug)
                if (0 == nswedit.status) {
                    if (isindef (lambda) && isindef (nswedit.lambda) == no) {
                        lambda = nswedit.lambda
                    } else {
                        if (dowarn) {
                            printlog ("NSAPPWAVE: Wavelength not from table \
                                for " // inimg // ".", l_logfile, l_verbose)
                        dowarn = no
                        }
                    }
                    if (isindef (delta) && isindef (nswedit.delta) == no) {
                        delta = nswedit.delta
                    } else {
                        if (dowarn) {
                            printlog ("NSAPPWAVE: Dispersion not from table \
                                for " // inimg // ".", l_logfile, l_verbose)
                        dowarn = no
                        }
                    }
                    if (isindef (radius) && isindef (nswedit.cradius) == no) {
                        radius = nswedit.cradius
                    } else {
                        if (dowarn) {
                            printlog ("NSAPPWAVE: cradius not from table \
                                for " // inimg // ".", l_logfile, l_verbose)
                        dowarn = no
                        }
                    }
                }
            }

            # If there are no values defined in the lookup table (i.e., lambda,
            # delta and radius are still INDEF), check if any values are
            # defined in the header
            if (isindef (lambda) && wave != 0.0) {
                if (f2)
                    # wave is already in Angstroms for FLAMINGOS-2 data
                    lambda = wave
                else
                    # Convert wave from microns to Angstroms 
                    lambda = 10000.0 * wave
                if (intdbg) print ("wavelength from header: " // lambda)
            }

            if (isindef (delta) && hdelta != 0.0) {
                if (f2)
                    # hdelta is already in Angstroms for FLAMINGOS-2 data
                    delta = hdelta
                else
                    # Convert hdelta from microns to Angstroms 
                    delta = 10000.0 * hdelta
                if (intdbg) print ("dispersion from header: " // delta)
            }

            # If there are no values defined in the header (i.e., lambda, delta
            # and radius are still INDEF), no valid value could be found. It's
            # ok for radius to not have a valid value though :)

            if (isindef (lambda)) {
                printlog ("ERROR - NSAPPWAVE: Unknown central wavelength.", \
                    l_logfile, verbose+) 
                printf ("  %10s %10s %10s %10s %10s                \
                    %2d %1d\n", camera, grating, filter, prism, fpmask, \
                    iorder, dispaxis) | scan (sline)
                printlog (sline, l_logfile, l_verbose) 
                goto clean
            }

            if (isindef (delta)) {
                printlog ("ERROR - NSAPPWAVE: Unknown dispersion.", \
                    l_logfile, verbose+) 
                printf ("  %10s %10s %10s %10s %10s                \
                    %2d %1d\n", camera, grating, filter, prism, fpmask, \
                    iorder, dispaxis) | scan (sline)
                printlog (sline, l_logfile, l_verbose) 
                goto clean
            }

            if (no == dowarn && first) {
                printlog ("           (further table messages \
                    supressed for this file).", l_logfile, l_verbose)
            }
        }

        printf ("  extension %2d: %7f %5f %5f\n", \
            version, lambda, delta, radius) | scan (sline)
        printlog (sline, l_logfile, l_verbose)

        if (intdbg) print (lambda // ", " // delta // l_verbose)

        # pixel position in the middle of the array
        extn = inimg // "[" // l_sci_ext // "," // version // "]"
        keypar (extn, "i_naxis" // dispaxis, silent+)
        if (keypar.found) {
            center = 0.5 * real (keypar.value)
        } else {
            printlog ("ERROR - NSAPPWAVE: Missing dimension data.", \
                l_logfile, verbose+) 
            goto clean
        }


        # Update headers

        if (intdbg) print ("updating headers")

        if (first) {
            outphu = outimg // "[0]"
            gemdate ()
            gemhedit (outphu, "NSAPPWAV", gemdate.outdate, 
                "UT Time stamp for NSAPPWAVE", delete-)
            gemhedit (outphu, "GEM-TLM", gemdate.outdate, 
                "UT Last modification with GEMINI IRAF", delete-)
        }

        for (i = 0; i <= 3; i = i + 1) {

            if (0 == i) {
                if (single)
                    extn = outimg // "[0]"
                else
                    next
            } else {
                #extn = outimg // "[" // extns[i] // "," // version // "]"
                #if (no == imaccess (extn)) next
                
                inextn = inimg//"["//extns[i]//","//version//"]"
                if (no == imaccess (inextn))
                    next
                
                # Create tmp file name used within this for-loop
                tmpextn = mktemp("tmpextn-")
                imcopy (inextn, tmpextn, verbose-)
                extn = tmpextn
            }
            if (intdbg) print ("headers: " // extn)

            # Reduce this to avoid mis-matches and errors
            # By this point in the process all images will be 2D,
            # although some may have started larger
            gemhedit (extn, "WCSDIM", 2, "", delete+)
            gemhedit (extn, "CD" // dispaxis // "_" // spatial, 0, "", delete+)
            gemhedit (extn, "CD" // spatial // "_" // dispaxis, 0, "", delete+)

            gemhedit (extn, "WAT" // dispaxis // "_000", "system=world", "", \
                delete-)
            gemhedit (extn, "WAT" // dispaxis // "_001", "wtype=linear", "", \
                delete-)
            gemhedit (extn, "WAT" // dispaxis // "_002", \
                "units=Angstroms label=Wavelength", "", delete-)
            gemhedit (extn, "CTYPE" // dispaxis, "LINEAR", \
                "coordinate type for the dispersion axis", delete-)

            if (intdbg) print ("CRPIX" // dispaxis // " = " // center)
            gemhedit (extn, "CRPIX" // dispaxis, center, \
                "coordinate of spectral reference pixel", delete-)
            gemhedit (extn, "CRVAL" // dispaxis, lambda, 
                "wavelength at spectral reference pixel", delete-)
            gemhedit (extn, "CD" // dispaxis // "_" // dispaxis, delta,
                "dispersion (Angstrom/pix)", delete-)

            if (no == isindef (cradius)) {
                gemhedit (extn, key_cradius, radius, \
                    "search radius for wavelength cal", delete-)
            }

            # The LTV1/2 keywords map the offset of the cut position to the
            # origin (in physical units) of the image from which it was cut. 
            # The presence of these keywords can prevent the correct
            # determination of the wavelength solution (due to this shift in
            # the WCS). Delete the physical offset (LTV) keywords so the WCS is
            # now in the frame of reference of the cut image and does not map
            # to the parent image coordinate system. This leaves LTM keywords
            # in place, as is. This is fine as they represent the
            # scale/rotation matrix which is still needed.
            gemhedit (extn, "LTV1", "", "", delete+)
            gemhedit (extn, "LTV2", "", "", delete+)

            gemhedit (extn, "WAT" // spatial // "_000", "system=physical", \
                "", delete-)
            gemhedit (extn, "WAT" // spatial // "_001", "wtype=linear", "", \
                delete-)
            gemhedit (extn, "CTYPE" // spatial, "LINEAR", \
                "coordinate type for the spatial axis", delete-)
            gemhedit (extn, "CRPIX" // spatial, 1., \
                "coordinate of the spatial reference pixel", delete-)
            gemhedit (extn, "CRVAL" // spatial, 1., \
                "spatial axis value at spatial reference pixel", delete-)
            gemhedit (extn, "CD" // spatial // "_" // spatial, 1., 
                "partial of spatial axis wrt spatial dim", delete-)

            gemhedit (extn, "DC-FLAG", 0, "", delete-)
            gemhedit (extn, "DISPAXIS", dispaxis, "Dispersion axis", delete-)

            if (no == single) {
                # append the tmp extension to the output file
                # Note: extn = tmpextn
                if (intdbg) print (outimg//"["//idxlastextn//"]")
                fxinsert (extn, outimg//".fits["//idxlastextn//"]", \
                    group="", verbose-)
                idxlastextn = idxlastextn + 1
                gemhedit (outimg//"["//idxlastextn//"]", "EXTNAME", extns[i],
                    "Extension name", delete-)
                gemhedit (outimg//"["//idxlastextn//"]", "EXTVER", version,
                    "Extension version", delete-)
                imdelete (tmpextn, verify-, >& "dev$null")
            }

        } # end for-loop through sci,var,dq
    } # end while-loop through input sci extns

    status = 0 # Everything ended OK


clean:

    scanin1 = ""
    scanin2 = ""
    scanin3 = ""

    if (intdbg)
        print ("delete: " // tmpin // ", " // tmpout // ", " \
            // tmpexpand // ", " // tmpall)

    delete (tmpin, verify-, >& "dev$null") 
    delete (tmpout, verify-, >& "dev$null") 
    delete (tmpexpand, verify-, >& "dev$null") 
    delete (tmpall, verify-, >& "dev$null") 

    printlog ("", l_logfile, l_verbose) 
    if (status == 0) 
        printlog ("NSAPPWAVE Exit status good", l_logfile, l_verbose) 

    printlog ("-------------------------------------------------------\
        -----------------------", l_logfile, l_verbose) 

end
