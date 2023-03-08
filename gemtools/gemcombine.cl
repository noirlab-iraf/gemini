# Copyright(c) 2002-2013 Association of Universities for Research in Astronomy, Inc.

procedure gemcombine (input, output)

# Combine MEF files
#
# Version Feb 28, 2002  BM  v1.3 release
#         Aug 10, 2002  BM,IJ  v1.4 release
#         Nov 20, 2002  IJ allow rejection with only 2 input images
#         Aug 12, 2003  KL IRAF2.12 - hedit: addonly- ; imcombine: headers,
#                           bpmasks,expmasks,outlimits ="", rejmask->rejmasks,
#                           plfile->nrejmasks
#         Nov 07, 2003  KL change default nkeep from 0 to 1

string input        {prompt = "Input MEF images"}
string output       {prompt = "Output MEF image"}

string title        {"", prompt = "Title for output SCI plane"}
string combine      {"average", enum = "average|median", prompt = "Combination operation"}
string reject       {"avsigclip", enum = "none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip", prompt = "Rejection algorithm"}
string offsets      {"none", prompt = "Input image offsets"}
string masktype     {"none", enum="none|goodvalue", prompt = "Mask type"}
real   maskvalue    {0., prompt = "Mask value"}
string scale        {"none", prompt = "Image scaling"}
string zero         {"none", prompt = "Image zeropoint offset"}
string weight       {"none", prompt = "Image weights"}
string statsec      {"[*,*]", prompt = "Statistics section"}
string expname      {"EXPTIME", prompt = "Exposure time header keyword"}
real   lthreshold   {INDEF, prompt = "Lower threshold"}
real   hthreshold   {INDEF, prompt = "Upper threshold"}
int    nlow         {1, min = 0, prompt = "minmax: Number of low pixels to reject"}
int    nhigh        {1, min = 0, prompt = "minmax: Number of high pixels to reject"}
int    nkeep        {1, prompt = "Minimum to keep or maximum to reject"}
bool   mclip        {yes, prompt = "Use median in sigma clipping algorithms?"}
real   lsigma       {3., min = 0, prompt = "Lower sigma clipping factor"}
real   hsigma       {3., min = 0, prompt = "Upper sigma clipping factor"}
string key_ron      {"RDNOISE", prompt = "Keyword for readout noise in e-"}
string key_gain     {"GAIN", prompt = "Keyword for gain in electrons/ADU"}
real   ron          {0.0, min = 0., prompt = "Readout noise rms in electrons"}
real   gain         {1.0, min = 0.00001, prompt = "Gain in e-/ADU"}
string snoise       {"0.0", prompt = "ccdclip: Sensitivity noise (electrons)"}
real   sigscale     {0.1, min = 0., prompt = "Tolerance for sigma clipping scaling correction"}
real   pclip        {-0.5, prompt = "pclip: Percentile clipping parameter"}
real   grow         {0.0, min = 0., prompt = "Radius (pixels) for neighbor rejection"}
string bpmfile      {"", prompt = "Name of bad pixel mask file or image."}
string nrejfile     {"", prompt = "Name of rejected pixel count image."}
string sci_ext      {"SCI", prompt = "Name(s) or number(s) of science extension"}
string var_ext      {"VAR", prompt = "Name(s) or number(s) of variance extension"}
string dq_ext       {"DQ", prompt = "Name(s) or number(s) of data quality extension"}
bool   fl_vardq     {no, prompt = "Make variance and data quality planes?"}
string logfile      {"gemcombine.log", prompt = "Log file"}
bool   fl_dqprop    {no, prompt = "Propagate all DQ values?"}
bool   verbose      {yes, prompt = "Verbose output?"}
int    status       {0, prompt = "Exit status (0=good)"}
struct  *flist      {"", prompt = "Internal use only"}
struct  *scanfile   {"", prompt = "Internal use only"}

begin

        string l_input, l_output, l_statsec, l_expname, l_bpmfile
        string l_logfile, l_offsets, l_masktype
        string l_combine, l_reject, l_scale, l_zero, l_weight
        string l_key_ron, l_key_gain, l_snoise, l_orig_statsec
        string l_sgain, l_sron
        real   l_ron, l_gain, l_lthreshold, l_hthreshold, l_lsigma, l_hsigma
        real   l_grow, l_sigscale, l_pclip, l_maskvalue
        int    l_nlow, l_nhigh, l_nkeep
        bool   l_verbose, s_verbose, l_mclip, l_fl_dqprop, l_fl_vardq, debug
        string l_sci_ext, l_var_ext, l_dq_ext, l_title, tsci, tvar, tdq
        string sextlist, vextlist, dextlist, extlist, sci[600], var[600]
        string dq[600], filelist, combout, combdq, combsig, img, dqsum, dqfits
        string dqsumold, suf, tmpdq, scilist, dqlist, plimg, tmplog
        string sect, infile[500], mdf, exphu, phuimg, bpmextn, bpm_extn
        real mean, epadu, rdnoise, gaineff, roneff
        int i, j, k, l, n, idum, nbad, len, nsci, nvar, ndq, junk
        int nimg, nextnd, secpos1, secpos2, nsciext, lcount
        bool useextver, usemdf, l_cleanbpm, isscale_list, bpmismef, bpmchecked
        struct sdate, line
        char l_nrejfile = ""
        char dqexp, mdqsum, sclist, scale_list, bpmname, tmpexprlist, dqpr_expr
        char statmods[3], l_statmod
        int nmodifiers, jj

        # Query parameters
        l_input = input; l_output = output
        l_title = title
        l_logfile = logfile
        l_combine = combine; l_reject = reject; l_offsets = offsets
        l_masktype= masktype
        l_maskvalue = maskvalue
        l_scale = scale; l_zero = zero; l_weight = weight
        l_statsec = statsec; l_expname = expname
        l_lthreshold = lthreshold; l_hthreshold = hthreshold
        l_nlow = nlow; l_nhigh = nhigh; l_nkeep = nkeep
        l_lsigma = lsigma; l_hsigma = hsigma
        l_verbose = verbose; l_mclip = mclip
        l_key_gain = key_gain; l_key_ron = key_ron
        l_gain = gain; l_ron = ron; l_snoise = snoise
        l_sgain = l_key_gain; l_sron = l_key_ron
        l_sigscale = sigscale; l_pclip = pclip; l_grow = grow
        l_sci_ext = sci_ext; l_dq_ext = dq_ext; l_var_ext = var_ext
        l_bpmfile = bpmfile; l_fl_dqprop = fl_dqprop; l_fl_vardq = fl_vardq
        junk = fscan (nrejfile, l_nrejfile)

        l_orig_statsec = l_statsec
        status = 0
        useextver = no
        usemdf = no
        bpmismef = no
        bpmchecked = no
        debug = no

        if (debug) print ("in gemcombine")

        # Keep imgets parameters from changing by outside world
        cache ("imgets", "hedit", "keypar", "gimverify", "gemextn", "gemdate")

        s_verbose = hedit.show
        hedit.show = l_verbose

        # Start logging to file
        if (l_logfile == "STDOUT") {
            l_logfile = ""
            l_verbose = yes
        }

        date | scan (sdate)
        printlog ("-----------------------------------------------------------\
            ---------------------", l_logfile, verbose = l_verbose)
        printlog ("GEMCOMBINE -- " // sdate, l_logfile, verbose = l_verbose)
        printlog (" ", l_logfile, verbose = l_verbose)
        printlog ("GEMCOMBINE: input = " // l_input, l_logfile, \
            verbose = l_verbose)
        printlog ("GEMCOMBINE: output = " // l_output, l_logfile, \
            verbose = l_verbose)
        printlog ("GEMCOMBINE: title = " // l_title, l_logfile, \
            verbose = l_verbose)
        printlog ("GEMCOMBINE: bpmfile = " // l_bpmfile, l_logfile, \
            verbose = l_verbose)
        printlog ("GEMCOMBINE: sci_ext = " // l_sci_ext, l_logfile, \
            verbose = l_verbose)
        printlog ("GEMCOMBINE: var_ext = " // l_var_ext, l_logfile, \
            verbose = l_verbose)
        printlog ("GEMCOMBINE: dq_ext = " // l_dq_ext, l_logfile, \
            verbose = l_verbose)
        printlog (" ", l_logfile, verbose = l_verbose)

        # Verify that the extension names are not empty, otherwise exit
        # gracefully
        if (l_sci_ext == "" || l_sci_ext == " ") {
            printlog ("ERROR - GEMCOMBINE: Extension name sci_ext is \
                missing", l_logfile, verbose+)
            goto error
        }
        if (l_dq_ext == "" || l_dq_ext == " ") {
            printlog ("ERROR - GEMCOMBINE: Extension name dq_ext is \
                missing", l_logfile, verbose+)
            goto error
        }

        # Define temporary files
        filelist = mktemp ("tmpfilelist")
        scilist = mktemp ("tmpscilist")
        dqlist = mktemp ("tmpdqlist")
        dqsum = mktemp ("tmpdqsum")
        tmpdq = mktemp ("tmpdq")
        tmplog = mktemp ("tmplog")
        sextlist = mktemp ("tmpsextlist")
        vextlist = mktemp ("tmpvextlist")
        dextlist = mktemp ("tmpdextlist")
        extlist = mktemp ("tmpextlist")
        mdf = mktemp ("tmpmdf")
        mdqsum = mktemp("tmpmdqsum")
        scale_list = mktemp ("tmpscale_list")

        # Modifiers of image section
        nmodifiers = 3
        statmods[1] = "input"
        statmods[2] = "output"
        statmods[3] = "overlap"

        # temporary files wannabes (values assigned later)
        combout = ""
        combdq = ""
        combsig = ""
        dqfits = ""

        if ("" != l_nrejfile) {
            if (no == l_fl_vardq) {
                printlog ("WARNING - GEMCOMBINE: nrejfile given, but \
                    fl_vardq-, so no output.", l_logfile, verbose+)
            } else {
                gemextn (l_nrejfile, check="absent", process="none", \
                    index="", extname="", extversion="", ikparams="", \
                    omit="", replace="", outfile="dev$null", logfile="",
                    glogpars="", verbose=l_verbose)
                if (gemextn.fail_count > 0 || 1 != gemextn.count) {
                    printlog ("ERROR - GEMCOMBINE: problem with nrejfile.", \
                        l_logfile, verbose+)
                    goto error
                }
                gemextn (l_nrejfile, check="absent", process="none", \
                    index="", extname="", extversion="", ikparams="", \
                    omit="exten", replace="", outfile="STDOUT", \
                    logfile="dev$null", glogpars="", verbose-) | scan (line)
                l_nrejfile = line
            }
        }

        #check that there are input files
        if (l_input == "" || l_input == " ") {
            printlog ("ERROR - GEMCOMBINE: input files not specified", \
                l_logfile, verbose+)
            goto error
        }

        # check existence of list file
        if (substr (l_input, 1, 1) == "@") {
            len = strlen (l_input)
            if (no==access (substr (l_input, 2, len) ) ) {
                printlog ("ERROR - GEMCOMBINE: " // substr (l_input, 2, len) \
                    // " does not exist.", l_logfile, verbose+)
                goto error
            }
        }

        # Check for an image section
        sect = ""
        secpos1 = stridx ("[", l_input)
        secpos2 = stridx ("]", l_input)
        if (secpos1 > 0 && secpos2 > 0) {
            sect = substr (l_input, secpos1, secpos2)
            l_input = substr (l_input, 1, secpos1-1)
        } else if (secpos1 > 0 || secpos2 > 0) {
            printlog ("ERROR - GEMCOMBINE: mismatched brackets in " // \
                l_input, l_logfile, verbose+)
            goto error
        }
        files (l_input, sort-, > filelist)

        #check that an output file is given
        if (l_output == "" || l_output == " ") {
            printlog ("ERROR - GEMCOMBINE: output files not specified", \
                l_logfile, verbose+)
            goto error
        }

        # Check that the output file does not already exist. If so, exit.
        if (imaccess (l_output) ) {
            printlog ("ERROR - GEMCOMBINE: Output file " // l_output // " \
                already exists.", l_logfile, verbose+)
            goto error
        }

        #put suffix on output image
        len = strlen (l_output)
        suf = substr (l_output, (len-4), len)
        if (suf != ".fits") {
            l_output = l_output // ".fits"
        }

        # Make list of sci/var/dq extensions
        files (l_sci_ext, sort-, > sextlist)
        count (sextlist) | scan (nsci)
        files (l_var_ext, sort-, > vextlist)
        count (vextlist) | scan (nvar)
        files (l_dq_ext, sort-, > dextlist)
        count (dextlist) | scan (ndq)
        if ( (nsci != nvar) || (nsci != ndq) ) {
            printlog ("ERROR - GEMCOMBINE: different numbers of SCI, VAR, or \
                DQ extensions.", l_logfile, verbose+)
            goto error
        }

        if (nsci > 600) {
            printlog ("ERROR - GEMCOMBINE: the number of science extensions \
                is larger than sci, var, and dq array lengths", l_logfile, \
                verbose+)
            goto error
        }

        joinlines (sextlist // "," // vextlist // "," // dextlist, "", \
            output=extlist)
        delete (sextlist // "," // vextlist // "," // dextlist, verify-, \
            >& "dev$null")
        flist = extlist
        i = 0
        while (fscan (flist, tsci, tvar, tdq) != EOF) {
            i = i+1
            sci[i] = tsci
            var[i] = tvar
            dq[i] = tdq
        }
        flist = ""
        delete (extlist, verify-, >& "dev$null")

        # Check that all images in the list exist and are MEF
        nbad = 0
        k = 0
        flist = filelist
        while (fscan (flist, img) != EOF) {
            secpos1 = stridx ("[", img)
            secpos2 = stridx ("]", img)
            if (secpos1 > 0 && secpos2 > 0) {
                sect = substr (img, secpos1, secpos2)
                img = substr (img, 1, secpos1-1)
            } else if (secpos1 > 0 || secpos2 > 0) {
                printlog ("ERROR - GEMCOMBINE: mismatched brackets in " // \
                    img, l_logfile, verbose+)
                goto error
            }
            gimverify (img)
            if (gimverify.status > 0) {
                nbad = nbad+1
            } else {
                k = k+1
                # name w/o suffix
                infile[k] = gimverify.outname // ".fits"
                if (k == 1) {
                    imgets (infile[k] // "[0]", "NSCIEXT", >& "dev$null")
                    if (debug)
                        print (infile[k] // "[0]" // " nsciext: " // \
                        imgets.value)
                    nsciext = int (imgets.value)
                    if (nsciext > 600) {
                        printlog ("ERROR - GEMCOMBINE: the number of science \
                            extensions is larger than sci, var, and dq array \
                            lengths", l_logfile, verbose+)
                        goto error
                    }
                    # check for MDF
                    imgets (infile[k] // "[0]", "NEXTEND", >& "dev$null")
                    if (debug) print (imgets.value)
                    idum = int (imgets.value)
                    for (l = 1; l <= idum; l += 1) {
                        keypar (infile[k] // "[" // l // "]", "EXTNAME", \
                            silent+)
                        if (keypar.value == "MDF" && keypar.found) {
                            tcopy (infile[k] // "[MDF]", mdf // ".fits", \
                                verbose-)
                            usemdf = yes
                        }
                    }
                    # check if EXTVER used
                    if (nsci == 1 && nsciext > 0) {
                        if (imaccess (infile[k] // "[" // sci[k] // "," // \
                            nsciext // "]") ) {
                            tsci = sci[1]
                            tvar = var[1]
                            tdq = dq[1]
                            for (i = 1; i <= nsciext; i += 1) {
                                sci[i] = tsci // "," // i
                                var[i] = tvar // "," // i
                                dq[i] = tdq // "," // i
                            }
                            nsci = nsciext
                            useextver = yes
                        } else {
                            if (debug) {
                                print ("strange test failed")
                                print (infile[k] // "[" // sci[k] // "," // \
                                    nsciext // "]")
                            }
                        }
                    }
                }
            }
        }
        #end of while loop
        flist = ""
        nimg = k

        # Exit if problems found
        if (nbad > 0) {
            printlog ("ERROR - GEMCOMBINE: " // nbad // " image(s) either do \
                not exist or are not MEF files.", l_logfile, verbose+)
            goto error
        }

        # Check offset file
        if (l_offsets != "none" && l_offsets != "wcs" && l_offsets != \
            "world" && !access (l_offsets) ) {
            printlog ("ERROR - GEMCOMBINE: " // l_offsets // " Offset file \
                not found", l_logfile, verbose+)
            goto error
        }

        # Check scale parameter for a comma list of files, list must match nsci
        # and within that list number of lines must match nimg
        isscale_list = no
        if (stridx(",",l_scale) > 0) {
            files (l_scale, sort-, > scale_list)
            scanfile = scale_list
            while (fscan(scanfile, sclist) != EOF) {
                if (!access(sclist)) {
                    printlog ("ERROR - GEMCOMBINE: Cannot access scale file "//
                        sclist, l_logfile, verbose+)
                    goto error
                } else {
                    count (sclist) | scan (lcount)
                    if (lcount != nimg) {
                        printlog ("ERROR - GEMCOMBINE: "//sclist//" doesn't \
                            contain the same number of\n"//\
                            "                    lines as input science \
                            extensions "//nimg//" "//lcount, \
                            l_logfile, verbose+)
                        goto error
                    }
                }
            }
            isscale_list = yes
        }

        # images ok, so continue
        l_cleanbpm = no # clean headers of BPM after combining
        if (debug) print ("nsci: " // nsci)
        for (i = 1; i <= nsci; i += 1) {
            # temporary files used within this loop
            combout = mktemp("tmpcombout")
            combdq = mktemp("tmpcombdq")
            combsig = mktemp("tmpcombsig")
            dqfits = mktemp("tmpdqfits")

            # If scaling parameters are in a list
            if (isscale_list) {
                if (i == 1) {
                    head (scale_list, nlines="1") | scan (l_scale)
                } else {
                    tail (scale_list, nlines=-(i-1)) | \
                        head ("STDIN", nlines="1") | scan (l_scale)
                }
                l_scale = "@"//l_scale
            }

            n = 0
            gaineff = 0.0
            roneff = 0.0
            flist = filelist
            while (fscan (flist, img) != EOF) {
                # Check if an image section is present
                secpos1 = stridx ("[", img)
                secpos2 = stridx ("]", img)
                if (secpos1 > 0) {
                    sect = substr (img, secpos1, secpos2)
                    img = substr (img, 1, secpos1-1)
                    len = secpos1-1
                } else {
                    len = strlen (img)
                }
                suf = substr (img, len-4, len)
                # Strip the suffix of the input file:
                if (suf != ".fits" && imaccess (img // ".fits") ) {
                    img = img // ".fits"
                }
                if (no==imaccess (img // "[" // sci[i] // "]") ) {
                    printlog ("ERROR - GEMCOMBINE: Could not access " // img \
                        // "[" // sci[i] // "]", l_logfile, verbose+)
                    goto error
                }
                n = n+1
                # if n=1, save the name for the output phu
                if (n == 1) {
                    phuimg = img

                    # Parse the image section for this image if required
                    l_statsec = l_orig_statsec
                    if (l_statsec != "[*,*]") {

                        # Check for imcombine modifiers
                        for (jj = 1; jj <= nmodifiers; jj += 1) {
                            l_statmod = statmods[jj]
                            if (strstr(l_statmod,l_statsec) > 0) {
                                 l_statsec = substr(l_statsec,\
                                     strlen(l_statmod)+1,\
                                     strlen(l_statsec))
                                 break
                            } else {
                                l_statmod = ""
                            }
                        }

                        # Parse the statsec for this image
                        gemsecchk (img//"["//sci[i]//"]", l_statsec, \
                            logfile=l_logfile, verbose=debug)

                        if (gemsecchk.status != 0) {
                            printlog ("ERROR - GEMCOMBINE: GEMSECCHK reurned \
                                a non-zero status", l_logfile, verbose+)
                            goto error
                        } else {
                            l_statsec = l_statmod//gemsecchk.out_imgsect
                        }
                    }
                } # End of working setting things up from n == 1

                # check gain and readnoise values
                if (l_key_gain != "" && l_key_gain != " ") {
                    if (l_sgain != str (l_gain) ) {
                        imgets (img // "[" // sci[i] // "]", l_key_gain, \
                            >& "dev$null")
                        if (imgets.value == "0") {
                            imgets (img // "[0]", l_key_gain, >& "dev$null")
                        }
                        if (imgets.value == "0") {
                            # only warn if it's going to be used
                            if ((l_reject == "ccdclip") || \
                                (l_reject == "crreject")) {
                                printlog ("WARNING - GEMCOMBINE: keyword " // \
                                    l_key_gain // " not found in " // img, \
                                    l_logfile, verbose+)
                                printlog ("Using gain = " // str (l_gain), \
                                    l_logfile, l_verbose)
                            }
                            l_sgain = str (l_gain)
                        } else {
                            l_sgain = imgets.value
                        }
                    }
                } else {
                    l_sgain = str (l_gain)
                }
                if (l_key_ron != "" && l_key_ron != " ") {
                    if (l_sron != str (l_ron) ) {
                        imgets (img // "[" // sci[i] // "]", l_key_ron, \
                            >& "dev$null")
                        if (imgets.value == "0") {
                            imgets (img // "[0]", l_key_ron, >& "dev$null")
                        }
                        if (imgets.value == "0") {
                            # only warn if it's going to be used
                            if ((l_reject == "ccdclip") || \
                                (l_reject == "crreject")) {
                                printlog ("WARNING - GEMCOMBINE: keyword " // \
                                    l_key_ron // " not found in " // img, \
                                    l_logfile, verbose+)
                                printlog ("Using ron = " // l_ron, l_logfile, \
                                    l_verbose)
                            }
                            l_sron = str (l_ron)
                        } else {
                            l_sron = imgets.value
                        }
                    }
                } else {
                    l_sron = str (l_ron)
                }
                gaineff = gaineff+real (l_sgain)
                roneff = roneff+real (l_sron) **2

                # Make sure exposure time is in the science header
                keypar (img // "[" // sci[i] // "]", l_expname, silent+)
                if (keypar.found == no) {
                    # check the phu
                    keypar (img // "[0]", l_expname, silent+)
                    if (keypar.found == no) {
                        printlog ("ERROR - GEMCOMBINE: " // l_expname // " \
                            not found", l_logfile, l_verbose)
                        goto error
                    } else {
                        gemhedit (img // "[" // sci[i] // "]", l_expname, \
                            keypar.value, "", delete-)
                    }
                }
                # science extension
                print (img // "[" // sci[i] // "]" // sect, >> scilist)

                # DQ extension
                if (imaccess (img // "[" // dq[i] // "]")) {

                    if (debug) {
                        printlog ("___Copying BPM...", l_logfile, verbose+)
                        date ("+%H:%M:%S.%N")
                    }

                    bpmname = tmpdq//"_"//n//".pl"
                    imcopy (img//"["//dq[i]//"]"//sect, bpmname, verbose-)

                    # Store BPM name for later use if propagating DQ values
                    if (l_fl_dqprop) {
                        print (bpmname, >> dqlist)
                    }

                    # Update the input image header with the BPM name
                    gemhedit (img // "[" // sci[i] // "]", "BPM", \
                        bpmname, "", delete-)

                    l_cleanbpm = yes

                    if (debug) {
                        printlog ("___    Done.", l_logfile, verbose+)
                        date ("+%H:%M:%S.%N")
                    }

                } else if (imaccess (l_bpmfile)) {

                    if (i == 1 && !bpmchecked) {
                        # Check whether the BPM file is MEF or not
                        fparse (l_bpmfile)

                        if (fparse.extension != ".pl") {
                            # It's not a mask file it could be a simple fits or
                            # MEF

                            # Check file is MEF
                            gemextn (l_bpmfile, check="exists", \
                                process="expand", \
                                index="", extname="", extversion="", \
                                ikparams="", omit="", replace="", \
                                outfile="dev$null", logfile=l_logfile, \
                                glogpars="", verbose=no)

                            if (gemextn.status != 0) {
                                printlog ("ERROR- GEMCOMBINE: GEMEXTN \
                                    returned a non-zero status. Exiting.", \
                                    l_logfile, verbose+)
                                goto error
                            } else if (gemextn.count != (nsci + 1) \
                                && gemextn.count > 1) {
                                # Check file for DQ planes
                                gemextn (l_bpmfile, check="exists", \
                                    process="expand", index="", \
                                    extname=l_dq_ext, extversion="1-"//nsci, \
                                    ikparams="", omit="", replace="", \
                                    outfile="dev$null", logfile=l_logfile, \
                                    glogpars="", verbose=no)

                                if (gemextn.status != 0) {
                                    printlog ("ERROR- GEMCOMBINE: GEMEXTN \
                                        returned a non-zero status. Exiting.",\
                                        l_logfile, verbose+)
                                    goto error
                                } else if (gemextn.count == nsci) {
                                    bpmismef = yes
                                    bpmextn = l_dq_ext//","
                                } else {
                                    printlog ("WARNING - GEMCOMBINE: "//\
                                        l_bpmfile//" is a MEF file but does \
                                        not\n                     contain \
                                        the same number of extensions \
                                        as the number of input science \
                                        extensions."//\
                                        "\n                      Nor does it \
                                        contain the same number of "//\
                                        l_dq_ext//" extensions as input \
                                        science extenions. Setting \
                                        bpmfile=\"\"", \
                                        l_logfile, verbose+)
                                    l_bpmfile = ""
                                    bpmchecked = yes
                                    goto BPMEND
                                }

                            } else if (gemextn.count != 1 \
                                && gemextn.count == (nsci + 1)) {

                                bpmismef = yes
                                if (imaccess(l_bpmfile//"["//l_dq_ext//",1]")){
                                    bpmextn = l_dq_ext//","
                                } else {
                                    bpmextn = ""
                                }
                            }
                        }
                        bpmchecked = yes
                    }

                    if (bpmismef) {
                        bpm_extn = "["//bpmextn//str(i)//"]"
                    } else {
                        bpm_extn = ""
                    }

                    bpmname = tmpdq//"_"//n//".pl"
                    imcopy (l_bpmfile//bpm_extn//sect, bpmname, verbose-)

                    gemhedit (img // "[" // sci[i] // "]", "BPM", bpmname, \
                        "", delete-)

                    # Store BPM name for later use if propagating DQ values
                    if (l_fl_dqprop) {
                        print (bpmname, >> dqlist)
                    }

                    l_cleanbpm = yes
BPMEND:
                }
            }
            flist = ""

            #combine images
            #special cases for low numbers of images
            if (n == 1) {
                printlog ("ERROR - GEMCOMBINE: only one image to combine.", \
                    l_logfile, verbose+)
                goto error
            }
            if (n <= 5) {
                printlog ("WARNING - GEMCOMBINE: five or less images to \
                    combine.", l_logfile, verbose+)
            }
            if (l_reject == "minmax" && n <= nhigh + nlow) {
                printlog ("ERROR - GEMCOMBINE: Too few images for minmax \
                    parameters (nhigh, nlow).", l_logfile, verbose+)
                goto error
            }
            if (l_fl_vardq == no) {
                combdq = ""
                combsig = ""
            }

            if (debug) {
                printlog ("___IMCOMBINE...", l_logfile, verbose+)
                date ("+%H:%M:%S.%N")
            }

            if (l_verbose) {
                if (debug) print ("calling imcombine")
                if (debug) print (l_scale)
                if (debug) print (l_zero)
                if (debug) print (l_weight)
                imcombine ("@"//scilist, combout, headers = "", bpmasks = "", \
                    rejmasks = "", nrejmasks = combdq, expmasks = "", \
                    sigmas = combsig, logfile = tmplog, combine = l_combine, \
                    reject = l_reject, project = no, outtype = "real", \
                    outlimits = "", offsets = l_offsets, \
                    masktype = l_masktype, maskvalue = l_maskvalue, \
                    blank = 0, scale = l_scale, zero = l_zero, \
                    weight = l_weight, statsec = l_statsec, \
                    expname = l_expname, lthreshold = l_lthreshold, \
                    hthreshold = l_hthreshold, nlow = l_nlow, \
                    nhigh = l_nhigh, nkeep = l_nkeep, mclip = l_mclip, \
                    lsigma = l_lsigma, hsigma = l_hsigma, rdnoise = l_sron, \
                    gain = l_sgain, snoise = l_snoise, sigscale = l_sigscale, \
                    pclip = l_pclip, grow = l_grow)
            } else {
                if (debug) print ("calling imcombine")
                imcombine ("@" // scilist, combout, headers = "", \
                    bpmasks = "", rejmasks = "", nrejmasks = combdq, \
                    expmasks = "", sigmas = combsig, logfile = tmplog, \
                    combine = l_combine, reject = l_reject, project = no, \
                    outtype = "real", outlimits = "", offsets = l_offsets, \
                    masktype = l_masktype, maskvalue = l_maskvalue, \
                    blank = 0, scale = l_scale, zero = l_zero, \
                    weight = l_weight, statsec = l_statsec, \
                    expname = l_expname, lthreshold = l_lthreshold, \
                    hthreshold = l_hthreshold, nlow = l_nlow, \
                    nhigh = l_nhigh, nkeep = l_nkeep, mclip = l_mclip, \
                    lsigma = l_lsigma, hsigma = l_hsigma, rdnoise = l_sron, \
                    gain = l_sgain, snoise = l_snoise, sigscale = l_sigscale, \
                    pclip = l_pclip, grow = l_grow, >& "dev$null")
            }

            if (debug) {
                printlog ("___    Done...", l_logfile, verbose+)
                date ("+%H:%M:%S.%N")
            }

            if (debug) print ("called imcombine")
            if (access (tmplog) ) {
                if (l_logfile != "" && l_logfile != " ") {
                    type (tmplog, >> l_logfile)
                }
                if (l_verbose) {
                    type (tmplog)
                }
                delete (tmplog, verify-, >& "dev$null")
            }

            # New effective gain and readnoise
            if (l_combine == "average") {
                roneff = sqrt (roneff)
            } else {
                gaineff = 2.*gaineff/3.
                roneff = sqrt (2.*roneff/3.)
            }

            # clean headers
            if (l_cleanbpm)
                hedit ("@" // scilist, "BPM", "", add-, addonly-, delete+, \
                    verify-, show-, update+)
            # update gain and readnoise
            imgets (combout, "GAINORIG", >& "dev$null")
            if (imgets.value == "0") {
                gemhedit (combout, "GAINORIG", (real (l_sgain)), \
                    "Input gain", delete-)
            }
            imgets (combout, "RONORIG", >& "dev$null")
            if (imgets.value == "0") {
                gemhedit (combout, "RONORIG", (real (l_sron)), \
                    "Input read-noise", delete-)
            }
            gemhedit (combout, l_key_gain, gaineff, "", delete-)
            gemhedit (combout, l_key_ron, roneff, "", delete-)

            # make variance image by squaring combsig
            #imarith(combsig,"*",combsig,combsig,verbose-)

            if (debug) {
                printlog ("___Creating OUT BPM...", l_logfile, verbose+)
                date ("+%H:%M:%S.%N")
            }

            # Use an "inverse" bad pixel mask (values are the number of pixels
            # used) to get the variance right.
            if (l_fl_vardq) {
                if ("" != l_nrejfile) {
                    imcopy (combdq // ".pl", \
                        l_nrejfile // "[" // i // ",append].fits", verbose-)
                }

                imcalc (combsig // ".fits," // combdq // ".pl", combsig // \
                    ".fits", "(im1**2/(" // n // "-im2))", pixtyp = "real", \
                    verbose-)

                imgets (combsig, "GAINORIG", >& "dev$null")
                if (imgets.value == "0") {
                    gemhedit (combsig, "GAINORIG", (real (l_sgain) **2), \
                        "Input gain", delete-)
                }
                imgets (combsig, "RONORIG", >& "dev$null")
                if (imgets.value == "0") {
                    gemhedit (combsig, "RONORIG", (real (l_sron) **2), \
                        "Input read-noise", delete-)
                }
                gemhedit (combsig, l_key_gain, (gaineff**2), "", delete-)
                gemhedit (combsig, l_key_ron, (roneff**2), "", delete-)

                # DQ file
                imcalc (combdq//".pl", combdq//".fits", "if (im1 == "//n// \
                    ") then 1 else 0", pixtype="ushort", verbose-)
                delete (combdq // ".pl", verify-)

                # propagate DQ values if requested
                if (access (dqlist) && l_fl_dqprop) {
                    printlog (" ", l_logfile, l_verbose)
                    printlog ("GEMCOMBINE: Propagating DQ values", \
                        l_logfile, l_verbose)

                    tmpexprlist = mktemp ("tmpexprlist")

                    dqexp = "im1"

                    # Loop over masks - to create expression - this gets
                    # printed to a file
                    for (j = 2; j <= n; j += 1) {
                        if (strlen(dqexp//" || im"//str(j)) > 64) {
                            print (dqexp, >> tmpexprlist)
                            dqexp = "|| im"//str(j)
                        } else {
                            dqexp = dqexp//" || im"//str(j)
                        }
                    }
                    print (dqexp, >> tmpexprlist)

                    # Add the input masks together
                    addmasks ("@"//dqlist, mdqsum // ".fits", \
                        "@"//tmpexprlist, flags = " ")

                    delete (tmpexprlist, verify-, >& "dev$null")

                    # Only propagate DQ values for pixels that are bad in
                    # output DQ plane from imcombine, if masking. Also, check
                    # for pixels that were all good on the way in but have
                    # somehow been rejected by imcombine. If not
                    # masking propagate all values. - MS
                    # Should really check the individual rejection masks output
                    # by imcomcine as part of the 3-d image, if DQ should be
                    # included for that image. - MS
                    # If masktype is ever anything other than none or goodvalue
                    # - i.e., bitvaues, then this will not work well. - MS
                    if (l_masktype == "none") {
                        dqpr_expr = "(a == 0) ? b : (16 | b)"
                    } else {
                        dqpr_expr = "(a != 0) ? ((b != 0) ? b : 16) : a"
                    }

                    imexpr (dqpr_expr, dqsum//".fits", combdq, \
                        mdqsum//".fits", dims="auto", intype="auto", \
                        outtype="ushort", refim="auto", bwidth=0, \
                        btype="nearest", bpixval=0., rangecheck=yes,\
                        verbose=no)

                    imdelete (combdq//".fits", verify-, >& "dev$null")
                    imrename (dqsum//".fits", combdq//".fits", verbose-)
                    imdelete (dqsum//","//mdqsum, \
                        verify=no, >& "dev$null")

                    delete ("@"//dqlist, verify-, >& "dev$null")
                    delete (dqlist, verify-, >& "dev$null")
                }
                # Headers
                gemhedit (combsig // "[0],", "BPM", "", "", delete+)
                gemhedit (combdq//"[0]", "BPM", "", "", delete+)
            }
            # update headers
            gemhedit (combout, "BPM", "", "", delete+)
            if (l_title != "" && l_title != " ") {
                gemhedit (combout, "i_title", l_title, "", delete-)
            }

            if (debug) {
                printlog ("___    Done.", l_logfile, verbose+)
                date ("+%H:%M:%S.%N")
            }

            if (debug) {
                printlog ("___Creating OUTPUT...", l_logfile, verbose+)
                date ("+%H:%M:%S.%N")
            }

            #make output MEF file
            files (sci[i], sort-) | scan (tsci)
            files (var[i], sort-) | scan (tvar)
            files (dq[i], sort-) | scan (tdq)
            if (i == 1) {
                if (l_fl_vardq) {
                    wmef (combout // "," // combsig // "," // combdq, \
                        l_output, extnames = tsci // "," // tvar // "," // \
                        tdq, phu = phuimg, verbose-, >& "dev$null")

                    if (useextver) {
                        gemhedit (l_output // "[1]," // l_output // "[2]," // \
                            l_output // "[3]", "EXTVER", i, "", delete-)
                    }
                } else {
                    if (debug) {
                        print ("calling wmef")
                        wmef (combout, l_output, extnames = tsci, phu = phuimg,
                            verbose+)
                        print ("wmef completed")
                    } else {
                        wmef (combout, l_output, extnames = tsci, phu = phuimg,
                            verbose-, >& "dev$null")
                    }
                    if (useextver) {
                        gemhedit (l_output // "[1]", "EXTVER", i, "", delete-)
                    }
                }
                if (wmef.status > 0) {
                    goto error
                }
                imgets (phuimg // "[0]", l_expname, >& "dev$null")
                exphu = imgets.value
                imgets (combout, l_expname, >& "dev$null")
                if (real (exphu) != real (imgets.value) ) {
                    gemhedit (l_output // "[0]", l_expname, imgets.value, "", \
                        delete-)
                }
                # update gain and readnoise in PHU?
                imgets (phuimg // "[0]", "GAIN", >& "dev$null")
                if (imgets.value != "0") {
                    imgets (phuimg // "[0]", "GAINORIG", >& "dev$null")
                    if (imgets.value == "0") {
                        gemhedit (l_output // "[0]", "GAINORIG", \
                            (real (l_sgain)), "Input gain", delete-)
                    }
                    imgets (phuimg // "[0]", "RONORIG", >& "dev$null")
                    if (imgets.value == "0") {
                        gemhedit (l_output // "[0]", "RONORIG", \
                            (real (l_sron)), "Input read-noise", delete-)
                    }
                    gemhedit (l_output // "[0]", l_key_gain, gaineff, "", \
                        delete-)
                    gemhedit (l_output // "[0]", l_key_ron, roneff, "", \
                        delete-)
                }
            } else {
                if (l_fl_vardq) {
                    imcopy (combdq, dqfits // ".fits", verbose-)
                    fxinsert (combout // ".fits," // combsig // ".fits," // \
                        dqfits // ".fits", l_output // "[" // ( (i-1) *3) // \
                        "]", "0", verbose-, >& "dev$null")
                    imgets (l_output // "[0]", "NEXTEND")
                    if (debug) print (imgets.value)
                    nextnd = int (imgets.value) +3
                    gemhedit (l_output // "[0]", "NEXTEND", nextnd, "", \
                        delete-)
                    gemhedit (l_output // "[" // ( (i-1) *3+1) // "]", \
                        "EXTNAME", tsci, "", delete-)
                    gemhedit (l_output // "[" // ( (i-1) *3+2) // "]", \
                        "EXTNAME", tvar, "", delete-)
                    gemhedit (l_output // "[" // ( (i-1) *3+3) // "]", \
                        "EXTNAME", tdq, "", delete-)
                    if (useextver) {
                        gemhedit (l_output // "[" // ( (i-1) *3+1) // "]", \
                            "EXTVER", i, "", delete-)
                        gemhedit (l_output // "[" // ( (i-1) *3+2) // "]", \
                            "EXTVER", i, "", delete-)
                        gemhedit (l_output // "[" // ( (i-1) *3+3) // "]", \
                            "EXTVER", i, "", delete-)
                    }
                    imdelete (dqfits // ".fits", verify-)
                } else {
                    if (debug) print ("calling fxinsert")
                    fxinsert (combout // ".fits", l_output // "[" // (i-1) // \
                        "]", "0", verbose-, >& "dev$null")
                    if (debug) print ("fxinsert completed")
                    imgets (l_output // "[0]", "NEXTEND")
                    if (debug) print (imgets.value)
                    nextnd = int (imgets.value) +1
                    gemhedit (l_output // "[0]", "NEXTEND", nextnd, "", \
                        delete-)
                    gemhedit (l_output // "[" // i // "]", "EXTNAME", tsci, \
                        "", delete-)
                    if (useextver) {
                        gemhedit (l_output // "[" // i // "]", "EXTVER", i, \
                            "", delete-)
                    }
                }
            }
            if (debug) {
                printlog ("___    Done.", l_logfile, verbose+)
                date ("+%H:%M:%S.%N")
            }

            imdelete (combout // "," // combsig // "," // combdq // "," // \
                dqfits, verify-, >& "dev$null")
            imdelete (tmpdq // "*.pl", verify-, >& "dev$null")
            delete (scilist, verify-, >& "dev$null")
        }
        # attach MDF
        if (usemdf) {
            if (l_verbose) {
                printlog ("GEMCOMBINE: Attaching the MDF", l_logfile, \
                    l_verbose)
            }
            fxinsert (mdf // ".fits", l_output // "[" // (i-1) // "]", "1", \
                verbose-, >& "dev$null")
            imgets (l_output // "[0]", "NEXTEND")
            if (debug) print (imgets.value)
            nextnd = int (imgets.value) +1
            gemhedit (l_output // "[0]", "NEXTEND", nextnd, "", delete-)
        }

        # Update PHU
        gemdate ()
        gemhedit (l_output // "[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
        gemhedit (l_output // "[0]", "GEMCOMB", gemdate.outdate,
            "UT Time stamp for GEMCOMBINE", delete-)
        # clean up
        goto clean

error:
        {
            status = 1
            goto clean
        }

clean:
        {
            delete (filelist // "," // scilist, verify-, >& "dev$null")
            imdelete (combout // "," // combsig // "," // combdq // ".fits", \
                verify-, >& "dev$null")
            imdelete (dqsum // "," // mdf // ".fits", verify-, >& "dev$null")
            imdelete (tmpdq // "*.pl," // combdq // ".pl", verify-, \
                >& "dev$null")
            delete (sextlist // "," // vextlist // "," // dextlist // "," // \
                extlist // "," // scale_list // "," // dqlist, \
                verify-, >& "dev$null")
            hedit.show = s_verbose
            # close log file
            printlog (" ", l_logfile, l_verbose)

            date | scan (sdate)
            printlog ("GEMCOMBINE -- Finished: "//sdate, \
                l_logfile, verbose = l_verbose)
            printlog (" ", l_logfile, l_verbose)
            if (status == 0) {
                printlog ("GEMCOMBINE -- Exit status: GOOD", \
                    l_logfile, l_verbose)
            } else {
                printlog ("GEMCOMBINE -- Exit status: ERROR", \
                    l_logfile, l_verbose)
            }
            printlog ("-------------------------------------------------------\
                -------------------------", l_logfile, l_verbose)
            flpr
        }

end
