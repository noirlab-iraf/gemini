# Copyright(c) 2010-2013 Association of Universities for Research in Astronomy, Inc.

procedure gafastsky(inimages)

# Make a quick-and-dirty median sky image for GSAOI images.

char    inimages     {prompt="GSAOI sky images to combine"}
char    outimages    {"",prompt="Output sky images"}
#char    outtitle     {"default",prompt="Title for output image"}
char    combine      {"default",enum="default|median|average",prompt="Type of combine operation (default|median|average)"}
char    reject       {"minmax",enum="none|minmax",prompt="Type of rejection (none|minmax)"}
int     nlow         {0,min=0,prompt="minmax: Number of low pixels to reject"}
int     nhigh        {1,min=0,prompt="minmax: Number of high pixels to reject"}
char    statsec      {"[5%]",prompt="Statistics section"}
char    masktype     {"goodvalue",enum="none|goodvalue",prompt="Bad Pixel Mask type (none|goodvalue)"}
real    maskvalue    {0.,prompt="Good pixel value in the BPM"}
char    badpix       {"gsaoi$data/gsaoibpm_high_full.fits",prompt="Static Bad Pixel Mask"}
char    key_exptime  {"EXPTIME",prompt="Keyword for exposure time"}
char    key_ron      {"RDNOISE",prompt="Header keyword for read noise (e-)"}
char    key_gain     {"GAIN",prompt="Header keyword for gain (e-/ADU)"}
bool    fl_vardq     {no,prompt="Create variance and data quality frames in output?"}
bool    fl_dqprop    {no,prompt="Propagate data quality information?"}
char    sci_ext      {"SCI",prompt="Name for science extensions"}
char    var_ext      {"VAR",prompt="Name for variance extensions"}
char    dq_ext       {"DQ",prompt="Name for data quality extensions"}
char    logfile      {"",prompt="Name of log file"}
bool    verbose      {yes,prompt="Verbose?"}
int     status       {0,prompt="Exit status (0=good)"}
struct  *scanfile    {"",prompt="Internal use only"}
struct  *scanfile2   {"",prompt="Internal use only"}
struct  *scanfile3   {"",prompt="Internal use only"}

######################################################

begin

    ########
    # Declare and set local variables; set default values; initiate temporary
    # files

    ####
    # Variable declaration
    char outgaimchk_list, l_inimages, l_outimages, l_combine, l_reject
    char l_logfile, l_key_exptime, l_key_filter, l_statsec
    char l_sci_ext, l_var_ext, l_dq_ext, l_key_ron, l_key_gain
    char l_rootname, l_rawpath, l_obs_allowed, l_obj_allowed, l_gaprep_pref
    char l_garedux_pref, l_key_allowed, tmpfile, l_redirect, filelist
    char tmplist, tmpconflist, conflist, tmpoutlist, tmpimgout
    char one_ron, one_gain, separator, orig_1_ron, orig_1_gain
    char tmpout_gaimchk, inimg, inphu, insci, inname, metaconf
    char l_separator, tmp_skipped, pr_string, outphu
    char l_title, l_masktype, l_scale, l_zero, l_weight
    char l_mclip, l_snoise, l_badpix, l_nrejfile, l_orig_combine
    char tmpbpm, bpmgemcomb, area_bpm, refname

    int total_num_created, num_created, nconf, i, nsci, lcount, n_out, nn
    int last_ron_char, last_gain_char, l_test, l_nkeep, l_nlow, l_nhigh, nfiles

    real l_maskvalue, l_lthreshold, l_hthreshold, l_lsigma, l_hsigma
    real l_sigscale, l_pclip, l_grow, l_expone, l_ron, l_gain

    bool l_verbose, l_fl_dqprop, l_fl_vardq, debug

    # Set local variables
    l_inimages = inimages
    l_outimages = outimages
    l_key_exptime = key_exptime
#    l_outtitle = outtitle
    l_orig_combine = combine
    l_reject = reject
    l_nlow = nlow
    l_nhigh = nhigh
    l_masktype = masktype
    l_maskvalue = maskvalue
    l_badpix = badpix
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_var_ext = var_ext
    l_fl_dqprop = fl_dqprop
    l_fl_vardq = fl_vardq
    l_key_ron = key_ron
    l_key_gain = key_gain
    l_statsec = statsec
    l_logfile = logfile
    l_verbose = verbose

    # Make temporary files
    outgaimchk_list = mktemp("tmpoutgaimchk")
    filelist = mktemp ("tmpfilelist")
    tmplist = mktemp("tmplist")
    tmpconflist = mktemp ("tmpconflist")
    tmpoutlist = mktemp ("tmpoutlist")
    tmpout_gaimchk = mktemp ("tmpoutchk")
    tmp_skipped = mktemp ("tmpskip")
    tmpbpm = mktemp ("tmpbpm")

    # Set default values for gaimchk
    l_rootname = ""
    l_rawpath = ""
    l_obs_allowed = "OBJECT"
    l_obj_allowed = ""
    l_gaprep_pref = "g"
    l_garedux_pref = "r"

    # Set other default values
    status=0
    debug = no

    if (l_verbose) {
        l_redirect = "STDOUT"
    } else {
        l_redirect = "dev$null"
    }

    #------------------------------------------------------------------------
    # Fun starts here

    # Keep gemcombine parameters from changing by outside world
    cache("gemcombine")

    # Test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = gsaoi.logfile
       if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
          l_logfile = "gsaoi.log"
          printlog ("WARNING - GAFASTSKY: Both gaprepare.logfile and \
              gsaoi.logfile fields are empty", l_logfile, l_verbose)
          printlog ("                     Using default file gsaoi.log", \
              l_logfile, l_verbose)
       }
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)

    # Print start time
    printlog ("", l_logfile, l_verbose)
    gemdate (zone="local")
    printlog ("GAFASTSKY -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    # Record all input parameters
    printlog ("GAFASTSKY: Input parameters...", l_logfile, l_verbose)
    printlog ("    inimages    = "//l_inimages, l_logfile, l_verbose)
    printlog ("    outimage    = "//l_outimages, l_logfile, l_verbose)
#    printlog ("    outtitile   = "//l_outtitle, l_logfile, l_verbose)
    printlog ("    key_exptime = "//l_key_exptime, l_logfile, l_verbose)
    printlog ("    combine     = "//l_orig_combine, l_logfile, l_verbose)
    printlog ("    statsec     = "//l_statsec, l_logfile, l_verbose)
    printlog ("    reject      = "//l_reject, l_logfile, l_verbose)
    printlog ("    nlow        = "//l_nlow, l_logfile, l_verbose)
    printlog ("    nhigh       = "//l_nhigh, l_logfile, l_verbose)
    printlog ("    badpix      = "//l_badpix, l_logfile, l_verbose)
    printlog ("    sci_ext     = "//l_sci_ext, l_logfile, l_verbose)
    printlog ("    var_ext     = "//l_var_ext, l_logfile, l_verbose)
    printlog ("    dq_ext      = "//l_dq_ext, l_logfile, l_verbose)
    printlog ("    key_ron     = "//l_key_ron, l_logfile, l_verbose)
    printlog ("    key_gain    = "//l_key_gain, l_logfile, l_verbose)
    printlog ("    fl_vardq    = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("    fl_dqprop   = "//l_fl_dqprop, l_logfile, l_verbose)
    printlog ("    logfile     = "//l_logfile, l_logfile, l_verbose)
    printlog ("    verbose     = "//l_verbose, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    # Check the l_badpix value
    if ((l_badpix == "") || (stridx(" ", l_badpix) > 0)) {
        # This assumes that the bpm covers all of the areas readout in each
        # detector in each image
        l_badpix = "none"
    } else if (!imaccess(l_badpix)) {
        printlog ("ERROR - GAFASTSKY: cannot access bpm "//l_badpix//\
            ". Exiting", l_logfile, verbose+)
    } else {
        # Check it hasn't been mosaiced
        keypar (l_badpix//"[0]", "GAMOSAIC", silent+)
        if (keypar.found) {
            printlog ("ERROR - GAFASTSKY: "//l_badpix//" has been mosaiced."//\
                " Please supply a static BPM that has not been mosaiced.", \
                l_logfile, verbose)
            goto crash
        }
    }

    # Verify that the extension names are not empty, otherwise exit gracefully
    if (l_sci_ext == "" || stridx(" ",l_sci_ext) > 0) {
        printlog("ERROR - GAFASTSKY: extension name sci_ext is missing.", \
            l_logfile,verbose+)
        status=1
        goto clean
    }

    if (l_fl_vardq) {
        if (l_dq_ext == "" || (stridx(" ",l_dq_ext) > 0) || \
            l_var_ext == "" || (stridx(" ",l_var_ext) > 0)) {

            printlog("WARNING - GAFASTSKY: var_ext or dq_ext have not been \
                set.\n                     Output image will not have VAR or \
                DQ planes.", l_logfile,verbose+)
            l_fl_vardq=no
        }
    }

    if (l_fl_dqprop && !l_fl_vardq) {
        printlog ("WARNING - GAFASTSKY: Cannot propagate DQ information with \
            fl_vardq=no. Setting fl_dqprop=no", l_logfile, verbose+)
            l_fl_dqprop = no
    }
    printlog ("GAFASTSKY: Calling GAIMCHK to check inputs...", \
        l_logfile, l_verbose)

    if (debug) gemdate (verbose+)

    # Uncomment if having to allow reduced files...
    # (and add another loop below) - MS
    # They ashould all be in the same state... can't mix reduce and prepared
    # This could get ugly...
    #    l_key_allowed="GAREDUCE"
    l_key_allowed=""

    # Call gaimchk to perform input checks
    gaimchk (inimages=l_inimages, rawpath=l_rawpath, rootname=l_rootname, \
        obstype_allowed=l_obs_allowed, object_allowed=l_obj_allowed, \
        key_allowed=l_key_allowed, key_forbidden="GAFSTSKY", \
        key_exists="GAPREPAR", fl_prep_check=yes, gaprep_pref=l_gaprep_pref, \
        fl_redux_check=no, garedux_pref=l_garedux_pref, fl_fail=no, \
        fl_vardq_check=l_fl_vardq, sci_ext=l_sci_ext, var_ext=l_var_ext, \
        dq_ext=l_dq_ext, outlist=outgaimchk_list, logfile=l_logfile, \
        verbose=no)

    if (gaimchk.status != 0) {
        printlog ("ERROR - GAFASTSKY: GAIMCHK returned a non-zero status. \
            Exiting.", l_logfile, verbose+)
        goto crash
    } else if (!access(outgaimchk_list)) {
        printlog ("ERROR - GAFASTSKY: Cannot access output list from GAIMCHK",\
            l_logfile, verbose+)
        goto crash
    }

    # Files that have been processed by GALFAT
    tmpfile = ""
    match ("tmpGAFSTSKY", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        printlog ("WARNING - GAFASTSKY: The following files have already been \
            processed by GAFASTSKY. Ignoring these files...", \
            l_logfile, verbose+)
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no) | \
            tee (l_logfile, out_type="text", append+)
    }

    # Files that are not prepared
    tmpfile = ""
    match ("tmpNotFound", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        printlog ("WARNING - GAFASTSKY: The following files have not been \
            processed by GAPREPRE. Ignoring these files...", \
            l_logfile, verbose+)
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no) | \
            tee (l_logfile, out_type="text", append+)
    }

    # Files that have been prepared
    tmpfile = ""
    match ("tmpGAPREPAR", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
                    print_file_names=no, >> tmplist)
    }

    lcount = 0
    if (access(tmplist)) {
        count (tmplist) | scan (lcount)
    }

    if (lcount == 0) {
        printlog ("ERROR - GAFASTSKY: No unique metaconfigurations values "//\
            "found.\n                   Exiting", l_logfile, verbose+)
        goto crash
    }

    if (debug) gemdate (verbose+)
    # Call gaimchk to perform output checks
    n_out = 0
    if (l_outimages != "" || stridx(" ",l_outimages) > 0) {

        printlog ("GAFASTSKY: Calling GAIMCHK to check outputs...", \
            l_logfile, l_verbose)

        if (debug) gemdate (verbose+)

        gaimchk (inimages=l_outimages, rawpath="", rootname="", \
            obstype_allowed="", object_allowed="", key_allowed="", \
            key_forbidden="", key_exists="",  fl_prep_check=no, \
            gaprep_pref=l_gaprep_pref, fl_redux_check=no, \
            garedux_pref=l_garedux_pref, fl_fail=no, \
            fl_out_check=yes, fl_vardq_check=no, sci_ext=l_sci_ext, \
            var_ext=l_var_ext, dq_ext=l_dq_ext, outlist=tmpout_gaimchk, \
            logfile=l_logfile, verbose=no)

        if (gaimchk.status != 0) {
            printlog ("ERROR - GAFASTSKY: GAIMCHK returned a non-zero status. \
                Exiting.", l_logfile, verbose+)
            goto crash
        } else if (!access(tmpout_gaimchk)) {
            printlog ("ERROR - GAFASTSKY: Cannot access output list from \
                GAIMCHK",\
                l_logfile, verbose+)
            goto crash
        }

        # Files that need to be prepared
        tmpoutlist = ""
        match ("tmpNotFound", tmpout_gaimchk, stop=no, print_file_name=no, \
            metacharacter=no) | scan (tmpoutlist)

        if (tmpoutlist == "") {
            printlog ("ERROR - GAFASTSKY: Output list does not exist. \
                Exiting", l_logfile, verbose+)
            goto crash
        } else {
            count (tmpoutlist) | scan (n_out)
            if (n_out == 0) {
                printlog ("ERROR - GAFASTSKY: Output list is empty. Exiting", \
                    l_logfile, verbose+)
                goto crash
            } else {
                scanfile3 = tmpoutlist
            }
        }

        if (debug) gemdate (verbose+)
    } else {
        l_outimages = "default"
    }

    # Loop over nextlist to obtain the METACONFIG keyword. This contains
    # information that will group files such that a) they are in the same
    # configuration and b) that they have been reduced in the same way.
    printlog ("\nGAFASTSKY: Accessing METACONF keyword values", \
        l_logfile, l_verbose)

    if (debug) gemdate (verbose+)

    scanfile = tmplist
    while (fscan(scanfile, inimg) != EOF) {

        inphu = inimg//"[0]"

        # Double check that you can access the file...
        if (!imaccess(inimg)) {
            printlog ("WARNING - GAFASTSKY: Image "//inimg//\
                " does not exist. Skipping this image.", l_logfile, verbose)

        } else {
            # Read the METACONFIG keyword
            keypar (inphu, "METACONF", silent+)
            if (!keypar.found) {
                printlog ("WARNING - GAFASTSKY: METACONF keyword not \
                    found in "//inimg//" Skipping this image.", \
                    l_logfile, verbose+)

            } else {
                metaconf = keypar.value

                # Check it's valid
                if (metaconf == "UNKNOWN" || metaconf == "" || \
                    (stridx(" ",metaconf) > 0)) {
                    printlog ("WARNING - GAFASTSKY: METACONF keyword is "//\
                        "\"UNKNOWN\" or bad in"//inimg//".\n"//\
                        "                     Skipping this image.", \
                        l_logfile, verbose+)

                } else {

                    # Check it's an OBJECT!
                    if (strstr("OBJ",metaconf) == 0) {
                        printlog ("WARNING - GAFASTSKY: Image \""//inimg//\
                            "\" is a calibration image. Skipping this image.",\
                            l_logfile, verbose+)
                    } else {
                        print (inimg//" "//metaconf, >> filelist)
                    }
                }
            }
        }

    } # End of loop to read METACONF keyword

    delete (tmplist, verify-, >& "dev$null")

    # Check if there are any images in the newly created filelist
    lcount = 0
    if (access(filelist)) {
        count (filelist) | scan (lcount)
    }

    if (lcount == 0) {
        printlog ("ERROR - GAFASTSKY: No images with a valid METACONF keyword \
            value. Please try again.", l_logfile, verbose+)
        goto crash
    }

    if (debug) gemdate (verbose+)

    # Check metaconf keywords
    # Determine the unique METACONF values for the file that made it
    printlog ("GAFASTSKY: Determining unique METACONFIG values", \
        l_logfile, l_verbose)

    # Make sure the reverse_sort is true to have closed GCAL flats after open
    # ones - MS
    fields (filelist, fields=2, lines="1-", quit_if_miss=no, \
        print_file_names=no) | sort ("STDIN", column=0, ignore_white=no, \
        numeric_sort=no, reverse_sort=yes) | unique ("STDIN", > tmpconflist)

    # Check if there are any images in the newly created tmpconflist
    lcount = 0
    if (access(tmpconflist)) {
        count (tmpconflist) | scan (lcount)
    }

    if (lcount == 0) {
        printlog ("ERROR - GAFASTSKY: No unique METACONF keywords found. \
            Please try again.", l_logfile, verbose+)
        goto crash
    } else if (l_outimages != "default" && lcount != n_out) {
        printlog ("ERROR - GAFASTSKY: Number of supplied output image \
            names does"//\
            "\n                   not match the number of unique \
            METACONFIGURATION keywaord values. Exiting", l_logfile, verbose+)
        goto crash
    }

    ####
    # Create the FLATs...

    total_num_created = 0

    # This becomes a maze of lists! Reader beware.

    if (debug) gemdate (verbose+)

    # Create separate lists for each unique METACONF keyword
    printlog ("GAFASTSKY: Creating METACONFIG lists and processing images...",\
        l_logfile, l_verbose)

    nconf = 1
    scanfile = tmpconflist
    while (fscan(scanfile, metaconf) != EOF) {

        num_created = 0

        printlog ("\n----\nGAFASTSKY: Processing configuration: "//\
            metaconf//"\n", l_logfile, l_verbose)

        if (debug) gemdate (verbose+)

        # Set up temporary configuration list
        conflist = tmplist//"_conf"//nconf//".lis"

        # Extract the files that have the same METACONF keyword value
        match (metaconf, filelist, stop=no, print_file_name=no, \
            metacharacter=no) | fields ("STDIN", fields=1, lines="1-", \
            quit_if_miss=no, print_file_names=no, > conflist)


        # Check the number of inputs and set up input into gemcombine
        # appropriately - if only one file skip
        lcount = 0
        if (access(conflist)) {
            count (conflist) | scan (lcount)
        } else {
            printlog ("ERROR - GAFASTSKY: Cannot access input \
                configuration list", l_logfile, verbose+)
            goto crash
        }
        if (lcount == 0) {
            printlog ("ERROR - GAFASTSKY: Configuration list for "//metaconf//\
                " is empty. Exiting", l_logfile, verbose+)
            goto crash
        } else if (lcount == 1) {
            printlog ("WARING - GAFASTSKY: Only one file to combine for "//\
                metaconf//". Skipping this configuration.",_logfile, verbose+)
            goto NEXT_CONFIGURATION
        }
        nfiles = lcount

        printlog ("GAFASTSKY: Using input files...\n", l_logfile, l_verbose)
        type (conflist) | \
            tee (l_logfile, out_type="text", append+, > l_redirect)

        if (debug) gemdate (verbose+)

        # Loop over inputs
        l_expone = 0.0
        nn = 0
        scanfile2 = conflist
        while (fscan(scanfile2, inimg) != EOF) {

            inphu = inimg//"[0]"
            fparse (inimg, verbose-)
            inname = fparse.root//fparse.extension

            if (nn == 0) {
                l_title = ""
                # Read the OBJECT for the first frame
                keypar (inphu, "OBJECT", silent+)
                if (keypar.found) {
                    l_title = keypar.value
                }
            }

            if (l_outimages == "default" && nn == 0) {
                tmpimgout = fparse.root//"_fsky.fits"
                if (imaccess(tmpimgout)) {
                    printlog ("\nERROR - GAFASTSKY: "//tmpimgout//\
                        " already exists. Exiting.", l_logfile, verbose+)
                    goto crash
                }
            } else if (l_outimages != "default" && nn == 0) {
                tmpimgout = ""
                l_test = fscan (scanfile3, tmpimgout)
                if (tmpimgout == "") {
                    printlog ("\nERROR - Cannot read output file name from "//\
                        tmpoutimg//". Exiting.", l_logfile, verbose+)
                    goto crash
                }
            }

            # check exposure time - Though this should already be done by the
            # unique metaconfig extraction
            keypar (inphu, l_key_exptime, silent+)
            if (!keypar.found) {
                printlog("ERROR - GAFASTSKY: "//l_key_exptime//" keyword not \
                    found in "//inname, l_logfile, verbose+)
                goto crash
            } else if (nn == 0) {
                l_expone = real(keypar.value)
            } else if (abs(real(keypar.value) - l_expone) > 0.1) {
                printlog("WARNING - GAFASTSKY: exposure times are \
                    significantly different. Continuing.", l_logfile,verbose+)
            }

            keypar (inphu, "NSCIEXT", silent+)
            if (!keypar.found) {
                printlog("ERROR - GAFASTSKY: NSCIEXT keyword not \
                    found in "//inname, l_logfile, verbose+)
                goto crash
            } else {
                nsci = int (keypar.value)
            }

            # Check readnoise and gain values are consistent
            for (i = 1; i <= nsci; i += 1) {

                insci = "["//l_sci_ext//","//i//"]"
                inname = inname//insci
                insci = inimg//insci

                # Check the dimensions of the image against the requested
                # statsec
                gadimschk (insci, section=l_statsec, chkimage="", \
                    key_check="CCDSEC", logfile=l_logfile, verbose-)

                if (gadimschk.status != 0) {
                    printlog ("ERROR - GAFASTSKY: GADIMSCHK returned a \
                        non-zero status. Exiting.", l_logfile, verbose-)
                    goto crash
                }

                if (nn != 0) {
                    if (i == 1) {
                        one_ron = orig_1_ron
                        one_gain = orig_1_gain
                    }

                    if (i < nsci) {
                        last_ron_char = stridx(",",one_ron) - 1
                        last_gain_char = stridx(",",one_gain) - 1
                    } else {
                        last_ron_char = strlen(one_ron)
                        last_gain_char = strlen(one_gain)
                    }

                    l_ron = real(substr(one_ron,1,last_ron_char))
                    l_gain = real(substr(one_gain,1,last_gain_char))

                    if (i < nsci) {
                        one_ron = substr(one_ron,last_ron_char+2,\
                            strlen(one_ron))
                        one_gain = substr(one_gain,last_gain_char+2,\
                            strlen(one_gain))
                    }

                } else if (i == 1) {
                    l_separator = ""
                    one_ron = ""
                    one_gain = ""
                } else {
                    l_separator = ","
                }

                # check read noise
                keypar (insci, l_key_ron, silent+)
                if (!keypar.found) {
                    printlog("ERROR - GAFASTSKY: "//l_key_ron//" keyword not \
                        found in "//inname, l_logfile, verbose+)
                    goto crash
                } else if (nn == 0) {
                    one_ron = one_ron//l_separator//keypar.value

                    if (i == nsci) {
                        orig_1_ron = one_ron
                    }
                } else if (abs(real(keypar.value) - l_ron) > 1.) {
                    printlog("WARNING - GAFASTSKY: read noise values are \
                        different.  Continuing, but the", l_logfile, verbose+)
                    printlog("                     read noise in the output \
                        header will be wrong.", l_logfile, verbose+)
                }

                keypar (insci, l_key_gain, silent+)
                if (!keypar.found) {
                    printlog("ERROR - GAFASTSKY: "//l_key_gain//" keyword not \
                        found in "//inname, l_logfile, verbose+)
                    goto crash
                } else if (nn == 0) {
                    one_gain = one_gain//l_separator//keypar.value

                    if (i == nsci) {
                        orig_1_gain = one_gain
                    }
                } else if (abs(real(keypar.value)- l_gain) > 0.5) {
                    printlog("WARNING - GAFASTSKY: read noise values are \
                        different.  Continuing, but the", l_logfile, verbose+)
                    printlog("                     read noise in the output \
                        header will be wrong.", l_logfile, verbose+)
                }
            } # End of for loop over sci

            nn += 1
        } # End of while of current configurations input list
        scanfile2 = ""

        if (debug) gemdate (verbose+)


        if (nfiles != nn) {
            printlog ("ERROR - GAFASTSKY: Mismatch of number of files in \
                current configuration list", l_logfile, verbose+)
            goto crash
        }

        l_title = l_title//" FSKY IMAGE from gemini.gsaoi.gafastsky"

        # Gemcombine parameter values
        ##M Open these up to the user?
        l_scale = "median"
        l_zero = "none"
        l_weight = "none"
        l_nkeep = 0
        l_mclip = yes
        l_lthreshold = INDEF
        l_hthreshold = INDEF
        l_lsigma = 3.0
        l_hsigma = 3.0
        l_snoise = "0.0"
        l_sigscale = 0.1
        l_pclip = -0.5
        l_grow = 0.
        l_nrejfile = ""

        if (l_orig_combine == "default") {

            if (l_reject != "minmax") {
                printlog ("\nWARNING - GAFASTSKY: combine is \"default\" \
                    setting reject=\"minmax\".", l_logfile, verbose+)
                l_reject = "minmax"
            }

            l_combine = "median"
            l_nhigh = 1
            l_nlow = 1

            if (nfiles < 5) {
                l_combine = "average"
                l_nlow -= 1
                printlog("\nWARNING - GAFASTSKY: Averaging 4 or fewer images \
                    with no low pixels rejected"//\
                    "\n                     and only 1 high pixel \
                    rejected", l_logfile, verbose+)
            } else if (nfiles >= 8) {
                l_nhigh += 1
            }

        } else {

            if (nfiles < 5) {
                printlog("\nWARNING - GAFASTSKY: //"//l_combine//" combining \
                    4 or fewer images.", l_logfile,verbose+)
                if (l_reject == "none") {
                    printlog("                      with no pixels rejected.",\
                        l_logfile, verbose+)
                } else {
                    printlog("                      with "//l_nlow//\
                        " low and "//l_nhigh//" high pixels rejected.",
                        l_logfile, verbose+)
                }
            } # end if(nfiles < 5)

            if ((nfiles <= (l_nlow + l_nhigh)) && (l_reject == "minmax")) {
                printlog("\nERROR - GAFASTSKY: Cannot reject more pixels than \
                    the number of images.", l_logfile,verbose+)
                status=1
                goto crash
            }
        } # end not-default section

        pr_string = strupr(substr(l_combine,1,1))//\
            strlwr(substr(l_combine,2,strlen(l_combine)))

        printlog ("\nGAFASTSKY: "//pr_string//" combining sky frames...", \
            l_logfile, l_verbose)

        if (l_reject == "minmax") {
            pr_string = " with "//l_nlow//" low and "//l_nhigh//\
                " high values rejected"
        } else {
            pr_string = ""
        }

        printlog ("           Rejection type is "//l_reject//pr_string, \
            l_logfile, l_verbose)

        bpmgemcomb = ""
        # Check the BPM dimensions against the first image in the list
        if (l_badpix != "none" && l_masktype != "none") {

            # Use the first file in the list - the rest should all be the same
            # dimensions - MS
            refname = ""
            head (conflist, nlines=1) | scan (refname)

            keypar (refname//"[0]", "NSCIEXT", silent+)
            nsci = int(keypar.value)
            for (i = 1; i <= nsci; i += 1) {
                gadimschk (refname//"["//l_sci_ext//","//i//"]", \
                    section="", \
                    chkimage=l_badpix//"["//l_dq_ext//","//i//"]", \
                    key_check="CCDSEC", logfile=l_logfile, verbose-)

                if (gadimschk.status != 0) {
                    printlog ("ERROR - GAFLAT: GADIMSCHK returned a non-zero \
                        status. Exiting.", l_logfile, verbose-)
                    goto crash
                } else {
                    # Output from gadimschk has always a sections appended
                    area_bpm = substr(gadimschk.out_chkimage,\
                        strldx("[",gadimschk.out_chkimage),\
                        strlen(gadimschk.out_chkimage)) == "[*,*]"

                    if (area_bpm == "[*,*]") {
                        bpmgemcomb = l_badpix
                    } else if (i == 1) {
                        imcopy (l_badpix//"[0]", tmpbpm, verbose-)
                        imcopy (gadimschk.out_chkimage, tmpbpm//"[append]", \
                            verbose-)
                        bpmgemcomb = tmpbpm
                    } else if (!imaccess(tmpbpm) && i > 1) {
                        printlog ("ERROR - GAFLAT: The array dimensions are \
                            different sizes for input image. Exiting.", \
                            l_logfile, verbose+)
                    } else {
                        imcopy (gadimschk.out_chkimage, tmpbpm//"[append]", \
                            verbose-)
                    }
                }
            } # End of loop over sci extensions
            refname = ""
        } # End of bpm dimensions check

        gemcombine ("@"//conflist, output=tmpimgout, title=l_title, \
            combine=l_combine, reject=l_reject, offsets="none", \
            masktype=l_masktype, maskvalue=l_maskvalue, scale=l_scale, \
            zero=l_zero, weight=l_weight, statsec=l_statsec, \
            expname=l_key_exptime, lthreshold=l_lthreshold, \
            hthreshold=l_hthreshold, nlow=l_nlow, nhigh=l_nhigh, \
            nkeep=l_nkeep, mclip=yes, lsigma=l_lsigma, hsigma=l_hsigma, \
            key_ron=l_key_ron, key_gain=l_key_gain, ron=0., gain=1., \
            snoise=l_snoise, sigscale=l_sigscale, pclip=l_pclip, \
            grow=l_grow, bpmfile=bpmgemcomb, nrejfile=l_nrejfile, \
            sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
            fl_vardq=l_fl_vardq, logfile=l_logfile, \
            fl_dqprop=l_fl_dqprop, verbose=l_verbose)

        # Check the output status of gemcombine
        if (gemcombine.status != 0) {
            printlog ("ERROR - GAFASTSKY: gemcombine returned a \
                non-zero status. Exiting.", l_logfile, verbose+)
            goto crash
        } else {
            printlog ("GAFASTSKY: Returned from GEMCOMBINE.\n", \
               l_logfile, l_verbose)
        }

        if (!imaccess(tmpimgout)) {
            printlog ("\nERROR - GAFASTSKY: Outout image \""//tmpimgout//\
                "\" not created. Exiting", l_logfile, verbose+)
            goto crash
        }

        if (debug) gemdate (verbose+)

        outphu = tmpimgout//"[0]"

        # Update the gain and units of the variance planes and remove
        # BUNIT keyword from DQ planes!
        if (l_fl_vardq) {
            keypar (outphu, "NSCIEXT", silent+)
            nsci = int(keypar.value)

            for (i = 1; i <= nsci; i += 1) {
                # Update the units of the variance plane
                keypar (tmpimgout//"["//l_var_ext//","//i//"]", \
                    "BUNIT", silent+)

                if (stridx("*",str(keypar.value)) == 0) {
                    gemhedit (tmpimgout//"["//l_var_ext//","//i//"]", \
                        "BUNIT", keypar.value//"*"//keypar.value, "", \
                        delete-)
                }

                # Delete the units of the DQ plane
                gemhedit (tmpimgout//"["//l_dq_ext//","//i//"]", \
                    "BUNIT", "", "", delete+)

                # Update the title of the DQ plane
                gemhedit (tmpimgout//"["//l_dq_ext//","//i//"]", \
                    "i_title", "Data Quality", "", delete-)
            }
        }

        # Update the title of the whole image
        gemhedit (outphu, "i_title", l_title, "", delete-)

        # Write statsec to PHU
        gemhedit (outphu, "GAFSTSTA", l_statsec, \
            "Statistics region used by GAFASTSKY", delete-)

        # Write combination type to PHU
        gemhedit (outphu, "GAFSTCOM", l_combine, \
            "Type of combine used by GAFASTSKY", delete-)

        # Write rejection method to PHU
        gemhedit (outphu, "GAFSTREJ", l_reject, \
            "Type of rejection used by GAFASTSKY", delete-)

        if (l_reject == "minmax") {
            # Write nigh and nlow to PHU
            gemhedit (outphu, "GAFSTNLO", l_nlow, \
                "Low pixels rejected (minmax)", delete-)
            gemhedit (outphu, "GAFSTNHI", l_nhigh, \
                "High pixels rejected (minmax)", delete-)
        }

        if (debug) gemdate (verbose+)

        # Write input files to PHU
        nn = 1
        scanfile2 = conflist
        while (fscan(scanfile2, inimg) != EOF) {
            ##M Bug here if more than 100 files!
            gemhedit (outphu, "GAFSIM"//str(nn), inimg, \
                "Input image combined with GAFASTSKY", delete-)
            nn += 1
        }

        gemdate (zone="UT")

        gemhedit (outphu, "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI", delete-)

        gemhedit (outphu, "GAFSTSKY", gemdate.outdate, \
            "UT Time stamp for GAFASTSKY", delete-)

        printlog ("GAFASTSKY: Created "//tmpimgout, l_logfile, l_verbose)

        if (debug) gemdate (verbose+)

        printlog ("\nGAFASTSKY: Finished processing configuration "//\
            metaconf//"\n----", l_logfile, l_verbose)

        # Add to the tally of files created
        num_created += 1

NEXT_CONFIGURATION:

        imdelete (tmpbpm, verify-, >& "dev$null")

        # Counter for number of configurations for creating lists
        nconf += 1

        # Keep track if any output flats were created
        if (num_created == 0) {
            printf ("%22s"//metaconf//"\n", " ", >> tmp_skipped)
        }

        delete (conflist, verify-, >& "dev$null")

        # Track the total number of images created
        total_num_created += num_created

    } # End of loop over unique configurations

    goto clean
    #-----------------------------------------------------------------------

crash:

    status = 1

clean:

    scanfile = ""
    scanfile2 = ""
    scanfile3 = ""

    # Check for configurations that were skipped and if any fils were created
    if (status == 0) {
        if (access(tmp_skipped)) {
            printlog ("\nWARNING - GAFASTSKY: The following configurations \
                were not processed - \n", l_logfile, verbose+)
            type (tmp_skipped) | tee (l_logfile, out_type="text", append+)
            delete (tmp_skipped, verify-, >& "dev$null")

            if (total_num_created == 0) {
                printlog ("\nERROR - GAFASTSKY: No output files created", \
                    l_logfile, verbose+)
                status = 1
            }
        }
    }

    # delete output list from gaimchk output checks
    if (access(tmpout_gaimchk)) {
        delete ("@"//tmpout_gaimchk//","//tmpout_gaimchk, verify-, \
            >& "dev$null")
    }

    # delete output list from gaimchk input checks
    if (access(outgaimchk_list)) {
        delete ("@"//outgaimchk_list//","//outgaimchk_list, verify-, \
            >& "dev$null")
    }

    imdelete (tmpbpm, verify-, >& "dev$null")

    delete (filelist//","//tmplist//","//tmpconflist, verify-, >& "dev$null")

    # Print finish time
    gemdate(zone="local")
    printlog ("\nGAFASTSKY -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("\nGAFASTSKY -- Exit status: GOOD", l_logfile, l_verbose)

    } else if (status == 1) {
        printlog ("\nGAFASTSKY -- Exit status: ERROR", l_logfile, l_verbose)

    }

    if (status != 0) {
        printlog ("          -- Please read the logfile \""//l_logfile//\
                "\" for more information.", l_logfile, l_verbose)
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)

end

