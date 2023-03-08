# Copyright(c) 2010-2015 Association of Universities for Research in Astronomy, Inc.

procedure gasky(inimages)

# Create sky frame and object masks from input. Inputs must be gaprepared! They
# must also be none-flats or darks!

# Source detection is performed before combining.

##M Add ccdsec and datasec parameters
##M Check for parked ODGWs
#    l_key_froz_odgw = "SKYFODGW"
##M Remove guiding ODGW CCDSEC keyword if not frozen in outfile
##M Remove creation of expression database and related files / variables

######################################################

char    inimages     {prompt="GSAOI sky images to combine"}
char    outimages    {"",prompt="Output sky images"}
char    outsufx      {"sky",prompt="Output suffix if outimages=\"\""}
char    combine      {"default",prompt="Type of combine operation (default|median|average)",enum="default|median|average"}
char    reject       {"minmax",prompt="Type of rejection (none|minmax|avsigclip)",enum="none|minmax|avsigclip"}
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
bool    fl_keepmasks {no,prompt="Keep object masks for each input image?"}
char    masksuffix   {"msk",prompt="Mask name suffix"}
real    threshold    {3.,min=1.5,prompt="Threshold in sigma for object detection"}
int     ngrow        {3,prompt="Number of iterations to grow objects into the wings"}
real    agrow        {3.,prompt="Area limit for growing objects into the wings"}
int     minpix       {6,prompt="Minimum number of pixels to be identified as an object"}
bool    fl_mask      {yes,prompt="Mask non-good pixels during source detection & calculation of statistics?"}
bool    qreduce      {yes,prompt="Quickly reduce images, with sky from GAFASTSKY for object masking?"}
char    qred_msktype {"goodvalue",enum="none|goodvalue",prompt="Bad Pixel Mask type (none|goodvalue) for GAFASTSKY"}
char    flatimg      {"",prompt="Flat field image, for use when qreduce=yes (\"\"|<filename>)"}
real    minflat_val  {1E-6,prompt="Minimum allowed pixel value in flat image; only for qreduce=yes"}
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
    char l_inimages, l_outimages, l_combine, l_reject, l_logfile
    char l_key_exptime, l_key_filter, l_statsec, outgaimchk_list, l_outsufx
    char l_sci_ext, l_var_ext, l_dq_ext, l_key_ron, l_key_gain
    char l_rootname, l_rawpath, l_obs_allowed, l_obj_allowed, l_gaprep_pref
    char l_garedux_pref, l_key_allowed, tmpfile, l_redirect, filelist
    char tmplist, tmpconflist, conflist, tmpoutlist, tmpimgout
    char one_ron, one_gain, separator, orig_1_ron, orig_1_gain
    char tmpout_gaimchk, inimg, inphu, insci, inname, inPrname, metaconf
    char l_separator, tmp_skipped, pr_string, outphu, l_title, l_masktype
    char l_scale, l_zero, l_weight, l_snoise, l_badpix, l_nrejfile
    char l_orig_combine, tmpfastsky, tmpflat, tmpgemcomb, tmpmasks
    char l_outtype, l_dqouttype, maskroot, maskname, l_sci_expr, l_var_expr
    char l_dq_expr, img_root, l_masksuffix, l_stattype, l_statextn, statvalue
    char stdvalue, sciextn, sciappend, dqextn, dqappend, l_flatimg, inflat
    char tmpinimg, inobjimg, dqfile, tmpgemfile, tmpmaskname, out_metaconf
    char tmpexprdb, cal_flat, todelete, tmpbpm, bpmgemcomb, area_bpm, refname
    char sci_end, var_end, tstring, tmpfcal, fpath, reffile, l_qred_msktype
    char scirep_val, varrep_val, dqrep_val

    int total_num_created, num_created, nconf, i, nsci, lcount, n_out, l_nkeep
    int last_ron_char, last_gain_char, l_nclip, l_test, l_test2, last_char
    int last2_char, l_ngrow, l_minpix, src_hthreshold, src_lthreshold
    int l_nlow, l_nhigh, nn, nfiles

    real l_lsigma, l_hsigma, stat_num, std_num, l_agrow
    real l_lthreshold, l_hthreshold, l_maskvalue
    real l_sigscale, l_pclip, l_grow, l_minflat_val
    real l_expone, l_ron, l_gain

    bool debug, process, reduce, noflatin, notflat_divided, l_fl_keepmasks
    bool l_fl_mask, l_mclip, write_to_header, can_use_one_file, del_masks
    bool l_verbose, l_fl_dqprop, l_fl_vardq

    ####
    # Set local variables
    l_inimages = inimages
    l_outimages = outimages
    l_outsufx = outsufx
    l_key_exptime = key_exptime
    l_orig_combine = combine
    l_reject = reject
    l_nlow = nlow
    l_nhigh = nhigh
    l_masktype = masktype
    l_maskvalue = maskvalue
    l_badpix = badpix
    l_flatimg = flatimg
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
    l_fl_keepmasks = fl_keepmasks
    l_masksuffix = masksuffix
    reduce = qreduce
    l_fl_mask = fl_mask
    l_agrow = agrow
    l_ngrow = ngrow
    l_minpix = minpix
    l_qred_msktype = qred_msktype

    ####
    # Make temporary files
    outgaimchk_list = mktemp ("tmpoutgaimchk")
    filelist = mktemp ("tmpfilelist")
    tmplist = mktemp ("tmplist")
    tmpconflist = mktemp ("tmpconflist")
    tmpoutlist = mktemp ("tmpoutlist")
    tmpout_gaimchk = mktemp ("tmpoutchk")
    tmp_skipped = mktemp ("tmpskip")
    tmpexprdb = mktemp ("tmpexprdb")
    todelete = mktemp ("tmptodelete")
    tmpfastsky = mktemp ("tmpfastsky")
    tmpflat = mktemp ("tmpflat")
    tmpgemcomb = mktemp ("tmpgemcombin")
    tmpmasks = mktemp ("tmpmasks")
    tmpbpm = mktemp ("tmpbpm")
    tmpfcal = mktemp ("tmpfcals")//".fits"

    ####
    # Set default values
    status = 0
    debug = no
    process = no
    can_use_one_file = no
    del_masks = yes # To delete existsing masks that do not have the correct
                    # number of extensions

    # Default values for source detection
    src_hthreshold = threshold
    src_lthreshold = 10

    # Used in the gaimchk call
    l_rootname = ""
    l_rawpath = ""
    l_obs_allowed = "OBJECT"
    l_obj_allowed = ""
    l_gaprep_pref = "g"
    l_garedux_pref = "r" # No longer used but required for call.

    # Pixel types for output images
    l_outtype = "real"
    l_dqouttype = "ushort"

    # Set this flag to yes if writing the coming step values to the PHU
    write_to_header = yes

    # dqrep_val and varrep_val for gemexpr call
#    varrep_val = "0"
#    dqrep_val = "256"
    varrep_val = "a[VAR]"
    dqrep_val = "a[DQ]"


    # This is for tee'd output
    if (l_verbose) {
        l_redirect = "STDOUT"
    } else {
        l_redirect = "dev$null"
    }

    #------------------------------------------------------------------------
    # Fun starts here

    # Cache in memory, tasks that are used
    cache("gaimchk", "gadimschk", "gemexpr", "gemcombine", "gafastsky")

    # Test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = gsaoi.logfile
       if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
          l_logfile = "gsaoi.log"
          printlog ("WARNING - GASKY: Both gaprepare.logfile and \
              gsaoi.logfile fields are empty", l_logfile, l_verbose)
          printlog ("                     Using default file gsaoi.log", \
              l_logfile, l_verbose)
       }
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)

    # Print start time
    printlog ("", l_logfile, verbose+)
    gemdate (zone="local")
    printlog ("GASKY -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    # Record all input parameters
    printlog ("GASKY: Input parameters...", l_logfile, l_verbose)
    printlog ("    inimages    = "//l_inimages, l_logfile, l_verbose)
    printlog ("    outimage    = "//l_outimages, l_logfile, l_verbose)
    printlog ("    outsufx     = "//l_outsufx, l_logfile, l_verbose)
    printlog ("    key_exptime = "//l_key_exptime, l_logfile, l_verbose)
    printlog ("    combine     = "//l_orig_combine, l_logfile, l_verbose)
    printlog ("    statsec     = "//l_statsec, l_logfile, l_verbose)
    printlog ("    reject      = "//l_reject, l_logfile, l_verbose)
    printlog ("    nlow        = "//l_nlow, l_logfile, l_verbose)
    printlog ("    nhigh       = "//l_nhigh, l_logfile, l_verbose)
    printlog ("    masktype    = "//l_masktype, l_logfile, l_verbose)
    printlog ("    maskvalue   = "//l_maskvalue, l_logfile, l_verbose)
    printlog ("    masksuffix  = "//l_masksuffix, l_logfile, l_verbose)
    printlog ("    badpix      = "//l_badpix, l_logfile, l_verbose)
    printlog ("    qreduce     = "//reduce, l_logfile, l_verbose)
    printlog ("    flatimg     = "//l_flatimg, l_logfile, l_verbose)
    printlog ("    sci_ext     = "//l_sci_ext, l_logfile, l_verbose)
    printlog ("    var_ext     = "//l_var_ext, l_logfile, l_verbose)
    printlog ("    dq_ext      = "//l_dq_ext, l_logfile, l_verbose)
    printlog ("    key_ron     = "//l_key_ron, l_logfile, l_verbose)
    printlog ("    key_gain    = "//l_key_gain, l_logfile, l_verbose)
    printlog ("    fl_vardq    = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("    fl_dqprop   = "//l_fl_dqprop, l_logfile, l_verbose)
    printlog ("    logfile     = "//l_logfile, l_logfile, l_verbose)

    # Create the expr database used by imexpr
    if (!isindef(minflat_val)) {
        tstring = str(minflat_val)
        l_test = fscan (tstring,l_minflat_val)
        if (l_test != 1) {
            printlog ("ERROR - GASKY: minflat_val is not a number. \
                Exiting", l_logfile, verbose+)
            goto crash
        } else {
            # Create a temporary expr databae for imexpr
            # Basically checking for near zero numbers
            print ("chkval(a)        ((abs(a) < "//l_minflat_val//") ? "//\
                "((a < 0) ? (-"//l_minflat_val//") : "//l_minflat_val//\
                ") : a)", > tmpexprdb)
            printlog ("    minflat_val = "//l_minflat_val, \
                l_logfile, l_verbose)
            scirep_val = str(l_minflat_val)
        }
    } else {
        l_minflat_val = minflat_val
        printlog ("    minflat_val = INDEF", \
            l_logfile, l_verbose)

        # Create a temporary expr databae for imexpr
        print ("chkval(a)        (a)", > tmpexprdb)

        scirep_val = "0"
    }
    printlog ("", l_logfile, l_verbose)

    if (debug)
        type (tmpexprdb) | tee (l_logfile, out_type="text", append+)

    # Check the flatimg setting though only if processing
    if (reduce && l_flatimg != "") {
        if (stridx(" ",l_flatimg) > 0) {
            printlog ("ERROR - GASKY: qreduce=yes but flatimg contains \
                spaces. Either supply a flatimg name or a null string \"\"", \
                l_logfile, verbose+)
            goto crash
        } else if (!imaccess(l_flatimg)) {
            printlog ("ERROR - GASKY: Cannot access flatimg \""//l_flatimg//\
                "\"", l_logfile, verbose+)
            goto crash
        } else {
            # Check it's not been mosaiced
            keypar (l_flatimg//"[0]", "GAMOSAIC", silent+)
            if (keypar.found) {
                printlog ("ERROR - GASKY: flatimg has been mosaiced with \
                    GAMOSIC. Please use a non-mosaiced version.", \
                    l_logfile, verbose+)
                goto crash
            } else {
                # Check it's actually a flat image
                keypar (l_flatimg//"[0]", "GAFLAT", silent+)
                if (!keypar.found) {
                    printlog ("ERROR - GASKY: flatimg \""//l_flatimg//\
                        "was not created by GAFLAT.", \
                        l_logfile, verbose+)
                    got crash
                }
            }
        }

        fparse (l_flatimg)
        l_flatimg = fparse.root//fparse.extension
        fpath = fparse.directory
        if (fpath == "") {
            fpath = "./"
        }

        # Send it off to gacalfind for checking later
        gacalfind (inimages=l_flatimg, calpath=fpath, \
            caltable=tmpfcal, fl_calrun=yes, fl_find=no, caltype="", \
            sciimg="", maxtime=INDEF, ignore_nlc=no, \
            sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
            logfile=l_logfile, verbose=no)

        if (gacalfind.status != 0) {
            printlog ("ERROR - GASKY: GACALFIND returned a non-zero \
                status. Exiting.", l_logfile, verbose+)
            goto crash
        } else if (!access(tmpfcal)) {
            printlog ("ERROR - GASKY: Cannot access \""//tmpfcal//\
                "\". Exiting.", l_logfile, verbose+)
            goto crash
        }
    } # End of flatimg check

    # Check the l_badpix value
    if ((l_badpix == "") || (stridx(" ", l_badpix) > 0)) {
        # This assumes that the bpm covers all of the areas readout in each
        # detector in each image
        l_badpix = "none"
    } else if (!imaccess(l_badpix)) {
        printlog ("ERROR - GASKY: cannot access bpm "//l_badpix//\
            ". Exiting", l_logfile, verbose+)
    } else {
        # Check it hasn't been mosaiced
        keypar (l_badpix//"[0]", "GAMOSAIC", silent+)
        if (keypar.found) {
            printlog ("ERROR - GASKY: "//l_badpix//" has been mosaiced."//\
                " Please supply a static BPM that has not been mosaiced.", \
                l_logfile, verbose)
            goto crash
        }
    }

    # Verify that the extension names are not empty, otherwise exit gracefully
    if (l_sci_ext == "" || stridx(" ",l_sci_ext) > 0) {
        printlog("ERROR - GASKY: extension name sci_ext is missing.", \
            l_logfile, verbose+)
        status=1
        goto clean
    }

    if (l_fl_vardq) {
        if (l_dq_ext == "" || (stridx(" ",l_dq_ext) > 0) || \
            l_var_ext == "" || (stridx(" ",l_var_ext) > 0)) {

            printlog("WARNING - GASKY: var_ext or dq_ext have not been \
                set.\n                     Output image will not have VAR or \
                DQ planes.", l_logfile, verbose+)
            l_fl_vardq=no
        }
    }

    if (l_fl_dqprop && !l_fl_vardq) {
        printlog ("WARNING - GASKY: Cannot propagate DQ information with \
            fl_vardq=no. Setting fl_dqprop=no", l_logfile, verbose+)
            l_fl_dqprop = no
    }

    printlog ("GASKY: Calling GAIMCHK to check inputs...", \
        l_logfile, l_verbose)

    if (debug) gemdate (verbose+)

    # Call gaimchk to perform input checks
    l_key_allowed=""

    gaimchk (inimages=l_inimages, rawpath=l_rawpath, rootname=l_rootname, \
        obstype_allowed=l_obs_allowed, object_allowed=l_obj_allowed, \
        key_allowed=l_key_allowed, key_forbidden="GASKY,GAFASTSKY", \
        key_exists="GAPREPAR", fl_prep_check=yes, gaprep_pref=l_gaprep_pref, \
        fl_redux_check=no, garedux_pref=l_garedux_pref, fl_fail=no, \
        fl_vardq_check=l_fl_vardq, sci_ext=l_sci_ext, var_ext=l_var_ext, \
        dq_ext=l_dq_ext, outlist=outgaimchk_list, logfile=l_logfile, \
        verbose=no)

    if (gaimchk.status != 0) {
        printlog ("ERROR - GASKY: GAIMCHK returned a non-zero status. \
            Exiting.", l_logfile, verbose+)
        goto crash
    } else if (!access(outgaimchk_list)) {
        printlog ("ERROR - GASKY: Cannot access output list from GAIMCHK", \
            l_logfile, verbose+)
        goto crash
    } else {
        printlog ("GASKY: Returned from GAIMCHK", l_logfile, l_verbose)
    }

    # Files that have been processed by GASKY
    tmpfile = ""
    match ("tmpGASKY", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        printlog ("\nWARNING - GASKY: The following files have already been \
            processed by GASKY. Ignoring these files...\n", \
            l_logfile, verbose+)
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no) | \
            tee (l_logfile, out_type="text", append+)
        printlog ("", l_logfile, verbose+)
    }

    # Files that were prepared
    tmpfile = ""
    match ("tmpGAPREPAR", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> tmplist)
    }

    # Files that were not prepared
    tmpfile = ""
    match ("tmpNotFound", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        printlog ("\nWARNING - GASKY: The following files have not been \
            processed by GAPREPARE. Ignoring these files...\n", \
            l_logfile, verbose+)
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no) | \
            tee (l_logfile, out_type="text", append+)
        printlog ("", l_logfile, verbose+)
    }

    if (access(tmplist)) {
        count (tmplist) | scan (lcount)
    }

    if (lcount == 0) {
        printlog ("ERROR - GASKY: No unique metaconfigurations values "//\
            "found.\n               Exiting", l_logfile, verbose+)
        goto crash
    }

    if (debug) gemdate (verbose+)

    # Call gaimchk to perform output checks
    n_out = 0
    if (l_outimages != "" && stridx(" ",l_outimages) == 0) {

        printlog ("GASKY: Calling GAIMCHK to check outputs...", \
            l_logfile, l_verbose)

        if (debug) gemdate (verbose+)

        gaimchk (inimages=l_outimages, rawpath="", rootname="", \
            obstype_allowed="", object_allowed="", key_allowed="", \
            key_forbidden="", key_exists="", fl_prep_check=no, \
            gaprep_pref=l_gaprep_pref, fl_redux_check=no, \
            garedux_pref=l_garedux_pref, fl_fail=no, \
            fl_out_check=yes, fl_vardq_check=no, sci_ext=l_sci_ext, \
            var_ext=l_var_ext, dq_ext=l_dq_ext, outlist=tmpout_gaimchk, \
            logfile=l_logfile, verbose=no)

        if (gaimchk.status != 0) {
            printlog ("ERROR - GASKY: GAIMCHK returned a non-zero status. \
                Exiting.", l_logfile, verbose+)
            goto crash
        } else if (!access(tmpout_gaimchk)) {
            printlog ("ERROR - GASKY: Cannot access output list from \
                GAIMCHK", l_logfile, verbose+)
            goto crash
        }

        # Files that need to be prepared
        tmpoutlist = ""
        match ("tmpNotFound", tmpout_gaimchk, stop=no, print_file_name=no, \
            metacharacter=no) | scan (tmpoutlist)

        if (tmpoutlist == "") {
            printlog ("ERROR - GASKY: Output list does not exist. \
                Exiting", l_logfile, verbose+)
            goto crash
        } else {
            count (tmpoutlist) | scan (n_out)
            if (n_out == 0) {
                printlog ("ERROR - GASKY: Output list is empty. Exiting", \
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

    if (l_outimages == "default" \
        && (l_outsufx == "" || stridx(" ", l_outsufx) > 0)) {
        printlog ("ERROR - GASKY: outimages is \"\" but outsufx is not \
            set properly. Exiting", l_logfile, verbose+)
        goto crash
    }

    # Loop over nextlist to obtain the METACONFIG keyword. This contains
    # information that will group files such that a) they are in the same
    # configuration and b) that they have been reduced in the same way.
    printlog ("\nGASKY: Accessing METACONF keyword values", \
        l_logfile, l_verbose)

    if (debug) gemdate (verbose+)

    scanfile = tmplist

    while (fscan(scanfile, inimg) != EOF) {

        inphu = inimg//"[0]"

        # Double check that you can access the file...
        if (!imaccess(inimg)) {
            printlog ("WARNING - GASKY: Image "//inimg//\
                " does not exist. Skipping this image.", l_logfile, verbose)

        } else {
            # Read the METACONFIG keyword
            keypar (inphu, "METACONF", silent+)
            if (!keypar.found) {
                printlog ("WARNING - GASKY: METACONF keyword not \
                    found in "//inimg//" Skipping this image.", \
                    l_logfile, verbose+)

            } else {
                metaconf = keypar.value

                # Check it's valid
                if (metaconf == "UNKNOWN" || metaconf == "" \
                    || (stridx(" ",metaconf) > 0)) {
                    printlog ("WARNING - GASKY: METACONF keyword is "//\
                        "\"UNKNOWN\" or bad in"//inimg//".\n"//\
                        "                     Skipping this image.", \
                        l_logfile, verbose+)

                } else {

                    # Check it's an OBJECT!
                    if (strstr("OBJ",metaconf) == 0) {
                        printlog ("WARNING - GASKY: Image \""//inimg//\
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
        printlog ("ERROR - GASKY: No images with a valid METACONF keyword \
            value. Please try again.", l_logfile, verbose+)
        goto crash
    }

    if (debug) gemdate (verbose+)

    # Check metaconf keywords
    # Determine the unique METACONF values for the file that made it
    printlog ("GASKY: Determining unique METACONFIG values", \
        l_logfile, l_verbose)

    # Reverse sort - MS
    fields (filelist, fields=2, lines="1-", quit_if_miss=no, \
        print_file_names=no) | sort ("STDIN", column=0, ignore_white=no, \
        numeric_sort=no, reverse_sort=yes) | unique ("STDIN", > tmpconflist)

    # Check if there are any images in the newly created tmpconflist
    lcount = 0
    if (access(tmpconflist)) {
        count (tmpconflist) | scan (lcount)
    }

    if (lcount == 0) {
        printlog ("ERROR - GASKY: No unique METACONF keywords found. \
            Please try again.", l_logfile, verbose+)
        goto crash
    } else if (l_outimages != "default" && lcount != n_out) {
        printlog ("ERROR - GASKY: Number of supplied output image \
            names does"//\
            "\n                   not match the number of unique \
            METACONFIGURATION keyword values. Exiting", l_logfile, verbose+)
        goto crash
    }

    ####
    # Create the FLATs...

    total_num_created = 0

    # This becomes a maze of lists! Reader beware.

    if (debug) gemdate (verbose+)

    # Create separate lists for each unique METACONF keyword
    printlog ("GASKY: Creating METACONFIG lists and processing images...",\
        l_logfile, l_verbose)

    nconf = 1
    scanfile = tmpconflist
    while (fscan(scanfile, metaconf) != EOF) {

        num_created = 0

        printlog ("\n----\nGASKY: Processing configuration: "//\
            metaconf//"\n", l_logfile, l_verbose)

        if (debug) gemdate (verbose+)

        # Set up temporary configuration list
        conflist = tmplist//"_conf"//nconf//".lis"

        print (conflist, >> todelete)

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
            printlog ("ERROR - GASKY: Cannot access input \
                configuration list", l_logfile, verbose+)
            goto crash
        }
        if (lcount == 0) {
            printlog ("ERROR - GASKY: Configuration list for "//metaconf//\
                " is empty. Exiting", l_logfile, verbose+)
            goto crash
        } else if (lcount == 1) {
            printlog ("WARNING - GASKY: Only one file to combine for "//\
                metaconf//".",  l_logfile, verbose+)

            if (!can_use_one_file) {
                printlog ("                 Skipping this configuration.", \
                    l_logfile, verbose+)
                goto NEXT_CONFIGURATION
            }
        }
        nfiles = lcount

        printlog ("GASKY: Using input files...\n", l_logfile, l_verbose)
        type (conflist) | \
            tee (l_logfile, out_type="text", append+, >> l_redirect)

        if (debug) gemdate (verbose+)

        # Loop over inputs
        l_expone = 0.0
        nn = 0
        scanfile2 = conflist
        while (fscan(scanfile2, inimg) != EOF) {

            inphu = inimg//"[0]"
            fparse (inimg, verbose-)
            img_root = fparse.root
            inname = img_root//fparse.extension

            if (nn == 0) {
                l_title = ""
                # Read the OBJECT for the first frame
                keypar (inphu, "OBJECT", silent+)
                if (keypar.found) {
                    l_title = keypar.value
                }
            }

            if (l_outimages == "default" && nn == 0) {
                tmpimgout = fparse.root//"_"//l_outsufx//".fits"
                if (imaccess(tmpimgout)) {
                    printlog ("\nERROR - GASKY: "//tmpimgout//\
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

            # Set output PHU
            outphu = tmpimgout//"[0]"

            if (nfiles == 1 && can_use_one_file) {
                printlog ("GASKY: Just copying input and upadting headers", \
                    l_logfile, l_verbose)
                fxcopy (inimg, tmpimgout, groups="", new_file=yes, verbose-)
                goto UPDATE_ONLY
            }

            # check exposure time - Though this should already be done by the
            # unique metaconfig extraction
            keypar (inphu, l_key_exptime, silent+)
            if (!keypar.found) {
                printlog("ERROR - GASKY: "//l_key_exptime//" keyword not \
                    found in "//inname, l_logfile, verbose+)
                goto crash
            } else if (nn == 0) {
                l_expone = real(keypar.value)
            } else if (abs(real(keypar.value) - l_expone) > 0.1) {
                printlog("WARNING - GASKY: exposure times are \
                    significantly different. Continuing.", l_logfile,verbose+)
            }

            keypar (inphu, "NSCIEXT", silent+)
            if (!keypar.found) {
                printlog("ERROR - GASKY: NSCIEXT keyword not \
                    found in "//inname, l_logfile, verbose+)
                goto crash
            } else {
                nsci = int (keypar.value)
            }

            # Check readnoise and gain values are consistent
            for (i = 1; i <= nsci; i += 1) {

                insci = "["//l_sci_ext//","//i//"]"
                inPrname = inname//insci
                insci = inimg//insci

                # Check the dimensions of the image against the requested
                # statsec
                gadimschk (insci, section=l_statsec, chkimage="", \
                    key_check="CCDSEC", logfile=l_logfile, verbose-)

                if (gadimschk.status != 0) {
                    printlog ("ERROR - GASKY: GADIMSCHK returned a non-zero \
                        status. Exiting.", l_logfile, verbose-)
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
                        one_ron = substr(one_ron,last_ron_char+2, \
                            strlen(one_ron))
                        one_gain = substr(one_gain,last_gain_char+2, \
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
                    printlog("ERROR - GASKY: "//l_key_ron//" keyword not \
                        found in "//inPrname, l_logfile, verbose+)
                    goto crash
                } else if (nn == 0) {
                    one_ron = one_ron//l_separator//keypar.value

                    if (i == nsci) {
                        orig_1_ron = one_ron
                    }
                } else if (abs(real(keypar.value) - l_ron) > 1.) {
                    printlog("WARNING - GASKY: read noise values are \
                        different.  Continuing, but the", l_logfile, verbose+)
                    printlog("                     read noise in the output \
                        header will be wrong.", l_logfile, verbose+)
                }

                keypar (insci, l_key_gain, silent+)
                if (!keypar.found) {
                    printlog("ERROR - GASKY: "//l_key_gain//" keyword not \
                        found in "//inPrname, l_logfile, verbose+)
                    goto crash
                } else if (nn == 0) {
                    one_gain = one_gain//l_separator//keypar.value

                    if (i == nsci) {
                        orig_1_gain = one_gain
                    }
                } else if (abs(real(keypar.value)- l_gain) > 0.5) {
                    printlog("WARNING - GASKY: read noise values are \
                        different.  Continuing, but the", l_logfile, verbose+)
                    printlog("                     read noise in the output \
                        header will be wrong.", l_logfile, verbose+)
                }
            } # End of for loop over sci

            # Check for the presence of an object mask if requested to
            # reduce only reduce if a mask is not found
            maskroot = img_root//"_"//l_masksuffix//".fits"

            # Check mask exists and has the correct number of extensions
            if (imaccess(maskroot)) {
                # Check file is MEF
                gemextn (inimg, check="mef", process="expand", index="", \
                    extname=l_dq_ext, extversion="1-", ikparams="", omit="", \
                    replace="", outfile="dev$null", logfile=l_logfile, \
                    glogpars="", verbose=debug)

                if (gemextn.status != 0) {
                    printlog ("ERROR - GASKY: GEMEXTN returned a non-zero \
                        status. Exiting.", l_logfile, verbose+)
                    goto crash
                } else if (gemextn.fail_count > 0 || gemextn.count != nsci) {

                    printlog ("WARNING - GASKY: "//maskroot//" exists "//\
                        "but does not"//\
                        "\n                 have the appropriate number of "//\
                        "extensions", l_logfile, verbose+)

                    # Delete the masks if vinternal variable is set above
                    if (del_masks) {
                        printlog ("                 Deleting "//maskroot, \
                            l_logfile, verbose+)
                        imdelete (maskroot, verify-, >& "dev$null")
                    } else {
                        goto crash
                    }
                }
            }

            # Double check it's existance
            if (reduce && !access(maskroot)) {
                process = yes
            }

            nn += 1
        } # End of while of current configurations input list
        scanfile2 = ""

        if (nfiles != nn) {
            printlog ("ERROR - GASKY: Mismatch of number of files in \
                current configuration list", l_logfile, verbose+)
            goto crash
        }

        if (debug) gemdate (verbose+)

        if (process) {
            printlog ("\nGASKY: calling GAFASTSKY", l_logfile, l_verbose)

            gafastsky ("@"//conflist, outimage=tmpfastsky, \
                key_exptime=l_key_exptime, combine="default", \
                statsec=l_statsec, reject="minmax", \
                nlow=0, nhigh=1, \
                masktype=l_qred_msktype, maskvalue=l_maskvalue, \
                badpix=l_badpix, sci_ext=l_sci_ext, var_ext=l_var_ext, \
                dq_ext=l_dq_ext, key_ron=l_key_ron, key_gain=l_key_gain, \
                fl_vardq=l_fl_vardq, fl_dqprop=l_fl_dqprop, \
                logfile=l_logfile, verbose=no)

            if (gafastsky.status != 0) {
                printlog("ERROR - GASKY: GAFASTSKY returned a non-zero \
                    status. Exiting.", l_logfile, verbose+)
                goto crash
            } else {
                printlog ("\nGASKY: Returned from GAFASTSKY.", \
                    l_logfile, l_verbose)
            }

            # Call gastat to extract the mean
            l_stattype = "mean"
            l_statextn = "ARRAY"
            l_lsigma = 4.
            l_hsigma = 4.
            l_nclip = 1

            gastat (tmpfastsky, stattype=l_stattype, statsec=l_statsec, \
                statextn=l_statextn, fl_mask=l_fl_mask, \
                badpix=l_badpix, calc_stddev=yes, lower=INDEF, \
                upper=INDEF, nclip=l_nclip, lsigma=l_lsigma, \
                usigma=l_hsigma, sci_ext=l_sci_ext, dq_ext=l_dq_ext, \
                logfile=l_logfile, verbose-)

            if (gastat.status != 0) {
                printlog ("ERROR - GASKY: GASTAT returned a \
                    non-zero status. Exiting.", \
                    l_logfile, verbose+)
                goto crash
            }

            statvalue = gastat.outstat
            stdvalue = gastat.stddev

            if (debug) printlog ("stat_value: "//statvalue//", stdvalue: "//\
                stdvalue, l_logfile, verbose+)

            l_test = 1
            l_test2 = 1
            j = 0

            while (l_test > 0 && l_test2 > 0) {
                j += 1

                l_test = stridx (",",statvalue)
                l_test2 = stridx (",",stdvalue)

                if (l_test == 0) {
                    last_char = strlen(statvalue)
                } else {
                    last_char = l_test - 1
                }

                if (l_test2 == 0) {
                    last2_char = strlen(stdvalue)
                } else {
                    last2_char = l_test2 - 1
                }

                # Parse the current value
                stat_num = real(substr(statvalue,1,last_char))
                std_num = real(substr(stdvalue,1,\
                    last2_char))

                if (debug) printlog ("stat_num: "//stat_num//", std_num: "//\
                    std_num, l_logfile, verbose+)

                gemhedit (tmpfastsky//"["//l_sci_ext//","//j//"]", \
                    "GASTATVA", stat_num, l_stattype//" 4 sigma clip", \
                    delete-)

                gemhedit (tmpfastsky//"["//l_sci_ext//","//j//"]", \
                    "GASTATSD", std_num, l_stattype//" std deviation", \
                    delete-)

                # Rest the output statistics for next value
                statvalue = substr(statvalue,l_test+1,\
                    strlen(statvalue))

                stdvalue = substr(stdvalue,\
                    l_test2+1, strlen(stdvalue))

                if (debug) printlog ("stat_value: "//statvalue//\
                    ", stdvalue: "//stdvalue, l_logfile, verbose+)
            }

            # Check if images were flat divided
            notflat_divided = yes
            keypar (tmpfastsky//"[0]", "FLATIMG", silent+)
            if (keypar.found) {
                notflat_divided = no
                printlog ("GASKY: Images are alredy flat divided. Not \
                   creating temporary flat field from GAFASTSKY output", \
                   l_logfile, l_verbose)
            }

            # Create the temporary flat if needed
            if (l_flatimg == "" && notflat_divided) {

                printlog ("\nGASKY: Creating temporary flat field image from \
                    GAFASTSKY output", l_logfile, l_verbose)

                # Create the temporary flat

                sci_end = "(chkval((a / a.GASTATVA),"//scirep_val//"))"
                var_end = "(chkvval((a[sci] / a[sci].GASTATVA),"//scirep_val//\
                    ",a[VAR],"//varrep_val//"))"

                l_sci_expr = "(((a > (a.GASTATVA + (4 * a.GASTATSD))) || "//\
                    "(a < (a.GASTATVA - (4 * a.GASTATSD)))) ? "//\
                    "1 : "//sci_end//")"

                l_var_expr = "(((a[sci] > (a[sci].GASTATVA + "//\
                    "(4 * a[sci].GASTATSD))) || "//\
                    "(a[sci] < (a[sci].GASTATVA - "//\
                    "(4 * a[sci].GASTATSD)))) ? "//\
                    "0 : "//var_end//")"

                l_dq_expr = "chkdval((a[sci] / a[sci].GASTATVA),"//\
                    scirep_val//",a[DQ],"//dqrep_val//")"

                gemexpr (l_sci_expr, tmpflat, tmpfastsky, l_minflat_val, \
                    var_expr=l_var_expr, dq_expr=l_dq_expr, \
                    mdf_ext="MDF", fl_vardq=l_fl_vardq, \
                    dims="default", intype="default", outtype="ref", \
                    refim="a", rangecheck=yes, verbose=no, \
                    exprdb="gsaoi$data/gsaoiEDB.dat", logfile=l_logfile)

                if (gemexpr.status != 0) {
                    printlog ("ERROR - GASKY: GEMEXPR returned a non-zero \
                        status. Exiting.", l_logfile, verbose+)
                    goto crash
                }

                inflat = tmpflat
            } else if (!notflat_divided) {
                # Double negative - they have been flat divided
                inflat = ""
            } else {
                # All files in conflist have the same METACONF - need to
                # the input flat is appropriate - MS
                head (conflist) | scan (reffile)

                if (debug) printlog ("____Checking flatimg", \
                    l_logfile, verbose+)

                # Call gacalfind to find the correct / check the flat
                gacalfind (l_flatimg, calpath=fpath, \
                    caltable=tmpfcal, \
                    fl_calrun=no, fl_find=yes, caltype="FLAT", \
                    sciimg=reffile, maxtime=INDEF,
                    ignore_nlc=no, sci_ext=l_sci_ext, \
                    var_ext=l_var_ext, dq_ext=l_dq_ext, \
                    logfile=l_logfile, verbose=no)

                if (gacalfind.status != 0) {
                    printlog ("ERROR - GASKY: GACALFIND returned a \
                        non-zero status. Exiting.", l_logfile, verbose+)
                    goto crash
                } else if (gacalfind.out_calimg == "") {
                    printlog ("WARNING - GASKY: GACALFIND could not \
                        find an appropriate FLAT IMAGE for "//inPrname//\
                        ".\n                    Setting fl_flat = no for \
                        this image.", l_logfile, verbose+)
                    inflat = ""
                } else {
                    inflat = gacalfind.out_calimg
                }

                if (debug) printlog ("____FLATIMG: "//inflat, \
                    l_logfile, verbose+)
            }
        } # End of if process loop

        # Loop over input images to create the masks...
        scanfile2 = conflist
        while (fscan(scanfile2,inimg) != EOF) {

            inphu = inimg//"[0]"
            fparse (inimg, verbose-)
            img_root = fparse.root
            inname = img_root//fparse.extension
            maskroot = img_root//"_"//l_masksuffix//".fits"

            printlog ("\nGASKY: Checking for existing object masks for "//\
                inname, l_logfile, l_verbose)

            tmpinimg = img_root//"_"//mktemp("tmpin")//".fits"

##M TODO
            # Check if ODGW (if present were frozen) - if all or some fraction
            # are and that they have the same coordinates set flag to write
            # keyword to phu of outfile - l_key_froz_odgw - MS

            ##M Check if the mask exists and check that it's valid

            # Check for any existing masks
            if (!access(maskroot)) {
                if (process) {
                    # If not present and process

                    # Determine the scaling factor
                    l_stattype = "midpt"
                    l_statextn = "ARRAY"
                    l_lsigma = 4.
                    l_hsigma = 4.
                    l_nclip = 1

                    gastat (inimg, stattype=l_stattype, statsec=l_statsec, \
                        statextn=l_statextn, fl_mask=l_fl_mask, \
                        badpix=l_badpix, calc_stddev=no, lower=INDEF, \
                        upper=INDEF, nclip=l_nclip, lsigma=l_lsigma, \
                        usigma=l_hsigma, \
                        sci_ext=l_sci_ext, dq_ext=l_dq_ext, \
                        logfile=l_logfile, verbose-)

                    if (gastat.status != 0) {
                        printlog ("ERROR - GASKY: GASTAT returned a \
                            non-zero status. Exiting.", \
                            l_logfile, verbose+)
                        goto crash
                    }

                    statvalue = gastat.outstat

                    l_test = 1
                    j = 0

                    while (l_test > 0) {
                        j += 1

                        sciextn = "["//l_sci_ext//","//j//"]"
                        sciappend = "["//l_sci_ext//","//j//",append]"
                        dqextn = "["//l_dq_ext//","//j//"]"
                        dqappend = "["//l_dq_ext//","//j//",append]"

                        l_test = stridx (",",statvalue)

                        if (l_test == 0) {
                            last_char = strlen(statvalue)
                        } else {
                            last_char = l_test - 1
                        }

                        # Parse the current value
                        stat_num = real(substr(statvalue,1,last_char))

                        # Rest the output statistics for next value
                        statvalue = substr(statvalue,l_test+1,\
                            strlen(statvalue))

                        if (j == 1) {
                            imcopy (inphu, tmpinimg, verbose-)

                            printlog ("GASKY: processing image "//\
                                inname//" into \""//tmpinimg//"\"", \
                                l_logfile, l_verbose)
                        }

                        cal_flat = inflat//sciextn

                        if (debug)
                            printlog ("____cal_flat(1): "//cal_flat, \
                                l_logfile, verbose+)

                        if (inflat == "") {
                            cal_flat = 1
                        } else if (!imaccess(cal_flat)) {
                            printlog ("WARNING - GASKY: Cannot access "//\
                                cal_flat//" setting flat image to a \
                                value of 1", l_logfile, verbose+)
                            cal_flat = 1
                        } else {

                            if (debug)
                                printlog ("____cal_flat(2): "//cal_flat, \
                                    l_logfile, verbose+)

                            if (j == 1) {
                                printlog ("GASKY: Will use "//inflat//\
                                    " as flatfield", l_logfile, l_verbose)
                            }

                            # Check the dimensions of the images
                            gadimschk (inimg//sciextn, section="", \
                                chkimage=cal_flat, \
                                key_check="CCDSEC", logfile=l_logfile, \
                                verbose-)

                            if (gadimschk.status !=0) {
                                printlog ("ERROR = GADIMSCHK returned a \
                                    non-zero status. Exiting.", \
                                    l_logfile, verbose+)
                            } else {
                                cal_flat = gadimschk.out_chkimage
                            }
                        }

                        if (debug) printlog ("____cal_flat(3): "//cal_flat, \
                            l_logfile, verbose+)

                        imexpr ("((a - (b * (c / b.GASTATVA))) / "//\
                            "chkval(d,"//scirep_val//"))"//" + c", \
                            tmpinimg//sciappend, inimg//sciextn, \
                            tmpfastsky//sciextn, \
                            stat_num, cal_flat, dims="auto", intype="auto", \
                            outtype=l_outtype, refim="auto", btype="nearest", \
                            bpixval=0., rangecheck=yes, verbose=no, \
                            exprdb="gsaoi$data/gsaoiEDB.dat")

                        # Need to copy over DQ plane too!
                        if (imaccess(inimg//dqextn)) {
                            if (imaccess(tmpfastsky//dqextn)) {
                                imexpr ("(a | b)", tmpinimg//dqappend, \
                                    inimg//dqextn, tmpfastsky//dqextn, \
                                    dims="auto", intype="auto", \
                                    outtype=l_dqouttype, refim="auto", \
                                    btype="nearest", bpixval=0., \
                                    rangecheck=yes, verbose=no, exprdb="none")
                            } else {
                                imcopy (inimg//dqextn, tmpinimg//dqappend, \
                                    verbose-)
                            }
                        }

                        inobjimg = tmpinimg
                    }
                } else {
                    inobjimg = inimg
                }

                imcopy (inphu, maskroot, verbose-)

                # Create the masks
                for (i = 1; i <= nsci; i += 1) {
                    sciextn = "["//l_sci_ext//","//i//"]"
                    dqextn = "["//l_dq_ext//","//i//"]"
                    dqappend = "["//l_dq_ext//","//i//",append]"

                    if (debug) printlog ("____running objmask on "//sciextn, \
                        l_logfile, verbose+)

                    maskname = maskroot//dqappend
                    tmpmaskname = mktemp ("tmpmaskroot")//".pl"

                    if (imaccess(inobjimg//dqextn) && l_fl_mask) {
                        if (i == 1) {
                            printlog ("GASKY: Supplying DQ plane as mask"//\
                                " for "//inobjimg//\
                                "\n       during source detection.",\
                                l_logfile, l_verbose)
                        }
                        dqfile = inobjimg//dqextn
                    } else {
                        dqfile = ""
                    }

                    nproto.objmasks (inobjimg//sciextn, \
                        tmpmaskname, \
                        omtype="numbers", skys="", sigmas="", \
                        masks=dqfile, extnames="", logfiles=l_logfile,
                        blkstep=1, blksize=-10, convolve="block 3 3", \
                        hsigma=src_hthreshold, lsigma=src_lthreshold, \
                        hdetect=yes, ldetect=no, \
                        neighbors=8, minpix=l_minpix, ngrow=l_ngrow, \
                        agrow=l_agrow) | tee (l_logfile, out_type="text", \
                        append+, >& l_redirect)

                    if (!imaccess(tmpmaskname)) {
                        printlog ("ERROR - GASKY: nproto.objmasks did not \
                            produce a mask for "//inname//sciextn, \
                            l_logfile, verbose+)

                        printlog ("               There are likely \
                            extreme pixel values in the input image.", \
                            l_logfile, verbose+)

                        if (process) {
                            printlog ("               Try increasing the \
                                'minflat_val' value.", \
                                l_logfile, verbose+)
                        }

                        imdelete (maskroot//","//tmpinimg, verify-, \
                            >& "dev$null")
                        goto crash

                    } else if (dqfile != "") {
                        imexpr ("a != 0 ? a : ((b > 10) ? 128 : b)", maskname,\
                            dqfile, tmpmaskname, dims="auto", intype="auto", \
                            outtype=l_dqouttype, refim="auto", \
                            btype="nearest", bpixval=0., rangecheck=yes, \
                            verbose=no, exprdb="none")

                    } else {
                        # Reset all of the numbers to 128
                        imexpr ("a != 0 ? 128 : 0", maskname,\
                            tmpmaskname, dims="auto", intype="auto", \
                            outtype=l_dqouttype, refim="auto", \
                            btype="nearest", bpixval=0., rangecheck=yes, \
                            verbose=no, exprdb="none")
                    }
                    imdelete (tmpmaskname, verify-, >& "dev$null")
                }

                # Update header
                gemhedit (maskroot//"[0]", "i_title", "Object mask for "//\
                    img_root, "", delete-)

                gemhedit (maskroot//"[0]", "NEXTEND", nsci, "", delete-, \
                    upfile="")

                gemdate (zone="UT")
                gemhedit (maskroot//"[0]", "GEM-TLM", gemdate.outdate, \
                    "UT Last modification with GEMINI", delete-)

                gemhedit (maskroot//"[0]", "GASKY", gemdate.outdate, \
                    "UT Time stamp for GASKY", delete-)

                if (!l_fl_keepmasks) {
                    print (maskroot, >> tmpmasks)
                }

                if (process) {
                    imdelete (tmpinimg, verify-, >& "dev$null")
                }

            } else if (access(maskroot)) {
                printlog ("WARNING - Using existing mask "//maskroot, \
                    l_logfile, verbose+)
            }

            # Create inputs for gemcombine - need the masks to be part of the
            # input file to gemcombine
            tmpgemfile = img_root//"_"//mktemp("tmp")//".fits"
            imcopy (inphu, tmpgemfile, verbose-)

            for (i = 1; i <= nsci; i += 1) {
                sciextn = "["//l_sci_ext//","//i//"]"
                sciappend = "["//l_sci_ext//","//i//",append]"
                dqextn = "["//l_dq_ext//","//i//"]"
                dqappend = "["//l_dq_ext//","//i//",append]"

                imcopy (inimg//sciextn, tmpgemfile//sciappend, verbose-)
                imcopy (maskroot//dqextn, tmpgemfile//dqappend, verbose-)
            }

            # Update the number of extensions so gemcombine doesn't complain
            gemhedit (tmpgemfile//"[0]", "NEXTEND", (2 * nsci), "", delete-)

            print (tmpgemfile, >> tmpgemcomb)

        } # End of loop over conflist to create object masks

        if (process) {
            imdelete (tmpfastsky, verify-, >& "dev$null")
            if (imaccess(tmpflat)) {
                imdelete (tmpflat, verify-, >& "dev$null")
            }
        }

        l_title = l_title//" SKY IMAGE from gemini.gsaoi.gasky"

        # Default gemcombine parameters
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

            printlog ("\nGASKY: combine is set to default. Setting default"//\
                "\n       GEMCOMBINE combine and rejection parameters.", \
                l_logfile, l_verbose)

            l_combine = "average"
            l_reject = "avsigclip"

            if (nfiles < 5) {
                printlog ("\nWARNING - GASKY: Averaging 4 or fewer images", \
                    l_logfile, verbose+)

                if (nfiles < 3) {
                    l_reject = "minmax"

                    l_nlow = 0
                    l_nhigh = 1

                    printlog ("\n                 with no low pixels \
                        rejected and only 1 high pixel rejected", \
                        l_logfile, verbose+)
                }
            }
        } else {

            l_combine = l_orig_combine

            if (nfiles < 5) {
                printlog("\nWARNING - GASKY: "//l_combine//" combining \
                    4 or fewer images", l_logfile,verbose+)
                if (l_reject == "none") {
                    printlog("                 with no pixels rejected.",\
                        l_logfile, verbose+)
                } else if (l_reject == "minmax") {
                    printlog("                 with "//l_nlow//\
                        " low and "//l_nhigh//" high pixels rejected.",
                        l_logfile, verbose+)
                }
            } # end if(nfiles < 5)

            if ((nfiles <= (l_nlow + l_nhigh)) && (l_reject == "minmax")) {
                printlog("\nERROR - GASKY: Cannot reject more pixels than \
                    the number of images.", l_logfile,verbose+)
                status=1
                goto crash
            }
        } # end not-default section

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

        pr_string = strupr(substr(l_combine,1,1))//\
            strlwr(substr(l_combine,2,strlen(l_combine)))

        printlog ("\nGASKY: "//pr_string//" combining sky frames...", \
            l_logfile, l_verbose)

        if (l_reject == "minmax") {
            pr_string = " with "//l_nlow//" low and "//l_nhigh//\
                " high values rejected"
        } else {
            pr_string = ""
        }

        printlog ("       Rejection type is "//l_reject//pr_string, \
            l_logfile, l_verbose)

        gemcombine ("@"//tmpgemcomb, output=tmpimgout, title=l_title, \
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
            printlog ("ERROR - GASKY: gemcombine returned a \
                non-zero status. Exiting.", l_logfile, verbose+)
            goto crash
        }

        imdelete ("@"//tmpgemcomb, verify-, >& "dev$null")
        delete (tmpgemcomb, verify-, >& "dev$null")

        printlog ("GASKY: Returned from GEMCOMBINE.\n", l_logfile, l_verbose)

        if (!imaccess(tmpimgout)) {
            printlog ("\nERROR - GASKY: Outout image \""//tmpimgout//\
                "\" not created. Exiting", l_logfile, verbose+)
            goto crash
        }

        if (debug) gemdate (verbose+)

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

        if (write_to_header) {
            # Write statsec to PHU
            gemhedit (outphu, "GASKSTA", l_statsec, \
                "Statistics region used by GASKY", delete-)

            # Write combination type to PHU
            gemhedit (outphu, "GASKCOM", l_combine, \
                "Type of combine used by GASKY", delete-)

            # Write rejection method to PHU
            gemhedit (outphu, "GASKREJ", l_reject, \
                "Type of rejection used by GASKY", delete-)

            if (l_reject == "minmax") {
                # Write nigh and nlow to PHU
                gemhedit (outphu, "GASKNLO", l_nlow, \
                    "Low pixels rejected (minmax)", delete-)
                gemhedit (outphu, "GASKNHI", l_nhigh, \
                    "High pixels rejected (minmax)", delete-)
            }
        }

        # Print that information to log file always
        printf("%-8s= \'%-20s\' / %-s\n", "GASKSTA", l_statsec, \
            "Statistics region used by GASKY", >> l_logfile)

        # Write combination type to logfile
        printf("%-8s= \'%-20s\' / %-s\n", "GASKCOM", l_combine, \
            "Type of combine used by GASKY", >> l_logfile)

        # Write rejection method to logfile
        printf("%-8s= \'%-20s\' / %-s\n", "GASKREJ", l_reject, \
            "Type of rejection used by GASKY", >> l_logfile)

        if (l_reject == "minmax") {
            # Write nigh and nlow to logfile
            printf("%-8s= %22.0f / %-s\n", "GASKNLO", l_nlow, \
                "Low pixels rejected (minmax)", >> l_logfile)
            printf("%-8s= %22.0f / %-s\n", "GASKNHI", l_nhigh, \
                "High pixels rejected (minmax)", >> l_logfile)
        }

        if (debug) gemdate (verbose+)

        # Write input files to PHU
        nn = 1
        scanfile2 = conflist
        while (fscan(scanfile2, inimg) != EOF) {
            if (write_to_header) {
                ##M Bug here if more than 100 files!
                gemhedit (outphu, "GASKIM"//str(nn), inimg, \
                    "Input image combined with GASKY", delete-)
            }
            printf("%-8s= \'%-20s\' / %-s\n","GASKIM"//str(nn), inimg, \
                "Input image combined with GASKY", >> l_logfile)
            nn += 1
        }

UPDATE_ONLY:

        # Update the METACONF keyword to read SKY instead of OBJ
        out_metaconf = substr(metaconf,1,strstr("OBJ",metaconf)-1)//"SKY"//\
            substr(metaconf,strstr("OBJ",metaconf)+strlen("OBJ"),\
            strlen(metaconf))

        gemhedit (outphu, "METACONF", out_metaconf, "", delete-)

        gemdate (zone="UT")
        gemhedit (outphu, "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI", delete-)

        gemhedit (outphu, "GASKY", gemdate.outdate, \
            "UT Time stamp for GASKY", delete-)

        printlog ("GASKY: Created "//tmpimgout, l_logfile, l_verbose)

        if (debug) gemdate (verbose+)

        printlog ("\nGASKY: Finished processing configuration "//\
            metaconf//"\n----", l_logfile, l_verbose)

        # Add to the tally of files created
        num_created += 1

NEXT_CONFIGURATION:

        imdelete (tmpbpm, verify-, >& "dev$null")

        # Counter for number of configurations for creating lists
        nconf += 1

        # Keep track if any output sky were created
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
            printlog ("\nWARNING - GASKY: The following configurations \
                were not processed - \n", l_logfile, verbose+)
            type (tmp_skipped) | tee (l_logfile, out_type="text", append+)
            delete (tmp_skipped, verify-, >& "dev$null")

            if (total_num_created == 0) {
                printlog ("\nERROR - GASKY: No output files created", \
                    l_logfile, verbose+)
                status = 1
            }
        }
    }

    imdelete (tmpfastsky, verify-, >& "dev$null")

    if (access(todelete)) {
        delete ("@"//todelete//","//todelete, verify-, >& "dev$null")
    }

    if (!l_fl_keepmasks && access(tmpmasks)) {
        imdelete ("@"//tmpmasks, verify-, >& "dev$null")
    }

    # delete output list from gaimchk output checks
    if (access(tmpout_gaimchk)) {
        delete ("@"//tmpout_gaimchk//","//tmpout_gaimchk, verify-, \
            >& "dev$null")
    }

    imdelete (tmpbpm, verify-, >& "dev$null")

    # delete output list from gaimchk input checks
    if (access(outgaimchk_list)) {
        delete ("@"//outgaimchk_list//","//outgaimchk_list, verify-, \
            >& "dev$null")
    }

    delete (filelist//","//tmplist//","//tmpconflist//","//tmpmasks//","//\
        tmpexprdb//","//tmpfcal, verify-, >& "dev$null")

    # Print finish time
    gemdate(zone="local")
    printlog ("\nGASKY -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("\nGASKY -- Exit status: GOOD", l_logfile, l_verbose)

    } else if (status == 1) {
        printlog ("\nGASKY -- Exit status: ERROR", l_logfile, l_verbose)

    }

    if (status != 0) {
        printlog ("      -- Please read the logfile \""//l_logfile//\
                "\" for more information.", l_logfile, l_verbose)
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)

end
