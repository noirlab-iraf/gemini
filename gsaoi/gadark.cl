# Copyright(c) 2010-2015 Association of Universities for Research in Astronomy, Inc.

procedure gadark (inimages)

# TODO
# - optimize parameters for combine.
# Make sure BPM is of the correct size for gemcombine call

# INPUTS:
#     GSAOI raw or gaprepare'd dark data
#     This task is smart and will find the DARKs in the given location if
#     *.fits is supplied.

# This task will:
#     combine darks per extension within a certain time interval from the first
#     one in the stack
#     only checks the exposure characteristics. All other configuration
#     information is discarded

# Version 2010 May 20 V0.0 CW - created
#         2010 Jul 06 V1.0 CW - implemented with vardq from gemcombine
#         2011 Apr 20 V1.2 CW - combine parameters depending on number of
#                               images allowing to process a single dark (may
#                               be required for twilights).
#         2011 Sep 28 V1.3 CW - small modification to handle ROIs (check on
#                               image sizes)
#         2011 Dec 07 V1.4 CW - and changed it back with the new metaconf
#
# Descriptions of future revisions are recorded in CVS logs

######################################################
# input and image selection
char    inimages    {prompt="Input GSAOI dark images; raw or prepared"}
char    rawpath     {"",prompt="Path for input images"}
char    outsufx     {"dark",prompt="Suffix for output dark"}
char    rootname    {"",prompt="Root name if supplying image number(s); blank for today's UT"}
int     mindark     {5,min=1,prompt="Minimum number of dark images to combine"}
real    maxtime     {3600.,prompt="Maximum time interval from first image in the list"}
char    datename    {"DATE-OBS",prompt="Date header keyword"}
char    timename    {"UT",prompt="Time stamp header keyword"}
bool    fl_vardq    {no,prompt="Create variance and data quality frames"}
char    badpix      {"gsaoi$data/gsaoibpm_high_full.fits",prompt="Static (MEF) Bad Pixel Mask"}

# gaprepare
char    gaprep_pref {"g",prompt="Prefix for output images"}
bool    fl_trim     {yes,prompt="Trim the images?"}
bool    fl_nlc      {yes,prompt="Apply non-linear correction to data?"}
bool    fl_sat      {yes,prompt="Include non-linear and saturated pixels in data quality planes"}
char    arraysdb    {"gsaoi$data/gsaoiAMPS.dat",prompt="Database file containing array information"}
char    non_lcdb    {"gsaoi$data/gsaoiNLC.dat",prompt="Database file containing non-linearity correction coefficients"}

# gemcombine parameters
char    combine     {"default",enum="default|average|median",prompt="Combination operation"}
char    reject      {"avsigclip",enum="none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip",prompt="Rejection algorithm"}
char    masktype    {"goodvalue",enum="none|goodvalue",prompt="Bad pixel mask type"}
real    maskvalue   {0.,prompt="Good pixel value in the BPM"}
char    scale       {"none",prompt="Image scaling (none|mode|median|mean|exposure|@<file>|!<keyword>)"}
char    zero        {"none",prompt="Image zeropoint offset (none|mode|median|mean|@<file>|!<keyword>)"}
char    weight      {"none",prompt="Image weights (none|mode|median|mean|exposure|@<file>|!<keyword>)"}
char    statsec     {"[*,*]",prompt="Statistics section for combining"}
char    expname     {"EXPTIME",prompt="Exposure time header keyword"}
real    lthreshold  {INDEF,prompt="Lower threshold"}
real    hthreshold  {INDEF,prompt="Upper threshold"}
int     nlow        {1,min=0,prompt="minmax: Number of low pixels to reject"}
int     nhigh       {1,min=0,prompt="minmax: Number of high pixels to reject"}
int     nkeep       {1,prompt="Minimum to keep or maximum to reject"}
bool    mclip       {yes,prompt="Use median in sigma clipping algorithms?"}
real    lsigma      {3.,min=0,prompt="Lower sigma clipping factor"}
real    hsigma      {3.,min=0,prompt="Upper sigma clipping factor"}
char    key_ron     {"RDNOISE",prompt="Keyword for readout noise in e-"}
char    key_gain    {"GAIN",prompt="Keyword for gain in electrons/ADU"}
real    ron         {0.0,min=0.,prompt="Readout noise rms in electrons"}
real    gain        {1.0,min=0.00001,prompt="Gain in e-/ADU"}
char    snoise      {"0.0",prompt="ccdclip: Sensitivity noise (electrons)"}
real    sigscale    {0.1,min=0.,prompt="Tolerance for sigma clipping scaling correction"}
real    pclip       {-0.5,prompt="pclip: Percentile clipping parameter"}
real    grow        {0.0,min=0.,prompt="Radius (pixels) for neighbor rejection"}
bool    fl_dqprop   {no,prompt="Propagate all DQ values?"}

# other stuff
char    sci_ext     {"SCI",prompt="Name of science extensions"}
char    var_ext     {"VAR",prompt="Name of variance extensions"}
char    dq_ext      {"DQ",prompt="Name of data quality extensions"}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

######################################################

begin

    ########
    # Declare and set local variables; set default values; initiate temporary
    # files

    ####
    # Variable declaration
    char l_inimages, l_rawpath, l_outsufx, l_logfile, l_rootname
    char l_badpix, l_combine, l_reject, l_masktype
    char l_scale, l_zero, l_weight, l_statsec, l_expname
    char l_key_ron, l_key_gain, l_snoise, l_timename, l_datename
    char l_temp, filelist, utdate, inimg, l_dir, inphu, metaconf
    char onlist, inlist, offlist, darkname, refname, l_redirect
    char tmplist, inimg_root, inname, tmpnotproc, l_orig_combine
    char l_gaprep_pref, l_sci_ext, l_var_ext, l_dq_ext, todelete, tmpgaprep
    char nextlist, tmpconflist, conflist, tmpinsert, new_metaconf
    char outgaimchk_list, tmpfile, l_obs_allowed, l_obj_allowed
    char l_arraysdb, l_non_lcdb, pr_string, def_combine, def_reject
    char tmpbpm, bpmgemcomb, area_bpm

    int i, l_mindark, l_nlow, l_nhigh, l_nkeep, num_in, nmatch, comma_pos
    int lcount, num_created, total_num_created, nconf, nsci

    real l_maskvalue, l_lthreshold, l_hthreshold, l_lsigma, l_hsigma
    real l_ron, l_gain, l_sigscale, l_pclip, l_grow, l_maxtime, l_distance

    bool l_fl_vardq, l_mclip, l_fl_dqprop, l_verbose
    bool l_fl_trim, l_fl_nlc, l_fl_sat, skipped, processing
    bool debug

    ####
    ##M Make sure these are in the same order as the parameters
    # Set local variables
    l_inimages = inimages
    l_rawpath = rawpath
    l_outsufx = outsufx
    l_logfile = logfile
    l_rootname = rootname
    l_badpix = badpix
    l_combine = combine
    l_reject = reject
    l_masktype = masktype
    l_scale = scale
    l_zero = zero
    l_weight = weight
    l_statsec = statsec
    l_expname = expname
    l_key_ron = key_ron
    l_key_gain = key_gain
    l_snoise = snoise
    l_timename = timename
    l_datename = datename
    l_fl_vardq = fl_vardq
    l_mclip = mclip
    l_verbose = verbose
    l_mindark = mindark
    l_nlow = nlow
    l_nhigh = nhigh
    l_nkeep = nkeep
    l_maskvalue = maskvalue
    l_lthreshold = lthreshold
    l_hthreshold = hthreshold
    l_lsigma = lsigma
    l_hsigma = hsigma
    l_ron = ron
    l_gain = gain
    l_sigscale = sigscale
    l_pclip = pclip
    l_grow = grow
    l_maxtime = maxtime
    l_gaprep_pref = gaprep_pref
    l_fl_dqprop = fl_dqprop
    l_fl_trim = fl_trim
    l_fl_nlc = fl_nlc
    l_fl_sat = fl_sat
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_arraysdb = arraysdb
    l_non_lcdb = non_lcdb

    ####
    # Set temporary files
    filelist = mktemp ("tmpfile")
    tmpnotproc = mktemp ("tmpnotproc")
    tmpgaprep = mktemp ("tmpgaprep")
    nextlist = mktemp ("tmpnextlist")
    tmplist = mktemp ("tmplist")
    todelete = mktemp ("tmptodelete")
    outgaimchk_list = mktemp ("tmpgaimchk")
    tmpbpm = mktemp ("tmpbpm")//".fits"

    ####
    # Set default values
    l_distance = INDEF
    status = 0
    nconf = 0
    inlist=""
    offlist=""
    total_num_created = 0
    l_orig_combine = l_combine
    debug = no

    if (l_verbose) {
        l_redirect = "STDOUT"
    } else {
        l_redirect = "dev$null"
    }

    # Default combine parameters - there ar edefault nhigh/nlow parameters
    # based on number of files too in the l_combine=default loop. Also, these
    # parameters get overridden if number of files are 4 or less. - MS
    def_combine = "average"
    def_reject = "avsigclip"

    ########
    # Here is where the actual work starts

    # Test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = gsaoi.logfile
       if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
          l_logfile = "gsaoi.log"
          printlog ("WARNING - GADARK: Both gaprepare.logfile and \
              gsaoi.logfile fields are empty", l_logfile, l_verbose)
          printlog ("                     Using default file gsaoi.log", \
              l_logfile, l_verbose)
       }
    }

    # Print start time
    gemdate (zone="local")
    printlog ("\nGADARK -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    # Record all input parameters relavent to this task only - other tasks will
    # print their inputs to log
    printlog ("GADARK: Input parameters...", l_logfile, l_verbose)
    printlog ("    inimages   = "//l_inimages, l_logfile, l_verbose)
    printlog ("    rawpath    = "//l_rawpath, l_logfile, l_verbose)
    printlog ("    outsufx    = "//l_outsufx, l_logfile, l_verbose)
    printlog ("    mindark    = "//l_mindark, l_logfile, l_verbose)
    printlog ("    maxtime    = "//l_maxtime, l_logfile, l_verbose)
    printlog ("    datename   = "//l_datename, l_logfile, l_verbose)
    printlog ("    timename   = "//l_timename, l_logfile, l_verbose)
    printlog ("    fl_vardq   = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("    badpix     = "//l_badpix, l_logfile, l_verbose)
    printlog ("    logfile    = "//l_logfile, l_logfile, l_verbose)
    printlog ("    verbose    = "//l_verbose, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    ####
    # Start user input checks and data validation

    # Check the l_badpix value
    if ((l_badpix == "") || (stridx(" ", l_badpix) > 0)) {
        # This assumes that the bpm covers all of the areas readout in each
        # detector in each image
        l_badpix = "none"
    } else if (!imaccess(l_badpix)) {
        printlog ("ERROR - GADARK: cannot access bpm "//l_badpix//\
            ". Exiting", l_logfile, verbose+)
    } else {
        # Check it hasn't been mosaiced
        keypar (l_badpix//"[0]", "GAMOSAIC", silent+)
        if (keypar.found) {
            printlog ("ERROR - GADARK: "//l_badpix//" has been mosaiced."//\
                " Please supply a static BPM that has not been mosaiced.", \
                l_logfile, verbose)
            goto crash
        }
    }
    printlog ("GADARK: badpix value after checks is \""//l_badpix//"\"\n", \
        l_logfile, l_verbose)

    # The checks in gaimchk convert these strings to lower case!
    l_obs_allowed = "DARK"
    l_obj_allowed = ""

    # Call gaimchk to perform inout checks
    gaimchk (inimages=l_inimages, rawpath=l_rawpath, rootname=l_rootname, \
        obstype_allowed=l_obs_allowed, object_allowed=l_obj_allowed, \
        key_allowed="GAPREPAR", key_forbidden="GADARK,GMOSAIC", \
        key_exists="", fl_prep_check=yes, gaprep_pref=l_gaprep_pref, \
        fl_redux_check=no, garedux_pref="r", fl_fail=no, fl_out_check=no, \
        fl_vardq_check=l_fl_vardq, sci_ext=l_sci_ext, var_ext=l_var_ext, \
        dq_ext=l_dq_ext, outlist=outgaimchk_list, logfile=l_logfile, \
        verbose=no)

    if (gaimchk.status != 0) {
        printlog ("ERROR - GADARK: GAIMCHK returned a non-zero status. \
            Exiting.", l_logfile, verbose+)
        goto crash
    } else if (!access(outgaimchk_list)) {
        printlog ("ERROR - GADARK: Cannot access output list from GAIMCHK", \
            l_logfile, verbose+)
        goto crash
    }

    # Files that have been processed by GADARK
    tmpfile = ""
    match ("tmpGADARK", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        printlog ("WARNING - GADARK: The following files have already been \
            processed by GADARK. Ignoring these files.", l_logfile, verbose+)
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no) | \
            tee (l_logfile, out_type="text", append+)
    }

    # Prepared files
    tmpfile = ""
    match ("tmpGAPREPAR", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        type (tmpfile, >> nextlist)
    }

    # Files that need to be prepared
    tmpfile = ""
    match ("tmpNotFound", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {

        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> tmpgaprep)

        printlog ("GADARK: Calling GAPREPARE to process unprepared input \
            files", l_logfile, l_verbose)

        gaprepare ("@"//tmpgaprep, rawpath="",outpref=l_gaprep_pref, \
            rootname="", fl_trim=l_fl_trim, fl_nlc=l_fl_nlc, \
            fl_vardq=l_fl_vardq, fl_sat=l_fl_sat, badpix=l_badpix, \
            sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
            arraysdb=l_arraysdb, non_lcdb=l_non_lcdb, logfile=l_logfile, \
            verbose=l_verbose)

        if (gaprepare.status != 0) {
            printlog ("ERROR - GADARK: GAPREPARE returned an non-zero status. \
                Exiting.", l_logfile, verbose+)
            goto crash
        } else {
            printlog ("GADARK: Returned from GAPREPARE\n", \
                l_logfile, l_verbose)
        }

        fields (tmpfile, fields=2, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> tmplist)
        sections (l_gaprep_pref//"//@"//tmplist, \
            option="fullname", >> nextlist)

        delete (tmplist, verify-, >& "dev$null")
        delete (tmpgaprep, verify=no, >& "dev$null")

    } # End of prepare loop

    # Check if there are any images in the next list list
    if (!access(nextlist)) {
        printlog ("ERROR - GADARK: No input images can be used. \
            Please try again.", l_logfile, verbose+)
        goto crash
    } else {
        print (nextlist, >> todelete)
    }

    # Check if there are any images in the next list list
    if (!access(nextlist)) {
        printlog ("ERROR - GADARK: No input images can be used. \
            Please try again.", l_logfile, verbose+)
        goto crash
    }

    # Loop over nextlist to obtain the METACONFIG keyword. This contains
    # information that will group files such that a) they are in the same
    # configuration and b) that they have been reduced in the same way.
    printlog ("GADARK: Accessing METACONF keyword values", \
        l_logfile, l_verbose)

    scanfile = nextlist
    while (fscan(scanfile, inimg, inname) != EOF) {

        inphu = inimg//"[0]"

        # Double check that you can access the file...
        if (!imaccess(inimg)) {
            printlog ("WARNING - GADARK: Image "//inname//\
                " does not exist. Skipping this image.", l_logfile, verbose)
            print (inimg, >> tmpnotproc)
        } else {
            # Read the METACONFIG keyword
            keypar (inphu, "METACONF", silent+)
            if (!keypar.found) {
                printlog ("WARNING - GADARK: METACONF keyword not found in "//\
                    inname//" Skipping this image.", l_logfile, verbose+)
                print (inimg, >> tmpnotproc)
            } else {
                metaconf = keypar.value

                # Check it's valid
                if (metaconf == "UNKOWN" || metaconf == "" || \
                    (stridx(" ",metaconf) > 0)) {
                    printlog ("WARNING - GADARK: METACONF keyword is "//\
                        "\"UNKNOWN\" or bad in"//inname//".\n"//\
                        "                  Skipping this image.", \
                        l_logfile, verbose+)
                    print (inimg, >> tmpnotproc)
                } else {
                    print (inimg//" "//metaconf, >> filelist)
                }
            }
        }

    } # End of loop to read METACONF keyword

    # Check if there are any images in the newly created filelist
    if (!access(filelist)) {
        printlog ("ERROR - GADARK: No images with a valid METACONF keyword \
            value. Please try again.", l_logfile, verbose+)
        goto crash
    }

    # Determine the unique METACONF values for the file that made it
    printlog ("GADARK: Determining unique METACONFIG values", \
        l_logfile, l_verbose)

    tmpconflist = mktemp ("tmpconflist")

    fields (filelist, fields=2, lines="1-", quit_if_miss=no, \
        print_file_names=no) | sort ("STDIN", column=0, ignore_white=no, \
        numeric_sort=no, reverse_sort=no) | unique ("STDIN", > tmpconflist)

    # Keep track of files to delete at the end
    print (tmpconflist, >> todelete)

    ####
    # Create the DARKs...

    # Create separate lists for each unique METACONF keyword
    printlog ("GADARK: Creating METACONFIG lists and processing images...", \
        l_logfile, l_verbose)

    scanfile = tmpconflist
    while (fscan(scanfile, metaconf) != EOF) {

        printlog ("\n----", l_logfile, l_verbose)
        printlog ("GADARK: Processing configuration "//metaconf//"\n", \
                    l_logfile, l_verbose)

        # Set up temporary configuration list
        conflist = tmplist//"_conf"//nconf//".lis"

        # Extract the files that have the same METACONF keyword value
        match (metaconf, filelist, stop=no, print_file_name=no, \
            metacharacter=no) | fields ("STDIN", fields=1, lines="1-", \
            quit_if_miss=no, print_file_names=no) | \
            sort ("STDIN", column=0, ignore_white=no, \
            numeric_sort=no, reverse_sort=no, > conflist)

        # Keep track of files to delete at the end
        print (conflist, >> todelete)

        # Default BPM for gemcombine
        bpmgemcomb = ""

        # Check the BPM dimensions against the first image in the list
        if (l_badpix != "none" && l_masktype != "none") {

            # Do this as it's possible (outside chance) that some files don't
            # have VAR/DQ whilst others do - MS

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
                    printlog ("ERROR - GADARK: GADIMSCHK returned a non-zero \
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
                        printlog ("ERROR - GADARK: The array dimensions are \
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

        # Now for each configuration list, check the time offset. The first
        # dark in the list is the reference and it parses the list for all
        # others that met the time constraint. Then it takes what is left and
        # repeats until none is left this can use gemoffsetlist both ways and
        # also to check the mindark constraint

        # Set flags, temporary files to be used inside the next loop
        processing = yes
        onlist = conflist
        inlist = mktemp("tmpinlist")
        offlist = mktemp("tmpofflist")
        skipped = no # This says whether the current configuation list has been
                     # processed more than once

        # Keep track of files to delete at the end
        print (inlist, >> todelete)
        print (offlist, >> todelete)

        if (debug) {
            printlog ("onlist -- "//conflist, l_logfile, verbose+)
            type (onlist) | tee (l_logfile, out_type="text", append+)
        }

        # Reset the num_created counter
        num_created = 0

        # Stay in this loop until all files have been exhausted
        while (processing) {

            # Determine the number of files in onlist
            count (onlist) | scan (num_in)

            # Only continue if the number of files in onlist are greater than
            # the minimum requested number of darks required to create a master
            # dark
            if (num_in < mindark) {

                # Can get here imediately or after processing the list already.
                # As such adjust the print statements accordingly

                if (skipped) {
                    tmpinsert = "now "
                } else {
                    tmpinsert = ""
                }

                printlog ("\nWARNING - GADARK: Number of input files for \
                    configuration "//metaconf//\
                    "\n                  is "//tmpinsert//\
                    "less than the minimum "//\
                    "requested "//mindark//" files.", \
                    l_logfile, verbose+)

                if (!skipped) {
                    printlog ("                  Skipping this \
                        configuration", l_logfile, verbose+)
                } else {
                    printlog ("                  Moving on from this \
                        configuration", l_logfile, verbose+)
                }

                # Print warning if no darks created
                type (onlist, >> tmpnotproc)
                processing = no

                if (num_created == 0) {
                    printlog ("\nWARNING - GADARK: No master darks created \
                        for configuration "//metaconf, l_logfile, verbose+)
                }

                goto MISSCONFIG
            }

            if (debug) {
                printlog ("HERE 1", l_logfile, verbose+)
            }

            # Read the reference file
            head (onlist, nlines=1) | scan (refname)

            if (debug) {
                printlog ("Refname -- "//refname, l_logfile, verbose+)
            }

            # Find any matching files by time
            # l_distance is INDEF so only seperating by time
            gemoffsetlist ("@"//onlist, refname, l_distance, l_maxtime, \
                inlist, offlist, fl_nearer=yes, direction=3, \
                fl_younger=yes, fl_noref=no, wcs_source="none", \
                key_xoff="XOFFSET",  key_yoff="YOFFSET", \
                key_date=l_datename, key_time=l_timename, \
                logfile=l_logfile, verbose=no, force=no, >& "dev$null")

            # Check the output status
            if (gemoffsetlist.status != 0) {
                printlog ("ERROR - GADARK: GEMOFFSETLIST returned a \
                    non-zero status. Exiting.", l_logfile, verbose+)
                goto crash
            }

            # Number of matches found -- will always contain the reference file
            nmatch = gemoffsetlist.count

            if (debug) {
                printlog ("NMATCH -- "//nmatch, l_logfile, verbose+)

                printlog ("HERE 3 matched list...", l_logfile, verbose+)

                type (inlist) | tee (l_logfile, out_type="text", append+)

                if (access(offlist)) {
                    printlog ("HERE 3b offlist...", l_logfile, verbose+)
                    type (offlist) | tee (l_logfile, out_type="text", append+)
                }
            }

            # Check the number of matches against the requested minimum
            # number of darks to combine
            if (nmatch < l_mindark) {
                # If not enough, remove the first and use the next as
                # reference
                printlog ("WARNING - GADARK: number of images in"//\
                    " list within "//l_maxtime//" seconds of", \
                    l_logfile, l_verbose)
                printlog ("                  reference "//refname//\
                    " less than minimum.", l_logfile, l_verbose)

                # Remove the temporary files
                delete (offlist//","//inlist, verify-, >& "dev$null")

                # Print the refernce file to a list to be printed out at
                # the end
                print (refname, >> tmpnotproc)

                # Remove the first file (reffile) from onlist and reset
                # onlist
                tail (onlist, nl=-1, >> offlist)

                # Rename offlist to be the new onlist
                delete (onlist, verify-, >& "dev$null")
                rename (offlist, onlist, field="all")

                if (debug) {
                    printlog ("HERE 3c new onlist...", l_logfile, verbose+)
                    type (onlist) | tee (l_logfile, out_type="text", append+)
                }

                # Set the skipped flag
                skipped = yes

            } else {
                # Enough were selected to combine [there may be files that
                # were not selected too]

                if (debug) {
                    printlog ("HERE 4", l_logfile, verbose+)
                }

                # Now the header is from the first  image, need the name to
                # be as well - so reset darkname and set up output name
                head (inlist, nlines=1) | scan (refname)
                darkname = refname
                fparse (darkname, verbose-)
                darkname = fparse.root//"_"//l_outsufx//".fits"

                # Check if the output image exists, if it does, exit
                # gracefully
                if (imaccess(darkname)) {
                    printlog ("ERROR - GADARK: output image "//darkname//\
                        " already exists. Exiting.", l_logfile, l_verbose)
                    goto crash
                }

                printlog ("GADARK: Creating output dark - "//darkname, \
                    l_logfile, l_verbose)

                # Check for a signal image to be combined, only copy if so. Can
                # only get here is user sets min_dark == 1 - MS
                if (nmatch == 1) {
                    printlog ("\nWARNING - GADARK: only a single image has \
                        been selected to be combined. Just copying image.", \
                        l_logfile, verbose+)

                    fxcopy (refname, darkname, groups="", new_file=yes, \
                        verbose=no)

                } else {

                    # Combine the darks

                    # Set the default gemcombine parameters if required
                    # Regardless of flat type this will always be what the
                    # user sets it to for each configuration. - MS
                    if (l_orig_combine == "default") {

                        printlog ("\nGADARK: Combine parameter is \
                            \"default\". \
                            \n        Setting default gemcombine combine \
                            and rejection parameters.", \
                            l_logfile, l_verbose)

                        # Refine by type - defined by location_type - defaults
                        # set at top of script
                        l_combine = def_combine
                        l_reject = def_reject

                        l_nhigh = 1
                        l_nlow = 1

                        # This overides the default settings if the number of
                        # inputs is too small
                        if (nmatch < 5) {
                            l_combine = "average"
                            l_reject = "minmax"
                            printlog("\nWARNING - GADARK: Averaging 4 or \
                                fewer images with 1 low pixels rejected"//\
                                "\n                  and only 1 high pixel \
                                rejected", l_logfile, verbose+)
                        }
                    } # End of default section

                    # Inform the user of the combination and rejection methods
                    pr_string = strupr(substr(l_combine,1,1))//\
                        strlwr(substr(l_combine,2,strlen(l_combine)))

                    printlog ("\nGADARK: "//pr_string//\
                        " combining dark frames...", \
                        l_logfile, l_verbose)

                    if (l_reject == "minmax") {
                        pr_string = " with "//l_nlow//" low and "//l_nhigh//\
                        " high values rejected"
                    } else {
                        pr_string = ""
                    }

                    printlog ("        Rejection type is "//l_reject//\
                        pr_string, l_logfile, l_verbose)


                    printlog ("\nGADARK: Input files are...\n", \
                        l_logfile, l_verbose)

                    type (inlist) | tee (l_logfile, out_type="text", append+, \
                        >> l_redirect)

                    # If fl_vardq=yes, the VAR plane is the variance of the
                    # SCI, and the DQ from imcombine unless fl_dqprop is yes.
                    # If no DQ planes exist, the BPM will be used in the future
                    # -MS
                    # NOTE: If the VAR should be built from the individual VAR
                    # planes, the code needs to be changed - cannot use
                    # GEMCOMBINE (or needs another loop for the VAR/DQ only)
                    # Doing this type of variance propagation is currently not
                    # suppoerted and may never be - MS

                    printlog ("", l_logfile, l_verbose)

                    gemcombine ("@"//inlist, darkname, \
                        title="Combined Dark", \
                        combine=l_combine, reject=l_reject, \
                        logfile=l_logfile, offsets="none", \
                        masktype=l_masktype, maskvalue=l_maskvalue, \
                        scale=l_scale, zero=l_zero, weight=l_weight, \
                        statsec=l_statsec, expname=l_expname, \
                        lthreshold=l_lthreshold, hthreshold=l_hthreshold, \
                        nlow=l_nlow, nhigh=l_nhigh, \
                        nkeep=l_nkeep, mclip=l_mclip, lsigma=l_lsigma, \
                        hsigma=l_hsigma, key_ron=l_key_ron, \
                        key_gain=l_key_gain, ron=l_ron, gain=l_gain, \
                        snoise=l_snoise, sigscale=l_sigscale,
                        pclip=l_pclip, \
                        grow=l_grow, bpmfile=bpmgemcomb, nrejfile="", \
                        sci_ext=l_sci_ext, var_ext=l_var_ext, \
                        dq_ext=l_dq_ext, fl_vardq=l_fl_vardq, \
                        fl_dqprop=l_fl_dqprop, \
                        verbose=l_verbose)

                    # Check the output status of gemcombine
                    if (gemcombine.status != 0) {
                        printlog ("ERROR - GADARK: gemcombine returned a \
                            non-zero status. Exiting.", l_logfile, verbose+)
                        goto crash
                    } else {
                        printlog ("GADARK: Returned from GEMCOMBINE\n", \
                            l_logfile, l_verbose)
                    }
                }

                # Check a master dark was created
                if (imaccess(darkname)) {
                    num_created += 1
                    total_num_created +=1
                } else {
                    printlog ("ERROR - GADARK: Cannot access "//darkname, \
                        l_logfile, verbose+)
                    goto crash
                }

                inphu = darkname//"[0]"

                # Update the gain and units of the variance planes and remove
                # BUNIT keyword from DQ planes!
                if (l_fl_vardq) {
                    keypar (inphu, "NSCIEXT", silent+)
                    nsci = int(keypar.value)

                    for (i = 1; i <= nsci; i += 1) {
                        # Update the units of the variance plane
                        keypar (darkname//"["//l_var_ext//","//i//"]", \
                            "BUNIT", silent+)

                        if (stridx("*",str(keypar.value)) == 0) {
                            gemhedit (darkname//"["//l_var_ext//","//i//"]", \
                                "BUNIT", keypar.value//"*"//keypar.value, "", \
                                delete-)
                        }

                        # Delete the units of the DQ plane
                        gemhedit (darkname//"["//l_dq_ext//","//i//"]", \
                            "BUNIT", "", "", delete+)
                    }
                }

                # Edit the header of the final flat to add a keyword
                # GADARK and update GEM-TLM both with the time stamp.
                # GEMCOMBINE will update GAIN and READNOISE
                # Update the METACONF keyword too

                new_metaconf = substr(metaconf,1,stridx("+",metaconf))//\
                    "DARK+"//\
                    substr(metaconf,stridx("+",metaconf)+1,strlen(metaconf))
                gemhedit (inphu, "METACONF", new_metaconf, "", delete-)

                gemdate ()
                gemhedit (inphu, "GEM-TLM", gemdate.outdate, \
                    "", delete-)
                gemhedit (inphu, "GADARK", gemdate.outdate, \
                    "UT Time stamp for GADARK", delete-)

                # If there's a list of files that didn't get selected
                # rename it and set it back to be the next onlist list
                if (access(offlist)) {
                    if (debug) {
                        printlog ("HERE 6", l_logfile, verbose+)
                    }

                    delete (onlist, verify-, >& "dev$null")
                    rename (offlist, onlist, field="all")
                    delete (inlist, verify-, >& "dev$null")
                    skipped = yes

                } else {
                    processing = no
                }
            } # End of else for if nmatch < mindark

MISSCONFIG:
        } # End of while processing loop

        imdelete (tmpbpm, verify-, >& "dev$null")

        printlog ("GADARK: Finished processing configuration "//\
            metaconf, l_logfile, l_verbose)
        printlog ("----\n", l_logfile, l_verbose)

        delete (conflist, verify-, >& "dev$null")

        nconf += 1

    } # End of loop over configuation values

    goto clean

    #--------------------------------------------------------------------------
crash:
    # Exit with error subroutine
    status = 1

clean:

    # Print the list of files that were not processed for some reason
    if (access(tmpnotproc)) {
        printlog ("WARNING - GADARK: The following files were not processed. \
            \n                  Check the logfile \""//l_logfile//\
            "\" for more information.\n",\
            l_logfile, verbose+)
        type (tmpnotproc) | tee (l_logfile, out_type="text", append+)
        delete (tmpnotproc, verify=no, >& "dev$null")
    }

    imdelete (tmpbpm, verify-, >& "dev$null")

    # Clean up todelete list if it exists
    if (access(todelete)) {
        delete("@"//todelete, verify=no, >& "dev$null")
        delete(todelete, verify=no, >& "dev$null")
    }

    if (access(outgaimchk_list)) {
        delete ("@"//outgaimchk_list, verify=no, >& "dev$null")
    }

    delete (filelist//","//nextlist//","//outgaimchk_list, \
        verify-, >& "dev$null")

    scanfile = ""

    # Check that some master biases were created - else change status to error
    if (status == 0 && total_num_created == 0) {
        printlog ("\nERROR - GADARK: No master darks were created.", \
            l_logfile, verbose+)
        status = 1
    }

    # Print finish time
    gemdate(zone="local")
    printlog ("GADARK -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("   ", l_logfile, l_verbose)
        printlog ("GADARK -- Exit status: GOOD", l_logfile, l_verbose)
    } else {
        printlog ("   ", l_logfile, l_verbose)
        printlog ("GADARK -- Exit status: ERROR", l_logfile, verbose+)
    }

    if (status != 0) {
        printlog ("       -- Please read the logfile \""//l_logfile//\
                "\" for more information.", l_logfile, verbose+)
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)
end

