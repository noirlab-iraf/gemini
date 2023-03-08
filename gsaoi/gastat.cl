# Copyright(c) 2012-2015 Association of Universities for Research in Astronomy, Inc.

procedure gastat (inimage)

# This calculates the requested statistic for the given input file
# Only handles gsaoi images - one at a time
# BPM must be a MEF with the same number of extensions as the number of science
#     extensions in the input
# Cannot have sci_ext == dq_ext including ""
# If fl_mask, if no dq_ext present in input then looks for badpix. If badpix
#     not found fl_mask gets switched off

##M Add binwidth option!

char    inimage     {prompt="Input GSAOI image"}
char    stattype    {"mean",enum="mean|mode|midpt",prompt="Type of statistic to compute"}
char    statsec     {"[*,*]",prompt="Statistics section to use (relative to an array)"}
char    statextn    {"DETECTOR",prompt="How to apply statsec. (DETECTOR|ARRAY|<extnname,extnver>|<index>)"}
bool    fl_mask     {yes,prompt="Mask non-good pixels when calculating statistics?"}
char    badpix      {"gsaoi$data/gsaoibpm_high_full.fits", prompt="Static Bad Pixel Mask - not mosaic"}
bool    calc_stddev {no,prompt="Calculate standard deviation in calculation"}
real    lower       {INDEF,prompt="Lower good data value."}
real    upper       {INDEF,prompt="Maximum good data value."}
int     nclip       {0,prompt="Number of clipping iterations"}
real    lsigma      {INDEF,prompt="Lower side clipping factor in sigma."}
real    usigma      {INDEF,prompt="Upper side clipping factor in sigma."}
char    outstat     {INDEF,prompt="Output statistic (comma separated list for each extension)"}
char    stddev      {INDEF,prompt="Output standard deviation (comma separated list for each extension)"}
char    sci_ext     {"SCI",prompt="Extension name used for science data frames"}
char    dq_ext      {"DQ",prompt="Extension name used for data quality frames"}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}

begin

    ########
    # Declare and set local variables; set default values; initiate temporary
    # files
    char l_inimage, l_stattype, l_statsec, l_statextn, l_badpix, l_outstat
    char l_stddev, l_sci_ext, l_dq_ext, l_logfile
    char inphu, l_extnversion, l_index, scilist, dqlist, tmpdqlist, inextn
    char pr_upper, pr_lower, pr_lsigma, pr_usigma, dqextn, img, dqimg
    char tile_list, dqtile_list, stat_file, dq_file, dq_statsec, ccdsec
    char dqccdsec, pr_string, tstring, imgsec

    int test, extn, comma1_location, comma2_location, last_charpos
    int first_charpos, nsci, inextver, naxis1, naxis2, i, j
    int ccdx1, ccdx2, ccdy1, ccdy2, dqccdx1, dqccdx2, dqccdy1, dqccdy2
    int dqdatax1, dqdatax2, dqdatay1, dqdatay2, l_nclip
    int statx1, statx2, staty1, staty2, inpix, npix

    real l_upper, l_lower, statval, tmpstddev, num_pix_threshold, pix_ratio
    real l_lsigma, l_usigma

    struct val1

    bool l_fl_mask, l_calc_stddev, l_verbose, dotile, allcalc, isbad
    bool checkdq, docalc, dimscheck, sxcheck, sycheck
    bool debug

    ####
    # Variable declaration
    l_inimage = inimage
    l_stattype = stattype
    l_statsec = statsec
    l_statextn = statextn
    l_fl_mask = fl_mask
    l_nclip = nclip
    l_badpix = badpix
    l_calc_stddev = calc_stddev
    l_sci_ext = sci_ext
    l_dq_ext = dq_ext
    l_logfile = logfile
    l_verbose = verbose

    #### Set defaults
    debug = no
    status = 0
    isbad = no
    statx1 = 0
    statx2 = 0
    staty1 = 0
    staty2 = 0
    num_pix_threshold = 0.5 # If ratio of number of pixels used to calculate
                            # static and number of input pixels is less than
                            # this value a warning is printed
    sxcheck = yes
    sycheck = yes
    dqlist = ""

    ####
    # Temporary Files
    scilist = mktemp ("tmpscilist")
    tmpdqlist = mktemp ("tmpdqlist")

    ########
    # Here is where the actual work starts

    # Test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = gsaoi.logfile
       if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
          l_logfile = "gsaoi.log"
          printlog ("WARNING - GASTAT: Both gastat.logfile and \
              gsaoi.logfile fields are empty", l_logfile, l_verbose)
          printlog ("                     Using default file gsaoi.log", \
              l_logfile, l_verbose)
       }
    }

    printlog ("", l_logfile, l_verbose)
    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)
    gemdate (zone="local")
    printlog ("GASTAT -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    # Check for the same extnames being supplied for sci_ext and dq_ext
    if (l_sci_ext == l_dq_ext) {
        printlog ("ERROR - GASTAT: sci_ext and dq_ext cannot be the same",\
            l_logfile, verbose+)
        goto crash
    }

    # Check upper input
    if (!isindef(upper)) {

        tstring = str(upper)
        test = fscan(tstring,l_upper)

        if (test == 0) {
            pr_string = "upper"
            isbad = yes
        }
        pr_upper = str(l_upper)
    } else {
        l_upper = upper
        pr_upper = "INDEF"
    }

    # Check lower input
    if (!isindef(lower)) {

        tstring = str(lower)
        test = fscan(tstring,l_lower)

        if (test == 0) {
            pr_string = "lower"
            isbad = yes
        }
        pr_lower = str(l_lower)
    } else {
        l_lower = lower
        pr_lower = "INDEF"
    }

    # Check usigma input
    if (!isindef(usigma)) {

        tstring = str(usigma)
        test = fscan(tstring,l_usigma)

        if (test == 0) {
            pr_string = "usigma"
            isbad = yes
        } else if (l_nclip == 0) {
            printlog ("WARING - GASTAT: usigma is not INDEF and nclip \
                is zero.\n                 No sigma clipping will take place",\
                l_logfile, verbose+)
        }
        pr_usigma = str(l_usigma)
    } else {
        l_usigma = usigma
        pr_usigma = "INDEF"
    }

    # Check lsigma input
    if (!isindef(lsigma)) {

        tstring = str(lsigma)
        test = fscan(tstring,l_lsigma)

        if (test == 0) {
            pr_string = "lsigma"
            isbad = yes
        } else if (l_nclip == 0) {
            printlog ("WARING - GASTAT: usigma is not INDEF and nclip \
                is zero.\n                 No sigma clipping will take place",\
                l_logfile, verbose+)
        }
        pr_lsigma = str(l_lsigma)
    } else {
        l_lsigma = lsigma
        pr_lsigma = "INDEF"
    }

    # Record all input parameters relevant to this task only - other tasks will
    # print their inputs to log
    printlog ("GASTAT: Input parameters...", l_logfile, l_verbose)
    printlog ("    inimage     = "//l_inimage, l_logfile, l_verbose)
    printlog ("    stattype    = "//l_stattype, l_logfile, l_verbose)
    printlog ("    statsec     = "//l_statsec, l_logfile, l_verbose)
    printlog ("    statextn    = "//l_statextn, l_logfile, l_verbose)
    printlog ("    fl_mask     = "//l_fl_mask, l_logfile, l_verbose)
    printlog ("    badpix      = "//l_badpix, l_logfile, l_verbose)
    printlog ("    calc_stddev = "//l_calc_stddev, l_logfile, l_verbose)
    printlog ("    lower       = "//pr_lower, l_logfile, l_verbose)
    printlog ("    upper       = "//pr_upper, l_logfile, l_verbose)
    printlog ("    lsigma      = "//pr_lsigma, l_logfile, l_verbose)
    printlog ("    usigma      = "//pr_usigma, l_logfile, l_verbose)
    printlog ("    sci_ext     = "//l_sci_ext, l_logfile, l_verbose)
    printlog ("    dq_ext      = "//l_dq_ext, l_logfile, l_verbose)
    printlog ("    logfile     = "//l_logfile, l_logfile, l_verbose)
    printlog ("    verbose     = "//l_verbose, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    if (isbad) {
        printlog ("ERROR - GASTAT: Bad "//pr_string//" value supplied", \
            l_logfile, verbose+)
        goto crash
    }

    # Only one image allowed at time
    # Check for sections
    fparse (l_inimage)
    l_inimage = fparse.directory//fparse.root//fparse.extension
    if (fparse.ksection != "") {
        printlog ("WARNING - GASTAT: Setting statextn to extension \
           supplied as part of inimage.", \
           l_logfile, verbose+)
        l_statextn = fparse.ksection
    }

    # Check the input file exists and is MEF
    gemextn (l_inimage, check="exists,mef", process="expand", index="", \
        extname="", extversion="", ikparams="", omit="", replace="", \
        outfile="dev$null", logfile=l_logfile, glogpars="", verbose=no)

    if (gemextn.status != 0) {
        printlog ("ERROR- GASTAT: GEMEXTN returned a non-zero status"//\
            ". Exiting.", l_logfile, verbose+)
        goto crash
    } else {
        if (gemextn.fail_count > 0 || gemextn.count == 0) {
            # This is to catch an input error to gemextn
            printlog ("ERROR - GASTAT: GEMEXTN returned a count of 0 \
                and a fail_count of 0. Files not found. \
                Exiting", l_logfile, verbose+)
            goto crash
        }
    }

    # Check the input sci_ext and change parameters appropriately
    if (l_sci_ext == "" || stridx(" ",l_sci_ext) > 0) {
        l_sci_ext = ""
        l_extnversion = ""
        l_index = "1-"
    } else {
        l_extnversion = "1-"
        l_index = ""
    }

    # Determine the number of science data extensions
    gemextn (l_inimage, check="exists", process="expand", index=l_index, \
        extname=l_sci_ext, extversion=l_extnversion, ikparams="", omit="", \
        replace="", outfile=scilist, logfile=l_logfile, glogpars="", \
        verbose=no)

    # Check the output status
    if (gemextn.status != 0) {
        printlog ("ERROR- GASTAT: GEMEXTN returned a non-zero \
            status. Exiting.", l_logfile, verbose+)
        goto crash
    } else if (gemextn.count == 0 && l_sci_ext != "") {
        printlog ("ERROR- GASTAT: GEMEXTN returned a value of zero when \
            counting the number of "//\
            "\n               \""//l_sci_ext//"\" extensions in "//l_inimage,\
            l_logfile, verbose+)
        goto crash

    } else {
        nsci = gemextn.count
    }

    if (debug) printlog ("___nsci = "//nsci, l_logfile, verbose+)

    # If masking check the presence of DQ planes / BPM
    if (l_fl_mask) {

        # Cannot use dq_ext == "" to inspect the input files for DQ planes
        checkdq = yes
        if (l_dq_ext == "" || stridx(" ",l_dq_ext) > 0) {
            l_dq_ext = ""
            l_extnversion = ""
            l_index = "1-"
            checkdq = no
        } else {
            l_extnversion = "1-"
            l_index = ""
        }

        # DQ planes - if "", then not prepared - can't have DQ
        if (checkdq) {

            if (debug) printlog ("____Checking for DQs", l_logfile, verbose+)

            # Determine the number of science data extensions
            gemextn (l_inimage, check="exists", process="expand", \
                index=l_index, extname=l_dq_ext, extversion=l_extnversion, \
                ikparams="", \
                omit="params,section", replace="", \
                outfile=tmpdqlist, logfile=l_logfile, \
                glogpars="", verbose=no)

            # Check the output status
            if (gemextn.status != 0) {
                printlog ("ERROR- GASTAT: GEMEXTN returned a non-zero \
                    status. Exiting.", l_logfile, verbose+)
                goto crash
            } else if (gemextn.count != nsci && gemextn.count != 0) {
                printlog ("ERROR - GASTAT: The number of "//l_dq_ext//\
                    " planes "//gemextn.count//" does not match the number"//\
                    " of "//l_sci_ext//" planes in \""//l_inimage//"\".", \
                    l_logfile, verbose+)
                goto crash
            } else if (gemextn.count == 0) {

                printlog ("GASTAT: Cannot find any \'"//l_dq_ext//"\" \
                    extensions in \""//l_inimage//"\". Will check BPM "//\
                    "\n        \""//l_badpix//"\"", l_logfile, l_verbose)

            } else {
                printlog ("GASTAT: Will use "//l_dq_ext//" extensions "//\
                    "from \""//l_inimage//"\" for masking.", \
                    l_logfile, l_verbose)
                dqlist = tmpdqlist
            }
        }

        # Test the BPM
        if (dqlist == "") {

            if (debug) printlog ("____Checking BPM", l_logfile, verbose+)

            # Delete the temporary file from gemextn
            delete (tmpdqlist, verify-, >& "dev$null")

            # All BPMs should have a EXTNAMES of DQ; ordq_ext should be set
            if (l_dq_ext == "") {
                l_dq_ext = "DQ"
            }
            l_index = ""
            l_extnversion = "1-"

            if (l_badpix == "" || stridx (" ",l_badpix) > 0) {
                # Check for a bad BPM input
                printlog ("WARNING - GASTAT: badpix is either an empty string \
                    or contains whitespace"//\
                    "\n                  Switching off masking", \
                    l_logfile, verbose+)
                l_fl_mask = no
            } else if (strlwr(l_badpix) == "none") {
                # Check for "none" too
                printlog ("WARNING - GASTAT: badpix is \""//l_badpix//\
                    "\".\n                  Switching masking off", \
                    l_logfile, verbose+)
                l_fl_mask = no

            } else {

                # Determine if it exists and then the number of extensions
                gemextn (l_badpix, check="exists,mef", process="expand", \
                    index=l_index, extname=l_dq_ext, extversion=l_extnversion,\
                    ikparams="", omit="params,section", replace="", \
                    outfile=tmpdqlist, logfile=l_logfile, glogpars="", \
                    verbose=no)

                if (debug) printlog ("___l_badpix gemextn.count = "//\
                    gemextn.count, l_logfile, verbose+)

                # Check the output status
                if (gemextn.status != 0) {
                    printlog ("ERROR- GASTAT: GEMEXTN returned a non-zero \
                        status. Exiting.", l_logfile, verbose+)
                    goto crash
                } else if (gemextn.count != nsci && gemextn.count != 0) {
                     print (nsci//" "//gemextn.count)
                    printlog ("WARNING - GASTAT: The number of "//\
                        l_dq_ext//" planes "//gemextn.count//\
                        " does not match the number"//\
                        " of "//l_sci_ext//" planes in "//\
                        "\n                  \""//l_inimage//"\"."//\
                        "\n                  Switching off masking", \
                        l_logfile, verbose+)
                    l_fl_mask = no

                    # Delete the temporary file from gemextn
                    delete (tmpdqlist, verify-, >& "dev$null")

                } else if (gemextn.count == 0) {
                    printlog ("WARNING - GASTAT: \""//l_badpix//\
                        "\" does not exist"//\
                        "\n                  or is not an MEF or does not "//\
                        "contain "//l_dq_ext//" extensions."//\
                        "\n                  Switching off masking", \
                        l_logfile, verbose+)
                    l_fl_mask = no

                } else {
                    printlog ("GASTAT: Will use \""//l_dq_ext//\
                        "\" extensions "//\
                        "from \""//l_badpix//"\" for masking.", \
                        l_logfile, l_verbose)
                    dqlist = tmpdqlist
                }
            }
        } # End of use_badpix loop

        if (dqlist == "") {
            # Delete the temporary file from gemextn
            delete (tmpdqlist, verify-, >& "dev$null")
        }

    } # End of checking for DQ/BPM


    # Two lists one containing the science extensions and the masking
    # extensions. The latter may not exist

    # Determine / check the way the statistic is to be calculated
    # l_statextn - allowed values are "DETECTOR", "ARRAY",
    # "<l_sci_ext>,<extver>" and "index"
    # "DETECTOR means we have to tile.
    # Loop over all extensions for all types (to create output) but do things
    # appropriately so only the minimum gets done.

    l_stattype = "npix,"//l_stattype

    if (l_calc_stddev) {
        l_stattype = l_stattype//",stddev"
    }

    dotile = no
    allcalc = yes
    extn = 0
    isbad = no
    inextver = 1

    dq_file = ""
    dq_statsec = ""

    if (l_statextn == "DETECTOR") {
        dotile = yes
        allcalc = no

    } else if (l_statextn != "ARRAY") {

        allcalc = no

        # Check for "[" & "]"
        if (substr(l_statextn,1,1) == "[") {
            l_statextn = substr(l_statextn,2,strlen(l_statextn))
        }

        if (substr(l_statextn,strlen(l_statextn),strlen(l_statextn)) == "]") {
            l_statextn = substr(l_statextn,1,strlen(l_statextn)-1)
        }

        if (debug) printlog ("___statextn (1) =  "//l_statextn, \
            l_logfile, verbose+)

        # Parse l_statextn
        test = fscan (l_statextn, extn)
        if (test != 1) {
            # Check for extension name - must be the same as l_sci_ext
            first_charpos = strstr(l_sci_ext,l_statextn)

            if (debug) printlog ("___first_charpos = "//first_charpos, \
                l_logfile, verbose+)

            if (first_charpos > 0) {
                comma1_location = stridx(",",l_statextn)
                comma2_location = strldx(",",l_statextn)

                if (comma1_location == comma2_location) {
                    last_charpos = strlen(l_statextn)
                } else {
                    last_charpos = comma2_location - 1
                }

                first_charpos += strlen(l_sci_ext) + 1

                inextn = substr(l_statextn,first_charpos,last_charpos)

                if (debug) printlog ("___inextn (1) = "//inextn, \
                    l_logfile, verbose+)

                test = fscan (inextn,extn)
                if (test != 1) {
                    isbad = yes
                }

                inextver = int(inextn)

                if (l_sci_ext != "") {
                    inextn = l_sci_ext//","//inextn
                }

            } else {
                inextn = ""
                isbad = yes
            }
        } else {
            inextn = str(extn)
        }

        if (debug) printlog ("___inextn (2) = "//inextn, l_logfile, verbose+)

        if (l_fl_mask) {
            comma1_location = stridx(",",inextn)
            if (l_dq_ext == "") {
                comma1_location += 1
            }

            dqextn = "["//l_dq_ext//\
                substr(inextn,comma1_location,strlen(inextn))//"]"
        } else {
            dqextn = ""
        }

        inextn = "["//inextn//"]"

        if (debug) printlog ("___inextn (3) = "//inextn, l_logfile, verbose+)
        if (debug) printlog ("___dqextn = "//dqextn, l_logfile, verbose+)
    }

    if (isbad) {
        printlog ("ERROR - GASTAT: Bad statextn supplied", \
            l_logfile, verbose+)
        goto crash
    }
    isbad = no

    # Check the statsec input - have to do final check on statsec when about to
    # read the actual section - MS

    docalc = no
    tile_list = ""
    dqtile_list = ""
    dimscheck = yes

    printlog ("", l_logfile, l_verbose)

    l_outstat = ""
    l_stddev = ""

    # Loop over nsci
    for (i = 1; i <= nsci; i += 1) {

        # This assumes that the files are in order and that the extension
        # versions of the two lists are the same
        fields (scilist, fields="1", lines=i, quit_if_miss=no, \
            print_file_n=no) | scan (img)

        # Parse the statistics section for the current extension
        gemsecchk (img, l_statsec, logfile=l_logfile, verbose=debug)

        if (gemsecchk.status != 0) {
            printlog ("ERROR - GASTAT: GEMSECCHK reurned a non-zero status", \
                l_logfile, verbose+)
            goto crash
        } else {
            imgsec = gemsecchk.out_imgsect
        }

        dqimg = ""
        dq_statsec = ""
        if (l_fl_mask) {
            fields (dqlist, fields="1", lines=i, quit_if_miss=no, \
                print_file_n=no) | scan (dqimg)
        }

        # Check what needs to be done
        if (!dotile && !allcalc) {
            if (strstr(inextn,img) > 0) {
                docalc = yes
                stat_file = img
                dq_file = dqimg
            } else if (i == nsci) {
                printlog ("ERROR - GASTAT: Could not find extension "//\
                    inextn//" in "//l_inimage, l_logfile, verbose+)
                goto crash
            } else {
                dimscheck = no
            }
        } else if (dotile) {

            tile_list = tile_list//","//img
            if (l_fl_mask) {
                dqtile_list = dqtile_list//","//dqimg
            }

        } else {
            docalc = yes
            stat_file = img
            dq_file = dqimg
        }

        # Check the dimensions of img against statsec
        if (dimscheck) {

            gadimschk (img, section=imgsec, chkimage=dqimg, \
                key_check="CCDSEC", logfile=l_logfile, verbose-)

            if (gadimschk.status != 0) {
                printlog ("ERROR - GASTAT: GADIMSCHK returned a non-zero \
                    status. Exiting.", l_logfile, verbose+)
                goto crash
            } else if (l_fl_mask) {
                dq_statsec = substr(gadimschk.out_chkimage,\
                    strlen(dqimg)+1,strlen(gadimschk.out_chkimage))
            }

            # Adjust tiling expression
            if (dotile) {
                dqtile_list = dqtile_list//dq_statsec
            }

            # Adjust tiling expression
            if (dotile) {
                tile_list = tile_list//imgsec
            }

        } # End of dimension checks

        # Tile when needed
        if (i == nsci && dotile) {
            docalc = yes

            fparse (img, verbose-)
            stat_file = fparse.root//"_tmptile"//mktemp("1")//".fits"

            tile_list = substr(tile_list,2,strlen(tile_list))

            if (debug) {
                printlog ("____tile_list: "//tile_list, \
                    l_logfile, verbose+)
            }

            imjoin (tile_list, stat_file, join_dimension=1, pixtype="r", \
                verbose=no)

            if (!imaccess(stat_file)) {
                printlog ("ERROR - Cannot access temporary tiled image", \
                    l_logfile, verbose+)
                goto crash
            }

            # Now that the appropriate sections have been extracted and joined
            # reset sections to *,*
            imgsec = "[*,*]"

            if (l_fl_mask) {

                dq_file = fparse.root//"_tmpdqtile"//mktemp("1")//".fits"

                dqtile_list = substr(dqtile_list,2,strlen(dqtile_list))

                if (debug) {
                    printlog ("____dqtile_list: "//dqtile_list, \
                        l_logfile, verbose+)
                }

                imjoin (dqtile_list, dq_file, \
                    join_dimension=1, pixtype="r", verbose=no)

                if (!imaccess(dq_file)) {
                    printlog ("ERROR - Cannot access temporary tiled \
                        dq image", l_logfile, verbose+)
                    goto crash
                }

                # Now that the appropriate sections have been extracted and
                # joined reset sections to *,*
                dq_statsec = "[*,*]"
            }

        }

        if (docalc) {

            npix = 0

            # Determine the area of the input image
            # The naxis keywords should have been checked already
            keypar (stat_file, "i_naxis1", silent+)
            naxis1 = int(keypar.value)

            keypar (stat_file, "i_naxis2", silent+)
            naxis2 = int(keypar.value)

            # The number of input pixels in the image
            inpix = naxis1 * naxis2

            if (debug) {
                printlog ("____stat_file//imgsec: "//stat_file//imgsec, \
                    l_logfile, verbose+)
                printlog ("____dq_file//dq_statsec: "//dq_file//dq_statsec, \
                    l_logfile, verbose+)
            }

            # Determine the statistic for this image
            statval = -99999.0
            tmpstddev = -99999.0

            # Do the calculation
            mimstatistics (stat_file//imgsec, imasks=dq_file//dq_statsec,\
                field=l_stattype, lower=l_lower, upper=l_upper, \
                lsigma=l_lsigma, usigma=l_usigma, nclip=l_nclip, \
                binwidth=0.1, format-, cache-) | scan (val1)

            if (debug) {
                printlog ("____stattype: "//l_stattype//" val1: "//\
                    val1, l_logfile, verbose+)

                mimstatistics (stat_file//imgsec, \
                    imasks=dq_file//dq_statsec,\
                    field=l_stattype, lower=l_lower, upper=l_upper, \
                    lsigma=l_lsigma, usigma=l_usigma, nclip=l_nclip, \
                    binwidth=0.1, format-, cache-)
            }

            if (l_calc_stddev) {
                print (val1) | scan (npix, statval, tmpstddev)

                # Check a value was returned
                if (tmpstddev == -99999.0) {
                    printlog ("ERROR - GASTAT: Could not "//\
                        "determine a value for the "//\
                        "stddev of "//stat_file, \
                        l_logfile, verbose+)
                    goto crash
                } else {
                    l_stddev = l_stddev//","//str(tmpstddev)
                }

            } else {
                print (val1) | scan (npix, statval)
            }

            # Calculate the ratio of pixels used to number input
            pix_ratio = ((1.0 * npix) / inpix)

            if (debug) {
                printlog ("____npix: "//npix//" inpix: "//inpix//" ratio: "//\
                    pix_ratio, l_logfile, verbose+)
            }

            if (npix == 0) {
                printlog ("ERROR - GASTAT: All pixels masked out. Please try \
                    again with fl_mask=no.", l_logfile, verbose+)
                goto crash
            } else if (pix_ratio < num_pix_threshold) {
                printlog ("WARNING - GASTAT: Less than "//\
                    str(int(100*num_pix_threshold))//\
                    "% of the input pixels were used to calculate "//\
                    "statistic\n", l_logfile, verbose+)
            }

            # Check a value was returned
            if (statval == -99999.0) {
                printlog ("ERROR - GASTAT: Could not "//\
                    "determine a value for the "//\
                    l_stattype//" of "//stat_file, \
                    l_logfile, verbose+)
                goto crash
            } else {
                printlog ("GASTAT: Calculated "//\
                    substr(l_stattype,stridx(",",l_stattype)+1,\
                    strlen(l_stattype))//\
                    " for "//stat_file//\
                    ":\n            "//\
                    substr(val1,stridx(" ",val1)+1,strlen(val1)), \
                    l_logfile, l_verbose)
                l_outstat = l_outstat//","//str(statval)
            }

            # Create correct output if tiling or specific extension requested
            if (!allcalc) {
                # Loop over nsci - 1
                for (j = 1; j <= (nsci - 1); j += 1) {
                    l_outstat = l_outstat//","//str(statval)
                    if (l_calc_stddev) {
                        l_stddev = l_stddev//","//str(tmpstddev)
                    }
                }

                if (dotile) {
                    imdelete (stat_file//","//dq_file, verify-, >& "dev$null")
                } else {
                    break
                }
            }
        }

    } # End of for loop

    # Remove the leading comma
    l_outstat = substr(l_outstat,2,strlen(l_outstat))
    printlog ("\nGASTAT: Output statistic - \""//l_outstat//"\"", \
        l_logfile, verbose-)

    if (l_stddev != "") {
        l_stddev = substr(l_stddev,2,strlen(l_stddev))
        printlog ("GASTAT: Output stddev - \""//l_stddev//"\"", \
            l_logfile, verbose-)
    }

    outstat = l_outstat
    stddev = l_stddev
    goto clean

    #--------------------------------------------------------------------------

crash:

    # Exit with error subroutine
    status = 1
    outstat = ""
    stddev = ""

clean:

    delete (scilist//","//dqlist//","//tmpdqlist, verify-, >& "dev$null")

    # Print finish time
    gemdate(zone="local")
    printlog ("\nGASTAT -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("\nGASTAT -- Exit status: GOOD", l_logfile, l_verbose)

    } else {
        printlog ("\nGASTAT -- Exit status: ERROR", l_logfile, l_verbose)
    }

    if (status != 0) {
        printlog ("       -- Please read the logfile \""//l_logfile//\
                "\" for more information.", l_logfile, l_verbose)
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)

end
