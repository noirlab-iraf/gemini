# Copyright(c) 2012-2013 Association of Universities for Research in Astronomy, Inc.

procedure gaimchk (inimages, outlist)

# This tasks does all of the standard input checks that are done in a script
# but all in one place. The output list contains a list of temporary files that
# contain files names that match the requested criteria. If not matches are
# found for a particular file then a tmpNotFound list is created.

# The checks performed are dependent on the flags set and the keywords
# supplied. See the scripts that call this for examples of its usage.

# It also does the smart parsing of input parameters for the inimages in all
# tasks. It will accept a comma separated lists (of lists too), lists, image
# number ranges, wild cards and files names.

# key_allowed is hierarchical!

# Currently only works for GSAOI - I don't think it would be much effort to
# expand to others. Need to, at least, change the check for the number of
# extensions. In addition will need to update print statements!

######################################################

char    inimages        {prompt="Input GSAOI images or list"}
char    outlist         {prompt="Output list containing output lists"}
char    rawpath         {"",prompt="Path for input raw images"}
char    rootname        {"",prompt="Root name if supplying image number; blank for today's UT"}
#char    instrument      {"GSAOI",enum="GSAOI",prompt="Required INSTRUME value"}
char    obstype_allowed {"",prompt="Allowed OBSTYPE keyword values"}
char    object_allowed  {"",prompt="Allowed OBJECT keyword values"}
char    key_allowed     {"",prompt="Check PHU for the existence of these keywords"}
char    key_forbidden   {"",prompt="PHU must not contain this keyword"}
char    key_exists      {"",prompt="PHU must contain this keyword"}
bool    fl_prep_check   {yes,prompt="Check for an already prepared version of the input file?"}
char    gaprep_pref     {"g",prompt="GAPREPARE prefix"}
bool    fl_redux_check  {no,prompt="Check for an already GAREDUCEd version of the input file?"}
char    garedux_pref    {"r",prompt="GAREDUCE prefix"}
bool    fl_fail         {no,prompt="Fail if fl_prep_check or fl_redux_check is true?"}
bool    fl_out_check    {no,prompt="Fail if file exists already?"}
bool    fl_vardq_check  {no,prompt="Check for VAR / DQ extensions?"}
char    sci_ext         {"SCI",prompt="Extension name for science data planes"}
char    var_ext         {"VAR",prompt="Extension name for variance planes"}
char    dq_ext          {"DQ",prompt="Extension name for data quality planes"}
char    logfile         {"",prompt="Logfile"}
bool    verbose         {yes,prompt="Verbose?"}
int     status          {0,prompt="Exit status (0=good)"}
struct  *scanfile       {"",prompt="Internal use only"}

######################################################

begin

    ########
    # Declare and set local variables; set default values; initiate temporary
    # files

    ####
    # Variable declaration

    char l_inimages, l_rawpath, l_logfile, l_rootname, orig_l_rawpath
    char l_key_allowed, l_key_forbidden, l_key_exists, l_obstype_allowed
    char l_object_allowed, l_instrument, l_key_prepared
    char l_gaprep_pref, l_garedux_pref, l_sci_ext, l_var_ext, l_dq_ext
    char pathtest, rphend, inlist, t_string, filelist, utdate, tmplist
    char l_range, t2_string, t3_string, img_list
    char inimg, inname, inphu, word_list, inword_list, tmp2list, envvar
    char l_outlist, l_orig_outlist, current_outlist, not_found_list, nextlist
    char l_key_mosaiced

    int i, atlocation, comma_pos, lcount, last_char, test, test2, t_num
    int dash_test, dash_ltest, xtest, xltest, t_length, allowed_counter
    int dotpos1, dotpos2, val_sci, nsci

    bool l_fl_gaprep_check, l_fl_redux_check, l_verbose, l_fl_vardq_check
    bool must_have, isallowed, not_allowed, wr_NF_to_outlist, write_to_outlist
    bool skip_vardq, l_fl_fail, l_fl_out_check, debug, empty_allowed, dotrim
    bool can_be_mosaiced, bad_num_extensions

    ####
    # Set local variables
    l_inimages = inimages
    l_rawpath = rawpath
    l_rootname = rootname
#    l_instrument = instrument # uncomment this andthe instrument parameter if
    # to be used with another instrument
    l_instrument = "GSAOI"
    l_obstype_allowed = obstype_allowed
    l_object_allowed = object_allowed
    l_key_allowed = key_allowed
    l_key_forbidden = key_forbidden
    l_key_exists = key_exists
    l_fl_gaprep_check = fl_prep_check
    l_gaprep_pref = gaprep_pref
    l_fl_redux_check = fl_redux_check
    l_garedux_pref = garedux_pref
    l_fl_fail = fl_fail
    l_fl_out_check = fl_out_check
    l_fl_vardq_check = fl_vardq_check
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_orig_outlist = outlist
    l_logfile = logfile
    l_verbose = verbose

    ####
    # Set default values
    debug = no
    status = 0
    scanfile = ""
    can_be_mosaiced = no
    # Set the valid number of science extensions
    ##M If expanding to use with other instruments set default to 1 and change
    ##M according to instrument. For GMOS a little more complicated (possibly
    ##M allow use to set? Something similar for l_key_prepared

    if (l_instrument == "GSAOI") {
        l_key_prepared = "GAPREPAR"
        val_sci = 4
        can_be_mosaiced = yes
        l_key_mosaiced = "GAMOSAIC"
    }

    ####
    # Set temporary files
    filelist = mktemp("tmpfilelist")
    tmplist = mktemp ("tmplist")
    tmp2list = mktemp ("tmp2list")
    l_outlist = mktemp ("tmpOUTlist")
    nextlist = mktemp("tmpnextlist")

    ########
    # Here is where the actual work starts

    # Test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = gsaoi.logfile
       if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
          l_logfile = "gsaoi.log"
          printlog ("WARNING - GAIMCHK: Both gaflat.logfile and \
              gsaoi.logfile fields are empty", l_logfile, l_verbose)
          printlog ("                     Using default file gsaoi.log", \
              l_logfile, l_verbose)
       }
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)
    gemdate (zone="local")
    printlog ("GAIMCHK -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    # Record all input parameters relavent to this task only - other tasks will
    # print their inputs to log
    printlog ("GAIMCHK: Input parameters...", l_logfile, l_verbose)
    printlog ("    inimages        = "//l_inimages, l_logfile, l_verbose)
    printlog ("    rawpath         = "//l_rawpath, l_logfile, l_verbose)
    printlog ("    rootname        = "//l_rootname, l_logfile, l_verbose)
    printlog ("    obstype_allowed = "//l_obstype_allowed, \
        l_logfile, l_verbose)
    printlog ("    object_allowed  = "//l_object_allowed, l_logfile, l_verbose)
    printlog ("    key_allowed     = "//l_key_allowed, l_logfile, l_verbose)
    printlog ("    key_forbidden   = "//l_key_forbidden, l_logfile, l_verbose)
    printlog ("    key_exists      = "//l_key_exists, l_logfile, l_verbose)
    printlog ("    fl_gaprep_check = "//l_fl_gaprep_check, \
        l_logfile, l_verbose)
    printlog ("    gaprep_pref     = "//l_gaprep_pref, l_logfile, l_verbose)
    printlog ("    fl_redux_check  = "//l_fl_redux_check, l_logfile, l_verbose)
    printlog ("    garedux_pref    = "//l_garedux_pref, l_logfile, l_verbose)
    printlog ("    fl_fail         = "//l_fl_fail, l_logfile, l_verbose)
    printlog ("    l_sci_ext       = "//l_sci_ext, l_logfile, l_verbose)
    printlog ("    l_var_ext       = "//l_var_ext, l_logfile, l_verbose)
    printlog ("    l_dq_ext        = "//l_dq_ext, l_logfile, l_verbose)
    printlog ("    outlist         = "//l_orig_outlist, l_logfile, l_verbose)
    printlog ("    logfile         = "//l_logfile, l_logfile, l_verbose)
    printlog ("    verbose         = "//l_verbose, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    if (l_orig_outlist == "" || stridx (" ",l_orig_outlist) > 0) {
        printlog ("ERROR - GAIMCHK: outlist is not specified", \
            l_logfile, verbose+)
        goto crash
    } else if (access(l_orig_outlist)) {
        printlog ("ERROR - GAIMCHK: outlist \""//l_orig_outlist//"\" already \
            exists. Exiting.", l_logfile, verbose+)
        goto crash
    }

    if (l_obstype_allowed == "" ) {
        l_obstype_allowed = "any"
    } else {
        l_obstype_allowed = strlwr(l_obstype_allowed)
    }

    if (l_object_allowed == "" ) {
        l_object_allowed = "any"
    } else {
        l_object_allowed = strlwr(l_object_allowed)
    }

    if (l_fl_out_check && (l_fl_redux_check || fl_prep_check)) {
        printlog ("WARNING - GAIMCHK: Conflicting checks set. Setting \
            fl_out_check=no.", l_logfile, verbose+)
        l_fl_out_check = no
    }

    if (l_fl_fail && !l_fl_redux_check && !fl_prep_check) {
        printlog ("WARNING - GAIMCHK: Conflicting checks set. Setting \
            fl_fail=no.", l_logfile, verbose+)
        l_fl_fail = no
    }

    # Check that the rawpath has a trailing slash and is a valid entry
    if (l_rawpath != "") {
        rphend = substr(l_rawpath,strlen(l_rawpath),strlen(l_rawpath))

        if (rphend == "$") {
            envvar = substr(l_rawpath,1,strlen(l_rawpath)-1)
            if (defvar (envvar)) {
                show (substr(l_rawpath,1,strlen(l_rawpath)-1)) | \
                    scan (pathtest)
                rphend = substr(pathtest,strlen(pathtest),strlen(pathtest))
            }
        }

        if (rphend != "/") {
            l_rawpath = l_rawpath//"/"
        }

        if (!access(l_rawpath)) {
            printlog ("ERROR - GAIMCHK: Cannot access rawpath: "//l_rawpath, \
                l_logfile, verbose+)
            goto crash
        }
    }
    printlog ("GAIMCHK: rawpath after checks is \""//l_rawpath//"\"\n", \
        l_logfile, l_verbose)

    # Require the path on disk - this makes life easier later on with section
    # calls , store original for later on output
    orig_l_rawpath = l_rawpath
    l_rawpath = osfn(orig_l_rawpath)

    # Load up input name list: @list, wildcards (including in comma separated),
    #                          comma-separated and image number ranges

    # The logic for this will be:
    #     if "@" is present, this is an @list
    #     if "*" or "?" are present, check for comma separated, then expand the
    #         wild cards
    #     if no wild cards, comma separated check firt entry to see if it's a
    #         file and it exists. If it exists it's a comma separated list. If
    #         it doesn't exit it's a range of file numbers -  form frame name
    #         using l_rootname (UT if not specified)
    #     if no @list, wild cards, or comma separated lists, check if entry is
    #         a number. If so, form frame name using l_rootname (UT if not
    #         specified). If not a number assume it's a file name

    if ((l_inimages == "") || (stridx(" ", l_inimages) > 0)) {
        printlog ("ERROR - GAIMCHK: No input files specified", \
            l_logfile, verbose+)
        goto crash
    }

    # Parse any commas from l_inimages
    comma_pos = 1
    l_range = ""
    img_list = l_inimages
    i = 0

    while (comma_pos != 0) {
        i += 1

        comma_pos = stridx (",", img_list)
        if (comma_pos == 0) {
            last_char = strlen(img_list)
        } else {
            last_char = comma_pos - 1
        }

        t_string = substr(img_list, 1, last_char)
        t_length = strlen(t_string)

        if (comma_pos != 0) {
            img_list = substr(img_list, comma_pos+1, strlen(img_list))
        }

        if (debug) printlog  ("____"//i//" t_string: "//t_string, \
            l_logfile, verbose+)
        if (debug) printlog  ("____"//i//" img_list: "//img_list, \
            l_logfile, verbose+)

        atlocation = stridx("@",t_string)
        if (atlocation > 0) {
            if (debug) printlog ("____In at list", l_logfile, verbose+)

            # Test for an @filelist
            inlist = substr(t_string,atlocation+1,strlen(t_string))
            if (!access(inlist)) {
                printlog ("ERROR - GAIMCHK: Input list "//inlist//\
                    " does not exist",l_logfile, yes)
                goto crash
            }
            sections (l_rawpath//"//"//t_string, option="fullname", >> \
                filelist)

        } else if (stridx("*?",t_string) > 0) {
            if (debug) printlog ("____In wildcard loop", l_logfile, verbose+)
            # Test for wild cards

            if (t_length >= 6) {
                if (substr(t_string,t_length - 4, t_length) != ".fits") {
                    # Add ".fits" otherwise sections will not find it
                    t_string = t_string//".fits"
                }

            } else {
                # File extension not supplied!
                t_string = t_string//".fits"
            }
            if (debug) printlog  ("____t_string updated: "//t_string, \
                l_logfile, verbose+)

            sections (l_rawpath//t_string, option="fullname", >> \
                filelist)
        } else {
            if (debug) printlog ("___Checking for ranges/files", \
                l_logfile, verbose+)

            # Check for everything else
            t_num = -99
            test = fscan (t_string, t_num)

            if (test == 1 && (str(t_num) == t_string)) {
                # It's a number - append it to l_range for later
                l_range = l_range//","//t_num
                if (debug) printlog ("____t_string is a number", \
                    l_logfile, verbose+)

            } else {
                if (debug) printlog ("____t_string is NOT a number - testing \
                    for range", l_logfile, verbose+)

                dash_test = stridx("-",t_string)
                dash_ltest = strldx("-",t_string)
                # Check if it's a range still!
                # For gemlist to work must contain only one "-" for a range,
                # can also contain an "x"
                if ((dash_test == dash_ltest) && (dash_test > 0)) {
                    # It contains only 1 "-"

                    if (dash_test == 1) {
                        t2_string = substr(t_string,2,strlen(t_string))
                    } else if (dash_test != strlen(t_string)) {
                        t2_string = substr(t_string,1,dash_test-1)//\
                            substr(t_string,dash_test+1,strlen(t_string))
                    } else {
                        # This means it's not a valid range
                        t2_string = substr(t_string,1,dash_test-1)
                    }
                    if (debug) printlog ("____t2_string: "//t2_string, \
                        l_logfile, verbose+)

                    xtest = stridx("x",t2_string)
                    xltest = strldx("x",t2_string)
                    # Check if it contains only 1 x
                    if ((xtest == xltest) && (xtest > 0)) {
                        if (xtest == 1) {
                            t3_string = substr(t2_string,2,strlen(t2_string))
                        } else if (xtest != strlen(t_string)) {
                            t3_string = substr(t2_string,1,xtest-1)//\
                                substr(t2_string,xtest+1,strlen(t2_string))
                        } else {
                            # It's not a valid range
                            t3_string = substr(t2_string,1,xtest-1)
                        }
                    } else {
                        t3_string = t2_string
                    }

                    if (debug) printlog ("____t3_string: "//t3_string, \
                        l_logfile, verbose+)

                    # Test if it's a number
                    test2 = fscan (t3_string, t_num)
                    if (test2 == 1) {
                        l_range = l_range//","//t_string
                        if (debug) printlog ("____t3_string is a number", \
                            l_logfile, verbose+)
                    } else {
                        # It's a string that isn't a range...
                        if (debug) printlog ("____t3_string is NOT a number", \
                            l_logfile, verbose+)

                        sections (l_rawpath//"//"//t3_string, \
                            option="fullname", >> filelist)
                    }
                } else {
                    if (debug) printlog ("____t_string is NOT a range", \
                        l_logfile, verbose+)

                    # It's a string that isn't a range...
                    sections (l_rawpath//"//"//t_string, \
                        option="fullname", >> filelist)
                }
            } # End of test if t_string is a number
        } # End of else checking for filenames / ranges
    } # End of while comma_pos > 0

    if (l_range != "") {
        l_range = substr(l_range,2,strlen(l_range))
        if (debug) printlog ("____l_range: "//l_range, l_logfile, verbose+)

        if ((l_rootname == "") || (stridx(" ", l_rootname) > 0)) {
            getfakeUT()
            utdate = getfakeUT.fakeUT
            printlog ("GAIMCHK: Determining today's UT date: "//\
                utdate, l_logfile, l_verbose)
            l_rootname = "S"//utdate//"S"
        }

        printlog ("GAIMCHK: rootname is "//l_rootname, \
            l_logfile, l_verbose)

        l_rootname = l_rawpath//l_rootname
        gemlist (root=l_rootname, range=l_range, >> filelist)
    }

    # Check if there are any images in the list
    lcount = 0
    if (access(filelist)) {
        count (filelist) | scan (lcount)
    }

    if (lcount == 0) {
        printlog ("ERROR - GAIMCHK: No input images supplied. \
            Please try again.", l_logfile, verbose+)
        goto crash
    }

    # Create the input for the next set of checks
    sort (filelist, column=1, ignore_white=no, numeric_sort=no, \
        reverse_sort=no, >> tmplist)

    # Loop over input file and perform validation checks
    scanfile = tmplist
    while (fscan(scanfile, inimg) != EOF) {

        # Parse the file name, set input name and set output name
        fparse (inimg, verbose-)
        inname = fparse.root//fparse.extension

        # Reset inimg to contain input path not the full system path
        inimg = orig_l_rawpath//substr(inimg,strlen(l_rawpath)+1,strlen(inimg))
        inphu = inimg//"[0]"

        printlog ("GAIMCHK: Inspecting status of file "//inname//".", \
            l_logfile, verbose=debug)

        # Do the checks to see if prepare or reduced files already exist...
        if (l_fl_redux_check && (imaccess(l_garedux_pref//inname) || \
            imaccess(l_garedux_pref//l_gaprep_pref//inname))) {

            if (imaccess(l_garedux_pref//inname)) {
                inimg = l_garedux_pref//inname
            } else {
                inimg = l_garedux_pref//l_gaprep_pref//inname
            }

            if (l_fl_fail) {
                printlog ("ERROR - GAIMCHK: GAREDUCEd version "//\
                    inimg//" already exists", \
                    l_logfile, verbose+)
                goto crash
            }

            # Double check this is truely the case
            keypar (inimg//"[0]", "GAREDUCE", silent+)
            if (keypar.found) {
                printlog ("WARNING - GAIMCHK: Using / checking the already \
                    GAREDUCEd version "//inimg, \
                    l_logfile, verbose+)
            }

        } else if (l_fl_gaprep_check && imaccess(l_gaprep_pref//inname)) {

            if (l_fl_fail) {
                printlog ("ERROR - GAIMCHK: GAPREPAREd version "//\
                    l_gaprep_pref//inname//" already exists", \
                    l_logfile, verbose+)
                goto crash
            }

            # Double check this is truely the case
            keypar (l_gaprep_pref//inname//"[0]", "GAPREPAR", silent+)
            if (keypar.found) {

                printlog ("WARNING - GAIMCHK: Using / checking the already \
                    GAPREPAREd version "//l_gaprep_pref//inname, \
                    l_logfile, verbose+)

                inimg = l_gaprep_pref//inname
            }

        } else {

            # Check file exists
            if (!imaccess(inimg)) {
                if (!l_fl_out_check) {
                    printlog ("WARNING - GAIMCHK: Image "//inname//\
                        " does not exist. Skipping this image.", \
                        l_logfile, l_verbose)
                    goto NEXTIMAGE
                } else {
                    print (inimg, >> tmp2list)
                    goto NEXTIMAGE
                }

            } else if (l_fl_out_check) {
                printlog ("ERROR - GAIMCHK: Image "//inname//\
                    " already exists. Exiting.", \
                    l_logfile, verbose+)
                goto crash
            }

            # Check file is MEF
            gemextn (inimg, check="mef", process="expand", index="1-", \
                extname="", extversion="", ikparams="", omit="", replace="", \
                outfile="dev$null", logfile=l_logfile, glogpars="", \
                verbose=debug)

            nsci = gemextn.count

            if (gemextn.status != 0) {
                printlog ("ERROR- GAIMCHK: GEMEXTN returned a non-zero \
                    status. Exiting.", l_logfile, verbose+)
                goto crash
            } else {
                if (gemextn.fail_count > 0) {
                    printlog ("WARNING - GAIMCHK: Image "//inname//\
                        " is not an MEF. Skipping this image.", \
                        l_logfile, verbose+)
                    goto NEXTIMAGE
                } else if (gemextn.count == 0) {
                    # This is to catch an input error to gemextn
                    printlog ("ERROR - GAIMCHK: GEMEXTN returned a count of 0 \
                        and a fail_count of 0. Exiting", l_logfile, verbose+)
                    goto crash
                }
            }

            # Check it's the correct instrument
            keypar (inphu, "INSTRUME", silent+)
            if (keypar.value != l_instrument) {
                printlog ("WARNING - GAIMCHK: \""//inname//"\" is not a "//\
                    l_instrument//" image. Skipping this image", \
                    l_logfile, verbose+)
                goto NEXTIMAGE
            } else {

                # Check it has the correct number of extensions
                keypar (inphu, l_key_prepared, silent+)
                if (keypar.found) {
                    keypar (inphu, "NSCIEXT", silent+)
                    if (!keypar.found) {
                        printlog ("WARNING GAIMCK: \""//inname//\
                            " is prepared but doesn't contain the NSCIEXT"//\
                            " keyword. Skipping this image", \
                            l_logfile, verbose+)

                        goto NEXTIMAGE
                    } else {
                        nsci = int (keypar.value)
                    }
                }

                if (nsci != val_sci) {
                    bad_num_extensions = yes
                    # Check in case it's been mosaiced
                    if (can_be_mosaiced) {
                        keypar (inphu, l_key_mosaiced, silent+)

                        if (keypar.found) {
                            if (nsci == 1) {
                                bad_num_extensions = no
                            }
                        }
                    }
                } else {
                    bad_num_extensions = no
                }

                if (bad_num_extensions) {
                    printlog ("WARNING - GAIMCHK: \""//inname//"\" is a "//\
                        l_instrument//" image but does not contain the\n"//\
                        "                   correct number of extensions. \
                        Skipping this image", \
                        l_logfile, verbose+)
                    goto NEXTIMAGE
                }
            }

        } # End of check for prepare / reduced / inputs

        # Print the names into a new list
        print (inimg, >> tmp2list)

NEXTIMAGE:
    } # End of loop over full list of input images detected

    lcount = 0
    if (access(tmp2list)) {
        count (tmp2list) | scan (lcount)
    }

    if (lcount == 0) {
        printlog ("ERROR - GAIMCHK: No input images pass first round of \
            checks. Please try again.", l_logfile, verbose+)
        goto crash
    }

    # Find only the unique names in the list - possible if prepared versions
    # etc already exist and check for a previously prepared version
    sort (tmp2list, column=1, ignore_white=no, numeric_sort=no, \
        reverse_sort=no) | unique ("STDIN", >> nextlist)

    delete (tmplist//","//tmp2list, verify-, >& "dev$null")

    #Set up default output list
    not_found_list = "tmpNotFound_"//tmplist//".lis"

    # Loop over unique input file list and perform keyword checks/vardq checks
    scanfile = nextlist
    while (fscan(scanfile, inimg) != EOF) {

        # Parse the file name, set input name and set output name
        fparse (inimg, verbose-)
        inname = fparse.root//fparse.extension
        inphu = inimg//"[0]"

        if (l_fl_out_check) {
            if (!access(not_found_list)) {
               print (not_found_list, >> l_outlist)
            }
            print (inimg//" "//inname, >> not_found_list)
            goto NEXT_IMAGE
        }

        printlog ("GAIMCHK: Inspecting keywords in file "//inname//".", \
            l_logfile, verbose-)

        # Check for any allowed OBSTYPE keywpords / OBJECT keyword values
        if (l_obstype_allowed != "any") {
            keypar (inphu, "OBSTYPE", silent+)
            if (strstr(strlwr(keypar.value),l_obstype_allowed) == 0) {
                printlog ("WARNING - GAIMCHK: OBSTYPE value for \""//inname//\
                    "\" is not in the allowed OBSTYPE list."//\
                    "\n                   Skipping this image.", \
                    l_logfile, verbose+)
                goto NEXT_IMAGE
            }
        }

        # Check for any allowed OBSTYPE keywpords / OBJECT keyword values
        if (l_object_allowed != "any") {
            keypar (inphu, "OBJECT", silent+)
            if (strstr(strlwr(keypar.value),l_object_allowed) == 0) {
                printlog ("WARNING - GAIMCHK: OBJECT value for \""//inname//\
                    "\" is not in the allowed OBJECT list."//\
                    "\n                   Skipping this image.", \
                    l_logfile, verbose+)
                goto NEXT_IMAGE
            }
        }

        # Form a long list of the three keyword lists to loop over
        inword_list = ""
        if (l_key_forbidden != "") {
            inword_list = l_key_forbidden
        }

        if (l_key_exists != "") {
            inword_list = inword_list//","//l_key_exists
        }

        if (l_key_allowed != "") {
            inword_list = inword_list//","//l_key_allowed
            empty_allowed = no
        } else {
            empty_allowed = yes
        }

        if (debug) printlog ("____inword_list: "//inword_list, \
            l_logfile, verbose+)

        # Only perform these checks if neccessary
        if (inword_list != "") {

            if (substr(inword_list,1,1) == ",") {
                inword_list = substr(inword_list,2,strlen(inword_list))
            }

            # Output list when all 'allowed' keywords not found
            if (access(not_found_list)) {
                wr_NF_to_outlist = no
            } else {
                wr_NF_to_outlist = yes
            }

            allowed_counter = 0
            i = 0
            comma_pos = 1
            word_list = inword_list
            while (comma_pos != 0) {

                must_have = no
                isallowed = no
                not_allowed = no
                skip_vardq = no

                i += 1
                comma_pos = stridx (",", word_list)

                if (comma_pos == 0) {
                    last_char = strlen(word_list)
                } else {
                    last_char = comma_pos - 1
                }

                t_string = substr(word_list, 1, last_char)
                t_length = strlen(t_string)

                if (comma_pos != 0) {
                    word_list = substr(word_list,comma_pos+1,strlen(word_list))
                }

                # Set up output list for this keyword
                current_outlist = "tmp"//t_string//"_"//tmplist//".lis"

                # Check for access to the current_file set flag to write to
                # out_list or not
                if (access(current_outlist)) {
                    write_to_outlist = no
                } else {
                    write_to_outlist = yes
                }

                # Set flags according to which list they came from
                if (strstr(t_string,l_key_forbidden) > 0) {
                    not_allowed = yes
                }

                if (strstr(t_string,l_key_exists) > 0) {
                    must_have = yes
                }

                if (strstr(t_string,l_key_allowed) > 0) {
                    isallowed = yes
                }

                if (debug) printlog  ("____"//i//" t_string: \""//\
                    t_string//"\"", l_logfile, verbose+)
                if (debug) printlog  ("____"//i//" word_list: "//word_list, \
                    l_logfile, verbose+)
                if (debug) printlog  ("____"//i//" isallowed: "//isallowed//\
                    " must_have: "//must_have//" not_allowed: "//not_allowed//\
                    " empty_allowed: "//empty_allowed, \
                    l_logfile, verbose+)

                # Check flags - checking for multiple inputs of the same
                # keyword
                if (not_allowed && (must_have || isallowed)) {
                    printlog ("ERROR - GAIMGCHK: Inconsistant \
                        keywords set. Exiting.", \
                        l_logfile, verbose+)
                    goto crash
                }

                if (!must_have && !not_allowed && !isallowed) {
                    printlog ("ERROR - GAIMGCHK: Inconsistant key_exists, \
                        key_allowed and key_forbidden keywords set. Exiting.",\
                            l_logfile, verbose+)
                        goto crash
                }

                # Check the current keyword - act appropriately
                keypar (inphu, t_string, silent+)
                if (keypar.found) {

                    if (not_allowed) {
                        if (debug) printlog ("____found and not allowed", \
                            l_logfile, verbose+)

                        printlog ("GAIMCHK: Image "//inname//\
                            " already contains the "//t_string//" keyword."//\
                            "\n         Skipping this image.", \
                            l_logfile, l_verbose)
                        skip_vardq = yes
                        goto KEY_NEXTIMAGE

                    } else if (must_have && empty_allowed && comma_pos == 0) {
                        if (debug) printlog ("____found must_have and \
                            empty_allowed", \
                            l_logfile, verbose+)

                        print (inimg//" "//inname, >> current_outlist)

                    } else if (isallowed && allowed_counter == 0) {
                        if (debug) printlog ("____found isallowed and \
                            allowed_counter=0", \
                            l_logfile, verbose+)

                        allowed_counter += 1
                        print (inimg//" "//inname, >> current_outlist)

                    } else if (comma_pos == 0 && allowed_counter == 0 \
                        && !empty_allowed) {

                        if (debug) printlog ("____found coma_pos=0, \
                            allowed_counter=0 and not empty_allowed", \
                            l_logfile, verbose+)

                        print (inimg//" "//inname, >> not_found_list)
                    }
                } else {
                    if (must_have) {

                        if (debug) printlog ("____not found and must have", \
                            l_logfile, verbose+)

                        printlog ("GAIMCHK: Image "//inname//\
                            " does not contain the "//t_string//" keyword."//\
                            "\n         Skipping this image.", \
                            l_logfile, l_verbose)
                        skip_vardq = yes
                        goto KEY_NEXTIMAGE

                    } else if (not_allowed && empty_allowed \
                        && comma_pos == 0) {
                        if (debug) printlog ("____not found, not allowed and \
                            empty_allowed", \
                            l_logfile, verbose+)

                        print (inimg//" "//inname, >> not_found_list)

                    } else if (comma_pos == 0 && allowed_counter == 0 \
                        && !empty_allowed) {

                        if (debug) printlog ("____not found, comma_pos=0, \
                            allowed_counter = 0 and not empty_allowed", \
                            l_logfile, verbose+)

                        print (inimg//" "//inname, >> not_found_list)
                    }
                }

KEY_NEXTIMAGE:
                if (access(current_outlist) && write_to_outlist) {
                    print (current_outlist, >> l_outlist)
                }

                if (skip_vardq) {
                    goto NEXT_IMAGE
                }
            }

            if (access(not_found_list) && wr_NF_to_outlist) {
                print (not_found_list, >> l_outlist)
            }
            # End of checking inword_list
        } else {
            if (access(not_found_list)) {
                wr_NF_to_outlist = no
            } else {
                wr_NF_to_outlist = yes
            }

            # Output list
            print (inimg//" "//inname, >> not_found_list)

            # Write to Output Output list when appropriate
            if (access(not_found_list) && wr_NF_to_outlist) {
                print (not_found_list, >> l_outlist)
            }
        }

        # Check for VAR and DQ planes - it's a fail if they don't exist and it
        # is prepared
        # For GSAOI it must be prepared!
        keypar (inphu, "GAPREPAR", silent+)
        if (keypar.found && l_fl_vardq_check) {

            # Check for the existance of VAR and DQ planes
            gemextn (inimg, check="exists", process="expand", index="", \
                extname=l_var_ext//","//l_dq_ext, extversion="1-", \
                ikparams="", omit="", replace="", outfile="dev$null", \
                logfile=l_logfile, glogpars="", verbose=no)

            # Check the output status
            if (gemextn.status != 0) {
                printlog ("ERROR- GAIMCHK: GEMEXTN returned a non-zero \
                    status. Exiting.", l_logfile, verbose+)
                goto crash
            } else {

                # Read the number of l_sci_ext - written by gaprepare
                keypar (inphu, "NSCIEXT", silent+)

                # Check that the number returned by gemextn is twice that of
                # the number of l_sci_ext
                if (gemextn.count != (2 * int (keypar.value))) {
                    printlog ("ERROR - GAIMCHK: Image "//inimg//" GAPREPAREd \
                        without "//l_var_ext//" and "//l_dq_ext//" planes "//\
                        "\n                 and fl_vardq=yes. Either supply \
                        input images with "//l_var_ext//" and "//\
                        "\n                 "//l_dq_ext//" planes or set \
                        fl_vardq=no. Exiting.", l_logfile, verbose+)
                    goto crash
                }
            } # End of loop checking gemextn results
        }

NEXT_IMAGE:
    }

    lcount = 0
    if (access(l_outlist)) {
        count (l_outlist) | scan (lcount)
    }

    if (lcount != 0) {
        # Rename the output
        rename (l_outlist, l_orig_outlist, field="all")
    } else {
        printlog ("\nERROR - GAIMCHK: No input files to return. None match \
            the critia requested by parent task.", l_logfile, verbose+)
        goto crash
    }

    goto clean

    #--------------------------------------------------------------------------
crash:

    # Exit with error subroutine
    status = 1
    if (access(l_outlist)) {
        delete ("@"//l_outlist//","//l_outlist, verify-, >& "dev$null")
    }

clean:

    delete (filelist//","//tmplist//","//tmp2list//","//nextlist, \
        verify-, >& "dev$null")
    scanfile = ""

    # Print finish time
    gemdate(zone="local")
    printlog ("\nGAIMCHK -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("\nGAIMCHK -- Exit status: GOOD", l_logfile, l_verbose)
    } else {
        printlog ("\nGAIMCHK -- Exit status: ERROR", l_logfile, l_verbose)
    }

    if (status != 0) {
        printlog ("        -- Please read the logfile \""//l_logfile//\
                "\" for more information.", l_logfile, l_verbose)
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)



end
