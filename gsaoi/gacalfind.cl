# Copyright(c) 2010-2013 Association of Universities for Research in Astronomy, Inc.

procedure gacalfind (inimages)

# The table created is stored in the PWD. However, it stores the absolute path
# to the calpath.

# Takes GSAOI processed flats and darks (4 extensions) and generates
# a table with information to match calibrations when reducing science
# data.

# It returns the appropriate calibration when requested.

# This is an ancillary task called by gareduce, but can be run on its own
# The output is a table with image name, caltype (dark or flat), and the
# bits of information from the metaconf keyword required to match with the
# science (for flats the filters, for darks the detector information).
#
# The task has an overwrite flag. If set to yes, will update the table
# (actually delete it and create a new one). If set to no, it checks if the
# table exists and run anyway if it does not. This way, the flag when
# called from gareduce is set to no to speed up execution.
#
# Can be slow if many images present. It select the candidates
# by INSTRUME=GSAOI, then goes through the list selecting the keywords
# GAFLAT/GADARK, OBSTYPE and METACONF. If OBSTYPE = object or flat, then is
# a flat, if obstype = dark, then is a dark.
# Parses metaconf correctly depending on obstype, and writes the table

# Version 2010 Oct 19 V0.0 CW - created
#         2010 Nov 17 V0.1 CW - mostly complete, no error handling
#
# More information on updates in CVS logs

######################################################

# input images and directory, output table
char    inimages   {prompt="Input GSAOI calibration images or list"}
char    calpath    {"./", prompt="Path for input calibration images"}
char    caltable   {"gsaoical.fits", prompt="Name of output table with calibration information"}
bool    fl_calrun  {no, prompt="When looking for calibrations, first re-create existing table"}
bool    fl_find    {no,prompt="Find the calibration specified by caltype for sciimg"}
char    caltype    {"",prompt="Type of calibration to find (""|DARK|FLAT<{DOME|TWLT|GCAL}>)"}
char    sciimg     {"",prompt="Image to find appropriate calibration for"}
char    out_calimg {"",prompt="Output calibration image"}
int     maxtime    {INDEF,prompt="Maximum allowed time for association with input image"}
bool    ignore_nlc {no,prompt="Ignore whether the calibration has been NLC?"}
char    sci_ext    {"SCI",prompt="Name of science extensions"}
char    var_ext    {"VAR",prompt="Name of variance extensions"}
char    dq_ext     {"DQ",prompt="Name of data quality extensions"}
char    logfile    {"",prompt="Logfile"}
bool    verbose    {yes,prompt="Verbose?"}
int     status     {0,prompt="Exit status (0=good)"}
struct  *scanfile  {"",prompt="Internal use only"}

######################################################

begin

    # variable declaration

    char l_inimages, l_calpath, l_caltable, l_logfile
    char filelist, tmpdata, tmpcdfile, l_temp, img, inimg
    char l_dir, sciphu, obskey, intype, ttype
    char metakey, metabit1, metabit2, incaltype, objct, flattype, img_root
    char l_datename, l_timename, uttime, datevalue, path_dir, path_machine
    char rel_path, l_pwd, out_dir, new_path, test_char, cal_table, tabphu
    char tabext, l_sciimg, l_caltype, envvar, outgaimchk_list, tmpfile
    char l_sci_ext, l_var_ext, l_dq_ext, rphend, pathtest, cal_requested
    char selexpr, proc_test, tmpproc_test, testmeta1, tmpcal_caltab
    char calfile, metaconf, inphu, l_calimg, l_tab_ext, rel_ptest[2]
    char type_requested

    string full_path

    int i, ii, l_status, l_test, ex_mark_pos, path_length, name_len_max
    int junk, obs_time, num_dir_strings, last_char, l_maxtime, row_id, nrows
    int num_counter, nsci, p_location, tvalmin, npathcols, max_meta1_len
    int max_meta2_len, len_rel_ptest, brk_pos1, brk_pos2, nlines
    int im_type_len

    bool l_fl_calrun, l_verbose, l_fl_find, l_create_table, chk_rel_path
    bool l_ignore_nlc, processing, redo_selection, dofind

    # Allowed input image types (the ones allowed to find calibrations for
    # The last one is used for selecting the requested flat type
    int n_allowed=6
    char allowed_types[6]="OBJ","SKY","DOME","TWLT","GCAL","DOMEON-OFF"

    # temporary files

    filelist  = mktemp("tmpfile")
    tmpdata = mktemp("tmpdata")
    tmpcdfile = mktemp("tmpcd")
    outgaimchk_list = mktemp("tmpoutgaimchk_list")
    tmpcal_caltab = mktemp("tmpcal_caltab")//".fits"

    # local variables

    l_inimages = inimages
    l_calpath = calpath
    l_caltable = caltable
    l_logfile = logfile
    l_status = status
    l_fl_find = fl_find
    l_caltype = caltype
    l_sciimg = sciimg
    l_fl_calrun = fl_calrun
    l_maxtime = maxtime
    l_ignore_nlc = ignore_nlc
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_verbose = verbose

    # Default values
    rel_path = ""
    status = 0
    out_calimg = ""
    type_requested = ""

    # Maybe in teh future append to an exiting calibratiuon file use this
    # parameter? and update check later on. - MS
    l_tab_ext = "CALINFO"

    # test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = gsaoi.logfile
       if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
          l_logfile = "gsaoi.log"
          printlog ("WARNING - GACALFIND: Both gareduce.logfile and \
                         gsaoi.logfile fields are empty", l_logfile, l_verbose)
          printlog ("                    Using default file gsaoi.log", \
              l_logfile, l_verbose)
       }
    }

    # Print start time
    printlog ("", l_logfile, l_verbose)
    gemdate (zone="local")
    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)
    printlog ("GACALFIND -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    printlog ("GACALFIND: Input parameters...", l_logfile, l_verbose)
    printlog ("    inimages    = "//l_inimages, l_logfile, l_verbose)
    printlog ("    calpath     = "//l_calpath, l_logfile, l_verbose)
    printlog ("    caltable    = "//l_caltable, l_logfile, l_verbose)
    printlog ("    fl_calrun   = "//l_fl_calrun, l_logfile, l_verbose)
    printlog ("    caltype     = "//l_caltype, l_logfile, l_verbose)
    printlog ("    sciimg      = "//l_sciimg, l_logfile, l_verbose)
    if (isindef(l_maxtime)) {
        printlog ("    maxtime     = INDEF", l_logfile, l_verbose)
    } else {
        printlog ("    maxtime     = "//l_maxtime, l_logfile, l_verbose)
    }
    printlog ("    ignore_nlc  = "//l_ignore_nlc, l_logfile, l_verbose)
    printlog ("    sci_ext     = "//l_sci_ext, l_logfile, l_verbose)
    printlog ("    var_ext     = "//l_var_ext, l_logfile, l_verbose)
    printlog ("    dq_ext      = "//l_dq_ext, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    # Check the tab_ext paramemter
    if (l_tab_ext == "" || stridx(" ",l_tab_ext) > 0) {
        printlog ("ERROR - GACALFIND: tab_ext parameter has not be set.", \
            l_logfile, verbose+)
        goto crash
    }

    # Check that the calpath has a trailing slash and is a valid entry
    if (l_calpath != "" && stridx(" ",l_calpath) == 0) {
        rphend = substr(l_calpath,strlen(l_calpath),strlen(l_calpath))

        if (rphend == "$") {
            envvar = substr(l_calpath,1,strlen(l_calpath)-1)
            if (defvar (envvar)) {
                show (substr(l_calpath,1,strlen(l_calpath)-1)) | \
                    scan (pathtest)
                rphend = substr(pathtest,strlen(pathtest),strlen(pathtest))

                if (rphend != "/") {
                    l_calpath = l_calpath//"/"
                }
            }
        }

        if (!access(l_calpath)) {
            printlog ("ERROR - GACALFIND: Cannot access calpath: "//\
                l_calpath, l_logfile, verbose+)
            goto crash
        }
    } else {
        printlog ("ERROR - GACALFIND: calpath is set incorrectly", \
            l_logfile, verbose+)
        goto crash
    }

    printlog ("GACALFIND: calpath after checks is \""//l_calpath//"\"\n", \
        l_logfile, l_verbose)

    if (l_caltable == "" || stridx(" ",l_caltable) > 0) {
        printlog ("ERROR - GACALFIND: caltable parameter is not set.", \
            l_logfile,, verbose+)
        goto crash
    } else {
        fparse (l_caltable)
        if (fparse.extension != ".fits") {
            printlog ("WARNING - GACALFIND: Setting output table extension to \
                \".fits\"", l_logfile, verbose+)
            l_caltable = fparse.directory//fparse.root//".fits"
        }
    }

    # Check for the exitance of the table - depending on overwrite flag act
    # accordingly
    l_create_table = yes
    if (!l_fl_calrun) {
        if (access (l_caltable)) {

            printlog ("GACALFIND: Found calibration table "//l_caltable,
                l_logfile, l_verbose)
            l_create_table = no
        } else {
            printlog ("GACALFIND: Calibration table "//l_caltable//\
                " does not exist. Creating it.", l_logfile, verbose)
        }
    } else {
        if (access (l_caltable)) {

            printlog ("GACALFIND: Calibration table "//l_caltable//\
                " already exists. Deleting and re-creating.", \
                l_logfile, verbose+)
            tdelete (l_caltable, verify-, >& "dev$null")
        } else {
            printlog ("GACALFIND: Calibration table "//l_caltable//\
                " does not exist. Creating it.", l_logfile, l_verbose)
        }
    }

    # Perform checks on input images when requested to find calibrations
    cal_requested = ""
    if (l_fl_find) {
        if (l_sciimg == "" || stridx(" ",l_sciimg) > 0) {
            printlog ("ERROR - GACALFIND: sciimg is not set properly", \
                l_logfile, verbose+)
            goto crash
        } else if (!imaccess(l_sciimg)) {
            printlog ("ERROR - GACALFIND: Cannot access \""//l_sciimg//"\".",\
                l_logfile, verbose+)
            goto crash
        } else {
            keypar (l_sciimg//"[0]", "GAPREPAR", silent+)
            if (!keypar.found) {
                if (!l_create_table) {
                    printlog ("ERROR - GACALFIND: \""//l_sciimg//"\" is not \
                        prepared. Exiting.", l_logfile, verbose+)
                    goto crash
                } else {
                    printlog ("WARNING - GACALFIND: \""//l_sciimg//"\" is not \
                        prepared. Setting fl_find=no", l_logfile, verbose+)
                    l_fl_find = no
                }
            } else {
                keypar (l_sciimg//"[0]", "METACONF", silent+)
                if (!keypar.found || keypar.value == "UNKNOWN") {

                    if (!l_create_table) {
                        printlog ("ERROR - GACALFIND: \""//l_sciimg//"\" "//\
                            "METACONF keyword is UNKNOWN or does not exist.", \
                            l_logfile, verbose+)
                        goto crash
                    } else {
                        printlog ("WARNING - GACALFIND: \""//l_sciimg//"\ "//\
                            "METACONF keyword is UNKNOWN or does not exist."//\
                            " Setting fl_find=no", l_logfile, verbose+)
                        l_fl_find = no
                    }
                } else {
                    metaconf = keypar.value
                    keypar (l_sciimg//"[0]", "NSCIEXT", silent+)
                    nsci = int(keypar.value)
                }
            }
        }
    } # End of first l_fl_find loop

    # Check the other setting for finding calibrations. Done again here as
    # l_fl_find can be switched off
    if (l_fl_find) {
        if (l_caltype == "" || stridx(" ",l_caltype) > 0) {
            printlog ("ERROR - GACALFIND: caltype is not set properly", \
                l_logfile, verbose+)
            goto crash
        } else {
            cal_requested = strupr (l_caltype)

            # Test for requested FLAT type
            brk_pos1 = stridx("{",l_caltype)
            brk_pos2 = stridx("}",l_caltype)
            if (brk_pos1 > 0 || brk_pos2  > 0) {
                if (brk_pos2 > brk_pos1 && brk_pos2 == strlen(cal_requested)) {
                    type_requested = substr(cal_requested,brk_pos1,\
                                         strlen(cal_requested))
                    cal_requested = substr(cal_requested,1,brk_pos1-1)
                } else {
                    printlog ("ERROR - GACALFIND: FLAT kind is not set \
                        properly")
                    goto crash
                }
            }

            if (cal_requested != "DARK" && cal_requested != "FLAT") {
                printlog ("ERROR - GACALFIND: caltype can only be FLAT or \
                    DARK. Please Try again.", l_logfile, verbose+)
                goto crash
            } else if (cal_requested == "FLAT" && type_requested != "") {
                for (i = 1; i <= n_allowed; i += 1) {
                    ttype = allowed_types[i]
                    if (strstr(ttype,type_requested) > 0) {
                        type_requested = ttype
                        break
                    } else if (i == n_allowed) {
                         printlog ("ERROR - GACALFIND: Unable to determine \
                             requested FLAT type type", l_logfile, verbose+)
                         goto crash
                    }
                }
            }
        }
    } # End of second l_fl_find loop

    # Check if anything needs to be done skip as required
    if (!l_create_table) {
        if (!l_fl_find) {
            printlog ("GACALFIND: There is nothing for GACALFIND \
                to do.", l_logfile, verbose+)
            goto clean
        } else {
            goto RETURN_CAL
        }
    }

    # Check input images...
    printlog ("GACALFIND: Calling GAIMCHK to check input files...", \
        l_logfile, l_verbose)

    # Call gaimchk to perform input checks
    gaimchk (inimages=l_inimages, rawpath=l_calpath, rootname="", \
        obstype_allowed="", object_allowed="", \
        key_allowed="GADARK,GAFLAT", key_forbidden="GAMOSAIC", \
        key_exists="GAPREPAR", \
        fl_prep_check=no, gaprep_pref="g", fl_redux_check=no, \
        garedux_pref="r", fl_fail=no, fl_out_check=no, \
        fl_vardq_check=no, sci_ext=l_sci_ext, var_ext=l_var_ext, \
        dq_ext=l_dq_ext, outlist=outgaimchk_list, logfile=l_logfile, \
        verbose=no)

    if (gaimchk.status != 0) {
        printlog ("ERROR - GACALFIND: GAIMCHK returned a non-zero status. \
            Exiting.", l_logfile, verbose+)
        goto crash
    } else if (!access(outgaimchk_list)) {
        printlog ("ERROR - GACALFIND: Cannot access output list from \
            GAIMCHK", l_logfile, verbose+)
        goto crash
    }
    printlog ("GACALFIND: Returned from GAIMCHK.", l_logfile, l_verbose)

    # Files that have have already been processed by GADARK
    tmpfile = ""
    match ("tmpGADARK", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> filelist)
    }

    # Files that have have already been processed by GAFLAT
    tmpfile = ""
    match ("tmpGAFLAT", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> filelist)
    }

    nlines = 0
    if (access(filelist)) {
        count (filelist) | scan (nlines)
    }
    if (nlines < 1) {
        printlog ("ERROR - GACALFIND: No appropriate input files found", \
            l_logfile, verbose+)
        goto crash
    }

    # Break the path up into more managable strings
    # Computer name
    pathnames (l_calpath) | scan (full_path)
    ex_mark_pos = stridx("!",full_path)
    path_machine = substr(full_path,1,ex_mark_pos)
    path_length = strlen(full_path) - ex_mark_pos
    path_dir = substr(full_path,ex_mark_pos+1,strlen(full_path))
    out_dir = ""
    num_dir_strings = 1

    # tabpar cannot handle > 64 charaters!
    if (path_length > 64) {
        l_test = 1
        new_path = substr(path_dir,1,64)

        while (l_test > 0) {
           # Find the last / before 64 characters
           last_char = strlen(new_path)

           test_char = substr(new_path,last_char,last_char)
           if (test_char == "/") {
               l_test = 0
               out_dir = out_dir//new_path

               path_dir = substr(path_dir,\
                   strstr(new_path,path_dir)+strlen(new_path), \
                   strlen(path_dir))

               new_path = path_dir

               if (new_path != "") {
                   out_dir = out_dir//"    "
                   num_dir_strings += 1
                   if (strlen(new_path) > 64) {
                       l_test = 1
                       new_path = substr(new_path,1,64)
                   } else {
                       out_dir = out_dir//new_path
                       l_test = 0
                   }
               }
           } else {
               new_path = substr(new_path,1,last_char-1)
           }
       }
    } else {
        out_dir = path_dir
    }

    # Here the main loop starts.

    scanfile = ""
    scanfile = filelist
    name_len_max = 0
    max_meta1_len = 0
    max_meta2_len = 0

    while (fscan(scanfile, img) != EOF) {
        fparse (img, verbose-)
        inimg = fparse.root//".fits"
        img_root = inimg
        l_dir = fparse.directory
        inimg = l_dir//inimg
        sciphu = inimg//"[0]"
        rel_path = l_dir

        len_rel_ptest = 2
        rel_ptest[1] = l_calpath
        rel_ptest[2] = "./"

        # Remove all "./" from beginning and check for l_calpath
        chk_rel_path = yes
        while (chk_rel_path) {

            # Check for l_calpath and "./" at the beginning of the strings
            for (ii = 1; ii <= len_rel_ptest; ii += 1) {
                if (strstr(rel_ptest[ii],rel_path) == 1) {
                    if (strlen(rel_path) == strlen(rel_ptest[ii])) {
                        rel_path = ""
                    } else {
                        rel_path = substr(rel_path,strlen(rel_ptest[ii])+1,\
                            strlen(rel_path))
                    }
                } else {
                    chk_rel_path = no
                }
            }
        }

        if (strlen(img_root) > name_len_max) {
            name_len_max = strlen(img_root)
        }

        printlog ("GACALFIND: processing file "//inimg//".", \
            l_logfile, verbose-)

        # get the metaconf and parses accordingly
        keypar (sciphu,"METACONF", silent+)
        if (!keypar.found || keypar.value == "UNKNOWN") {
            printlog ("ERROR - GACALFIND: METACONF keyword is either missing \
                or is UNKOWN in "//inname//". Exiting.", l_logfile, verbose+)
            goto crash
        }
        metakey = keypar.value

        l_timename = "UT"
        # The observation time of the flat
        keypar (sciphu, l_timename, silent+)
        if (!keypar.found) {
           printlog ("ERROR - GACALFIND: "//l_timename//" keyword not \
               found in "//off_flat//". Exiting", \
               l_logfile, verbose+)
           goto crash
        }
        uttime = keypar.value

        l_datename = "DATE-OBS"
        keypar (sciphu, l_datename, silent+)
        if (!keypar.found) {
           printlog ("ERROR - GACALFIND: "//l_datename//" keyword not \
               found in "//off_flat//". Exiting", \
               l_logfile, verbose+)
           goto crash
        }
        datevalue = keypar.value

        cnvtsec (datevalue, uttime) | scan (obs_time)

        # now if obstype = dark the metabit I need for GAREDUCE is the detector
        # configuration
        # which is the METACONF in the prepared Dark frame minus the ROI
        # otherwise, GAREDUCE needs the filter configuration for selecting the
        # flat
        # and the METACONF in the flat header is
        # FILTER1_FILTER2+ITIME_LNRS_COADDS+GCALSHUT+ROI
        # so the metabit I need is from the start to the first +
        # then from the last + to the end

        junk = stridx("+",metakey)
        metabit1 = substr(metakey,1,junk-1)
        metabit2 = substr(metakey,junk+1,strlen(metakey))

        # Check for DARK and FLAT in string - this is a very secondary check
        if ((strstr("DARK",metabit2) == 0) && (strstr("FLAT",metabit2) == 0)) {
            printlog ("WARNING - GACALFIND: The METACONF keyword value for "//\
                inimg//"\n                     does not contain a valid \
                value (\"DARK\" or \"FLAT\")."//\
                "\n                     Skipping image", \
                l_logfile, verbose+)
            goto NEXTIMG
        }

        junk = stridx("+",metabit2)
        incaltype = substr(metabit2,1,junk-1)
        flattype = "NULL"
        if (strstr("_FLAT",incaltype) > 0) {
            flattype = substr(incaltype,1,strstr("_FLAT",incaltype)-1)
            incaltype = "FLAT"
        }
        metabit2 = substr(metabit2,junk+1,strlen(metabit2))

        if (strlen(metabit1) > max_meta1_len) {
            max_meta1_len = strlen(metabit1)
        }
        if (strlen(metabit2) > max_meta2_len) {
            max_meta2_len = strlen(metabit2)
        }

        # write to temporary file IMAGENAME, OBSTYPE, METABIT
        print (img_root//"    "//incaltype//"    "//flattype//\
            "    "//metabit1//"    "//metabit2//"    "//\
            obs_time//"    "//path_machine//"   "//out_dir//"    "//\
            rel_path, >> tmpdata)

# next image
NEXTIMG:

    }

    # Setup configuration file for creating table
    print ("IMAGE_NAME   ch*"//name_len_max, >> tmpcdfile)
    print ("CALTYPE      ch*12", >> tmpcdfile)
    print ("FLATTYPE     ch*12", >> tmpcdfile)
    print ("METABIT1     ch*"//max_meta1_len, >> tmpcdfile)
    print ("METABIT2     ch*"//max_meta2_len, >> tmpcdfile)
    print ("OBS_TIME    int*12", >> tmpcdfile)

    # PATH1 is the OS path to the path supplied by the user
    # PATH2 is the path from that directory to the path containing the
    # calibrations
    print ("MACHINE     ch*"//ex_mark_pos, >> tmpcdfile)
    for (i = 1; i <= num_dir_strings; i += 1) {
        print ("PATH"//str(i)//"       ch*64", >> tmpcdfile)
    }
    print ("REL_PATH       ch*64", >> tmpcdfile)

    pwd | scan (l_pwd)

    # Now create the table
    tcreate (l_caltable, tmpcdfile, tmpdata, uparfile="", nskip=0, \
        nlines=1, nrows=0, tbltype="default")

    # Update the table header with some information
    cal_table = l_caltable
    tabphu = cal_table//"[0]"
    tabext = cal_table//"[1]"

    thedit (tabext, "EXTNAME", l_tab_ext, delete-, show=no)
    thedit (tabext, "EXTVER", 1, delete-, show=no)

    thedit (tabext, "i_title", "Cal info", delete-, show=no)
    thedit (tabext, "NPATHCOL", num_dir_strings, delete-, show=no)

    gemhedit (tabphu, "i_title", "Cal Table", \
        "GSAOI Cal Table from gsaoi.gacalfind", delete-)

    gemdate ()
    gemhedit (tabphu, "GACALFIN", gemdate.outdate, \
        "UT Time stamp for GACALFIND")

    gemhedit (tabphu, "GEM-TLM", gemdate.outdate, \
        "UT Last modification with GEMINI")

RETURN_CAL:

    # Skip the next bit if not finding the calibration to return.
    if (!l_fl_find) {
        goto clean
    }

    # Find the requested calibration for the input file
    l_caltable = l_caltable
#    l_caltable = l_calpath//l_caltable

    # Want to find an appropriate calibration
    # The default will be if maxtime is indef then an exact match of the
    # metaconfiguration keyword takes presidence.

    ##M Can thoretically switch off NLC too
    # If maxtime is not indef then if cannot find a complete match in the time
    # frame find the closest match in that time range!

    # Determine the UT in seconds of the file
    l_timename = "UT"
    inphu = l_sciimg//"[0]"

    # The observation time of the flat
    keypar (inphu, l_timename, silent+)
    if (!keypar.found) {
       printlog ("ERROR - GAREDUCE: "//l_timename//" keyword not \
           found in "//off_flat//". Exiting", \
           l_logfile, verbose+)
       goto crash
    }
    uttime = keypar.value

    l_datename = "DATE-OBS"
    keypar (inphu, l_datename, silent+)
    if (!keypar.found) {
       printlog ("ERROR - GAREDUCE: "//l_datename//" keyword not \
           found in "//off_flat//". Exiting", \
           l_logfile, verbose+)
       goto crash
    }
    datevalue = keypar.value

    cnvtsec (datevalue, uttime) | scan (obs_time)

    printlog ("GACALFIND: Metaconfiguration of input image is: "//metaconf, \
        l_logfile, l_verbose)

    # Determine which allowed intype it is
    for (i = 1; i <= n_allowed; i += 1) {
        ttype = allowed_types[i]
        if (strstr(ttype,metaconf) > 0) {
            intype = ttype
            break
        } else if (i == n_allowed) {
             printlog ("ERROR - GACALFIND: Unable to determine input image \
                 type", l_logfile, verbose+)
             goto crash
        }
    }

    if (cal_requested == "DARK") {
        testmeta1 = substr(metaconf,1,stridx("+",metaconf)-1)
    } else {
        testmeta1 = substr(metaconf,stridx("+",metaconf)+1,\
            strstr(intype,metaconf)-2)
    }

    im_type_len = strlen(intype)
    if (intype == "DOME" || intype == "GCAL") {
        im_type_len += 5
    }
    proc_test = substr(metaconf,strstr(intype,metaconf)+im_type_len,\
        strlen(metaconf))

    # In case the file has been processed in someway.
    if (strstr("_",proc_test) == 1) {
        proc_test = substr(proc_test,strstr("+",proc_test)+1,\
            strlen(proc_test))
    } else {
        proc_test = substr(proc_test,2,strlen(proc_test))
    }

    if (l_ignore_nlc) {
        proc_test = substr(proc_test,1,strstr("NLC",proc_test)-2)
    }

    tmpproc_test = proc_test

    ##M If user requests a certain type of flat add the select
    ##M expression

    num_counter = 0
    row_id = 1
    processing = yes
    while (processing) {
        calfile = ""
        redo_selection = no

        num_counter += 1

        # Read the table for the rquired information
        selexpr = "CALTYPE == '"//cal_requested//"' && "//\
            "METABIT1 == '"//testmeta1//"'"

        if (type_requested != "") {
            selexpr = selexpr//" && FLATTYPE == '"//type_requested//"'"
        }

        if (tmpproc_test != "") {
            selexpr = selexpr//" && "//"METABIT2 ?= '"//tmpproc_test//"'"
        }

        # Only do this twice - exact match then ROI match
        if (num_counter <= 2) {
            tselect (l_caltable, tmpcal_caltab, expr=selexpr)
            tinfo (tmpcal_caltab, ttout=no)
            nrows = tinfo.nrows
        }

        if (nrows != 0) {

            tcalc (tmpcal_caltab, "OBS_DIFF", \
                "abs(OBS_TIME - "//obs_time//")", datatype="int")

            tsort (tmpcal_caltab, columns="OBS_DIFF", ascend=yes, casesens=yes)

            tabpar (tmpcal_caltab, "OBS_DIFF", row_id)
            if (tabpar.undef) {
                printlog ("ERROR - Cannot read differnce in observation \
                    times", l_logfile, verbose+)
                goto crash
            } else {
                tvalmin = int(tabpar.value)
            }

            # Check if it's within the requested time if needed
            # Select output

            if (isindef(l_maxtime)) {
                dofind = yes
            } else if (tvalmin < l_maxtime) {
                dofind = yes
            } else {
                dofind = no
            }

            if (dofind) {

                tabpar (tmpcal_caltab, "IMAGE_NAME", row_id)

                if (tabpar.undef) {
                    printlog ("WARNING - GAREDUCE: Cannot access IMAGE_NAME \
                       in "//tmpcal_caltab//".", l_logfile, verbose+)

                } else {
                    l_calimg = tabpar.value

                    # Need to figure out how many paths there are!
                    calfile = ""
                    keypar (tmpcal_caltab, "NPATHCOL", silent+)
                    npathcols = int(keypar.value)
                    for (i = 1; i <= npathcols; i += 1) {
                        tabpar (tmpcal_caltab, "PATH"//str(i), row_id)
                        calfile = calfile//tabpar.value
                    }

                    tabpar (tmpcal_caltab, "REL_PATH", row_id)
                    # If it conatins " " then there is no relative path
                    if (stridx(" ",tabpar.value) == 0) {
                        calfile = calpath//tabpar.value
                    }

                    pwd | scan (l_pwd)
                    if (strstr(l_pwd,calfile) == 1) {
                        calfile = "./"//\
                            substr(calfile,\
                            strstr(l_pwd,calfile)+strlen(l_pwd)+1,\
                            strlen(calfile))
                    }

                    calfile = calfile//l_calimg

                    for (ii = 1; ii <= nsci; ii += 1) {
                        # Check the dimensions
                        gadimschk (l_sciimg//"["//l_sci_ext//","//ii//"]", \
                            section = "", \
                            chkimage=calfile//"["//l_sci_ext//","//ii//"]", \
                            key_check="CCDSEC", logfile=l_logfile, verbose=no)
                        if (gadimschk.status != 0) {
                            redo_selection = yes
                            break
                        } else if (ii == nsci) {
                            redo_selection = no
                        }
                    }
                }
            } else {
                redo_selection = yes
            }
        } else if (num_counter == 1) {
            redo_selection = yes
        } else {
            printlog ("WARNING - GACALFIND: Cannot find an appropriate "//\
                cal_requested//" for "//l_sciimg, l_logfile, verbose+)
            goto cal_clean
        }

        if (redo_selection) {
            processing = yes
            if (num_counter == 1) {
                tdelete (tmpcal_caltab, verify-, >& "dev$null")

                p_location = stridx("+",proc_test)

                if (p_location > 0) {
                    # Strip the
                    tmpproc_test = substr(proc_test,p_location+1,\
                        strlen(proc_test))
                } else {
                    tmpproc_test = ""
                }
            } else if (row_id == nrows) {
                printlog ("WARNING - GACALFIND: Cannot find an appropriate "//\
                    cal_requested//" for "//l_sciimg, l_logfile, verbose+)
                goto cal_clean
            } else {
                row_id += 1
            }
        } else if (calfile != "") {
            processing = no
            out_calimg = calfile
        }

    } # End of processing loop

    printlog ("GACALFIND: Output calibration image is: "//out_calimg, \
        l_logfile, l_verbose)

    goto clean

#--------------------------------------------------------------------------
crash:
    # Exit with error subroutine
    status = 1

cal_clean:

    out_calimg = ""

clean:

    delete (filelist, verify-, >& "dev$null")
    delete (tmpcdfile, verify-, >& "dev$null")
    delete (tmpdata, verify-, >& "dev$null")

    if (access(tmpcal_caltab)) {
        tdelete (tmpcal_caltab, verify-, >& "dev$null")
    }

    if (access(outgaimchk_list)) {
        delete ("@"//outgaimchk_list//","//outgaimchk_list, \
            verify-, >& "dev$null")
    }

    scanfile = ""

    gemdate (zone="local")
    printlog ("\nGACALFIND -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("\nGACALFIND -- Exit status: GOOD", l_logfile, l_verbose)
    } else {
        printlog ("\nGACALFIND -- Exit status: ERROR", l_logfile, l_verbose)
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)

end


