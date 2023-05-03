# Copyright(c) 2014-2015 Association of Universities for Research in Astronomy, Inc.

######################################################
procedure gsscatsub(inimages)

######################################################
# Notes:

# Subtraction of instrumental scattered light for GMOS LS and MOS observations
#     Mainly a MEF wrapper for apscatter
#
# Version Feb 21, 2011 BM First version
#         Mar 11, 2011 BM Added interactive flag

# Example call:
#    gsscatsub tcrgsN20020212S122.fits t_order=5 order1=11 sample1="4:192 299:466" order2=7


# TODO Can this be used by other instruments?
#     If, so move to gemtools and rename
#     Also, have user set dispersion axis: then update (internally) the
#     appropriate values / rename parameters
# TODO allow users to supply reference traces
# TODO open up database parameter

######################################################
# Parameters:

string  inimages  {prompt="Images for background subtraction"}
string  outimages {"", prompt="Output image with scattered light subtracted"}
string  prefix    {"b", prompt="Output image prefix"}
int     nfind     {1, prompt="Number of source regions between background"}
int     column    {INDEF, prompt="Column for finding sources (APFIND line par)"}
int     t_order   {3, prompt="Trace order"}
string  order1    {"11", prompt="Cross-dispersion fit order"}
string  sample1   {"*", prompt="Cross-dispersion sample points"}
string  order2    {"7", prompt="Order in dispersion direction"}
int     niterate2 {3, prompt="Clipping iterations along the dispersion"}
string  database  {"database", prompt="Database to store traces"}
bool    fl_inter  {no, prompt="Interactive?"}
bool    fl_display{yes, prompt="Display scattered light image to be removed?"}
bool    fl_vardq  {no, prompt="Propagate variance and data quality extensions?"}
string  logfile   {"", prompt="Logfile"}
bool    verbose   {no, prompt="Verbose?"}
int     status    {0, prompt="Exit status (0=good)"}
struct* scanfile  {"", prompt="Internal use only"}

######################################################
# Task:

begin

    ####
    # Variable declaration

    # Task messaging variables
    string task_name, err_msg, warn_msg, task_msg, mdelim, msg_types[4], tmpstr
    string pad_types[4], err_pad, warn_pad, task_pad, mpad_delim, log_wrap
    string pkg_name
    int    n_msg_types, lwrap_len
    struct l_date_struct

    # Declare variables not related to input parameters
    string tlist_file, tmpoutlist, testfile, output_list, inlist, outlist
    string test_list, tmpinlist, inphu, outfile
    string l_key_prepare, l_key_nsciext, l_key_task
    string bkgimg, tmporder,apsimg, sciextn
    string infile, l_sci_ext, l_var_ext, l_dq_ext, l_mdf_ext, tmplooplist
    string outextn, thisextn, exit_status[2]
    string varextn, outvarextn, dqextn, outdqextn
    string l_databasesave, l_logfilesave
    int    l_dispaxissave, l_ordersave
    int    atpos, nimages, nsci, ii, junk, l_xorder[60], l_yorder[60]
    int    n_order_vals, l_val, l_nsum, l_line, curr_nsci
    int    chk_nimages
    bool   debug, curr_fl_vardq, l_verbosesave, stored

    # Declare local parameter variables
    string l_inimages
    string l_outimages
    string l_prefix
    int    l_nfind
    int    l_column
    int    l_t_order
    string l_order1
    string l_sample1
    string l_order2
    int    l_niterate2
    string l_database
    bool   l_fl_inter
    bool   l_fl_display
    bool   l_fl_vardq
    string l_logfile
    bool   l_verbose
    int    l_status

    # Set local parameter variables
    # fscan strips spaces
    junk = fscan(inimages, l_inimages)
    l_outimages = outimages
    l_prefix = prefix
    l_nfind = nfind
    l_column = column
    l_t_order = t_order
    junk = fscan(order1, l_order1)
    junk = fscan(sample1, l_sample1)
    junk = fscan(order2, l_order2)
    l_niterate2 = niterate2
    junk = fscan(database, l_database)
    l_fl_display = fl_display
    l_fl_inter = fl_inter
    l_fl_vardq = fl_vardq
    l_logfile = logfile
    l_verbose = verbose

    ####
    # Script messages

    # Set task and package names
    task_name = strupr("gsscatsub")
    pkg_name = "gmos"

    # Set up message types
    err_msg = "ERROR - "//task_name//": "
    warn_msg = "WARNING - "//task_name//": "
    task_msg = task_name//": "
    mdelim = "____"//strlwr(task_name)//" -- "

    # Store the message types
    n_msg_types = 4
    msg_types[1] = err_msg
    msg_types[2] = warn_msg
    msg_types[3] = task_msg
    msg_types[4] = mdelim

    # Set up the whitespace padding
    for (i = 1; i <= n_msg_types; i += 1) {
        tmpstr = ""
        for (j = 1; j <= strlen(msg_types[i]); j += 1) {
            tmpstr += " "
        }
        pad_types[i] = tmpstr
    }

    # Sore the whitespace padding for a given message
    err_pad = pad_types[1]
    warn_pad = pad_types[2]
    task_pad = pad_types[3]
    mpad_delim = pad_types[4]

    # Add newline to mdelim
    mdelim = "\n"//mdelim

    # Set up the log file wrapper
    # lwrap_len is the number of "-"
    lwrap_len = 60
    log_wrap = ""
    for (i = 1; i <= lwrap_len; i += 1) {
       log_wrap += "-"
    }

    ####
    # Default values

    # Assume task will fail
    l_status = 1

    debug = no
    if (debug) {
        l_verbose = yes
    }

    l_key_task = substr(task_name, 1, 8)
    l_key_nsciext = "NSCIEXT"
    l_key_prepare = "GPREPARE"
    l_sci_ext = "SCI"
    l_var_ext = "VAR"
    l_dq_ext = "DQ"
    l_mdf_ext = "MDF"

    ####
    # Temporary files
    inlist = mktemp("tmpinlist")
    outlist = mktemp("tmpoutlist")
    tmpinlist = mktemp("tmptmpinlist")
    tmpoutlist = mktemp("tmptmpoutlist")
    tmporder = mktemp("tmporder")
    tmplooplist = mktemp("tmplooplist")
    bkgimg = mktemp("tmpbkgimg")
    apsimg = mktemp("tmpapsimg")

    ####
    # Caching
    cache("imgets", "gimverify", "specred", "aptrace", "apfind")

    ########

    # Start task: User input and log file checking.

    ####
    # Check log file
    if (l_logfile == "" || stridx(" ",l_logfile) > 0) {
        ## Cannot easily use pkg_name to determine the logfile to use
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile) > 0) {
            l_logfile = pkg_name//".log"

            printlog ("\n"//warn_msg//"both "//task_name//".logfile and "//\
                pkg_name//".logfile are empty.\n"//\
                warn_pad//"Using default file "//pkg_name//".log.",
                l_logfile, l_verbose)
        }
    } # End of logfile check

    ####
    # Start log / show time
    date | scan (l_date_struct)
    printlog ("\n"//log_wrap, l_logfile, l_verbose)
    printlog (task_name//" Started -- "//l_date_struct//"\n", \
        l_logfile, l_verbose)

    ####
    # Record all input parameters
    printlog (task_msg//"Input parameters...\n", l_logfile, l_verbose)
    printlog ("    inimages   = "//l_inimages, l_logfile, l_verbose)
    printlog ("    outimages  = "//l_outimages, l_logfile, l_verbose)
    printlog ("    prefix     = "//l_prefix, l_logfile, l_verbose)
    printlog ("    nfind      = "//l_nfind, l_logfile, l_verbose)
    printlog ("    column     = "//l_column, l_logfile, l_verbose)
    printlog ("    t_order    = "//l_t_order, l_logfile, l_verbose)
    printlog ("    order1     = "//l_order1, l_logfile, l_verbose)
    printlog ("    sample1    = "//l_sample1, l_logfile, l_verbose)
    printlog ("    order2     = "//l_order2, l_logfile, l_verbose)
    printlog ("    fl_display = "//l_fl_display, l_logfile, l_verbose)
    printlog ("    fl_inter   = "//l_fl_inter, l_logfile, l_verbose)
    printlog ("    fl_vardq   = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("    logfile    = "//l_logfile, l_logfile, l_verbose)
    printlog ("    verbose    = "//l_verbose, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    ####
    # User input checks

    # Inputs
    test_list = l_inimages
    atpos = stridx("@", test_list)
    if (atpos > 0) {
        tlist_file = substr(test_list, atpos+1, strlen(test_list))
        if (!access(tlist_file)) {
            printlog (err_msg//"Cannot access list file: "//tlist_file, \
                l_logfile, verbose+)
            goto error
        }
    }

    sections (l_inimages, option="fullname", >> tmpinlist)
    scanfile = tmpinlist
    ## Needed?
    nsci = -1
    output_list = inlist
    while(fscan(scanfile, testfile) != EOF) {
        gimverify (testfile)
        if (gimverify.status != 0) {
            printlog (err_msg//"input image \""//testfile//\
                "\" not a valid MEF or does not exist", l_logfile, verbose+)
            goto error
        }
        testfile = gimverify.outname//".fits"

        inphu = testfile//"[0]"

        # Check task not already ran on this file
        keypar (inphu, l_key_task, silent+)
        if (keypar.found) {
            printlog (err_msg//"image "//testfile//" has already been "//\
                "processed with "//task_name, l_logfile, verbose+)
            goto error
        }

        # Check prepared
        keypar (inphu, l_key_prepare, silent+)
        if (!keypar.found) {
            printlog (err_msg//"image "//testfile//" has not been prepared", \
                l_logfile, verbose+)
            goto error
        }

        # Check NSCIEXT exists
        keypar (inphu, l_key_nsciext, silent+)
        if (!keypar.found) {
            printlog (err_msg//"Cannot find "//l_key_nsciext//\
                " keyword in image "//testfile, l_logfile, verbose+)
            goto error
        }

        # Require all images to have the same number of science extensions
        ##TODO possible bug for MOS Data?
        curr_nsci = int(keypar.value)
        if (nsci == -1) {
            nsci = curr_nsci
        } else if (curr_nsci != nsci) {
            printlog (err_msg//"Mismatch of "//l_key_nsciext//\
                " value for image \""//testfile//"\"", l_logfile, verbose+)
            goto error
        }

        ##TODO Test for var and dq inputs if fl_vardq+; reset variable if they
        ##TODO don't exist
        print (testfile, >> output_list)
    } # End of loop checking inputs
    scanfile = ""

    # Check images have actually been supplied and count them for later use
    if (access(output_list)) {
        count(output_list) | scan (nimages)
    } else {
        nimages = 0
    }

    if (nimages == 0) {
        printlog (err_msg//"No input images supplied", l_logfile, verbose+)
        goto error
    }

    # Output checks
    if (l_outimages == "" || stridx (" ",l_outimages) > 0) {
        printlog (task_msg//"Using out prefix \""//l_prefix//"\"", \
            l_logfile, l_verbose)
        # Set out image to prefixed inimages
        l_outimages = l_prefix//"//@"//inlist
    }

    test_list = l_outimages
    atpos = stridx("@", test_list)
    if (atpos > 0) {
        tlist_file = substr(test_list, atpos+1, strlen(test_list))
        if (!access(tlist_file)) {
            printlog (err_msg//"Cannot access list file: "//tlist_file, \
                l_logfile, verbose+)
            goto error
        }
        # Check that number of files supplied match the number of input images
    }

    sections (l_outimages, option="fullname", >> tmpoutlist)
    scanfile = tmpoutlist

    # Check output files do not exist
    output_list = outlist
    while(fscan(scanfile, testfile) != EOF) {
        gimverify (testfile)
        if (gimverify.status != 1) {
            printlog (err_msg//"output image \""//testfile//\
                "\" already exists", l_logfile, verbose+)
            goto error
        }
        testfile = gimverify.outname//".fits"
        print (testfile, >> output_list)
    }
    scanfile = ""

    # Check that the number of outfiles match the input files
    if (access(output_list)) {
        count(output_list) | scan (chk_nimages)
    } else {
        chk_nimages = 0
    }

    if (chk_nimages != nimages) {
        printlog (err_msg//"Number of ouput images supplied ("//\
            chk_nimages//") does not match number of supplied "//\
            "input images ("//nimages//")", l_logfile, verbose+)
        goto error
    }

    ## End of file checks

    # Other input checks and setups
    files(l_order2, sort-, > tmporder)
    count(tmporder) | scan(n_order_vals)
    if (n_order_vals > nsci) {
        printlog (err_msg//"Too many order1 values supplied (>"//\
            max_extensions//")", l_logfile, verbose+)
        goto error
    }

    if (n_order_vals == 1) {
        for (ii = 1; ii <= nsci; ii += 1) {
            l_xorder[ii] = int(l_order2)
        }
    } else if (n_order_vals == nsci) {
        ii = 0
        scanfile = tmporder
        while(fscan(scanfile, l_val) != EOF) {
            ii = ii + 1
            l_xorder[ii] = l_val
        }
        scanfile = ""
    } else {
        print(err_msg//"Number of order2 entries must be one or nsci ("//\
            nsci//")", l_logfile, verbose+)
        goto error
    }
    delete(tmporder, verify-, >& "dev$null")

    files(l_order1, sort-, > tmporder)
    count(tmporder) | scan(n_order_vals)
    if (n_order_vals > nsci) {
        printlog (err_msg//"Too many order1 values supplied (>"//\
            nsci//")", l_logfile, verbose+)
        goto error
    }

    if (n_order_vals == 1) {
        for (ii = 1; ii <= nsci; ii += 1) {
            l_yorder[ii] = int(l_order1)
        }
    } else if (n_order_vals == nsci) {
        ii = 0
        scanfile = tmporder
        while(fscan(scanfile, l_val) != EOF) {
            ii = ii + 1
            l_yorder[ii] = dum
        }
        scanfile = ""
    } else {
        print(err_msg//"Number of order1 entries must be one or nsci ("//\
            nsci//")", l_logfile, verbose+)
        goto error
    }
    delete(tmporder, verify-, >& "dev$null")

    # Checks complete

    ####
    # Hard work part

    joinlines (inlist, outlist, output=tmplooplist, delim=" ",
        missing="MISSING", maxchars=161, shortest=yes, verbose=no)

    l_dispaxissave = specred.dispaxis
    l_verbosesave = specred.verbose
    l_databasesave = specred.database
    l_logfilesave = specred.logfile
    l_ordersave = aptrace.order
    stored = yes

    specred.dispaxis = 1
    specred.verbose = l_verbose
    specred.database = l_database
    specred.logfile = l_logfile
    aptrace.order = l_t_order

    # TODO open these up to users
    l_nsum = 10 # TODO change this depending on type of spectra
    l_line = INDEF

    # Loop over images
    scanfile = tmplooplist
    while (fscan(scanfile, infile, outfile) != EOF) {
        printlog (task_msg//"Working on file: "//infile, l_logfile, l_verbose)

        curr_fl_vardq = l_fl_vardq

        imcopy (infile//"[0]", outfile, verbose=debug)
        # Initiate tinfo.tabletype
        tinfo.tbltype = ""

        # read the MDF table imformation
        tinfo (infile//"["//l_mdf_ext//"]", ttout=no, >& "dev$null")

        # Test the tbltype
        if (tinfo.tbltype == "fits") {
            # Copy out MDF
            nttools.tcopy (intable=infile//"["//l_mdf_ext//"]", \
                outtable=outfile//"["//l_mdf_ext//"]", \
                verbose-)
        }

        # Loop over extensions
        for (ii = 1; ii <= nsci; ii += 1) {
            printlog (task_pad//"Extension version: "//ii, \
                l_logfile, l_verbose)

            sciextn = "["//l_sci_ext//", "//ii
            outextn = sciextn//",append+]"
            sciextn = sciextn//"]"

            thisextn = infile//sciextn
            outextn = outfile//outextn

            varextn = "["//l_var_ext//","//ii
            outvarextn = outfile//varextn//",append+]"
            varextn = infile//varextn//"]"

            dqextn = "["//l_dq_ext//","//ii
            outdqextn = outfile//dqextn//",append+]"
            dqextn = infile//dqextn//"]"

            ## TODO do we need this to happen for each extension
            bkgimg = mktemp("tmpbkgimg")//".fits"
            apsimg = mktemp("tmpapsimg")//".fits"

            # TODO Open more parameters to user
            # Call apscatter
            # Do not subtract in place!
            # First call apfind to create a trace in the database, otherwise
            # apscatter calls that task without its obligatory nfind parameter,
            # causing it to prompt the user even in non-interactive mode.
            apfind(thisextn, nfind=l_nfind, apertures="", references="", \
                interactive=l_fl_inter, find+, recenter-, resize-, edit+, \
                line=l_column, nsum=l_nsum)
                #line=l_column, nsum=l_nsum, apfind.order="increasing")
            apscatter(thisextn, apsimg, apertures="", \
                scatter=bkgimg, interactive=l_fl_inter, find=no, \
                recenter=no, resize=yes, edit=no, trace=yes, fittrace=yes, \
                subtract=yes, smooth=yes, fitscatter=yes, line=l_column, \
                nsum=l_nsum, buffer=2.0, apscat1.order=l_yorder[ii], \
                apscat2.order=l_xorder[ii], \
                apscat1.sample=l_sample1, apscat2.niterate=l_niterate2)

            # TODO Check if it worked - check for outputs?
            #      Delete database file for input img if wasn't there already?

            if (l_fl_display) {
                # TODO Fix this call
                printlog (task_msg//"Displaying scattered light image", \
                    l_logfile, l_verbose)
                display(bkgimg, 1)
            }

            # Apply the correction
            imarith (thisextn, "-", bkgimg, outextn, verbose=debug)

            # TODO write information to header if required

            # TODO If not propagating VARDQ update NEXTEND
            if (curr_fl_vardq) {
                if (imaccess(varextn) && imaccess(dqextn)) {
                    # TODO Use all parameters
                    # The surface fit is approximately noiseless since it
                    # should have a low order and be highly overconstrained,
                    # so propagate the input variance unmodified; the change in
                    # apparent S/N is reflected in the signal level instead.
                    imcopy (varextn, outvarextn, verbose-)
                    imcopy (dqextn, outdqextn, verbose-)
                } else if (ii == 1) {
                    printlog (warn_msg//l_var_ext//" and / or "//l_dq_ext//\
                        "not present in "//infile, \
                        l_logfile, verbose+)
                    curr_fl_vardq = no
                }
            }
            imdelete (bkgimg//","//apsimg, verify-, >& "dev$null")

        } # End of loop over extensions

        gemdate()
        gemhedit (outfile//"[0]", l_key_task, gemdate.outdate,
            "UT Time stamp for "//task_name, delete-)
        gemhedit (outfile//"[0]", "GEM-TLM", gemdate.outdate,
            "Last modification with GEMINI", delete-)

    } # End of looping over images
    delete (tmplooplist, verify=no, >& "dev$null")

    # If we made it here all is good
    l_status = 0
    goto clean

error:
;

clean:

    if (stored) {
        specred.dispaxis = l_dispaxissave
        specred.verbose = l_verbosesave
        specred.database = l_databasesave
        specred.logfile = l_logfilesave
        aptrace.order = l_ordersave
    }

    # Clean up a bit
    delete (tmpinlist//","//inlist//","//tmpoutlist//","//\
        outlist//","//tmporder//","//tmplooplist//","//bkgimg//","//apsimg, \
        verify-, >& "dev$null")

    scanfile = ""
    status = l_status

    # Finish time
    date | scan (l_date_struct)
    printlog ("\n"//task_name//" Finished -- "//l_date_struct//"\n", \
        l_logfile, l_verbose)

    exit_status[1] = "GOOD"
    exit_status[2] = "ERROR"

    printlog (task_name//" exit status:  "//exit_status[l_status + 1]//\
        ".", l_logfile, l_verbose)

    printlog (log_wrap//"\n", l_logfile, l_verbose)

end
