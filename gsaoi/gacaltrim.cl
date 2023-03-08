# Copyright(c) 2013 Association of Universities for Research in Astronomy, Inc.

# Write a script to trim calibration data to the science science size, if
# required, updating any rquired information as desired. Output is a file with
# a name. Not sure how wrapper script will handle the interface but I think
# this script will save space in gareduce.

# If no trimming is required the out name parameter will be set to the input
# calimg

######################################################

procedure gacaltrim (sciimg, calimg)

# Input / output parameters
char    sciimg      {prompt="Input GSAOI image"}
char    calimg      {prompt="Input GSAOI calibration image to check"}
char    outimg      {"",prompt="Outname of trimmed calibration file"}
char    key_check   {"CCDSEC",enum="CCDSEC|DETSEC",prompt="Keyword to check dimensions against"}

# Other stuff
char    sci_ext     {"SCI",prompt="Name of science extensions"}
char    var_ext     {"VAR",prompt="Name of variance extensions"}
char    dq_ext      {"DQ",prompt="Name of data quality extensions"}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="For internal use only"}

######################################################

begin

    ########
    ####
    # Local variables

    # Local versions of input parameters
    char l_sciimg, l_calimg, l_outname, l_key_check, l_logfile
    char l_sci_ext, l_var_ext, l_dq_ext
    bool l_verbose
    int  l_status

    # Logging / debugging
    char task_name, err_msg, warn_msg, task_msg, mdelim, msg_types[4], tmpstr
    char err_pad, warn_pad, mpad_delim, log_wrap, exit_status[2], pad_types[4]
    char task_pad
    int  num_msg_types
    int  lwrap_len
    bool debug

    # General
    char ext_names[3], sec_to_copy, input, output, tmp_tocopy, sciphu, calphu
    char tmpout_calfile, secvalue, inextn, calextn
    int  nsciext, i, j, next_names
    bool doscichk, docalchk, have_to_copy, update_nextend
    struct l_date_struct

    ####
    # Assign inputs to local variables
    l_sciimg = sciimg
    l_calimg = calimg
    l_key_check = key_check
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_logfile = logfile
    l_verbose = verbose

    ####
    # Set task and package names
    task_name = strupr("gacaltrim")

    ####
    # Default values
    l_status = 1 # Assume it will fail reset just before end if not
    debug = no

    # Reset verbosity if debug set to true
    if (debug)
        l_verbose = yes

    ####
    # Temporay files
    tmpout_calfile = mktemp("tmpcal")//".fits"
    tmp_tocopy = mktemp("tmpcopylist")

    ####
    # Set up message types
    err_msg = "ERROR - "//task_name//": "
    warn_msg = "WARNING - "//task_name//": "
    task_msg = task_name//": "
    mdelim = "____"//strlwr(task_name)//" -- "

    # Store the message types
    num_msg_types = 4
    msg_types[1] = err_msg
    msg_types[2] = warn_msg
    msg_types[3] = task_msg
    msg_types[4] = mdelim

    # Set up the whitespace padding
    for (i = 1; i <= num_msg_types; i += 1) {
        tmpstr = ""
        for (j = 1; j <= strlen(msg_types[i]); j += 1) {
            tmpstr += " "
        }
        pad_types[i] = tmpstr
    }

    # Store the whitespace padding for a given message
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

    ########

    # Start task: Input and log file cheking.

    ####
    # Check log file
    if (l_logfile == "" || stridx(" ",l_logfile) > 0) {
        l_logfile = gsaoi.logfile
        if (l_logfile == "" || stridx(" ",l_logfile) > 0) {
            l_logfile = "gsaoi.log"

            printlog ("\n"//warn_msg//"both "//task_name//".logfile and "//\
                "gsaoi.logfile are empty.\n"//\
                warn_pad//"Using default file gsaoi.log.",
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
    printlog ("    sciimg    = "//l_sciimg, l_logfile, l_verbose)
    printlog ("    calimg    = "//l_calimg, l_logfile, l_verbose)
    printlog ("    key_check = "//l_key_check, l_logfile, l_verbose)
    printlog ("    logfile   = "//l_logfile, l_logfile, l_verbose)
    printlog ("    verbose   = "//l_verbose, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    # Check input sciimg
    if (l_sciimg == "" || stridx (" ",l_sciimg) > 0) {
        printlog (err_msg//"sciimg is not set properly", \
            l_logfile, verbose+)
        goto crash
    } else if (!imaccess(l_sciimg)) {
        printlog (err_msg//l_sciimg//" does not exist", \
            l_logfile, verbose+)
        goto crash
    }

    # Check input calimg
    if (l_calimg == "" || stridx (" ",l_calimg) > 0) {
        printlog (err_msg//"calimg is not set properly", \
            l_logfile, verbose+)
        goto crash
    } else if (!imaccess(l_calimg)) {
        printlog (err_msg//l_calimg//" does not exist", \
            l_logfile, verbose+)
        goto crash
    }

    ########

    # Start the actual work

    ####
    # Loop over the input science image and get

    sciphu = l_sciimg//"[0]"
    calphu = l_calimg//"[0]"

    # NSCIEXT should be present if prepared
    nsciext = INDEF
    keypar(sciphu, "NSCIEXT", silent-)
    if (!keypar.found) {
        printlog (err_msg//"NSCIEXT keyword not found in "//sciphu, \
            l_logfile, verbose+)
        goto crash
    }
    nsciext = int(keypar.value)

    # Set extension names to copy
    next_names = 3
    ext_names[1] = l_sci_ext
    ext_names[2] = l_var_ext
    ext_names[3] = l_dq_ext

    # Flag to track whether calibration needs to be trimmed
    have_to_copy = no
    update_nextend = yes

    # Create a list to copy then copy... allows different science areas
    for (i = 1; i <= nsciext; i += 1) {
        # Get the name and version of the current extension
        for (j = 1; j <= next_names; j += 1) {

            inextn = l_sciimg//"["//ext_names[j]//","//i//"]"
            calextn = l_calimg//"["//ext_names[j]//","//i//"]"

            if (debug) {
                printlog (inextn//" "//calextn, l_logfile, verbose+)
            }

            # Flags for checking extensions
            doscichk = no
            docalchk = no

            if (imaccess(inextn)) {
                doscichk = yes
            }
            if (imaccess(calextn)) {
                docalchk = yes
            }

            # Check that the same extensions in teh science image exist in the
            # calibration image
            if (doscichk) {
                if (!docalchk) {
                    printlog (err_msg//"No "//ext_names[j]//\
                        " extensions in "//\
                        l_calimg, l_logfile, verbose+)
                    goto crash
                } else if (j > 1) {
                    # Need to keep track of updating NEXTEND
                    update_nextend = no
                }
            } else if (j == 1) {
                # Only fail if sci extensions don't exist
                printlog (err_msg//l_sciimg//" does not contain any "//\
                    ext_names[j]//" extensions", l_logfile, verbose+)
                goto crash
            } else {
                # No var or DQ extensions
                goto SKIP_NEXT
            }

            # Now check the dimensions; assume var and dq
            # dimensions are the same as the science extensions
            # This makes the assumption that the extension versions in each of
            # the files represent the same array.
            if (j == 1) {

                gadimschk (inextn, section="", chkimage=calextn, \
                    key_check=l_key_check, logfile=l_logfile, verbose=debug)

                if (gadimschk.status != 0) {
                    printlog (err_msg//"GADIMSCHK returned a non-zero "//\
                        "status for extensions: "//inextn//" "//calextn, \
                        l_logfile, verbose+)
                    goto crash
                } else if (strstr("[*,*]",gadimschk.out_chkimage) == 0) {
                    # sections do not match
                    have_to_copy = yes
                }

                # Store the section for printing to list for copying
                fparse (gadimschk.out_chkimage)
                sec_to_copy = fparse.section

                # Read the key_check to update output - should be there as
                # it's needed by gadimschk
                keypar (inextn, l_key_check, silent+)
                secvalue = keypar.value
            }
            # Print the file names to a tmp file in case a trimmed copy of the
            # calibration file need be created
            print (calextn//sec_to_copy//" "//tmpout_calfile//" "//\
                secvalue, >> tmp_tocopy)
SKIP_NEXT:
        } # End of loop over j, ext_names
    } # End of loop over i, input science image extensions

    # Copy if needed, a subsection, update output name too.
    if (have_to_copy) {
        imcopy (calphu, tmpout_calfile, verbose=debug)

        if (debug)
            type (tmp_tocopy)

        scanfile = tmp_tocopy
        i = 0
        while (fscan(scanfile, input, output, secvalue) != EOF) {
            i += 1
            imcopy (input, output//"[append]", verbose=debug)

            # Update the value of the output key_check value
            fparse (input)
            ##M BUG if DATASEC is to be used in the future
            gemhedit (output//fparse.ksection, l_key_check, \
                secvalue, "", delete-, upfile="")
        }
        if (update_nextend) {
            gemhedit (tmpout_calfile//"[0]", "NEXTEND", i, "", delete-, \
                upfile="")
        }
        outimg = tmpout_calfile
    } else {
        outimg = l_calimg
    }

    # All good to get here
    l_status = 0
    goto clean

crash:

    outimg = ""
    if (access(tmpout_calfile)) {
        imdelete (tmpout_calfile, verify-, >& "dev$null")
    }

clean:

    scanfile = ""
    delete (tmp_tocopy, verify-, >& "dev$null")

    # Finish time
    date | scan (l_date_struct)
    printlog ("\n"//task_name//" Finished -- "//l_date_struct//"\n", \
        l_logfile, l_verbose)

    status = l_status

    exit_status[1] = "GOOD"
    exit_status[2] = "ERROR"

    printlog (task_name//" exit status:  "//exit_status[status + 1]//\
        ".", l_logfile, verbose+)

    printlog (log_wrap//"\n", l_logfile, l_verbose)

end
