# Copyright(c) 2013 Association of Universities for Research in Astronomy, Inc.

######################################################
procedure gemdqexpand (inimage, outimage)

######################################################
# Notes:

# This task assumes that all DQ values are > 0. Therefore if using signed
# values only > 0 values will be returned and not quite all the possible bit
# information will be returned. For example: int32 will decompose values from
# 0 upto [(2**31) - 1]. Where as ushort will decompose values from 0 upto
# [(2**16)-1].

######################################################
# Parameters

char    inimage       {prompt="Input image"}
char    outimage      {prompt="Input image"}
char    logfile       {"",prompt="Logfile"}
char    dq_ext        {"DQ",prompt="Data quality extension name"}
bool    fl_real       {no,prompt="Convert output to real?"}
bool    verbose       {yes,prompt="Verbose?"}
int     status        {0,prompt="Task status. 0=GOOD, 1=ERROR"}

######################################################

begin

    ########
    ####
    # Local versions of input parameters
    char l_inimage, l_outimage, task_name, err_msg, warn_msg, task_msg, mdelim
    char msg_types[4], pad_types[4], tmpstr, log_wrap, l_logfile, l_dq_ext
    char err_pad, warn_pad, task_pad, mpad_delim, pixtype, exit_status[2]
    char tmpoutimg, curr_outext, curr_appoutext, zero_pixtype

    int n_msg_types, lwrap_len, l_status, intype, type_index
    int ii, max_value, curr_bit, bit_base, mvalue, ext_num, max_nbits

    bool debug, l_verbose, l_fl_real

    int    npixtypes=12
    char   l_pixtype[12] = "bool","char","short","int","long","real","double","complex","pointer","struct","ushort","ubyte"
    # l_pixtype integers -- > 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 (from iraf$unix/lib/iraf.h)

    # Allowed inputs anf their respctive maximum number of bits (note < is
    # issued later rather than <= when determining the maximum bit value)
    int    n_allowed_ptypes=5
    int    allowed_ptypes[5]=3,4,5,11,12
    int    nbits[5]=15,31,31,16,16

    struct l_date_struct

    ####
    # Assign inputs to local variables
    l_inimage = inimage
    l_outimage = outimage
    l_dq_ext = dq_ext
    l_logfile = logfile
    l_verbose = verbose
    l_fl_real = fl_real

    ####
    # Set task and package names
    task_name = strupr("gemdqexpand")

    ####
    # Default values
    l_status = 1 # Assume it will fail reset just before end if not
    debug = no
    bit_base = 2

    # Reset verbosity if debug set to true
    if (debug)
        l_verbose = yes

    ####
    # Temporray files
    tmpoutimg = mktemp ("tmpoutimg")//".fits"

    ####
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

    ########

    # Start task: Input and log file cheking.

    ####
    # Check log file
    if (l_logfile == "" || stridx(" ",l_logfile) > 0) {
        l_logfile = gemtools.logfile
        if (l_logfile == "" || stridx(" ",l_logfile) > 0) {
            l_logfile = "gemtools.log"

            printlog ("\n"//warn_msg//"both "//task_name//".logfile and "//\
                "gemtools.logfile are empty.\n"//\
                warn_pad//"Using default file gemtools.log.",
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
    printlog ("    inimage    = "//l_inimage, l_logfile, l_verbose)
    printlog ("    outimage   = "//l_outimage, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    ####
    # Check for output image
    if (l_outimage == "" || stridx (" ",l_outimage) > 0) {
        printlog (err_msg//"outimage is not set properly", \
            l_logfile, verbose+)
        goto crash
    } else if (imaccess(l_outimage)) {
        printlog (err_msg//l_outimage//" already exists", \
            l_logfile, verbose+)
        goto crash
    }

    # Check input image
    if (l_inimage == "" || stridx (" ",l_inimage) > 0) {
        printlog (err_msg//"inimage is not set properly", \
            l_logfile, verbose+)
        goto crash
    } else if (!imaccess(l_inimage)) {
        printlog (err_msg//l_inimage//" does not exist", \
            l_logfile, verbose+)
        goto crash
    }

    # Check that the image is an integer - obtain the max value
    keypar (l_inimage, "i_pixtype", silent+)
    if (!keypar.found) {
        printlog (err_msg//"Cannot access pixtype keyword for "//l_inimage, \
            l_logfile, verbose+)
        goto crash
    } else {
        intype = int(keypar.value)
        if (debug)
            printlog (mdelim//"intype: "//intype, l_logfile, verbose+)

        # Determine if it's a valid DQ value
        # Only allow
        for (ii = 1; ii <= n_allowed_ptypes; ii += 1) {
            if (intype == allowed_ptypes[ii]) {
                type_index = ii
                break
            } else if (ii == n_allowed_ptypes) {
                printlog (err_msg//l_inimage//" has a bitxpix of "//intype//\
                    " which is not an integer pixel type", \
                    l_logfile, verbose+)
                goto crash
            }
        }
    }

    # Check dq_ext
    if (l_dq_ext == "" || stridx(" ",l_dq_ext) > 0) {
        printlog (err_msg//"dq_ext is not set properly", \
            l_logfile, verbose+)
        goto crash
    }

    ########

    # Set the output pixtype
    pixtype = l_pixtype[intype]
    if (debug)
        printlog (mdelim//"Pixtype: "//pixtype, l_logfile, verbose+)

    # Find the highest value of the DQ plane
    max_value = INDEF

    imstatistics (l_inimage, fields="max", lower=INDEF, upper=INDEF, \
        nclip=0, lsigma=3., usigma=3., binwidth=0.1, format-, cache-) | \
        scan (max_value)

    if (isindef(max_value)) {
        printlog (err_msg//"Cannot determine maximum value in "//l_inimage, \
            l_logfile, verbose+)
        goto crash
    }

    ext_num = 1
    max_nbits = nbits[type_index]
    if (max_value > max_nbits && max_value <= max_nbits) {
        if (debug)
            printlog (mdelim//"Max value > than max_nits", \
                l_logfile, verbose+)
        max_nbits = 2 * max_nbits
    }

    if (debug) {
        printlog (mdelim//"Max values: "//max_nbits//" "//max_value, \
            l_logfile, verbose+)
    }

    zero_pixtype = pixtype
    if (l_fl_real) {
        zero_pixtype = "real"
    }
    # Create a zero bit mask for masking faked data
    imexpr ("a & b", l_outimage//"["//l_dq_ext//",0,append]", \
        l_inimage, 0, dims="auto", intype="auto", outtype=zero_pixtype, \
        refim="auto", bwidth=0, btype="constant", bpixval=0., \
        rangecheck=yes, verbose-, exprdb="none")

    # Loop over the bits that make up the max_value
    for (ii = 0; ii < max_nbits; ii += 1) {

        curr_bit = bit_base ** ii

        if (debug)
            printlog (mdelim//"ii: "//ii//" curr_bit = "//curr_bit, \
                l_logfile, verbose+)

        curr_outext = l_outimage//"["//l_dq_ext//","//curr_bit//"]"
        curr_appoutext = l_outimage//"["//l_dq_ext//","//curr_bit//",append]"

        if (curr_bit > max_value) {
            goto CREATED_OUTPUT
        }

        imexpr ("a & b", tmpoutimg, l_inimage, curr_bit, dims="auto", \
            intype="auto", outtype=pixtype, \
            refim="auto", bwidth=0, btype="constant", bpixval=0., \
            rangecheck=yes, verbose-, exprdb="none")

        # Check if there are any values in this extension
        # Find the highest value of the DQ plane
        mvalue = INDEF

        imstatistics (tmpoutimg, fields="max", lower=INDEF, \
            upper=INDEF, nclip=0, lsigma=3., usigma=3., binwidth=0.1, \
            format-, cache-) | scan (mvalue)

        if (isindef(mvalue)) {
            printlog (err_msg//"Cannot determine maximum value in "//\
                l_inimage,
                l_logfile, verbose+)
            goto crash
        } else if (mvalue != 0) {
            if (l_fl_real) {
                chpixtype (tmpoutimg, curr_appoutext, newpixtype="real", \
                    oldpixtype="all", verbose-)
            } else {
                imcopy (tmpoutimg, curr_appoutext, verbose-)
            }
            ext_num += 1
        }

        delete (tmpoutimg, verify-, >& "dev$null")
    }

CREATED_OUTPUT:

    # Update the output phu with the number of extensions
    gemhedit (l_outimage//"[0]", "NEXTEND", ext_num, "Number of extensions", \
        delete-, upfile="")

    # Final header updates
    gemdate (zone="UT")
    gemhedit (l_outimage//"[0]", substr(task_name,1,8), gemdate.outdate, \
        "UT Time stamp for "//task_name, delete-)
    gemhedit (l_outimage//"[0]", "GEM-TLM", gemdate.outdate, \
        "UT Last modification with GEMINI", delete-)

    l_status = 0
    goto clean

crash:
;
clean:

    delete (tmpoutimg, verify-, >& "dev$null")

    status = l_status

    # Finish time
    date | scan (l_date_struct)
    printlog ("\n"//task_name//" Finished -- "//l_date_struct//"\n", \
        l_logfile, l_verbose)

    status = l_status

    exit_status[1] = "GOOD"
    exit_status[2] = "ERROR"

    printlog (task_name//" exit status:  "//exit_status[status + 1]//\
        ".", l_logfile, l_verbose)

    printlog (log_wrap//"\n", l_logfile, l_verbose)

end
