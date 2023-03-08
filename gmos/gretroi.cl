# Copyright(c) 2012-2017 Association of Universities for Research in Astronomy, Inc.

procedure gretroi (inimage)

# This task determines which extensions of MEF GMOS image that contain part of
# the ROI specified in the PHU (a specific ROI can be requested); it also
# determines upon request the area of the data section / image that relate to
# the ROIs. The output can be:
#     (1) An image that just contains the untouched extension that contain part
#         of an ROI
#     (2) An image that contains extensions that make up the ROI(s), the
#         extensions themselves contain only data within the ROI(s).
#     (3) A list of image extensions, with an image section that relates to the
#         ROI(s), related DETSEC and CCDEC for this image section, binning,
#         the ROI that the extension relates too and CRPIX values.
#
#     In the case of 1 and 2, any var_ext or dq_ext extensions will be returned
#     too if fl_vardq+. Also, the PHU will have a timestamp. In addtion, three
#     keywords will be added to the PHU: GRETROI - timestamp; RETROI - which
#     ROI was returned (0=all, then numbered there after); and RETROIT - how
#     the extensions that are returned were returned (EXTENSION||ROI). Any
#     present MDFs will be copied too. PHU counters are updated too.
#
#     In case of 3, only the sci extensions will be returned (it's assumed the
#     calling script will deal with the VAR/DQ extensions itself under the
#     assumption that VAR and DQ planes relate to the SCI extensions.

# INPUTS cannot be GMOSAICed
# One image at a time!

char    inimage     {prompt="Input GMOS image"}
char    outimage    {"", prompt="Output image"}
char    outfile     {"STDOUT", prompt="Output text file"}
char    rawpath     {"", prompt="Path to inimage"}
int     req_roi     {0, min=0, prompt="Requested ROI: 0=all"}
bool    roi_only    {yes, prompt="Trim extensions to contain ROI data only?"}
bool    fl_vardq    {yes, prompt="Return var_ext and dq_ext extensions in image?"}
char    sci_ext     {"SCI", prompt="Name of science extension(s)"}
char    var_ext     {"VAR", prompt="Name of variance extension(s)"}
char    dq_ext      {"DQ", prompt="Name of data quality extension(s)"}
char    mdf_ext     {"MDF", prompt="Mask definition file extension name"}
char    key_detsec  {"DETSEC", prompt="Header keyword for detector section"}
char    key_ccdsec  {"CCDSEC", prompt="Header keyword for CCD section"}
char    key_datasec {"DATASEC", prompt="Header keyword for data section"}
char    key_biassec {"BIASSEC", prompt="Header keyword for overscan section"}
char    key_ccdsum  {"CCDSUM", prompt="Header keyword for CCD binning"}
char    logfile     {"", prompt="Name of logfile to use"}
bool    verbose     {no, prompt="Verbose?"}
int     status      {0, prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

begin

    # Setup local variables for input parameters
    char l_inimage = ""
    char l_outimage = ""
    char l_outfile = ""
    char l_rawpath = ""
    int  l_req_roi = 0
    bool l_roi_only = yes
    bool l_fl_vardq = yes
    char l_sci_ext = ""
    char l_var_ext = ""
    char l_dq_ext = ""
    char l_mdf_ext = ""
    char l_key_detsec = ""
    char l_key_ccdsec = ""
    char l_key_datasec = ""
    char l_key_biassec = ""
    char l_key_ccdsum = ""
    char l_logfile = ""
    bool l_verbose = yes
    int  l_status

    # Set up variables used in code
    char inimg, inphu
    char l_gemindex, l_extn, l_extver, l_omit, tmpsci_list
    char tmpsci_list2, inextn, chk_tasks[3], l_sciext, outextn
    char err_str, warn_str, prt_str, this_task, img_root, curr_task
    char l_key_task, l_key_retroi, l_key_typeret, l_key_extnroi
    char roi_key_detsec, roi_key_ccdsec, roi_key_datasec, tmphaeder
    char l_outsciext, tmpheader
    char ret_key_detsec, ret_key_ccdsec, ret_key_datasec, copy_datasec
    char l_inum_str, l_char_str, out_datasec, type_ret, img_extn_num
    char ext_kernel_info, ret_ccdsec_comm, ret_datasec_comm, ret_detsec_comm
    char roi_detsec_comm, roi_ccdsec_comm, roi_datasec_comm
    char detsec_comm, ccdsec_comm, datasec_comm, tstring
    char in_extname[3], out_extname[3], outphu, dettype

    int  nextend, ii, jj, j, atlocation, num_images, err_len, warn_len, prt_len
    int out_nsci, out_nextend, ntasks, nsci, row_removal_value
    int ndetroi, junk, xbin, ybin, iistart, iistop, curr_roi, num_extn
    int hdroix, hdroixs, hdroiy, hdroiys, comma_pos1, comma_pos2
    int x1offset, x2offset, x2x1offset, x1x2offset, nout_extensions
    int y1offset, y2offset, y2y1offset, y1y2offset

    real l_crpix1, l_crpix2

    bool yvalid, xvalid, missextn, sec_update, reset_extn, clean_header
    bool fl_useprefix, hasmdf, hasvar, hasdq, out_imgset, out_fileset
    bool check_sec, ishamamatsu, already_exists, istrimmed
    bool debug, profile

    struct sdate, pstruct

    # Max of 5 ROIs; max 12 extensions = 60
    int roidetx1[5], roidetx2[5], roidety1[5], roidety2[5]
    char l_detsec, l_ccdsec, l_datasec, detsec[60], ccdsec[60], datasec[60]
    char biassec[60], l_biassec
    int tx1, tx2, ty1, ty2, detx1[60], detx2[60], dety1[60], dety2[60]
    int ccdx1[60], ccdx2[60], ccdy1[60], ccdy2[60]
    int datax1[60], datax2[60], datay1[60], datay2[60]
    int biasx1[60], biasx2[60], biasy1[60], biasy2[60], bias_width[60]
    int data_width[60], data_height[60]

    # Read input parameters
    l_inimage = inimage
    l_outimage = outimage
    l_outfile = outfile
    l_rawpath = rawpath
    l_req_roi = req_roi
    l_roi_only = roi_only
    l_fl_vardq = fl_vardq
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_mdf_ext = mdf_ext
    l_key_detsec = key_detsec
    l_key_ccdsec = key_ccdsec
    l_key_datasec = key_datasec
    l_key_biassec = key_biassec
    l_key_ccdsum = key_ccdsum
    l_logfile = logfile
    l_verbose = verbose

    # Set default parameters
    debug = no
    profile = no
    if (debug) {
        profile = yes
        l_verbose = yes
    }
    l_status = 1 # Default is bad
    already_exists = no
    reset_extn = no
    check_sec = yes

    # General task print string
    this_task = "GRETROI"
    prt_str = this_task//": "
    prt_len = strlen(prt_str)
    prt_str = "\n"//prt_str
    l_inum_str = "%-8s %d \"%s\"\n"
    l_char_str = "%-8s \"%s\" \"%s\"\n"

    # Error string & warning strings
    err_str = "ERROR - "//this_task//": "
    err_len = strlen(err_str)
    err_str = "\nERROR - "//this_task//": "
    warn_str = "WARNING - "//this_task//": "
    warn_len = strlen(warn_str)
    warn_str = "\nWARNING - "//this_task//": "

    # Output file header keywords
    l_key_task = "GRETROI"
    l_key_retroi = "RETROI"
    l_key_typeret = "RETROITY"
    roi_key_detsec = "ROIDSEC"
    roi_key_ccdsec = "ROICSEC"
    roi_key_datasec = "ROIDASEC"
    l_key_extnroi = "EXTNROI"

    if (l_roi_only) {
        type_ret = "ROI"
    } else {
        type_ret = "EXTENSION"
    }

    roi_detsec_comm = "ROI DETSEC"
    roi_ccdsec_comm = "ROI CCDSEC"
    roi_datasec_comm = "ROI datasec"
    detsec_comm = "."
    ccdsec_comm = "."
    datasec_comm = "."

    # Fail tasks
    ntasks = 1
    ntasks = 3
    chk_tasks[1] = "GMOSAIC"
    chk_tasks[2] = "GTILE"
    chk_tasks[3] = this_task

    # Temporary files
    tmpsci_list = mktemp ("tmpscilist")
    tmpsci_list2 = mktemp ("tmpscilist2")
    tmpheader = mktemp ("tmpheader")

    # Store required tasks in memeory
    cache ("gimverify")

    ###########
    # Test the logfile here to start logging:
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
        l_logfile = gmos.logfile
        if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
            l_logfile = "gmos.log"
            printlog (warn_str//"Both gretroi.logfile and "//\
                "gmos.logfile fields are empty.", l_logfile, verbose+)
            printlog ("                   Using default file gmos.log", \
                l_logfile, verbose+)
        }
    }

    if (profile) {
        date ("+%H:%M:%S.%N%n") | scan (sdate)
        printlog ("____"//this_task//" Started: "//sdate, l_logfile, verbose+)
    }

    # Print start time
    gemdate ()
    printlog ("\n"//this_task//" - Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    ################
    # Start checks #
    ################

    #### Perform checks on input/output lists; setup filenames ####

    out_imgset = no
    out_fileset = no

    # Check the input exists
    gimverify (l_rawpath//l_inimage)
    if (gimverify.status != 0) {
        if (gimverify.status == 1) {
            printlog (err_str//l_inimage//" does not exist", \
                l_logfile, verbose+)
        } else {
            printlog (err_str//l_inimage//" is not a MEF", \
                l_logfile, verbose+)
        }
        goto crash
    }

    l_inimage = gimverify.outname//".fits"

    # Check for requested output
    if (l_outimage == "" && l_outfile == "") {
        printlog (err_str//"Either outimage or outfile must be set", \
            l_logfile, verbose+)
        goto crash
    }

    if (l_outimage != "") {
        out_imgset = yes
        if (stridx(" ",l_outimage) > 0) {
            printlog (err_str//"outimage cannot contain spaces", \
                l_logfile, verbose+)
            goto crash
        } else if (imaccess(l_outimage)) {
            printlog (err_str//"outimage "//l_outimage//" already exists", \
                l_logfile, verbose+)
            already_exists = yes
            goto crash
        }

        fparse (l_outimage)
        if (fparse.extension == "") {
            l_outimage = fparse.directory//fparse.root//".fits"
        }
    }

    if (l_outfile != "") {
        out_fileset = yes
        if (stridx(" ",l_outfile) > 0) {
            printlog (err_str//"outfile cannot contain spaces", \
                l_logfile, verbose+)
            goto crash
        } else if (access(l_outfile) && l_outfile != "STDOUT") {
            printlog (err_str//l_outfile//" already exists", \
                l_logfile, verbose+)
            already_exists = yes
            goto crash
        }
    }

    # only one output allowed
    if (out_imgset && out_fileset) {
        if (l_outfile != "STDOUT") {
            printlog (warn_str//"Both outimage and outfile are set. Only an \
                output image will be created", \
                l_logfile, verbose+)
        }
        out_fileset = no
    }

    # If file requested not image switch vardq off
    if (out_fileset && l_fl_vardq) {
        l_fl_vardq = no
    }

    #### Check extension names ####

    # Require that sci_ext is set for output image (if it doesn't exist in
    # image then it will be found...
    if (stridx(" ",l_sci_ext) > 0) {# && out_imgset) {
        printlog (err_str//"sci_ext is not set properly", \
                l_logfile, verbose+)
        goto crash
    }

    if (l_fl_vardq && out_imgset) {
        l_fl_vardq = no
        if (l_var_ext == "" || (stridx(" ",l_var_ext) > 0)) {
            printlog (warn_str//l_var_ext//" is not set correctly", \
                l_logfile, verbose+)
        } else if (l_dq_ext == "" || (stridx(" ",l_dq_ext) > 0)) {
            printlog (warn_str//l_dq_ext//" is not set correctly", \
                l_logfile, verbose+)
        } else {
            l_fl_vardq = yes
        }

        if (!l_fl_vardq) {
            printf ("%"//warn_len//"s%s\n", " ","Setting fl_vardq = no\n") | \
                scan (pstruct)
            printlog (pstruct, l_logfile, verbose+)
        }
    }

    #################
    # End of checks #
    #################

    ##########################
    # Start of the hard work #
    ##########################

    # Set flags
    hasvar = no
    hasdq = no

    inimg = l_inimage
    fparse (inimg)
    img_root = fparse.root//fparse.extension
    inphu = inimg//"[0]"

    printlog (prt_str//" Working on image - "//img_root, l_logfile, l_verbose)

    clean_header = no
    # Check if this file has been preprocessed before
    for (jj = 1; jj <= ntasks; jj += 1) {

        curr_task = chk_tasks[jj]

        keypar (inphu, curr_task, silent+)
        if (keypar.found) {
            if (jj == 1) {
                # GMOSAIC
                printlog (err_str//curr_task//" has already been ran on "//\
                    img_root, l_logfile, verbose+)
                goto crash
            } else if (jj == 2) {
                # GTILE
                check_sec = no
                break
            } else if (jj == 3) {
                # this_task - check requested inputs against work already done

                # Need to check if requesting extensions when already only ROI
                # avaliable && to check if headers require to cleaned up - MS
                keypar (inphu, l_key_typeret, silent+)
                if (keypar.found) {
                    if (keypar.value != type_ret && !l_roi_only) {
                        printlog (warn_str//img_root//" contains only ROI \
                            data but roi_only=no; reseting roi_only=yes.", \
                            l_logfile, verbose+)
                        l_roi_only = yes
                    } else if (keypar.value != type_ret && l_roi_only) {
                        clean_header = yes
                    }
                }

                # Check requested ROI was returned previously if not requesting
                # all ROIs
                if (l_req_roi != 0) {
                    keypar (inphu, l_key_retroi, silent+)
                    if (keypar.found) {
                        if (l_req_roi == int(keypar.value)) {
                            # They match reset l_req_roi to 0
                            l_req_roi = 0
                        } else {
                            printlog (err_str//"Requested ROI previously \
                                removed by "//this_task, l_logfile, verbose+)
                            goto crash
                        }
                    } else {
                        printlog (err_str//"Cannot access "//l_key_retroi//\
                            " keyword in "//inphu, l_logfile, verbose+)
                        goto crash
                    }
                }
            }
        }
    }

    l_gemindex = ""
    l_extn = l_sci_ext
    l_extver = "1-"
    l_omit = ",index"

    if (l_extn == "") {
        l_gemindex = "1-"
        l_extver = ""
        l_omit = ""
    }

    # Set nsci
    nsci = INDEF

    # Get a list of input sci extensions
    gemextn (inimg, check="exists,mef", index=l_gemindex, extname=l_extn, \
        extversion=l_extver, ikparams="", omit="section,params"//l_omit, \
        outfile=tmpsci_list, logfile=l_logfile, glogpars="", verbose=no)

    if (gemextn.status != 0) {
        printlog (err_str//"GEMEXTN retuned a non-zero status", \
            l_logfile, verbose+)
        goto crash
    } else if (l_extn != "" && gemextn.count == 0 &&  \
        gemextn.fail_count == 0) {

        # l_sci_ext doesn't exist in inimg
        delete (tmpsci_list, verify-, >& "dev$null")

        printlog (warn_str//l_extn//" extensions not found in "//\
            img_root//". Trying sci_ext=\"\"", \
            l_logfile, verbose+)

        if (debug) {
            printlog ("____sci_ext = "//l_sci_ext//": gemextn cannot \
                find any of these extensions in "//img_root, \
                l_logfile, verbose+)
        }

        # Re-perform the check

        l_gemindex = "1-"
        l_extn = ""
        l_extver = ""
        l_omit = ""

        gemextn (inimg, check="exists,mef", index=l_gemindex, \
            extname=l_extn, extversion=l_extver, ikparams="", \
            omit="section,params"//l_omit, \
            outfile=tmpsci_list2, logfile=l_logfile, glogpars="", \
            verbose=no)

        if (gemextn.status != 0) {
            printlog (err_str//"GEMEXTN retuned a non-zero status", \
                l_logfile, verbose+)
            goto crash
        } else if (gemextn.count == 0 && gemextn.fail_count == 0) {

            # Something is wrong
            printlog (err_str//img_root//" contains no data - GEMEXTN", \
                l_logfile, verbose+)
            goto crash
        }

        # Only want extensions with no extname / version
        match ("][", tmpsci_list2, stop=yes, print_file_name=no, \
            metacharacters=no, > tmpsci_list)

        # Check the number of extensions
        count (tmpsci_list) | scan (nsci)

        if (isindef(nsci) || nsci <= 0) {

            # Something is wrong
            printlog (err_str//img_root//" contains no data - COUNT", \
                l_logfile, verbose+)
            goto crash
        } else {
            printlog (warn_str//"Reseting sci_ext=\"\" internally", \
                l_logfile, verbose+)
            if (l_sci_ext == "" && out_imgset) {
                printf ("%"//warn_len//"s%s\n", " ", \
                    "Setting output science EXTNAME to SCI") | scan (pstruct)
                printlog (pstruct, l_logfile, verbose+)
                l_sci_ext = "SCI"
                reset_extn = yes
            }
        }

    } else if (l_extn == "" && gemextn.count == 0 &&  \
        gemextn.fail_count == 0) {

        # Something is wrong!

        printlog (err_str//img_root//" contains no data", \
            l_logfile, verbose+)
        goto crash
    } else if (l_extn == "") {
        # Strip out anymore extensions

        # Only want extensions with no extname / version
        match ("][", tmpsci_list, stop=yes, print_file_name=no, \
            metacharacters=no, > tmpsci_list2)

        # Check the number of extensions
        count (tmpsci_list) | scan (nsci)

        if (isindef(nsci) || nsci <= 0) {

            # Something is wrong
            printlog (err_str//img_root//" contains no data", \
                l_logfile, verbose+)
            goto crash
        }

        delete (tmpsci_list, verify-, >& "dev$null")

        # Reset the tmplist variable
        tmpsci_list = tmpsci_list2

        if (l_sci_ext == "" && out_imgset) {
            printlog (prt_str//"Setting output science EXTNAME to SCI", \
                l_logfile, l_verbose)
            l_sci_ext = "SCI"
            reset_extn = yes
        }

    } else {
        # Set the number of extensions (to perform check of var and dq exts
        nsci = gemextn.count
    }

    # Check nsci!
    if (isindef(nsci)) {
        printlog (err_str//"Cannot determine number of science extensions", \
            l_logfile, verbose+)
        goto crash
    }

    # Set the extension names to use for reading inputs and writing outputs
    l_outsciext = l_sci_ext
    if (l_sci_ext != "") {
        l_outsciext = l_outsciext//","
    }

    l_sciext = l_extn
    if (l_extn != "") {
        l_sciext = l_sciext//","
    }

    # Check if the sci_ext name has been reset
    if (l_fl_vardq) {

        if (l_sci_ext != "" && l_extn == "" && !reset_extn) {
            printf ("%"//warn_len//"s%s\n", " ", "Setting fl_vardq = no\n") | \
                scan (pstruct)
            printlog (pstruct, l_logfile, verbose+)
            l_fl_vardq = no
        }
    }

    if (l_fl_vardq) {
        # Check for var / dq extns
        l_gemindex = ""
        l_omit = ""
        l_extn = l_var_ext
        l_extver = "1-"

        for (j = 1; j <= 2; j += 1) {

            if (j == 2) {
                l_extn = l_dq_ext
            }

            gemextn (inimg, check="exists,mef", index=l_gemindex, \
                extname=l_extn, extversion=l_extver, ikparams="", \
                omit="section,params"//l_omit, \
                outfile="dev$null", logfile=l_logfile, glogpars="", \
                verbose=no)

            if (gemextn.status != 0) {
                printlog (err_str//"GEMEXTN retuned a non-zero status", \
                    l_logfile, verbose+)
                goto crash
            } else if (gemextn.count == 0 && gemextn.fail_count == 0) {

                # Something is wrong
                printlog (prt_str//img_root//" contains no "//l_extn//\
                    " data", l_logfile, l_verbose)

            } else if (nsci != gemextn.count && gemextn.count > 0) {
                printlog (err_str//"Number of "//l_var_ext//\
                    " not equal to the number of "//l_sci_ext, \
                    l_logfile, verbose+)
            } else {

                if (l_extn == l_var_ext) {
                    hasvar = yes
                } else if (l_extn == l_dq_ext) {
                    hasdq = yes
                }
            }
        } # End of loop over var / dq

        if (!hasvar && !hasdq) {
            printf ("%"//prt_len//"s%s\n", " ", "Setting fl_vardq=no") | \
                scan (pstruct)
            printlog (pstruct, l_logfile, l_verbose)
        }
    }  # End of var/dq checks

    nextend = 0

    # If creating output image copy out PHU
    if (out_imgset) {

        # Copy out PHU
        imcopy (inphu, l_outimage, verbose-)

        # Set this regardless here
        outphu = l_outimage//"[0]"

        # Check for an MDF
        # Check for var_exts
        l_gemindex = ""
        l_omit = ",index"
        l_extn = l_mdf_ext
        l_extver = "1-"

        gemextn (inimg, check="exists,mef", index=l_gemindex, \
            extname=l_extn, extversion=l_extver, ikparams="", \
            omit="section,params"//l_omit, \
            outfile="dev$null", logfile=l_logfile, glogpars="", \
            verbose=no)

        if (gemextn.status != 0) {
            printlog (err_str//"GEMEXTN retuned a non-zero status", \
                l_logfile, verbose+)
            goto crash
        } else if (gemextn.count > 0) {
            tcopy (l_inimage//"["//l_mdf_ext//"]", l_outimage//\
                "["//l_mdf_ext//"]", verbose-)
            nextend += 1
        }
    }

    #### PHU ROI INFO ####
    # Read the number of ROIs
    keypar (inphu, "DETNROI", silent+)
    if (!keypar.found) {
        printlog (err_str//"DETNROI keyword not found in "//inphu, \
            l_logfile, verbose+)
        goto clean
    }
    ndetroi = int(keypar.value)

    # Check if it's Hamamatsu data
    keypar (inphu, "DETTYPE", silent+)
    if (!keypar.found) {
        printlog (err_str//"DETTYPE keyword not found in "//inphu, \
            l_logfile, verbose+)
        goto clean
    }
    dettype = keypar.value

    ishamamatsu = no
    row_removal_value = 0
    if (dettype == "S10892" || dettype == "S10892-N") {
        ishamamatsu = yes
        row_removal_value = 48
    }

    # Check if the image has been trimmed
    istrimmed = no
    keypar (inphu, "TRIMMED", silent+)
    if (keypar.found) {
        istrimmed = yes
    }

    # Read the dimensions of the ROIs and convert to equivalent of DETSEC
    # Read binning from first sci extension
    keypar (inimg//"["//l_sciext//"1]", "CCDSUM", silent+)
    if (!keypar.found) {
        printlog (err_str//"Cannot access CCDSUM "//\
            "keyword in "//inimg//"["//l_sciext//"1]", \
             l_logfile, verbose+)
        goto clean
    } else {
        print (keypar.value) | scan (xbin, ybin)
    }

    # Used in multiple places!
    iistart = 1
    iistop = ndetroi

    if (l_req_roi != 0) {

        if (l_req_roi > ndetroi) {
            printlog (err_str//"Requested ROI is greater than the number of \
                available ROIs.", l_logfile, verbose+)
            goto crash
        }
        iistart = l_req_roi
        iistop = iistart
    }

    # Loop over desired ROIs
    for (ii = iistart; ii <= iistop; ii += 1) {

        hdroix = INDEF
        hdroixs = INDEF
        hdroiy = INDEF
        hdroiys = INDEF

        hselect (inphu, "DETRO"//ii//"X,DETRO"//ii//"XS,"//\
            "DETRO"//ii//"Y,DETRO"//ii//"YS", yes) | \
             scan (hdroix, hdroixs, hdroiy, hdroiys)

        if (isindef(hdroix)) {
            printlog (err_str//"Cannot access DETRO"//ii//"X "//\
                "keyword in "//inphu, \
                l_logfile, verbose+)
            goto clean
        } else {
            roidetx1[ii] = hdroix
        }

        if (isindef(hdroixs)) {
            printlog (err_str//"Cannot access DETRO"//ii//"XS "//\
                "keyword in "//inphu, \
                l_logfile, verbose+)
            goto clean
        } else {
            roidetx2[ii] = (hdroixs * xbin) + (roidetx1[ii] - 1)
        }

        if (isindef(hdroiy)) {
            printlog (err_str//"Cannot access DETRO"//ii//"Y "//\
                "keyword in "//inphu, \
                l_logfile, verbose+)
            goto clean
        } else if (ishamamatsu) {
            if  (l_roi_only && (hdroiy <= row_removal_value / ybin)) {
                # Adjust roidety1 for bottom unused rows - do not update the
                # extension sections, let the code determining if extension is
                # within an ROI or not do the trimming
                roidety1[ii] = row_removal_value + 1
            } else {
                roidety1[ii] = hdroiy
            }
        } else {
            roidety1[ii] = hdroiy
        }

        if (isindef(hdroiys)) {
            printlog (err_str//"Cannot access DETRO"//ii//"YS "//\
                "keyword in "//inphu, \
                l_logfile, verbose+)
            goto clean
        } else {
            roidety2[ii] = (hdroiys * ybin) + (roidety1[ii] - 1)
        }

    } # End of loop reading ROIs

    #### Work on image extensions ####

    # Reset nsci to 1
    nsci = 1

    # Loop over tempsci_list
    scanfile = tmpsci_list
    j = 0
    nout_extensions = 0
    while (fscan(scanfile,inextn) != EOF) {

       if (debug) {
            printlog ("\n%%%%"//inextn, l_logfile, verbose+)
       }

        # Parse the extension name / index
        fparse (inextn)
        if (fparse.ksection != "") {
            ext_kernel_info = substr(fparse.ksection,2,\
                strlen(fparse.ksection)-1)

            comma_pos1 = stridx(",",ext_kernel_info)
            comma_pos2 = strldx(",",ext_kernel_info)

            if (comma_pos1 != comma_pos2) {
                # There is extra information in the kernel section
                tstring = substr(ext_kernel_info,\
                    stridx(",",ext_kernel_info)+1,strlen(ext_kernel_info))

                ext_kernel_info = substr(ext_kernel_info,1,comma_pos1)//\
                    substr(tstring,1,stridx(",",tstring)-1)
            }

            img_extn_num = ext_kernel_info
            comma_pos1 = stridx(",",img_extn_num)
            if (comma_pos1 != 0) {
                img_extn_num = substr(img_extn_num,comma_pos1+1,\
                    strlen(img_extn_num))
            }
        } else {
            ext_kernel_info = str(fparse.cl_index)
            img_extn_num = ext_kernel_info
        }

        # Update output image index / recording information
        j += 1

        # Set default values for values to be read
        l_detsec = INDEF
        l_ccdsec = INDEF
        l_datasec = INDEF
        l_biassec = INDEF
        l_crpix1 = INDEF
        l_crpix2 = INDEF

        # Read all of the required information
        keypar (inextn, l_key_detsec, silent+)
        if (keypar.found) {
            l_detsec = keypar.value
        }

        keypar (inextn, l_key_ccdsec, silent+)
        if (keypar.found) {
            l_ccdsec = keypar.value
        }

        keypar (inextn, l_key_datasec, silent+)
        if (keypar.found) {
            l_datasec = keypar.value
        }

        keypar (inextn, l_key_biassec, silent+)
        if (keypar.found) {
            l_biassec = keypar.value
        }

        keypar (inextn, "CRPIX1", silent+)
        if (keypar.found) {
            l_crpix1 = real(keypar.value)
        }

        keypar (inextn, "CRPIX2", silent+)
        if (keypar.found) {
            l_crpix2 = real(keypar.value)
        }

        # DETSEC
        if (isindef(l_detsec)) {
            printlog (err_str//l_key_detsec//" not found. "//\
                "Exiting.", l_logfile, verbose+)
            goto crash
        }

        detsec[j] = l_detsec
        junk = fscanf (l_detsec, "[%d:%d,%d:%d]", tx1, tx2, ty1, ty2)
        detx1[j] = tx1
        detx2[j] = tx2
        dety1[j] = ty1
        dety2[j] = ty2

        # CCDSEC
        if (isindef(l_ccdsec)) {
            printlog (err_str//l_key_ccdsec//" not found. "//\
                "Exiting.", l_logfile, verbose+)
            goto crash
        }

        ccdsec[j] = l_ccdsec
        junk = fscanf (l_ccdsec, "[%d:%d,%d:%d]", tx1, tx2, ty1, ty2)
        ccdx1[j] = tx1
        ccdx2[j] = tx2
        ccdy1[j] = ty1
        ccdy2[j] = ty2

        # DATASEC
        if (isindef(l_datasec)) {
            printlog (err_str//l_key_datasec//" not found. "//\
                "Exiting.", l_logfile, verbose+)
            goto crash
        }

        datasec[j] = l_datasec
        junk = fscanf (l_datasec, "[%d:%d,%d:%d]", tx1, tx2, ty1, ty2)
        datax1[j] = tx1
        datax2[j] = tx2
        datay1[j] = ty1
        datay2[j] = ty2
        data_width[j] = datax2[j] - (datax1[j] - 1)
        data_height[j] = datay2[j] - (datay1[j] - 1)

        # Only need biassec if not trimmed (for crpix1 / hamamamtsu data)
        if (!istrimmed) {
            # BIASSEC
            if (isindef(l_biassec)) {
                printlog (warn_str//l_key_biassec//" not found. "//\
                    "Assuming it doesn't exist for "//inextn, \
                    l_logfile, verbose+)
                # Create fake biassec
                l_biassec = "[1:0,0:0]"
            }

            biassec[j] = l_biassec
            junk = fscanf (l_biassec, "[%d:%d,%d:%d]", tx1, tx2, ty1, ty2)
            biasx1[j] = tx1
            biasx2[j] = tx2
            biasy1[j] = ty1
            biasy2[j] = ty2
            bias_width[j] = biasx2[j] - (biasx1[j] - 1)
        }

        # CRPIX1
        if (isindef(l_crpix1)) {
            printlog (warn_str//"CRPIX1 not found. "//\
                "Seeting to 0. Output file WCS may be wrong.", \
                 l_logfile, verbose+)
            l_crpix1 = 0
        } else if (!istrimmed && biasx1[j] == 1 && l_crpix1 > 0 ) {
            # Update crpix1 if bias is on the left and crpix1 is > 0.
            l_crpix1 -= bias_width[j]
        }

        # CRPIX 2
        if (isindef(l_crpix2)) {
            printlog (warn_str//"CRPIX2 not found. "//\
                "Seeting to 0. Output file WCS may be wrong.", \
                 l_logfile, verbose+)
            l_crpix2 = 0
        }

        # If Hamamatsu data - check that DETSEC is above row_removal_value
        # Remove if neccessary
        ##M Removed this to let the next block of code checking for section
        ##M relation to ROIs trim the bottom 48 rows
#        if (ishamamatsu) {
#            if (dety1[j] <= row_removal_value) {
#                datay1[j] = datay1[j] + \
#                    ((row_removal_value - (dety1[j] - 1)) / ybin)
#                print (":::::_____(.) "//datay1[j])
#                data_height[j] = datay2[j] - (datay1[j] - 1)
#                dety1[j] = row_removal_value + 1
#                ccdy1[j] = row_removal_value + 1
#                l_crpix2 -= ((row_removal_value - dety1[j]) / ybin)
#            }
#        }

        # Default current ROI is 1
        curr_roi = 1

        if (!check_sec) {
            goto MISS_SEC_CHK
        }

        # Check if the array is in the ROIs
        for (ii = iistart; ii <= iistop; ii += 1) {

            if (debug) {
                printlog ("++++"//ii, l_logfile, verbose+)
            }

            curr_roi = ii

            sec_update = yes
            yvalid = no
            xvalid = no
            missextn = no

            x1offset = roidetx1[ii] - detx1[j]
            x2offset = roidetx2[ii] - detx2[j]
            x2x1offset = roidetx2[ii] - detx1[j]
            x1x2offset = roidetx1[ii] - detx2[j]

            y1offset = roidety1[ii] - dety1[j]
            y2offset = roidety2[ii] - dety2[j]
            y2y1offset = roidety2[ii] - dety1[j]
            y1y2offset = roidety1[ii] - dety2[j]

            # Check y values
            if (((y1offset <= 0 && y2offset <= 0) || \
                (y1offset >= 0 && y2offset >= 0) || \
                (y1offset >= 0 && y2offset <= 0) || \
                (y1offset <= 0 && y2offset >= 0)) && \
                y2y1offset >=0 && y1y2offset <= 0) {

                 yvalid = yes
             }

            # First pass of x values
            if (x2x1offset >= 0 && x1x2offset <=0) {
                xvalid = yes
            }

            if (debug) {
                printlog ("____"//\
                    detx1[j]//" "//roidetx1[ii]//" "//detx2[j]//\
                    " "//roidetx2[ii]//" "//dety1[j]//" "//roidety1[ii]//\
                    " "//dety2[j]//" "//roidety2[ii], l_logfile, verbose+)
                printlog (x1offset//" "//x2offset//" "//y1offset//\
                    " "//y2offset, l_logfile, verbose+)
                printlog (x2x1offset//" "//x1x2offset//" "//y2y1offset//\
                    " "//y1y2offset, l_logfile, verbose+)
                printlog ("====xvalid: "//xvalid//" yvalid: "//yvalid, \
                    l_logfile, verbose+)
            }


            if (!xvalid || !yvalid) {
                 missextn = yes
            } else if (x1offset <= 0 && x2offset >= 0 && y1offset <= 0 \
                    && y2offset >= 0) {
                # Extension is soley within ROI
                if (debug) {
                    printlog ("____Extension is soley within ROI", \
                        l_logfile, verbose+)
                }

                sec_update = no
                break

            } else if (x1offset >= 0 && x2offset <= 0 && y1offset >= 0 \
                    && y2offset <= 0) {

                # ROI soley in extension
                if (debug) {
                    printlog ("____ROI soley in extension", \
                        l_logfile, verbose+)
                }
                break

            } else if (x1offset >= 0 && x2offset >= 0 && yvalid) {

                if (debug) {
                    printlog ("____Case 1", \
                        l_logfile, verbose+)
                }
                # Case 1
                break

            } else if (x1offset <= 0 && x2offset <= 0  && yvalid) {

                if (debug) {
                    printlog ("____Case 2", \
                        l_logfile, verbose+)
                }

                # Case 2
                break

            } else if (x1offset <= 0 && x2offset >= 0 && yvalid) {

                if (debug) {
                    printlog ("____Case 3", \
                        l_logfile, verbose+)
                }

                # Case 3
                break

            } else if (x1offset >= 0 && x2offset <= 0 && yvalid) {

                if (debug) {
                    printlog ("____Case 4", \
                        l_logfile, verbose+)
                }

                # Case 4
                break

            } else {
                missextn = yes
            }

            if (missextn && ii == iistop) {

                if (debug) {
                    printlog ("____extn: "//inextn//" does not \
                        lie in a ROI. Skipping.", l_logfile, verbose+)
                }
                j -= 1
                goto NEXTEXTN
            }
        } # End of loop over requested ROIs - removing unneed extensions

        # Determine region of input extension datasec is part of ROI
        # Calculate updated CCDSEC and DETSEC for this region

        if (sec_update && l_roi_only) {
            x1offset = max (0,x1offset)
            x2offset = min (0,x2offset)
            y1offset = max (0,y1offset)
            y2offset = min (0,y2offset)

            if (debug) {
                printlog (x1offset//" "//x2offset//" "//y1offset//\
                    " "//y2offset, l_logfile, verbose+)
            }

            # DETSEC
            detx1[j] = detx1[j] + x1offset
            detx2[j] = detx2[j] + x2offset
            dety1[j] = dety1[j] + y1offset
            dety2[j] = dety2[j] + y2offset
            l_detsec = "["//detx1[j]//":"//detx2[j]//","//\
                dety1[j]//":"//dety2[j]//"]"

            # CCDSEC
            ccdx1[j] = ccdx1[j] + x1offset
            ccdx2[j] = ccdx2[j] + x2offset
            ccdy1[j] = ccdy1[j] + y1offset
            ccdy2[j] = ccdy2[j] + y2offset
            l_ccdsec = "["//ccdx1[j]//":"//ccdx2[j]//","//\
                ccdy1[j]//":"//ccdy2[j]//"]"

            # DATASEC
            datax1[j] = datax1[j] + (x1offset / xbin)
            datax2[j] = datax2[j] + (x2offset / xbin)
            datay1[j] = datay1[j] + (y1offset / ybin)
            datay2[j] = datay2[j] + (y2offset / ybin)
            l_datasec = "["//datax1[j]//":"//datax2[j]//","//\
                datay1[j]//":"//datay2[j]//"]"
            data_width[j] = datax2[j] - (datax1[j] - 1)
            data_height[j] = datay2[j] - (datay1[j] - 1)

            l_crpix1 -= x1offset
            l_crpix2 -= y1offset
        }

MISS_SEC_CHK:

        nout_extensions += 1

        # l_detsec, l_ccdsec and l_datasec now refer to the ROI section within
        # the extension

        # Set the output keywords
        if (l_roi_only) {
            copy_datasec = l_datasec
            out_datasec = "[1:"//data_width[j]//",1:"//data_height[j]//"]"
            ret_key_detsec = l_key_detsec
            ret_key_ccdsec = l_key_ccdsec
            ret_key_datasec = l_key_datasec
            ret_detsec_comm = detsec_comm
            ret_ccdsec_comm = ccdsec_comm
            ret_datasec_comm = datasec_comm
        } else {
            copy_datasec = ""
            out_datasec = l_datasec
            ret_key_detsec = roi_key_detsec
            ret_key_ccdsec = roi_key_ccdsec
            ret_key_datasec = roi_key_datasec
            ret_detsec_comm = roi_detsec_comm
            ret_ccdsec_comm = roi_ccdsec_comm
            ret_datasec_comm = roi_datasec_comm
        }

        if (debug) {
            printlog ("____copy_datasec: "//copy_datasec, \
                l_logfile, verbose+)
        }

        # Write the output
        if (out_imgset) {
            # Write to output image

            # Output extension
            outextn = l_outimage//"["//l_outsciext//j//",append]"

            num_extn = 1
            in_extname[num_extn] = inextn
            out_extname[num_extn] = outextn

            if (hasvar) {
                num_extn += 1
                in_extname[num_extn] = img_root//\
                    "["//l_var_ext//","//img_extn_num//",append]"
                out_extname[num_extn] = l_outimage//\
                    "["//l_var_ext//","//j//",append]"
            }

            if (hasdq) {
                num_extn += 1
                in_extname[num_extn] = img_root//\
                    "["//l_dq_ext//","//img_extn_num//",append]"
                out_extname[num_extn] = l_outimage//\
                    "["//l_dq_ext//","//j//",append]"
            }

            # Create file to update headers
            print ("default_pars update+ add+ show- verify- before=\"\" \
                after=\"\" delete-", >> tmpheader)

            if (clean_header) {
                # Delete previous hedaers written by this task
                print (roi_key_detsec//" delete+", >> tmpheader)
                print (roi_key_ccdsec//" delete+", >> tmpheader)
                print (roi_key_datasec//" delete+", >> tmpheader)
            }

            printf (l_char_str, ret_key_detsec, l_detsec, ret_detsec_comm, \
                >> tmpheader)

            printf (l_char_str, ret_key_ccdsec, l_ccdsec, ret_ccdsec_comm, \
                >> tmpheader)

            printf (l_char_str, ret_key_datasec, out_datasec, \
                ret_datasec_comm,
                >> tmpheader)

            printf (l_inum_str, l_key_extnroi, curr_roi, \
                "ROI extn belongs to", >> tmpheader)

            if (l_roi_only) {
                printf (l_char_str, "TRIMSEC", copy_datasec, \
                    "Trimmed section(s)", >> tmpheader)
            }

            # Loop over extensions to copy

            for (jj = 1; jj <= num_extn; jj += 1) {

                # Copy the requested data from input to output
                imcopy (in_extname[jj]//copy_datasec, out_extname[jj], \
                    verbose-)

                gemhedit (out_extname[jj], "", "", "", delete-, \
                    upfile=tmpheader)

                if (l_roi_only) {
                    # Reset the LTV1 and 2 keywords using wcsreset
                    # wcs=physicsal
                    # This is so logical and physical coordinates relate to
                    # the new file not any images they may have come from
                    wcsreset (out_extname[jj], wcs="physical",\
                        verbose-)
                }
            }

            delete (tmpheader, verify-, >& "dev$null")

        } else if (out_fileset) {
            # Write to text file

            if (j == 1 && l_outfile == "STDOUT") {

                printf ("\n#%-6s %-18s %-21s %-21s %3s %-12s %-12s\n", \
                    "EXTN", "DATASEC", "DETSEC", \
                    "CCDSEC", "ROI", "CRPIX1", "CRPIX2", >> l_outfile)
            }

            printf ("%-7s %18s %21s %21s %3d %.11f %.11f\n", \
                ext_kernel_info, l_datasec, \
                l_detsec, l_ccdsec, curr_roi, l_crpix1, l_crpix2, >> l_outfile)
        }

NEXTEXTN:

    } # End of loop over input extensions

    # Check for at least one extension lieing in an ROI
    if (nout_extensions == 0) {
        printlog ("\nERROR - GRETROI: No extensions lie within any of the \
            defined ROIs", l_logfile, verbose+)
        goto crash
    }

    # Update output counters
    out_nsci = j
    out_nextend = j
    if (hasvar) {
        out_nextend += j
    }
    if (hasdq) {
        out_nextend += j
    }

    # Update the PHU if rquired
    if (out_imgset) {
        # Create file to update headers
        print ("default_pars update+ add+ before=\"\" \
           after=\"\" delete-", >> tmpheader)

        if (l_req_roi != 0) {

            printf (l_inum_str, "DETNROI", 1, \
                "No. regions of interest after "//this_task, >> tmpheader)
            printf (l_inum_str, "DETRO1X", hdroix, \
                "ROI 1 X start after "//this_task, >> tmpheader)
            printf (l_inum_str, "DETRO1XS", hdroixs, \
                "ROI 1 X size after "//this_task, >> tmpheader)
            printf (l_inum_str, "DETRO1Y", hdroiy, \
                "ROI 1 Y start after "//this_task, >> tmpheader)
            printf (l_inum_str, "DETRO1YS", hdroiys, \
                "ROI 1 Y size after "//this_task, >> tmpheader)
        }

        printf (l_inum_str, l_key_retroi, l_req_roi, \
            "ROI returned by "//l_key_task//" \(0\=all\)", >> tmpheader)

        printf (l_inum_str, "NSCIEXT", out_nsci, \
            "Number of science extensions", >> tmpheader)

        printf (l_inum_str, "NEXTEND", out_nsci, \
            "Number of extensions", >> tmpheader)

        printf (l_inum_str, l_key_retroi, l_req_roi, \
            "ROI returned by "//l_key_task, >> tmpheader)

        printf (l_char_str, l_key_typeret, type_ret, "Data returned by "//\
            l_key_task, >> tmpheader)

        if (l_roi_only) {
            printf (l_char_str, "TRIMMED", "yes", "Overscan section trimmed", \
                >> tmpheader)
        }

        gemdate(zone="UT")

        printf (l_char_str, l_key_task, gemdate.outdate, \
            "UT Time stamp for "//l_key_task, >> tmpheader)

        printf (l_char_str, "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI", >> tmpheader)

        gemhedit (outphu, "", "", "", delete-, \
            upfile=tmpheader)

        delete (tmpheader, verify-, >& "dev$null")

    }

    # All must be good to get here
    l_status = 0

    goto clean

crash:

    l_status = 1
    if (!already_exists && out_imgset) {
        imdelete (l_outimage, verify-, >& "dev$null")
    }
    if (!already_exists && out_fileset) {
        delete (l_outfile, verify-, >& "dev$null")
    }
clean:

    # Reset the status
    status = l_status

    scanfile = ""

    delete (tmpsci_list//","//tmpsci_list2//","//tmpheader, \
        verify-, >& "dev$null")

    gemdate ()
    printlog (this_task//" - Finished: "//gemdate.outdate//"\n", \
        logfile, l_verbose)

    if (status == 0) {
        printlog (this_task//" - Exit status: GOOD\n", logfile, l_verbose)
    } else {
        printlog (this_task//" - Exit status: ERROR\n", logfile, verbose+)
    }

    if (profile) {
        date ("+%H:%M:%S.%N%n") | scan (sdate)
        printlog ("____"//this_task//" Finished: "//sdate, l_logfile, verbose+)
    }

end
