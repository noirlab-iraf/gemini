# Copyright(c) 2011-2017 Association of Universities for Research in Astronomy, Inc.

procedure gtile (inimages)

# Tile the amplifiers of the GMOS input images into their respective CCDs and
# return the input image as phu plus MDF if presenet, plus 3 (each CCD)
# versions each type of extension (e.g., SCI, VAR, DQ) with appropriate header
# keywords updated. All supplied CCDs are returned by default (out_ccds="all").
# An individual CCD science frame (with no phu) is returned when fl_stats_only
# is yes and out_ccds != "all". If fl_stats_only=yes and fl_tile_det=yes
# the entire detector tiled together can be returned. There are checks for
# flag clashes, where there are clashes flags get reset. - MS

# This task is hidden and designed for use when doing things
# like imstat calls on CCD2 in giflat and gemlocal tasks
# It makes the following assumptions:
#     The detector orientation is the same for North and South; i.e.,
#         CCD1 detsec X1 starts at 1 and CCD3 detsec X2 ends at 6144

# Changed the behaviour to allow tiling of all extensions into detector -
#    meaning fl_tile_det can work with fl_stats_only-

# Updated to handle multiple ROIs

char    inimages    {prompt="GMOS images to tile (single/list)"}
char    outimages   {"", prompt="Output image name or list of names"}
char    outpref     {"ti", prompt="prefix for output images"}
char    out_ccds    {"all", prompt="Which CCDs to return. (all|1|2|3)"}
bool    ret_roi     {yes, prompt="Return only ROI data"}
int     req_roi     {0, min=0, prompt="ROI to return. (0=all)"}
bool    fl_stats_only {no, prompt="Return only the science data for use with statistics tasks"}
bool    fl_tile_det {no, prompt="To tile all three CCDs"}
bool    fl_app_rois {no, prompt="Append tiled ROIs?"}
bool    fl_pad      {no, prompt="Include chip gaps in tiled detector image?"}
real    sci_padval  {0., prompt="Value of chipgap pixels when fl_pad=yes for sci_ext"}
real    var_padval  {0., prompt="Value of chipgap pixels when fl_pad=yes for var_ext"}
real    dq_padval   {16., prompt="Value of chipgap pixels when fl_pad=yes for dq_ext"}
real    sci_fakeval {0., prompt="Value of any faked pixles for sci_ext"}
real    var_fakeval {0., prompt="Value of any faked pixels for var_ext"}
real    dq_fakeval  {16., prompt="Value of any faked pixels for dq_ext"}
char    chipgap     {"default", prompt="Chipgap width in pixels"}
char    sci_ext     {"SCI", prompt="Name of science extension(s)"}
char    var_ext     {"VAR", prompt="Name of variance extension(s)"}
char    dq_ext      {"DQ", prompt="Name of data quality extension(s)"}
char    mdf_ext     {"MDF", prompt="Mask definition file extension name"}
char    key_detsec  {"DETSEC", prompt="Header keyword for detector section"}
char    key_ccdsec  {"CCDSEC", prompt="Header keyword for CCD section"}
char    key_datasec {"DATASEC", prompt="Header keyword for data section"}
char    key_biassec {"BIASSEC", prompt="Header keyword for overscan section"}
char    key_ccdsum  {"CCDSUM", prompt="Header keyword for CCD binning"}
char    rawpath     {"", prompt="Path to raw data file(s)"}
char    logfile     {"", prompt="Name of logfile to use"}
bool    fl_verbose  {no, prompt="Verbose"}
int     status      {0, prompt="Exit status (0=good)"}
struct  *scanfile   {prompt="For internal use only"}
struct  *scanfileII {prompt="For internal use only"}

begin

    # Setup local variables for input parameters
    char    l_inimages = ""
    char    l_outimages = ""
    char    l_outpref = ""
    char    l_out_ccds = ""
    bool    l_fl_stats_only
    bool    l_fl_tile_det
    bool    l_fl_pad
    real    l_sci_padval = 0.
    real    l_var_padval = 0.
    real    l_dq_padval = 0.
    real    l_sci_fakeval = 0.
    real    l_var_fakeval = 0.
    real    l_dq_fakeval = 0.
    char    l_chipgap = ""
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_mdf_ext = ""
    int     l_req_roi = 0
    bool    l_ret_roi = yes
    bool    l_fl_app_rois = no
    char    l_key_detsec = ""
    char    l_key_ccdsec = ""
    char    l_key_datasec = ""
    char    l_key_biassec = ""
    char    l_key_ccdsum = ""
    char    l_rawpath = ""
    char    l_logfile = ""
    bool    l_fl_verbose

    # Max 5 ROIs; max 12 amps = 6.
    char    filename, inimage_files[200], outimg_files[200]
    char    datasec[60], detsec[60], ccdsec[60], ccdsum
    char    l_detsec, l_ccdsec, l_datasec, extn
    char    extn_list, in_extn, imtile_list, tstring, imextn
    char    out_pixtype, inlist, outlist, next_inp, l_sciext
    char    out_ampsname, amplifiername[60]
    char    l_sec, l_key, in_tile_file, listfile, instrument, dettype
    char    rphend, pathtest, tmpinfile, tmpsec, default_gemh_pars, tmpheader
    char    tmpccd_file[3], name_chipgap[3]
    char    tmpdetsec, tmpccdsec, tmpdatasec, tmplist, tmplist2, tmpgretroi
    char    insciext, tmpamp_file, tmpamp_stem, ampsec, curr_extn
    char    l_key_typeret, l_key_extnroi, type_ret
    char    roi_key_detsec, roi_key_ccdsec, roi_key_datasec, l_gtcomm
    char    l_CCD_gttype, l_APP_gttype, l_DET_gttype, l_key_gttype, l_gttype
    int     tx1, tx2, ty1, ty2, ampsorder[60], amps2order[60]
    int     detx1[60], detx2[60], dety1[60], dety2[60]
    int     ccdx1[60], ccdx2[60], ccdy1[60], ccdy2[60]
    int     datax1[60], datax2[60], datay1[60], datay2[60]
    int     data_width[60], ccd_height[2,3]
    int     data_height[60], extn_roi[60], in_order[60]
    int     maxfiles, junk, in_pixtype, num_images, total_num_amplifiers
    int     number_of_ampsperccd[3,5] # Max 3 CCDs 5 ROIs
    int     i, j, k, ii, jj, kk, refextn[5], naxis1, tmpint, num_ccds[5]
    int     jjstart, jjstop, jstart, current_extn, height_iccd
    int     comma_pos, curr_roi, ref_loop_max, refextn_ccd[5]
    int     det_start[3], det_stop[3], ordered_refextn[5], num_ext_roi[5]
    int     num_sci_extns, num_extns, num_of_this_extension
    int     detroiy_start, row_removal_value, row_removal_value_binned
    int     xbin, ybin, image_width, image_height, site, iccd
    int     secx1, secx2, secy1, secy2, atlocation, raw_width[2,3], gapwidth
    int     total_num_sciext, iistart, iistop, fkccdy1, fkccdx1, fkccdx2
    int     ndetroi, ccddist, detdist, namps, ampwidth, y1det, pix_type_ext
    int     roidetx1, roidetx1s, roidety1, roidety1s, l_ncols, num_dets
    int     num_ext_last,outnum, wcs_ref, fake_width, num_ext_names, type_ext
    real    fake_values[3], pad_values[3]
    real    crpix1_in[60], l_crpix1, l_crpix2
    real    l_cd1_1, l_cd1_2, l_cd2_1, l_cd2_2
    real    crpix2_in[60], crpix1_out, crpix2_out, dety_start
    bool    ret_one_roi, tile_amps, ret_sec, clean_header
    bool    missextn, fl_useprefix, mdf_exists, dotile, ref_ext_used
    bool    update_rois, fake_ccd_created, ishamamatsu, chk_roi
    bool    one_key_found, debug, profile
    struct  sdate

    char    array_keys[3]
    char    l_pixtype[12] = "bool","char","short","int","long","real","double","complex","pointer","struct","ushort","ubyte"
    # l_pixtype integers -- > 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 (from iraf$unix/lib/iraf.h)

    # Set up parameters caches of the following tasks
    cache ("gimverify", "gemextn", "gemdate", "imjoin", "gretroi")

    # Set internal debug flag
    debug = no
    profile = no
    if (debug) {
        profile = yes
        l_fl_verbose = yes
    }

    if (profile) {
        date ("+%H:%M:%S.%N%n") | scan (sdate)
        printlog ("____GTILE Started: "//sdate, l_logfile, verbose+)
    }

    # Set default values
    maxfiles = 200
    total_num_amplifiers = 0
    total_num_sciext = 0
    status = 0
    l_sciext = ""
    ishamamatsu = no
    default_gemh_pars = "default_pars update+ add+ before=\"\" \
        after=\"\" delete-"

    # The total number of det_starts and stops
    num_dets = 3

    # Set det limits for start and stop of each CCD
    det_start[1] = 1
    det_start[2] = 2049
    det_start[3] = 4097

    det_stop[1] = 2048
    det_stop[2] = 4096
    det_stop[3] = 6144

    # Define the chip gaps for each site and detector in unbinned pixels
    # raw_width[site, iccd]
    # site: 1 = GMOS-N; 2 = GMOS-S
    # iccd: 1 = EEV; 2 = e2vDD; 3 = Hamamatsu
    raw_width[1,1] = 37
    raw_width[1,2] = 37
    raw_width[1,3] = 67  #KL check this
    raw_width[2,1] = 37
    raw_width[2,3] = 61 ##M CHIP_GAP

    # Define the CCD height for each detector type - always trimmed
    ccd_height[1,1] = 4608 # 4644 Physical size
    ccd_height[1,2] = 4608 # 4644 Physical size
#    ccd_height[1,3] = 4176 # 4224 Physical size
    ccd_height[1,3] = 4224 # 4224 Physical size  #KL check this
    ccd_height[2,1] = 4608 # 4644 Physical size
    ccd_height[2,3] = 4224 # 4224 Physical size

    # Read input parameters
    l_inimages = inimages
    l_outimages = outimages
    l_outpref = outpref
    l_fl_stats_only = fl_stats_only
    l_fl_tile_det = fl_tile_det
    l_fl_pad = fl_pad
    l_chipgap = chipgap
    junk = fscan (out_ccds, l_out_ccds)
    l_sci_ext = sci_ext
    junk += fscan (var_ext, l_var_ext)
    junk += fscan (dq_ext, l_dq_ext)
    junk += fscan (mdf_ext, l_mdf_ext)
    junk += fscan (key_detsec, l_key_detsec)
    junk += fscan (key_ccdsec, l_key_ccdsec)
    junk += fscan (key_datasec, l_key_datasec)
    junk += fscan (key_biassec, l_key_biassec)
    junk += fscan (key_ccdsum, l_key_ccdsum)
    l_sci_padval = sci_padval
    l_var_padval = var_padval
    l_dq_padval = dq_padval
    l_sci_fakeval = sci_fakeval
    l_var_fakeval = var_fakeval
    l_dq_fakeval = dq_fakeval
    l_rawpath = rawpath
    l_logfile = logfile
    l_fl_verbose = fl_verbose
    l_fl_app_rois = fl_app_rois
    l_ret_roi = ret_roi
    l_req_roi = req_roi

    tmpinfile = mktemp ("tmpinfile")//".fits"
    tmpamp_stem = mktemp ("tmpamp")//".fits"
    tmplist = mktemp("tmplist")
    tmplist2 = mktemp("tmp2list")
    tmpgretroi = mktemp("tmpgretroi")
    tmpheader = mktemp("tmpheader")

    # For updating headers later on
    array_keys[1] = l_key_datasec
    array_keys[2] = l_key_ccdsec
    array_keys[3] = l_key_detsec

    # Need to check if headeres need cleaning from previous runs through
    # GRETROI
    l_key_typeret = "RETROITY"
    l_key_extnroi = "EXTNROI"
    roi_key_detsec = "ROIDSEC"
    roi_key_ccdsec = "ROICSEC"
    roi_key_datasec = "ROIDASEC"

    # Keyword, values and comment for the type of tiling performed
    # This is required to allow a user to tile an image into it's constiuent
    # arrays and then be able to tile that image into it's detector - MS
    l_key_gttype = "GTILETYP"
    l_CCD_gttype = "ARRAY"
    l_APP_gttype = "ROI APPENDED"
    l_DET_gttype = "DETECTOR"
    l_gtcomm = "Type of gtiling performed"

    if (l_ret_roi) {
        type_ret = "ROI"
    } else {
        type_ret = "EXTENSION"
    }
    clean_header = no

    ###########
    # Test the logfile here to start logging:
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
        l_logfile = gmos.logfile
        if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
            l_logfile = "gmos.log"
            printlog ("GTILE - WARNING: Both gtile.logfile and "//\
                "gmos.logfile fields are empty.", l_logfile, verbose+)
            printlog ("                   Using default file gmos.log", \
                l_logfile, verbose+)
        }
    }

    # Print start time
    gemdate (zone="UT")
    printlog ("\nGTILE - Started: "//gemdate.outdate//"\n", \
        l_logfile, l_fl_verbose)

    ################
    # Start checks #
    ################

    # Check the fscanned inputs == 9
    if (junk != 9) {
        printlog ("ERROR - GTILE: Either out_ccds or one or more of the key \
            name parameters are an empty string", l_logfile, verbose+)
        goto crash
    }

    # Perform checks on input/output lists; setup filenames;

    # Create temporary inlist and outlist files
    inlist = mktemp("tmpinlist")
    outlist = mktemp("tmpoutlist")

    #### Rawpath ####
    if ((l_rawpath != "") && (stridx(" ", l_rawpath) == 0)) {
        # Check it has a trailing slash after expaning any environmental
        # variables
        rphend = substr(l_rawpath,strlen(l_rawpath),strlen(l_rawpath))
        if (rphend == "$") {
            show (substr(l_rawpath,1,strlen(l_rawpath)-1)) | scan (pathtest)
            rphend = substr(pathtest,strlen(pathtest),strlen(pathtest))
        }
        if (rphend != "/") {
            l_rawpath = l_rawpath//"/"
        }
        if (!access(l_rawpath)) {
            printlog ("ERROR - GIREDUCE: Cannot access rawpath: "//l_rawpath, \
                l_logfile, verbose+)
            goto crash
        }
    }

    atlocation = 0
    atlocation = stridx("@", l_inimages)
    if (atlocation > 0) {
        listfile = substr(l_inimages,(atlocation + 1),strlen(l_inimages))
        # Check to se if file exists
        if (!access(listfile)) {
            printlog ("ERROR - GTILE: Cannot access file: \""//listfile//\
                "\". Exiting",\
                l_logfile, verbose+)
            goto crash
        }
    }

    # Check for an input list
    if (atlocation >= 1) {
        # Input is a list
        sections (l_inimages, > inlist)
    } else {
        # Images not in an @list
        files (l_inimages, sort=no, > inlist)
    } # End of input list checks

    # Set up parameters for reading input files
    scanfile = inlist
    filename = ""
    num_images = 1

    # Loop over input list of images
    while (fscan(scanfile,filename) != EOF) {
        # Verify image
        gimverify (l_rawpath//filename)

        if (gimverify.status > 0) {
             # Image bad/doesn't exist goto to crash
             printlog ("\nERROR - GTILE: File: "//filename//\
                 " does not exist or is bad", l_logfile, verbose+)
             goto crash
        }
        inimage_files[num_images] = gimverify.outname//".fits"
        num_images += 1

    } # End of while loop verifying images and populating inimage_files[]

    # Fix the num_images counter to be correct
    num_images -= 1

    # Set default flag for use of the prefix
    fl_useprefix = yes

    # Check for an output list
    if (l_outimages != "" && l_outimages != " ") {
        if (stridx("@",l_inimages) >= 1) {
            # Output is a list
            sections (l_outimages, > outlist)
        } else {
            # Images not in an @list
            files (l_outimages, sort=no, > outlist)
        }

        # Check the number of output files is equal to the number of inputs
        sections ("@"//outlist, option="nolist")
        if (num_images == int(sections.nimages)) {
            fl_useprefix = no
        } else {
            delete (outlist, verify-, >& "dev$null")
        }
        # End of if statement checking number of in and out files

    } # End of output list checks

    # Check if prefix is to be used
    if (fl_useprefix) {
        # Set list to loop over to be the input list
        outlist = inlist
    }

    # Set up parameters for reading output file list
    scanfile = outlist
    filename = ""
    i = 1

    # Loop over outlist of images to populate outimg_files and check they
    # don't exist
    while (fscan(scanfile,filename) != EOF) {

        # Reset output filename if prefix is to be used
        if (fl_useprefix) {
            fparse (filename)
            filename = l_outpref//fparse.root
        }

        # Verify image
        gimverify (filename)
        if (gimverify.status != 1) {
             # If image doesn't exit
             printlog ("\nERROR - GTILE: File: "//filename//\
                 " already exists", l_logfile, verbose+)
             # Change status to bad
             goto crash
        }
        outimg_files[i] = gimverify.outname//".fits"
        i += 1
    } # End of while loop verifying and populating outimg_files[]

    ###########
    # Check some of the flags for clashes
    if (l_fl_stats_only && l_out_ccds == "all" && !l_fl_tile_det) {
        printlog ("\nGTILE - WARNING: stats_only flag set to yes, out_ccds "//\
            "is set to all", l_logfile, verbose+)
        printlog ("                 and the tile_det flag is set to no.", \
            l_logfile, verbose+)
        printlog ("                 **Resetting fl_tile_det to yes. Only "//\
            "the science frames", l_logfile, verbose+)
        printlog ("                  will be returned tiled into the "//\
            "the available size of the entire detector.\n", \
            l_logfile, verbose+)
        l_fl_tile_det = yes
    }

# Removed to allow tiling of all extnames into a detector - MS
#    if (!l_fl_stats_only && l_fl_tile_det) {
#        printlog ("\nGTILE - WARNING: fl_stats_only flag set to no and "//\
#            "fl_tile_det set to yes.", l_logfile, verbose+)
#        printlog ("                 **Resetting fl_stats_only to yes. "//\
#            "Only the science frames", l_logfile, verbose+)
#        printlog ("                 will be returned tiled into the "//\
#            "available size of the entire detector.", l_logfile, verbose+)
#        l_fl_stats_only = yes
#    }

    if (!l_fl_stats_only && l_out_ccds != "all" && !l_fl_tile_det) {
        printlog ("\nGTILE - WARNING: out_ccds not set to all, but "//\
            "fl_stats_only flag is no", l_logfile, verbose+)
        printlog ("                 and the tile_det flag is set to no.", \
            l_logfile, verbose+)
        printlog ("                 **Resetting out_ccds to all. All "//\
            "the extension types", l_logfile, verbose+)
        printlog ("                 will be returned tiled into the "//\
            "the available size of the entire detector.\n", \
            l_logfile, verbose+)
        l_out_ccds = "all"
    }

    if (l_out_ccds != "all" && l_fl_tile_det) {
        printlog ("\nGTILE - WARNING: out_ccds not set to all, but fl_tile"//\
            "_det flag is yes", l_logfile, verbose+)
        printlog ("                 **Resetting fl_tile_det to no.\n", \
            l_logfile, verbose+)
        l_fl_tile_det = no
    }

    # Check for conflicting flags
    if (l_fl_pad && (l_out_ccds != "all" && !l_fl_tile_det)) {

        printlog ("\nGTILE - WARNING: fl_pad is set to yes but one or more \
             of the following parameters ", \
            l_logfile, verbose+)
        printlog ("                 are not set correctly:", \
            l_logfile, verbose+)
        printlog ("                     fl_stats_only, out_ccds or \
            fl_tile_det", l_logfile, verbose+)
        printlog ("                 They are required to be set as follows"//\
            " to use the fl_pad=yes",
            l_logfile, verbose+)
        printlog ("                 option: ",
            l_logfile, verbose+)
        printlog ("                     fl_stats_only=yes, out_ccds=\"all\" \
            and fl_tile_det=yes", l_logfile, verbose+)
        printlog ("                 **Resetting fl_pad to no.", \
            l_logfile, verbose+)
        l_fl_pad = no
    }

    #### Appending ROIs ####
    # fl_app_rois has presidence over fl_tile_det - for puposes of tiling the
    # correct way
    if (l_fl_app_rois && !l_fl_tile_det) {
        l_fl_tile_det = yes
    }

    #### Set the value of the l_key_gttype keyword
    if (l_fl_app_rois) {
        l_gttype = l_APP_gttype
    } else if (l_fl_tile_det) {
        l_gttype = l_DET_gttype
    } else {
        l_gttype = l_CCD_gttype
    }

    #### Chipgap related ####
    # Check chipgap input
    if ((l_chipgap == "") || (stridx(" ",l_chipgap) > 0)) {
        printlog ("ERROR - GTILE: chipgap parameter not set. Exiting", \
            l_logfile, verbose+)
        goto crash
    } else if (l_chipgap != "default") {
        # Check the input is a number
        gemisnumber (l_chipgap, ttest="decimal", verbose=no)
        if (gemisnumber.status == 0) {
            if (gemisnumber.fl_istype) {
                gapwidth = int(l_chipgap)
            } else {
                printlog ("ERROR - GTILE: chipgap parameter is not a \
                    value. Exiting", l_logfile, verbose+)
                goto crash
            }
        } else {
            printlog ("ERROR - GTILE: gemisnumber returned a non-zero \
                status. Exiting", l_logfile, verbose+)
            goto crash
        }
    }

    #################
    # End of checks #
    #################

    ##########################
    # Start of the hard work #
    ##########################

    # Input files are in inimage_files[]
    # Output files are in outimages_files[]
    # Total number of input images is num_files [loop < num_files]

    # Printlog warning on usage
    printlog ("\nGTILE - WARNING: GTILE only tiles the data sections of "//\
        "each type of "//\
        "extension", l_logfile, l_fl_verbose)
    printlog ("                 No overscan sections will be included in "//\
        "the output image\n", l_logfile, l_fl_verbose)

    # Loop over input images
    for (i=1; i<=num_images; i+=1) {

        # Initiate the list of types of extensions to loop over
        extn_list = mktemp("tmpextn_list")

        # Read the instrument name
        keypar (inimage_files[i]//"[0]", "INSTRUME", silent+)
        instrument = keypar.value
        if (instrument == "GMOS-N" || instrument == "GMOS") {
            # GN
            site = 1
        } else if (instrument == "GMOS-S") {
            # GS
            site = 2
        } else {
            printlog ("ERROR - GTILE: Unrecognizable INSTRUME keyword value"//\
                instrument//" for image "//inimage_files[i]//". Exiting", \
            l_logfile, verbose+)
            goto crash
        }

        # Read the detector type
        keypar (inimage_files[i]//"[0]", "DETTYPE", silent+)
        dettype = keypar.value
        if (dettype == "S10892" || dettype == "S10892-N") {
            # Hamamatsu CCDs
            ##M Remove in the future?
            ishamamatsu = yes
            row_removal_value = 48
            iccd = 3
        } else if (dettype == "SDSU II e2v DD CCD42-90") {
            # e2vDD CCDs
            iccd = 2
        } else if (dettype == "SDSU II CCD") {
            # EEV CCDs
            iccd = 1
        } else {
            printlog ("ERROR - GTILE: Unrecognizable DETTYPE keyword value"//\
                dettype//" for image "//inimage_files[i]//". Exiting", \
            l_logfile, verbose+)
            goto crash
        }

        # Check if headers need cleaning up after previous GRETROI runs
        keypar (inimage_files[i]//"[0]", l_key_typeret, silent+)
        if (keypar.found) {
            if (keypar.value != type_ret && !l_ret_roi) {
                printlog ("WARNING - GTILE: "//inimage_files[i]//" only ROI \
                    data but roi_only=no; reseting roi_only=yes.", \
                    l_logfile, verbose+)
                    l_ret_roi = yes
            } else if (keypar.value != type_ret && l_ret_roi) {
                clean_header = yes
            }
        }

        # Set mdf_exists flag
        mdf_exists = no
        total_num_sciext = INDEF
        num_ext_names = 1
        fake_values[num_ext_names] = l_sci_fakeval
        pad_values[num_ext_names] = l_sci_padval

        gemextn (inimage_files[i], check="exists,mef", process="expand", \
            index="", extname=l_sci_ext, extversion="1-", ikparams="", \
            omit="", replace="", outfile="dev$null", logfile=l_logfile, \
            glogpars="", verbose=debug)

        if (gemextn.status != 0) {
            printlog ("ERROR - GTILE: GEMEXTN returned a non-zero status", \
                l_logfile, verbose+)
            goto crash
        } else if (gemextn.count > 0) {
            total_num_sciext = gemextn.count
        } else {

            l_sci_ext = ""

            # Determine the number of amplifiers
            gemextn (inimage_files[i], check="exists,mef", process="expand", \
                index="1-", extname=l_sci_ext, extversion="", ikparams="", \
                omit="section", replace="", outfile=tmplist, \
                logfile=l_logfile, glogpars="", verbose=debug)

            if (gemextn.status != 0) {
                printlog ("ERROR - GTILE: GEMEXTN returned a non-zero \
                    status", l_logfile, verbose+)
            } else if (gemextn.count > 0) {

                # Only want extensions with no extname / version
                match ("][", tmplist, stop=yes, print_file_name=no, \
                    metacharacters=no, > tmplist2)

                count (tmplist2) | scan (total_num_sciext)
            }

            if (isindef(total_num_sciext)) {
                printlog ("ERROR - GTILE: Cannot determine the number \
                    of extensions", l_logfile, verbose+)
                goto crash
            }
        }
        delete (tmplist//","//tmplist2, verify-, >& "dev$null")

        # Set up list of extensions to loop over
        l_sciext = l_sci_ext

        if (l_sciext != "") {
            l_sciext = l_sciext//","
        }
        print (l_sciext, > extn_list)

        # Check if stats_only, if not -  see if vardq planes and MDF exist
        # If so, set up the list to loop over properly
        if (!l_fl_stats_only) {

            ##M Make this check them separately - add all parameters
            # Check for VARDQ planes exist
            # Determine the number of amplifiers
            gemextn (inimage_files[i], check="exists,mef", process="expand", \
                index="", extname=l_var_ext, extver="1-", ikparams="", \
                omit="section", replace="", outfile="dev$null", \
                logfile=l_logfile, glogpars="", verbose=debug)

            if (gemextn.count > 0) {

                printlog ("GTILE - Variance planes found", \
                    l_logfile, l_fl_verbose)

                if (gemextn.count != total_num_sciext) {
                    printlog ("ERROR - GTILE: The number of \""//l_var_ext//\
                        "\" do not match the number of \""//l_sci_ext//\
                        "\" extnsions", \
                        l_logfile, verbose+)
                    goto crash
                }

                print (l_var_ext//",", >> extn_list)
                num_ext_names +=1
                fake_values[num_ext_names] = l_var_fakeval
                pad_values[num_ext_names] = l_var_padval
            }

            gemextn (inimage_files[i], check="exists,mef", process="expand", \
                index="", extname=l_dq_ext, extver="1-", ikparams="", \
                omit="section", replace="", outfile="dev$null", \
                logfile=l_logfile, glogpars="", verbose=debug)

            if (gemextn.count > 0) {
                printlog ("GTILE - Data quality planes found", \
                    l_logfile, l_fl_verbose)

                if (gemextn.count != total_num_sciext) {
                    printlog ("ERROR - GTILE: The number of \""//l_dq_ext//\
                        "\" do not match the number of \""//l_sci_ext//\
                        "\" extensions", \
                        l_logfile, verbose+)
                    goto crash
                }

                print (l_dq_ext//",", >> extn_list)
                num_ext_names +=1
                fake_values[num_ext_names] = l_dq_fakeval
                pad_values[num_ext_names] = l_dq_padval
            }

            # Find any MDFs
            gemextn (inimage_files[i], check="exists,mef", process="expand", \
                index="", extname=l_mdf_ext, extver="1-", ikparams="", \
                omit="section", replace="", outfile="dev$null", \
                logfile=l_logfile, glogpars="", verbose=debug)

            if (gemextn.count > 0) {
                printlog ("GTILE - MDF found", l_logfile, l_fl_verbose)
                mdf_exists = yes
            }
        } # End of if l_fl_stats_only flag check

        # Check for 0 amplifiers
        if (total_num_sciext < 1) {
            printlog ("ERROR - GTILE: The file "//inimage_files[i]//\
                " could not be tiled.", l_logfile, verbose+)
            printlog("                 No image extensions were found", \
                l_logfile, verbose+)
            printlog("                 Exiting.", l_logfile, verbose+)

            goto crash
        }

        # Call gretroi to create list of requested extensions
        gretroi (inimage_files[i], outimage="", outfile=tmpgretroi, \
            rawpath="", req_roi=l_req_roi, roi_only=l_ret_roi, \
            fl_vardq=no, sci_ext=l_sci_ext, var_ext=l_var_ext, \
            dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, \
            key_detsec=l_key_detsec, key_ccdsec=l_key_ccdsec, \
            key_datasec=l_key_datasec, key_ccdsum=l_key_ccdsum, \
            key_biassec=l_key_biassec, logfile=l_logfile, verbose=debug)

        if (gretroi.status != 0) {
            printlog ("ERROR - GTILE: GRETROI returned a non-zero status", \
                l_logfile, verbose+)
            goto crash
        } else if (!access(tmpgretroi)) {
            printlog ("ERROR - GTILE: Cannot access output from GRETROI", \
                l_logfile, verbose+)
            goto crash
        } else {
            type (tmpgretroi, >> l_logfile)
        }

        # Read the binning for this image
        keypar (inimage_files[i]//"["//l_sciext//"1]",
            l_key_ccdsum, silent+)
        if (keypar.found) {
            print (keypar.value) | scanf ("%d %d", xbin, ybin)
            if (debug) {
                printlog ("____XBIN: "//xbin//" YBIN: "//ybin, l_logfile, \
                    verbose+)
            }
        } else {
            xbin = 1
            ybin = 1
        }

        # Read the number of ROIs
        keypar (inimage_files[i]//"[0]", "DETNROI", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GTILE: DETNROI keyword not found in "//\
                inimage_files[i]//"[0]", l_logfile, verbose+)
            goto crash
        } else {
            ndetroi = int(keypar.value)

            # Check for previous gtile passes on image
            keypar (inimage_files[i]//"[0]", "GTILE", silent+)
            if (keypar.found) {

                # Check the previous type of tiling
                keypar (inimage_files[i]//"[0]", l_key_gttype, silent+)
                if (keypar.found) {

                    # If CCD and obviously previously GTILEd reset ndetroi to 1
                    # to allow tiling into a detector
                    if (keypar.value == l_CCD_gttype) {
                        ndetroi = 1
                    }
                } else {
                    printlog ("ERROR - GTILE: "//l_key_gttype//" not found "//\
                        "in "//inimage_files[i]//"[0]", l_logfile, verbose+)
                    goto crash
                }
            }

            if (debug) {
                printlog ("____NDETROI: "//ndetroi, l_logfile, verbose+)
            }
        }

        update_rois = no
        ret_one_roi = no
        # Read the ROIs if needed
        if (l_req_roi != 0) {
            ret_one_roi = yes
            if (ndetroi != 1) {
                # If ndetroi != 1 mode likely not been ran through gretroi
                # before. Set flag to update ROIs at end
                update_rois = yes
                ndetroi = 1
            }
        }

        if (l_fl_app_rois) {
            ref_loop_max = ndetroi
        } else {
            ref_loop_max = 1
        }

        if (update_rois) {

            ii = l_req_roi

            # Read the X start position
            keypar (inimage_files[i]//"[0]", "DETRO"//ii//"X", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GTILE: Cannot access DETRO"//ii//"X "//\
                    "keyword in "//inimage_files[i]//"[0]", \
                    l_logfile, verbose+)
                goto crash
            } else {
                roidetx1 = int(keypar.value)
            }

            # Read the X size - convert to coordinate
            keypar (inimage_files[i]//"[0]", "DETRO"//ii//"XS", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GTILE: Cannot access DETRO"//ii//"XS "//\
                    "keyword in "//inimage_files[i]//"[0]", \
                    l_logfile, verbose+)
                goto crash
            } else {
                roidetx1s = int(keypar.value)
            }

            # Read the Y start position
            keypar (inimage_files[i]//"[0]", "DETRO"//ii//"Y", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GTILE: Cannot access DETRO"//ii//"Y "//\
                    "keyword in "//inimage_files[i]//"[0]", \
                    l_logfile, verbose+)
                goto crash
            } else {
                roidety1 = int(keypar.value)
            }

            # Read the Y size - convert to coordinate
            keypar (inimage_files[i]//"[0]", "DETRO"//ii//"YS", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GTILE: Cannot access DETRO"//ii//"YS "//\
                    "keyword in "//inimage_files[i]//"[0]", \
                    l_logfile, verbose+)
                goto crash
            } else {
                roidety1s = int(keypar.value)
            }

        } # End of loop reading required ROIs

        # Need NAMPS if returning more than one ROI and not appending them
        if (!ret_one_roi && !l_fl_app_rois && ndetroi > 1) {
         # Read the number of ampliers per CCD
            keypar (inimage_files[i]//"[0]", "NAMPS", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GTILE: Cannot find NAMPS keyword in "//\
                    inimage_files[i]//"[0]", l_logfile, verbose+)
                goto crash
            }
            namps = int (keypar.value)
        }

        # Set image height - default if full CCD height - sections are used
        # later on
        image_height = int(ccd_height[site,iccd]/ybin)

        ##M This may require an update for GMOS-N updates to Hamamatsu CCDs,
        ##  i.e., upgraded from the e2vDD CCDs
        dety_start = 1
        if ((site == 1 || site == 2) && iccd == 3) {
            keypar (inimage_files[i]//"[0]", "TRIMMED", silent+)
            if (keypar.found || l_ret_roi) {
                image_height -= int (row_removal_value / ybin)
                dety_start = int (row_removal_value / ybin) + 1
            }
        }
        height_iccd = image_height

        if (debug) {
            printlog ("____image_height: "//image_height, \
                l_logfile, verbose+)
        }

        # To get to here with l_fl_pad true all of the correct flags must have
        # been set appropriately - create full-frame version then use section
        if (l_fl_pad) {
            # Set the gap width if default - already set if user defined
            if (l_chipgap == "default") {
                gapwidth = int(raw_width[site,iccd]/xbin)
            }

            for (ii = 1; ii <= num_ext_names; ii += 1) {
                if (debug)
                    print ("II="//ii//" "//gapwidth//" "//image_height)

                name_chipgap[ii] = mktemp ("tmpfile_of_chipgap")//".fits"

                # Create the image that will be the chip gap when tiled
                mkpattern (name_chipgap[ii], pattern="constant", \
                    option="replace", v1=pad_values[ii], v2=0, pixtype="real",\
                    ncols=gapwidth, nlines=image_height)
            }
        }

        # Initialise refextn and number of extensions per roi
        for (j = 1; j <= ref_loop_max; j += 1) {
            refextn[j] = 0
            num_ext_roi[j] = 0
        }

        #### Read the output from GRETROI ####

        # Initialise variables
        imextn = INDEF
        l_datasec = INDEF
        l_detsec = INDEF
        l_ccdsec = INDEF
        curr_roi = INDEF
        l_crpix1 = INDEF
        l_crpix2 = INDEF

        j = 0
        scanfileII = tmpgretroi
        while (fscan(scanfileII,imextn,l_datasec,l_detsec,l_ccdsec,curr_roi,\
            l_crpix1,l_crpix2) != EOF) {

            # Read the DETSEC, CCDSEC, CRPIX1, DATASEC and CRPIX values
            if (isindef(imextn) || isindef(l_datasec) || isindef(l_detsec) \
                || isindef(l_ccdsec) || isindef(curr_roi) || \
                isindef(l_crpix1) || isindef(l_crpix2)) {

                printlog ("ERROR - GTILE cannot read information from \
                    GRETROI output", l_logfile, verbose+)
                goto crash
            }

            # Set the extension properly
            comma_pos = stridx(",",imextn)
            tstring = imextn
            if (comma_pos > 0) {
                tstring = substr(tstring,comma_pos+1,strlen(tstring))
            }
            junk = fscan(tstring,kk)

            insciext = inimage_files[i]//"["//imextn//"]"
            j += 1

            if (debug) {
                printlog (insciext//" "//kk//" "//j, \
                    l_logfile, verbose+)
            }

            # DETSEC
            detsec[j] = l_detsec
            junk = fscanf (l_detsec, "[%d:%d,%d:%d]", tx1, tx2, ty1, ty2)
            detx1[j] = tx1
            detx2[j] = tx2
            dety1[j] = ty1
            dety2[j] = ty2

            # DATASEC
            datasec[j] = l_datasec
            junk = fscanf (l_datasec, "[%d:%d,%d:%d]", tx1, tx2, ty1, ty2)
            datax1[j] = tx1
            datax2[j] = tx2
            datay1[j] = ty1
            datay2[j] = ty2
            data_width[j] = tx2 - (tx1 -1)
            data_height[j] = datay2[j] - (datay1[j] - 1)

            if (debug) {
                printlog ("____Width of data section for "//\
                    inimage_files[i]//"["//l_sciext//j//"] is: "//\
                    data_width[j], logfile, verbose+)
            }

            # CCDSEC
            ccdsec[j] = l_ccdsec
            junk = fscanf (l_ccdsec, "[%d:%d,%d:%d]", tx1, tx2, ty1, ty2)
            ccdx1[j] = tx1
            ccdx2[j] = tx2
            ccdy1[j] = ty1
            ccdy2[j] = ty2

            if (debug) {
                printlog ("____"//l_detsec//" "//l_ccdsec//" "//l_datasec//\
                    " "//l_crpix1//" "//l_crpix2, \
                    logfile, verbose+)
            }

            # CRPIX1 & CRPIX2 (these have been adjust for biassec by this point
            crpix1_in[j] = l_crpix1
            crpix2_in[j] = l_crpix2

            if (crpix1_in[j] != 0 && crpix2_in[j] != 0) {
                # Determine the reference extension
                if (crpix1_in[j] > 0 && \
                    (crpix1_in[j] + datax1[j]) < datax2[j]) {
                    if (l_fl_app_rois && ref_loop_max > 1) {
                        refextn[curr_roi] = j
                    } else {
                        refextn[1] = j
                    }
                }
            }

            # initialize ampsorder to current non sorted amps order
            ampsorder[j] = j

            # This is used to store the extension index / version of the
            # extension that ampsorder is refering to
            in_order[j] = kk

            # Read the amplifier's name
            keypar (insciext, "AMPNAME", silent+)
            if (keypar.found) {
                amplifiername[j] = keypar.value
            } else {
                amplifiername[j] = "EXTN NUM "//imextn
            }

            # Store the roi
            extn_roi[j] = curr_roi

            # reset inputs
            imextn = INDEF
            l_datasec = INDEF
            l_detsec = INDEF
            l_ccdsec = INDEF
            curr_roi = INDEF
            l_crpix1 = INDEF
            l_crpix2 = INDEF

        } # End of loop over extensions obtaining information

        scanfileII = ""
        delete (tmpgretroi, verify-, >& "dev$null")

        # Reset total_num_sciextn
        total_num_sciext = j

        if (debug) {
            printlog ("____New total_num_sciext: "//total_num_sciext, \
                l_logfile, verbose+)
        }

        if (total_num_sciext == 0) {
            printlog ("ERROR - GTILE: There are no science extensions that \
                lie within the requested ROIs. Exiting.", \
                l_logfile, verbose+)
            goto crash
        }

        # Get the amplifiers into the correct order
        # Ordered by DETX1
        for (ii = 1; ii < total_num_sciext; ii += 1) {
            jj = ampsorder[ii]
            kk = ampsorder[ii+1]
            if (detx1[jj] > detx1[kk]) {
                ampsorder[ii] = kk
                ampsorder[ii+1] = jj
                if (ii != total_num_sciext) {
                    ii = 0
                }
            }
        }

        # Need to sort by ROIs too if appending ROIs
        if (l_fl_app_rois) {
            for (ii = 1; ii < total_num_sciext; ii += 1) {
                jj = ampsorder[ii]
                kk = ampsorder[ii+1]
                if (extn_roi[jj] > extn_roi[kk]) {
                    ampsorder[ii] = kk
                    ampsorder[ii+1] = jj
                    if (ii != total_num_sciext) {
                        ii = 0
                    }
                }
            }
        }

        # Need to know the number of extensions per ROI
        if (ref_loop_max == 1) {
            curr_roi = 1
            num_ext_roi[curr_roi] = total_num_sciext
        } else {
            for (ii = 1; ii <= total_num_sciext; ii += 1) {
                jj = ampsorder[ii]
                curr_roi = extn_roi[jj]
                num_ext_roi[curr_roi] = num_ext_roi[curr_roi] + 1
            }
        }

        if (debug) {
            printlog ("____ampsorder:", l_logfile, verbose+)
            for (j = 1; j <= total_num_sciext; j = j+1) {
                printlog ("____"//j//": "//ampsorder[j]//"; ROI: "//\
                    extn_roi[ampsorder[j]], \
                    l_logfile, verbose+)
            }
        }

        # Find number of CCDs used by looping of the amps and
        # comparing detsec and ccdsec also find number_of_amps_per_ccd
        # Inititate them - there will always be 1 ccd and 1 amp at least
        # by this point
        curr_roi = 1
        num_ccds[curr_roi] = 1
        number_of_ampsperccd[num_ccds[curr_roi],curr_roi] = 1

        # Check where the first extension lies if ndetroi > 1 and not
        # appending ROIs and not stats only
        if (!l_fl_app_rois && !l_fl_stats_only && \
            ndetroi > 1) {
            for (j = 2; j <= num_dets; j += 1) {
                if (detx1[ampsorder[1]] > det_start[j]) {
                    number_of_ampsperccd[num_ccds[curr_roi],curr_roi] = 0
                    num_ccds[curr_roi] = num_ccds[curr_roi] + 1
                    number_of_ampsperccd[num_ccds[curr_roi],curr_roi] = 1
                }
            }
        }

        for (j=2; j<=total_num_sciext; j+=1) {

            if (l_fl_app_rois) {
                chk_roi = yes
                # Reset curr_roi and j for next ROI
                while (chk_roi) {

                    if (extn_roi[ampsorder[j]] != extn_roi[ampsorder[j-1]]) {
                        if (debug) {
                            printlog ("____Different ROI: "//\
                                extn_roi[ampsorder[j]], l_logfile, verbose+)
                        }

                        curr_roi += 1
                        num_ccds[curr_roi] = 1
                        number_of_ampsperccd[num_ccds[curr_roi],curr_roi] = 1

                        j += 1

                        if (j > total_num_sciext) {
                            goto ROI_CHECKS_END
                        }
                    } else {
                        chk_roi = no
                    }
                }
            }

            ccddist = ccdx1[ampsorder[j]] - ccdx1[ampsorder[j-1]]
            detdist = detx1[ampsorder[j]] - detx1[ampsorder[j-1]]

            if (debug) {
                printlog ("____detdist: "//detdist//" ccddist: "//\
                    ccddist, l_logfile, verbose+)
            }

            if (detx1[ampsorder[j]] > detx1[ampsorder[j-1]] && \
                (ccdx1[ampsorder[j]] <= ccdx1[ampsorder[j-1]] || \
                (detdist > ccddist && ccddist >= 0))) {

                if (debug) {
                    printlog ("____new_ccd at j: "//j, \
                        l_logfile, verbose+)
                }

                num_ccds[curr_roi] = num_ccds[curr_roi] + 1

                number_of_ampsperccd[num_ccds[curr_roi],curr_roi] = 1

                # Need to determine if a CCD is missing... when not appending
                # ROIS and not fl_stats_only
                if (!l_fl_app_rois && !l_fl_stats_only) {

                    if (detdist > \
                        (det_stop[num_ccds[curr_roi]] - \
                        (det_start[num_ccds[curr_roi]] - 1)) && \
                        (detx1[ampsorder[j]] > \
                        det_stop[num_ccds[curr_roi]])) {

                        if (debug) {
                            printlog ("____extra CCD at j: "//j, \
                                l_logfile, verbose+)
                        }

                        number_of_ampsperccd[num_ccds[curr_roi],curr_roi] = 0
                        num_ccds[curr_roi] = num_ccds[curr_roi] + 1
                        number_of_ampsperccd[num_ccds[curr_roi],curr_roi] = 1
                    }
                }

            } else {

                number_of_ampsperccd[num_ccds[curr_roi],curr_roi] = \
                    number_of_ampsperccd[num_ccds[curr_roi],curr_roi] + 1
            }

ROI_CHECKS_END:
        }

        # Chcek if the number of CCDs is not equal to the num_dets in certain
        # circumstances - basically adding a CCD
        curr_roi = 1
        if (!l_fl_app_rois && !l_fl_stats_only && ndetroi > 1 && \
            num_ccds[curr_roi] > 1 && num_ccds[curr_roi] != num_dets) {
            for (j = (num_ccds[curr_roi] + 1); j <= num_dets; j += 1) {
                num_ccds[curr_roi] = num_ccds[curr_roi] + 1
                number_of_ampsperccd[num_ccds[curr_roi],curr_roi] = 0
            }
        }

        if (debug) {
            printlog ("____Number of ccds is: "//num_ccds, logfile, verbose+)
            printlog ("____num_extensions_per_ccd:", l_logfile, verbose+)

            for (ii = 1; ii <= ndetroi; ii += 1) {
                for (j = 1; j <= num_ccds[ii]; j = j+1) {
                    printlog ("____Number of ccds for ROI "//ii//" is: "//\
                        num_ccds[ii], logfile, verbose+)
                    printlog ("____num_extensions_per_ccd:", \
                        l_logfile, verbose+)
                    printlog ("____    "//j//": "//\
                        number_of_ampsperccd[j,ii], \
                        l_logfile, verbose+)
                }
            }
        }

        num_ext_last = 0
        # Reset REFEXTN if it is equal to 0 - set it to the LHS.
        for (curr_roi = 1; curr_roi <= ref_loop_max; curr_roi += 1) {

            if (refextn[curr_roi] == 0) {
                if (debug) {
                    printlog ("____refextn["//curr_roi//"]: "//\
                        refextn[curr_roi], l_logfile, verbose+)
                }
                refextn[curr_roi] = ampsorder[num_ext_last + 1]
            }

            num_ext_last += num_ext_roi[curr_roi]

            if (debug) {
                printlog ("____refextn["//curr_roi//"]: "//\
                    refextn[curr_roi], l_logfile, verbose+)
            }
        }

#########
#        ### Left this here in case the usage is to be changed
#        # If there is only amp per CCD, not tiling detector and not stats_only
#        # finish with this image
#        if (!l_fl_tile_det && !l_fl_stats_only && \
#            total_num_sciext == num_ccds && !l_fl_stats_only) {
#            printlog ("GTILE WARNING - Only one amplifier per CCD exiting")
#            goto crash
#        }
#########

        # Find the position of the reference extension in ampsorder
        for (curr_roi = 1; curr_roi <= ref_loop_max; curr_roi += 1) {
            ordered_refextn[curr_roi] = 1
            for (j=1; j<=total_num_sciext; j+=1) {
                if (ampsorder[j] == refextn[curr_roi]) {
                    ordered_refextn[curr_roi] = j
                }
            }

            if (debug) {
                printlog ("____The ordered reference extension is: "//\
                    ampsorder[ordered_refextn[curr_roi]]//\
                    " for ROI "//
                    extn_roi[ampsorder[ordered_refextn[curr_roi]]], \
                    logfile, verbose+)
            }
        }

        if (debug) {
            for (jj = 1 ; jj <= total_num_sciext; jj += 1) {
                printlog ("AMP: "//jj//"---------------", logfile, verbose+)
                printlog ("CRPIX1 is: "//\
                    crpix1_in[ampsorder[jj]], logfile, verbose+)
                printlog ("CRPIX2 is: "//\
                    crpix2_in[ampsorder[jj]], logfile, verbose+)
                printlog ("datasec is: "//\
                    datasec[ampsorder[jj]], logfile, verbose+)
                printlog ("ccdsec is: "//\
                    ccdsec[ampsorder[jj]], logfile, verbose+)
                printlog ("detsec is: "//\
                    detsec[ampsorder[jj]], logfile, verbose+)
            }
        }

        # Copy the phu to the outfile
        if (!l_fl_stats_only) {
            imcopy (inimage_files[i]//"[0]", outimg_files[i], verbose=debug)
            if (mdf_exists) {
                tcopy (inimage_files[i]//"["//l_mdf_ext//"]", \
                    outimg_files[i]//"[append,"//l_mdf_ext//"]", \
                    verbose=l_fl_verbose)
            }
        }

        # Check input file against request return ccd
        if (num_ccds[1] == 1 && l_out_ccds != "all" && !l_fl_app_rois) {
            if ((detx1[ampsorder[1]] < det_start[int(l_out_ccds)] || \
                detx1[ampsorder[1]] > det_stop[int(l_out_ccds)]) &&
                (detx2[ampsorder[number_of_ampsperccd[1,1]]] > \
                det_stop[int(l_out_ccds)] || \
                detx2[ampsorder[number_of_ampsperccd[1,1]]] < \
                det_start[int(l_out_ccds)])) {

                printlog ("\nGTILE - WARNING: CCD "//l_out_ccds//\
                    " requested but not no data with"//\
                    " any information from that CCD in it has been provided",\
                    logfile, verbose+)
                printlog ("                 Continuing anyway...", logfile, \
                    verbose+)
            }
        }

        # Tile the amplifiers into their respective CCDs or total detector
        scanfile = extn_list
        in_extn = ""
        dotile = yes
        ref_ext_used = no

        # Initialise counters
        num_extns = 0
        num_sci_extns = 0

        # If an MDF is present and not producing a statistics image
        # Update num_extns by 1
        if (mdf_exists && !l_fl_stats_only) {
            num_extns += 1
        }

        # Counter for input extension names
        ii = 0

        # Loop over the extensions to tile
        while (fscan(scanfile,in_extn) != EOF) {

            ii += 1
            fake_ccd_created = no

            if (debug) {
                printlog ("Working on extensions called: \""//in_extn//"\"", \
                    logfile, verbose+)
            }

            # Initiate num_of_this_extension
            num_of_this_extension = 0

            # Initiate jstart for use below in detemining which amps to use
            jstart = 1

            # Initiate out_ampsname
            out_ampsname = ""
            image_width = 0
            curr_roi = 1

            # loop over ccds
            for (j = 1; j <= num_ccds[curr_roi]; j += 1) {

                if (debug) {
                    printlog ("()()Number of amplifiers for index "//j//" "//\
                        "of j (num_ccds): "//num_ccds[curr_roi]//" "//\
                        number_of_ampsperccd[j,curr_roi]//\
                        " "//curr_roi, logfile, verbose+)
                }

                # Set output extension
                extn = ""
                outnum = INDEF

                # Reset output extension if not stats_only
                if (!l_fl_stats_only && !l_fl_tile_det) {
                    outnum = j
                } else if (l_fl_tile_det && l_fl_app_rois) {
                    outnum = curr_roi
                } else if (l_fl_tile_det && !l_fl_stats_only) {
                    outnum = 1
                }

                if (!isindef(outnum)) {
                    extn = "["//in_extn//outnum//",append]"
                } else {
                    # Stats only!
                    extn = ""
                }

                # Determine the amplifiers to use
                # Don't use the line below it assumes the same number of amps
                # per CCD
                # jjstart = j + ((j-1)*(number_of_ampsperccd[j]-1))
                # Use this instead (jstart gets updated below)
                jjstart = jstart
                jjstop = jjstart + \
                    number_of_ampsperccd[j,curr_roi]
                if (debug) {
                    printlog ("___jjstart: "//jjstart//" "//\
                        "jjstop: "//jjstop, l_logfile, verbose+)
                }

                # Determine output pixtype
                pix_type_ext = ampsorder[jjstart]

                # Have to adjust the extnsion if fake CCD
                if (number_of_ampsperccd[j,curr_roi] == 0) {
                    type_ext = max((jjstart - 1),1)
                    if (j != 1) {
                        pix_type_ext = ampsorder[type_ext]
                    }
                }

                hselect (inimage_files[i]//"["//in_extn//\
                    pix_type_ext//"]", \
                    "i_pixtype", "yes") | scan (in_pixtype)

                out_pixtype = substr(l_pixtype[in_pixtype],1,1)

                # Setup reference WCS extension
                wcs_ref = ampsorder[ordered_refextn[curr_roi]]

                if (l_fl_tile_det) {

                    # Initiate the input list for the task that will combine
                    # the amplifiers
                    if (j == 1) {
                        ref_ext_used = no
                        imtile_list = mktemp("tmpimtile_list")

                        # Set output crpix values
                        crpix1_out = crpix1_in[wcs_ref]
                        crpix2_out = crpix2_in[wcs_ref]
                    } else {

                        crpix1_out += image_width

                        # If the ref extension doesn't lie in CCD2 don't add
                        # chip gap
                        if (j < num_ccds[curr_roi] && l_fl_pad) {
                            if (!ref_ext_used && j == 2) {
                                crpix1_out += gapwidth
                            }
                        }
                    }
                } else {

                    # Initiate the imtile input list
                    imtile_list = mktemp("tmpimtile_list")

                    if (number_of_ampsperccd[j,curr_roi] != 0) {
                        # Set output crpix values for this ccd - use the far
                        # left one
                        crpix1_out = crpix1_in[ampsorder[jjstart]]
                        crpix2_out = crpix2_in[ampsorder[jjstart]]
                    }
                }

                # If tiling detector and middle CCD is missing create tmp file
                # and write it to tile list - should only happen for CCD2 - MS

                if (number_of_ampsperccd[j,curr_roi] == 0 && !l_fl_app_rois) {

                    l_ncols = ((det_stop[j] - (det_start[j] - 1)) / xbin)

                    if (!fake_ccd_created) {
                        tmpccd_file[ii] = mktemp ("tmpccd_file")//".fits"

                        mkpattern (tmpccd_file[ii], pattern="constant", \
                            option="replace", v1=fake_values[ii], v2=0, \
                            pixtype=out_pixtype, ncols=l_ncols, \
                            nlines=image_height, \
                            header=inimage_files[i]//"["//in_extn//\
                            in_order[wcs_ref]//"]")

                        fake_ccd_created = yes
                    }

                    if (!ref_ext_used) {
                        image_width += l_ncols
                    } else {
                        image_width = 0
                    }

                    tmpdetsec = "["//\
                        det_start[j]//":"//\
                        det_stop[j]//","//\
                        "1:"//(image_height * ybin)//"]"

                    tmpccdsec = "[1:"//(l_ncols * xbin)//\
                            ",1:"//(image_height * ybin)//"]"

                    tmpdatasec = "[1:"//l_ncols//",1:"//image_height//"]"

                    # Update the headers with various sections
                    print (default_gemh_pars, >> tmpheader)

                    print (l_key_detsec//" \""//tmpdetsec//\
                        "\" \"Detector section(s)\"", >> tmpheader)

                    print (l_key_ccdsec//" \""//tmpccdsec//\
                        "\" \"CCD section(s)\"", >> tmpheader)

                    print (l_key_datasec//" \""//tmpdatasec//\
                        "\" \"Data section(s)\"", >> tmpheader)

                    # I think this is a hack - MS
                    if (!l_fl_tile_det) {

                        # Assume the centre CCD is the reference
                        if (j == nint(num_ccds[curr_roi]/2.)) {

                            crpix1_out = (l_ncols * 0.5)

                        } else if (real(j) > (num_ccds[curr_roi]/2.)) {

                            crpix1_out = 0 - (l_ncols * ((j - 1 - \
                                nint((num_ccds[curr_roi]/2.))) + 0.5))

                        } else if (real(j) < (num_ccds[curr_roi]/2.)) {

                            crpix1_out = (l_ncols * (j + 0.5))
                        }

                        crpix2_out = (image_height/2)

                        printf ("CRPIX1 %.11f \"Ref pix of axis 1\"\n", \
                            crpix1_out, >> tmpheader)

                        printf ("CRPIX2 %.11f \"Ref pix of axis 2\"\n", \
                            crpix2_out, >> tmpheader)

                    }

                    gemhedit (tmpccd_file[ii], "", "", "", \
                        delete-, upfile=tmpheader)

                    delete (tmpheader, verify-, > "dev$null")

                    print (tmpccd_file[ii]//"[0]", >> imtile_list)

                    if (debug) {
                        printlog ("____Adding fake CCD", \
                            l_logfile, verbose+)
                    }

                    jjstop = 0
                    # Have to add pad parameter here as it may get missed \
                    # later on

                    if (l_fl_pad && j < num_ccds[curr_roi]) {
                        tmpsec = "[*,1:"//image_height//"]"
                        print (name_chipgap[ii]//"[0]"//tmpsec, >> imtile_list)
                    }

                    # Then skip
                    goto SKIP_CHECKS

                } else if (ndetroi > 1 && !l_fl_app_rois) {
                    # Set up the amplifiers to write all ROIs too

                    l_ncols = ((det_stop[j] - (det_start[j] - 1)) / xbin)

                    # Only tile this CCDs worth of amplifiers if
                    # requested - this is duplicated below for the
                    # ndetroi == 1 case
                    if (num_ccds[curr_roi] > 1 && l_fl_stats_only && \
                        l_out_ccds != "all") {
                        # Need to check boundaries and set a flag
                        print (det_start[j]//" "//det_start[int(l_out_ccds)])

                        if ((det_start[j] >= det_start[int(l_out_ccds)] &&\
                            det_start[j] <= det_stop[int(l_out_ccds)]) &&
                            (det_stop[j] <= \
                            det_stop[int(l_out_ccds)] && \
                            det_stop[j] >= \
                            det_start[int(l_out_ccds)])) {

                            tile_amps = yes
                        } else {
                            tile_amps = no
                            goto SKIP_AMPCREATION
                        }
                    } else {
                        tile_amps = yes
                    }

                    for (kk = 1; kk <= namps; kk += 1) {

                        tmpamp_file = "tmpamp_"//j//"_"//kk//"_"//tmpamp_stem

                        # Need to keep all of the pertinant header info so use
                        # the header from jjstart extension

                        fake_width = ((det_stop[j] - (det_start[j] - 1)) / \
                            (namps * xbin))

                        mkpattern (tmpamp_file, pattern="constant", \
                            option="replace", v1=fake_values[ii], v2=0, \
                            pixtype="real", \
                            ncols=fake_width, nlines=image_height, \
                            header=inimage_files[i]//"["//in_extn//\
                            ampsorder[jjstart]//"]")

                        ampwidth = det_stop[j] - (det_start[j] - 1) / namps

                        tmpdetsec = "["//\
                            (det_start[j] + ((kk - 1) * ampwidth))//":"//\
                            (det_start[j] + (kk * ampwidth) - 1)//","//\
                            dety_start//":"//(image_height * ybin)//"]"

                        tmpccdsec = "["//\
                            (1 + ((kk - 1) * ampwidth))//":"//\
                            (kk * ampwidth)//","//\
                            dety_start//":"//(image_height * ybin)//"]"

                        tmpdatasec = "[1:"//l_ncols//",1:"//image_height//"]"

                        # Update the headers with various sections

                        # Update the headers with various sections
                        print (default_gemh_pars, >> tmpheader)

                        print (l_key_detsec//" \""//tmpdetsec//\
                            "\" \"Detector section(s)\"", >> tmpheader)

                        print (l_key_ccdsec//" \""//tmpccdsec//\
                            "\" \"CCD section(s)\"", >> tmpheader)

                        print (l_key_datasec//" \""//tmpdatasec//\
                            "\" \"Data section(s)\"", >> tmpheader)

                        gemhedit (tmpamp_file, "", "", "", \
                            delete-, upfile=tmpheader)

                        delete (tmpheader, verify-, > "dev$null")

                        print (tmpamp_file//"[0]", >> imtile_list)
                    }
SKIP_AMPCREATION:
                }

                if (debug) {
                    printlog ("___outpixtype: "//out_pixtype, \
                        l_logfile, verbose+)
                }

                # Reset image_width
                image_width = 0

                # Loop over amps for this ccd and create input for the
                # combining task
                for (jj=jjstart; jj<jjstop; jj+=1) {

                    # Obtain the current amplifier
                    k = ampsorder[jj]
                    if (debug) {
                        printlog ("Now looking at amplifiers in order: "//\
                            "amp "//k, logfile, verbose+)
                    }

                    # Set flag for updating CRPIX1
                    if (k == wcs_ref) {
                        ref_ext_used = yes
                    }

                    curr_extn = inimage_files[i]//"["//\
                        in_extn//in_order[k]//"]"//\
                        datasec[k]

                    # If needed paste into amplifiers first
                    if (!l_fl_app_rois && ndetroi > 1) {

                        if (!tile_amps) {

                            dotile = no
                            goto NEXT_EXTENSION
                        } else {

                            dotile = yes
                        }

                        # Paste the extenions into one image
                        for (kk = 1; kk <= namps; kk += 1) {
                            tmpamp_file = "tmpamp_"//j//"_"//kk//"_"//\
                                tmpamp_stem//"[0,overwrite]"

                            if (debug) {
                                printlog ("____j: "//j//" kk: "//kk//\
                                    " "//tmpamp_file, l_logfile, verbose+)
                            }

                            detdist = (det_stop[j] - (det_start[j] - 1))
                            ccddist = detdist / namps

                            if (debug) {
                                printlog ("____ccddist: "//detdist//" "//\
                                    "ampdist: "//ccddist, l_logfile, verbose+)
                                printlog ("____ccdsec limits: "//\
                                    (1 + ((kk - 1) * ccddist))//" "//\
                                    (1 + ((kk) * ccddist)), \
                                     l_logfile, verbose+)
                                printlog ("____ccdx1: "//ccdx1[k]//\
                                    " ccdx2: "//ccdx2[k], l_logfile, verbose+)
                            }

                            fkccdx1 = (1 + ((kk - 1) * ccddist))
                            fkccdx2 = (1 + (kk * ccddist))

                            # No longer worrying about Hamamatsu data
                            fkccdy1 = 1

                            if ((ccdx1[k] >= fkccdx1) && (ccdx2[k] < fkccdx2)){

                                # Create binned version of CCDSEC relative to
                                # the current amplifier.
                                tx1 = (((ccdx1[k] - fkccdx1) / xbin) + 1)
                                tx2 = tx1 + \
                                    (((ccdx2[k] - (ccdx1[k] - 1)) / xbin) - 1)

                                ty1 = (((ccdy1[k] - fkccdy1) / ybin) + 1)
                                ty2 = ty1 + \
                                    (((ccdy2[k] - (ccdy1[k] - 1)) / ybin) - 1)

                                ampsec = "["//tx1//":"//tx2//","//ty1//":"//\
                                    ty2//"]"

                                if (debug) {
                                    printlog ("____ampsec is: "//\
                                        ampsec, logfile, verbose+)
                                }

                                imcopy (curr_extn, tmpamp_file//ampsec, \
                                    verbose=debug)
                                break
                            }
                        }

                        # Sort the CRPIX values out - another hack - MS
                        if (k == wcs_ref && l_fl_tile_det) {
                            # update image_width as is used to update crpix2
                            # when tiling
                            image_width += tx1 + \
                                (((kk - 1) + ((j - 1) * namps)) * fake_width)

                            crpix2_out += ty1

                        } else if (jj == jjstart && !l_fl_tile_det) {

                            crpix2_out += ty1

                            # Have to update crpix1 differently depending on
                            # it's sign
                            crpix1_out += tx1

                            # Need to add any additional amplifier widths if
                            # they are before this extension
                            if (kk != 1) {
                                if (j == nint(num_ccds[curr_roi]/2.)) {
                                    crpix1_out += ((kk - 1) * fake_width)

                                } else if (real(j) > (num_ccds[curr_roi]/2.)) {

                                    crpix1_out += (((kk - 1) + (j - \
                                        1 - nint(num_ccds[curr_roi]/2.)) * \
                                        namps) * fake_width)

                                } else if (real(j) < (num_ccds[curr_roi]/2.)) {
                                    crpix1_out += (((kk - 1) + \
                                        ((j - 1) * namps)) * fake_width)
                                }
                            }
                        }

                        goto NEXT_EXTENSION
                    } else {
                        if (!ref_ext_used) {
                            image_width += data_width[k]
                        }
                        image_height = data_height[k]
                    }

                    # Setup the input list for imtile

                    # Only tile this CCDs worth of amplifiers if requested
                    if (num_ccds[curr_roi] > 1 && l_fl_stats_only && \
                        l_out_ccds != "all") {

                        # reset dotile flag
                        dotile = no

                        # Need to check boundaries and set a flag
                        if ((detx1[k] >= det_start[int(l_out_ccds)] && \
                            detx1[k] <= det_stop[int(l_out_ccds)]) &&
                            (detx2[k] <= \
                            det_stop[int(l_out_ccds)] && \
                            detx2[k] >= \
                            det_start[int(l_out_ccds)])) {

                            # Printlog file name to imtile input list
                            print (curr_extn, >> imtile_list)
                            if (jj == jjstop - 1) {
                                # Reset the do tile flag
                                dotile = yes
                            }
                       }

                    } else { # Else do tiling regardless
                        print (curr_extn, >> imtile_list)
                    }

NEXT_EXTENSION:
                    # Need to reset the section used of the chip gap file
                    if (l_fl_pad && jj == (jjstop - 1) && \
                        j < num_ccds[curr_roi]) {

                        tmpsec = "[*,1:"//image_height//"]"

                        print (name_chipgap[ii]//"[0]"//tmpsec, >> imtile_list)
                    }

                    # Update the output amp name list only if it is to be used
                    if (dotile) {
                        out_ampsname = out_ampsname//","//amplifiername[k]
                    }

                    if (debug) {
                        printlog ("Current out_ampsname is: "//out_ampsname, \
                            logfile, verbose+)
                    }

                } # End of configuring paste variable (for loop)

SKIP_CHECKS:
                # Do the tiling!
                if ((!l_fl_tile_det && dotile) \
                    || (l_fl_tile_det && j==num_ccds[curr_roi])) {

                    if (debug) {
                        printlog ("In first section of the dotile loop", \
                            logfile, verbose+)
                    }

                    # Only call imjoin if more than one amp - else imcopy
                    if (!l_fl_tile_det && \
                        number_of_ampsperccd[j,curr_roi] == 1 && ndetroi == 1){
                        if (debug) {
                            printlog ("Only imcopying the file", logfile, \
                                verbose+)
                            printlog ("Output name is: "//outimg_files[i]//\
                                extn, logfile, verbose+)

                        }
                        scanfileII = imtile_list
                        while (fscan(scanfileII, in_tile_file) != EOF) {
                            imcopy (in_tile_file, outimg_files[i]//extn, \
                                verbose=debug)
                            # Update the extension version - only needed
                            # for raw frames when not producing a statistics
                            # image
                            num_of_this_extension += 1
                            if (in_extn == "") {
                                if (!l_fl_stats_only) {

                                    print (default_gemh_pars, >> tmpheader)

                                    print ("EXTVER "//num_of_this_extension, \
                                        >> tmpheader)

                                    # Have to update EXTNAME too to get arround
                                    # Pyraf not liking the imcopy /append
                                    # of raw data frames
                                    print ("EXTNAME \"TEST\"", \
                                        >> tmpheader)

                                    gemhedit (outimg_files[i]//extn, \
                                        "", "", "", \
                                        delete-, upfile=tmpheader)

                                    delete (tmpheader, verify-, > "dev$null")

                                    if (debug) {
                                        printlog ("Updated EXTVER to: "//\
                                            num_of_this_extension, \
                                            l_logfile, verbose+)
                                    }
                                }

                            }

                            # Slight code duplication to update DATASEC - MS

                            # Obtain the keyword value for datasec
                            keypar (outimg_files[i]//extn, l_key_datasec, \
                                silent+)
                            if (keypar.found) {
                                l_sec = keypar.value
                                junk = fscanf (l_sec, "[%d:%d,%d:%d]", tx1,\
                                    tx2, ty1, ty2)

                                # Extract the individual values
                                secx1 = tx1
                                secx2 = tx2
                                secy1 = ty1
                                secy2 = ty2

                                # Obtain the NAXIS1 value for this extension to
                                # update DATASEC, DETSEC and CCDSEC
                                keypar (outimg_files[i]//extn, "i_naxis1", \
                                    silent+)
                                if (keypar.found) {
                                    image_width = int(keypar.value)
                                } else {
                                    image_width = 0
                                }

                                # Reset to the new DATASEC
                                secx1 = 1
                                secx2 = image_width

                                # Define the section value
                                l_sec = "["//secx1//\
                                    ":"//secx2//","//\
                                    secy1//":"//secy2//"]"

                                # Update the keyword
                                print (default_gemh_pars, >> tmpheader)

                                print (l_key_datasec//" \""//l_sec//\
                                    "\" \"Data section(s)\"", >> tmpheader)

                                gemhedit (outimg_files[i]//extn, \
                                    "", "", "", \
                                    delete-, upfile=tmpheader)

                                delete (tmpheader, verify-, > "dev$null")

                            } else {
                                printlog ("WARNING - GTILE: "//l_key_datasec//\
                                    " not found", l_logfile, verbose+)
                                printlog ("                 Ignoring and "//\
                                    "going to next keyword", l_logfile, \
                                    verbose+)
                            }
                        }
                        scanfileII = ""
                    } else {

                        if (debug) {
                            type (imtile_list)
                            printlog ("Writing to "//outimg_files[i]//extn, \
                                logfile, verbose+)
                        }

                        # Doing the joining
                        imjoin (input="@"//imtile_list, \
                            output=outimg_files[i]//extn, join_dim=1,\
                            pixtype=out_pixtype, verbose=yes, >> l_logfile)

                        # Update the extension version - only needed
                        # for raw frames when not producing a statistics image
                        num_of_this_extension += 1
                        if (in_extn == "") {
                            if (!l_fl_stats_only) {

                               # Update the keyword
                                print (default_gemh_pars, >> tmpheader)

                                print ("EXTVER "//num_of_this_extension, \
                                    >> tmpheader)

                                gemhedit (outimg_files[i]//extn, \
                                    "", "", "", \
                                    delete-, upfile=tmpheader)

                                delete (tmpheader, verify-, > "dev$null")

                                if (debug) {
                                    printlog ("Updated EXTVER to: "//\
                                        num_of_this_extension, \
                                            l_logfile, verbose+)
                                }
                            }
                        }

                        # Update the amplifier name in output extension
                        if (substr(out_ampsname,1,1) == ",") {
                            out_ampsname=substr(out_ampsname,2, \
                                strlen(out_ampsname))
                        }

                        # Update AMPNAME, CRPIX, DATASEC, DETSEC and CCDSEC
                        # keywords here

                        # AMPNAME
                        # Only do this if not DQ plane
                        if (in_extn != l_dq_ext//",") {

                            # Update the keyword
                            print (default_gemh_pars, >> tmpheader)

                            # AMPNAME
                            print ("AMPNAME \""//out_ampsname//"\" \".\"", \
                                >> tmpheader)

                            # CRPIX1
                            printf ("CRPIX1 %.11f \".\"\n", crpix1_out,\
                                >> tmpheader)

                            # CRPIX2
                            printf ("CRPIX2 %.11f \".\"\n", crpix2_out,\
                                >> tmpheader)

                            gemhedit (outimg_files[i]//extn, \
                                "", "", "", \
                                delete-, upfile=tmpheader)

                            delete (tmpheader, verify-, > "dev$null")
                        }

                        # Obtain the NAXIS(1/2) value for this extension to
                        # update DATASEC, DETSEC and CCDSEC
                        keypar (outimg_files[i]//extn, "i_naxis1", silent+)
                        if (keypar.found) {
                            image_width = int(keypar.value)
                        } else {
                            image_width = 0
                        }

                        keypar (outimg_files[i]//extn, "i_naxis2", silent+)
                        if (keypar.found) {
                            image_height = int(keypar.value)
                        } else {
                            image_height = 0
                        }

                        # imjoin uses the header information from the 1st file
                        # it receives in that call as the header for the output

                        # Update the keyword
                        print (default_gemh_pars, >> tmpheader)
                        one_key_found = no

                        # Loop over array describing header keywords
                        for (kk=1; kk<=3; kk+=1) {

                            # Set the header keyword to read
                            l_key = array_keys[kk]

                            # Obtain the keyword value
                            keypar (outimg_files[i]//extn, l_key, silent+)
                            if (keypar.found) {
                                one_key_found = yes

                                l_sec = keypar.value
                                junk = fscanf (l_sec, "[%d:%d,%d:%d]", tx1,\
                                    tx2, ty1, ty2)

                                # Extract the individual values
                                secx1 = tx1
                                secx2 = tx2
                                secy1 = ty1
                                secy2 = ty2

                                # Update the key value width and height - even
                                # if they are not required to be updated;
                                # covers all cases then
                                if (l_key != l_key_datasec) {
                                    # Override secy1 value as they are the for
                                    # DETSEC and CCDSEC - This is a bit of a
                                    # hack for Hamamaatsu CCDs - MS
                                    if ((ndetroi > 1 && !l_fl_app_rois) || \
                                        (image_height == height_iccd)) {
                                        secy1 = dety_start
                                    }

                                    secx2 = (secx1 - 1) + (image_width * xbin)
                                    secy2 = (secy1 - 1) + (image_height * ybin)
                                } else {
                                    secx1 = 1
                                    secx2 = image_width
                                    secy1 = 1
                                    secy2 = image_height
                                }

                                # Define the section value
                                l_sec = "["//secx1//\
                                    ":"//secx2//","//\
                                    secy1//":"//secy2//"]"

                                # Update the keyword
                                print (l_key//" \""//l_sec//"\" \".\"", \
                                    >> tmpheader)

                            } else {
                                printlog ("WARNING - GTILE: "//l_key//" not"//\
                                    " found", l_logfile, verbose+)
                                printlog ("                 Ignoring and "//\
                                    "going to next keyword", \
                                    l_logfile, verbose+)
                            }
                        } # End of for loop updating keywords

                        if (one_key_found) {
                            gemhedit (outimg_files[i]//extn, \
                                "", "", "", \
                                delete-, upfile=tmpheader)
                        }
                        delete (tmpheader, verify-, > "dev$null")
                    } # End of imjoin / imcopy calls

                    if (l_fl_tile_det && dotile) {
                        # Read the cd?_? keywords from the wcs_ref extension
			l_cd1_1 = 1.
			l_cd1_2 = 1.
			l_cd2_1 = 1.
			l_cd2_2 = 1.

                        hselect (inimage_files[i]//"["//\
                            in_extn//in_order[wcs_ref]//"]", \
                            "CD1_1,CD1_2,CD2_1,CD2_2", yes) | \
                            scan (l_cd1_1, l_cd1_2, l_cd2_1, l_cd2_2)

                        print (default_gemh_pars, >> tmpheader)
                        print ("CD1_1 "//l_cd1_1//" \".\"", >> tmpheader)
                        print ("CD1_2 "//l_cd1_2//" \".\"", >> tmpheader)
                        print ("CD2_1 "//l_cd2_1//" \".\"", >> tmpheader)
                        print ("CD2_2 "//l_cd2_2//" \".\"", >> tmpheader)
                        gemhedit (outimg_files[i]//extn, "", "", "", \
                                delete-, upfile=tmpheader)
                        delete (tmpheader, verify-, > "dev$null")
                    }

                    # Clean header keywords from GRETROI if required
                    if (clean_header) {
                        print (default_gemh_pars, >> tmpheader)
                        print (l_key_extnroi//" delete+", >> tmpheader)
                        print (roi_key_detsec//" delete+", >> tmpheader)
                        print (roi_key_ccdsec//" delete+", >> tmpheader)
                        print (roi_key_datasec//" delete+", >> tmpheader)
                        gemhedit (outimg_files[i]//extn, "", "", "", \
                                delete-, upfile=tmpheader)
                        delete (tmpheader, verify-, > "dev$null")
                    }

                    # Reset out_ampsname
                    out_ampsname = ""

                    # Reset the LTV1 and 2 keywords using wcsreset
                    # wcs=physicsal
                    # This is so logical and physical coordinates relate to
                    # the new file not any images they may have come from
                    wcsreset (outimg_files[i]//extn, wcs="physical",\
                        verbose-)

                    # Count the total number of science extensions
                    if (in_extn == l_sciext) {
                        num_sci_extns += 1
                    }

                    # Count the total number of extensions in the output
                    num_extns += 1

                } # End of tiling

                # Update jstart to be correct and not expect the same number
                # of amps per ccd
                jstart += number_of_ampsperccd[j,curr_roi]

                # If appending ROIs reset the curr_roi and j
                if (j == num_ccds[curr_roi] && l_fl_app_rois) {

                    if (curr_roi != ndetroi) {
                        j = 0
                        curr_roi += 1
                    }
                }
            } # End of pasting and looping over ccds

            # Update NEXTEND and NSCIETN here
            print (default_gemh_pars, >> tmpheader)

            # NEXTEND
            print ("NEXTEND "//num_extns//" \"Number of extensions\"", \
                >> tmpheader)

            # NSCIEXT
            print ("NSCIEXT "//num_sci_extns//\
                " \"Number of science extensions\"", \
                >> tmpheader)

            gemhedit (outimg_files[i]//"[0]", \
                "", "", "", \
                delete-, upfile=tmpheader)

            delete (tmpheader, verify-, > "dev$null")

            # Due to PyRAF not liking to append imjoined raw frames to
            # the output, update EXTVER of the output frames to -1 here - MS
            # If input was raw should be no other extensions other than
            # those joined earlier; though will double check for MDFs - MS
            # Only needs to be done for tiled raw frame with non-statistics
            # output
            curr_roi = 1
            if (in_extn == "" && !l_fl_stats_only && !l_fl_tile_det) {

                for (j = 1; j <= num_ccds[curr_roi]; j += 1) {
                    print (default_gemh_pars, >> tmpheader)

                    current_extn = j
                    if (mdf_exists) {
                        current_extn += 1
                    }
                    if (debug) {
                        printlog ("current_extn for "//j//" is: "//\
                            current_extn, l_logfile, verbose+)
                    }

                    # EXTVER
                    print ("EXTVER -1", >> tmpheader)

                    # This test and subsequent update is here for
                    # the fact that Pyraf doesn't like ot append raw data
                    # especially when imcopying
                    keypar (outimg_files[i]//"["//current_extn//"]", \
                        "EXTNAME", silent+)
                    if (keypar.value == "TEST") {
                        # EXTVER
                        print ("EXTNAME \"\"", >> tmpheader)
                    }
                    gemhedit (outimg_files[i]//"["//current_extn//"]", \
                        "", "", "", \
                        delete-, upfile=tmpheader)

                    delete (tmpheader, verify-, > "dev$null")

                    if (j == num_ccds[curr_roi] && l_fl_app_rois) {
                        if (curr_roi != ndetroi) {
                            j = 1
                            curr_roi += 1
                        }
                    }
                } # End of for loop
            }

            # Can only delete these tmpimtile_list files here as they are
            # needed in case of l_fl_stats_only. Saves a number of if loops -
            # MS
            delete ("tmpimtile_list*", verify-, >& "dev$null")

        } # End of loop over extension types

        delete (extn_list, verify-, >& "dev$null")

        if (fake_ccd_created) {
            for (ii = 1; ii <= num_ext_names; ii += 1) {
                imdelete (tmpccd_file[ii], verify-, >& "dev$null")
            }
        }

        if (l_fl_pad) {
            for (ii = 1; ii <= num_ext_names; ii += 1) {
                imdelete (name_chipgap[ii], verify=no, >& "dev$null")
            }
        }

        # Write TRIMMED keyword to PHU
        print (default_gemh_pars, >> tmpheader)
        print ("TRIMMED yes \"Overscan section trimmed\"", \
            >> tmpheader)
        print (l_key_typeret//" \""//type_ret//"\" \"Data returned by \
            GRETROI\"", >> tmpheader)

        # Write updated ROIs to PHU if needed
        if (update_rois) {
            print ("DETNROI 1 \"No. regions of interest after GTILE\"", \
                >> tmpheader)
            print ("DETRO1X "//roidetx1//" \"ROI 1 X start after GTILE\"", \
                >> tmpheader)
            print ("DETRO1XS "//roidetx1s//"\"ROI 1 X size after GTILE\"", \
                >> tmpheader)
            print ("DETRO1Y "//roidety1//" \"ROI 1 Y start after GTILE\"", \
                >> tmpheader)
            print ("DETRO1YS "//roidety1s//" \"ROI 1 Y size after GTILE\"", \
                >> tmpheader)
        }

        # Type of tiling done
        print (l_key_gttype//" \""//l_gttype//"\" \""//l_gtcomm//"\"", \
            >> tmpheader)

        # Write GTILE time stmap to header
        gemdate (zone="UT")

        # GTILE
        print ("GTILE \""//gemdate.outdate//"\" \"UT Time stamp for GTILE\"", \
            >> tmpheader)

        # GEM-TLM
        # Added extra space due to old time stamps being really long!
        print ("\"GEM-TLM\" \""//gemdate.outdate//\
            "\" \"UT Last modification with GEMINI            \"", \
            >> tmpheader)

        gemhedit (outimg_files[i]//"[0]", \
            "", "", "", delete-, upfile=tmpheader)

        delete (tmpheader, verify-, > "dev$null")
    } # End of loop over input

    ####################
    # End of hard work #
    ####################

    goto clean

crash:

    status = 1

clean:

    # Clean up a bit - delete inlist and outlist
    scanfile = ""
    scanfileII = ""
    imdelete ("*"//tmpamp_stem, verify-, >& "dev$null")
    delete (inlist//","//outlist//","//tmplist//","//tmplist2, \
        verify-, >& "dev$null")

    # Printlog finish time
    gemdate (zone="UT")
    printlog ("\nGTILE - Finished: "//gemdate.outdate//"\n", \
        logfile, l_fl_verbose)

    if (status == 0) {
        printlog ("GTILE - Exit status: GOOD\n", logfile, l_fl_verbose)
    } else {
        printlog ("GTILE - Exit status: ERROR\n", logfile, l_fl_verbose)
    }

    if (profile) {
        date ("+%H:%M:%S.%N%n") | scan (sdate)
        printlog ("____GTILE Finished: "//sdate, l_logfile, verbose+)
    }
end

