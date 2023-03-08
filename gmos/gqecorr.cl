# Copyright(c) 2011-2017 Association of Universities for Research in Astronomy, Inc.

procedure gqecorr (inimages)

# This task will correct the QE discrepancies seen between the 3 GMOS CCDs such
# that the QE efficiencies of CCDs 1 and 3 match that of CCD 2. This can be done
# in all spectroscopic observing modes.

# This correction will be made available for all types of GMOS detectors
# though are not be available for earlier GMOS-S detectors 2003-04-01

# As a wavelength solution is required, a reference frame is needed. For a
# GSAPPWAVE solution, only a combined (possibly cut) flat is required to be
# ran through GSAPPWAVE.

# Inputs must be gprepared'd!

# Correction has to be applied to science and flat frames only but has to be
# done immediately after bias/dark correction (+ overscan correction?) before
# any other manipulation.

# Any changes in the naming convention of the correction images / data files
# must be reflected in gireduce

# Overview:
#
#    GQECORR will create QE correction frame(s) (spectroscopy only) from
#    supplied reference image(s) based upon the GSAPPWAVE solution (currently).
#    The QE information required other than a wavelength solution, is read
#    from a file in the gmos$data/ directory (see default parameter value for
#    filename). In addition the user can supply their own QE data in file form
#    with the same format; or the can put it in by hand in the format:
#    qecorr_data = "1.2,0.33,3.5E-3:1:0.89,0.32,1.5e-4"
#    where the information for each CCD is separated by a colon.
#
#    The same applies for imaging data, however, no correction image is
#    created. Instead a '.dat' file containing the CCD names and the
#    reciprocal of the QE factor is written. Currently imaging data is not
#    supported.
#
#    GQECORR can also correct the input data providing they are not GMOSAIC'ed
#    and have been BIAS/DARK subtracted.
#
#    Any existing correction files will be used instead of creating new ones.
#    The user can supply a correction image rather than a reference image if
#    correcting the data or wishing to create a N&S correction image based on
#    that correction image, as when reducing N&S data.

char    inimages    {prompt="GMOS image(s) to QE correct"}
char    outimages   {"", prompt="Output image(s)"}
char    outpref     {"q", prompt="Prefix for output image(s)"}
char    refimages   {"", prompt="Reference image(s) with wavelength solution(s) - spectroscopy only"}
bool    fl_correct  {yes, prompt="Correct the input image(s)?"}
bool    fl_keep     {no, prompt="Keep correction data?"}
char    corrimages  {"", prompt="Name for QE correction image(s)"}
char    corrimpref  {"qecorr", prompt="Prefix for QE correction data"}
int     ifu_preced  {1, min=1, max=2, prompt="IFU slit to take precedence in IFU-2 mode if there is overlap between slits"}
bool    fl_vardq    {yes, prompt="Propagate variance and data quality planes?"}
char    sci_ext     {"SCI", prompt="Name of science extension(s)"}
char    var_ext     {"VAR", prompt="Name of variance extension(s)"}
char    dq_ext      {"DQ", prompt="Name of data quality extension(s)"}
char    mdf_ext     {"MDF", prompt="Mask definition file extension name"}
char    key_detsec  {"DETSEC", prompt="Header keyword for detector section"}
char    key_ccdsec  {"CCDSEC", prompt="Header keyword for CCD section"}
char    key_datasec {"DATASEC", prompt="Header keyword for data section"}
char    key_biassec {"BIASSEC", prompt="Header keyword for overscan section"}
char    key_ccdsum  {"CCDSUM", prompt="Header keyword for CCD binning"}
char    qecorr_data {"gmos$data/gmosQEfactors.dat", prompt="Data file that contains the QE correction information."}
char    database    {"database/", prompt="Location of wavelength solution database"}
char    logfile     {"", prompt="Logfile to use"}
bool    verbose     {no, prompt="Verbose"}
int     status      {0, prompt="Exit status (0=good)"}
struct  *scanfile   {prompt="For internal use only"}
char    scanfile2   {"",prompt="For internal use only"}

# Note scanfile2 is only used in gireduce when gqecorr is called, such that a
# a tmp file can be written to disk so gireduce can easily access the
# correction data files
begin

    ################
    # Set up variables, default values and chache parameters of tasks used
    ################

    # Setup local variables for input parameters
    char    l_inimages = ""
    char    l_outimages = ""
    char    l_outpref = ""
    char    l_refimages = ""
    bool    l_correct_inimages
    bool    l_keep_corrimages
    char    l_corrimages = ""
    char    l_corrimpref = ""
    bool    l_fl_vardq
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_mdf_ext = ""
    char    l_key_detsec = ""
    char    l_key_ccdsec = ""
    char    l_key_datasec = ""
    char    l_key_biassec = ""
    char    l_key_ccdsum = ""
    char    l_qecorr_data = ""
    char    l_database = ""
    char    l_logfile = ""
    bool    l_verbose

    # Set up variables used in code
    char    filename, tmpfilename, inlist, outlist, reflist, corrimlist
    char    inimage_files[200], outimage_files[200], refimage_files[200]
    char    corrimage_files[200], listfile, test_string, errstr
    char    l_refimg, l_input, l_corrimg, l_output
    char    l_sciext, l_varext, l_dqext, extn_list[3]
    char    l_key_gprepare, l_key_gqecorr, l_key_overscan, l_key_biasimg
    char    l_key_darkimg, l_key_obsmode, l_key_gscut
    char    l_key_gfextract, l_key_gsappwav, l_key_gmosaic, l_key_gswave
    char    l_key_trimmed, l_key_maskname, l_key_dettype, l_key_detector
    char    l_key_instrument, l_key_nsciext, l_key_filter1
    char    l_key_filter2, l_key_noddistance, l_key_anodcount, l_key_bnodcount
    char    l_key_crpix1, l_key_refpix, l_key_imgfactor
    char    l_key_crval1, l_key_cd1_1, l_key_ccdname, l_key_mdfrow
    char    l_key_corrim, l_imgfactor, l_key_gqenod, l_key_gemtlm
    char    l_key_refimage, l_key_grating, l_key_grtilt, nodanswer
    char    l_key_date1, l_key_date2, l_task
    char    l_tabkey_slittilt, l_tabkey_slittilt_LS, secondop, inmaskname
    char    indettype, indetector, ininstrument, ingrating, ingrtilt
    char    inobsmode[200], refobsmode[200], thisobsmode
    char    tmpcorrim, tmppixels, tmponlypix, tmppiximg, tmp1dimg, tmp2dimg
    char    tmpoutimg, iminscoordds, tmp2dimgtilt
    char    thisfilter1, thisfilter2, lastfilter1, lastfilter2
    char    thisdettype, thisextn, thisinstrument, thisccd, thisccd_data
    char    thisdetector, db_end
    char    correction_type, tmpcoeffs, tmpimgcorrected, tmpqedata
    char    stored_dettype, currentarea, currentsection, lambda_expr
    char    iminscoords, corrapp, gemzone
    char    imexpr_char, calexpr, lasttmpcorrim, thismode
    char    ccdname[3], inextn, outextn, inccdname, dataccdname, tmpccdname
    char    corrccdname, thisextnappend, insertextn, outextnappend
    char    datasecin[200,12], ccdsecin[12], detsecin[12], imgpixtype
    char    somefilename, nextsomefilename, numstring, img_corrfact[200,3,2]
    char    im_expr_shift, tmpimgshift, tmpimg2shift, indatasec, file_extn
    char    ns_prefix, qefile_to_read, tmpqe2data, factor_description
    char    column_to_read, stopdatekey, obsdate_inkey_value
    char    obsdate_refkey_value, tmpifu, wavextn, wav_extname, cut_wavextn
    int     hh, ii, jj, kk, ll, nn, atlocation, junk, chk_lines
    int     num_inimages, num_outimages, num_refimages, num_corrimages
    int     num_sci, num_ref_sci, count_imgobsmode
    int     xinbin[200], yinbin[200], xrefbin[200], yrefbin[200]
    int     xbin, ybin, fullxdim, fullydim, num_ccds,  num_dettypes
    int     ccdwidth[3], det_info[4,2], chipgap, ampwidth, ampheight
    int     datax1[12], datax2[12], datay1[12], datay2[12], imgshift
    int     detx1, detx2, dety1, dety2, xstart, xstop, ystart, ystop
    int     detinx1[12], detinx2[12], detiny1[12], detiny2[12]
    int     xmax, ymax, thisccdstart, thisccdstop, thisorder
    int     ccdsecx1[12], ccdsecx2[12], ccdsecy1[12], ccdsecy2[12]
    int     ccdx1, ccdx2, ccdy1, ccdy2, fail_counter
    int     tmpint, ampsorder[12], number_of_ampsperccd[3]
    int     dummyx1, dummyx2, dummyy1, dummyy2, total_num_amplifiers
    int     row_to_read, wav_adjustment, ycent_int, user_qenum_count
    int     extns_to_correct, num_corrim_exists, num_ns_corrim_exists
    int     num_ns_infiles, num_corrimgs_to_use, num_incorrimages
    int     inobsdate, refobsdate, gmossccddate, gmossbiasdate, n_extentions
    int     ifuslit_precedence, first_stop, overlap
    real    refpix, refval, refxdelta, refydelta, thiscoeff, tilt
    real    pi, ang_to_use, ycentre, maskpa
    bool    l_fl_useprefix, isimaging, nocorrimages, alldone
    bool    norefimages, fl_filter_match, iscutLS, istilted, userdefined
    bool    userdefined_orig
    bool    corrim_exists[200], refim_gqenod[200], corrim_skip_flag
    struct  currentccdsec, l_struct_description
    char    det_names[4]="SDSU II CCD","SDSU II e2v DD CCD42-90","S10892",\
                         "S10892-N"

    # Hidden parameter to be moved futher down in the future / Setup for users
    bool    fl_fullcorr = no # Full wavelength solution; no = GSAPPWAVE
    bool    l_fl_fullcorr = yes
    bool    fl_relative_corr = yes # Type of QE correction to apply
    bool    l_fl_relative_corr = yes
    bool    fl_untransform_corrim = no # Untransform the data (ungmosaic)
    bool    l_fl_untransform_corrim #reset later on to fl_untransform_corrim
    char    l_key_refimg = "" # For reading refimages from a keyword
    bool    debug = no

    char tmpcoords, tmp2coords, lambda_fname, l_dbfile
    char l_wav_solution, ext_formated

    # Read input parameters
    junk = fscan (inimages, l_inimages)
    junk = fscan (refimages, l_refimages)
    junk = fscan (outimages, l_outimages)
    junk = fscan (outpref, l_outpref)
    l_correct_inimages = fl_correct
    l_keep_corrimages  = fl_keep
    junk = fscan (corrimages, l_corrimages)
    junk = fscan (corrimpref, l_corrimpref)
    ifuslit_precedence = ifu_preced
    l_fl_vardq = fl_vardq
    junk = fscan (sci_ext, l_sci_ext)
    junk = fscan (var_ext, l_var_ext)
    junk = fscan (dq_ext, l_dq_ext)
    junk = fscan (mdf_ext, l_mdf_ext)
    junk = fscan (key_detsec, l_key_detsec)
    junk = fscan (key_ccdsec, l_key_ccdsec)
    junk = fscan (key_datasec, l_key_datasec)
    junk = fscan (key_biassec, l_key_biassec)
    junk = fscan (key_ccdsum, l_key_ccdsum)
    junk = fscan (qecorr_data, l_qecorr_data)
    junk = fscan (database, l_database)
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    # Set local variables here that users are not to set
    # Used thoughout the code

    # Keywords
    l_key_gprepare   =  "GPREPARE" # GPREPARE header keyword - phu
    l_key_gqecorr    =  "GQECORR" # gqecorr timestamp - phu
    l_key_gqenod     =  "GQENOD" # gqecorr adjusted corrim for N&S - phu - bool
    l_key_gemtlm     =  "GEM-TLM" # Last update from gemini task - phu
    l_key_imgfactor  =  "QEFACTOR" # QE correction factor - imaging only
    l_key_corrim     =  "QECORRIM" # GQECORR correction image - phu
    l_key_refimage   =  "GQEREFIM" # GQECORR reference image - phu
    l_key_overscan   =  "OVERSEC" # keyword for overscan correction - sci
    l_key_trimmed    =  "Trimmed" # gireduce trimmed keyword - sci
    l_key_biasimg    =  "BIASIM" # Bias image - phu
    l_key_darkimg    =  "DARKIM" # Dark image - phu
    l_key_obsmode    =  "OBSMODE" # Observation mode keyword - phu
    l_key_gscut      =  "GSCUT" # gsappwave keyword - phu
    l_key_gfextract  =  "GFEXTRACT" # gsappwave keyword - phu
    l_key_gsappwav   =  "GSAPPWAV" # gsappwave keyword - phu
    l_key_gswave     =  "GSWAVELE" # gswavelength keyword - phu
    l_key_gmosaic    =  "GMOSAIC" # gmosaic keyword - phu
    l_key_maskname   =  "MASKNAME" # Mask name keyword - phu
    l_key_dettype    =  "DETTYPE" # Dettype keyword - phu
    l_key_detector   =  "DETECTOR" # Detector keyword - phu
    l_key_instrument =  "INSTRUME" # Instrument keyword - phu
    l_key_nsciext    =  "NSCIEXT" # Number of sci extensions - phu
    l_key_filter1    =  "FILTER1" # FILTER1 keyword - phu
    l_key_filter2    =  "FILTER2" # FILTER1 keyword - phu
    l_key_noddistance=  "NODPIX" # Keyword for N&S pixel shift - phu
    l_key_anodcount  =  "ANODCNT" # Keyword for N&S A NOD possition count - phu
    l_key_bnodcount  =  "BNODCNT" # Keyword for N&S A NOD possition count - phu
    l_key_crpix1     =  "CRPIX1" # Keyword for WCS X reference pixel - SCI
    l_key_refpix     =  "REFPIX" # WCS X reference pix set by GSAPPWAVE - SCI
    l_key_crval1     =  "CRVAL1" # Keyword for WCS X reference pixel value -SCI
    l_key_cd1_1      =  "CD1_1" # Keyword for change WCS along x - SCI
    l_key_ccdname    =  "CCDNAME" # Keyword for CCDNAME - SCI
    l_key_mdfrow     =  "MDFROW" # Corresponding row in MDF for SCI extn - SCI
    l_key_grating    =  "GRATING" # Says which grating was used - PHU
    l_key_grtilt     =  "GRTILT" # Grating tilt; detemines central wavlen - PHU
    l_key_date1      =  "DATE-OBS" # OBSDATE keyword attempt 1 - PHU
    l_key_date2      =  "DATE" # OBSDATE keyword attempt 2- PHU

    # Table Column Names
    l_tabkey_slittilt =  "SLITTILT" # Column name for slit tilt in MOS data
                                    # including NS data

    l_tabkey_slittilt_LS = "slittilt_m" # Column name for slit tilt in LS
                                        # data including NS data - These
                                        # are in the gmos/data directory
                                        # and they have no pa. It's the
                                        # slit tilt in the mask cutter
                                        # Cutter reference

    # Set default variable values here
    #
    status = 1 # Output status of task; 0 - Good; 1 - Bad; assume bad
    isimaging = no # yes for all imaging data only
    norefimages = no # Expects user to supply at least one
    l_sciext = l_sci_ext//"," # SCI extension
    l_varext = l_var_ext//"," # VAR extension
    l_dqext = l_dq_ext//"," # DQ extension
    num_refimages = 0 # Used as a check later on but may not get set by then
    fl_filter_match = yes # No if different filters bewteen imaging frames
    imgpixtype = "real" # Output pixel type for various tasks
    userdefined = no # Whether the user supplies the QE information themselves
    corrapp = "/" # How to apply the QE correction to the input data
    factor_description = "1"//corrapp # If corrapp == *, set this to ""; imaging
    gemzone = "UT" # gemdate time zone for header keyword values
    correction_type = "RELATIVE" # Type of correction to apply
    ns_prefix = "ns" # Prefix for existing corr files that need the ns applying
    num_corrimgs_to_use = -99 # Number of corrimg files to use (reset later on)
    stopdatekey = "old" # Used when parsing QE data from data file
    gmossbiasdate = 841449600 # Date in seconds GMOS-S had new board fitted
                              # 2006-08-31
    gmossccddate = 733622400 # Date GMOS-S in seconds had first new CCD fitted
                             # 2003-04-01
    wav_extname = "WAV"

    # Reassign correction type if an absolute correction is desired
    if (!l_fl_relative_corr) {
        correction_type = "ABSOLUTE"
    }

    # Set up / store default chipgaps/ variables that will be DETECTOR
    # dependant here so they can be assigned later on to another variable
    # allowing as generic as possible code

    # Set detector specific information according to DETTYPE
    num_dettypes = 4

    # det_names defined above: "SDSU II CCD","SDSU II e2v DD CCD42-90",
    #     "S10892", "S10892-N"
    # det_info[?,1] = unbinned chip gap
    # det_info[?,2] = binned chip gap

    # SDSU-II (EEV)
    det_info [1,1] = 37
    det_info [1,2] = 36
    # SDSU II e2v DD CCD42-90 (e2v)
    det_info [2,1] = 37 ##M
    det_info [2,2] = 36 ##M
    # S10892 (Hamamatsu SOUTH)
    det_info [3,1] = 61 ##M CHIP_GAP
    det_info [3,2] = 61 ##M CHIP_GAP
    # S10892-N  (Hamamatsu NORTH)
    det_info [4,1] = 67 #KL
    det_info [4,2] = 67 #KL
    
    # Set up parameter caches of the following tasks
    cache ("gimverify", "gemextn", "gemdate", "gemisnumber", \
        "mkpattern", "imexpr", "blkrep", "imshift")

    # Test the logfile here to start logging:
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
        l_logfile = gmos.logfile
        if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
            l_logfile = "gmos.log"
            printlog ("GQECORR WARNING - Both gqecorr.logfile and "//\
                "gmos.logfile fields are empty.", l_logfile, verbose+)
            printlog ("                   Using default file gmos.log", \
                l_logfile, verbose+)
        }
    }

    # Print start time
    gemdate (zone="local")
    printlog ("\nGQECORR - Started: "//gemdate.outdate//"\n", \
        l_logfile, verbose+)

    ########################
    # Check user set flags #
    ########################

    # Check the settings of l_correct_inimages and l_keep_corrimages
    if (!l_correct_inimages && !l_keep_corrimages) {
        # l_correct_inimages overrides l_keep_corrimages
        l_keep_corrimages = yes
    }

    #########################
    # Start standard checks #
    #########################

    # Perform checks on input/output/refimage/corrimage lists and
    # QE data file/user file/user input; setup required filenames.

    # Create temporary file lists
    inlist = mktemp("tmpinlist")
    outlist = mktemp("tmpoutlist")
    reflist = mktemp("tmpreflist")
    corrimlist = mktemp("tmpcorrimlist")

    ################
    # Check the QE DATA information

    ## Put checks in for other ways these numbers can be input...
    ## Set the userdefined flag? Need to override if put in by hand
    # Check the existance of the file that contains the QE correction factors
    # Check the user has supplied something or left it as the default value
    if ((l_qecorr_data != "") && (stridx(" ", l_qecorr_data) == 0)) {

        # Check to see if it's a file and hence if it exists
        if (!access(l_qecorr_data)) {

            # Make a temporary file
            tmpqedata = mktemp("tmpqedata")

            # File doesn't exist or user supplied numbers
            # Check the users input should be in the form of 1,2,3:4,5,6:8,9
            print (l_qecorr_data) | translit (infile="STDIN", \
                from_string=":,", to_string="\n", delete=no, > tmpqedata)

            # Initiate a counter for error message
            user_qenum_count = 0

            # Loop of the tmpqedata file and check that eash line is a number
            scanfile = tmpqedata
            while (fscan(scanfile,numstring) != EOF) {

                # I think "decimal" will pick all types of decmial numbers
                gemisnumber (instring=numstring, ttest="decimal", \
                    verbose=no)

                # Check status of gemisnumber
                if (gemisnumber.status != 0) {
                    printlog ("GQECORR ERROR - call to GEMISNUMBER returned"//\
                        " a non-zero status. Exiting.", l_logfile, \
                        verbose+)
                    goto clean
                } else {
                   # Check is string is a number
                   if (!gemisnumber.fl_istype) {
                       # It is not a number

                       if (user_qenum_count == 0) {
                           # User may have a supplied a file that couldn't be
                           # accessed
                           printlog ("GQECORR ERROR - "//numstring//\
                               " is neither a number or a file that could"//\
                               " be accssed. Exiting.", l_logfile, verbose+)
                       } else {
                           printlog ("GQECORR ERROR - "//numstring//\
                               "is not a number. Exiting.", \
                               l_logfile, verbose+)
                       }
                       goto clean
                   } # End of checking if input is a number

               } # End of gemisumber check

               user_qenum_count += 1
            } # End of while loop

            # resset the user defined flag
            userdefined = yes

            delete (tmpqedata, verify-, >& "dev$null")

            # Re-write number to tmp file for use later but all of a CCD's info
            # on one line
            print (l_qecorr_data) | translit (infile="STDIN", \
                from_string=":", to_string="\n", delete=no, > tmpqedata)

        } # End of checks to see if user has input numbers

    } else {
        # Exit due to no input
        printlog ("GQECORR ERROR - the qecorr_data parameter has not been"//\
            " set. Exiting", l_logfile, verbose+)
        goto clean
    } # End of qe_data checks

    ################
    # Standard INPUT checks

    # Check for a list file and see if it exists
    # Find the first occurance of '@' in l_input. value of zero means '@' not
    # found in test_string
    test_string = ""
    test_string = l_inimages
    atlocation = stridx("@", test_string)
    if (atlocation > 0) {
        listfile = substr(test_string,(atlocation + 1),strlen(test_string))
        # Check to se if file exists
        if (!access(listfile)) {
            printlog ("GQECORR ERROR - Cannot access file: \""//listfile//\
                "\". Exiting",\
                l_logfile, verbose+)
            goto clean
        }
    }

    # Expand input files into a file
    sections (l_inimages, option="fullname", > inlist)

    # Set up parameters for reading input files
    scanfile = inlist
    filename = ""
    num_inimages = 1

    # Loop over input list of images
    while (fscan(scanfile,filename) != EOF) {

        # Verify image
        gimverify (filename)
        if (gimverify.status > 0) {
             # Image bad/doesn't exist goto to clean
             printlog ("", l_logfile, l_verbose)
             printlog ("GQECORR ERROR - File: "//filename//\
                 " does not exist or is bad", l_logfile, verbose+)
             # Change status to bad
             goto clean
        }
        filename = gimverify.outname

        # Check if the inputs have been gprepare'd
        keypar (filename//"[0]", l_key_gprepare, silent+)
        if (!keypar.found) {
            printlog ("GQECORR EROOR - The input file "//filename//\
                " appears to be a raw file.", l_logfile, verbose+)
            printlog ("                If you wish to use this "//\
                "image, please prepare the data using GPREPARE.", \
                l_logfile, verbose+)
            goto clean
        }

        # Check if the inputs have been gmosaiced
        keypar (filename//"[0]", l_key_gmosaic, silent+)
        if (keypar.found) {
            printlog ("GQECORR EROOR - The input file "//filename//\
                " has been GMOSAIC'ed.", l_logfile, verbose+)
            printlog ("                GQECORR cannot correct GMOSAIC'ed "//\
                "data.", l_logfile, verbose+)
            goto clean
        }

        # If l_refimages == KEYWORD need to access it and populate reflist here
        if (l_refimages == "KEYWORD") { # refimage is keyword in phu of input
            keypar (filename//"[0]", l_key_refimg, silent+)
            if (keypar.found) {
                l_refimg = keypar.value

                # Make sure l_refimg has been set
                if (l_refimg != "" || (stridx(" ", l_refimg) > 0)) { # Not set
                    printlog ("GQECORR ERROR - "//l_key_refimg//" not set "//\
                        "in image "//filename//". Exiting", l_logfile, \
                        verbose+)
                    goto clean
                } else { # l_refimg set
                    # Write l_refimg to reflist for standard checking later on
                    print (l_refimg, >> reflist)
                    printlog ("Input: "//filename//" Refimage: "//l_refimg, \
                        l_logfile, l_verbose)
                    goto clean
                }
            } else {
                printlog ("GQECORR ERROR - "//l_key_refimg//" keyword not "//\
                    "found in image "//filename, l_logfile, verbose+)
                goto clean
            }
        }

        inimage_files[num_inimages] = filename//".fits"
        num_inimages += 1

    } # End of while loop verifying images and populating inimage_files[]

    # Fix the num_inimages counter to be correct
    num_inimages -= 1

    ################
    # Standard OUTPUT checks / set up

    # Set prefix usage flag
    l_fl_useprefix = yes

    # Check for an output list
    if (l_outimages != "" && (stridx(" ", l_outimages) == 0)) {

        # Check for a list file and see if it exists
        # Find the first occurance of '@' in l_input. value of 0 means '@' not
        # found in test_string
        test_string = ""
        test_string = l_outimages
        atlocation = stridx("@", test_string)
        if (atlocation > 0) {
            listfile = substr(test_string,(atlocation + 1),strlen(test_string))
            # Check to se if file exists
            if (!access(listfile)) {
                printlog ("GQECORR ERROR - Cannot access file: \""//listfile//\
                    "\". Exiting",\
                    l_logfile, verbose+)
                goto clean
                }
        }

        # Expand output files into a list
        sections (l_outimages, option="fullname", > outlist)

        # Check the number of output files is equal to the number of inputs
        sections ("@"//outlist, option="nolist")
        if (num_inimages == int(sections.nimages)) {
            l_fl_useprefix = no
        } # End of if statement checking number of in and out files

        # Delete the file if required to prefix the prefix
        if (l_fl_useprefix) {
            delete (outlist, verify-, >& "dev$null")
        }

    } # End of output list checks

    # Check if prefix is to be used
    if (l_fl_useprefix) {
        # Use prefix and input list
        scanfile = inlist
        while (fscan(scanfile, filename) != EOF) {
            fparse (filename)
            print (l_outpref//fparse.root//".fits", >> outlist)
        }
        scanfile = ""
    }

    # Set up parameters for reading output file list
    scanfile = outlist
    filename = ""
    num_outimages = 1

    # Loop over outlist of images to populate outimage_files and check they
    # don't exist (only if correcting the data exit)
    while (fscan(scanfile,filename) != EOF) {
        # Verify image
        gimverify (filename)
        if (gimverify.status != 1 && l_correct_inimages) {
             # If image doesn't exit
             printlog (" ",l_logfile, l_verbose)
             printlog ("GQECORR ERROR - File: "//filename//" already exists",\
                 l_logfile, verbose+)
             # Change status to bad
             goto clean
        }
        filename = gimverify.outname
        outimage_files[num_outimages] = filename//".fits"
        num_outimages += 1
    } # End of while loop verifying and populating outimage_files[]

    # Fix the num_outimages counter to be correct
    num_outimages -= 1

    ################
    # Standard REFIMAGES checks

    # Check the value of the refimages parameter
    if (l_refimages == "" || (stridx(" ", l_refimages) != 0)) {
        # Set flag for test later on (it is aloud for imaging data)
        norefimages = yes
    } else { # l_refimages set
        if (l_refimages != "KEYWORD") { # reflist will already exist if KEYWORD
           # Check for a list file and see if it exists
           # Find the first occurance of '@' in l_input. value of zero means
           # '@' not found in test_string
           test_string = ""
           test_string = l_refimages
           atlocation = stridx("@", test_string)
           if (atlocation > 0) {
               listfile = substr(test_string,(atlocation + 1), \
                   strlen(test_string))
               # Check to se if file exists
                if (!access(listfile)) {
                    printlog ("GQECORR ERROR - Cannot access file: \""//\
                        listfile//"\". Exiting",\
                        l_logfile, verbose+)
                    goto clean
                }
            }

            # Expand reference files into a list
            sections (l_refimages, option="fullname", > reflist)
        }

        # Set up parameters for reading input files
        scanfile = reflist
        filename = ""
        num_refimages = 1

        # Loop over input list of images
        while (fscan(scanfile,filename) != EOF) {
            # Verify image
            gimverify (filename)
            if (gimverify.status > 0) {
                 # Image bad/doesn't exist goto to clean
                 printlog ("", l_logfile, l_verbose)
                 printlog ("GQECORR ERROR - File: "//filename//\
                     " does not exit or is bad", l_logfile, verbose+)
                 # Change status to bad
                 goto clean
            }
            filename = gimverify.outname

            file_extn = substr(filename,(strlen(filename)-4),strlen(filename))
            if (file_extn == ".fits") {
                refimage_files[num_refimages] = \
                    substr(filename,1,(strlen(filename)-5))
            } else {
                refimage_files[num_refimages] = filename
            }
            num_refimages += 1
        } # End of while loop verifying images and populating refimage_files[]

        # Fix the num_refimages counter to be correct
        num_refimages -= 1

        # Check that the num_refimages is equal to num_inimages or 1
        if ((num_refimages != 1) && (num_refimages != num_inimages)) {
            printlog ("GQECORR ERROR - The number of reference files "//\
                "supplied is not one and is not equal to the number of "//\
                "input images", l_logfile, verbose+)
            goto clean
        }
    } # End of standard REFIMAGE checks

    ################
    # Standard CORRIMAGE checks / setup (Regardless of l_keep_corrimages)
    ## Need to check if they aleardy exist and are to be used?

    # Set prefix usage flag
    l_fl_useprefix = yes
    nocorrimages = no
    num_incorrimages = -99

    # Check for an output list
    if ((l_corrimages != "") && (stridx(" ", l_corrimages) == 0)) {

        # Check for a list file and see if it exists
        # Find the first occurance of '@' in l_input. value of zero means '@'
        # not found in test_string
        test_string = ""
        test_string = l_corrimages
        atlocation = stridx("@", test_string)
        if (atlocation > 0) {
            listfile = substr(test_string,(atlocation + 1),strlen(test_string))
            # Check to se if file exists
            if (!access(listfile)) {
                printlog ("GQECORR ERROR - Cannot access file: \""//listfile//\
                "\". Exiting",\
                l_logfile, verbose+)
                goto clean
            }
        }

        # Expand input correction files into corrimlist
        sections (l_corrimages, option="fullname", > corrimlist)

        # Check the number of corrim output files is equal to the
        # number of input reference files or the input files supplied
        # If neither are equal to the number of corrimages use prefix
        sections ("@"//corrimlist, option="nolist")
        num_incorrimages = section.nimages

        if ((num_inimages == num_incorrimages) || \
            (num_refimages == num_incorrimages) || \
                (num_incorrimages == 1)) {
                    # Don't use the prefix
                    l_fl_useprefix = no
        }

        # Delete the list so that the prefix can be added and written to
        # disk later on
        if (l_fl_useprefix) {
            delete (corrimlist, verify-, >& "dev$null")
        }

    } else {
        # nocorrimages flag for use as a check later
        nocorrimages = yes
    } # End of corrimage list checks

    # The tests later on are more robust and can, therefore, make a small
    # assumption here as only the correction information names are being
    # inititated. So, one can check the first image and see if it's
    # imaging data. Only all imaging or all spectroscopy is allowed.
    # Any failures in this call will be caught in the tests later on
    thisobsmode = ""
    keypar (inimage_files[1]//"[0]", l_key_obsmode, silent+)
    thisobsmode = keypar.value

    # Check for GMOS-N and GMOS-S data taken with EEV CCDs
    # in imaging mode - this cannot as yet be supported.
    # Read DETTYPE of first image too
    thisdettype = ""
    keypar (inimage_files[1]//"[0]", l_key_dettype, silent+)
    thisdettype = keypar.value
    # Quick way of exiting
#    if (thisobsmode == "IMAGE") {
    if (thisobsmode == "IMAGE" && 
        (thisdettype != det_names[3] && thisdettype != det_names[4])) {
#        printlog ("GQECORR ERROR - QE correction for GMOS imaging data "//\
#            "is not currently supported.\n"//\
#            "                Exiting.", \
#            l_logfile, verbose+)
        printlog ("GQECORR ERROR - QE correction for GMOS imaging data "//\
            "with a "//l_key_dettype//" of \""//det_names[1]//"\"\n"//\
            "                is not currently supported. Exiting.",
            l_logfile, verbose+)
    ## Delete before release
#        printlog ("\n***** This is only printed out currently as a "//\
#            "warning *****"//\
#            "\n      only one FILTER combination is accepted to "//\
#            "allow"//\
#            "\n      the imaging loop to be checked. "//\
#            "\n      See "//l_qecorr_data//" for the allowed filters\n", \
#            l_logfile, verbose+)
    ## Uncomment before release
        goto clean
    }

    # Check if prefix is to be used
    if (l_fl_useprefix) {
        # Need to get the naming convention correct. For imaging data, the
        # the correction name uses the name of the input file
        # For spectroscopy it uses the reference image.
        # If the user supplies the files then it is assumed they know what
        # to do / expect.

        if (norefimages || (thisobsmode == "IMAGE")) {
            # Use prefix and inlist
            scanfile = inlist
        } else {
            # Use prefix and reflist
            scanfile = reflist
        }
        while (fscan(scanfile, filename) != EOF) {
            fparse (filename)
            print (l_corrimpref//fparse.root//".fits", >> corrimlist)
        }
        scanfile = ""
    }

    # Set up parameters for reading output file list
    scanfile = corrimlist
    filename = ""
    num_corrimages = 1
    num_corrim_exists = 0

    # Loop over corrimlist of images to populate corrimage_files and check
    # they don't exist / do exisat and act accordingly
    while (fscan(scanfile,filename) != EOF) {

        # Strip the .fits file extension (appended
        # with .fits or .dat later on)
        file_extn = substr(filename,(strlen(filename)-4),strlen(filename))
        if (file_extn == ".fits") {
            corrimage_files[num_corrimages] = \
                substr(filename,1,(strlen(filename)-5))
        } else {
            corrimage_files[num_corrimages] = filename
        }

        # Initiate the flag that records is an correction file already exists
        corrim_exists[num_corrimages] = no

        # Verify image
        gimverify (filename)
        if (gimverify.status == 0) {
             # If image exits and is an MEF
             printlog ("GQECORR WARNING - Correction file "//filename//\
                 " already exists.", \
                  l_logfile, verbose+)

             # Check if created by GQECORR
             keypar (filename//"[0]", l_key_gqecorr, silent+)
             if (keypar.found) {
                 # Used later on to determine whether to obtain QE info or not
                 corrim_exists[num_corrimages] = yes
                 num_corrim_exists += 1
             } else {
                 # Can't use it
                 printlog ("                  However, it was not "//\
                     "creaeted by GQECORR.\n"//\
                     "                  Cannot use this image for QE "//\
                     "correction.\n"//\
                     "                  Continuing anyway.", \
                     l_logfile, verbose+)
             }

        } else if (gimverify.status == 1) {
            # Image doesn't exist

            # Check if a .dat file exits
            if (access(corrimage_files[num_corrimages]//".dat")) {
                printlog ("GQECORR WARNING - Correction file "//\
                    corrimage_files[num_corrimages]//".dat already exists.",  \
                    l_logfile, verbose+)

                 # Used later on to determine whether to obtain QE info or not
                 corrim_exists[num_corrimages] = yes
                 num_corrim_exists += 1
                 # Can't check validity of file, have to assume it's OK.
             }
        } else {
             # Image exists but isn't an MEF
             printlog (" ",l_logfile, l_verbose)
             printlog ("GQECORR ERROR - Correction file "//filename//\
                 " already exists and is not a MEF file", \
                 l_logfile, verbose+)
             goto clean
        }

        num_corrimages += 1
    } # End of while loop verifying and populating outimage_files[]

    # Fix the num_outimages counter to be correct
    num_corrimages -= 1

    if (debug) {
        printlog ("num_corrimages: "//num_corrimages//\
            "\nnum_corrim_exists: "//num_corrim_exists//\
            "\nnum_inimages: "//num_inimages//\
            "\nnum_incorrimages: "//num_incorrimages, \
            l_logfile, verbose+)
    }

    # Do a check of the number of existing correction files
    if (num_corrim_exists != 0) {
        # Has to be the same number as input or equal to 1
        if (((num_corrim_exists != num_inimages) && \
            (num_corrim_exists != 1))) {
                printlog ("GQECORR WARNING - The number of existing "//\
                    "correction files is not one and does not match the "//\
                    "\n                   number of input files. Will "//\
                    "continue anyway but will not "//\
                    "use all the correction "//\
                    "\n                   files that already exist.", \
                    l_logfile, verbose+)

            # Reset num_corrim_exists tp zero
            num_corrim_exists = 0
        } else {
            # switch refimage_files[] to corrimage_files[].
            # Only if thisobsmode set ealier is not equsl to IMAGE
            # Cannot run the refimage checks on IMAGE correction data
            if (thisobsmode != "IMAGE") {
                # Use correction files that already exist
                printlog ("GQECORR - Using the correction files that "//\
                    "already exist.", l_logfile, verbose+)

                # Reset the norefimages flag to no - so the standard checks
                # can be done and the number od re_images too
                norefimages = no
                num_refimages = num_corrimages

                # Repopulate the refimage_files[]
                for (ii = 1; ii <= num_corrimages; ii += 1) {
                    refimage_files[ii] = corrimage_files[ii]
                }
            }
        }
    }


    ##########################
    # End of standard checks #
    ##########################

    ########################################
    # Start of specific requirement checks #
    ########################################

    # Populate the binning and obsmode arrays for both refimages and inimages
    #     Check their validity and check that the modes, detector etc. match.
    # For the reference images check that specific tasks have been run on them
    #     MOS - GSCUT, IFU - GFEXTRACT
    # For the reference image, if doing a full correction check that the
    #     database files exist.

    ################
    # REFIMAGE checks

    # Only needs to be done if there are reference images defined
    if (!norefimages) {

        # Loop of the reference images
        for (ii = 1; ii<= num_refimages; ii += 1) {

            # Set the local input parameter to be the current image
            l_refimg = refimage_files[ii]//".fits"
            if (debug) {
                printlog ("REFIMAGE tests... working on: "//l_refimg, \
                    l_logfile, verbose+)
            }

            # Obtain the number of l_sci_ext frames in the reference image
            gemextn (l_refimg, process="expand", check="exists", \
                extver="1-", extname=l_sci_ext, outfile="dev$null", \
                logfile=l_logfile, verbose=l_verbose)

            # Check if GEMEXTN call was succesfull
            if (gemextn.status == 0) { # GEMEXTN was successful
                if (gemextn.count > 0) { # l_sci_ext exist
                    # Obtain the number of l_sci_ext extensions in l_refimg
                    num_ref_sci = gemextn.count
                } else { # l_sci_ext don't exist
                    # No science extensions or l_sci_ext is incorrectly set
                    printlog ("GQECORR ERROR - No "//l_sci_ext//\
                        " extensions found in image"//l_refimg//\
                        " (Error 011). Exiting", l_logfile, verbose+)
                    goto clean
                }
            }  else {
                # GEMEXTN returned an error
                printlog ("GQCORR ERROR - Call to GEMEXTN returned non-zero"//\
                    " status for image "//l_refimg//". Exiting. (Error "//\
                    "011)", logfile, verbose+)
                goto clean
            }

            # Compare num_ref_sci to NSCIEXT (so it can be used with
            #     confidence later
            keypar (l_refimg//"[0]", l_key_nsciext, silent+)
            if (keypar.found) {
                if (int(keypar.value) != num_ref_sci) {
                    printlog ("GQECORR ERROR - The number of "//l_sci_ext//\
                        " extensions does not match the value of "//\
                        l_key_nsciext//". Exiting.", l_logfile, verbose+)
                    goto clean
                }
            } else {
                printlog ("GQECORR ERROR - "//l_key_nsciext//" keyword not "//\
                    " found in reference image"//l_refimg//". Exiting.", \
                    l_logfile, verbose+)
                    goto clean
            }

            # Obtain the OBSMODE
            keypar (l_refimg//"[0]", l_key_obsmode, silent+)
            thisobsmode = ""
            if (keypar.found) {
                thisobsmode   = keypar.value
                refobsmode[ii] = thisobsmode

                # Unreadable obsmode
                if ((thisobsmode != "MOS") && (thisobsmode != "IFU") && \
                    (thisobsmode != "LONGSLIT") && (thisobsmode != "IMAGE")) {
                        printlog ("GQECORR ERROR - "//l_key_obsmode//" not "//\
                            " recognised in image "//l_refimg//". Exiting", \
                            l_logfile, verbose+)
                        goto clean
                } else if (thisobsmode == "IMAGE") {
                    # Imaging data
                    printlog ("GQECORR ERROR - Refimage "//l_refimg//\
                        " has an OBSMODE of IMAGE. "//l_key_obsmode//\
                        " must be LONGSLIT/MOS/IFU for the reference image."//\
                        " Reference images are not required for imaging"//\
                        " data", l_logfile, verbose+)
                    goto clean
                }

            } else {
                printlog ("GQECORR ERROR - "//l_key_obsmode//" keyword not "//\
                    " found in reference image"//l_refimg//". Exiting.", \
                    l_logfile, l_verbose)
            }

            # Refimage will be spectroscopic by this point

            ################
            # Specific keyword checks by mode

            if (corrim_exists[ii]) {
                # Skip the next few tests
                goto SKIP_TESTS
            }

            # Check the database exists if using fullcorr
            if (l_fl_fullcorr && ii == 1) {
                # Check access to database directory
                l_database = (osfn(l_database))

                if (!access(l_database)) {
                    printlog ("ERROR - GQECORR: Cannot access database "//\
                        "directory "//l_database, l_logfile, verbose+)
                    goto clean
                }
            }

            # For IFU check gfextract has been run
            if (thisobsmode == "IFU") {
                keypar (l_refimg//"[0]", l_key_gfextract, silent+)
                if (!keypar.found) {
                    # gscut keyword not found
                    printlog ("GQECORR ERROR - The reference image: "//\
                        l_refimg//" appears not to have been run through "//\
                        "GFEXTRACT. Exiting", l_logfile, verbose+)
                    goto clean
                }
                ## Double check the presence of MDF
            } # End of IFU checks

            # For MOS data check gscut has been run
            if (thisobsmode == "MOS") {
                keypar (l_refimg//"[0]", l_key_gscut, silent+)
                if (!keypar.found) {
                    # gscut keyword not found
                    printlog ("GQECORR ERROR - The reference image: "//\
                        l_refimg//" appears not to have been run through "//\
                        "GSCUT but it is MOS/NS data. Exiting", \
                        l_logfile, verbose+)
                    goto clean
                }

            } # End of MOS checks

SKIP_TESTS:

            # Double check the presence of MDF for MOS and IFU data
            if ((thisobsmode == "MOS") || (thisobsmode == "IFU")) {
                # Initiate tinfo.tabletype
                tinfo.tbltype = ""

                # read the MDF table imformation
                tinfo (l_refimg, ttout=no, >& "dev$null")

                # Test the tbltype
                if (tinfo.tbltype != "fits") {
                    printlog ("GQECOR ERROR - OBSMODE for "//l_refimg//\
                        " is "//thisobsmode//" and a MDF is required "//\
                        "but no MDF is present in "//l_refimg//". Exiting.", \
                        l_logfile, verbose+)
                    goto clean
                }
            }

            # Check if GQECORR has already applied a N&S shift to the reference
            # image - will only be the case if corrim_exists[] is true
            keypar (l_refimg//"[0]", l_key_gqenod, silent+)
            if (keypar.found) {
                refim_gqenod[ii] = yes
            } else {
                refim_gqenod[ii] = no
            }

            # Depending on desired correction -
            #     l_fl_fullcorr - dispcorr or fitcoord solution check the
            #         database files exist
            #     !l_fl_fullcorr - check gsappwav has been run
            if (!l_fl_fullcorr) {
                # Check that gsappwave has been run on the refimg
                keypar (l_refimg//"[0]", l_key_gsappwav, silent=yes)
                l_task = "GSAPPWAVE"

            } else {
                # Full correction

                # Check that gswavelength has been run on the refimg
                keypar (l_refimg//"[0]", l_key_gswave, silent=yes)
                l_task = "GSWAVELENGTH"
            }

            # If using a QECORR created correction image then
            # do not need to check for GSAPPWAV keyword
            if (!keypar.found && !corrim_exists[ii]) {
                # gsappwav keyword not found
                printlog ("GQECORR ERROR - The reference image: "//\
                    l_refimg//" appears not to have been run through "//\
                    l_task//". Exiting", l_logfile, verbose+)
                goto clean

            } # End of WCS checks

            # CCDSUM (Use first l_sciext and store it for later)
            xbin = 0
            ybin = 0
            keypar (l_refimg//"["//l_sciext//"1]", l_key_ccdsum, silent+)
            if (keypar.found) {
                print (keypar.value) | scanf ("%d %d", xbin, ybin)
                xrefbin[ii] = xbin
                yrefbin[ii] = ybin
            } else {
                printlog ("GQECORR ERROR - "//l_key_ccdsum//" keyword not "//\
                    "found in image "//l_refimg, l_logfile, verbose+)
            } # End of binning checks

        } # End of ii loop over refimages
    } # End of refimages specific checks

    ################
    # INPUT checks

    kk = 1
    count_imgobsmode = 0
    num_ns_corrim_exists = 0
    num_ns_infiles = 0

    for (ii = 1; ii<= num_inimages; ii += 1) {

        # Set the local input parameter to be the current image
        l_input = inimage_files[ii]

        if (debug) {
            printlog ("Doing checks on: "//l_input, l_logfile, verbose+)
        }

        # Only do the following check if correcting the output image
        if (l_correct_inimages) {

            # Initiate the fail counter
            fail_counter = 0

            # Check if image has been QE corrected before (no need to check
            # if not correcting the images)
            keypar (l_input//"[0]", l_key_gqecorr, silent+)
            if (keypar.found){
                printlog ("GQECORR ERROR - Image "//l_input//" has already"//\
                    " been QE correcteed by GQECORR. Exiting.", l_logfile, \
                    verbose+)
                goto clean
            } # End of check for previous GQECORR correction

            # Check if image has been QE corrected before - GIREDUCE Spectra
            keypar (l_input//"[0]", l_key_corrim, silent+)
            if (keypar.found){
                printlog ("GQECORR ERROR - Image "//l_input//" has already"//\
                    " been QE correcteed by GIREDUCE. Exiting.", l_logfile, \
                    verbose+)
                goto clean
            } # End of check for previous QE spectroscopy correction

            # Check if image has been QE corrected before - GIREDUCE Imaging
            keypar (l_input//"["//l_sciext//"1]", l_key_imgfactor, silent+)
            if (keypar.found){
                printlog ("GQECORR ERROR - Image "//l_input//" has already"//\
                    " been QE correcteed by GIREDUCE. Exiting.", l_logfile, \
                    verbose+)
                goto clean
            } # End of check for previous QE spectroscopy correction

            # Check if trimmed corrected
#            keypar (l_input//"[0]", l_key_trimmed, silent+)
#            if (!keypar.found){
#                # Image not overscan corrected - exit with error
#                printlog ("GQECORR ERROR - Image "//l_input//\
#                    " has not been TRIMMED.", logfile, verbose+)
#                fail_counter += 1
#
#            } # End of trimmed check

            # Image should be BIAS or DARK or OVERSCAN subtracted!
            # Check if bias corrected
            keypar (l_input//"[0]", l_key_biasimg, silent+)
            if (!keypar.found){
                # If not, check if dark corrected
                keypar (l_input//"[0]", l_key_darkimg, silent+)
                if (!keypar.found) {
                    keypar (l_input//"["//l_sciext//"1]", l_key_overscan, \
                        silent+)
                    if (!keypar.found) {
                        # Image not bias or dark subtracted or overscan
                        # corrected - exit with error
                        printlog ("GQECORR ERROR - Image "//\
                            l_input//\
                            " has not been BIAS or DARK subtracted or "//\
                            " overscan corrected.", \
                            logfile, verbose+)
                        fail_counter += 1
                    }
                }
            } # End of BIAS/DARK subtraction checks

            # If there is a failure print this reminder then fail
            if (fail_counter > 0) {
                printlog ("QCECORR ERROR - Images must be BIAS or "//\
                    "DARK subtracted or overscan corrected to apply QE "//\
                    "correction. Exiting.", logfile, verbose+)
                goto clean
            }
        } # End of checks if trying to correct the input images

        # Obtain the number of l_sci_ext frames in the input image
        gemextn (l_input, process="expand", check="exists", \
            extver="1-", extname=l_sci_ext, outfile="dev$null", \
            logfile=l_logfile, verbose=l_verbose)

        # Check if GEMEXTN call was succesfull
        if (gemextn.status == 0) { # GEMEXTN was successful
            if (gemextn.count > 0) { # l_sci_ext exist
                # Obtain the number of l_sci_ext extensions in l_refimg
                num_sci = gemextn.count
            } else { # l_sci_ext don't exist
                # No science extensions or l_sci_ext is incorrectly set
                printlog ("GQECORR ERROR - No "//l_sci_ext//\
                    " extensions found in image "//l_input//\
                    " (Error 012). Exiting", l_logfile, verbose+)
                goto clean
            }
        }  else { # GEMEXTN returned an error
            printlog ("GQCORR ERROR - Call to GEMEXTN returned non-zero"//\
                " status for image "//l_input//". Exiting. (Error "//\
                "012)", logfile, verbose+)
            goto clean
        } # End of GEMEXTN checks

        # Compare num_ref_sci to NSCIEXT (so it can be used with
        #     confidence later)
        keypar (l_input//"[0]", l_key_nsciext, silent+)
        if (keypar.found) {
            if (int(keypar.value) != num_sci) {
                printlog ("GQECORR ERROR - The number of "//l_sci_ext//\
                    " extensions does not match the value of "//\
                    l_key_nsciext//". Exiting", l_logfile, verbose+)
                goto clean
            }
        } else {
            printlog ("GQECORR ERROR - "//l_key_nsciext//" keyword not "//\
                " found in input image"//l_inout//". Exiting.", \
                l_logfile, l_verbose)
            goto clean
        } # End of NSCIEXT check

        # CCDSUM (Use first l_sciext and store for later use)
        xbin = 0
        ybin = 0
        keypar (l_input//"["//l_sciext//"1]", l_key_ccdsum, silent+)
        if (keypar.found) {
            print (keypar.value) | scanf ("%d %d", xbin, ybin)
            xinbin[ii] = xbin
            yinbin[ii] = ybin
        } else {
            printlog ("GQECORR ERROR - "//l_key_ccdsum//" keyword not "//\
                "found in image "//l_input, l_logfile, verbose+)
        } # End of chceking and storing CCDSUM

        # OBSMODE
        inobsmode[ii] = ""
        keypar (l_input//"[0]", l_key_obsmode, silent+)
        thisobsmode = ""
        if (keypar.found) {
            thisobsmode   = keypar.value
            inobsmode[ii] = thisobsmode

            # Unreadable obsmode
            if ((thisobsmode != "MOS") && (thisobsmode != "IFU") && \
                (thisobsmode != "LONGSLIT") && (thisobsmode != "IMAGE")) {
                    printlog ("GQECORR ERROR - "//l_key_obsmode//" not "//\
                        "recognised in image "//l_input//". Exiting", \
                        l_logfile, verbose+)
                    goto clean
            } else if (thisobsmode != "IMAGE") {

                #######
                # Spectroscopic data checks

                # Check here for a reference image || correection images
                if (norefimages && nocorrimages) {
                    # refimages and corrimages not supplied
                    printlog ("GQECORR ERROR - Sepctroscopic data input ("//\
                        l_input//"). However, refimages\n"//\
                        "                and corrimages "//\
                        "parameters not set. Exiting",\
                        l_logfile, verbose+)
                    goto clean
                } else if (norefimages && !nocorrimages && \
                    (num_corrim_exists == 0)) {
                        # refimages not supplied, corrimages supplied
                     printlog ("GQECORR ERROR - Sepctroscopic data input ("//\
                        l_input//"). However, refimages\n"//\
                        "                parameter not set and the supplied"//\
                        " corrimages do not exist or were not\n"//\
                        "                created by QECORR. Exiting",\
                        l_logfile, verbose+)
                    goto clean
                } else { # Reference images supplied

                    # Set the refence image index kk - default is 1
                    if (num_refimages != 1) {
                        # If a reference image for every input image -
                        # confirmed in standard checks
                        kk = ii
                    }

                    # Assign reference image
                    l_refimg = refimage_files[kk]

                    # Check the obsmode against the refimage obsmode
                    if (thisobsmode != refobsmode[kk]) {
                        printlog ("GQECORR ERROR - OBSMODE mismatch between"//\
                            " "//l_input//" and "//l_refimg//". Exiting", \
                            l_logfile, verbose+)
                        goto clean
                    }
                } # End of check for REFIMG and comparison of OBSMODE

                # Check MASKNAME, DETTYPE, DETECTOR and INSTRUME
                # against l_refimage
                # the keywords are read for l_refimage too
                indettype = ""
                indetector = ""
                ininstrument = ""

                # Obtain the MASKNAME of the input image
                keypar (l_input//"[0]", l_key_maskname, silent+)
                if (keypar.found) {
                    inmaskname = keypar.value
                    # Obtain the MASKNAME of the refimage
                    keypar (l_refimg//"[0]", l_key_maskname, silent+)
                    if (keypar.found) {
                        # Check if they are the same
                        if (keypar.value != inmaskname) {
                            if (thisobsmode == "LONGSLIT") {
                                errstr="WARNING"
                            } else {
                                errstr="ERROR"
                            }
                            printlog ("GQECORR "//errstr//" - "//\
                                l_key_maskname//" mismatch between "//\
                                l_input, l_logfile, verbose+)
                            printlog ("  and "//l_refimg, l_logfile, verbose+)
                            if (thisobsmode != "LONGSLIT") goto clean
                        }
                    } else {
                        printlog ("GQECORR ERROR - "//l_key_maskname//\
                            " keyword not found in image "//l_refimg, \
                            l_logfile, verbose+)
                        goto clean
                    }
                } else {
                    printlog ("GQECORR ERROR - "//l_key_maskname//" keyword"//\
                        " not found in image "//l_input, l_logfile, verbose+)
                    goto clean
                } # End of MASKNAME checks

                # Obtain the DETTYPE of the input image
                keypar (l_input//"[0]", l_key_dettype, silent+)
                if (keypar.found) {
                    indettype = keypar.value
                    # Obtain the DETTYPE of the refimage
                    keypar (l_refimg//"[0]", l_key_dettype, silent+)
                    if (keypar.found) {
                        # Check if they are the same
                        if (keypar.value != indettype) {
                            printlog ("QECORR ERROR - "//l_key_dettype//\
                                " mismatch between "//l_input//" and "//\
                                l_refimg, l_logfile, verbose+)
                            goto clean
                        }
                    } else {
                        printlog ("GQECORR ERROR - "//l_key_dettype//\
                            " keyword not found in image "//l_refimg, \
                            l_logfile, verbose+)
                        goto clean
                    }
                } else {
                    printlog ("GQECORR ERROR - "//l_key_dettype//" keyword "//\
                        "not found in image "//l_input, l_logfile, verbose+)
                    goto clean
                } # End of DETTYPE checks

                # Obtain the DETECTOR of the input image
                keypar (l_input//"[0]", l_key_detector, silent+)
                if (keypar.found) {
                    indetector = keypar.value
                    # Obtain the DETTYPE of the refimage
                    keypar (l_refimg//"[0]", l_key_detector, silent+)
                    if (keypar.found) {
                        # Check if they are the same
                        if (keypar.value != indetector) {
                            printlog ("GQECORR ERROR - "//l_key_detector//\
                                " mismatch between "//l_input//" and "//\
                                l_refimg, l_logfile, verbose+)
                            goto clean
                        }
                    } else {
                        printlog ("GQECORR ERROR - "//l_key_detector//\
                            " keyword not found in image "//l_refimg, \
                            l_logfile, verbose+)
                        goto clean
                    }
                } else {
                    printlog ("GQECORR ERROR - "//l_key_detector//" keyword"//\
                        " not found in image "//l_input, l_logfile, verbose+)
                    goto clean
                } # End of DETECTOR checks

                # Obtain the INSTRUMENT of the input image
                keypar (l_input//"[0]", l_key_instrument, silent+)
                if (keypar.found) {
                    ininstrument = keypar.value
                    # Obtain the INSTRUMENT of the refimage
                    keypar (l_refimg//"[0]", l_key_instrument, silent+)
                    if (keypar.found) {
                        # Check if they are the same
                        if (keypar.value != ininstrument) {
                            printlog ("QECORR ERROR - "//l_key_instrument//\
                                " mismatch between "//l_input//" and "//\
                                l_refimg, l_logfile, verbose+)
                            goto clean
                        }
                    } else {
                        printlog ("GQECORR ERROR - "//l_key_instrument//\
                            " keyword not found in image "//l_refimg, \
                            l_logfile, verbose+)
                        goto clean
                    }
                } else {
                    printlog ("GQECORR ERROR - "//l_key_instrument//\
                        " keyword not found in image "//l_input, \
                        l_logfile, verbose+)
                    goto clean
                } # End of INSTRUMENT checks


                # Check that grating the same for input and refimage
                # Obtain the GRATING of the refimage
                keypar (l_input//"[0]", l_key_grating, silent+)
                if (keypar.found) {
                    ingrating = keypar.value
                    # Obtain the GRATING of the refimage
                    keypar (l_refimg//"[0]", l_key_grating, silent+)
                    if (keypar.found) {
                        # Check if they are the same
                        if (keypar.value != ingrating) {
                            printlog ("QECORR ERROR - "//l_key_grating//\
                                " mismatch between "//l_input//" and "//\
                                l_refimg, l_logfile, verbose+)
                            goto clean
                        }
                    } else {
                        printlog ("GQECORR ERROR - "//l_key_grating//\
                            " keyword not found in image "//l_refimg, \
                            l_logfile, verbose+)
                        goto clean
                    }
                } else {
                    printlog ("GQECORR ERROR - "//l_key_grating//\
                        " keyword not found in image "//l_input, \
                        l_logfile, verbose+)
                    goto clean
                } # End of GRATING checks

                # Check that grating tilt are the same for input and refimage
                # Obtain the GRTILT of the refimage
                keypar (l_input//"[0]", l_key_grtilt, silent+)
                if (keypar.found) {
                    ingrtilt = keypar.value
                    # Obtain the GRTILT of the refimage
                    keypar (l_refimg//"[0]", l_key_grtilt, silent+)
                    if (keypar.found) {
                        # Check if they are the same
                        if (keypar.value != ingrtilt) {
                            # This was previously a fatal error but that
                            # prevented using an arc at a slightly different
                            # wavelength (eg. dither pos) for quick-look
                            # reduction, notably in gmosexamples MOS.
                            printlog ("QECORR WARNING - "//l_key_grtilt//\
                                " mismatch between "//l_input//" and "//\
                                l_refimg, l_logfile, verbose+)
                            # goto clean
                        }
                    } else {
                        printlog ("GQECORR ERROR - "//l_key_grtilt//\
                            " keyword not found in image "//l_refimg, \
                            l_logfile, verbose+)
                        goto clean
                    }
                } else {
                    printlog ("GQECORR ERROR - "//l_key_grtilt//\
                        " keyword not found in image "//l_input, \
                        l_logfile, verbose+)
                    goto clean
                } # End of GRTILT checks

                # OBSERVATION date checks - Only needed really for GMOS-S
                obsdate_inkey_value = ""
                keypar (l_input//"[0]", l_key_date1, silent+)
                if (keypar.found) {
                    obsdate_inkey_value = keypar.value
                } else {
                    keypar (l_input//"[0]", l_key_date2, silent+)
                    if (keypar.found) {
                        obsdate_inkey_value = keypar.value
                    } else {
                        printlog ("GQECORR ERROR - Cannot find "//\
                            l_key_date1//" or "//l_key_date2//" in file:"//\
                            l_input//". Exiting", l_logfile, verbose+)
                        goto clean
                    }
                }

                # Only need to do actuall OBSDATE checks if GMOS-S
                if (obsdate_inkey_value != "" && ininstrument == "GMOS-S") {
                    # Convert OBS date to second for easier before-after
                    # date comparison
                    cnvtsec (obsdate_inkey_value, "00:00:00") | \
                        scan (inobsdate)

                    # OBSERVATION refimage date checks
                    obsdate_refkey_value = ""
                    keypar (l_refimg//"[0]", l_key_date1, silent+)
                    if (keypar.found) {
                        obsdate_refkey_value = keypar.value
                    } else {
                        keypar (l_refimg//"[0]", l_key_date2, silent+)
                        if (keypar.found) {
                            obsdate_refkey_value = keypar.value
                        } else {
                            printlog ("GQECORR ERROR - Cannot find "//\
                                l_key_date1//" or "//l_key_date2//\
                                " in file:"//\
                                l_refimg//". Exiting", l_logfile, verbose+)
                            goto clean
                        }
                    }

                    # Compare the in file and the ref images observation dates
                    # They need to be on the same side of (31-08-2006) for
                    # GMOS-S only
                    if (obsdate_refkey_value != "") {
                        # Convert OBS date to second for easier before-after
                        # date comparison
                        cnvtsec (obsdate_refkey_value, "00:00:00") | \
                            scan (refobsdate)

                        # Read the OBSERVATION date from the reference image
                        if (!(((inobsdate > gmossbiasdate) && \
                            (refobsdate > gmossbiasdate)) || \
                            ((inobsdate < gmossbiasdate) && \
                            (refobsdate < gmossbiasdate)))){
                            printlog ("GQECORR ERROR - OBSSERVATION dates "//\
                                "of "//l_input//" and "//l_refimg//" are "//\
                                "not on the same side of 841449600 "//\
                                "(31-08-2006). Exiting.", \
                                l_logfile, verbose+)
                        } else {
                            # Only needed for EEV detector and
                            # DETTYPE SDSU II CCD
                            # Set some flags for later on
                            # stopdatekey is used to stop unwanted lines when
                            # parsing QE data - default is old
                            if ((inobsdate < gmossbiasdate) && \
                                (indettype == det_names[1])) {
                                    stopdatekey = "newer"
                                    #Check for original GMOS-S CCDs
                                    # Not supported
                                    if (inobsdate < gmossccddate) {
                                        printlog ("GQECOR ERROR - GMOS-S "//\
                                            "observations before 2003-04-01"//\
                                            " are not supported. Exiting",
                                            l_logfile, verbose+)
                                    }
                            }
                        }
                    } else {
                        printlog ("GQECORR ERROR - Keyword "//\
                            obskey_to_useref//\
                            " found in file: "//l_refimg//\
                            " but it is an empty string. Exiting.", \
                            l_logfile, verbose+)
                        goto clean
                    }
                } else if (ininstrument == "GMOS-S") {
                    printlog ("GQECORR ERROR - Keyword "//l_key_date1//\
                        " or "//l_key_date2//" found in file: "//l_input//\
                        " but it is an empty "//\
                        "string. Exiting.", l_logfile, verbose+)
                    goto clean
                } # End of OBSDATE checks

                # Check bining against l_refimage
                if ((xinbin[ii] != xrefbin[kk]) || \
                    (yinbin[ii] != yrefbin[kk])) {
                        printlog ("GQECORR ERROR - Binning mismatch between"//\
                            " "//l_input//" and "//l_refimg//". Exiting", \
                            l_logfile, verbose+)
                    goto clean
                } # End of binning check against refimg

                # Check for the existance of N&S NODPIX keyword
                # but don't check against refimg as they will not exist in them
                imgshift = 0
                keypar (l_input//"[0]", l_key_noddistance, silent+)
                if (!keypar.found) {
                    # Check for the existance of some of the other N&S keywords
                    #     to double check that l_key_noddistance isn't just
                    #     missing

                    # Check for ANODCNT
                    keypar (l_input//"[0]", l_key_anodcount, silent+)
                    if (!keypar.found) {
                        # Check for BNODCNT
                        keypar (l_input//"[0]", l_key_anodcount, silent+)
                        if (keypar.found) {
                            # N&S keywords exist but some missing
                            printlog ("GQECORR ERROR - There are Nod-and-"//\
                                "Shuffle keywords present in "//l_input//\
                                " but the required "//l_key_noddistance//\
                                " is not present in the phu. Exiting.",
                                l_logfile, verbose+)
                            goto clean
                        }
                    } else {
                        # N&S keywords exist but some missing
                        printlog ("GQECORR ERROR - There are Nod-and-"//\
                            "Shuffle keywords present in "//l_input//\
                            " but the required "//l_key_noddistance//\
                            " is not present in the phu. Exiting.",
                            l_logfile, verbose+)
                            goto clean
                    }

                    if (refim_gqenod[kk]) {
                        printlog ("GQECORR ERROR - Reference/correction"//\
                            "file "//l_refimg//" has been corrected\n"//\
                            "                for nod-and-shuffle but "//\
                            l_input//" has not been shuffled.\n"//\
                            "                A new reference or correction "//\
                            "image is required.\n"//\
                            "                Exiting.", \
                            l_logfile, verbose+)
                        goto clean
                    }
                } else if (corrim_exists[kk]) {

                    # Assign imgshift
                    imgshift = int(keypar.value)

                    # Counter for the number of N&S data files input
                    num_ns_infiles += 1

                    # Printlog saying input is N&S
                    printlog ("GQECORR - Nod-and-Shuffle data found. "//\
                        l_input//" has a shift of "//imgshift//" pixels.", \
                        l_logfile, l_verbose)

                    # If using a GQECORR correction previously made check
                    # if GQECORR applied the nod and shuffle offest to the
                    # correction image
                    # If using a GSAPPWAV solution, no need to apply correction
                    # for GMOS-S N&S IFU data. Only needs to be done for a full
                    # wavelength olution
                    if (refim_gqenod[kk]) {
                        # To get to here NODPIX should exist
                        keypar (l_refimg//"[0]", l_key_noddistance, silent+)
                        # Check the NODPIX values match
                        if (int(keypar.value) != imgshift) {
                            printlog ("GQECORR ERROR - "//l_key_noddistance//\
                                " values do not match between "//l_input//\
                                " and "//l_refimg//". Exiting.", \
                                l_logfile, verbose+)
                            goto clean
                        }

                        # Counter for the number of pre-existing
                        # N&S corrimages
                        num_ns_corrim_exists += 1
                    } else {
                        # Check if the ns prefixed l_refimg exists
                        if (imaccess(ns_prefix//l_refimg)) {

                            # Double check for the l_key_gqenod keywprd
                            keypar (ns_prefix//l_refimg//"[0]", l_key_gqenod, \
                                silent+)
                            if (keypar.found) {
                                printlog ("GQECORR WARNING - A nod-and-"//\
                                    "shuffled "//\
                                    "version of "//l_refimg//" has been "//\
                                    "found on disk \n                  "//\
                                    "(prefixed with 'ns'). "//\
                                    "Using 'ns' version instead.", \
                                    l_logfile, verbose+)
                                # Reset refimages name
                                refimage_files[kk] = ns_prefix//l_refimg
                                # Reset the flag saying the existing corrim
                                # has been shuffled
                                refim_gqenod[kk] = yes
                                # Counter for the number of pre-existing
                                # N&S corrimages
                                num_ns_corrim_exists += 1
                            }
                        }
                    }
                }

            # End of input against refimg checks for spectroscopic data

            } else if (thisobsmode == "IMAGE") {
                #######
                # Imaging data

                lastfilter1 = ""
                lastfilter2 = ""

                # If more than one image store last image's FILTER keywords
                if (ii > 1) {
                    lastfilter1 = thisfilter1
                    lastfilter2 = thisfilter2
                }

                # Check FILTER1 and FILTER2 have entries
                # FILTER1
                keypar (l_input//"[0]", l_key_filter1, silent+)
                if (keypar.found) {
                    thisfilter1 = keypar.value
                    if (thisfilter1 == "" || (stridx(" ", thisfilter1) > 0)) {\
                        printlog ("GQECORR ERROR - "//l_key_filter1//\
                            " keyword has no value in image"//l_input//\
                            ". Exiting", l_logfile, verbose+)
                        goto clean
                    }
               } else {
                    printlog ("GQECORR ERROR - "//l_key_filter1//\
                        " keyword not found in image "//l_input//". Exiting", \
                        l_logfile, verbose+)
                    goto clean
                }

                # Make the open filter1 string equal to open
                if (substr(thisfilter1,1,4) == "open") {
                    thisfilter1 = "open"
                }

                # FILTER2
                keypar (l_input//"[0]", l_key_filter2, silent+)
                if (keypar.found) {
                    thisfilter2 = keypar.value
                    if (thisfilter2 == "" || (stridx(" ", thisfilter2) > 0)) {\
                        printlog ("GQECORR ERROR - "//l_key_filter1//\
                            " keyword has no value in image"//l_input//\
                            ". Exiting", l_logfile, verbose+)
                        goto clean
                    }
                } else {
                    printlog ("GQECORR ERROR - "//l_key_filter2//\
                        " keyword not found in image "//l_input, \
                        l_logfile, verbose+)
                    goto clean
                }

                # Make the open filter2 string equal to open
                if (substr(thisfilter2,1,4) == "open") {
                    thisfilter2 = "open"
                }

                # If more than one image compare FILTER values to last image
                # Check incase the filter names got swapped in keywords
                if (ii > 1) {
                    if (((thisfilter1 != lastfilter1) || \
                            (thisfilter2 != lastfilter2)) && \
                                ((thisfilter1 != lastfilter2) || \
                                    (thisfilter2 != lastfilter1))) {

                         # If they don't match change matching filter to false
                         fl_filter_match = no
                     }
                }

                # Count the number of (bsmode == IMAGE) files
                count_imgobsmode += 1

            } # End of input imaging checks

        } else { # OBSMODE keyword not found in input image
            printlog ("GQECORR ERROR - "//l_key_obsmode//" keyword not "//\
                " found in input image"//l_input//". Exiting.", \
                l_logfile, l_verbose)
            goto clean
        } # End of specfic OBSMODE checks

    } # End of for ii loop checking specific requirements for INPUT IMAGES

    ##################
    # Checks to see if anything else is to be and to set the isimaging flag

    # Set the imagimg data only flag - checks for mixing of spectroscopic and
    # imaging data
    if (count_imgobsmode == num_inimages && count_imgobsmode != 0) {
        #Imaging data only
        isimaging = yes

        printlog ("GQECORR - All files are imaging data.", l_logfile, \
            l_verbose)

        # Check for userdefined (by hand) numbers
        if (userdefined) {

            # If relative correcion and not 3 numbers
            if (l_fl_relative_corr && user_qenum_count != 3) {
                # Exit
                printlog ("GQECORR ERROR - Relative correction requested "//\
                    "but user has not input 3 QE factors. Exiting", \
                    l_logfile, verbose+)
                goto clean
            # If absolute has to be one or three
            } else if (!l_fl_relative_corr && ((user_qenum_count != 3) || \
                (user_qenum_count != 1))) {
                    # Exit
                    printlog ("GQECORR ERROR - Absolute correction "//\
                        "requested but user has not input 1 or 3 QE factors"//\
                        ". Exiting", \
                        l_logfile, verbose+)
                    goto clean
                }
        }

    } else if (count_imgobsmode != 0) {
        # Mixed imaging and spectroscopic data
        printlog ("GQECORR ERROR - Imaging and spectroscopic input files "//\
            "supplied. Cannot mix imaging and \n"//\
            "                spectroscopic data in GQECORR. "//\
            "Different spectroscopic data can \n"//\
            "                be supplied, so long as each "//\
            "frame has a reference file with matching\n"//\
            "                OBSMODE supplied. Exiting.", \
            l_logfile, verbose+)
        goto clean
    } # End of checks for imaging only / mixed imaging and spectroscopic data

    # If not correcting the data check if anything else needs to be done
    # Initiate flag to escape
    alldone = no

    # Spectroscopic data
    # If all the ns correction images already exist exit or
    # If no N&S data present and all correction images exist exit
    if (!isimaging && (((num_corrim_exists == num_ns_corrim_exists || \
        num_corrim_exists == 1) && \
        (num_ns_corrim_exists != 0)) || \
            ((num_ns_infiles == 0) && ((num_corrim_exists == 1 || \
                num_corrim_exists == num_inimages))))) {
        alldone = yes
    }

    if (debug) {
        printlog ("num_corrim_exists: "//num_corrim_exists//\
            "\n"//"num_inimages: "//num_inimages//\
            "\n"//"num_ns_infiles: "//num_ns_infiles//\
            "\n"//"num_ns_corrim_exists: "//num_ns_corrim_exists, \
            l_logfile, verbose+)
        printlog ("Alldone is: "//alldone, l_logfile, verbose+)
    }

    # Imagaing
    if (isimaging) {
        # All the same filter
        if (fl_filter_match) {

            printlog ("          Matching filters were used.", \
                l_logfile, l_verbose)

            # Either one correction file exists or there are teh same
            # number as input files
            if ((num_corrim_exists == 1) || \
                (num_corrim_exists == num_inimages)) {

                    alldone = yes
            }
        } else {
            printlog ("          More than one set of filters were used.", \
                l_logfile, l_verbose)
            if (num_corrim_exists == num_inimages) {
                # More than one pair of filters used
                # All input images must have a correction file
                alldone = yes
            }
        }

    }

    # Set this here in case of nothing else to do and is being called
    # by gireduce.
    num_corrimgs_to_use = num_corrim_exists

    if (debug) {
        printlog ("Setting num_corrimgs_to_use to: "//num_corrimgs_to_use,\
            l_logfile, verbose+)
        printlog ("Alldone: "//alldone//" l_correct_images: "//\
            l_correct_inimages, l_logfile, verbose+)
    }

    if (alldone && !l_correct_inimages) {
        # No need to create correction data and no need to correct inputs
        printlog ("GQECORR - All the correction images already exist "//\
            "and fl_correct=no.\n"//\
            "          There is nothing left for GQECORR to do.\n"//\
            "          Exiting.", \
            l_logfile, verbose+)
        # Exit
        goto FILE_FOR_GIREDUCE
    } else if (alldone && l_correct_inimages) {
        # No need to create correction data but need to correct inputs
        # Cannot skip the next section of code as it does other checks that
        # are needed.
        printlog ("GQECORR - All the correction images already exist "//\
            "and fl_correct=yes.\n", l_logfile, l_verbose)
    } else if (!alldone) {
        # Still correction images to create
        printlog ("GQECORR - Creating any required/missing correction data.", \
            l_logfile, l_verbose)
    }

    ######################################
    # End of specific requirement checks #
    ######################################

    ##########################
    # Start of the hard work #
    ##########################

    # Input files are in inimage_files[]
    # Output files are in outimage_files[]
    # Refernce files are in refimage_files[] if supplied
    # Corrimage files are in corrimage_files[]

    # All check have been completed that should halt the script bar a couple
    #    that should liekly have been picked up by other tasks to get this far

    ################
    # Create the correction images (spectroscopy) or obtain correction factors
    # (imaging)

    # It will be a mosaiced image but then cut according to input images later

    # Do everything with respect to DETSEC to minimise the amount of checks
    #    later on in the correction loop

    # For imaging data, only one factor per CCD is applied to the data set -
    #     so, no correction images are created.

    # Create the correction image
    if (!isimaging) { # Spectroscopic data

        if (debug) {
            printlog ("num_refimages is: "//num_refimages, l_logfile, verbose+)
        }

        # Loop of reference images
        for (ii = 1; ii<= num_refimages; ii += 1) {

            # Set the refimg variable
            l_refimg = refimage_files[ii]//".fits"

            # This is used when fl_fullcorr is yes
            l_wav_solution = refimage_files[ii]

            printlog ("GQECORR - Refimage: "//l_refimg//\
                      "; Wavelength solution: "//l_wav_solution, \
                      l_logfile, l_verbose)

            # Set associted input file
            l_input = inimage_files[ii]

            printlog ("GQECORR - Input image: "//l_input, \
                      l_logfile, l_verbose)

            # Assign correction image output name
            l_corrimg = corrimage_files[ii]
            l_corrimg = l_corrimg//".fits"

            printlog ("GQECORR - Corrimg: "//l_corrimg, l_logfile, l_verbose)

            # Store l_corrimg again
            corrimage_files[ii] = l_corrimg

            # Read the number of science frames in the reference image
            keypar (l_refimg//"[0]", l_key_nsciext, silent+)
            num_ref_sci = int(keypar.value)

            # By this point the binning is the same for l_input and l_refimg
            xbin = xrefbin[ii]
            ybin = yrefbin[ii]

            # Initiate the cut LS data flag
            iscutLS = no

            # Set l_fl_untransform_corrim to user input request
            # may get reset later on...
            l_fl_untransform_corrim = fl_untransform_corrim

            # Initiate the N&S pixel shift in case required
            imgshift = 0

            # Read the OBSMODE
            keypar (l_refimg//"[0]", l_key_obsmode, silent+)
            thisobsmode = keypar.value
            # Check if LS then check for gscut
            if (thisobsmode == "LONGSLIT") {
                #Check if gscut
                keypar (l_refimg//"[0]", l_key_gscut, silent+)
                if (keypar.found) {
                    iscutLS = yes
                }
            } else {
                # Check for N&S nodpix value - MOS / IFU
                # This keyword only exists in data that has been N&S'd, so not
                # flats or arcs - need to use input
                keypar (l_input//"[0]", l_key_noddistance, silent+)
                if (keypar.found) {
                    # NOD distance is in unbinned pixels.
                    imgshift = int(keypar.value)

                    # Correct for binning
                    imgshift = int(keypar.value) / ybin
                }
            }

            corrim_skip_flag = no
            # Check if the correction image exists for this file exists
            if (corrim_exists[ii]) {
                # If N&S data but refimg has not been corrected for N&S
                # create a new correction image

                # If using a GSAPPWAV solution, no need to apply correction
                # for GMOS-S N&S IFU data. Only needs to be done for a full
                # wavelength solution - though due to imexpr call to apply
                # shifted image currently there is no need to check the OBSMODE
                # and requested correction - can just continue.
                if ((imgshift > 0) && (!refim_gqenod[ii])) {
                #&& (thisobsmode != "IFU" && !l_fl_fullcorr)) {

                    printlog ("GQECORR - Nod-and-shuffle input data.\n"//\
                        "          Supplied correction images have not "//\
                        "been adjusted for the shuffle.", \
                        l_logfile, verbose+)

                    # Prefix ns to the front of the correction image
                    l_corrimg = ns_prefix//l_corrimg

                    printlog ("          Creating "//l_corrimg//" that will"//\
                        " have the shuffle applied.", l_logfile, verbose+)

                    # Set flag to skip things and t odo N&S correction
                    corrim_skip_flag = yes

                    # Skip some of the checks
                    goto CORRIM_EXISTS_NS

                } else {
                    # No need to do anything
                    goto NEXT_REFIMAGE
                }
            }

            ########
            # Need to access the input images for sizes, extensions and CCD
            # names - to get to here there are the same number of refimages
            # as there are inout images or there is only one ref image
            # and they are all the same dettype, detetctor, same mode etc.
            # So OK to do this here - it should have failed by now if not.

            # Obtain the number of science extensions - checked already
            keypar (l_input//"[0]", l_key_nsciext, silent+)
            num_sci = int(keypar.value)

            # Acess DETTYPE and any other keywords needed to ditinguish
            # between detectors with different pixel dimensions etc.
            # DETTYPE already checked earlier on
            keypar (l_input//"[0]", l_key_dettype, silent+)
            thisdettype = keypar.value

            # DETECTOR - already checked
            keypar (l_input//"[0]", l_key_detector, silent+)
            thisdetector = keypar.value

            # INSTRUMENT - already checked
            keypar (l_input//"[0]", l_key_instrument, silent+)
            thisinstrument = keypar.value

            # Read the chip gaps from stored information at the start of
            # the script. num_dettypes hardcoded at top too
            for (kk = 1 ; kk <= num_dettypes; kk += 1) {
                stored_dettype = det_names[kk]
                if (stored_dettype == thisdettype) {
                    if (xbin == 1) {
                        chipgap = det_info[kk,1] # Unbinned chipgap
                    } else {
                        chipgap = det_info[kk,2] # Binned chipgap
                    }
                    break
                }
            } # End of kk loop of dettypes to obtain chip gaps

            # Read the dimensions (DATASEC/CCDSEC/DETSEC) of the input file
            # Allows for absolute CCD2 only corrections plus other setups that
            # may not be common. Assumes a 1-D row of 2-D CCDs

            # Loop over science extensions of input image
            for (jj = 1; jj <= num_sci; jj += 1) {

                # Set SCI extension to read
                thisextn = "["//l_sciext//jj//"]"

                # Access DATASEC keyword and update fullxdim
                keypar (l_input//thisextn, l_key_datasec, silent+)
                # Parse data section values
                # Keep for correcting data later
                datasecin[ii,jj] = keypar.value
                junk = fscanf(keypar.value, "[%d:%d,%d:%d]", dummyx1, \
                    dummyx2, dummyy1, dummyy2)

                datax1[jj] = dummyx1
                datax2[jj] = dummyx2
                datay1[jj] = dummyy1
                datay2[jj] = dummyy2

                # Access DETSEC keyword and update fullxdim
                keypar (l_input//thisextn, l_key_detsec, silent+)
                # Parse data section values
                detsecin[jj] = keypar.value
                junk = fscanf(keypar.value, "[%d:%d,%d:%d]", dummyx1, \
                    dummyx2, dummyy1, dummyy2)

                detinx1[jj] = dummyx1
                detinx2[jj] = dummyx2
                detiny1[jj] = dummyy1
                detiny2[jj] = dummyy2

                # Access CCDSEC keyword and update fullxdim
                keypar (l_input//thisextn, l_key_ccdsec, silent+)
                # Parse data section values
                ccdsecin[jj] = keypar.value
                junk = fscanf(keypar.value, "[%d:%d,%d:%d]", dummyx1, \
                    dummyx2, dummyy1, dummyy2)

                ccdsecx1[jj] = dummyx1
                ccdsecx2[jj] = dummyx2
                ccdsecy1[jj] = dummyy1
                ccdsecy2[jj] = dummyy2

                # Initiate the ampsorder array
                ampsorder[jj] = jj

            } # End of jj loop over input SCI extensions to store DETSEC etc.

            # Get the amplifiers into the correct order
            for (j=2; j<= num_sci; j = j+1) {
                for (jj=1; jj<j; jj=jj+1) {
                    if (detinx1[j] < detinx1[jj]) {
                        tmpint = ampsorder[jj]
                        ampsorder[jj] = ampsorder[j]
                        ampsorder[j] = tmpint
                    }
                }
            } # End of j and jj loop to determine amps order

            # Find number of CCDs used by looping of the amps and
            # comparing detsec and ccdsec also find number_of_ampsperccd
            # Inititate them - there will always be 1 ccd and 1 amp at least
            # by this point. Also, initiate mosaiced correction output
            # dimensions.
            fullxdim = 0
            fullydim = 0
            num_ccds = 1
            number_of_ampsperccd[num_ccds] = 1

            # Loop over ordered amplifiers
            for (j = 2; j <= num_sci; j += 1) {

                # Initiate ccdwidth (binned) and ampheight and CCDNAME using
                # j - 1 = 1 (the first amplifier)
                if (j==2) {
                    # Calaculate the ampwidth
                    ccdwidth[num_ccds] = (datax2[ampsorder[j-1]] - \
                        (datax1[ampsorder[j-1]] - 1))

                    # Calaculate the first amp's width
                    ampheight = (datax2[ampsorder[j-1]] - \
                        (datax1[ampsorder[j-1]] - 1))
                    fullxdim = ampheight

                    # Calaculate the first amp's height
                    ampheight = (datay2[ampsorder[j-1]] - \
                        (datay1[ampsorder[j-1]] - 1))
                    fullydim = ampheight

                    # Read CCDNAME
                    keypar (l_input//"["//l_sciext//str(ampsorder[j-1])//"]", \
                        l_key_ccdname, silent+)
                    ccdname[num_ccds] = str(keypar.value)
                }

                # Uses largest y dimension - only valid for 1row of 2D arrays
                if ((datay2[ampsorder[j]] - (datay1[ampsorder[j]] - 1)) > \
                    fullydim) {
                        ampheight = datay2[ampsorder[j]] - \
                            (datay1[ampsorder[j]] - 1)
                        fullydim = ampheight
                }

                # Calaculate the current ampwidth
                ampwidth = (datax2[ampsorder[j]] - \
                    (datax1[ampsorder[j]] - 1))

                # Check for a new CCD
                if (detinx1[ampsorder[j]] > detinx1[ampsorder[j-1]] && \
                    ccdsecx1[ampsorder[j]] <= ccdsecx1[ampsorder[j-1]]) {
                    # New CCD found

                    # Set extension to read - for CCDNAME
                    thisextn = "["//l_sciext//jj//"]"

                    # Increment the number of ccds
                    num_ccds += 1
                    number_of_ampsperccd[num_ccds] = 1

                    # Initialise this CCD's width
                    ccdwidth[num_ccds] = ampwidth

                    # Add the chip gap to fullxdim
                    fullxdim += (chipgap / xbin)

                    # Read CCDNAME
                    keypar (l_input//"["//l_sciext//str(ampsorder[j])//"]", \
                        l_key_ccdname, silent+)
                    ccdname[num_ccds] = str(keypar.value)

                    # Check the ccdname has actually changed
                    if (ccdname[num_ccds] == ccdname[num_ccds-1]) {
                        printlog ("GQECORR ERROR - When counting amplifiers"//\
                           " and CCDs, a new CCD was counted but the "//\
                           l_key_ccdname//" keyword did not change. "//\
                           "Exiting.", l_logfile, verbose+)
                        goto clean
                    }

                } else {
                    # Not a new CCD

                    # Increment the number of amplifiers per CCD
                    number_of_ampsperccd[num_ccds] = \
                        number_of_ampsperccd[num_ccds] + 1

                    # Calculate the ccdwidth
                    ccdwidth[num_ccds] = ccdwidth[num_ccds] + ampwidth
                }

                # Update fullxdim
                fullxdim += ampwidth

            } # End of j loop counting num_ccds, aamps per CCD, full(x/y)dim

            ############
            # Some very final checks on the user requested parameters and
            # the input data

            # Check if only CCD2 data - assumes only 1 CCD means CCD2
            if (num_ccds == 1) {
                # Check if trying to perform a relative correction
                if (l_fl_relative_corr) {
                    # If so fail
                    printlog ("GQECORR ERROR - Cannot apply a relative "//\
                        "correction to CCD2 data only. Exiting.",
                        l_logfile, verbose+)
                    goto clean
                }

                # Check if trying to unmosaic only CCD2 data
                if (l_fl_untransform_corrim) {
                    # Reset unmosaic flag to no
                    l_fl_untransform_corrim = no
                    printlog ("QECORR WARNING - cannot untransform CCD2 "//\
                        "only data. Swithing fl_untransform_corrim to no "//\
                        "for "//l_input//".", l_logfile, verbose+)
                }
            } # End of CCD2 only checks

            ## Check that the number of CCDs worth of coeefficents supplied
            ## by user (if done) are the same as the number of CCDs
            ## if (user_qenum_count != num_ccds)

            #############
            # Creation of the correction image

            # Set tmp correction image name
            tmpcorrim = mktemp("tmpcorrim")//".fits"

            if (debug) {
                printlog ("GQECORR - Dimensions of the input file mosaiced "//\
                    "are: x="//fullxdim//" y="//fullydim, \
                    l_logfile, l_verbose)
            }

            # Create the temporary correction image here to the fullxdim and
            # fullydim determined above (binned) (all pixels have a value of 1)
            mkpattern (tmpcorrim, output="", option="replace", \
                pattern="constant", v1=0.0, v2=0, pixtype=imgpixtype, \
                title="QE Correction Image created using "//l_refimg, ndim=2, \
                ncols=fullxdim, nlines=fullydim, header="")

CORRIM_EXISTS_NS:

            # Initiate the output correction file
            # Do here to allow for pre-made correction images to be shuffled
            # Copy the PHU of the input to the output correction file
            # Use input if creating image from scratch otherwise use refimg
            if (!corrim_skip_flag) {
                imcopy (l_input//"[0]", l_corrimg, verbose=debug)
            } else {
                imcopy (l_refimg//"[0]", l_corrimg, verbose=debug)
            }

            # Copy refimg MDF to output too
            # Send the output to dev$null because longslit data may not have
            # a MDF. Always use the refimg MDF
            ttools.tcopy (l_refimg//"["//l_mdf_ext//"]", \
                l_corrimg//"["//l_mdf_ext//"]", verbose=no, \
                >& "dev$null")

            # Now loop over science extensions of reference image
            # to start to populate the correction image with pixels, who's
            # values are set to the wavelength of that pixel
            for (kk = 1; kk <= num_ref_sci; kk += 1) {

                # Set extension to read
                thisextn = "["//l_sciext//kk//"]"
                wavextn = "["//wav_extname//","//kk//"]"
                thisextnappend = "["//l_sciext//kk//",append]"

                printlog ("GQECORR - Working on refimg extension: "//thisextn, \
                          l_logfile, l_verbose)

                # Correction image exists just need to correct for N&S
                if (corrim_skip_flag) {
                    # Set correct file name
                    somefilename = l_refimg//thisextn

                    # Skip some more things
                    goto CORRIM_EXISTS_APPLY_NS
                }

                # Need to define the region to replace in the correction image
                # with the created the lambda image, for the given science
                # frame. Read the detector section if iscutLS is false. If
                # true, set it to fullxdim and fullydim
                # Applies to already created refimages being shuffled
                if ((!iscutLS && thisobsmode == "LONGSLIT") || \
                    corrim_skip_flag) {
                    xstart = 1
                    xstop  = fullxdim
                    ystart = 1
                    ystop  = fullydim
                } else {
                    # Read DETSEC of reference SCI extension
                    keypar (l_refimg//thisextn, l_key_detsec, silent+)
                    # Parse det section values
                    junk = fscanf(keypar.value, "[%d:%d,%d:%d]", detx1, \
                        detx2, dety1, dety2)

                    # Adjust for binning - now in reference of the output
                    # correcion image (chip gaps are already taken care of
                    # by GSCUT)
                    xstart = nint(((detx1 - 1) / xbin) + 1)
                    xstop  = nint(detx2 / xbin)
                    ystart = nint(((dety1 - 1) / ybin) + 1)
                    ystop  = nint(dety2 / ybin)

                }

                #############
                # Start the creation of the wavelength lamba imagine for this
                # extnsion - the pixel values will be that of the wavelength of
                # the pixel

                # Set dimensions for width and height of images
                xmax = xstop - (xstart - 1)
                ymax = ystop - (ystart - 1)

                # Initiate temporary files
                #########
                # These few have been left in, in case of implementing full
                # correction
                # Text file to contain the list of pixels&values from ref_image
                tmppixels = mktemp("tmppixels")
                # Text file to just contain pixel number
                #tmponlypix = mktemp("tmponlypix")
                #########
                # FITS image of tmppixels
                tmppiximg = mktemp("tmppiximg")//".fits"
                # FITS image of correct lambda of tmmpiximg
                tmp1dimg = mktemp("tmp1dimg")//".fits"
                # FITS image of 2D blkrep imgae of tmp1dimg
                tmp2dimg = mktemp("tmp2dimg")//".fits"
                # Temporay output file of SCI extension - to be replaced in
                # actual output
                tmpoutimg = mktemp("tmpoutimg")

                # Type of wavelength solution to use
                if (!l_fl_fullcorr) { # GSAPPWAVE solution

                    # Read WCS
                    # REFPIX - in reference of datasec
                    keypar (l_refimg//thisextn, l_key_refpix, silent+)
                    if (keypar.found) {
                        refpix = real(keypar.value)
                    } else {
                        # CRPIX1 - in reference of datasec
                        keypar (l_refimg//thisextn, l_key_crpix1, silent+)
                        if (keypar.found) {
                            refpix = real(keypar.value)
                        } else {
                            printlog ("GQECORR ERROR - No WCS "//\
                                l_key_refpix//" or "//l_key_crpix1//\
                                " keywords found in "//l_refimg//thisextn//\
                                " . Exiting.", l_logfile, verbose+)
                            goto clean
                        }
                    }

                    # CRVAL1
                    keypar (l_refimg//thisextn, l_key_crval1, silent+)
                    if (keypar.found) {
                        refval = real(keypar.value )
                    } else {
                        printlog ("GQECORR ERROR - No WCS "//\
                            l_key_cdval1//" keyword found in "//\
                            l_refimg//thisextn//\
                            ". Exiting.", l_logfile, verbose+)
                        goto clean
                    }

                    # CD1_1
                    keypar (l_refimg//thisextn, l_key_cd1_1, silent+)
                    if (keypar.found) {
                        refxdelta = real(keypar.value)
                    } else {
                        printlog ("GQECORR ERROR - No WCS "//\
                            l_key_cd1_1//" keyword found in "//\
                            l_refimg//thisextn//\
                            ". Exiting.", l_logfile, verbose+)
                        goto clean
                    }

                    # Initate tilted slit flag
                    istilted = no
                    tilt = 0

                    # If MOS check MDF for tilted slits (this check indludes
                    # NS LS files)
                    if (thisobsmode == "MOS") {

                        # Read the mask name
                        keypar (l_refimg//"[0]", l_key_maskname, silent+)
                        if (keypar.found) {
                            # To catch NS LS data
                            if (substr(keypar.value,1,2) == "NS") {
                                column_to_read = l_tabkey_slittilt_LS

                            # For NS LS data this will return a tilt of 90
                            # degrees due to the definition of slittilt_m
                            # no need to check if it's tilted, it's LS data
                            goto NS_LS_SKIP
                            } else {
                                # All MOS data
                                column_to_read = l_tabkey_slittilt
                            }
                        } else {
                            printlog ("GQECORR ERROR - keyword "//\
                                l_key_maskname//" is not present in "//\
                                l_refimg//". Exiting", \
                                l_logfile, verbose+)
                            goto clean
                        }

                        # Obtain the MDF row number for this SCI extension
                        keypar (l_refimg//thisextn, l_key_mdfrow, silent+)
                        if (keypar.found) {
                            row_to_read = int(keypar.value)
                            if (debug) {
                                printlog ("GQECORR - MDF row for "//\
                                    l_refimg//thisextn//" is: "//\
                                    row_to_read//".", \
                                    l_logfile, l_verbose)
                            }
                        } else {
                            printlog ("GQECORR ERROR - keyword "//\
                                l_key_mdfrow//" is not present in "//\
                                l_refimg//thisextn//". Exiting", \
                                l_logfile, verbose+)
                            goto clean
                        }

                        # Read the slittilt (degrees; relative to MASK_PA)
                        tprint (l_refimg//"["//l_mdf_ext//"]", prparam=no, \
                            prdata=yes, \
                            showrow=no, showhdr=no, showunits=no, \
                            col=column_to_read, rows=row_to_read, \
                            pwidth=80) | scan (tilt)

                        # Tilted slit
                        if (tilt != 0.0) {
                            istilted = yes
                            printlog ("GQECORR - Tilted slit found. "//\
                                l_refimg//thisextn//" has a tilt of "//tilt//\
                                " degrees.", l_logfile, l_verbose)

                            # Have to calculate change in wavelength w.r.t. y
                            # Crude way is to use MASK_PA, slittilt, CD1_1
                            # and simple trig to calculate the cange in
                            # wavelength w.r.t. y.

                            maskpa = 0

                            # Read the MASK_PA from MDF [in degrees]
                            # Only in MOS masks not in package masks
                            keypar (l_refimg//"["//l_mdf_ext//"]", \
                                l_key_maskpa, silent+)
                            if (keypar.found) {
                                maskpa = real(keypar.value)
                            } else {
                                printlog ("GQECORR ERROR - "//l_key_maskpa//\
                                    " missing in "//l_refimg//"["//\
                                    l_mdf_ext//"]", \
                                    l_logfile, verbose+)
                                goto clean
                            }

                            # Need to convert degrees to radians to use
                            # tangent in iraf - define PI
                            pi = (asin(1) * 2.0)

                            # Calculate the angle in radians
                            ang_to_use = (maskpa + tilt) * (pi / 180.0)

                            # refydelta is the change in wavelength for a 1
                            # pixel shift in positive y
                            refydelta = (refxdelta / tan(ang_to_use))

                        } # End of MOS tilted slit calculations

                    } # End of MOS tilted slits check

NS_LS_SKIP:

                    # Make a 1D image who's pixel values are the same as their
                    # X coordinate - in DATASEC frame of reference of SCI extn
                    mkpattern (input=tmppiximg, output="", \
                        pattern="coordinates", option="replace", \
                        pixtype="real", ndim=2, ncols=xmax, nlines=1, \
                        header="")

                    # Use imexpression to calculate the correct wavelength
                    # for each pixel
                    # lambda_expr = ((imagepixel - refpix) * refxdelta) + \
                    #     refval
                    lambda_expr = "(((a - b) * c) + d)"

                    # Don't specify refim parameter in imexpr it gets confused
                    # refim="auto"
                    imexpr (lambda_expr, tmp1dimg, tmppiximg, refpix, \
                        refxdelta, refval, dims="auto", intype="auto", \
                        outtype="auto", bwidth=0, \
                        btype="nearest", bpixval=0, verbose=no)

                    # Block replicate the image to the desired size
                    blkrep (tmp1dimg, tmp2dimg, b1=1, b2=ymax)

                    somefilename = tmp2dimg

                    # Correct the wavelength for tilted slits
                    if (istilted) {
                        # Is a tilted slit

                        # make a temporary file again
                        tmp2dimgtilt = mktemp("tmp2dimgtilt")

                        # Use imexpr to do the correction
                        # Assume the middle of the slit has the most
                        # accurate wavelength form GSAPPWAVE
                        # ycentre needs to be a non-int value hence the 1.0
                        ycentre = 1.0 * ymax / 2
                        # Using IRAFs way of rounding down when casting real to
                        # int, to test for even values later
                        ycen_int = int(ycentre)

                        # Will need to adjust for even number of pixels in y
                        # Set wav_adjustment accordingly
                        wav_adjustment = 0
                        if ((ycentre - ycen_int) == 0.0) {
                            wav_adjustment = 1
                        }

                        # J is the imexpression y coordinate of input
                        # a is the blkrep image
                        lambda_expr = "J == c ? a + (b * 0.5 * d) : "//\
                            "a + ((J - c + (b * 0.5)) * d)"

                        # Don't specify refim parameter in imexpr it gets
                        # confused - refim="auto"
                        imexpr (lambda_expr, tmp2dimgtilt, tmp2dimg, \
                            wav_adjustment, ycent_int, \
                            refydelta, dims="auto", intype="auto", \
                            outtype="auto", bwidth=0, \
                            btype="nearest", bpixval=0, verbose=no)

                        somefilename = tmp2dimgtilt
                    }

                } else { # Full identify or fitscoord solution
                    currentarea = ""

                    # Create text file containing coordinates from input image
                    tmpcoords = mktemp ("tmpcoords")

                    # Format the extension number for accessing the database
                    if (thisobsmode == "MOS") {
                        # Obtain the MDF row number for this SCI extension
                        keypar (l_refimg//thisextn, l_key_mdfrow, silent+)
                        printf ("%03d", int(keypar.value)) | \
                            scan (ext_formated)

                    } else if (thisobsmode == "LONGSLIT") {
                        printf ("%03d", kk) | scan (ext_formated)
                    }

                    if (thisobsmode == "MOS" || thisobsmode == "LONGSLIT") {
                        fparse (l_wav_solution)
                        lambda_fname = fparse.root//"_"//ext_formated
                        l_dbfile = l_database//"fc"//lambda_fname
                        if (!access(l_dbfile)) {
                            printlog ("ERROR - GQECORR: Cannot access "//\
                                      "fitcoords database: "//l_dbfile, \
                                      l_logfile, verbose+)
                            goto clean
                        }

                        # List the pixels and determine their wavelength using
                        # fceval
                        ##M Check the wcs used?
                        ##M Do ARCs need to be cut for LS?
                        ##M Are the WCS correct?
                        ##M Assume the entire image has been fitcoord -
                        ##M MOS (N&S?)/LS

                        listpix (l_refimg//thisextn, wcs="logical", \
                            formats="", verbose=no, >> tmppixels)

                        fceval (tmppixels, tmpcoords, fitnames=lambda_fname, \
                                database=l_dbfile)

                        tmp2coords = mktemp("tmp2coords")
                        fields (tmpcoords, fields="3", >> tmp2coords)

                        # Create the wavelength image
                        rtextimage (tmp2coords, output=tmp2dimg, otype="real", \
                            header=no, pixels=yes, nskip=0, dim=xmax//","//ymax)

                        delete (tmppixels//","//tmpcoords//","//tmp2coords, \
                                verify-, >& "dev$null")
                        somefilename = tmp2dimg

                    } else if (thisobsmode == "IFU") {
                        # Call gfunexwl SPP code
                        if (kk == 1) {
                            tmpifu = mktemp ("tmpifu_wav_image")//".fits"

                            # Needs only be called once - handle naming
                            # afterwards
                            gfunexwl (l_wav_solution, tmpifu, ncols=11, \
                                      xorder=4, yorder=3, database=l_database,
                                      verbose=debug)

                            # No exit status for gfunexwl - may hard crash
                            if (!imaccess(tmpifu)) {
                                printlog ("ERROR - GQECORR: Cannot access "//\
                                          "output from GFUNEXWL", l_logfile,
                                          verbose+)
                                goto clean
                            }
                        }

                        # Allow user to determine which slit takes precedence
                        # in the IFU-2 mode

                        # This syntax assume the red slit is in extension 1 and
                        # the blue in extension 2

                        ##M Maybe change for r / b syntax
                        ##M Overlaping slits not tested thoroughly

                        # Assumes 1 is default and at most two slits
                        cut_wavextn = ""
                        if (ifuslit_precedence == 1) {
                            if (kk == 1) {
                                # Requires to know information from this slit
                                # to adjust position and size of next slit
                                first_stop = xstop
                            } else {
                                # Adjust insert position and extension string
                                # to slice image
                                if (xstart <= first_stop) {
                                    # In the detecor reference
                                    overlap = first_stop - (xstart - 1)
                                    xstart += overlap

                                    # In the second slit reference frame
                                    cut_wavextn = "["//(overlap + 1)//":"//\
                                                    xmax//",*]"
                                } # Else no overlap - nothing to do
                            }
                        } # Else just insert the second slit over the first

                        imcopy (tmpifu//wavextn//cut_wavextn, tmp2dimg, \
                                verbose=debug)
                        somefilename = tmp2dimg
                        if (kk == num_ref_sci) {
                            imdelete (tmpifu, verify-, >& "dev$null")
                        }
                    }
                }

CORRIM_EXISTS_APPLY_NS:

                # Check if non-IFU N&S - and not using a non-shuffled existing
                # correction image
                if ((thisobsmode != "IFU") && (imgshift != 0) && \
                    (!corrim_skip_flag)) {

                    # Adjust ystart - Only this needs to be changed due to the
                    # way iminsert works (it only requires the bottom left
                    # corner of replacement region in image to be inserted into
                    if (imgshift > 0) { # sky below object (-ve means above)
                        ystart -= imgshift

                        printlog ("GQECORR - Applying Nod and Shuffle "//\
                            "offset for "//l_refimg//thisextn, \
                            l_logfile, l_verbose)
                    }

                    # Blockreplicate the image to double its size in y
                    blkrep (somefilename, tmpoutimg, b1=1, b2=2)

                    nextsomefilename = tmpoutimg

                # If N&S IFU or existing correction image which needs shuffling
                } else if ((thisobsmode == "IFU" && imgshift != 0) || \
                     corrim_skip_flag) {

                     # N&S IFU GMOS South only
                     # Use imshift to correct for N&S using imgshift

                     # Make a temporary file names
                     tmpimgshift = mktemp("tmpimgshift")
                     tmpimg2shift = mktemp("tmpimg2shift")

                     # Copy somefilename to tmpimgshift
                     imcopy (somefilename, tmpimgshift, verbose=debug)

                     # Shift tmpimgshift by imgshift
                     imshift (input=tmpimgshift, output=tmpimg2shift, \
                         xshift=0, yshift=-imgshift, shifts_file="", \
                         interp_type="nearest", boundary_type="constant", \
                         constant=1)

                     # Create temporay filename for nextsomefilename
                     nextsomefilename = mktemp("tmpnextsomefilename")

                     # Set up imexpr expression for adding shifted image
                     im_expr_shift = "a == 1 ? (a + b - 1) : a"

                     # Add the somefilename and tmpimg2shift together where
                     # appropriatete.
                     # Don't specify refim parameter in imexpr it gets confused
                     # refim="auto"
                     imexpr (im_expr_shift, nextsomefilename, \
                         somefilename, tmpimg2shift, dims="auto", \
                         intype="auto", \
                         outtype=imgpixtype, bwidth=0, \
                         btype="nearest", bpixval=0, verbose=no)

                     # Do not need tmpimg?shift anymore
                     imdelete (tmpimgshift//", "//tmpimg2shift, verify-, \
                         >& "dev$null")

                     # if working on a pre-made correction image
                     if (corrim_skip_flag) {
                         # Copy nextsomefilename to correction image
                         imcopy (nextsomefilename, l_corrimg//thisextnappend, \
                             verbose=debug)

                         # Skip to next extension
                         goto NEXT_REF_EXTENSION
                     }

                } else {
                    # Standard spectroscopy no N&S
                    # Assign next part's input name
                    nextsomefilename = somefilename

                } # End of N&S adjustments

                if (debug) {
                    fxhead (nextsomefilename)
                }

                # Write coordinates of where the bottom left corner of the
                # current correction SCI extension being created
                # will go in the temporray mosaiced correction file, to a
                # text file
                iminscoords = mktemp("tmpiminscoords")
                print (xstart//" "//ystart, > iminscoords)

                # Insert this science extension into temporary correction image
                iminsert (tmpcorrim, nextsomefilename, tmpcorrim, "replace", \
                    coordfile=iminscoords, offset1=0, offset2=0, xcol="c1", \
                    ycol="c2")

                # Clean up these files
                delete (iminscoords, \
                    verify-, >& "dev$null")

                imdelete (tmppiximg//", "//tmp1dimg//", "//tmp2dimg//", "//\
                    tmpoutimg, verify-, >& "dev$null")
NEXT_REF_EXTENSION:

                imdelete (somefilename//", "//nextsomefilename, \
                    verify-, >& "dev$null")

            } # End of for kk loop over science extensions in reference image

            # If just shuffling a pre-made correction image - all is complete
            if (corrim_skip_flag) {
                # Skip the rest of the file creation
                goto UPDATE_CORRIM_PHU
            }

            #########
            # For a given CCD calculate the QE correction and apply it
            # Creates lasttmpcorrim - contains phu and 3 SCI extensions (CCDs)

            # Make another temporary file - the last one before the
            # correction image is finished
            lasttmpcorrim = mktemp("tmplasttmpcorrim")//".fits"

            # Initiate lasttmpcorrim by copy input phu to it
            imcopy (l_input//"[0]", lasttmpcorrim, verbose=debug)

            # Make a temporary text file
            tmpcoeffs = mktemp("tmpcoeffs")

            # Initiate thisccdstart and thisccdstop
            thisccdstart = 1
            thisccdstop = 0

            # Loop over the CCDs
            for (ll = 1; ll <= num_ccds; ll += 1) {

                # Assign the CCDNAME (ccdname is order due to ordering amps)
                thisccd = ccdname[ll]

                # Due to the extreme number of white space in the CCDNAMEs
                # the CCDNAMEs in the data file are stored with the spaces
                # replaced with underscores. So, need an _ version of CCDNAME
                print (thisccd) | translit ("STDIN", from_string=" ", \
                    to_string="_", delete=no, collapse=no) | \
                        scan(thisccd_data)

                # Set up extensions to be used in this loop
                thisextn = "["//ll//"]"
                thisextnappend = "["//ll//",append]"
                insertextn = "["//str(ll-1)//"]"

                # Set up the section of tmpcorrimg to correct
                if (ll > 1) {
                    # Only add the chip gap if more than one CCD
                    thisccdstart = thisccdstop + (chipgap / xbin)
                }
                thisccdstop = thisccdstart + (ccdwidth[ll] - 1)
                currentsection = "["//thisccdstart//":"//thisccdstop//\
                    ",*]"

                # Create temporary image to write correected CCD to
                tmpimgcorrected = mktemp("tmpimgcorrected")

                ########
                # Correct the current CCDs worth of information

                if (debug) {
                    printlog ("GQECORR - Current CCD is: "//thisccd_data//\
                        " and correction "//\
                        "type is: "//correction_type, l_logfile, verbose+)
                }

                # Apply corrections
                # Don't do anything to CCD2 if correction is relative
                if ((ll == 2) && (l_fl_relative_corr)) {

                    if (debug) {
                        printlog ("GQECORR - Correction is "//\
                            correction_type//" and current CCD is CCD2 - "//\
                            " not calculating QE"//\
                            " correction.", l_logfile, verbose+)
                    }

                    # This needs to be an extension of all 1's again!
                    mkpattern (tmpimgcorrected, output="", option="replace", \
                        pattern="constant", v1=1.0, v2=0, \
                        pixtype=imgpixtype, \
                        title="QE Correction Image created using "//\
                        l_refimg, ndim=2, \
                        ncols=ccdwidth[ll], nlines=fullydim, header="")

                } else {

                    # Use default parameters / user supplied file which
                    # is set out like the default one (the layout of the user's
                    # file is not checked!
                    if (!userdefined) {
                        # Set mode to parse
                        thismode = "SPEC"

                        if (debug) {
                            print (thismode)
                        }

                        # Parse the QE correction coefficients to tmptext file
                        # The fields are correct as off 2011-09-05
                        match (pattern=stopdatekey, files=l_qecorr_data, \
                            stop=yes, \
                            print_file_names=no) | \
                            match (pattern=thisinstrument, stop=no, \
                                print_file_names=no) | \
                            match (pattern=thisdettype, stop=no, \
                                print_file_names=no) | \
                            match (pattern=thisdetector//"#|", stop=no,\
                                print_file_names=no) | \
                            match (pattern=correction_type, stop=no, \
                                print_file_names=no) | \
                            match (pattern=thisccd_data, stop=no, \
                                print_file_names=no) | \
                            match (pattern=thismode, stop=no, \
                                print_file_names=no) | \
                            fields (fields="1,5,7", \
                                lines="1-", quit_if_missing=no, > tmpcoeffs)

                        count (tmpcoeffs) | scan (chk_lines)
                        if (chk_lines <= 0) {
                            printlog ("GQECORR ERROR - "//\
                                thisdetector//\
                                ": no coefficients found. "//
                                "Exiting ", \
                                l_logfile, verbose+)
                            goto clean
                        }

                        if (debug) {
                            type (tmpcoeffs)
                        }

                        # Scan in the coefficients creating the imexpr
                        # expression
                        scanfile = tmpcoeffs
                        calexpr = ""
                        imexpr_char = "a"

                        # Read the coeeficients
                        while (fscanf(scanfile, "%s %d %g", dataccdname, \
                            thisorder,thiscoeff) != EOF) {

                            # Remove the " "s and add _s for comparison
                            print (ccdname[ll]) | translit ("STDIN", \
                                 from_string=" ", to_string="_", delete=no, \
                                 collapse=no) | scan(tmpccdname)

                            if (debug) {
                                printlog (dataccdname//" "//tmpccdname, \
                                    l_logfile, verbose+)
                            }

                            # Double check the CCDNAMES
                            if (tmpccdname != dataccdname) {
                                printlog ("GQECORR ERROR - "//\
                                    l_key_ccdname//\
                                    " mismatch between "//l_input//" and "//\
                                    "information in "//l_qecorr_data//\
                                    ". Exiting ", \
                                    l_logfile, verbose+)
                                goto clean
                            }
                            calexpr = calexpr//"+(("//str(thiscoeff)//")*("//\
                                imexpr_char//"**("//str(thisorder)//")))"
                        } # End of while loop reading coefficients from a file

                        if (debug) {
                            print (calexpr)
                        }
                    } else {
                        # user supplied numbers on command line
                        # Assumes the user ordered the numbers according to
                        # CCD, i.e., 1-3, then by coefficient order, 0-n

                        # One line per CCD - assumed ot be in order 1-3
                        # and each line will be one string
                        fields (tmpqedata, fields="1", lines=ll, \
                            quit_if_missing=no, print_file_names=no) | \
                            translit (infile="STDIN", from_string=",", \
                            to_string="\n", delete=no, > tmpcoeffs)

                        # Unfortunately had to duplicate this code -
                        # leaving the it's use up to the peril of the user

                        # Scan in the coefficients creating the imexpr
                        # expression
                        scanfile = tmpcoeffs
                        calexpr = ""
                        imexpr_char = "a"
                        thisorder = 0

                        # Read the coeeficients
                        while (fscanf(scanfile, "%g", thiscoeff) != EOF) {

                            calexpr = calexpr//"+("//str(thiscoeff)//"*"//\
                                imexpr_char//"**"//str(thisorder)//")"

                            thisorder += 1

                        } # End of while loop reading coeffs from user input

                    } # End of reading of QE coefficients


                    # Strip the leading "+"
                    if (substr(calexpr,1,1) == "+") {
                        calexpr = substr(calexpr,2,strlen(calexpr))
                    }

                    # Only perform calculation on area of image that
                    # that are no longer 0 (a is the correction image)
                    calexpr = "a == 0.0 ? 1.0 : "//calexpr

                    if (debug) {
                        printlog ("GQECORR - Imexpr expression for QE "//\
                            "correction for CCD "//thisccd//" is: ", \
                            l_logfile, verbose+)
                        printlog ("  "//calexpr, l_logfile, verbose+)
                    }

                    # Apply the correction using imexpression
                    # Don't specify refim parameter in imexpr it gets confused
                    # refim="auto"
                    imexpr (calexpr, tmpimgcorrected, \
                        tmpcorrim//currentsection, dims="auto", \
                        intype="double", \
                        outtype=imgpixtype, bwidth=0, \
                        btype="nearest", bpixval=0, verbose=no)

                } # End of if statement to perform the correction

                # See if the images are to be unmosaiced - don't unmosaic CCD2
                # Instances of only CCD2 data being supplied and unmosaicing
                # being requested are caught earlier on
                if ((ll != 2) && (l_fl_untransform_corrim)) {
                    # Ungmosaic the images

                } else {
                    # insert the frame into the output image
                    fxinsert (input=tmpimgcorrected, \
                        output=lasttmpcorrim//insertextn, groups="", \
                        verbose=debug)
                } # End of unmosaicing/fxinserting into lastcorrection image

                # Update current lastcorrction extension's
                # CCDNAME, EXTNAME and EXTVER as used later on
                gemhedit (lasttmpcorrim//thisextn, l_key_ccdname, thisccd, "",\
                    delete-)
                gemhedit (lasttmpcorrim//thisextn, "EXTNAME", l_sci_ext, "",\
                    delete-)
                gemhedit (lasttmpcorrim//thisextn, "EXTVER", ll, "",\
                    delete-)

                # Clean up
                imdelete (tmpimgcorrected, \
                    verify-, >& "dev$null")
                delete (tmpcoeffs, verify-, >& "dev$null")
                scanfile = ""

            } # End of ll loop over CCDs to create lasttmpcorrim

            # Delete tmpcoeffs coeffficients
            delete (tmpcoeffs, verify-, >& "dev$null")
            imdelete (tmpcorrim, verify-, >& "dev$null")

            ##########
            # Create the cut up final output corerction image l_corrimg

            # Need to match the input extensions to a CCD and cut using CCDSEC

            # Loop over SCI extensions in input image
            for (kk = 1; kk <= num_sci; kk += 1) {

                # Assign some extnsion names to use throughout loop
                inextn = "["//l_sciext//kk//"]"
                outextnappend = "["//l_sciext//kk//",append]"
                outextn = "["//l_sciext//kk//"]"

                if (debug) {
                    printlog ("GQECORR - Looking at input file extension: "//\
                        l_input//inextn, l_logfile, verbose+)
                }

                # Read CCDNAME
                keypar (l_input//inextn, l_key_ccdname, silent+)
                inccdname = keypar.value

                if (debug) {
                    printlog ("inccdname is: "//inccdname, \
                        l_logfile, verbose+)
                }

                # Read CCDSEC from ccdsec arrays collected earlier
                # Correct for binning
                ccdx1 = ((ccdsecx1[kk] - 1) / xbin) + 1
                ccdx2 = ccdsecx2[kk] / xbin
                ccdy1 = ((ccdsecy1[kk] - 1) / ybin) + 1
                ccdy2 = ccdsecy2[kk] / ybin
                # Assign the current extensions binned CCDSEC
#                currentccdsec = "["//ccdx1//":"//ccdx2//","//ccdy1//\
#                    ":"//ccdy2//"]"
                currentccdsec = "["//ccdx1//":"//ccdx2//",*]"

                if (debug) {
                    printlog ("GQECORR - binned ccdsec for current "//\
                        "extension is: "//currentccdsec, \
                        l_logfile, verbose+)
                }

                # Loop over lastcorrimage SCI extensions
                for (ll = 1; ll <= num_ccds; ll +=1 ) {

                    # Assign extension to look at
                    thisextn = "["//l_sciext//ll//"]"

                    # Read CCDNAME
                    keypar (lasttmpcorrim//thisextn, l_key_ccdname, silent+)
                    corrccdname = keypar.value

                    # This CCDNAME shouldmatch the stored one
                    if (corrccdname != ccdname[ll]) {
                        printlog ("GQECORR WARNING - corrccdname and "//\
                            "ccdname["//ll//"] don't match. Exiting.", \
                            l_logfile, verbose+)
                        goto clean
                    }

                    if (inccdname == corrccdname) {
                        if (debug) {
                            printlog ("GQECORR - current extensions "//\
                                "CCDNAME is: "//inccdname, \
                                l_logfile, verbose+)
                        }
                        break
                    }
                } # End of loop over correction image CCDs

                # Imcopy out the input CCDSEC from the lasttmpcorr image to
                # l_corrimg
                imcopy (lasttmpcorrim//thisextn//currentccdsec, \
                    l_corrimg//outextnappend, verbose=debug)

                # Update header keywords

                # DETSEC
                gemhedit (l_corrimg//outextn, l_key_detsec, detsecin[kk], \
                    "Detector section(s)", delete-)
                # CCDSEC
                gemhedit (l_corrimg//outextn, l_key_ccdsec, ccdsecin[kk], \
                    "CCD section(s)", delete-)
                # DATASEC
                gemhedit (l_corrimg//outextn, l_key_datasec, \
                    "[1:"//str(ccdx2 - (ccdx1 - 1))//",1:"//\
                    str(ccdy2 - (ccdy1 - 1))//"]", \
                    "Data section(s)", delete-)
                # CCDSUM
                gemhedit (l_corrimg//outextn, l_key_ccdsum, xbin//" "//ybin, \
                    "CCD sum", delete-)

            } # End of loop over input SCI extensions

            # Clean up a bit
            imdelete (lasttmpcorrim, verify-, >& "dev$null")

UPDATE_CORRIM_PHU:

            # Update some PHU header keywords

            # REFIMG - Reference image used
            gemhedit (l_corrimg//"[0]", l_key_refimage, l_refimg,
                "Reference image used to create QE correction", \
                 delete-)
            # GQENOD - Whether the image has been corrected for N&S
            if (imgshift > 0) {
                nodanswer = "YES"
                gemhedit (l_corrimg//"[0]", l_key_gqenod, nodanswer,
                    "Correction image been corrected for N&S offset", \
                    delete-)
                # Likelihood is that reference image will not have NODPIX
                gemhedit (l_corrimg//"[0]", l_key_noddistance, imgshift*ybin,
                    "Number of rows shuffled", delete-)

            }

            # Obtain current time in UT
            gemdate(zone=gemzone)
            # Update GQECORR timestamp
            gemhedit (l_corrimg//"[0]", l_key_gqecorr, gemdate.outdate, \
                "UT Time stamp for GQECORR", delete-)
            # Updaet GEM-TLM
            gemhedit (l_corrimg//"[0]", l_key_gemtlm, gemdate.outdate, \
                "Last modification with GEMINI", delete-)

NEXT_REFIMAGE:

        } # End of for ii loop over refimages creating correction images

        # Reset num_corrimgs to the number of ref images - may need to be
        # made more robust
        num_corrimgs_to_use = num_refimages

        # End of loop creating correction images for spectroscopic data

    } else { # All imaging data

        # Reset the corrimage files to use counter to zero
        num_corrimgs_to_use = 0

        for (ii = 1; ii<= num_inimages; ii += 1) {

            # Assign correction data file output name
            l_corrimg = corrimage_files[ii]

            # Append .dat to output name
            l_corrimg = l_corrimg//".dat"

            # Store l_corrimg again
            corrimage_files[ii] = l_corrimg

            # Inititate corrim_skip_flag
            corrim_skip_flag = no

            if (debug) {
                printlog ("In imaging: "//ii//", "//l_corrimg//", "//\
                    corrim_exists[ii], \
                    l_logfile, verbose+)
            }

            # Check if the correction image exists for this file
            # If not correcting skip to the next inimage
            # If it doesn't exist reset l_corrimg
            if (corrim_exists[ii] && !l_correct_inimages) {
                num_corrimgs_to_use += 1
                goto NEXT_INIMAGE

            } else if (corrim_exists[ii] && l_correct_inimages) {
                # Correction image exists and correcting the inputs
                # Need to read the QE factors in, reset corrim_skip_flag
                corrim_skip_flag = yes

                # Set tmp name for file
                tmpqe2data = l_corrimg

                # Set it up to keep origiinal version of userdefined flag
                userdefined_orig = userdefined
                userdefined = yes
            }

            # Set input file
            l_input = inimage_files[ii]

            # Obtain the number of science extensions - checked already
            keypar (l_input//"[0]", l_key_nsciext, silent+)
            num_sci = int(keypar.value)

            # Acess DETTYPE and any other keywords needed to ditinguish
            # between detectors with different pixel dimensions etc.
            # DETTYPE already checked earlier on
            keypar (l_input//"[0]", l_key_dettype, silent+)
            thisdettype = keypar.value

            # DETECTOR - already checked
            keypar (l_input//"[0]", l_key_detector, silent+)
            thisdetector = keypar.value

            # INSTRUMENT - already checked
            keypar (l_input//"[0]", l_key_instrument, silent+)
            thisinstrument = keypar.value

            # Initiate CCDNAME counter
            num_ccds = 1

            # Loop over the science frames of the input image and read CCDNAME
            for (kk = 1; kk <= num_sci; kk += 1) {

                # Set up the extension to read
                thisextn = "["//l_sciext//kk//"]"

                # Read the CCDNAME
                keypar (l_input//thisextn, l_key_ccdname, silent+)
                thisccd = keypar.value

                if (kk == 1) {
                    # Assign the first CCNAME array value
                    ccdname[num_ccds] = thisccd
                } else {
                    # Compare the current CCDNAME to the stored values
                    # Then store if required
                    for (ll = 1; ll <= num_ccds; ll += 1) {
                        if ((ccdname[num_ccds] != thisccd) && \
                            (ll == num_ccds)) {
                                # Store the CCDNAME
                                num_ccds += 1
                                # There should only be 3 CCDs
                                if (num_ccds > 3) {
                                    # Exit
                                    printlog ("GQECORR ERROR - More than 3 "//\
                                        l_key_ccdname//"s found. Exiting.", \
                                        l_logfile, verbose+)
                                    goto clean
                                } else {
                                    ccdname[num_ccds] = thisccd
                                }
                         } # End of if comparison of CCDNAME and stored CCDNAME
                    } # End of ll loop through stored CCDNAMES
                } # End of CCDNAME checks / storage
            } # End of kk loop over l_input SCI extensions

            # Break if only 1 CCDs worth (assumed to be CCD2) and only
            # doing a relative correction
            if (num_ccds == 1 && l_fl_relative_corr) {
                printlog ("GQECORR ERROR - Only one CCDs worth of "//\
                    "supplied (assuming CCD2) and a relative correction "//\
                    " was requested. This is not possible. Exiting", \
                    l_logfile, verbose+)
            } # End of CCD2 and relative corr request check

            # Read the filter values if first image or not all the images have
            # the same filters
            if ((ii == 1) || (!fl_filter_match)) {

                # FILTER 1 - already checked
                keypar (l_input//"[0]", l_key_filter1, silent+)
                thisfilter1 = keypar.value

                # Make the open filter1 string equal to open
                if (substr(thisfilter1,1,4) == "open") {
                    thisfilter1 = "open"
                }

                # FILTER 2 - already checked
                keypar (l_input//"[0]", l_key_filter2, silent+)
                thisfilter2 = keypar.value

                # Make the open filter2 string equal to open
                if (substr(thisfilter2,1,4) == "open") {
                    thisfilter2 = "open"
                }

                # Loop over CCDNAMEs
                for (ll = 1; ll <= num_ccds; ll += 1) {

                    # Set CCDNAME
                    thisccd = ccdname[ll]

                    # Due to the extreme number of white space in the CCDNAMEs
                    # the CCDNAMEs in the data file are stored with the spaces
                    # replaced with underscores. So, need an _ version of
                    # CCDNAME
                    print (thisccd) | translit ("STDIN", from_string=" ", \
                        to_string="_", delete=no, collapse=no) | \
                            scan(thisccd_data)

                    # Set obsmode
                    thismode = "IMAGE"
                    if (debug) {
                        printlog (thisinstrument//" "//thisdettype//" "//\
                            thisdetector//" "//correction_type//" "//thisccd//\
                            " "//thismode//", "//thisccd_data,\
                            l_logfile, verbose+)
                    }

                    # Parse the QE correction factors
                    if (!userdefined) {

                        # Parse the QE correction to tmptext file
                        # The fields are correct as off 2011-09-05
                        match (pattern=stopdatekey, files=l_qecorr_data, \
                            stop=yes, \
                            print_file_names=no) | \
                            match (pattern=thisinstrument, stop=no, \
                                print_file_names=no) | \
                            match (pattern=thisdettype, stop=no, \
                                print_file_names=no) | \
                            match (pattern=thisdetector, stop=no,\
                                print_file_names=no) | \
                            match (pattern=correction_type, stop=no, \
                                print_file_names=no) | \
                            match (pattern=thisccd_data, stop=no, \
                                print_file_names=no) | \
                            match (pattern=thismode, stop=no, \
                                print_file_names=no) | \
                            match (pattern=thisfilter1, stop=no, \
                                print_file_names=no) | \
                            match (pattern=thisfilter2, stop=no, \
                                print_file_names=no) | \
                            fields (fields="1,5", \
                                lines="1-", quit_if_missing=no) |\
                            scan (dataccdname, thiscoeff)

                            # Remove the " "s and add __s for comparison
                            print (ccdname[ll]) | translit ("STDIN", \
                                 from_string=" ", to_string="_", delete=no, \
                                 collapse=no) | scan(tmpccdname)

                            # Double check the CCDNAMES
                            if (tmpccdname != dataccdname) {
                                printlog ("GQECORR ERROR - "//\
                                    l_key_ccdname//\
                                    " mismatch between "//l_input//" and "//\
                                    "information in "//l_qecorr_data//\
                                    ". Exiting", \
                                    l_logfile, verbose+)
                                goto clean
                            }

                    } else {
                        # user supplied numbers on command line
                        # Assumes the user ordered the numbers according to
                        # CCD, i.e., 1-3, then by coefficient order, 0-n

                        # Set up file to read appropriately
                        if (corrim_skip_flag) {
                            # Correction data file already exists
                            qefile_to_read = tmpqe2data
                            userdefined = userdefined_orig
                        } else {
                            # User input
                            qefile_to_read = tmpqedata
                        }

                        # One line per CCD - assumed to be in order 1-3
                        # and each line will be one string
                        fields (qefile_to_read, fields="1", lines=ll, \
                            quit_if_missing=no, print_file_names=no) | \
                            translit (infile="STDIN", from_string=",", \
                            to_string="\n", delete=no)| \
                                scan (thiscoeff)
                    }

                    # Store or write information to disk
                    if (l_correct_inimages) {
                        # Store for later as correcting the inout images
                        img_corrfact[ii,ll,1] = thisccd
                        img_corrfact[ii,ll,2] = thiscoeff

                    }

                    if (l_keep_corrimages) {
                        # Write to disk
                        # Print to file to be parsed by another task or user
                        print (thiscoeff//" "//thisccd_data, >> l_corrimg)

                    } # End of storing or writing to disk

                } # End of ll loop over CCDNAMEs to parse QE

                # Increment the number of corr images to use by one
                num_corrimgs_to_use += 1

                # End of if loop reading filters, parsing QE info and
                # storage of info
            } else {
                # ii > 1 and fl_filter match is true - no need to continue
                break
            } # End of if i==1 or !fl_filter_match

NEXT_INIMAGE:
        } # End of ii loop over input images
    } # End of imaging data

    # End of QE correction image creation / gathering of QE imaging factors

    ################

    ################
    # Input image correction

    # Correct input images if desired
    if (l_correct_inimages) {

        # Initiate extensions to work on
        extn_list[1] = l_sciext
        extn_list[2] = l_varext
        extn_list[3] = l_dqext

        # Loop over input images
        for (ii = 1; ii<= num_inimages; ii += 1) {

            # Inititate extn to loop over read
            extns_to_correct = 1

            # Set input file
            l_input = inimage_files[ii]

            printlog ("GQECORR - Correcting input image "//l_input, \
                l_logfile, l_verbose)

            # Set output file
            l_output = outimage_files[ii]

            # Read the number of science extensions from l_input
            keypar (l_input//"[0]", l_key_nsciext, silent+)
            num_sci = int(keypar.value)

            # Set correction image
            if (num_corrimgs_to_use == 1) {
                hh = 1
            } else {
                hh = ii
            }

            if (l_fl_vardq) {
                # Check for VARDQ planes exist
                gemextn (l_input,process="expand", \
                    check="exists", extver="1-", \
                    extname=l_var_ext//","//l_dq_ext, \
                    outfile="dev$null")

                if (gemextn.status == 0) { # GEMEXTN was successful
                    if ((gemextn.count/2) == num_sci) { # l_var_ext exist
                        printlog ("          "//l_var_ext//" planes found"//\
                            " in "//l_input, l_logfile, l_verbose)

                        # Increase the number of extensions to correct by 2
                        extns_to_correct += 2
                    } else {
                        printlog ("GQECORR WARNING - No "//l_var_ext//" or "//\
                            l_dq_ext//" planes "//\
                            "found or their number do not match the",\
                            l_logfile, l_verbose)
                        printlog ("                  number of "//l_sci_ext//\
                            " extensions in "//l_input//".", \
                            l_logfile, l_verbose)
                        printlog ("                  Cannot propagate the "//\
                            "variance for this image.", l_logfile, l_verbose)
                    }
                }  else { # GEMEXTN returned an error
                    printlog ("GQCORR ERROR - Call to GEMEXTN returned "//\
                        "a non-zero status for image "//l_input//\
                        ". Exiting (Error 013)", logfile, verbose+)
                        goto clean
                } # End of GEMEXTN checks
            } # End of VARDQ checks

            # Copy the phu of the inut file to the output file
            imcopy (l_input//"[0]", l_output, verbose=debug)

            # Loop over the types of extensions to correct
            for (jj = 1; jj <= extns_to_correct; jj += 1) {

                # Loop over the number of SCI extensions
                for (kk = 1; kk <=num_sci; kk += 1) {

                    # Set the extn to read
                    thisextn = "["//extn_list[jj]//kk//"]"

                    # Set extension to write to
                    outextn = "["//extn_list[jj]//kk//",append]"

                    # imcopy out the DQ planes
                    if (extn_list[jj] == l_dqext) {
                        imcopy (l_input//thisextn, \
                            l_output//outextn, verbose=debug)
                    # Use imarith to perform calculations
                    } else {

                        # Set the DATASEC of the input image
                        keypar (l_input//thisextn, l_key_datasec, silent+)
                        indatasec = str(keypar.value)

                        if (isimaging) {

                            # Read the CCDNAME from l_input
                            keypar (l_input//thisextn, l_key_ccdname, silent+)
                            thisccd = str(keypar.value)

                            # Set the default thiscoeff to ""
                            secondop = ""
                            if (debug) {
                                printlog (secondop, l_logfile, verbose+)
                                printlog ("----------*******-----------", \
                                    l_logfile, verbose+)
                            }

                            # Loop over stored coeeficient information
                            # Making the assumption that there is a max of
                            # 3 CCDs
                            for (ll = 1; ll <= 3; ll += 1) {
                                if (debug) {
                                    printlog (img_corrfact[hh,ll,1]//", "//\
                                        img_corrfact[hh,ll,2], \
                                        l_logfile, verbose+)
                                }

                                if (thisccd == img_corrfact[hh,ll,1]) {
                                    secondop = img_corrfact[hh,ll,2]
                                    break
                                }
                            } # End of ll loop over stored datat for l_input

                            if (secondop == "") {
                                printlog ("GQECORR ERROR - cannot "//\
                                    "match CCD names with stored CCD "//\
                                    "names for "//l_input//". Exiting.", \
                                    l_logfile, verbose+)
                                goto clean
                            }

                            # Set up printlog statement for imaging QE factor
                            l_struct_description = ""
                            print ("          QE factor for extension "//\
                                thisextn//" is "//factor_description) | \
                                    scan (l_struct_description)
                            if (extn_list[jj] != l_varext) {
                                printlog (l_struct_description//secondop, \
                                    l_logfile, l_verbose)
                            } else {
                                printlog (l_struct_description//"("//\
                                    secondop//")**2", l_logfile, l_verbose)
                            }
                        } else {
                            # Spectroscopic data

                            # Assign correction image
                            l_corrimg = corrimage_files[hh]
                            l_refimg = refimage_files[hh]

                            if (kk == 1 && extn_list[jj] == l_sciext) {
                                printlog ("          Using correction "//\
                                    "image "//l_corrimg, \
                                    l_logfile, l_verbose)
                            }

                            # Set the secondop to be the correct SCI extn in
                            # the correction image
                            # Corrorction image will only have SCI extns
                            secondop = l_corrimg//"["//l_sciext//kk//"]"
                        } # End of testing whether imaging or spectroscopy

                        # Set up imexpression depending on extension
                        calexpr = ""
                        if (extn_list[jj] == l_sciext) {
                            calexpr = "(a)"//corrapp//"(b)"
                        } else {
                            calexpr = "(a)"//corrapp//"(b**2)"
                        }

                        # Correct the data
                        imexpr (calexpr, l_output//outextn,
                            l_input//thisextn//indatasec, \
                            secondop, dims="auto", \
                            intype="auto", \
                            outtype=imgpixtype, bwidth=0, \
                            btype="nearest", bpixval=0, verbose=no)

                        # Update header. Imaging only.
                        if (isimaging) {
                            # Write QE factor to header
                            gemhedit (l_output//thisextn, l_key_imgfactor, \
                                secondop, "QE correction factor", delete-)
                        }

                    } # End of test of type of extn
                } # End of kk loop over the num_sci (extver)
            } # End of jj loop over extns_to_correct

            # If spectroscopic data
            if (!isimaging) {
                # Copy out MDF

                # Initiate tinfo.tabletype
                tinfo.tbltype = ""

                # read the MDF table imformation
                tinfo (l_input//"["//l_mdf_ext//"]", ttout=no, >& "dev$null")

                # Test the tbltype
                if (tinfo.tbltype == "fits") {
                    # Copy out MDF
                    ttools.tcopy (intable=l_input//"["//l_mdf_ext//"]", \
                        outtable=l_output//"["//l_mdf_ext//"]", \
                        verbose-)
                }

                # Write either the correction image or the reference frame
                # name to the header
                gemhedit (l_output//"[0]", l_key_corrim, l_corrimg, \
                    "QE correction image", delete-)

                # Read the reference image from the correction image
                keypar (l_corrimg//"[0]", l_key_refimage, silent+)
                # Write it to the header
                gemhedit (l_output//"[0]", l_key_refimage, \
                    keypar.value,
                    "Reference image used to create QE correction", \
                    delete-)
            }

            # Update the number of extensions properly
            gemextn (l_output, check="exists", process="expand", index="1-", \
                extname="", extversion="", ikparams="", omit="", replace="", \
                outfile="dev$null", logfile=l_logfile, glogpars="", \
                verbose=no)

            if (gemextn.status != 0) {# GEMEXTN returned an error
                printlog ("GQCORR ERROR - Call to GEMEXTN returned "//\
                    "a non-zero status for image "//l_input//\
                    ". Exiting (Error 033)", logfile, verbose+)
                goto clean
            } else {
                n_extentions = gemextn.count
            } # End of GEMEXTN checks

            # Update the number of extensions present
            gemhedit (l_output//"[0]", "NEXTEND", n_extentions, \
                "", delete-, upfile="")

            # Obtain current time in UT
            gemdate(zone=gemzone)
            # Update GQECORR timestamp
            gemhedit (l_output//"[0]", l_key_gqecorr, gemdate.outdate, \
                "UT Time stamp for GQECORR", delete-)
            # Updaet GEM-TLM
            gemhedit (l_output//"[0]", l_key_gemtlm, gemdate.outdate, \
                "Last modification with GEMINI", delete-)

        } # End of for ii loop over num_inimages
    } # End of correcting inimages

    # Clean up the correction images if needed.

    # If not keeping the correction image - spectroscopic data only
    # (Imaging no corr files are created)
    if (!l_keep_corrimages) {

        # Loop over correction files and delete only if they
        # didn't already exist
        for (ii = 1; ii <= num_corrimages; ii += 1) {
            if (!corrim_exists[ii]) {
                delete (corrimage_files[ii], verify-, >& "dev$null")
            }
        }

    } # End of cleaning up correction images

    ########################
    # End of the hard work #
    ########################

FILE_FOR_GIREDUCE:

    if (debug) {
        printlog ("Num of CORRIMGS is: "//num_corrimages, \
            l_logfile, l_verbose)
        printlog ("Num of CORRIMGS to use is: "//num_corrimgs_to_use, \
            l_logfile, l_verbose)
    }

    # To allow easy use in gireduce
    # Write the correction images to a specified file
    if (scanfile2 != "") {
        for (ii = 1; ii <= num_corrimgs_to_use; ii += 1) {
            print (corrimage_files[ii]//" "//corrim_exists[ii], >> \
                scanfile2)
            if (debug) {
                printlog (corrimage_files[ii]//" "//corrim_exists[ii], \
                    l_logfile, verbose+)
            }
        }
        scanfile2 = ""
    }

    # If we got to here everything went OK
    status = 0

clean:

    # Clean up a bit - delete inlist, outlist, reflist, corimlist etc.
    scanfile = ""
    delete (inlist//", "//outlist, verify-, >& "dev$null")

    # The following may not exist
    delete (reflist//", "//corrimlist, verify-, >& "dev$null")

    if (userdefined) {
        delete (tmpqedata, verify-, >& "dev$null")
    }

    # Print finish time
    gemdate (zone="local")
    printlog ("\nGQECORR - Finished: "//gemdate.outdate, \
        l_logfile, verbose+)

    if (status == 0) {
        printlog ("GQECORR - Exit status: GOOD\n", l_logfile, verbose+)
    } else {
        printlog ("GQECORR - Exit status: ERROR\n", l_logfile, verbose+)
    }

end
