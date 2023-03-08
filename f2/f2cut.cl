# Copyright(c) 2012-2017 Association of Universities for Research in Astronomy, Inc.

# F2CUT - Cut spectroscopic FLAMINGOS-2 data
# Full details of this task are given in the help file

procedure f2cut(inimages)

char    inimages    {prompt = "Input FLAMINGOS-2 image(s)"}
char    outimages   {"", prompt = "Output image(s)"}
char    outprefix   {"c", prompt = "Prefix for output image(s)"}
bool    fl_vardq    {yes, prompt = "Propagate variance and data quality planes"}
char    gradimage   {"", prompt = "Image to use for finding slit edges using the gradient method (MOS)"}
char    outslitim   {"slits.fits", prompt = "Output image to be used to determine s-distortion correction (MOS)"}
char    minsep      {"default", prompt = "Minimum separation between spectra in pixels (MOS, apfind)"}
real    edgediff    {2.0, prompt = "Allowed difference in pixels when associating slit edges to a slit (MOS)"}
char    refimage    {"", prompt = "Reference image for slit positions (MOS)"}
char    database    {"", prompt = "Directory for files containing slit trace (MOS)"}
char    logfile     {"", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose"}
int     status      {0, prompt = "Exit status (0=good)"}
struct  *scanfile1  {"", prompt = "Internal use only"}
struct  *scanfile2  {"", prompt = "Internal use only"}
struct  *scanfile3  {"", prompt = "Internal use only"}

begin

    char    l_inimages = ""
    char    l_outimages = ""
    char    l_outprefix = ""
    bool    l_fl_vardq
    char    l_gradimage = ""
    char    l_outslitim = ""
    char    l_minsep = ""
    real    l_edgediff
    char    l_refimage = ""
    char    l_database = ""
    char    l_logfile = ""
    bool    l_verbose

    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_pixscale = ""
    char    l_key_grating = ""
    char    l_key_filter = ""
    char    l_key_slit = ""
    char    l_key_cut_section = ""
    char    l_key_delta = ""
    char    l_key_order = ""

    int     junk, input_count, output_count, detid, slit_num, last_row
    int     nx, ny, yoffset, slit_width, half_slit_width, width_even
    int     slit_height, half_slit_height, height_even, ix1, ix2, iy1, iy2, iy
    int     atlocation, slits, l_dispaxissave, priority, extver, count
    int     appos, apneg, obs_date_seconds
    real    slitpos_mx, slitpos_my, slitsize_mx, slitsize_my, x_ccd, y_ccd
    real    slitsize_x, slitsize_y, x1, x2, y1, y2, dx, dy, y_cen, slitpos_y
    real    pixscale, filter_lower, filter_upper, filter_width
    real    orig_dispersion, dispersion, asecmm, grad_left_edge, mdf_left_edge 
    real    final_left_edge, grad_right_edge, mdf_right_edge, final_right_edge
    real    slit_offset_left, slit_offset_right, tolerance, xmin, xmax, ymin
    real    ymax, coeff1, coeff2, coeff3,x,n, y, tmpvalue, aplow, aphigh
    real    lsltedge, rsltedge, bpmx1, bpmx2, bpmy1, bpmy2, min_width
    real    apfind_minsep, edge[2048]
    char    tmpinlist, tmpinput, tmpoutput, tmpgradlist, tmpreflist, tmptab
    char    tmpcenter, tmpapnum, tmpdb, tmpgradpos, tmpgradneg, tmpgradleft
    char    tmpmdfleft, tmpgradright, tmpmdfright, tmptmatch, tmpbpmfile
    char    tmpbpm, tmpdq, tmptodel, l_logfilesave
    char    l_databasesave, badhdr, input_file, grad_file, ref_file, inimg
    char    line, input, output, phu, mdf, sci, var, dq, outsci, filterdb
    char    offsetdb, nsappwavedb, grism, filter, slit, section, mdftest
    char    filter1, filter2, lyot
    char    tmpgradim, tmprefim, word, curve, func, ord, params[7], maskexpr
    char    obs_date
    bool    debug, seeplot, longslit, mos, useref, usegrad
    bool    usetable, l_verbosesave, stored, getvalues, dolefttrace
    bool    dorighttrace, maskoverlap, maketmpbpm, updatedq, first, noleftedge
    struct  l_struct

    status = 1
    debug = no
    seeplot = no
    useref = no
    usegrad = no
    usetable = no
    maskoverlap = no
    maketmpbpm = no
    updatedq = no
    noleftedge = no

    # asecmm is the plate scale (the scale of images in the focal plane)
    # measured in arcsec/mm (see section 3.9 of the f/16 optical design summary
    # http://www.gemini.edu/documentation/webdocs/rpt/rpt-o-g0047.pdf). For an
    # 8m f/16 telescope, where the focal distance = 8m / (1/16) = 128000mm:
    #     206264.8 [arcsec/radian] * (1 / 128000) [radian/mm] 
    #         = 1.611444 [arcsec/mm]
    asecmm = 1.611444

    # Cache tasks used throughout the script
    cache ("keypar", "gemdate")

    # Set the local variables
    l_inimages = inimages
    l_outimages = outimages
    l_outprefix = outprefix
    l_fl_vardq = fl_vardq
    l_gradimage = gradimage
    l_outslitim = outslitim
    l_minsep = minsep
    l_edgediff = edgediff
    l_refimage = refimage
    l_database = database
    l_logfile = logfile
    l_verbose = verbose

    # Shared definitions, define parameters from nsheaders
    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_slit, l_key_slit)
    if ("" == l_key_slit) badhdr = badhdr + " key_slit"
    junk = fscan (nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"
    junk = fscan (nsheaders.key_grating, l_key_grating)
    if ("" == l_key_grating) badhdr = badhdr + " key_grating"
    junk = fscan (nsheaders.key_filter, l_key_filter)
    if ("" == l_key_filter) badhdr = badhdr + " key_filter"
    junk = fscan (nsheaders.key_cut_section, l_key_cut_section)
    if ("" == l_key_cut_section) badhdr = badhdr + " key_cut_section"
    junk = fscan (nsheaders.key_delta, l_key_delta)
    if ("" == l_key_delta) badhdr = badhdr + " key_delta"
    junk = fscan (nsheaders.key_order, l_key_order)
    if ("" == l_key_order) badhdr = badhdr + " key_order"

    # Make temp files:
    tmpinlist = mktemp ("tmpinlist") 
    tmpgradlist = mktemp ("tmpgradlist") 
    tmpreflist = mktemp ("tmpreflist") 
    tmpinput = mktemp ("tmpinput")
    tmpoutput = mktemp ("tmpoutput")
    tmpgradpos = mktemp ("tmpgradpos")
    tmpgradneg = mktemp ("tmpgradneg")
    tmptab = mktemp ("tmptab")//".fits"
    tmpcenter = mktemp("tmpcenter")
    tmpapnum = mktemp("tmpapnum")
    tmpdb = mktemp("tmpdb")
    tmpgradleft = mktemp("tmpgradleft")//".fits"
    tmpmdfleft = mktemp("tmpmdfleft")
    tmpgradright = mktemp("tmpgradright")//".fits"
    tmpmdfright = mktemp("tmpmdfright")
    tmptmatch = mktemp("tmptmatch")//".fits"
    tmpbpmfile = mktemp ("tmpbpmfile")
    tmpbpm = mktemp ("tmpbpm")
    tmpdq = mktemp ("tmpdq")
    tmptodel = mktemp("tmptodel")

    phu = "[0]"
    mdf = ".fits[MDF]"
    sci = "[" // l_sci_ext // ",1]"
    var = "[" // l_var_ext // ",1]"
    dq = "[" // l_dq_ext // ",1]"

    # Define the look-up tables used to calculate the section to cut for 
    # FLAMINGOS-2 data
    filterdb = "f2$data/f2filters.fits"
    nsappwavedb = "gnirs$data/nsappwave.fits"

    # Test that name of logfile makes sense
    if (l_logfile == "" || l_logfile == " ") {
        junk = fscan (f2.logfile, l_logfile) 
        if (l_logfile == "" || l_logfile == " ") {
            l_logfile = "f2.log"
            printlog ("WARNING - F2CUT: Both f2cut.logfile and \
                f2.logfile are undefined.", l_logfile, verbose+)
            printlog ("                  Using " // l_logfile, l_logfile, 
                verbose+)
        }
    }

    # Open logfile
    date | scan (l_struct)
    printlog ("--------------------------------------------------------------\
        --------------", l_logfile, l_verbose)
    printlog ("F2CUT -- " // l_struct, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    if ("" != badhdr) {
        printlog ("ERROR - F2CUT: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+)
        goto clean
    }

    if (l_database == "") {
        junk = fscan (f2.database, l_database)
        if (l_database == "") {
            l_database = "database"
            printlog ("WARNING - F2CUT: Both f2cut.database and \
                f2.database are", l_logfile, verbose+)
            printlog ("                   undefined. Using " // l_database, \
                l_logfile, verbose+)
        }
    }

    # Remove any trailing slashes from the database name
    if (strlen (l_database) > 1 && \
        substr (l_database, strlen (l_database), strlen (l_database)) \
        == "/") {
        l_database = substr (l_database, 1, strlen (l_database) - 1)
    }

    # Store default specred values
    l_logfilesave = specred.logfile
    l_databasesave = specred.database
    l_dispaxissave = specred.dispaxis
    l_verbosesave = specred.verbose
    stored = yes

    specred.logfile = l_logfile
    specred.database = l_database
    specred.dispaxis = 2
    specred.verbose = l_verbose

    if (debug) {
        print ("TIMESTAMP: starting checks")
        date "+%H:%M:%S.%N"
    }

    # Check that list file exists
    atlocation = stridx ("@", l_inimages)
    if (atlocation > 0) {
        input_file = substr (l_inimages, atlocation + 1, strlen (l_inimages))
        if (!access (input_file)) {
            printlog ("ERROR - F2CUT: Input file " // input_file // " not \
                found.", l_logfile, verbose+)
            goto clean
        }
    }

    # Check that all input images exist and have a [SCI,1] extension
    if (debug) print ("exist check")
    gemextn (l_inimages, check="exists,mef,image", process="append", \
        index="", extname=l_sci_ext, extversion="1", ikparams="", \
        omit="", replace="", outfile="dev$null", logfile="", glogpars="",
        verbose=l_verbose)
    if (gemextn.fail_count != 0 || gemextn.count == 0) {
        printlog ("ERROR - F2CUT: One or more of the input images are missing \
            science data.", l_logfile, verbose+) 
        goto clean
    }
    if (l_fl_vardq) {

        # Check that all input images have a [VAR,1] and a [DQ,1] extensions
        gemextn (l_inimages, check="exists,mef,image", process="append", \
            index="", extname=l_var_ext, extversion="1", ikparams="", \
            omit="", replace="", outfile="dev$null", logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - F2CUT: One or more of the input images are \
                missing variance data.", l_logfile, verbose+) 
            goto clean
        }
        gemextn (l_inimages, check="exists,mef,image", process="append", \
            index="", extname=l_dq_ext, extversion="1", ikparams="", \
            omit="", replace="", outfile="dev$null", logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - F2CUT: One or more of the input images are \
                missing data quality information.", l_logfile, verbose+) 
            goto clean
        }
    }

    # Create a list of input file names
    gemextn (l_inimages, check="", process="append", index="", 
        extname=l_sci_ext, extversion="1", ikparams="", omit="kernel,exten",
        replace="", outfile=tmpinlist, logfile="", glogpars="",
        verbose=l_verbose)
    input_count = gemextn.count

    scanfile1 = tmpinlist
    longslit = no
    mos = no
    while (fscan (scanfile1, inimg) != EOF) {

        if (debug) print ("inimg: " // inimg)

        # Determine whether the input images are longslit or MOS
        keypar (inimg // phu, l_key_slit, silent+)
        if (keypar.found) {
            slit = str(keypar.value)
        } else {
            printlog ("ERROR - F2CUT: No header value for " \
                // l_key_slit // " in " //inimg // ".", l_logfile, verbose+)
            goto clean
        }
        if (debug) print ("slit: " // slit)

        # Determine whether the input images have been cut previously
        if (debug) print ("f2cut check")
        keypar (inimg // phu, "F2CUT", silent+) 
        if (keypar.found) {
            printlog ("WARNING - F2CUT: File " // inimg // " has already been \
                run through F2CUT. No further processing will be performed.", \
                l_logfile, verbose+) 
        } else {
            if (substr (slit, 1, (strlen (slit) - 1)) == "mos") {
  	        mos = yes
		slit = "mos"
	    }

            if (substr (slit, strlen (slit) - 3, strlen (slit)) == "slit")
	        longslit = yes

	    if (mos && longslit) {
                printlog ("ERROR - F2CUT: Mixture of longslit and MOS data \
                    in input images.", l_logfile, verbose+)
                goto clean
            }

            print (inimg, >> tmpinput) 
        }

        # Determine whether the input images have been prepared
        if (debug) print ("prepare check")
        keypar (inimg // phu, "PREPARE", silent+)
        if (keypar.found == no) {
            printlog ("ERROR - F2CUT: Image " // inimg // " has not been \
                processed by F2PREPARE.", l_logfile, verbose+)
            goto clean
        }
    }
    scanfile1 = ""

    if (!access (tmpinput)) {
        printlog ("ERROR - F2CUT: no input files to process.", l_logfile, \
            verbose+)
        goto clean
    }

    # For MOS data, check that all input files have an [MDF] extension
    if (mos) {
        gemextn ("@" // tmpinput, check="exists,mef,table", process="append", \
            index="", extname="MDF", extversion="", ikparams="", omit="", \
            replace="", outfile="dev$null", logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - F2CUT: One or more input images do not have an \
                MDF attached.", l_logfile, verbose+)
            goto clean
        }
    }

    # Check the gradimage and refimage parameters
    if (l_gradimage != "" && l_refimage == "") {

        # Use the gradient image to detect the slit edges
        if (debug) print ("gradient check")

        # Check the list file exists
        atlocation = stridx ("@", l_gradimage)
        if (atlocation > 0) {
            grad_file = substr (l_gradimage, atlocation + 1, \
                strlen (l_gradimage))
            if (!access (grad_file)) {
                printlog ("ERROR - F2CUT: Gradient file " // grad_file // " \
                    not found.", l_logfile, verbose+)
                goto clean
            }
        }

        # Check that the gradient image exists and create a list of gradient
        # file names 
        gemextn (l_gradimage, check="exists", process="append", index="", 
            extname=l_sci_ext, extversion="1", ikparams="", 
            omit="kernel,exten", replace="", outfile=tmpgradlist, logfile="",
            glogpars="", verbose=l_verbose)

        # Only one gradient image should be supplied
        if (gemextn.count == 1) {
            type (tmpgradlist) | scan (tmpgradim)
            printlog ("F2CUT: Using gradient image " // tmpgradim // " to \
                detect the slit edges.", l_logfile, l_verbose)
        } else {
            printlog ("ERROR - F2CUT: A single image must be supplied as the \
                gradient image.", l_logfile, verbose+)
            goto clean
        }

        # Check that the gradient image has a [SCI,1] extension
        gemextn (tmpgradim, check="exists,mef,image", process="append", \
            index="", extname=l_sci_ext, extversion="1", ikparams="", \
            omit="", replace="", outfile="dev$null", logfile="", glogpars="", \
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - F2CUT: Gradient image " // tmpgradim // " is \
                missing science data.", l_logfile, verbose+)
            goto clean
        }

        # Check that the gradient image has an [MDF] extension
        gemextn (tmpgradim, check="exists,mef,table", process="append", \
            index="", extname="MDF", extversion="", ikparams="", omit="",
            replace="", outfile="dev$null", logfile="", glogpars="",
            verbose=l_verbose) 
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - F2CUT: Gradient image " // tmpgradim // " does \
                not have an MDF attached.", l_logfile, verbose+)
            goto clean
        }

        # If the areas contaminated by overlap are to be masked using the
        # information from the trace of the gradient image, fl_vardq must be
        # set to yes
        if (l_fl_vardq) {
            printlog ("F2CUT: Areas contaminated by overlap from other slits \
                will be masked and", l_logfile, l_verbose) 
            printlog ("       recorded in the output DQ extension.", \
                l_logfile, l_verbose) 
            maskoverlap = yes
            maketmpbpm = yes
            updatedq = yes
        } else {
            printlog ("WARNING - F2CUT: Areas contaminated by overlap from \
                other slits will not be masked,", l_logfile, verbose+)
            printlog ("                 since fl_vardq-.", l_logfile, \
                verbose+)
        }

        # Compute the gradient images
        gradient (tmpgradim // sci, tmpgradpos, "0", boundary="nearest", \
            constant=0.) 
        imarith (tmpgradpos, "*", -1., tmpgradneg)

        # Set pixels with values less than zero to zero to help apfind find the
        # correct slit edges and determine an appropriate aperture (sometimes
        # the aperture is resized to be much larger than it needs to be due to 
        # the negative peak in the image from the adjacent slit)
        imreplace (tmpgradpos, value=0.0, imaginary=0.0, lower=INDEF, \
            upper=0.0, radius=0.0)
        imreplace (tmpgradneg, value=0.0, imaginary=0.0, lower=INDEF, \
            upper=0.0, radius=0.0)

        # If this point is reached, the gradient file has been fully checked
        # and validated
        usegrad = yes

    } else if (l_gradimage == "" && l_refimage != "") {

        # Use the reference image to cut input images
        if (debug) print ("reference check")

        # Check that list file exists
        atlocation = stridx ("@", l_refimage)
        if (atlocation > 0) {
            ref_file = substr (l_refimage, atlocation + 1, \
                strlen (l_refimage))
            if (!access (ref_file)) {
                printlog ("ERROR - F2CUT: Reference file " // ref_file // " \
                    not found.", l_logfile, verbose+)
                goto clean
            }
        }

        # Check that the reference image exists and create a list of reference
        # file names 
        gemextn (l_refimage, check="exists", process="append", index="", 
            extname=l_sci_ext, extversion="1", ikparams="", 
            omit="kernel,exten", replace="", outfile=tmpreflist, logfile="",
            glogpars="", verbose=l_verbose)

        # Only one reference image should be supplied
        if (gemextn.count == 1) {
            type (tmpreflist) | scan (tmprefim)
            printlog ("F2CUT: Using the slit position defined in " // \
                tmprefim // " to cut the", l_logfile, l_verbose)
            printlog ("       input images.", l_logfile, l_verbose)
        } else {
            printlog ("ERROR - F2CUT: A single image must be supplied as the \
                reference image.", l_logfile, verbose+)
            goto clean
        }

        # Check that the reference image has a [SCI,1] extension
        gemextn (tmprefim, check="exists,mef,image", process="append",
            index="", extname=l_sci_ext, extversion="1", ikparams="", omit="",
            replace="", outfile="dev$null", logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - F2CUT: Reference image " // tmprefim // " is \
                missing science data.", l_logfile, verbose+)
            goto clean
        } 

        # Check that the reference image has an [MDF] extension
        gemextn (tmprefim, check="exists,mef,table", process="append",
            index="", extname="MDF", extversion="", ikparams="", omit="",
            replace="", outfile="dev$null", logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - F2CUT: Reference image " // tmprefim // " does \
                not have an MDF attached.", l_logfile, verbose+)
            goto clean
        }
 
        # If the areas contaminated by overlap are to be masked using the
        # information contained in the DQ extension of the reference image,
        # fl_vardq must be set to yes and the reference image must have a DQ
        # extension 
        gemextn (tmprefim, check="exists,mef,image", process="append", \
            index="", extname=l_dq_ext, extversion="1", ikparams="", \
            omit="", replace="", outfile="dev$null", logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("WARNING - F2CUT: Reference image " // tmprefim // " is \
                missing data quality information.", l_logfile, verbose+) 
            printlog ("                 Areas contaminated by overlap from \
                other slits will not be masked.", l_logfile, verbose+)
        } else {
            maskoverlap = yes
            if (l_fl_vardq) {
                printlog ("F2CUT: Areas contaminated by overlap from other \
                    slits will be masked and", l_logfile, l_verbose) 
                printlog ("       recorded in the output DQ extension.", \
                    l_logfile, l_verbose) 
                updatedq = yes
            } else {
                printlog ("F2CUT: Areas contaminated by overlap from other \
                    slits will be masked.", l_logfile, l_verbose) 
            }
        }

        # Check if the slit edge information exists in the MDF 
        mdftest = ""
        tlcol (tmprefim // mdf, nlist=4) | match ("SECX2", stop-, \
            print_file_n+, metacharacte+) | scan (mdftest)
        if (mdftest == "") {
            printlog ("ERROR - F2CUT: The MDF in reference image " // \
                tmprefim // " does not have slit edge info.", l_logfile, \
                verbose+)
            goto clean
        }

        # If this point is reached, the reference file has been fully checked
        # and validated
        useref = yes

    } else if (l_gradimage != "" && l_refimage != "") {
        if (longslit) {
            printlog ("ERROR - F2CUT: Neither gradimage nor refimage should \
                be set for longslit data.", l_logfile, verbose+) 
            goto clean
        }
        if (mos) {
            printlog ("ERROR - F2CUT: Either gradimage or refimage should be \
                set, not both.", l_logfile, verbose+) 
            goto clean
        }

    } else {
        # l_gradimage == "" && l_refimage == ""
        # This is expected for longslit data, but for MOS data, one of either
        # gradimage or refimage should be set
        if (mos) {
            printlog ("ERROR - F2CUT: Either gradimage or refimage should be \
                set.", l_logfile, verbose+) 
            goto clean
        }
    }

    # Generate the output images
    if (debug) print ("output check")
    gemextn (l_outimages, check="absent", process="none", index="", \
        extname="", extversion="", ikparams="", omit="kernel,exten", 
        replace="", outfile=tmpoutput, logfile="", glogpars="",
        verbose=l_verbose)
    if (gemextn.fail_count != 0) {
        printlog ("ERROR - F2CUT: Existing or incorrectly formatted output \
            files.", l_logfile, verbose+) 
        goto clean
    }
    output_count = gemextn.count

    # If tmpoutput is empty, the output files names should be created with the
    # prefix parameter 
    if (output_count == 0) {
        if (debug) print ("using prefix")
        gemextn ("%^%" // l_outprefix // "%" // "@" // tmpinput, \
            check="absent", process="none", index="", extname="", \
            extversion="", ikparams="", omit="kernel,exten", replace="", \
            outfile=tmpoutput, logfile="", glogpars="", verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - F2CUT: Existing or incorrectly formatted \
                output files.", l_logfile, verbose+)
            goto clean
        }
    }
    output_count = gemextn.count

    if (input_count != output_count) {
        printlog ("ERROR - F2CUT: Number of input and output files differ.", \
            l_logfile, verbose+) 
        goto clean
    }

    # The outslitim parameter is only used if a gradient image is used to
    # detect the slit edges 
    if (usegrad) {

        # Check that the output files do not already exist
        if (debug) print ("output slit check")
        gemextn (l_outslitim, check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", omit="kernel,exten", 
            replace="", outfile="dev$null", logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0) {
            printlog ("ERROR - F2CUT: Slit image " // l_outslitim // \
                " exists.", l_logfile, verbose+) 
            goto clean
        }
    }

    if (debug) {
        print ("TIMESTAMP: finished checks")
        date "+%H:%M:%S.%N"
    }

    scanfile1 = tmpinput
    scanfile2 = tmpoutput
    first = yes
    while (fscan (scanfile1, input) != EOF && \
        fscan (scanfile2, output) != EOF) {

        if (debug) {
            print ("input: " // input)
            print ("output: " // output)
            print ("TIMESTAMP: obtain header information")
            date "+%H:%M:%S.%N"
        }

        # Determine which offset lookup table to use (all data taken before
        # 2013 requires the use of an engineering lookup table)
        keypar (input // phu, "DATE-OBS", silent+)
        if (keypar.found) {
            obs_date = str(keypar.value)

            # Convert the date to seconds
            cnvtsec (obs_date, "00:00:00") | scan (obs_date_seconds)

            # cnvtsec ("2013-01-01", "00:00:00") returns 1041465600 seconds
            if (obs_date_seconds > 1041465600) {
                offsetdb = "f2$data/f2offsets.fits"
            } else {
                offsetdb = "f2$data/f2offsets_eng.fits"
            }
        } else {
            offsetdb = "f2$data/f2offsets_eng.fits"
        }

        # Obtain information from the headers
        keypar (input // sci, "i_naxis1", silent+)
        if (keypar.found) {
            nx = int (keypar.value)
        } else {
            printlog ("ERROR - F2CUT: No header value for i_naxis1 in " \
                // input // sci // ".", l_logfile, verbose+)
            goto clean
        }

        keypar (input // sci, "i_naxis2", silent+)
        if (keypar.found) {
            ny = int (keypar.value)
        } else {
            printlog ("ERROR - F2CUT: No header value for i_naxis2 in " \
                // input // sci // ".", l_logfile, verbose+)
            goto clean
        }
        if (debug) print ("image size: " // nx // " x " // ny)

        keypar (input // phu, l_key_pixscale, silent+)
        if (keypar.found) {
            pixscale = real (keypar.value)
        } else {
            printlog ("ERROR - F2CUT: No header value for " \
                // l_key_pixscale // " in " // input // phu // ".", \
                l_logfile, verbose+)
            goto clean
        }
        if (first)
            printlog ("\nPixel scale: " // pixscale, l_logfile, l_verbose)
 
        keypar (input // phu, l_key_grating, silent+)
        if (keypar.found) {
            grism = keypar.value
        } else {
            printlog ("ERROR - F2CUT: No header value for " \
                // l_key_grating // " in " // input // phu // ".", l_logfile, \
                verbose+)
            goto clean
        }
        if (first)
            printlog ("Grism: " // grism, l_logfile, l_verbose)

        keypar (input // phu, l_key_filter, silent+)
        if (keypar.found) {
            filter = keypar.value
        } else {
            printlog ("ERROR - F2CUT: No header value for " \
                // l_key_filter // " in " // input // phu // ".", l_logfile, \
                verbose+)
            goto clean
        }

        # KL 20220216
        # Now get back that component ID that was removed in f2prepare...
        # Also, check the LYOT keyword, they are now putting filters in there.
        
        keypar (input // phu, l_key_filter//"1", silent+)
        if (keypar.found) {
            filter1 = keypar.value
        } else { # early ENG data
            keypar (input // phu, "FILT1POS", silent+)
            if (keypar.found) {
                filter1 = keypar.value
            }
        }
        keypar (input // phu, l_key_filter//"2", silent+)
        if (keypar.found) {
            filter2 = keypar.value
        } else {  #  early ENG data
            keypar (input // phu, "FILT2POS", silent+)
            if (keypar.found) {
                filter2 = keypar.value
            }
        }
        keypar (input // phu, "LYOT", silent+)
        if (keypar.found) {
            lyot = keypar.value
        } else { # early ENG data
            keypar (input // phu, "LYOTPOS", silent+)
            if (keypar.found) {
                lyot = keypar.value
            }
        }
        
        if ( (filter == "Open") && \
             (  (substr(lyot, 1, 1) != "f") && \
                (substr(lyot, 1, 4) != "GEMS") && \
                (substr(lyot, 1, 4) != "Hart") \
              ) ) {
            filter = lyot
        } else if ( (substr(filter1, 1, (strlen (filter1) - 6)) == filter) || \
             (substr(filter1, 1, (strlen (filter1) - 11)) == filter) ) {
            filter = filter1
        } else if ( (substr(filter2, 1, (strlen (filter2) - 6)) == filter) || \
             (substr(filter2, 1, (strlen (filter2) - 11)) == filter) ) {
            filter = filter2
        } 
        if (first)
            printlog ("Filter: " // filter // "\n", l_logfile, l_verbose)

        keypar (input // phu, "DETID", silent+)
        if (keypar.found) {
            detid = int(keypar.value)
        } else {
            printlog ("ERROR - F2CUT: No header value for DETID in " \
                // input // phu // ".", l_logfile, verbose+)
            goto clean
        }
        if (debug) print ("detid: " // str(detid))

        # Get the yoffset from the lookup table, which is used to define which 
        # pixel on the array the central wavelength falls on (w.r.t. the centre
        # of the array, i.e., y = 1024)
        tselect (offsetdb, tmptab, "grism == '" // grism // "' && filter == \
            '" // filter // "' && slit == '" // slit // "'")

        tinfo (tmptab, ttout-)
        if (tinfo.nrows == 1)
            usetable = yes

        if (debug) print (tinfo.nrows // ", " // usetable)

        if (usetable) {
            tabpar (table=tmptab, col="yoffset", row=1, format-)
            if (tabpar.undef) {
                printlog ("WARNING - F2CUT: No offset value found in \
                    " // offsetdb // " for", l_logfile, verbose+)
                printlog ("                 grism == " // grism // ", \
                    filter == " // filter // " and slit == " // slit // ".", \
                    l_logfile, verbose+)
                printlog ("                 Using full y range.", \
                    l_logfile, verbose+)
                yoffset = INDEF
            } else {
                yoffset = int(tabpar.value)
            }

            # Get the x1 and x2 values from the lookup table, which are used
            # for cutting longslit data only
            if (longslit) {
                x1 = 30
                x2 = 1525
                tabpar (table=tmptab, col="x1", row=1, format-)
                if (tabpar.undef) {
                    printlog ("WARNING - F2CUT: No value for x1 found in \
                        " // offsetdb // " for", l_logfile, verbose+)
                    printlog ("                 grism == " // grism // ", \
                        filter == " // filter // " and slit == " // slit // \
                        ".", l_logfile, verbose+)
                    printlog ("                 Using default value \
                        (x1 = " // x1 // ").", l_logfile, verbose+)
                    ix1 = int(x1)
                } else {
                    ix1 = int(tabpar.value)
                }
                tabpar (table=tmptab, col="x2", row=1, format-)
                if (tabpar.undef) {
                    printlog ("WARNING - F2CUT: No value for x2 found in \
                        " // offsetdb // " for", l_logfile, verbose+)
                    printlog ("                 grism == " // grism // ", \
                        filter == " // filter // " and slit == " // slit // \
                        ".", l_logfile, verbose+)
                    printlog ("                 Using default value \
                        (x2 = " // x2 // ").", l_logfile, verbose+)
                    ix2 = int(x2)
                } else {
                    ix2 = int(tabpar.value)
                }
            usetable = no
            }
        } else {
            printlog ("WARNING - F2CUT: No entry found in " // offsetdb // " \
                for grism == " // grism // ", filter == " // filter // " and \
                slit == " // slit // ".", l_logfile, verbose+)
            if (longslit) {
                x1 = 30
                x2 = 1525
                ix1 = int(x1)
        		ix2 = int(x2)
                printlog ("                 Using x1 = " // x1 // ", x2 = " \
                    // x2 // " and full y range to cut longslit data.", \
                    l_logfile, verbose+)
            } else {
                printlog ("                 Using full y range.", l_logfile, \
                    verbose+)
            }
            yoffset = INDEF
        }
        tdelete (tmptab, yes, verify-, >& "dev$null")

        # If there is no yoffset available, no need to get the filter width
        # or dispersion, since the full y range will be used
        if (!isindef(yoffset)) {

            # Get the filter width from the lookup table
            tselect (filterdb, tmptab, "filter == '" // filter // "'")

            tinfo (tmptab, ttout-)
            if (tinfo.nrows == 1)
                usetable = yes
            if (debug) print (tinfo.nrows // ", " // usetable)

            if (usetable) {
                # Use the 50% cut-on and cut-off wavelengths to define the
                # section to cut
                tabpar (table=tmptab, col="cuton50", row=1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2CUT: No value for filter width \
                        found in " // filterdb // ".", l_logfile, verbose+)
                    goto clean
                } else {
                    filter_lower = real(tabpar.value)
                }
                tabpar (table=tmptab, col="cutoff50", row=1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2CUT: No value for filter width \
                        found in " // filterdb // ".", l_logfile, verbose+)
                    goto clean
                } else {
                    filter_upper = real(tabpar.value)
                }
                usetable = no
                filter_width = filter_upper - filter_lower
            } else {
                printlog ("ERROR - F2CUT: No value for filter width found \
                    in " // filterdb // " for filter == " // filter // ".", \
                    l_logfile, verbose+)
                goto clean
            }
            tdelete (tmptab, yes, verify-, >& "dev$null")

            # Get the dispersion from the lookup table
            tselect (nsappwavedb, tmptab, "grating == '" // grism // "' && \
                filter == '" // filter // "'")

            tinfo (table=tmptab, >& "dev$null")
            if (tinfo.nrows > 1)
                usetable = yes

            if (debug) print (tinfo.nrows // ", " // usetable)

            if (usetable) {
                tabpar (table=tmptab, col="DELTA", row=1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2CUT: No dispersion value found \
                        in " // nsappwavedb // ".", l_logfile, verbose+)
                    goto clean
                } else {
                    # The dispersion is required in microns, but the value
                    # in the table is in (negative) Angstroms
                    orig_dispersion = real(tabpar.value)
                    dispersion = orig_dispersion / -10000
                }
                usetable = no
            } else {
                printlog ("ERROR - F2CUT: No dispersion value found in \
                    " // nsappwavedb // " for grism == " // grism // ", \
                    filter == " // filter // " and slit == " // slit // \
                    ".", l_logfile, verbose+)
                goto clean
            }
            tdelete (tmptab, yes, verify-, >& "dev$null")
        } else {
	    orig_dispersion = INDEF
            dispersion = INDEF
        }
 
        if (debug) {
            print ("yoffset: " // str(yoffset))
            print ("filter_width: " // str(filter_width))
            print ("dispersion: " // str(dispersion))
        }

        if (debug) {
            print ("TIMESTAMP: finished obtaining header information")
            date "+%H:%M:%S.%N"
        }

        # Determine the number of slits to cut; for longslit data, only one
        # slit is cut
        if (longslit) {
            slits = 1

            # Create the output image
            imcopy (input // phu, output, verbose-, >& "dev$null")
        }

        if (mos) {
            tinfo (input // mdf, ttout-)
            slits = tinfo.nrows

            # Copy the MDF to the output image
            wmef (input // mdf, output, extnames="MDF", phu=input // phu, \
                verbose-, >& "dev$null")
            if (wmef.status != 0) {
                printlog ("ERROR - F2CUT: Could not write final MEF file \
                    (WMEF).", l_logfile, verbose+)
                goto clean
            }

            # Add section columns to the MDF
            tcalc (output // mdf, "SECX1", 1, datatype="int")
            tcalc (output // mdf, "SECX2", 1, datatype="int")
            tcalc (output // mdf, "SECY1", 1, datatype="int")
            tcalc (output // mdf, "SECY2", 1, datatype="int")
            tcalc (output // mdf, "EXTVER", 0, datatype="int")

	    # Do some prep work from the MDF
	    if (usegrad) {

                # Determine the minimum slit length from the MDF (corresponding
                # to the minimum spectral width) for use with apfind (if minsep
                # = "default")
                tstat (input // mdf, column="slitsize_mx", outtable="", \
                    lowlim=INDEF, highlim=INDEF, rows="-", n_tab="table", \
                    n_nam="column", n_nrows="nrows", n_mean="mean", \
                    n_stddev="stddev", n_median="median", n_min="min", \
                    n_max="max")
                min_width = tstat.vmin

                # Use the minimum slit width (by default) from the MDF to
                # determine the value to use for apfind.minsep (the minimum
                # separation between spectra)
                if (l_minsep == "default") {
                    apfind_minsep = min_width * asecmm / pixscale
                } else {
                    apfind_minsep = real(l_minsep)
                }
	    }
        }

        print (output, >> tmptodel)

        if (debug) print ("number of slits: " // slits)

        printlog (input // " --> " // output, l_logfile, l_verbose)

        extver = 0
        for (slit_num = 1; slit_num <= slits; slit_num += 1) {

            dolefttrace = yes
            dorighttrace = yes
            slit_offset_left = INDEF
            slit_offset_right = INDEF

            if (usegrad) {

                # Retrieve slit information from MDF
                tabpar (table=input // mdf, col="slitpos_mx", row=slit_num, \
                    format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2CUT: No data found in " // input // \
                        mdf, l_logfile, verbose+)
                    goto clean
                } else {
                    slitpos_mx = real (tabpar.value)
                }

                tabpar (table=input // mdf, col="slitpos_my", row=slit_num, 
                    format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2CUT: No data found in " // input // \
                        mdf, l_logfile, verbose+)
                    goto clean
                } else {
                    slitpos_my = real (tabpar.value)
                }

                tabpar (table=input // mdf, col="slitsize_mx", row=slit_num, 
                    format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2CUT: No data found in " // input // \
                        mdf, l_logfile, verbose+)
                    goto clean
                } else {
                    slitsize_mx = real (tabpar.value)
                }

                tabpar (table=input // mdf, col="slitsize_my", row=slit_num, 
                    format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2CUT: No data found in " // input // \
                        mdf, l_logfile, verbose+)
                    goto clean
                } else {
                    slitsize_my = real (tabpar.value)
                }

                tabpar (table=input // mdf, col="priority", row=slit_num, 
                    format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2CUT: No data found in " // input // \
                        mdf, l_logfile, verbose+)
                    goto clean
                } else {
                    priority = int (tabpar.value)
                }

                if (debug)
                    print ("Information for slit " // slit_num // " retrieved")

                if (priority != 0) {

                    # Apply a distortion correction to slitpos_mx so that,
                    # after conversion from the mask plane to the array, it
                    # matches with the actual centre of the spectrum as seen on
                    # the array. The distortion correction was determined using
                    # the x pixel position of the left slit edge on the array
                    # at a y pixel position of y_ccd (determined below) and
                    # comparing that value with x_ccd. Both slitpos_mx and
                    # x_ccd are in mm.
                    x_ccd = (slitpos_mx * 1.00035894) - 0.86698091

		    # Based on June 2017 data
		    #x_ccd = (slitpos_mx * 0.99548) + 3.75662

                    # Convert x_ccd and y_ccd from mm to pixels. asecmm is
                    # measured in arcsec/mm and pixscale is measured in
                    # arcsec/pixel
                    x_ccd = (x_ccd * asecmm / pixscale) + (nx / 2)
                    y_ccd = (slitpos_my * asecmm / pixscale) + (ny / 2)

                    if (debug)
                        printf ("centre of slit %d from MDF: %.2f, %.2f\n",
                            slit_num, x_ccd, y_ccd)

                    slitsize_x = slitsize_mx * asecmm / pixscale
                    slitsize_y = slitsize_my * asecmm / pixscale

                    mdf_left_edge = x_ccd - (slitsize_x / 2)
                    mdf_right_edge = x_ccd + (slitsize_x / 2)

                    y1 = y_ccd + yoffset - (filter_width / dispersion) / 2
                    y2 = y_ccd + yoffset + (filter_width / dispersion) / 2

                    dy = y2 - y1
                    y_cen = y1 + (dy / 2)
                    iy1 = int(y1)
                    iy2 = int(y2)

                    if (debug) {
                        print ("TIMESTAMP: calling apfind for left slit edges")
                        date "+%H:%M:%S.%N"
                    }

                    # Locate the slit edges in the gradient image
                    # Find the left slit edges. Since the expected position of
                    # the slit after applying the distortion correction should
                    # be correct at y_ccd, detect the edges at this y pixel
                    # position
                    if (debug) printf ("detect edges at y = %.2f\n", y_cen)
                    apfind (tmpgradpos, nfind=slits, apertures="", \
                        references="", interactive=seeplot, find+, recenter-, \
                        resize-, edit+, line=y_cen, nsum=20., \
                        minsep=apfind_minsep, maxsep=1000., \
                        order="increasing", >& "dev$null")
                    match ("center", l_database // "/ap"//tmpgradpos, stop-, \
                        > tmpcenter)
                    match ("begin", l_database // "/ap"//tmpgradpos, \
                        stop+) | match ("aperture", "STDIN", stop-, > tmpapnum)
                    tmerge (tmpcenter // "," // tmpapnum, tmpgradleft, \
                        option="merge", allcols+, tbltype="default", \
                        allrows=1, extracol=0)
                    delete (tmpcenter // "," // tmpapnum, verify-, \
                        >& "dev$null")

                    # tmpgradleft is a file that contains three columns
                    # ("center" x_edge y_edge) for each left slit edge found in
                    # the gradient image. Find the left slit edge found using
                    # the gradient image (the second column in tmpgradleft)
                    # that matches with the slit slit_num in the MDF 
                    print (mdf_left_edge, > tmpmdfleft)
                    if (debug)
                        printf ("left slit edge from MDF: %.2f\n", 
                            mdf_left_edge)
                    tolerance = 0.75 * slitsize_x
                    tmatch (tmpgradleft, tmpmdfleft, tmptmatch, "c2", "c1", \
                        tolerance, >& "dev$null")

                    if (debug) tprint (tmptmatch)
                    tinfo (tmptmatch, ttout-)

                    if (tinfo.nrows > 0) {
                        # tmptmatch is a table containing seven columns (row,
                        # c1_1, c2, c3, c4, c5, c1_2), where c1_1 is "center",
                        # c2 is the x_edge (from tmpgradleft) of the left slit
                        # edge found using the gradient image, c3 is the y_edge
                        # (from tmpgradleft), c4 is "aperture", c5 is the
                        # aperture number and c1_2 is the left slit edge from
                        # the MDF. Find the left slit edge (c2) that is closest
                        # to the slit in the MDF (c1_2) 
                        tcalc (tmptmatch, "c6", "abs(c2-c1_2)", \
                            datatype="real", colunits="", colfmt="f6.3")
                        tcalc (tmptmatch, "c7", "c2-c1_2", datatype="real", \
                            colunits="", colfmt="f6.3")
                        tsort (tmptmatch, "c6", ascend+, casesens+)

                        tabpar (table=tmptmatch, column="c2", row=1, format-)
                        if (tabpar.undef) {
                            printlog ("ERROR - F2CUT: No data found in " // \
                                tmptmatchmdf, l_logfile, verbose+)
                            goto clean
                        } else {
                            grad_left_edge = real (tabpar.value)
                            if (debug)
                                printf ("left slit edge from grad image: \
                                    %.2f\n", grad_left_edge)
                        }

                        tabpar (table=tmptmatch, column="c5", row=1, format-)
                        if (tabpar.undef) {
                            printlog ("ERROR - F2CUT: No data found in " // \
                                tmptmatchmdf, l_logfile, verbose+)
                            goto clean
                        } else {
                            appos = int (tabpar.value)
                            if (debug)
                                printf ("aperture: %d\n", appos)
                        }

                        tabpar (table=tmptmatch, column="c7", row=1, format-)
                        if (tabpar.undef) {
                            printlog ("ERROR - F2CUT: No data found in " // \
                                tmptmatchmdf, l_logfile, verbose+)
                            goto clean
                        } else {
                            slit_offset_left = real (tabpar.value)
                            if (debug)
                                printf ("gradimage edge - mdf edge: %.2f\n",
                                    (slit_offset_left))
                        }

                        final_left_edge = grad_left_edge

                    } else {
                        # If the left slit edge isn't found, wait and see if
                        # the right slit edge is found
                        printlog ("WARNING - F2CUT: No left slit edge found \
                            in gradimage within " // tolerance, l_logfile, \
                            verbose+)
                        printlog ("                 pixels of slit edge found \
                            using MDF.", l_logfile, verbose+)

                        noleftedge = yes
                        grad_left_edge = INDEF

                        # No left edge was found, so no tracing can be done
                        dolefttrace = no
                    }
                    delete (tmptmatch, verify-, >& "dev$null")

                    if (debug) {
                        print ("TIMESTAMP: calling apfind for right slit \
                            edges")
                        date "+%H:%M:%S.%N"
                    }

                    # Find the right slit edges
                    apfind (tmpgradneg, nfind=slits, apertures="", \
                        references="", interactive=seeplot, find+, recenter-, \
                        resize-, edit+, line=y_cen, nsum=20., \
                        minsep=apfind_minsep, maxsep=1000., \
                        order="increasing", >& "dev$null")
                    match ("center", l_database // "/ap"//tmpgradneg, stop-, \
                        > tmpcenter)
                    match ("begin", l_database // "/ap"//tmpgradneg, \
                        stop+) | match ("aperture", "STDIN", stop-, > tmpapnum)
                    tmerge (tmpcenter // "," // tmpapnum, tmpgradright, \
                        option="merge", allcols+, tbltype="default", \
                        allrows=1, extracol=0)
                    delete (tmpcenter // "," // tmpapnum, verify-, \
                        >& "dev$null")

                    # tmpgradright is a file that contains three columns
                    # ("center" x_edge y_edge) for each right slit edge found
                    # in the gradient image. Find the right slit edge found
                    # using the gradient image (the second column in
                    # tmpgradright) that matches with the slit slit_num in the
                    # MDF
                    print (mdf_right_edge, > tmpmdfright)
                    if (debug)
                        printf ("right slit edge from MDF: %.2f\n",
                            mdf_right_edge)
                    tmatch (tmpgradright, tmpmdfright, tmptmatch, "c2", "c1", \
                        tolerance, >& "dev$null")

                    if (debug) tprint (tmptmatch)
                    tinfo (tmptmatch, ttout-)

                    if (tinfo.nrows > 0) {

                        # tmptmatch is a table containing seven columns (row,
                        # c1_1, c2, c3, c4, c5, c1_2), where c1_1 is "center",
                        # c2 is the x_edge (from tmpgradright) of the right
                        # slit edge found using the gradient image, c3 is the
                        # y_edge (from tmpgradright), c4 is "aperture", c5 is
                        # the aperture number and c1_2 is the right slit edge
                        # from the MDF. Find the right slit edge (c2) that is 
                        # closest to the slit in the MDF (c1_2)
                        tcalc (tmptmatch, "c6", "abs(c2-c1_2)", \
                            datatype="real", colunits="", colfmt="f6.3")
                        tcalc (tmptmatch, "c7", "c2-c1_2", datatype="real", \
                            colunits="", colfmt="f6.3")
                        tsort (tmptmatch, "c6", ascend+, casesens+)

                        tabpar (table=tmptmatch, column="c2", row=1, format-)
                        if (tabpar.undef) {
                            printlog ("ERROR - F2CUT: No data found in " // \
                                tmptmatchmdf, l_logfile, verbose+)
                            goto clean
                        } else {
                            grad_right_edge = real (tabpar.value)
                            if (debug)
                                printf ("right slit edge from grad image: \
                                    %.2f\n", grad_right_edge)
                        }

                        tabpar (table=tmptmatch, column="c5", row=1, format-)
                        if (tabpar.undef) {
                            printlog ("ERROR - F2CUT: No data found in " // \
                                tmptmatchmdf, l_logfile, verbose+)
                            goto clean
                        } else {
                            apneg = int (tabpar.value)
                            if (debug)
                                printf ("aperture: %d\n", apneg)
                        }

                        tabpar (table=tmptmatch, column="c7", row=1, format-)
                        if (tabpar.undef) {
                            printlog ("ERROR - F2CUT: No data found in " // \
                                tmptmatchmdf, l_logfile, verbose+)
                            goto clean
                        } else {
                            slit_offset_right = real (tabpar.value)
                            if (debug)
                                printf ("gradimage edge - mdf edge: %.2f\n",
                                    (slit_offset_right))
                        }

                        final_right_edge = grad_right_edge

                    } else {
                        if (noleftedge) {
                            # Neither the left slit edge nor the right slit
                            # edge were found, so move on to the next slit.
                            printlog ("WARNING - F2CUT: No right slit edge \
                                found in gradimage within " // tolerance, \
                                l_logfile, verbose+)
                            printlog ("                 pixels of slit edge \
                                found using MDF. Skipping slit " // \
                                slit_num // ".", l_logfile, verbose+)
                            section = ""
                            delete (tmptmatch, verify-, >& "dev$null")
                            goto NEXTSLIT
                        } else {
                            # If the right slit edge isn't found, use the left
                            # slit edge plus the slit width
                            printlog ("WARNING - F2CUT: No right slit edge \
                                found in gradimage within " // tolerance, \
                                l_logfile, verbose+)
                            printlog ("                 pixels of slit edge \
                                found using MDF. Using left slit", \
                                l_logfile, verbose+)
                            printlog ("                  edge plus MDF slit \
                                width to get right slit edge.", l_logfile, \
                                verbose+)

                            grad_right_edge = INDEF
                            final_right_edge = final_left_edge + slitsize_x

                            # No right edge was found, so no tracing can be
                            # done 
                            dorighttrace = no
                        }
                    }
                    delete (tmptmatch, verify-, >& "dev$null")

                    # Check that the offset between the left slit edge found
                    # using the gradient image and the left slit edge
                    # determined from the MDF is within the number of pixels
                    # given by the parameter edgediff of the offset between the
                    # right slit edge found using the gradient image and the
                    # right slit edge determined from the MDF (i.e., make sure
                    # the edges found using the gradient image belong to the
                    # same slit). Typically, the difference between the offsets
                    # is less than 0.1 pixels. Troublesome slits will have a
                    # difference greater than this.
                    if ((!isindef(grad_left_edge)) && \
                        (!isindef(grad_right_edge))) {

                        if (abs (slit_offset_left - slit_offset_right) > \
                            l_edgediff) {

                            printlog ("WARNING - F2CUT: Left and right slit \
                                edges determined from gradimage are", \
                                l_logfile, l_verbose)
                            printlog ("                 not associated with \
                                the same slit within " // l_edgediff // " \
                                pixels", l_logfile, l_verbose)
                            printlog ("                 (difference = " // \
                                abs (slit_offset_left - slit_offset_right) // \
                                ").", l_logfile, l_verbose)

                            # Determine which slit edge found using the
                            # gradient image is closest to the slit edge
                            # determined from the MDF and use that slit edge
                            # for tracing
                            if (abs (slit_offset_left) < \
                                abs (slit_offset_right)) {

                                printlog ("                 Using left slit \
                                    edge plus MDF slit width to get right", \
                                    l_logfile, l_verbose)
                                printlog ("                 slit edge.", \
                                    l_logfile, l_verbose)

                                final_right_edge = final_left_edge + slitsize_x

                                dorighttrace = no
                            }

                            if (abs (slit_offset_right) < \
                                abs (slit_offset_left)) {

                                printlog ("                 Using right slit \
                                    edge minus MDF slit width to get left", \
                                    l_logfile, l_verbose)
                                printlog ("                 slit edge.", \
                                    l_logfile, l_verbose)

                                final_left_edge = final_right_edge - slitsize_x

                                dolefttrace = no
                            }
                        }
                    }

                    if (isindef(grad_left_edge) && !isindef(grad_right_edge)) {
                        printlog ("WARNING - F2CUT: Using right slit edge \
                            minus MDF slit width to get left", l_logfile, \
                            verbose+)
                        printlog ("                 slit edge.", l_logfile, \
                            verbose+)

                        final_left_edge = final_right_edge - slitsize_x
                    }

                    # To ensure that the cut section contains the entirety of
                    # each slit, trace the edges of the slits to determine the
                    # full extent of each slit (distortions means that simple
                    # assumptions such as the bottom (i.e., at low y) of each
                    # slit have the lowest x are not true). Use order=3 to fit
                    # the slits appropriately. Set resize+ to update the
                    # aperture used for each peak in the database file. First,
                    # trace the left slit edge.
                    if (dolefttrace) {

                        if (debug) {
                            print ("TIMESTAMP: calling aptrace for left slit \
                                edge")
                            date "+%H:%M:%S.%N"
                        }

                        aptrace (tmpgradpos, apertures=appos, \
                            references=tmpgradpos, interactive=seeplot, \
                            find-, recenter-, resize+, edit-, trace+, \
                            fittrace-, line=y_cen, nsum=20, step=10, nlost=3, \
                            function="legendre", order=3, \
                            sample=iy1 // ":" // iy2, naverage=1, niterate=0, \
                            low_reject=3., high_reject=3., grow=0., \
                            >& "dev$null")

                        # Determine the minimum and maximum x pixel value from
                        # the fit from aptrace. From the help file for aptrace,
                        # the following equation is used for an order 3
                        # legendre (i.e., i = 3):
                        #   x = (coeff1 * z_1)+(coeff2 * z_2)+(coeff3 * z_3)
                        # where:
                        #   z_1 = 1
                        #   z_2 = n
                        #   z_3 = ((2 * 3 - 3) * n * z_2 - 1 * z_1) / 2
                        #       = (3 * n * n - 1) / 2
                        # and:
                        #   n = (2 * y - (ymax + ymin)) / (ymax - ymin)
                        # So:
                        #   x = (coeff1)+(coeff2 * n)+(coeff3 * (3n^2 - 1)/2)

                        if (debug) {
                            print ("TIMESTAMP: parse database file")
                            date "+%H:%M:%S.%N"
                        }

                        # First match the correct entry in the database file.
                        # Use the aperture number determined from the apfind
                        # match.
                        scanfile3 = l_database // "/ap"//tmpgradpos
                        count = 0
                        while (fscan (scanfile3, word, tmpvalue) != EOF) {
                            count += 1
                            #if (debug) print ("word: " // word)
                            if (word == "aperture") {
                                if (tmpvalue == appos) {
                                    if (debug) {
                                        print ("aperture: " // tmpvalue)
                                    }
                                    break
                                }
                            }
                        }
                        scanfile3 = ""

                        if (debug) {
                            print ("line number: " // count)
                            print ("TIMESTAMP: parse updated database file")
                            date "+%H:%M:%S.%N"
                        }

                        # Extract the database information for the current slit
                        # edge. The "aperture" line is line 3 of the 27 lines
                        # that need to be extracted 
                        tail (l_database // "/ap"//tmpgradpos, \
                            nlines=-(count-3)) | head ("STDIN", \
                            nlines=27, > tmpdb)
                        if (debug) type (tmpdb)

                        # Get the appropriate information for the current slit
                        scanfile3 = tmpdb
                        getvalues = no
                        count = -1

                        # Reset the params array before setting new values
                        for (i = 1; i <= 7; i += 1) {
                            params[i] = 0
                        }

                        while (fscan (scanfile3, word, tmpvalue) != EOF) {
                            if (word == "low")
                                aplow = tmpvalue
                            if (word == "high")
                                aphigh = tmpvalue
                            if (word == "curve")
                                getvalues = yes
                            if (getvalues) {
                                # Skip the next entry (that's the word
                                # "curve"), then retrieve the next 7 entries:
                                # 2. (function, 2 = Legendre), 3. (order),
                                # ymin, ymax, coeff1, coeff2, coeff3
                                count += 1
                                if (count > 0) {
                                    if (debug) print ("word: " // word)
                                    params[count] = word
                                }
                            }
                        }
                        scanfile3 = ""
                        delete (tmpdb, verify-, >& "dev$null")

                        func = params[1]
                        ord = params[2]
                        ymin = real(params[3])
                        ymax = real(params[4])
                        coeff1 = real(params[5])
                        coeff2 = real(params[6])
                        coeff3 = real(params[7])

                        # Actually use y1 and y2 for ymin and ymax
			# Can't do this because equation won't be correct
                        #ymin = iy1
                        #ymax = iy2

                        if (debug)
                            print (aplow // " " // aphigh // " " // func // \
                                " " // ord // " " // ymin // " " // ymax // \
                                " " // coeff1 // " " // coeff2 // " " // \
                                coeff3)

			xmin = 10
			xmax = -10
			for (iy = iy1; iy <= iy2; iy +=1 ){
			    n = (2*iy - (ymax+ymin)) / (ymax-ymin)
			    edge[iy] = coeff1 + coeff2*n + 0.5*coeff3*(3*n*n-1)
			    if (edge[iy] < xmin) xmin = edge[iy]
			    if (edge[iy] > xmax) xmax = edge[iy]
			}

                        if (debug) {
                            print ("xmin: " // xmin)
                            print ("xmax: " // xmax)
                        }

                        # xmin and xmax define the extent of the *peaks* of the
                        # spectrum. Since the whole peak needs to be included
                        # in the slit image (so that the s-distortion
                        # correction can be reliably determined), use the width
                        # of the peak from the database file and include it in
                        # the calculation of x1 and x2. In the case that the
                        # right slit edge cannot be traced, make the assumption
                        # that the right edge trace will be identical to the
                        # left edge trace.
                        x1 = final_left_edge + xmin + aplow
                        ix1 = int (x1 + 0.5)

                        # Use the equation defining the slit to create a bad
                        # pixel file that can be used later to blank out the
                        # areas beyond the edge of the slit (i.e., mask other
                        # slits that may be within the cut area). The bad pixel
                        # file consists of lines describing rectangular regions
                        # of the image in the format x1 x2 y1 y2
                        if (debug) {
                            print ("TIMESTAMP: create bpm pixel file for left \
                                slit edge")
                            date "+%H:%M:%S.%N"
                        }

                        if (debug) {
                            print ("final left edge = " // final_left_edge)
                            print ("aplow = " // aplow)
                            print ("ix1 = " // ix1)
                            print ("iy1 = " // iy1)
                        }

                        # The following entry in the bpm pixel file will mask
                        # everything to the left of the left slit edge in the
                        # cut slit frame
			for (iy = iy1; iy <= iy2; iy += 1) {
			    x2 = edge[iy] + final_left_edge + aplow - ix1 + 1
			    print (1, int(x2), (iy-iy1+1), (iy-iy1+1),
			        >> tmpbpmfile)
			}
                    }

                    if (!dorighttrace) {
                        printlog ("WARNING - F2CUT: Assuming right edge trace \
                            will be identical to left edge", l_logfile, \
                            verbose+)
                        printlog ("                 trace.", l_logfile, \
                            verbose+)
                        x2 = final_right_edge + xmax + aphigh
                        ix2 = int (x2 + 0.5)

                        if (debug) {
                            print ("TIMESTAMP: update bpm pixel file for \
                                right slit edge")
                            date "+%H:%M:%S.%N"
                        }

                        if (debug) {
                            print ("final right edge = " // final_right_edge)
                            print ("aphigh = " // aphigh)
                            print ("ix1 = " // ix1)
                            print ("iy1 = " // iy1)
                        }

                        # The following entry in the bpm pixel file will mask
                        # everything to the right of the right slit edge in the
                        # cut slit frame 
			for (iy = iy1; iy <= iy2; iy += 1) {
			    x1 = edge[iy] + final_right_edge + aphigh - ix1 + 1
			    print (int(x1), (ix2-ix1+1), (iy-iy1+1), (iy-iy1+1), \
			        >> tmpbpmfile)
			}
                    }

                    # Now trace the right slit edge
                    if (dorighttrace) {

                        if (debug) {
                            print ("TIMESTAMP: calling aptrace for right slit \
                                edge")
                            date "+%H:%M:%S.%N"
                        }

                        aptrace (tmpgradneg, apertures=apneg, \
                            references=tmpgradneg, interactive=seeplot, \
                            find-, recenter-, resize+, edit-, trace+, \
                            fittrace-, line=y_cen, nsum=20, step=10, nlost=3, \
                            function="legendre", order=3, \
                            sample=iy1 // ":" // iy2, naverage=1, niterate=0, \
                            low_reject=3., high_reject=3., grow=0., \
                            >& "dev$null")

                        if (debug) {
                            print ("TIMESTAMP: parse database file")
                            date "+%H:%M:%S.%N"
                        }

                        # First match the correct entry in the database file.
                        # Use the aperture number determined from the apfind
                        # match
                        scanfile3 = l_database // "/ap"//tmpgradneg
                        count = 0
                        while (fscan (scanfile3, word, tmpvalue) != EOF) {
                            count += 1
                            #if (debug) print ("word: " // word)
                            if (word == "aperture") {
                                if (tmpvalue == apneg) {
                                    if (debug) {
                                        print ("aperture: " // tmpvalue)
                                    }
                                    break
                                }
                            }
                        }
                        scanfile3 = ""

                        if (debug) {
                            print ("line number: " // count)
                            print ("TIMESTAMP: parse updated database file")
                            date "+%H:%M:%S.%N"
                        }

                        # Extract the database information for the current slit
                        # edge. The "aperture" line is line 3 of the 27 lines
                        # that need to be extracted 
                        tail (l_database // "/ap"//tmpgradneg, \
                            nlines=-(count-3)) | head ("STDIN", \
                            nlines=27, > tmpdb)
                        if (debug) type (tmpdb)

                        # Get the appropriate information for the current slit
                        scanfile3 = tmpdb
                        getvalues = no
                        count = -1

                        # Reset the params array before setting new values
                        for (i = 1; i <= 7; i += 1) {
                            params[i] = 0
                        }

                        while (fscan (scanfile3, word, tmpvalue) != EOF) {
                            if (word == "low")
                                aplow = tmpvalue
                            if (word == "high")
                                aphigh = tmpvalue
                            if (word == "curve")
                                getvalues = yes
                            if (getvalues) {
                                # Skip the next entry (that's the word 
                                # "curve"), then retrieve the next 7 entries:
                                # 2. (function, 2 = Legendre), 3. (order), 
                                # ymin, ymax, coeff1, coeff2, coeff3
                                count += 1
                                if (count > 0) {
                                    if (debug) print ("word: " // word)
                                    params[count] = word
                                }
                            }
                        }
                        scanfile3 = ""
                        delete (tmpdb, verify-, >& "dev$null")

                        func = params[1]
                        ord = params[2]
                        ymin = real(params[3])
                        ymax = real(params[4])
                        coeff1 = real(params[5])
                        coeff2 = real(params[6])
                        coeff3 = real(params[7])

                        if (debug)
                            print (aplow // " " // aphigh // " " // func // \
                                " " // ord // " " // ymin // " " // ymax // \
                                " " // coeff1 // " " // coeff2 // " " // \
                                coeff3)

			xmin = 10
			xmax = -10
			for (iy = iy1; iy <= iy2; iy +=1 ){
			    n = (2*iy - (ymax+ymin)) / (ymax-ymin)
			    edge[iy] = coeff1 + coeff2*n + 0.5*coeff3*(3*n*n-1)
			    if (edge[iy] < xmin) xmin = edge[iy]
			    if (edge[iy] > xmax) xmax = edge[iy]
			}

                        if (debug) {
                            print ("xmin: " // xmin)
                            print ("xmax: " // xmax)
                        }

                        # xmin and xmax define the extent of the *peaks* of the
                        # spectrum. Since the whole peak needs to be included
                        # in the slit image (so that the s-distortion
                        # correction can be reliably determined), use the width
                        # of the peak from the database file and include it in
                        # the calculation of x1 and x2. In the case that the
                        # left slit edge cannot be traced, make the assumption
                        # that the left edge trace will be identical to the
                        # right edge trace.
                        x2 = final_right_edge + xmax + aphigh
                        ix2 = int (x2 + 0.5)

			# Need to define left edge of cut slit if we couldn't
			# trace the left edge
			if (!dolefttrace) {
			   x1 = final_left_edge + xmin + aplow
			   ix1 = int (x1 + 0.5)
			}

                        # Use the equation defining the slit to create a bad
                        # pixel file that can be used later to blank out the
                        # areas beyond the edge of the slit (i.e., mask other
                        # slits that may be within the cut area). The bad pixel
                        # file consists of lines describing rectangular regions
                        # of the image in the format x1 x2 y1 y2
                        if (debug) {
                            print ("TIMESTAMP: update bpm pixel file for \
                                right slit edge")
                            date "+%H:%M:%S.%N"
                        }

                        if (debug) {
                            print ("final right edge = " // final_right_edge)
                            print ("aphigh = " // aphigh)
                            print ("ix1 = " // ix1)
                            print ("iy1 = " // iy1)
                        }

                        # The following entry in the bpm pixel file will mask
                        # everything to the right of the right slit edge in the
                        # cut slit frame 
			for (iy = iy1; iy <= iy2; iy += 1) {
			    x1 = edge[iy]+final_right_edge+aphigh-ix1+1
			    print (int(x1), (ix2-ix1+1), (iy-iy1+1), (iy-iy1+1), \
			        >> tmpbpmfile)
			}
                    }

                    if (!dolefttrace) {
                        printlog ("WARNING - F2CUT: Assuming left edge trace \
                            will be identical to right edge", l_logfile, \
                            verbose+)
                        printlog ("                 trace.", l_logfile, \
                            verbose+)

                        if (debug) {
                            print ("TIMESTAMP: update bpm pixel file for left \
                                slit edge")
                            date "+%H:%M:%S.%N"
                        }

                        if (debug) {
                            print ("final left edge = " // final_left_edge)
                            print ("aplow = " // aplow)
                            print ("ix1 = " // ix1)
                            print ("iy1 = " // iy1)
                        }

                        # The following entry in the bpm pixel file will mask
                        # everything to the left of the left slit edge in the
                        # cut slit frame 
			for (iy = iy1; iy <= iy2; iy += 1) {
			    x2 = edge[iy]+final_left_edge+aplow-ix1+1
			    print (1, int(x2), (iy-iy1+1), (iy-iy1+1), \
			        >> tmpbpmfile)
			}
                    }

                    if (debug) {
                        print ("TIMESTAMP: determine section")
                        date "+%H:%M:%S.%N"
                    }

                    # Determine the value of the slit edges in the cut image
                    # and write them to the header.
                    lsltedge = final_left_edge - ix1 + 1
                    rsltedge = final_right_edge - ix1 + 1

                    if (debug) {
                        print ("left slit edge in cut image: " // \
                            lsltedge)
                        print ("right slit edge in cut image: " // \
                            rsltedge)
                    }

		    if (ix1 < 1) {
                        printlog ("WARNING - F2CUT: Left edge beyond detector \
				 limit.", l_logfile, verbose+)
		        ix1 = 1
		    }
		    if (ix2 > nx) {
                        printlog ("WARNING - F2CUT: Right edge beyond detector \
				 limit.", l_logfile, verbose+)
                        ix2 = nx
		    }
                    section = "[" // ix1 // ":" // ix2 // "," \
                        // iy1 // ":" // iy2 // "]"

                    # Put image section info in MDF columns
		    extver += 1
                    partab (ix1, output // mdf, "SECX1", slit_num)
                    partab (ix2, output // mdf, "SECX2", slit_num)
                    partab (iy1, output // mdf, "SECY1", slit_num)
                    partab (iy2, output // mdf, "SECY2", slit_num)
                    partab (extver, output // mdf, "EXTVER", slit_num)

                } else {
                    section = ""
                }

NEXTSLIT:
                if (debug) print ("section: " // section)

                # Delete the temporary files used in this loop
                delete (tmpgradleft // "," // tmpmdfleft // "," // \
                    tmpgradright // "," // tmpmdfright, verify-, \
                    >& "dev$null")
            }

            if (useref) {

                tabpar (table=tmprefim // mdf, col="priority", \
                    row=slit_num, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2CUT: No data found in " // \
                        tmprefim, l_logfile, verbose+)
                    goto clean
                } else {
                    priority = int (tabpar.value)
                }

                tabpar (table=tmprefim // mdf, col="EXTVER", \
                    row=slit_num, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - F2CUT: No data found in " // \
                        tmprefim, l_logfile, verbose+)
                    goto clean
                } else {
                    extver = int (tabpar.value)
                }

                if (priority != 0 && extver != 0) {

                    # Obtain the section to cut from the MDF
                    tabpar (table=tmprefim // mdf, col="SECX1", \
                        row=slit_num, format-)
                    if (tabpar.undef) {
                        printlog ("ERROR - F2CUT: No data found in " // \
                            tmprefim, l_logfile, verbose+)
                        goto clean
                    } else {
                        ix1 = int (tabpar.value)
                    }

                    tabpar (table=tmprefim // mdf, col="SECX2", \
                        row=slit_num, format-)
                    if (tabpar.undef) {
                        printlog ("ERROR - F2CUT: No data found in " // \
                            tmprefim, l_logfile, verbose+)
                        goto clean
                    } else {
                        ix2 = int (tabpar.value)
                    }

                    tabpar (table=tmprefim // mdf, col="SECY1", \
                        row=slit_num, format-)
                    if (tabpar.undef) {
                        printlog ("ERROR - F2CUT: No data found in " // \
                            tmprefim, l_logfile, verbose+)
                        goto clean
                    } else {
                        iy1 = int (tabpar.value)
                    }

                    tabpar (table=tmprefim // mdf, col="SECY2", \
                        row=slit_num, format-)
                    if (tabpar.undef) {
                        printlog ("ERROR - F2CUT: No data found in " // \
                            tmprefim, l_logfile, verbose+)
                        goto clean
                    } else {
                        iy2 = int (tabpar.value)
                    }

                    if (debug)
                        print ("Information for slit " // slit_num // " \
                            retrieved")

                    section = "[" // ix1 // ":" // ix2 // "," \
                        // iy1 // ":" // iy2 // "]"

                    if (maskoverlap) {
                        # Obtain the appropriate DQ extension from the
                        # reference image so that the areas contaminated by
                        # overlap can be masked.
                        tmpbpm = tmprefim // "[" // l_dq_ext // "," // \
                            extver // "]"

                    }

                    # Put image section info in MDF columns
                    partab (ix1, output // mdf, "SECX1", slit_num)
                    partab (ix2, output // mdf, "SECX2", slit_num)
                    partab (iy1, output // mdf, "SECY1", slit_num)
                    partab (iy2, output // mdf, "SECY2", slit_num)
                    partab (extver, output // mdf, "EXTVER", slit_num)

                } else {
                    section = ""
                }

                if (debug) print ("section: " // section)
            }

            if (longslit) {

                extver += 1

                # If yoffset is INDEF, use the full y range
                if (isindef(yoffset)) {
                    y1 = 1
                    y2 = 2048
                } else {
                    slitpos_y = 1024
                    y1 = slitpos_y + yoffset - (filter_width / dispersion) / 2
                    y2 = slitpos_y + yoffset + (filter_width / dispersion) / 2
                }

                if (y1 < 1)
                    y1 = 1
                if (y2 > ny)
                    y2 = ny

                iy1 = int(y1)
                iy2 = int(y2)

                section = "[" // ix1 // ":" // ix2 // "," \
                    // iy1 // ":" // iy2 // "]"

                if (debug) print ("section: " // section)
            }

            # Copy the section to the output file
            if (section != "") {

                if (debug) {
                    print ("TIMESTAMP: create output file")
                    date "+%H:%M:%S.%N"
                }

                outsci = "[" // l_sci_ext // "," // extver // "]"

                if (first)
                    printlog (" " // sci // section // " --> " // outsci, \
                        l_logfile, l_verbose)

                if (l_fl_vardq) {
                    imcopy (input // dq // section, tmpdq, verbose-, \
                        >& "dev$null")
                }

                if (maketmpbpm) {

                    # Use the information in the bad pixel file to create a
                    # temporary bpm file. Set the value of the bad pixels equal
                    # to 32 to indicate that they are contaminated by overlap.
                    badpiximage (fixfile=tmpbpmfile, template=tmpdq, \
                        image=tmpbpm, goodvalue=0, badvalue=32)
                    delete (tmpbpmfile, verify-, >& "dev$null")
                }

                if (debug) {
                    print ("TIMESTAMP: create sci ext")
                    date "+%H:%M:%S.%N"
                }

                if (maskoverlap) {

                    # Replace the pixels in both the science extension and the
                    # variance extension in the output image that are
                    # contaminated by overlap (i.e., have a value of 32 in the 
                    # temporary bpm file) with a value of 0 and update the DQ
                    # extension.
                    imexpr ("b == 32 ? 0 : a", output // "[" // l_sci_ext // \
                        "," // extver // ",append]", input // sci // section, \
                        tmpbpm, dims="auto", intype="auto", outtype="auto", \
                        refim="auto", bwidth=0, btype="nearest", bpixval=0.0, \
                        rangecheck+, verbose-, exprdb="none")

                    if (l_fl_vardq) {
                        imexpr ("b == 32 ? 0 : a", output // "[" // \
                            l_var_ext // "," // extver // ",append]", input \
                            // var // section, tmpbpm, dims="auto", \
                            intype="auto", outtype="auto", refim="auto", \
                            bwidth=0, btype="nearest", bpixval=0.0, \
                            rangecheck+, verbose-, exprdb="none")
                    }
                } else {
                    imcopy (input // sci // section, output // "[" // \
                        l_sci_ext // "," // extver // ",append]", verbose-, \
                        >& "dev$null")

                    if (l_fl_vardq) {
                        imcopy (input // var // section, output // "[" // \
                            l_var_ext // "," // extver // ",append]", \
                            verbose-, >& "dev$null")
                    }
                }

                if (debug) {
                    print ("TIMESTAMP: update sci ext headers")
                    date "+%H:%M:%S.%N"
                }

                # Update headers
                gemhedit (output // outsci, l_key_cut_section, section, \
                    "Region extracted by F2CUT", delete-)
                gemhedit (output // outsci, l_key_order, 1, "Spectral order", \
                    delete-)

                if (l_fl_vardq) {
                    if (updatedq) {
                        imexpr ("b == 32 ? b : a", output // "[" // l_dq_ext \
                            // "," // extver // ",append]", tmpdq, tmpbpm, \
                            dims="auto", intype="ushort", outtype="ushort", \
                            refim="auto", bwidth=0, btype="nearest",
                            bpixval=0.0, rangecheck+, verbose-, exprdb="none")
                    } else {
                        imcopy (tmpdq, output // "[" // l_dq_ext // "," // \
                            extver // ",append]", verbose-, >& "dev$null")
                    }
                    imdelete (tmpdq, verify-, >& "dev$null")

                    # Update headers
                    gemhedit (output // "[" // l_var_ext // "," // extver // \
                        "]", l_key_cut_section, section, \
                        "Region extracted by F2CUT", delete-)
                    gemhedit (output // "[" // l_var_ext // "," // extver // \
                        "]", l_key_order, 1, "Spectral order", \
                        delete-)
                    gemhedit (output // "[" // l_dq_ext // "," // extver // \
                        "]", l_key_cut_section, section, \
                        "Region extracted by F2CUT", delete-)
                    gemhedit (output // "[" // l_dq_ext // "," // extver // \
                        "]", l_key_order, 1, "Spectral order", \
                        delete-)
                }

                # If a gradient image was used, create the output slit image,
                # which will be used to determine the s-distortion correction
                if (usegrad) {

                    # Create the output slit image
                    if (extver == 1) {
                        imcopy (tmpgradim // phu, l_outslitim, verbose-, \
                            >& "dev$null")
                    }

		    # Zero pixels outside this slit to avoid contamination
                    imexpr ("c == 32 ? 0 : a+b", l_outslitim // "[" // \
                        l_sci_ext // "," // extver // ",append]", \
			tmpgradpos // section, tmpgradneg // section, tmpbpm, \
			verbose-)

                    #imarith (tmpgradpos // section, "+", tmpgradneg // \
                    #    section, l_outslitim // "[" // l_sci_ext // "," // \
                    #    extver // ",append]", title="", divzero=0.0, \
                    #    hparams="", pixtype="", calctype="", verbose-, noact-)

                    # Update headers
                    gemhedit (l_outslitim // outsci, l_key_cut_section, \
                        section, "Region extracted by F2CUT", delete-)

                    gemhedit (l_outslitim // outsci, l_key_order, 1, \
                        "Spectral order", delete-)

                    gemhedit (l_outslitim // outsci, "LSLTEDGE", lsltedge, \
                        "x pixel position of the left slit edge", delete-)

                    gemhedit (l_outslitim // outsci, "RSLTEDGE", rsltedge, \
                        "x pixel position of the right slit edge", \
                        delete-)
                }
                imdelete (tmpbpm, verify-, >& "dev$null")

                if (debug) {
                    print ("TIMESTAMP: finished slit")
                    date "+%H:%M:%S.%N"
                }
            }
        }

        if (debug) {
            print ("TIMESTAMP: update headers")
            date "+%H:%M:%S.%N"
        }

        # Write the original image size to the header (required by nsflat) 
        gemhedit (output // phu, "ORIGXSIZ", nx, "Original size in X", delete-)
        gemhedit (output // phu, "ORIGYSIZ", ny, "Original size in Y", delete-)
        # Update headers
        gemhedit (output // phu, "NSCIEXT", extver, \
            "Number of science extensions", delete-)
        # Write the dispersion to the header for nsappwave
	if (!isindef(orig_dispersion)) {
            gemhedit (output // phu, l_key_delta, orig_dispersion, \
                "Dispersion written by F2CUT", delete-)
        }
        gemdate ()
        gemhedit (output // phu, "F2CUT", gemdate.outdate, \
            "UT Time stamp for F2CUT", delete-)
        gemhedit (output // phu, "NSCUT", gemdate.outdate, \
            "Dummy keyword written by F2CUT", delete-)
        gemhedit (output // phu, "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI", delete-)

        if (usegrad) {
            gemhedit (l_outslitim // phu, "NSCIEXT", extver, \
                "Number of science extensions", delete-)
            gemhedit (l_outslitim // phu, "F2CUT", gemdate.outdate, \
                "UT Time stamp for F2CUT", delete-)
            gemhedit (l_outslitim // phu, "NSCUT", gemdate.outdate, \
                "Dummy keyword written by F2CUT", delete-)
            gemhedit (l_outslitim // phu, "GEM-TLM", gemdate.outdate, \
                "UT Last modification with GEMINI", delete-)
        }
        first = no
    }

    # If we have managed to arrive here, everything has worked correctly
    status = 0

clean:
    # Clean up
    scanfile1 = ""
    scanfile2 = ""
    scanfile3 = ""

    if (status == 1) {
        if (access(tmptodel)) {
            # Delete the output files if the task failed
            imdelete ("@" // tmptodel, verify-, >& "dev$null")
        }
    }

    delete (tmpinlist // "," // tmpgradlist // "," // tmpreflist // "," // \
        tmpinput // "," // tmpoutput // "," // tmptab // "," // tmpcenter // \
        "," // tmpapnum // "," // tmpdb // "," // tmpgradleft // "," // \
        tmpmdfleft // "," // tmpgradright // "," // tmpmdfright // "," // \
        tmptmatch // "," // \
        tmpbpmfile // "," // tmptodel, verify-, >& "dev$null")
    imdelete (tmpgradpos // "," // tmpgradneg // "," // "," // tmpbpm // "," \
        // tmpdq, verify-, >& "dev$null")

    # Restore default specred values
    if (stored) {
        specred.logfile = l_logfilesave
        specred.database = l_databasesave
        specred.dispaxis = l_dispaxissave
        specred.verbose = l_verbosesave
    }

    printlog (" ", l_logfile, l_verbose)
    if (status == 0) {
        printlog ("F2CUT exit status: good.", l_logfile, l_verbose)
    } else {
        printlog ("F2CUT exit status: failed.", l_logfile, l_verbose)
    }

    printlog ("--------------------------------------------------------------\
        --------------", l_logfile, l_verbose)

end
