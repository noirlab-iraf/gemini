# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gscut (inimage)
    
    
# Cut GMOS MOS spectra into 2D pieces
#
# Version Feb 28, 2002  BM  v1.3 release
#         May 07, 2002  BM  add simple x distortion correction
#         Aug 15, 2002  BM  reintroduce x/yoffset which override
#                           database values if != INDEF, don't cut acq stars
#         Aug 18, 2002  BM  fix writing NSCIEXT
#         Aug 20, 2002  BM  removed call to GMOSAIC, inimage must be
#                           GMOSAICed if outimage desired
#         Aug 20, 2002  IJ  priority=1 for NS longslit that is handled as MOS
#         Aug 26, 2002  IJ  bugfix for longslit spectra, logfile for *offset=INDEF
#         Sept 20, 2002  IJ added w2 since this does in fact not work in the far red
#         Sept 20, 2002     v1.4 release
#         Oct 14, 2002  IJ  don't use j,i from within scripts
#         Feb 14, 2003  IJ  bugfix - had overlooked one of those i's
#         Mar 20, 2003  BM  pixscales for both instruments
#         May  9, 2003  IJ  change in instrument logic to support old GMOS-N data
#         May 24, 2003  BM  Remove y-distortion correction for GMOS-S
#         Aug 20, 2003  MB  added a parameter for extra slit width in pixels for MOS  
#         Aug 26, 2003  KL  IRAF2.12 - new parameter, addonly, in hedit
#                             DQ header, REFPIX1: 'hedit' replaced for 
#                                                 'gemhedit' as it should be. 
#         Feb 29, 2004  BM  Add support for two filters
#         Apr 08, 2004  BM  Fix problem in getting wmn2 for filter2
#         Sep 09, 2004  BM,IJ  Edge detection of slits, code cleaning
#         Mar 12, 2013  BM   Improved logic for the use of y2 from the gradient image
#         Jun 03, 2013  BM   Catch cases that caused crashes (ry1/2 not defined)

char    inimage     {prompt="Input image"}
char    outimage    {"",prompt="Output image"}
char    secfile     {"",prompt="Output file for image sections"}
bool    fl_update   {no,prompt="Update inimage header if no output image produced"}
bool    fl_vardq    {no,prompt="Propagate variance and data quality planes"}
bool    fl_oversize {yes,prompt="Use 1.05x slit length to accommodate distortion?"}
char    gratingdb   {"gmos$data/GMOSgratings.dat",prompt="Gratings database file"}
char    filterdb    {"gmos$data/GMOSfilters.dat",prompt="Filters database file"}
char    bpmfile     {"gmos$data/chipgaps.dat",prompt="Info on location of chip gaps"}
char    gradimage   {"",prompt="Image to use for finding slit edges using the gradient method"}
char    refimage    {"",prompt="Reference image for slit positions"}
real    xoffset     {INDEF,prompt="X offset in wavelength [nm]"}
real    yoffset     {INDEF,prompt="Y offset in unbinned pixels"}
real    yadd        {0,prompt="Additional unbinned pixels to add to each end of MOS slitlet lengths"}
real    w2          {INDEF,prompt="Upper wavelength limit of cut spectra (nm)"}
char    sci_ext     {"SCI",prompt="Name of science extension"}
char    var_ext     {"VAR",prompt="Name of variance extension"}
char    dq_ext      {"DQ",prompt="Name of data quality extension"}
char    key_gain    {"GAIN",prompt="Header keyword for gain (e-/ADU)"}   
char    key_ron     {"RDNOISE",prompt="Header keyword for readout noise"}  
real    ron         {3.5,prompt="Readout noise in electrons"}                   
real    gain        {2.2,prompt="Gain in e-/ADU"}                               
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose"}
int     status      {0,prompt="Exit status (0=good)"}

begin

    # local variables for task parameters
    char    l_inimage
    char    l_outimage = ""
    char    l_secfile = ""
    char    l_gradimage = "" 
    char    l_refimage = ""
    char    l_logfile = "" 
    char    l_key_gain = ""
    char    l_key_ron = ""
    char    l_gratingdb, l_filterdb, l_bpmfile
    char    l_sci_ext, l_var_ext, l_dq_ext
    real    l_yoff, l_uxoffset, l_uyoffset, l_uw2
    real    l_yadd, l_gain, l_ron
    bool    l_verbose, l_fl_vardq, l_fl_update, l_fl_oversize

    # other variables
    char  sec, datsec, detsec, mdf, impaste
    char  sci, var, dq, modobs
    char  grating, filter[2], gname, fname, ffile, outroot, slittype
    char  gradpos, gradneg, tmpcoo, tmppos, tmpneg, tmptmatch
    real    cwave, grule, gblaze, gR, gcoverage, nmppx, xscale, yscale
    real    gwave1, gwave2, fwave1, fwave2, wave1, wave2
    real    a, tilt, gtilt, greq, pi
    real    xccd, yccd, slitx, slity, slitwid, slitlen
    real    pixcwave, refpix, asecmm, speccen, apcol
    real    ssize_mx, ssize_my, ssize_mr, ssize_mw, stilt_m, spos_mx, spos_my
    real    wavoffset, dx, y, pixscale[2,3] # pixscale[inst,iccd]
    real    wmn1, wmn2, wmx1, wmx2, ry1, ry2
    int     detector_upper_spec_limit[3] # detector_upper_spec_limit[iccd]
    int     xbin, ybin, speclen, specwid, center, xcen, ycen, priority
    int     slitsepp, objid, nxpix, nypix, inst, iccd, apnsum, aphsum
    int     nslit, j1, id1, nr, x1, y1, x2, y2, n_j, nextnd, n_i, nsciext
    int		grorder
    bool    found, foundup, ishamamatsu, debug, fullframe
    struct  sdate
    int     s_dispaxis, junk, ngap, gapl[2], gaph[2]
    char    s_logfile, s_database, tmpoff, mdftest

    # Query parameters
    junk = fscan (inimage, l_inimage) 
    junk = fscan (outimage, l_outimage)
    junk = fscan (secfile, l_secfile)
    l_fl_update=fl_update
    l_fl_vardq = fl_vardq
    l_fl_oversize = fl_oversize
    junk = fscan (gratingdb, l_gratingdb)
    junk = fscan (filterdb, l_filterdb)
    junk = fscan (bpmfile, l_bpmfile)
    junk = fscan (gradimage, l_gradimage) 
    junk = fscan (refimage, l_refimage) 
    l_uxoffset = xoffset
    l_uyoffset = yoffset
    l_yadd = yadd
    l_uw2 = w2
    junk = fscan (sci_ext, l_sci_ext)
    junk = fscan (var_ext, l_var_ext)
    junk = fscan (dq_ext, l_dq_ext)
    junk = fscan (key_gain, l_key_gain)
    junk = fscan (key_ron, l_key_ron)
    l_ron = ron
    l_gain = gain
    junk = fscan (logfile, l_logfile)
    l_verbose=verbose

    # Parameter for debugging
    debug = no

    # Define some deafult / lookup values
    status = 0
    ishamamatsu = no
    pi = 3.14159265
    asecmm = 1.611444
    
    grorder = 1

    # pixelscale[inst,iccd]
    pixscale[1,1] = 0.0727 # GMOS-N EEV pixscale
    pixscale[1,2] = 0.07288 # GMOS-N e2vDD pixscale
    pixscale[1,3] = 0.0807 # GMOS-N Hamamatsu pixscale ##M
    pixscale[2,1] = 0.073 # GMOS-S EEV pixscale 
    pixscale[2,3] = 0.0800 # GMOS-S Hamamatsu pixscale ##M PIXEL_SCALE

    # Define the spectral cut-off limit (red limit) according to the iccd
    # (detector type). Value is in nm. If needed, this can be changed to also
    # accomodate diffrenet values for GMOS-N and GMOS-S, see pixscale
    detector_upper_spec_limit[1] = 1025 # EEV CCDs 
    detector_upper_spec_limit[2] = 1050 # e2vDD CCDs ##M
    detector_upper_spec_limit[3] = 1080 # Hamamatsu CCDs ##M
    
    # speccen: this is the default value for MOS with refimage="" and 
    # gradimage="".  I (KL) am not sure what it should be set to.
    speccen = 0
    mdftest = ""

    # Keep imgets parameters from changing by outside world
    cache ("imgets", "gemhedit", "tinfo", "fparse", "tabpar", "gmos")
    cache ("gimverify", "specred", "gemdate")
    
    # Save user parameters
    s_dispaxis = specred.dispaxis ; specred.dispaxis = 1
    s_logfile = specred.logfile ; specred.logfile = ""
    s_database = specred.database ; specred.database = "database"

    # Define temporary files
    mdf = mktemp("tmpmdf")
    impaste = mktemp("tmppaste")
    gradpos = mktemp("tmpgradpos") ; gradneg=mktemp("tmpgradneg")
    tmpoff = mktemp("tmpoff")

    # Assign variables to be used later
    sci = ""
    var = ""
    dq = ""
    tmppos = "" ; tmpneg = "" ; tmpcoo = "" ; tmptmatch = ""

    # Start logging to file
    if (l_logfile == "STDOUT") {
        l_logfile = ""
        l_verbose = yes
    }

    # Test the logfile:
    if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GSCUT: Both gscut.logfile and gmos.logfile \
                fields are empty", l_logfile, verbose+)
            printlog ("                    Using default file gmos.log",
                l_logfile, verbose+)
        }
    }

    date | scan(sdate)
    printlog ("-----------------------------------------------------------\
        ---------------------", l_logfile, verbose=l_verbose)
    printlog ("GSCUT -- "//sdate, l_logfile, verbose=l_verbose)
    printlog (" ", l_logfile, verbose=l_verbose)
    printlog ("inimage     = "//l_inimage, l_logfile, verbose=l_verbose)
    printlog ("outimage    = "//l_outimage, l_logfile, verbose=l_verbose)
    printlog ("secfile     = "//l_secfile, l_logfile, verbose=l_verbose)
    printlog ("gratingdb   = "//l_gratingdb, l_logfile, verbose=l_verbose)
    printlog ("filterdb    = "//l_filterdb, l_logfile, verbose=l_verbose)
    printlog ("bpmfile     = "//l_bpmfile, l_logfile, verbose=l_verbose)
    printlog ("gradimage   = "//l_gradimage, l_logfile, verbose=l_verbose)
    printlog ("refimage    = "//l_refimage, l_logfile, verbose=l_verbose)
    printlog ("fl_update   = "//l_fl_update, l_logfile, l_verbose)
    printlog ("fl_oversize = "//l_fl_oversize, l_logfile, l_verbose)
    if (l_uxoffset!=INDEF)
        printlog ("xoffset     = "//str(l_uxoffset), l_logfile, 
            verbose=l_verbose)
    else
        printlog ("xoffset     = INDEF", l_logfile, verbose=l_verbose)
    if (l_uyoffset!=INDEF)
        printlog ("yoffset     = "//str(l_uyoffset), l_logfile, 
            verbose=l_verbose)
    else
        printlog ("yoffset     = INDEF", l_logfile, verbose=l_verbose)
    if (l_uw2!=INDEF)
        printlog ("w2          = "//str(l_uw2), l_logfile, verbose=l_verbose)
    else
        printlog ("w2          = INDEF", l_logfile, verbose=l_verbose)
    printlog (" ", l_logfile, verbose=l_verbose)

    #check that there are input files
    if (l_inimage == "" || l_inimage == " "){
        printlog ("ERROR - GSCUT: input files not specified", l_logfile, 
            verbose+)
        goto error
    }

    # Check input image
    gimverify (l_inimage)
    if (gimverify.status!=0) {
        printlog ("ERROR - GSCUT: input file doesn't exist or is not a MEF",
            l_logfile, verbose+)
        goto error
    }
    l_inimage = gimverify.outname//".fits"

    # check existence of list file
    if (substr(l_inimage,1,1) == "@") {
        printlog ("ERROR - GSCUT: lists are currently not supported",
            l_logfile, verbose+)
        goto error
    }

    # check output name
    if (l_outimage != "") {
        fparse (l_outimage,verbose-)
        outroot = fparse.root
        if (fparse.extension == "")
            l_outimage = l_outimage//".fits"

        # Check that the output file does not already exist. If so, exit.
        if (imaccess(l_outimage)) {
            printlog ("ERROR - GSCUT: Output file "//l_outimage//\
                " already exists.", l_logfile, verbose+)
            goto error
        }
    } else if (!l_fl_update && l_secfile=="") {
            printlog ("ERROR - GSCUT: No output file specified, "//\
                "fl_update=no and secfile=\"\". Exiting.", l_logfile, verbose+)
            goto error
    }

    if (l_secfile != "" && access(l_secfile)) {
        printlog ("ERROR - GSCUT: Section file "//l_secfile//" already \
            exists.", l_logfile, verbose+)
        goto error
    }

    # Check that input image has been GMOSAICed
    imgets (l_inimage//"[0]","GMOSAIC", >& "dev$null")
    if (imgets.value == "0" && l_outimage!="") {
        printlog ("ERROR - GSCUT: inimage must be GMOSAICed", l_logfile, 
            verbose+)
        goto error
    } else
        impaste = l_inimage

    # get header information
    # Which instrument?
    imgets (impaste//"[0]","INSTRUME", >& "dev$null")
    if (imgets.value == "0") {
        printlog ("ERROR - GSCUT: Instrument keyword not found.",
            l_logfile, verbose+)
        goto error
    }

    inst = 1 # default is GMOS-N
    if (imgets.value == "GMOS-S")
        inst = 2

    # Set Hamamatsu flag
    imgets (impaste//"[0]", "DETTYPE", >& "dev$null")
    if (imgets.value == "S10892" || imgets.value == "S10892-N") {# Hamamatsu
        ishamamatsu = yes
        iccd = 3
    } else if (imgets.value == "SDSU II e2v DD CCD42-90") {# e2vDD CCDs
        iccd = 2
    } else if (imgets.value == "SDSU II CCD") {# EEV CCDs
        iccd = 1
    } else {
        printlog ("ERROR - GSCUT: DETTYPE keyword not known or not found",
            l_logfile, verbose+)
        goto clean
    }

    # check reference image
    if (l_refimage != "") {
        gimverify(l_refimage)
        if (gimverify.status!=0) {
            printlog ("ERROR - GSCUT: reference doesn't exist or is not a MEF",
                l_logfile, verbose+)
            goto error
        }
        l_refimage = gimverify.outname//".fits"
    }

    # Check tha the MDF exists & contains the appropriate keywords
    if (l_refimage != "")
        tinfo (l_refimage//"[MDF]", ttout-, >& "dev$null")
    else
        tinfo (l_inimage//"[MDF]", ttout-, >& "dev$null")
    if (tinfo.tbltype!="fits") {
        printlog ("ERROR - GSCUT: MDF file does not exist", l_logfile, 
            verbose+)
        goto error
    }
    if (l_refimage != "") {
        tlcol (l_refimage//"[MDF]", nlist=4) | match ("SECX2", stop-, \
                print_file_n+, metach+) |scan (mdftest)
        if (mdftest == "") {
            printlog ("ERROR - GSCUT: reference "//l_refimage//" MDF does ",
                       l_logfile, verbose+)
            printlog ("               not have slit edge info.", l_logfile, \
                verbose+)
            goto error
        }
    }

    # check gradient image
    if (l_gradimage != "" && l_refimage == "") {
        gimverify (l_gradimage)
        if (gimverify.status!=0) {
            printlog ("ERROR - GSCUT: gradient image doesn't exist or is \
                not a MEF", l_logfile, verbose+)
            goto error
        } else {
            l_gradimage = gimverify.outname//".fits"
            imgets (l_gradimage//"[0]", "GMOSAIC", >& "dev$null")
            if (imgets.value == "0") {
                printlog ("ERROR - GSCUT: gradimage must be GMOSAICed",
                    l_logfile, verbose+)
                goto error
            }
        }
    }
    # Compute gradient images
    if (l_gradimage != "" && l_refimage == "") {
        gradient (l_gradimage//"["//l_sci_ext//",1]", gradpos, "90",
            boundary="nearest", constant=0.)
        imarith (gradpos, "*", -1., gradneg)
    }
    
    # Print warning if creating a new variance plane and check GAIN/READNOISE
    if (fl_vardq && !imaccess(impaste//"["//l_var_ext//",1]")) {
        printlog ("WARNING - GSCUT: Image "//impaste//" does not have \
            a variance plane.", l_logfile, l_verbose)
        printlog ("Creating new variance planes for each cut using VAR = \
            (READNOISE)^2 + SCI/GAIN.", l_logfile, l_verbose)
        printlog ("Please be advised that this is an approximation of the \
             variance.", l_logfile, l_verbose)
        if (l_key_ron == "") {
            printlog ("WARNING - GSCUT: Parameter key_ron not found. Using ", \
                l_logfile, l_verbose)
            printlog ("default keyword, 'READNOISE.'", l_logfile, l_verbose)
            l_key_ron = "READNOISE"
        }
        if (l_key_gain == "") {
            printlog ("WARNING - GSCUT: Parameter key_gain not found. Using", \
                l_logfile, verbose+)
            printlog ("default keyword, 'GAIN.'", l_logfile, l_verbose)
                l_key_gain = "GAIN"
        }
        imgets (impaste//"[0]", l_key_gain, >& "dev$null")
        if (imgets.value == "0") {
            printlog ("WARNING - GSCUT: Keyword \
                "//l_key_gain//" not found in \
                header using a gain", l_logfile, verbose+)
            printlog ("of "//l_gain//" electrons per \
                ADU.", l_logfile, verbose)
        } else {
            l_gain = real(imgets.value)
        }
        imgets (impaste//"[0]", l_key_ron, >& "dev$null")
        if (imgets.value == "") {
            printlog ("WARNING - GSCUT: Keyword \
                "//l_key_ron//"not found in \
                header using a read noise", l_logfile, verbose+)
            printlog ("of "//l_ron//" electrons.", l_logfile, verbose+)
        } else { 
            l_ron = real(imgets.value)
        }
    }
    
    # Find slits
    # PyRAF doesn't like hselect/translit, prefer keypar
    #hselect (impaste//"["//l_sci_ext//",1]", "CCDSUM", yes) | \
    #    translit ("STDIN", '"', "", delete+, collapse-) | scan (xbin, ybin)
    keypar(impaste//"["//l_sci_ext//",1]", "CCDSUM")
    print(keypar.value) | scan(xbin,ybin)
    xscale = pixscale[inst,iccd]*real(xbin)
    yscale = pixscale[inst,iccd]*real(ybin)

    hselect (impaste//"[0]", "GRATING,FILTER1,FILTER2,GRWLEN,GRTILT,GRORDER", yes) | \
        scan (grating, filter[1], filter[2], cwave, tilt, grorder)
    tilt = tilt*pi/180.
    hselect (impaste//"["//l_sci_ext//",1]", "i_naxis1,i_naxis2", yes) |\
        scan (nxpix, nypix)

    # according to Mark, the input from gmosaic is either full frame (with
    # ROIs pasted into it if applicable) or CCD2
    if (nxpix > 2048/xbin) fullframe = yes
    else fullframe = no

    # get chip gap locations, if needed to avoid finding grad peaks in them;
    if (l_gradimage != "" && fullframe) {
        if (access(l_bpmfile)) {
            if (ishamamatsu) {
                if (inst == 2) {  #GMOS-S
                    fields (l_bpmfile, "1,2", lines="5") | scan (gapl[1], gaph[1])
                    fields (l_bpmfile, "1,2", lines="6") | scan (gapl[2], gaph[2])
                } else {   #GMOS-N
                    fields (l_bpmfile, "1,2", lines="8") | scan (gapl[1], gaph[1])
                    fields (l_bpmfile, "1,2", lines="9") | scan (gapl[2], gaph[2])
                }
            } else {
                fields (l_bpmfile, "1,2", lines="2") | scan (gapl[1], gaph[1])
                fields (l_bpmfile, "1,2", lines="3") | scan (gapl[2], gaph[2])
            }
            if (debug)
                printlog ("Gaps "//gapl[1]//"-"//gaph[1]//", "//\
                    gapl[2]//"-"//gaph[2], l_logfile, verbose+)
        } else {
            printlog ("ERROR - GSCUT: can't open "//l_bpmfile, l_logfile, \
                verbose+)
            goto error
        }
    }

    # get grating information
    match (grating, l_gratingdb, stop-, print+, meta+) | \
        scan (gname, grule, gblaze, gR, gcoverage, gwave1, gwave2,
        wavoffset, l_yoff)

    if (l_uxoffset != INDEF)
        wavoffset = l_uxoffset

    if (l_uyoffset != INDEF)
        l_yoff = l_uyoffset

    l_yoff = l_yoff/ybin
    print (((cwave*grule)/1.e6)) | \
        interp ("gmos$data/gratingeq.dat","STDIN", int_mode="spline", 
        curve_gen-) | scan (greq, gtilt)
    
    ##M Needs checking
    #greq=(cwave*grule)/1.e6
    gtilt = gtilt * pi/180.
    a = sin(gtilt+0.872665) / sin(gtilt)
    gR = 206265. * greq/(0.5*81.0*sin(gtilt))
    nmppx = a*xscale*cwave*81.0*sin(gtilt)/(206265.*greq)
    wave1 = gwave1
    wave2 = gwave2
    
    if (grorder == 2) {
    	gR = gR*4.
    	nmppx = nmppx/2.
    }

    printlog ("Central wavelength "//cwave, l_logfile, l_verbose)
    printlog ("Grating: "//gname, l_logfile, l_verbose)
    printlog ("Resolution (0.5'' slit): "//gR, l_logfile, l_verbose)
    printlog ("Anamorphic factor: "//a, l_logfile, l_verbose)
    printlog ("Grating tilt (header): "//(tilt*180./pi), l_logfile, l_verbose)
    printlog ("Calculated tilt: "//(gtilt*180./pi), l_logfile, l_verbose)
    printlog ("nm/pix = "//nmppx, l_logfile, l_verbose)
    printlog ("Wavelength offset [nm] ="//wavoffset, l_logfile, l_verbose)
    printlog ("Y offset (binned pixels) = "//l_yoff, l_logfile, l_verbose)

    # get filter information
    fwave1 = 0.0     ; wmn1 = 0.0     ; wmn2 = 0.0
    fwave2 = 99999.0 ; wmx1 = 99999.0 ; wmx2 = 99999.0
    if (filter[1] != "" && substr(filter[1],1,4) != "open")
        match (filter[1], l_filterdb, stop-, print+, meta+) | \
            scan (fname, wmn1, wmx1, ffile)
    if (filter[2] != "" && substr(filter[2],1,4) != "open")
        match (filter[2], l_filterdb, stop-, print+, meta+) | \
            scan (fname, wmn2, wmx2, ffile)

    fwave1 = max(wmn1,wmn2)
    fwave2 = min(wmx1,wmx2)

    # determine whether filter or grating limits wavelength coverage
    wave1 = max(wave1,fwave1)
    wave2 = min(wave2,fwave2)

    # This sets the hard red limit according to detector type if user doesn't
    # supply an upper limit
    if (wave2 > detector_upper_spec_limit[iccd])
        wave2 = detector_upper_spec_limit[iccd]

    # reset wave2 if user defined
    if (l_uw2 != INDEF) {
        # Check supplied value is greater than wave1
        if (l_uw2 > wave1) {
            wave2 = l_uw2
        } else {
            printlog ("WARNING - GSCUT: User supplied upper wavength limit "//\
                "is less than the calculated\n"//\
                "                 lower wavelength limit.\n"//\
                "                 Reseting the upper wavelngth to the "//\
                "calulated value ("//wave2//" [nm]).", l_logfile, verbose+)
        }
    }

    # in pixels
    speclen = nint((wave2-wave1)/nmppx)
    printlog ("Filter1: "//filter[1], l_logfile, l_verbose)
    printlog ("Filter2: "//filter[2], l_logfile, l_verbose)

    # pixel value of central wavelength from left (red) end of spectrum
    pixcwave = speclen - (cwave-wave1)/nmppx

    if (debug) {
        printlog ("wave1: "//wave1//"[nm] wave2: "//wave2//" [nm]", \
            l_logfile, verbose+)
        printlog ("speclen: "//speclen//" [pix] pixcwave: "//pixcwave//\
            " [pix]", l_logfile, verbose+)
    }

    # Copy MDF to output image
    if (l_outimage != "") {
        # slice out the MDF
        if (l_refimage != "") {
            tcopy (l_refimage//"[MDF]", mdf//".fits", verbose-)
        } else { 
            tcopy (l_inimage//"[MDF]", mdf//".fits", verbose-)
            # Add section and select columns to MDF
            tcalc (mdf//".fits", "EXTVER", 0, datatype="int")
            tcalc (mdf//".fits", "REFPIX1", 1, datatype="int")
            tcalc (mdf//".fits", "SECX1", 1, datatype="int")
            tcalc (mdf//".fits", "SECX2", 1, datatype="int")
            tcalc (mdf//".fits", "SECY1", 1, datatype="int")
            tcalc (mdf//".fits", "SECY2", 1, datatype="int")
            tcalc (mdf//".fits", "SPECCEN", 0., colfmt="f6.3")
        }
        wmef (mdf//".fits", l_outimage, extnames="MDF", verbose-, 
            phu=l_inimage, >& "dev$null")
        tcalc (l_outimage//"[MDF]", "SELECT", 1, datatype="int")
    } else if (l_fl_update) {
        tcalc (l_inimage//"[MDF]", "SELECT", 1, datatype="int")
        tcalc (l_inimage//"[MDF]", "EXTVER", 0, datatype="int")
        tcalc (l_inimage//"[MDF]", "REFPIX1", 1, datatype="int")
        tcalc (l_inimage//"[MDF]", "SECX1", 1, datatype="int")
        tcalc (l_inimage//"[MDF]", "SECX2", 1, datatype="int")
        tcalc (l_inimage//"[MDF]", "SECY1", 1, datatype="int")
        tcalc (l_inimage//"[MDF]", "SECY2", 1, datatype="int")
        tcalc (l_inimage//"[MDF]", "SPECCEN", 0., colfmt="f6.3")
    }

    # How many slits to extract? If LONGSLIT nr=1.
    hselect (l_inimage//"[0]", "OBSMODE", yes) | scan (modobs)

    if (modobs == "LONGSLIT") {
        nr = 1
    } else if (modobs == "MOS") {
        tinfo (l_inimage//"[MDF]", ttout-)
        nr = tinfo.nrows
    }
    
    # Loop over slits
    slittype = "" ; slitlen = 0.0
    xcen = real(nxpix)/2.
    ycen = real(nypix)/2.

    if (debug) {
        printlog ("xcen: "//xcen//" ycen: "//ycen, l_logfile, verbose+)
        printlog ("Number of rows in MDF = "//nr, l_logfile, verbose+)
    }

    n_j = 0
    nsciext = 0
    for (n_i=1; n_i<=nr; n_i+=1) {
        # Create tmp FITS file names used within this loop
        sci = mktemp("tmpsci")
        var = mktemp("tmpvar")
        dq = mktemp("tmpdq")
        tmppos = mktemp("tmppos")
        tmpneg = mktemp("tmpneg")
        tmpcoo = mktemp("tmpcoo")
        tmptmatch = mktemp("tmptmatch")

        if (debug) {
            printlog ("++++ Reading MDF row "//n_i//" ++++", \
                l_logfile, verbose+)
        }
        
        if (modobs == "MOS") {
            if (l_refimage != "") {
                tprint (l_refimage//"[MDF]", prparam-, prdata+, showrow-,
                    showhdr-, showunits-,
                    col="SECX1,SECX2,SECY1,SECY2,EXTVER,REFPIX1,priority,\
                    SPECCEN", rows=n_i,pwidth=160) | \
                    scan (x1, x2, y1, y2, n_j, refpix, priority, speccen)
                printlog ("Read parameters for slit "//n_i,
                    l_logfile, l_verbose)
            } else {
                # Make sure we have a priority for NS longslit
                priority = 1

                # Get parameters for slit from MDF (for MOS only)
                tprint (l_inimage//"[MDF]", prparam-, prdata+, pwidth=160,
                    plength=0, showrow-, orig_row+, showhdr-, showunits-, 
                    col="slitid,slitpos_mx,slitpos_my,slitsize_mx,slitsize_my,\
                    slitsize_mr,slitsize_mw,slittilt_m,SLITTYPE,priority", 
                    rows=n_i, option="plain", align=yes, sp_col="", lgroup=0) \
                    | scan (objid, spos_mx, spos_my, ssize_mx, ssize_my, \
                    ssize_mr, ssize_mw, stilt_m, slittype, priority)

                if (debug) {
                    tprint (l_inimage//"[MDF]", prparam-, prdata+, pwidth=160,
                        plength=0, showrow-, orig_row+, showhdr-, showunits-, 
                        col="slitid,x_ccd,slitpos_mx,slitpos_my,slitsize_mx,\
                            slitsize_my,slitsize_mr,slitsize_mw,slittilt_m,\
                            SLITTYPE,priority", \
                            rows=n_i, option="plain", align=yes, sp_col="", \
                            lgroup=0, >> l_logfile) 
                    tprint (l_inimage//"[MDF]", prparam-, prdata+, pwidth=160,
                        plength=0, showrow-, orig_row+, showhdr-, showunits-, 
                        col="slitid,x_ccd,slitpos_mx,slitpos_my,slitsize_mx,\
                            slitsize_my,slitsize_mr,slitsize_mw,slittilt_m,\
                            SLITTYPE,priority", \
                            rows=n_i, option="plain", align=yes, sp_col="", \
                            lgroup=0) 
                }

                if (priority != 0) {
                    # Convert from mask to pixel coordinates and correct for 
                    # binning in both directions
                    xccd = spos_mx * asecmm/xscale

                    # distortion correction in y
                    if (inst==2) {
                        #yccd=spos_my*asecmm/yscale
                        # Not only there is a y-offset (85) but there is also 
                        # a distortion.  The solution below is not perfect but 
                        # it is already much better than the first order 
                        # solution [Kathleen Labrie]
                        yccd = 0.99911*spos_my - 1.7465E-5*spos_my**2 + \
                            3.0494E-7*spos_my**3
                        yccd = yccd * asecmm/yscale
                    } else {
                        ##M Needs checking - if need put check of iccd in.
                        yccd = 0.99591859227*spos_my + \
                            5.3042211333437E-8*spos_my**2 + \
                            1.7447902551997E-7*spos_my**3
                        yccd = yccd * asecmm/yscale
                    }

                    if (debug) {
                        printlog ("First xccd: "//xccd//" yccd: "//yccd//\
                            " relative to array centre", \
                            l_logfile, verbose+) 
                    }

                    slitwid = ssize_mx*asecmm
                    slitlen = ssize_my*asecmm

                    # set slit length if the aperture is a circle, 
                    # radius=slitwid
                    if (slittype=="circle")
                        slitlen=2.*slitwid

                    xccd = xcen+xccd
                    yccd = ycen+yccd

                    # simple correction for distortion in x
                    y = (yccd/real(nypix)-0.5)
                    dx = real(nxpix) * (0.0014*y - 0.0167*y**2)
                    #print(yccd," ",y," ",dx)

                    if (debug) {
                        printlog ("Second xccd: "//xccd//" yccd: "//yccd//\
                            " relative to 1,1 coordiinate", \
                             l_logfile, verbose+) 
                        printlog ("slitwid: "//slitwid//" slitlen: "//slitlen,\
                            l_logfile, verbose+) 
                        printlog ("y: "//y//" dx: "//dx, \
                            l_logfile, verbose+) 
                    }

                    # slit height
                    # The old 1.05x for distortions is no longer needed when
                    # using gradient images but is kept as the default for
                    # backwards compatibility and to limit RTF breakage:
                    if (l_fl_oversize) {
                        specwid = nint(1.05*slitlen/yscale)
                    } else {
                        specwid = nint(slitlen/yscale)  # arg. already real
                    }
                    if (gradimage == "") {
                        printlog ("Slit "//n_i//": slit height = "//specwid,
                            l_logfile, l_verbose)
                    } else {
                        printlog ("Slit "//n_i//":", \
                            l_logfile, l_verbose)
                    }

                    center = specwid/2

                    if (debug) {
                        printlog ("specwid (height): "//specwid//" centre: "//\
                            center, \
                            l_logfile, verbose+) 
                    }
                }
            }
        } else if (modobs == "LONGSLIT") {
            dx = 0. # the distortion correction, IJ
            specwid = real(nypix)
            center = specwid/2.
            xccd = xcen
            yccd = center
            priority = 1
        }
        if (priority != 0) {
            if (l_refimage == "") {
                n_j = n_j+1
                refpix = pixcwave

                # Position of object, take into account that lambda decreases 
                # with x
                x1 = nint(xcen-(xcen-xccd)/a-pixcwave) + wavoffset/nmppx + dx
                #+0.000189*(yccd-2304./ybin) - 1.044231E-5*(yccd-2304./ybin)**2
                x2 = x1 + speclen-1
                y1 = nint(yccd-center+l_yoff)
                y2 = y1 + specwid-1

                if (debug) {
                    printlog ("Pre-checks x1: "//x1//" x2: "//x2//" y1: "//\
                        y1//" y2: "//y2, l_logfile, verbose+)
                    printlog ("    REFPIX1 "//refpix, l_logfile, l_verbose)
                }

                # check spectrum isn't off chip
                if (x1 < 1) {
                    refpix = refpix+real(x1-1)
                    x1 = 1
                }

                if (x2 > nxpix)
                    x2 = nxpix
                if (y1 < 1)
                    y1 = 1
                if (y2 > nypix)
                    y2 = nypix

                if (debug) {
                    printlog ("After-checks x1: "//x1//" x2: "//x2//\
                        " y1: "//y1//" y2: "//y2, l_logfile, verbose+)
                    printlog ("    REFPIX1 "//refpix, l_logfile, l_verbose)
                }

                if (modobs == "MOS" && l_gradimage != "") {

                    # This is a real number to be sure of maintaining the
                    # behaviour prior to 1.41, though int makes more sense:
                    apcol = x1+refpix-1
                    if (apcol < 1.0)
                        apcol = 1.0
                    else if (apcol > nxpix)
                        apcol = nxpix

                    if (debug)
                        printlog ("Initial apcol "//apcol, l_logfile, \
                            l_verbose)

                    # Columns to sum over for apfind:
                    apnsum = 20
                    aphsum = apnsum/2

                    # If the apfind column happens to fall in a chip gap,
                    # adjust it back onto the nearest CCD (adapted from IJ,
                    # 2014jun24 version):
                    if (fullframe) {

                        # Which gap (if any) is it in?
                        for (ngap=2; ngap >= 0; ngap-=1) {
                            if (ngap==0) break
                            if ((apcol+aphsum)*xbin >= gapl[ngap] && \
                                (apcol-aphsum)*xbin <= gaph[ngap]) break
                        }
                        if (debug)
                            printlog ("Falls in gap "//ngap, l_logfile, \
                                l_verbose)

                        # Move to whichever side of the gap is nearer:
			if (ngap > 0) {
                            if (gaph[ngap]-apcol*xbin < apcol*xbin-gapl[ngap])
                                apcol = (gaph[ngap]+aphsum*xbin+1)/xbin
                            else
                                apcol = (gapl[ngap]-aphsum*xbin)/xbin
                        }
                    }

                    if (debug) {
                        printlog ("MOS mode and gradimage = "//l_gradimage, \
                            l_logfile, verbose+)
                        printlog ("line number for apfind: "//apcol,\
                            l_logfile, verbose+)
                    }

                    # Find slit parameters from gradient image
                    # Find lower edges
                    apfind (gradpos, nfind=nr, line=apcol, nsum=apnsum,
                        inter-, minsep=3., > "dev$null")
                    match ("center", "database/ap"//gradpos, stop-, > tmppos)
                    # Find upper edges
                    apfind (gradneg, nfind=nr, line=apcol, nsum=apnsum,
                        inter-, minsep=3., > "dev$null")
                    match ("center","database/ap"//gradneg, stop-, > tmpneg)


                    print (y1, > tmpcoo)

                    if (debug) {
                        printlog ("y1 again: "//y1//" max(0.75*specwid,5.)"//\
                        " = "//max(0.75*specwid,5.), l_logfile, l_verbose)
                    }

                    tmatch (tmppos, tmpcoo, tmptmatch, "c3", "c1",
                        max(0.75*specwid,5.), >& "dev$null")
                    tinfo (tmptmatch, ttout-)

                    if (debug) {
                        printlog ("Passed lower edge tmatch and tinfo", \
                            l_logfile, l_verbose)
                    }

                    found = no
                    foundup = no
                    if(tinfo.nrows!=0) {
                        tcalc (tmptmatch, "c4", "abs(c3-c1_2)", colfmt="f6.3")
                        tcalc (tmptmatch, "c5", "c3-c1_2", colfmt="f6.3")
                        tsort (tmptmatch, "c4")
                        tprint (tmptmatch, col="c5", row="1", showr-, showh-,
                            >> tmpoff)
                        tabpar (tmptmatch,"c3",1)
                        y1 = nint(real(tabpar.value))
                        ry1 = real(tabpar.value)
                        found = yes
                    } else {
                        ry1=real(y1)
                    }
                    printlog ("MDF row="//n_i//" lower edge of slit \
                        found ="//found, l_logfile, l_verbose)
                    delete (tmpcoo//","//tmptmatch, verify-)
                    print (y2, > tmpcoo)

                    tmatch (tmpneg, tmpcoo, tmptmatch, "c3", "c1",
                        max(0.75*specwid,5.), >& "dev$null")
                    tinfo (tmptmatch, ttout-)

                    if (debug) {
                        printlog ("Passed upper edge tmatch and tinfo", \
                            l_logfile, l_verbose)
                    }

                    if (tinfo.nrows == 0) {
                        if (found) {
                            y2 = y1 + specwid-1
                            ry2 = ry1 + specwid-1
                        } else {
                            ry2=real(y2)
                        }
                    } else {
                        tcalc (tmptmatch, "c4", "abs(c3-c1_2)", colfmt="f6.3")
                        tcalc (tmptmatch, "c5", "c3-c1_2", colfmt="f6.3")
                        tsort (tmptmatch, "c4")
                        tprint (tmptmatch, col="c5", row="1", showr-, showh-,
                            >> tmpoff)
                        tabpar (tmptmatch, "c3", 1)
                        y2 = nint(real(tabpar.value))
                        ry2 = real(tabpar.value)
                        foundup = yes
                        if (!found) {
                            y1 = y2 - specwid+1
                            ry1 = ry2 - specwid+1
                        }
                    }

                    if (debug) {
                        printlog ("Calculated (int)  y1: "//y1//" and y2: "//\
                            y2, l_logfile, verbose+)
                        printlog ("Calculated (real) y1: "//ry1//\
                            " and y2: "//ry2, l_logfile, verbose+)
                    }

                    printlog ("MDF row="//n_i//" upper edge of slit \
                        found ="//foundup, l_logfile, l_verbose)

                    if ( abs(y2-y1-specwid+1)>3. ) {
                        printlog ("WARNING - GSCUT: Slitlength from gradimage \
                            more than 3 pixels different from MDF info",
                            l_logfile, verbose+)
                        printlog ("                 Using lower edge + MDF \
                            info to get upper edge", l_logfile, verbose+)
                        y2 = y1 + specwid-1
                        ry2 = ry1 + specwid-1
                    }
                    speccen = (ry2 + ry1) / 2.
                    if (y1 < 1)
                        y1 = 1
                    if (y2 > nypix)
                        y2 = nypix

                    printlog ("slit height = "//(y2-(y1-1)),
                        l_logfile, l_verbose)

                    if (debug) {
                        printlog ("Re-calcualted (int)  y1: "//y1//\
                            " and y2: "//y2, \
                            l_logfile, verbose+)
                        printlog ("Re-calcualted (real) y1: "//ry1//\
                            " and y2: "//ry2, \
                            l_logfile, verbose+)
                    }

                    # clean
                    delete ("database/ap"//gradpos, verify-, >& "dev$null")
                    delete ("database/ap"//gradneg, verify-, >& "dev$null")
                    delete (tmpcoo//","//tmppos//","//tmpneg//","//tmptmatch,
                        verify-, >& "dev$null")
                }
            }
            
            # Add the extra y here if needed.  Ensures it is added in all cases
            y1 = nint(y1 - l_yadd/ybin)
            y2 = nint(y2 + l_yadd/ybin)

            if (l_yadd != 0) {
                printlog ("Added "//nint(l_yadd/ybin)//\
                    " binned pixels to both "//\
                    "ends of slit height.", l_logfile, l_verbose)
                printlog ("Slight height is now: "//(y2-(y1-1)), \
                    l_logfile, l_verbose)
            }

            if (debug) {
                printlog ("l_y add: "//l_yadd//" ybin: "//ybin, \
                    l_logfile, verbose+)
                printlog ("Re-re-calcualted (int)  y1: "//y1//\
                    " and y2: "//y2, \
                    l_logfile, verbose+)
            }

            sec = "["//x1//":"//x2//","//y1//":"//y2//"]"
            datsec = "[1:"//(x2-x1+1)//",1:"//(y2-y1+1)//"]"
            # detsec="["//x1//":"//x2//","//y1//":"//y2//"]"
            # Correct detsec for binning, such that gdisplay will work on 
            # output
            detsec = "["//str((x1-1)*xbin+1)//":"//str(x2*xbin)//","//\
                str((y1-1)*ybin+1)//":"//str(y2*ybin)//"]"

            if (debug) {
                printlog ("Section of input image to be cut:\n    "//sec,\
                    l_logfile, verbose+)
                printlog ("Output data section:\n    "//datsec,\
                    l_logfile, verbose+)
                printlog ("Detector section of cut slit:\n    "//detsec,\
                    l_logfile, verbose+)
            }

            if (l_secfile != "") 
                print (sec, >> l_secfile)

            # Put image section info in MDF columns
            if (l_outimage == "" && l_fl_update) {
                partab (x1, l_inimage//"[MDF]", "SECX1", n_i)
                partab (x2, l_inimage//"[MDF]", "SECX2", n_i)
                partab (y1, l_inimage//"[MDF]", "SECY1", n_i)
                partab (y2, l_inimage//"[MDF]", "SECY2", n_i)
                partab (speccen, l_inimage//"[MDF]", "SPECCEN", n_i)
                partab (n_j, l_inimage//"[MDF]", "EXTVER", n_i)
                partab (refpix, l_inimage//"[MDF]", "REFPIX1", n_i)
            }

            if (l_outimage != "") {
                # Put image section info in MDF columns
                if (modobs == "MOS") {
                    partab (x1, l_outimage//"[MDF]", "SECX1", n_i)
                    partab (x2, l_outimage//"[MDF]", "SECX2", n_i)
                    partab (y1, l_outimage//"[MDF]", "SECY1", n_i)
                    partab (y2, l_outimage//"[MDF]", "SECY2" ,n_i)
                    partab (speccen, l_outimage//"[MDF]", "SPECCEN", n_i)
                    partab (n_j, l_outimage//"[MDF]", "EXTVER", n_i)
                    partab (refpix, l_outimage//"[MDF]", "REFPIX1", n_i)
                }

                # Cut this slit and insert in output file
                imcopy (impaste//"["//l_sci_ext//",1]"//sec, sci, verbose-)

                # Delete LTV1/2 keywords as, LTV1 messes up gswavelength when 
                # SECX1 is not 1. LTV1 maps the offset of the cut position
                # to the origin (in physical units) of the image from which
                # it was cut. - MS
                # This is not done here as the WCS etc., still refers to the
                # original image and therefore LTV1/2 should only be updated
                # when we update the WCS in GSAPPWAVE - MS

                hselect (l_outimage//"[0]", "NEXTEND", yes) | scan (nextnd)
                fxinsert (sci//".fits", l_outimage//"["//nextnd//"]","0",
                    verbose-, >& "dev$null")
                nextnd = nextnd+1
                gemhedit (l_outimage//"[0]", "NEXTEND", nextnd, "", delete-)
                gemhedit (l_outimage//"["//nextnd//"]", "EXTNAME", l_sci_ext,
                    "", delete-)
                gemhedit (l_outimage//"["//nextnd//"]", "EXTVER", n_j, "", \
                    delete-)
                printlog ("Inserted slit as extension [SCI,"//str(n_j)//"]",
                    l_logfile, l_verbose)
                gemhedit (l_outimage//"["//nextnd//"]", "MDFROW", n_i,
                    "Corresponding row in MDF", delete-)
                gemhedit (l_outimage//"["//nextnd//"]", "DETSEC", detsec, \
                    "Detector section(s)", delete-)
                gemhedit (l_outimage//"["//nextnd//"]", "DATASEC", datsec, \
                    "Data section(s)", delete-)
                gemhedit (l_outimage//"["//nextnd//"]", "CCDSUM", \
                    xbin//" "//ybin, "", delete-)

                #update header with wavelength information
                gemhedit (l_outimage//"["//nextnd//"]", "REFPIX1", refpix,
                    "Pixel of central wavelength", delete-)
                # Variance and DQ planes
                if (l_fl_vardq) {
                    if (imaccess(impaste//"["//l_var_ext//",1]")) {
                        imcopy (impaste//"["//l_var_ext//",1]"//sec,
                            var, verbose-)
                    } else {

                        # Do the math, create variance plane using
                        # gain, readnoise, and sci image
                        
                        imexpr ("a**2 + b/c", var, l_ron, sci, l_gain, \
                            dims="auto", intype="auto", outtype="real", \
                            refim="auto", bwidth=0, btype="nearest", \
                            bpixval=0., rangecheck=yes, verbose=no,
                            exprdb="none", lastout="")
                    }

                    if (imaccess(impaste//"["//l_dq_ext//",1]")) {
                        imcopy (impaste//"["//l_dq_ext//",1]"//sec,
                            dq, verbose-)
                    } else {
                        imarith (sci, "*", 0, dq, pixtype="int",
                            calctype="int", verbose-)
                    }
                    # Insert
                    fxinsert (var//".fits,"//dq//".fits",
                        l_outimage//"["//nextnd//"]", "0", verbose-,
                        >& "dev$null")
                    gemhedit (l_outimage//"[0]", "NEXTEND", (nextnd+2), "",
                        delete-)
                    # Variance header
                    gemhedit (l_outimage//"["//(nextnd+1)//"]", "EXTNAME",
                        l_var_ext, "", delete-)
                    gemhedit (l_outimage//"["//(nextnd+1)//"]", "EXTVER", n_j,
                        "", delete-)
                    gemhedit (l_outimage//"["//(nextnd+1)//"]", "MDFROW", n_i,
                        "Corresponding row in MDF", delete-)
                    gemhedit (l_outimage//"["//(nextnd+1)//"]", "REFPIX1",
                        refpix, "Reference pixel of central wavelength",
                        delete-)
                    gemhedit (l_outimage//"["//(nextnd+1)//"]", "DETSEC", 
                        detsec, "Detector section(s)", delete-)
                    gemhedit (l_outimage//"["//(nextnd+1)//"]", "DATASEC", 
                        datsec, "Data section(s)", delete-)
                    gemhedit (l_outimage//"["//(nextnd+1)//"]", "CCDSUM",
                        xbin//" "//ybin, "", delete-)

                    # DQ header
                    gemhedit (l_outimage//"["//(nextnd+2)//"]", "EXTNAME",
                        l_dq_ext, "", delete-)
                    gemhedit (l_outimage//"["//(nextnd+2)//"]", "EXTVER", n_j,
                        "", delete-)
                    gemhedit (l_outimage//"["//(nextnd+2)//"]", "MDFROW", n_i,
                        "Corresponding row in MDF", delete-)
                    gemhedit (l_outimage//"["//(nextnd+2)//"]", "REFPIX1",
                        refpix, "Reference pixel of central wavelength",
                        delete-)
                    gemhedit (l_outimage//"["//(nextnd+2)//"]", "DETSEC", 
                        detsec, "Detector section(s)", delete-)
                    gemhedit (l_outimage//"["//(nextnd+2)//"]", "DATASEC", 
                        datsec, "Data section(s)", delete-)
                    gemhedit (l_outimage//"["//(nextnd+2)//"]", "CCDSUM",
                        xbin//" "//ybin, "", delete-)
                }
            } #end outimage

            # In cases where there is an entry at the end of the MDF with a
            # priority equal to 0, the n_j variable is overwritten with a value
            # of 0. Since the number of science extensions (== n_j) keyword is
            # written outside of this loop, store this value in a different
            # variable.
            nsciext = n_j

        } # should be end priority
        imdelete (sci//","//var//","//dq, verify-, >& "dev$null")

    } #end loop over slits

    # Median yoffset
    if (access(tmpoff)) {
        tstat (tmpoff, "c1", outtable="STDOUT", rows="-", lowlim=INDEF,
            highlim=INDEF, >& "dev$null")
        printlog ("Median Y offset for automatically found slits "//\
            str(tstat.median+l_yoff), l_logfile, l_verbose)
        delete (tmpoff, verify-, >& "dev$null")
    }

    # final header update
    if (l_outimage != "") {
        gemhedit (l_outimage//"[0]", "NSCIEXT", nsciext, "", delete-)
            
        gemdate ()
        gemhedit (l_outimage//"[0]", "GSCUT", gemdate.outdate, 
            "UT Time stamp for GSCUT", delete-)
        gemhedit (l_outimage//"[0]", "GEM-TLM", gemdate.outdate, 
            "UT Last edit with GEMINI", delete-)
    } else if (l_fl_update) {
        gemdate ()
        
        #  TODO
        #  The image has not been cut, only the MDF has been updated.
        #  Now, I believe there should still be a GSCUT time stamp for
        #  that.  However, I believe that the gmos package uses the GSCUT
        #  keyword to check if the image has actually been cut.
        #  I would suggest using another keyword, i.e. CUT (??) to indicate
        #  that the data have been cut and keep GSCUT as the time stamp.
        
        #gemhedit (l_inimage//"[0]", "GSCUT", gemdate.outdate, 
        #    "UT Time stamp for GSCUT", delete-)
        gemhedit (l_inimage//"[0]", "GEM-TLM", gemdate.outdate, 
            "UT Last edit with GEMINI", delete-)
    }

    # clean up
    goto clean

error:
    status = 1
    goto clean

clean:
    specred.dispaxis = s_dispaxis ; specred.logfile = s_logfile
    specred.database = s_database
    
    delete ("tmpslit*.fits,"//mdf//".fits", verify-, >& "dev$null")
    imdelete (sci//","//var//","//dq, verify-, >& "dev$null")
    imdelete (gradpos//","//gradneg, verify-, >& "dev$null")
    # close log file
    printlog (" ", l_logfile, l_verbose)
    if (status == 0)
        printlog ("GSCUT done", l_logfile, l_verbose)
    else
        printlog ("GSCUT failed", l_logfile, l_verbose)
    printlog ("---------------------------------------------------------\
        -----------------------", l_logfile, l_verbose)

end
