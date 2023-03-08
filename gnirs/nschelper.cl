# Copyright(c) 2001-2011 Association of Universities for Research in Astronomy, Inc.

# combined combiner code - called by nscombine and nsstack

procedure nschelper (inimages) 

char    inimages    {prompt = "Input images or spectra to shift/group and combine"}			# OLDP-1-primary-combine-suffix=_add
real    tolerance   {0.5, min = 0., prompt = "Maximum offset for grouping the frames (arcsec)"}	# OLDP-2
char    outimages   {"", prompt = "File name(s) for output shifted/grouped and combined spectra"}	# OLDP-1-output
char    output_suffix   {"_comb", prompt = "Suffix for output shifted and combined image"}		# OLDP-3
bool    fl_shift    {yes, prompt = "Shift?  Otherwise, group."}					# OLDP-1
char    bpm         {"", prompt = "Name of bad pixel mask."}					# OLDP-1
int     dispaxis    {1, prompt = "Dispersion axis (if not in header)"}				# OLDP-3
real    pixscale    {1., prompt = "Pixel scale in arcsec/pixel (if not in header)"}			# OLDP-3
bool    fl_cross    {no, prompt = "Update WCS from cross-correlation?"}				# OLDP-3
bool    fl_keepshift    {no, prompt = "Keep shifted images?"}             				# OLDP-3
bool    fl_shiftint {yes, prompt = "Shift frames by integer pixels?"}				# OLDP-2
char    interptype  {"linear", prompt = "Interpolation type for shifting",  enum = "nearest|linear|poly3|poly5|spline3|sinc|drizzle"}						# OLDP-2
char    boundary    {"nearest", prompt = "Boundary type for shifting", enum = "constant|nearest|reflect|wrap"}									# OLDP-2
real    constant    {0., prompt = "Constant value for boundary extension when shifting"}		# OLDP-2
char    combtype    {"average", prompt = "Combination operation", enum = "average|median"}		# OLDP-2
char    rejtype     {"sigclip", prompt = "Rejection algorithm when combining", enum = "none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip"}					# OLDP-2
char    masktype    {"goodvalue", enum="none|goodvalue", prompt="Mask type when combining"} # OLDP-3
real    maskvalue   {0., prompt="Mask value when combining"}                # OLDP-3
char    statsec     {"[*,*]", prompt = "Image section to be used for statistics"}			# OLDP-2
char    scale       {"none", prompt = "Image scaling"}						# OLDP-2
char    zero        {"none", prompt = "Image zeropoint offset"}					# OLDP-2
char    weight      {"none", prompt = "Image weights"}						# OLDP-2
real    lthreshold  {INDEF, prompt = "Lower threshold"}						# OLDP-2
real    hthreshold  {INDEF, prompt = "Upper threshold"}						# OLDP-2
int     nlow        {1, min = 0, prompt = "minmax: Number of low pixels to reject"}			# OLDP-2
int     nhigh       {1, min = 0, prompt = "minmax: Number of high pixels to reject"}		# OLDP-2
int     nkeep       {0, prompt = "Minimum to keep or maximum to reject"}  				# OLDP-2
bool    mclip       {yes, prompt = "Use median in sigma clipping algorithms?"}			# OLDP-2
real    lsigma      {5., prompt = "Lower sigma clipping factor"}					# OLDP-2
real    hsigma      {5., prompt = "Upper sigma clipping factor"}					# OLDP-2
real    ron         {0.0, min = 0., prompt = "Readout noise rms in electrons"}			# OLDP-2
real    gain        {1.0, min = 0.00001, prompt = "Gain in e-/ADU"}					# OLDP-2
char    snoise      {"0.0", prompt = "ccdclip: Sensitivity noise (electrons)"}			# OLDP-2
real    sigscale    {0.1, prompt = "Tolerance for sigma clipping scaling correction"}		# OLDP-2
real    pclip       {-0.5, prompt = "pclip: Percentile clipping parameter"}				# OLDP-2
real    grow        {0.0, prompt = "Radius (pixels) for neighbor rejection"}			# OLDP-2
char    nrejfile    {"", prompt = "Name of rejected pixel count image."}				# OLDP-3
bool    fl_vardq    {yes, prompt = "Create output variance and data quality frames?"}		# OLDP-2
bool    fl_inter    {no, prompt = "Measure cross-correlation peak interactively?"}			# OLDP-4
char    logfile     {"", prompt = "Logfile"}							# OLDP-1
char    logname     {"NSCHELPER", prompt = "Name printed in messages"}				# OLDP-4
bool    verbose     {yes, prompt = "Verbose output?"}						# OLDP-2
bool    debug       {no, prompt = "Very verbose output?"}						# OLDP-2
bool    force       {no, prompt = "Force use with earlier IRAF versions?"}				# OLDP-3
int     status      {0, prompt = "Exit status (0=good)"}						# OLDP-4
struct  *scanfile1  {"", prompt = "Internal use only"} 						# OLDP-4
struct  *scanfile2  {"", prompt = "Internal use only"} 						# OLDP-4
struct  *scanfile3  {"", prompt = "Internal use only"} 						# OLDP-4
struct  *scanfile4  {"", prompt = "Internal use only"} 						# OLDP-4

begin
        
    char    l_inimages = ""
    real    l_tolerance
    char    l_output = ""
    bool    l_fl_shift
    char    l_output_suffix = ""
    bool    l_shift
    char    l_bpm = ""
    int     l_dispaxis
    real    l_pixscale
    bool    l_fl_cross
    bool    l_fl_keepshift
    bool    l_fl_shiftint
    char    l_interptype = ""
    char    l_boundary = ""
    real    l_constant
    char    l_combtype = ""
    char    l_rejtype = ""
    char    l_masktype = ""
    real    l_maskvalue
    char    l_statsec = ""
    char    l_scale = ""
    char    l_zero = ""
    char    l_weight = ""
    real    l_lthreshold
    real    l_hthreshold
    int     l_nlow
    int     l_nhigh
    int     l_nkeep
    bool    l_mclip
    real    l_lsigma
    real    l_hsigma
    real    l_ron
    real    l_gain
    char    l_snoise = ""
    real    l_sigscale
    real    l_pclip
    real    l_grow
    char    l_nrejfile = ""
    bool    l_fl_vardq
    bool    l_fl_inter
    char    l_logfile = ""
    char    l_logname = ""
    bool    l_verbose
    bool    l_debug
    bool    l_force

    char    l_key_pixscale = ""
    char    l_key_xoff = ""
    char    l_key_yoff = ""
    char    l_key_exptime = ""
    char    l_key_ron = ""
    char    l_key_gain = ""
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_dispaxis = ""

    char    shft_suffix
    int     varext, dqext, index, nbad, i, junk, nsci
    char    badhdr, tmp, tmpout, tmpdel, phu, imgin, extname, target
    bool    intdbg, dynout, first, shfted
    char    imgout, imgnrej, imgname, keywdname
    struct  sdate, line
    char    tmpimages, tmpokimages, tmptarget, tmpoffset, tmpshiftcoo
    char    tmpshift, img, keyfound, tmpnrej, tmpshiftout2
    char    tmpshiftvar, tmpshiftvarout, tmpshiftdq, tmpshiftdqout
    char    dummy, tmphead, tmpcomb, outimg, nrejimg, tmpshifted
    char    combfile, tmpshiftout, tmpshiftin, imgfrom, imgto, tmpcopy
    char    tmpdel2, tmpcross
    int     nfiles, nn, ntarg, noff, nshift, disp_axis, ninput
    real    xoff, yoff, xrefoff, yrefoff, dispoff, pix_scale

    junk = fscan (  inimages, l_inimages)
    l_tolerance =   tolerance
    junk = fscan (  outimages, l_output)
    l_fl_shift =    fl_shift
    junk = fscan (  output_suffix, l_output_suffix)
    junk = fscan (  bpm, l_bpm)
    l_dispaxis =    dispaxis
    l_pixscale =    pixscale
    l_fl_cross =    fl_cross
    l_fl_keepshift =    fl_keepshift
    l_fl_shiftint = fl_shiftint
    junk = fscan (  interptype, l_interptype)
    junk = fscan (  boundary, l_boundary)
    l_constant =    constant
    junk = fscan (  combtype, l_combtype)
    junk = fscan (  rejtype, l_rejtype)
    junk = fscan (  masktype, l_masktype)
    l_maskvalue = maskvalue
    l_statsec =     statsec		# may contain spaces
    junk = fscan (  scale, l_scale)
    junk = fscan (  zero, l_zero)
    junk = fscan (  weight, l_weight)
    l_lthreshold =  lthreshold
    l_hthreshold =  hthreshold
    l_nlow =    nlow
    l_nhigh =   nhigh
    l_nkeep =   nkeep
    l_mclip =   mclip
    l_lsigma =  lsigma
    l_hsigma =  hsigma
    l_ron =     ron
    l_gain =    gain
    junk = fscan (  snoise,	l_snoise)
    l_sigscale =    sigscale
    l_pclip =       pclip
    l_grow =        grow
    junk = fscan (  nrejfile, l_nrejfile)
    l_fl_vardq =    fl_vardq
    l_fl_inter =    fl_inter
    junk = fscan (  logfile, l_logfile)
    junk = fscan (  logname, l_logname)
    l_verbose =     verbose
    l_debug =       debug
    l_force =       force

    cache ("keypar", "gemcombine", "nsoffset", "gnirs", "gemextn", "gemdate") 


    shft_suffix = "_shift" 

    status = 1
    intdbg = no

    badhdr = ""
    junk = fscan (  nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"
    junk = fscan (  nsheaders.key_xoff, l_key_xoff)
    if ("" == l_key_xoff) badhdr = badhdr + " key_xoff"
    junk = fscan (  nsheaders.key_yoff, l_key_yoff)
    if ("" == l_key_yoff) badhdr = badhdr + " key_yoff"
    junk = fscan (  nsheaders.key_exptime, l_key_exptime)
    if ("" == l_key_exptime) badhdr = badhdr + " key_exptime"
    junk = fscan (  nsheaders.key_ron, l_key_ron)
    if ("" == l_key_ron) badhdr = badhdr + " key_ron"
    junk = fscan (  nsheaders.key_gain, l_key_gain)
    if ("" == l_key_gain) badhdr = badhdr + " key_gain"
    junk = fscan (  nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (  nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (  nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (  nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"



    tmpimages = mktemp ("tmpimages") 
    tmpokimages = mktemp ("tmpokimages") 
    tmpoffset = mktemp ("tmpoffset") 
    tmpshiftcoo = mktemp ("tmpshiftcoo") 
    tmpshift = mktemp ("tmpshift") 
    tmpshiftvar = mktemp ("tmpshiftvar") 
    tmpshiftvarout = mktemp ("tmpshiftvarout") 
    tmpshiftdq = mktemp ("tmpshiftdq") 
    tmpshiftdqout = mktemp ("tmpshiftdqout") 
    tmphead = mktemp ("tmphead") 
    tmpnrej = mktemp ("tmpnrej")
    tmpout = mktemp ("tmpout")
    tmpcomb = mktemp ("tmpcomb")
    tmpdel = mktemp ("tmpdel")
    tmpdel2 = mktemp ("tmpdel2")
    tmpcross = mktemp ("tmpcross")
    tmptarget = mktemp ("tmptarget") 

    extname = l_sci_ext
    if (l_fl_vardq)
        extname = extname // "," // l_var_ext // "," // l_dq_ext


    if (l_logname == "") l_logname = "NSCHELPER"

    if (l_logfile == "") {
        junk = fscan (gnirs.logfile, l_logfile) 
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - "//l_logname//": Both nscombine.logfile \
                and gnirs.logfile are empty.", l_logfile, verbose+) 
            printlog ("                     Using default file " \
                // l_logfile, l_logfile, verbose+) 
        }
    }

    date | scan (sdate) 
    printlog ("---------------------------------------------------------\
        -----------------------", l_logfile, l_verbose) 
    printlog (l_logname//": " // sdate, l_logfile, verbose = l_verbose) 
    printlog (" ", l_logfile, verbose = l_verbose) 
    if (l_fl_shift) {
        printlog ("Shifting and combining to single file.", \
            l_logfile, l_verbose) 
    } else {
        printlog ("Grouping and combining by position.", \
            l_logfile, l_verbose) 
    }
    printlog ("Input images = " // l_inimages, l_logfile, l_verbose) 
    if (l_verbose)
        gemextn (l_inimages, check="", process="none", index="", extname="",
            extversion="", ikparams="", omit="extension,kernel", replace="",
            outfile="STDOUT", logfile="dev$null", glogpars="", verbose-)
    gemextn (l_inimages, check="", process="none", index="", \
        extname="", extversion="", ikparams="", omit="extension,kernel", \
        replace="", outfile="STDOUT", logfile="dev$null", glogpars="",
        verbose-, >> l_logfile) 
    printlog ("Bad pixel mask = " // l_bpm, l_logfile, l_verbose) 
    printlog ("Offset tolerance = " // l_tolerance, l_logfile, l_verbose) 
    printlog ("Statistics section = " // l_statsec, l_logfile, l_verbose) 
    printlog ("Combine algorithm = " // l_combtype, l_logfile, l_verbose) 
    printlog ("Rejection algorithm = " // l_rejtype, l_logfile, l_verbose) 
    if (fl_shift) {
        printlog ("Use integer pixel shifts? " // l_fl_shiftint, \
            l_logfile, l_verbose) 
        printlog ("Boundary extension  = " // l_boundary, \
            l_logfile, l_verbose) 
    }

    if ("" != badhdr) {
        printlog ("ERROR - NCHELPER: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }


    if (intdbg) print ("checking stats section")
    junk = fscan (l_statsec, tmp)	# (may contain spaces)
    if (tmp == "") {
        printlog ("WARNING - "//l_logname//": No statsec defined.  \
            Using entire array.", l_logfile, verbose+) 
        l_statsec = "[*,*]"
    }


    if (intdbg) print ("checking pixel scale")
    if (l_key_pixscale == "") {
        if (l_pixscale != 0) {
            printlog ("WARNING - "//l_logname//": key_pixscale not set, \
                using pixscale value = " // l_pixscale, \
                l_logfile, verbose+) 
        } else {
            printlog ("ERROR - "//l_logname//": No pixscale.", \
                l_logfile, verbose+) 
            goto clean
        }
    }


    if (intdbg) print ("checking input files")
    gemextn (l_inimages, check="exists,mef", process="none", index="", \
        extname="", extversion="", ikparams="", omit="extension,kernel", \
        replace="", outfile=tmpokimages, logfile="", glogpars="",
        verbose=l_verbose)
    if (0 < gemextn.fail_count) {
        printlog ("ERROR - "//l_logname//": Bad or missing file in \
            inimages.", l_logfile, verbose+) 
        goto clean
    }
    ninput = gemextn.count

    # If only one file to combine - exit.
    # TODO - why is this an error?  seems ok to me.
    if (ninput == 1) {
        printlog ("ERROR - "//l_logname//": Only one image to process.", \
            l_logfile, verbose+) 
        goto clean
    }


    if (intdbg) print ("checking for nsprepare etc")
    nbad = 0
    scanfile1 = tmpokimages
    while (fscan (scanfile1, imgin) != EOF) {
        phu = imgin // "[0]"
        keyfound=""
        hselect(phu, "*PREPAR*", yes) | scan(keyfound)
        if (keyfound == "") {
            printlog ("ERROR - "//l_logname//": Image " // imgin \
                // " not PREPAREd.", l_logfile, l_verbose) 
            nbad += 1
        }
    }

    if (nbad > 0) {
        printlog ("ERROR - "//l_logname//": " // nbad // " image(s) \
            have not been run through *PREPARE", l_logfile, verbose+) 
        goto clean
    }


    if (intdbg) print ("checking var/dq")
    if (l_fl_vardq) {
        gemextn (l_inimages, check="exists,mef", process="expand", \
            index="", extname=l_sci_ext, extversion="1-", ikparams="", \
            omit="extension,kernel", replace="", outfile="dev$null", \
            logfile="", glogpars="", verbose=l_verbose)
        nsci = gemextn.count
        gemextn (l_inimages, check="exists,mef", process="expand", \
            index="", extname=l_var_ext, extversion="1-", ikparams="", \
            omit="extension,kernel", replace="", outfile="dev$null", \
            logfile="", glogpars="", verbose=l_verbose)
        if (nsci != gemextn.count) {
            printlog ("WARNING - "//l_logname//": Some " // l_var_ext \
                // " extensions missing.", l_logfile, verbose+) 
            printlog ("                    Setting fl_vardq=no", \
                l_logfile, verbose+) 
            l_fl_vardq = no
        }
        gemextn (l_inimages, check="exists,mef", process="expand", \
            index="", extname=l_dq_ext, extversion="1-", ikparams="", \
            omit="extension,kernel", replace="", outfile="dev$null", \
            logfile="", glogpars="", verbose=l_verbose)
        if (l_fl_vardq && nsci != gemextn.count) {
            printlog ("WARNING - "//l_logname//": Some " // l_dq_ext \
                // " extensions missing.", l_logfile, verbose+) 
            printlog ("                    Setting fl_vardq=no", \
                l_logfile, verbose+) 
            l_fl_vardq = no
        }
    }


    if (intdbg) print ("checking output")

    # - either from a list (which we know must have a single entry 
    # if fl_shift+) or dynamically by adding postfix to input

    gemextn (l_output, check="", process="none", index="", \
        extname="", extversion="", ikparams="", omit="extension,kernel", \
        replace="", outfile=tmpout, logfile="", glogpars="", 
        verbose=l_verbose)
    if (gemextn.count > 1 && fl_shift) {
        printlog ("ERROR - "//l_logname//": Several output files, but \
            fl_shift+.", l_logfile, verbose+) 
        goto clean
    } else if (gemextn.count > 0) {
        dynout = no
        if (intdbg) type (tmpout)
    } else {
        dynout = yes
        touch (tmpout)
    }
    if (intdbg) print ("dynamic? " // dynout)

    if (no == dynout) {
        gemextn ("@" // tmpout, check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", \
            omit="extension,kernel", replace="", outfile="dev$null", \
            logfile="", glogpars="", verbose=l_verbose)
        if (0 != gemextn.fail_count || 1 != gemextn.count) {
            printlog ("ERROR - "//l_logname//": Output file exists.", \
                l_logfile, verbose+)
            goto clean
        }
    }


    if (intdbg) print ("checking nrejfile")
    if ("" != l_nrejfile) {
        if (no == l_fl_vardq) {
            printlog ("WARNING - "//l_logname//": nrejfile given, but \
                fl_vardq-, so no output.", l_logfile, verbose+)
            l_nrejfile = ""
        } else {
            gemextn (l_nrejfile, check="absent", process="none", \
                index="", extname="", extversion="", ikparams="", \
                omit="", replace="", outfile=tmpnrej, logfile="",
                glogpars="", verbose=l_verbose)
            if (gemextn.fail_count > 0 || 0 == gemextn.count) {
                printlog ("ERROR - "//l_logname//": Problem with \
                    nrejfile.", l_logfile, verbose+)
                goto clean
            }
        }
    }
    if ("" != l_nrejfile) touch (tmpnrej)



    # First processing step - group with gemoffsetlist and shift
    # if required        

    head (tmpokimages, nlines = 1) | scan (target) 

    first = yes
    while (ninput > 0) {

        head (tmpokimages, nlines = 1) | scan (img) 
        tmptarget = mktemp ("tmptarget") 

        if (intdbg) print ("grouping with gemoffsetlist")
        gemoffsetlist ("@" // tmpokimages, img, \
            distance = l_tolerance, \
            age = INDEF, fl_younger+, fl_nearer+, fl_noref-, \
            wcs_source = "phu", targetlist = tmptarget, \
            direction = 3, offsetlist = tmpoffset, key_xoff = l_key_xoff, \
            key_yoff = l_key_yoff, logfile = l_logfile, \
            verbose = l_debug, key_date="", key_time = "", force = l_force) 
        if (gemoffsetlist.status == 1) {
            printlog ("ERROR - "//l_logname//": gemoffsetlist could not \
                determine the offsets.", l_logfile, verbose+) 
            goto clean
        }

        ntarg = 0; noff = 0
        count (tmptarget) | scan (ntarg) 
        ninput = ninput - ntarg
        delete (tmpokimages, ver-, >& "dev$null")
        if (access (tmpoffset)) rename (tmpoffset, tmpokimages, field="all")


        # Sharpen WCS using cross-correlation if required

        if (l_fl_cross && ! first) {

            tmpcross = mktemp (tmpcross)
            touch (tmpcross)
            scanfile1 = tmptarget
            while (fscan (scanfile1, img) != EOF) {
                img = mktemp (img // "-")
                print (img, >> tmpcross)
                print (img, >> tmpdel2)
            }

            nsoffset (target, "@" // tmptarget, \
                outimages = "@" // tmpcross, outprefix = "", \
                shift = INDEF, fl_measure+, axis = INDEF, \
                pixscale = INDEF, fl_apply+, fl_all-, \
                fl_inter = l_fl_inter, logfile = l_logfile, \
                verbose = l_verbose)
            if (0 != nsoffset.status) {
                printlog ("ERROR - "//l_logname//": nsoffset failed.", \
                    l_logfile, verbose+)
                goto clean
            }

            delete (tmptarget, verify-, >& "dev$null")
            tmptarget = tmpcross
        }


        # Shift images if required

        if (l_fl_shift) {

            if (no == first) {
                if (intdbg) print ("measuring with gemoffsetlist")
                gemoffsetlist ("@" // tmptarget, target, \
                    distance = l_tolerance, \
                    age = INDEF, fl_younger+, fl_nearer+, fl_noref+, \
                    wcs_source = "phu", targetlist = "dev$null", \
                    direction = 3, offsetlist = "dev$null", \
                    key_xoff = l_key_xoff, \
                    key_yoff = l_key_yoff, logfile = l_logfile, \
                    verbose = l_debug, key_date="", key_time = "",
                    force = l_force) 
                if (gemoffsetlist.status == 1) {
                    printlog ("ERROR - "//l_logname//": gemoffsetlist \
                        could not determine the offsets.", l_logfile, \
                        verbose+) 
                    goto clean
                }
            }

            head (tmptarget, nlines = 1) | scan (img)
            phu = img // "[0]"
            keypar (phu, l_key_xoff, silent+) 
            xoff = real (keypar.value) 
            keypar (phu, l_key_yoff, silent+) 
            yoff = real (keypar.value) 

            if (first) {

                printlog (l_logname//": Unshifted images:", \
                    l_logfile, l_verbose) 
                if (l_verbose) {
                    type (tmptarget)
                    type (tmptarget, >> l_logfile)
                }

                if (intdbg) print ("getting reference offsets")
                
                xrefoff = xoff; yrefoff = yoff
                if (intdbg) print ("ref: " // xrefoff // ", " // yrefoff)

                keypar (phu, l_key_dispaxis, silent+)
                disp_axis = l_dispaxis
                if (keypar.found) disp_axis = int (keypar.value)
                if (intdbg) print ("disp_axis: " // disp_axis)

                keypar (phu, l_key_pixscale, silent+) 
                pix_scale = l_pixscale
                if (keypar.found) pix_scale = real (keypar.value) 
                if (intdbg) print ("pix_scale: " // pix_scale)

                tmpshifted = mktemp ("tmpshifted")
                print (tmpshifted, >>& tmpcomb)
                type (tmptarget, >>& tmpshifted)

            } else {

                if (intdbg) {
                    print ("shifting input")
                    print ("from " // xoff // "," // yoff)
                }
                xoff = xoff - xrefoff; yoff = yoff - yrefoff
                if (intdbg) print ("to " // xoff // "," // yoff)
                xoff = xoff / pix_scale;   yoff = yoff / pix_scale
                if (intdbg) print ("to " // xoff // "," // yoff)
                

                if ((1 == disp_axis && 0.001 < abs (xoff)) || \
                    (2 == disp_axis && 0.001 < abs (yoff))) {

                    if (1 == disp_axis) dispoff = xoff
                    else dispoff = yoff

                    printlog ("WARNING - "//l_logname//": Non-zero shift \
                        in dispersion direction: " // dispoff, \
                        l_logfile, verbose+)
                    printlog ("                     (setting to zero)", \
                        l_logfile, verbose+)

                    if (1 == disp_axis) xoff = 0.0
                    else yoff = 0.0

                }		

                if (l_fl_shiftint) {
                    xoff = nint (xoff)
                    yoff = nint (yoff)
                }


                if (intdbg) print ("generating shift output")
                tmpshiftout2 = mktemp ("tmpshiftout2")
                gemextn ("@" // tmptarget, check="", process="append", 
                    index="0", extname="", extversion="", ikparams="",
                    omit="", replace="%\[0\]%"//shft_suffix//"%", 
                    outfile=tmpshiftout2, logfile="", glogpars="",
                    verbose=l_verbose)
                type (tmpshiftout2, >>& tmpshifted)
                type (tmpshiftout2, >>& tmpdel)


                printlog (l_logname//": Images shifted by " // xoff \
                    // "," // yoff // ":", l_logfile, l_verbose) 
                if (l_verbose) {
                    type (tmptarget)
                    type (tmptarget, >> l_logfile)
                }


                if (intdbg) print ("copy PHUs")

                if (intdbg) print ("generating copy output")
                tmpshiftout = mktemp ("tmpshiftout")
                gemextn ("@" // tmpshiftout2, check="absent", 
                    process="none", index="", extname="", extversion="",
                    ikparams="", omit="", replace="", 
                    outfile=tmpshiftout, logfile="", glogpars="",
                    verbose=l_verbose)
                if (0 != gemextn.fail_count) {
                    printlog ("ERROR - "//l_logname//": Shift PHU files \
                        already exist.", l_logfile, verbose+) 
                    goto clean
                }

                if (intdbg) print ("generating copy input")
                tmpshiftin = mktemp ("tmpshiftin")
                gemextn ("@" // tmptarget, check="exists", process="expand",
                    index="0", extname="", extversion="", ikparams="", 
                    omit="", replace="", outfile=tmpshiftin, logfile="",
                    glogpars="", verbose=l_verbose)
                if (0 != gemextn.fail_count) {
                    printlog ("ERROR - "//l_logname//": Missing input \
                        files!.", l_logfile, verbose+) 
                    goto clean
                }

                tmpcopy = mktemp ("tmpcopy")
                joinlines (tmpshiftin // "," // tmpshiftout, "", \
                    output=tmpcopy, delim=" ", missing="none", \
                    maxchars=1000, shortest+, verbose-)
                scanfile1 = tmpcopy
                while (fscan (scanfile1, imgfrom, imgto) != EOF) {
                    if (intdbg) print (imgfrom // " --> " // imgto)
                    imcopy (imgfrom, imgto, verbose-, >& "dev$null")
                    gemextn (imgto, proc="none", check="", omit="params", 
                        outfile="STDOUT", logfile="", glogpars="",
                        verbose-) | scan (imgto)
                    if (intdbg) print ("setting in " // imgto)
                    gemhedit (imgto, "NSCHLX", xoff, "x shift (pix)")
                    gemhedit (imgto, "NSCHLY", yoff, "y shift (pix)")
                }
                delete (tmpcopy, verify-, >& "dev$null")


                if (intdbg) print ("copy MDFs")

                if (intdbg) print ("generating copy output")
                delete (tmpshiftout, verify-, >& "dev$null")
                tmpshiftout = mktemp ("tmpshiftout")
                gemextn ("@" // tmpshiftout2, check="absent", 
                    process="append", index="", extname="MDF", 
                    extversion="", ikparams="append", omit="", replace="", 
                    outfile=tmpshiftout, logfile="", glogpars="",
                    verbose=l_verbose)
                if (0 != gemextn.fail_count) {
                    printlog ("ERROR - "//l_logname//": Shift MDF files \
                        already exist.", l_logfile, verbose+) 
                    goto clean
                }

                if (intdbg) print ("generating copy input")
                delete (tmpshiftin, verify-, >& "dev$null")
                tmpshiftin = mktemp ("tmpshiftin")
                gemextn ("@" // tmptarget, check="exists", process="expand",
                    index="", extname="MDF", extversion="", ikparams="", 
                    omit="", replace="", outfile=tmpshiftin, logfile="",
                    glogpars="", verbose=l_verbose)
                if (0 != gemextn.fail_count) {
                    printlog ("ERROR - "//l_logname//": Missing input \
                        files!.", l_logfile, verbose+) 
                    goto clean
                }

                tmpcopy = mktemp ("tmpcopy")
                joinlines (tmpshiftin // "," // tmpshiftout, "", \
                    output=tmpcopy, delim=" ", missing="none", \
                    maxchars=1000, shortest+, verbose-)
                if (intdbg) type (tmpcopy)
                scanfile1 = tmpcopy
                while (fscan (scanfile1, imgfrom, imgto) != EOF) {
                    if (intdbg) print (imgfrom // " --> " // imgto)
                    tcopy (imgfrom, imgto, verbose-, >& "dev$null")
                }
                delete (tmpcopy, verify-, >& "dev$null")


                if (intdbg) print ("generating file names for shift")
                delete (tmpshiftout, verify-, >& "dev$null")
                tmpshiftout = mktemp ("tmpshiftout")
                gemextn ("@" // tmptarget, check="absent", process="expand",
                    index="", extname=extname, extversion="1-", 
                    ikparams="append", omit="", 
                    replace="%\[%"//shft_suffix//"[%", outfile=tmpshiftout,
                    logfile="", glogpars="", verbose=l_verbose)

                delete (tmpshiftin, verify-, >& "dev$null")
                tmpshiftin = mktemp ("tmpshiftin")
                gemextn ("@" // tmptarget, check="", process="expand", \
                    index="", extname=extname, extversion="1-", \
                    ikparams="", omit="", replace="", outfile=tmpshiftin, \
                    logfile="", glogpars="", verbose=l_verbose)


                if (intdbg) {
                    print ("shifting...")
                    type (tmpshiftin)
                    print ("by " // xoff // "," // yoff // " to...")
                    type (tmpshiftout)
                }
                imshift ("@" // tmpshiftin, "@" // tmpshiftout, \
                    xshift=xoff, yshift=yoff, shifts_file="", \
                    interp_type=l_interptype, boundary_type=l_boundary, \
                    constant=l_constant)


                delete (tmpshiftin, verify-, >& "dev$null")
                delete (tmpshiftout, verify-, >& "dev$null")


                if (no == l_fl_shiftint) {

                    if (intdbg) print ("peaking dq")
                    tmpshiftin = mktemp ("tmpshiftin")
                    gemextn ("@" // tmptarget, check="exists", \
                        process="expand", index="", extname=l_dq_ext, \
                        extversion="1-", ikparams="", omit="", \
                        replace="%\[%"//shft_suffix//"[%", \
                        outfile=tmpshiftin, logfile="dev$null", glogpars="",
                        verbose-)
                    tmpshiftout = mktemp ("tmpshiftout")
                    gemextn ("@" // tmptarget, check="exists", \
                        process="expand", index="", extname=l_dq_ext, \
                        extversion="1-", ikparams="overwrite", omit="", \
                        replace="%\[%"//shft_suffix//"[%", \
                        outfile=tmpshiftout, logfile="dev$null",
                        glogpars="", verbose-)

                    scanfile1 = tmpshiftin
                    scanfile2 = tmpshiftout
                    while (fscan (scanfile1, imgin) != EOF) {
                        junk = fscan (scanfile2, imgout)

                        if (intdbg) print (imgin // " --> " // imgout)
                        imexpr ("a > 0 ? 1 : 0", imgout, imgin, \
                            dims="auto", intype="auto", outtype="short", \
                            refim="auto", bwidth=0, btype="nearest", \
                            bpixval=0, rangecheck+, verbose+)
                    }

                    delete (tmpshiftin, verify-, >& "dev$null")
                    delete (tmpshiftout, verify-, >& "dev$null")
                }

                delete (tmpshiftout2, verify-, >& "dev$null")
            }

            delete (tmptarget, verify-, >& "dev$null")

        } else {
            print (tmptarget, >>& tmpcomb)
        }

        first = no
    }


    # Second processing step - group with gemcombine using o/p 
    # files as required

    scanfile1 = tmpcomb
    scanfile2 = tmpout
    scanfile3 = tmpnrej

    while (fscan (scanfile1, combfile) != EOF) {

        if (intdbg) {
            print ("will now combine:")
            type (combfile)
        }


        if (intdbg) print ("checking output file")
        if (dynout) {
            head (combfile, nlines = 1) | scan (imgout)
            imgout = imgout // output_suffix
        } else {
            if (fscan (scanfile2, imgout) == EOF) {
                printlog ("ERROR - "//l_logname//": Missing output \
                    file.", l_logfile, verbose+) 
                goto clean
            }
        }
        if (intdbg) print (imgout)
        gemextn (imgout, check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", omit="", replace="",
            outfile="dev$null", logfile="", glogpars="", verbose=l_verbose)
        if (0 != gemextn.fail_count || 1 != gemextn.count) {
            printlog ("ERROR - "//l_logname//": Bad/existing output \
                file.", l_logfile, verbose+) 
            goto clean
        }


        if (intdbg) print ("checking nrejfile")
        imgnrej = ""
        if ("" != l_nrejfile) {
            if (fscan (scanfile3, imgnrej) == EOF) {
                printlog ("WARNING - "//l_logname//": Missing nrej \
                    file.", l_logfile, verbose+) 
            } else {
                if (intdbg) print (imgnrej)
                gemextn (imgnrej, check="absent", process="none", \
                    index="", extname="", extversion="", ikparams="", \
                    omit="", replace="", outfile="dev$null", \
                    logfile="", glogpars="", verbose=l_verbose)
                if (0 != gemextn.fail_count || 1 != gemextn.count) {
                    printlog ("ERROR - "//l_logname//": Bad/existing nrej \
                        file.", l_logfile, verbose+) 
                    goto clean
                }
            }
        }


        if (no == l_fl_shift) {
            printlog (l_logname//": Grouping:", l_logfile, l_verbose) 
            if (l_verbose) {
                type (combfile)
                type (combfile, >> l_logfile)
            }
        }


        if (intdbg) 	print ("combining...")
        gemcombine ("@" // combfile, imgout, title = "", \
            logfile = l_logfile, combine = l_combtype, \
            reject = l_rejtype, masktype=l_masktype, maskvalue=l_maskvalue, \
            scale = l_scale, zero = l_zero, \
            weight = l_weight, statsec = l_statsec, \
            expname = l_key_exptime, lthreshold = l_lthreshold, \
            hthreshold = l_hthreshold, \
            nlow = l_nlow, nhigh = l_nhigh, nkeep = l_nkeep, \
            mclip = l_mclip, lsigma = l_lsigma, \
            hsigma = l_hsigma, key_ron = l_key_ron, \
            key_gain = l_key_gain, ron = l_ron, gain = l_gain, \
            snoise = l_snoise, sigscal = l_sigscale, \
            pclip = l_pclip, grow = l_grow, bpmfile = l_bpm, \
            nrejfile = imgnrej, sci_ext = l_sci_ext, \
            var_ext = l_var_ext, dq_ext = l_dq_ext, fl_vardq = l_fl_vardq, 
            offsets = "none", verbose = l_debug) 

        if (0 != gemcombine.status) {
            printlog ("ERROR - "//l_logname//": Error in gemcombine.", \
                l_logfile, verbose+) 
            goto clean
        }



        if (intdbg) print ("fixing header stuff")

        gemdate ()
        phu = imgout // "[0]"
        gemhedit (phu, "GEM-TLM", gemdate.outdate, 
            "UT Last modification with GEMINI") 

        # hedit(phu,"GXOFF",add-,addonly-,delete+,verify-,show-,update+)
        # hedit(phu,"GYOFF",add-,addonly-,delete+,verify-,show-,update+)
        # hedit(phu,"GOFFREF",add-,addonly-,delete+,verify-,show-,update+)

        keywdname = substr(l_logname, 1, 8)
        keywdname = strupr (keywdname)
        printf ("%-8s= \'%-18s\' / %-s\n", keywdname, gemdate.outdate, \
            "UT Time stamp for "//l_logname//"", >> tmphead)
        printf ("%-8s= %20.3f / %-s\n", "NSCHLTOL", l_tolerance, \
            "Spatial tolerance for "//l_logname//"", >> tmphead) 
        printf ("%-8s= \'%-18s\' / %-s\n", "NSCHLSTA", l_statsec, \
            "Statistics section used by "//l_logname//"", >> tmphead) 
        printf ("%-8s= \'%-18b\' / %-s\n", "NSCHLSFT", l_fl_shiftint, \
            "Use integer pixel shifts?", >> tmphead) 
        if (no ==l_fl_shiftint) {
            printf ("%-8s= \'%-18s\' / %-s\n", "NSCHLINT", l_interptype, \
                "Pixel interpolation method", >> tmphead) 
        }
        printf ("%-8s= \'%-18s\' / %-s\n", "NSCHLBND", l_boundary, \
            "Boundary type for shifting", >> tmphead) 
        if (l_boundary == "constant") {
            printf ("%-8s= %20.3f / %-s\n", "NSCHLCON", l_constant, \
                "Constant for boundary extension for shifting", \
                >> tmphead) 
        }
        printf ("%-8s= \'%-18s\' / %-s\n", "NSCHLCOM", l_combtype, \
            "Combine method used by "//l_logname//"", >> tmphead) 
        if (l_zero != "none") {
            printf ("%-8s= \'%-18s\' / %-s\n", "NSCHLZER", l_zero, \
                "Statistic used to compute additive offsets", >> tmphead) 
        }
        if (l_scale != "none") {
            printf ("%-8s= \'%-18s\' / %-s\n", "NSCHLSCL", l_scale, \
                "Statistic used to compute scale factors", >> tmphead) 
        }
        if (l_weight != "none") {
            printf ("%-8s= \'%-18s\' / %-s\n", "NSCHLWEI", l_weight, \
                "Statistic used to compute relative weights", >> tmphead) 
        }
        if (l_rejtype != "none") {
            printf ("%-8s= \'%-18s\' / %-s\n", "NSCHLREJ", l_rejtype, \
                "Rejection algorithm used by "//l_logname//"", >> tmphead) 
            printf ("%-8s= \'%-18s\' / %-s\n", "NSCHLLTH", str(l_lthreshold), \
                "Lower threshold before combining", >> tmphead) 
            printf ("%-8s= \'%-18s\' / %-s\n", "NSCHLHTH", str(l_hthreshold), \
                "Upper threshold before combining", >> tmphead) 
            printf ("%-8s= %20.2f / %-s\n", "NSCHLGRW", l_grow, \
                "Radius for additional pixel rejection", >> tmphead) 
        }
        if (l_rejtype == "minmax") {
            printf ("%-8s= %20.0f / %-s\n", "NSCHLNLO", l_nlow, \
                "Low pixels rejected (minmax)", >> tmphead) 
            printf ("%-8s= %20.0f / %-s\n", "NSCHLNHI", l_nhigh, \
                "High pixels rejected (minmax)", >> tmphead) 
        }
        if (l_rejtype == "sigclip" || l_rejtype == "avsigclip" \
            || l_rejtype == "ccdclip" || l_rejtype == "crreject" \
            || l_rejtype == "pclip") {
            
            printf ("%-8s= %20.2f / %-s\n", "NSCHLLSI", l_lsigma, \
                "Lower sigma for rejection", >> tmphead) 
            printf ("%-8s= %20.2f / %-s\n", "NSCHLHSI", l_hsigma, \
                "Upper sigma for rejection", >> tmphead) 
            printf ("%-8s= %20.0f / %-s\n", "NSCHLNKE", l_nkeep, \
                "Min(max) number of pixels to keep(reject)", >> tmphead) 
            printf ("%-8s= \'%-18b\' / %-s\n", "NSCHLMCL", l_mclip, \
                "Use median in clipping algorithms", >> tmphead) 
        } 
        if (l_rejtype == "sigclip" || l_rejtype == "avsigclip" \
            || l_rejtype == "ccdclip" || l_rejtype == "crreject") {
            
            printf ("%-8s= %20.2f / %-s\n", "NSCHLSSC", l_sigscale, \
                "Tolerance for sigma clip scaling correction", >> tmphead) 
        }
        if (l_rejtype == "pclip") {
            printf ("%-8s= %20.2f / %-s\n", "NSCHLPCL", l_pclip, \
                "Percentile clipping factor used by pclip", >> tmphead) 
        }     
        if (l_rejtype == "ccdclip") {
            printf ("%-8s= \'%-18s\' / %-s\n", "NSCHLSNO", l_snoise, \
                "Sensitivity noise (e) used by ccdclip", >> tmphead) 
        }

        index = 1
        scanfile4 = combfile
        while (fscan (scanfile4, img) != EOF) {
            imgname = img
            shfted = \
                substr (img, strlen (img) -5, strlen (img)) == "_shift"
            if (shfted) imgname = substr (img, 1, strlen (img) - 6) 
            printf ("%-5s%03d= \'%-18s\' / %-s\n", "NSCHL", index, imgname, \
                "Image " // index // " combined with "//l_logname//"", \
                >> tmphead) 
            if (l_fl_shift) {
                xoff = 0
                yoff = 0
                if (shfted) {
                    if (intdbg) print ("trying " // img // "[0]")
                    hselect (img // "[0]", "NSCHLX", yes) | scan (xoff)
                    hselect (img // "[0]", "NSCHLY", yes) | scan (yoff)
                }
                printf ("%-6s%02d= %20.2f / %-s\n", "NSCHLX", index, xoff, \
                    "x shift (pix) for input image " // index, >> tmphead) 
                printf ("%-6s%02d= %20.2f / %-s\n", "NSCHLY", index, yoff, \
                    "y shift (pix) for input image " // index, >> tmphead)
            }
            index += 1
        }

        if (intdbg) type (tmphead)
        mkheader (phu, tmphead, append+, verbose-)
        delete (tmphead, verify-, >& "dev$null")
        tmphead = mktemp ("tmphead")

    }

    status = 0

clean:

    if (access (tmpcomb)) {
        scanfile1 = tmpcomb
        while (fscan (scanfile1, combfile) != EOF) {
            delete (combfile, verify-, >& "dev$null")
        }
    }

    if (access (tmpdel)) {
        if (no == l_fl_keepshift)
            imdelete ("@" // tmpdel, verify-, >& "dev$null")
        delete (tmpdel, verify-, >& "dev$null")
    }
    if (access (tmpdel2)) {
        imdelete ("@" // tmpdel2, verify-, >& "dev$null")
        delete (tmpdel2, verify-, >& "dev$null")
    }
    delete (tmpcomb, verify-, >& "dev$null")
    delete (tmpout, verify-, >& "dev$null")
    delete (tmptarget, verify-, >& "dev$null")

    scanfile1 = ""
    scanfile2 = ""
    scanfile3 = ""
    scanfile4 = ""

    if (status == 0) 
        printlog (l_logname//": Exit status good.", l_logfile, l_verbose) 
    else 
        printlog (l_logname//": Exit status error.", l_logfile, l_verbose) 
    printlog ("-------------------------------------------------------\
        -------------------------", l_logfile, l_verbose) 

end
