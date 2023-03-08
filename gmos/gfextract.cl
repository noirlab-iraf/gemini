# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gfextract (inimage)

# Extract GMOS IFU spectra
#
# Version   Sept 20, 2002 BM v1.4 release
#           Aug 25, 2003  KL IRAF2.12 - new parameter, addonly, in hedit
#           Nov 12, 2003  BM generalized for GMOS-S
#           Feb 29, 2004  BM added thresh param, fix determination of
#                            slit mode for N&S, add gnsskysub option
#           Mar 19, 2004  BM Significant improvement to algorithm for finding
#                            image sections to extract, params added
#           Mar 25, 2004  BM Fix refpix for blue one-slit mode
#           Apr 06, 2004  BM Call gnsskysub before gmosaic
#           Apr 07, 2004  BM fix locations of mktemp
#           Apr 08, 2004  BM generalize matching of open filters
#           Sep 25, 2008  JT solve problem calculating wavelength regions and
#                            add the exslit option for odd setups
#           Mar 15, 2013  JT Propagate VAR/DQ instead of using apall stats
#           Sep 27, 2013  JT/MS Handle multiple DQ planes correctly
#           Jun 02, 2014  JT Option to re-interpolate & grow chip gaps

string inimage   {prompt="Input image"}
string outimage  {"",prompt="Output file"}
string outpref   {"e",prompt="Prefix for output image"}
string title     {"",prompt="Title for output SCI plane"}
string reference {"",prompt="Reference file"}
string response  {"",prompt="Fiber response file"}
string exslits   {"*",enum="red|blue|*",prompt="Which slit(s) to extract (red,blue,*)"}
int    line      {INDEF,prompt="Line/column for finding apertures"}
int    nsum      {10,min=1,prompt="Number of columns to use for finding apertures"}
bool   trace     {yes,prompt="Trace spectra?"}
bool   recenter  {yes,prompt="Recenter apertures?"}
real   thresh    {200.,prompt="Detection threshold for profile centering"}
string function  {"chebyshev",enum="chebyshev|spline1|spline3|legendre",prompt="Function for trace"}
int    order     {5,min=1,prompt="Order of trace fit"}
int    t_nsum    {10,min=1,prompt="Number of dispersion lines to sum for the trace"}
string weights   {"variance",enum="variance|none",prompt="Weighting during extraction"}
string bpmfile   {"gmos$data/chipgaps.dat",prompt="Bad pixel mask for column interpolation"}
real   grow      {1.0, min=0, prompt="Gap growth radius in pixels"}
string  gaindb  {"default", prompt="Database with gain data"}
string gratingdb  {"gmos$data/GMOSgratings.dat",prompt="Gratings database file"}
string filterdb   {"gmos$data/GMOSfilters.dat",prompt="Filters database file"}
real   xoffset    {INDEF,prompt="X offset in wavelength [nm]"}
real   perovlap   {10.,prompt="Percentage by which to shrink overlapping spectra"}
string sci_ext    {"SCI",prompt="Name of science extension"}
string var_ext    {"VAR",prompt="Name of variance extension"}
string dq_ext     {"DQ",prompt="Name of data quality extension"}
bool   fl_inter   {no,prompt="Run APALL interactively?"}
bool   fl_vardq   {no,prompt="Extract variance and DQ planes?"}
bool   fl_novlap  {yes,prompt="Avoid spectral overlap?"}
bool   fl_gnsskysub {no,prompt="Run GNSSKYSUB?"}
bool   fl_fixnc     {no,prompt="Auto-correct for nod count mismatch?"}
bool   fl_fixgaps {no,prompt="Re-interpolate chip gaps after extraction?"}
bool   fl_gsappwave {yes,prompt="Run GSAPPWAVE?"}
bool   fl_fulldq    {no, prompt="Decompose DQ during gmosaic transformation; apply correct handling?"}            # OLDP-3
real   dqthresh     {0.1, min=0.0, max=0.5, prompt="Threshold applied to DQ when fl_fulldq=yes"}
string logfile    {"",prompt="Logfile"}
bool   verbose    {yes,prompt="Verbose output?"}
int    status     {0,prompt="Exit status (0=good)"}

begin

    string l_inimage, l_outimage, l_prefix, l_title, l_logfile, l_response
    string l_exslits, l_grating, l_filter[2], l_gratingdb, l_filterdb
    string l_reference, l_bpmfile, l_gaindb, l_function, l_weights
    string l_sci_ext, l_var_ext, l_dq_ext
    string imgsec, varsec, dqsec, slitreg[2], l_detsec[2], mdf, impaste
    string sci, var, dq, aptable, mdfslit, nsskysub, dqtype, tmpmef1, tmpmef2
    string gname, fname, ffile, outroot, refroot, outlog
    string sr_logfile, l_key_qecorrim
    string tmpin, tmpslit, tmptmpslit, tmpvar, tmpdq, tmpdqsec, tmpplane
    string tmpplane2, tmpaccum
    int    l_xbin, l_ybin, l_line, l_nsum, l_t_nsum, l_order, l_naxis2
    int    lmm, sx, sy, wx, width, center, xccd[2], yccd[2], slitnum
    int    sposy, slitsepp, apid, apno, detx1, detx2, dety1, dety2
    real   l_wave, blambda, R, coverage, nmppx, scale[2,3], asecmm, l_xoffset
    real   gwave1, gwave2, fwave1, fwave2, wave1, wave2, wmin, wmax, wavoffset
    real   wmn1, wmn2, wmx1, wmx2, l_perovlap, slitl1, slitl2, l_grow
    real   tilt, pi, a, greq, refpix[2], rempix[2], refmin, remmin, l_thresh
    real   detector_upper_spec_limit[3], l_dqthresh, a_lower, a_upper, a_width
    real   a_radius, a_minsep, a_maxsep, sfact
    int    nslit, j1, id1, nr, x1[2], x2[2], i, j, k, len, nx[2], n
    int    nextnd, inst, nxmax, novlap, chipgap, chipgap_bin, iccd
    int    dqmin, dqmax, ndq, ymin, ymax
    bool   l_verbose, l_inter, l_recenter, l_trace, delpaste, l_fl_vardq
    bool   l_fl_gsappwave, l_fl_fixpix, l_fl_gnsskysub, l_fl_fixnc, l_fl_fixgaps
    bool   l_fl_novlap, in_qestate, resp_qestate, l_fl_fulldq
    struct sdate, msg

    # Query parameters
    l_inimage=inimage ; l_outimage=outimage ; l_title=title ; l_prefix=outpref
    l_response=response ; l_reference=reference ; l_exslits=exslits 
    l_bpmfile=bpmfile ; l_grow=grow ; l_line=line; l_logfile=logfile
    l_verbose=verbose ; l_inter=fl_inter ; l_fl_vardq=fl_vardq
    l_gratingdb=gratingdb ; l_filterdb=filterdb ; l_trace=trace; l_nsum=nsum
    l_t_nsum=t_nsum ; l_order=order ; l_function=function; l_weights=weights
    l_fl_gsappwave=fl_gsappwave ; l_recenter=recenter ; l_thresh=thresh
    l_sci_ext=sci_ext ; l_var_ext=var_ext ; l_dq_ext=dq_ext
    l_xoffset=xoffset ; l_fl_gnsskysub=fl_gnsskysub ; l_fl_fixnc=fl_fixnc
    l_fl_fixgaps=fl_fixgaps; l_fl_novlap=fl_novlap ; l_perovlap=perovlap/100.
    l_gaindb=gaindb

    # DQ handling
    l_fl_fulldq = fl_fulldq
    l_dqthresh = dqthresh

    status=0
    pi=3.1415927
    asecmm=1.611444
    delpaste=yes
    # Pixel scales scale[inst,iccd]
    # inst: 1 - GMOS-N, 2 - GMOS-S
    # iccd: 1 - Current CCDs (2011-07), 2 - New EEV2 CCDs, 3 - Hamamatsu CCDs
    scale[1,1]=0.0727
    scale[1,2]=0.07288 ##M
    scale[1,3]=0.0807
    scale[2,1]=0.073
    scale[2,3]=0.0800 ##M PIXEL_SCALE
    # Upper (red) limit of the given detector type in [nm].
    # detector_upper_spec_limit[iccd] 
    detector_upper_spec_limit[1] = 1025.0 # EEV CCDs
    detector_upper_spec_limit[2] = 1050.0 # e2vDD CCDs ##M
    detector_upper_spec_limit[3] = 1080.0 # Hamamatsu CCDs ## KL

    # Set QE correction image keyword
    l_key_qecorrim = "QECORRIM"

    # Keep imgets parameters from changing by outside world
    cache("imgets", "fparse", "tinfo", "tabpar", "specred", "gemdate")

    # Define temporary files
    #imgsec=mktemp("tmpimgsec") #tmp name not used
    mdf=mktemp("tmpmdf")
    impaste=mktemp("tmppaste")
    nsskysub=mktemp("tmpnsskysub")
    aptable=mktemp("tmpaptable")
    mdfslit=mktemp("tmpmdfslit")//".fits"
    tmpin = mktemp ("tmpin")
    tmpdqsec = mktemp ("tmpdqsec")//".fits"
    tmpplane = mktemp ("tmpplane")//".fits"
    tmpplane2 = mktemp ("tmpplane2")//".fits"
    tmpaccum = mktemp ("tmpaccum")//".fits"
    # Temporary files to be defined later
    sci=""
    var=""
    dq=""

    # Start logging to file
    if (l_logfile == "STDOUT") {
        l_logfile = ""
        l_verbose=yes
    }

    # Test the logfile:
    gemlogname (logpar=l_logfile, package="gmos")
    if (gemlogname.status != 0)
        goto error
    l_logfile=gemlogname.logname

    # Start logfile
    date | scan(sdate)
    printlog("---------------------------------------------------------\
        -----------------------", l_logfile, verbose=l_verbose)
    printlog("GFEXTRACT -- "//sdate, l_logfile, verbose=l_verbose)
    printlog(" ", l_logfile, verbose=l_verbose)
    printlog("inimage   = "//l_inimage, l_logfile, verbose=l_verbose)
    printlog("outimage  = "//l_outimage, l_logfile, verbose=l_verbose)
    printlog("outpref   = "//l_prefix, l_logfile, verbose=l_verbose)
    printlog("title     = "//l_title, l_logfile, verbose=l_verbose)
    printlog("response  = "//l_response, l_logfile, verbose=l_verbose)
    printlog("reference = "//l_reference, l_logfile, verbose=l_verbose)
    printlog("exslits   = "//l_exslits, l_logfile, verbose=l_verbose)
    printlog("trace     = "//l_trace, l_logfile, verbose=l_verbose)
    printlog("function  = "//l_function, l_logfile, verbose=l_verbose)
    printlog("order     = "//l_order, l_logfile, verbose=l_verbose)
    printlog("weights   = "//l_weights, l_logfile, verbose=l_verbose)
    printlog("bpmfile   = "//l_bpmfile, l_logfile, verbose=l_verbose)
    printlog("gratingdb = "//l_gratingdb, l_logfile, verbose=l_verbose)
    printlog("filterdb  = "//l_filterdb, l_logfile, verbose=l_verbose)
    if (l_xoffset != INDEF) {
        printlog("xoffset   = "//str(l_xoffset), l_logfile, verbose=l_verbose)
    } else {
        printlog("xoffset   = INDEF", l_logfile, verbose=l_verbose)
    }
    printlog("fl_vardq  = "//l_fl_vardq, l_logfile, verbose=l_verbose)
    printlog(" ", l_logfile, verbose=l_verbose)

    #check that there are input files
    if (l_inimage == "" || l_inimage == " "){
        printlog("ERROR - GFEXTRACT: input files not specified", l_logfile, \
            verbose+)
        goto error
    }

    # check existence of list file
    if (substr(l_inimage, 1, 1) == "@") {
        printlog("ERROR - GFEXTRACT: lists are currently not supported",
            l_logfile, verbose+)
        goto error
    }

    # check existence of input file
    gimverify(l_inimage)
    if (gimverify.status != 0) {
        printlog("ERROR - GFEXTRACT: "//l_inimage//" does not exist or is \
            not a MEF", l_logfile, verbose+)
        goto error
    }
    l_inimage=gimverify.outname//".fits"

    #check that an output file is given
    if ((l_outimage == "" || l_outimage == " ") && (l_prefix=="" || \
        l_prefix==" ")) {
        printlog("ERROR - GFEXTRACT: output file not specified", l_logfile, \
            verbose+)
        goto error
    }

    # Check if using prefix
    if (l_outimage == "" || l_outimage == " ")
        l_outimage=l_prefix//l_inimage

    #check output name
    fparse(l_outimage, verbose-)
    outroot=fparse.root
    if (fparse.extension == "")
        l_outimage=l_outimage//".fits"

    # check response file if not ""
    if (l_response != "" && l_response != " ") {
        gimverify(l_response)
        if (gimverify.status != 0) {
            printlog("ERROR - GFEXTRACT: response image does not exist or is \
                not a MEF.", l_logfile, verbose+)
            goto error
        }
        l_response=gimverify.outname//".fits"
        imgets(l_response//"[0]", "GFRESPON", >>& "dev$null")
        if (imgets.value == "0") {
            printlog("ERROR - GFEXTRACT: response image has not been \
                processed by GFRESPONSE.", l_logfile, verbose+)
            goto error
        }

        # Check the QE correction state of the response file
        keypar (l_inimage//"[0]", l_key_qecorrim, silent+)
        if (keypar.found) {
            resp_qestate = yes
        } else {
            resp_qestate = no
        }

    } else {
        l_response=""
    }

    # check bpmfile file if not ""
    if (l_bpmfile != "" && l_bpmfile != " ") {
        if (!access(l_bpmfile)) {
            printlog("ERROR - GFEXTRACT: bpmfile image does not exist.", 
                l_logfile, verbose+)
            goto error
        }
        l_fl_fixpix=yes
    } else {
        l_fl_fixpix=no
    }

    # reference
    fparse(l_reference, verbose-)
    refroot=fparse.root

    # Check that the output file does not already exist. If so, exit.
    if (imaccess(l_outimage)) {
        printlog("ERROR - GFEXTRACT: Output file "//l_outimage//" already \
            exists.", l_logfile, verbose+)
        goto error
    }

    # Parse the specified slit into a number to match existing logic below
    # (names were used for consistency with gfreduce)
    if (l_exslits=="red") 
        slitnum=1
    else if (l_exslits=="blue") 
        slitnum=2
    else 
        slitnum=0

    # Determine which instrument
    imgets(l_inimage//"[0]", "INSTRUME", >>& "dev$null")
    if (imgets.value == "0") {
        printlog ("ERROR - GFEXTRACT: Instrument keyword not found.", 
            l_logfile, verbose+)
        goto error
    }
    inst=1 # Default is GMOS-N, support for old data
    if (imgets.value == "GMOS-S")
        inst=2

    # Deteremine the dectector type used
    imgets(l_inimage//"[0]", "DETTYPE", >>& "dev$null")
    if (imgets.value == "0") {
        printlog ("ERROR - GFEXTRACT: Instrument keyword not found.", 
            l_logfile, verbose+)
        goto error
    }
    iccd=1 # Default is old EEV2 for GMOS-N and GMOS-S
    # Chipgaps in pixels
    chipgap = 37 # Default is old EEV2 for GMOS-N and GMOS-S
    chipgap_bin = 36 # Default is old EEV2 for GMOS-N and GMOS-S. Binned data
    if (imgets.value == "SDSU II e2v DD CCD42-90") { # New e2vDD CCDs
        iccd = 2 
        chipgap = 37
        chipgap_bin = 36 # For use with binned data
    } else if (imgets.value == "S10892") { # Hamamatsu CCDs SOUTH
        iccd = 3 
        chipgap = 61 ##M CHIP_GAP
        chipgap_bin = 60 # For use with binned data ##M CHIP_GAP
    } else if (imgets.value == "S10892-N") { # Hamamatsu NORTH
        iccd = 3 
        chipgap = 67 ## KL CHIP_GAP
        chipgap_bin = 67 # For use with binned data ## KL CHIP_GAP
    }

    # slice out the MDF
    tcopy(l_inimage//"[MDF]", mdf//".fits", verbose-)
    if (!access(mdf//".fits")) {
        printlog("ERROR - GFEXTRACT: Input image does not contain a MDF", 
            l_logfile, verbose+)
        goto error
    }

    # Find slits
    # get header information
    imgets(l_inimage//"[0]", "GRATING")
    l_grating=imgets.value
    if (l_grating=="MIRROR") {
        printlog("ERROR - GFEXTRACT: input image is not dispersed.", 
            l_logfile, l_verbose)
        goto error
    }

    # Obtain the NAXIS2 value - used to write the DETSEC keyword
    # which is the area of the mosaiced image that the slits are extracted from
    imgets(l_inimage//"[sci,1]", "i_naxis2")
    print(imgets.value) | scan(l_naxis2)

    imgets(l_inimage//"[sci,1]", "CCDSUM")
    print(imgets.value) | scan(l_xbin, l_ybin)
    
    imgets(l_inimage//"[0]", "FILTER1")
    l_filter[1]=imgets.value

    imgets(l_inimage//"[0]", "FILTER2")
    l_filter[2]=imgets.value

    imgets(l_inimage//"[0]", "GRWLEN")
    l_wave=real(imgets.value)

    # Read the QE correction state of the input image
    keypar (l_inimage//"[0]", l_key_qecorrim, silent+)
    if (keypar.found) {
        in_qestate = yes
    } else {
        in_qestate = no
    }

    # get grating information
    match(l_grating, l_gratingdb, stop-, print+, meta+) | \
        scan(gname, lmm, blambda, R, coverage, gwave1, gwave2, wavoffset)

    if (l_xoffset != INDEF) 
        wavoffset=l_xoffset

    #l_wave=l_wave+wavoffset
    print(((l_wave*lmm)/1.e6)) | interp("gmos$data/gratingeq.dat", "STDIN", 
        int_mode="spline", curve_gen-) | scan(greq, tilt)

    tilt=tilt*pi/180.
    a=sin(tilt+0.872665)/sin(tilt)
    R=206265.*greq/(0.35*81.0*sin(tilt))
    nmppx=a*scale[inst,iccd]*real(l_xbin)*l_wave/(R*0.35)
    coverage=nmppx*(6144.+real(2*chipgap))/real(l_xbin)
    slitsepp=175.*asecmm/(scale[inst,iccd]*a)

    wave1=gwave1
    wave2=gwave2

    printlog("Grating: "//gname, l_logfile, l_verbose)
    printlog("Resolution (0.35'' slit) = "//R, l_logfile, l_verbose)
    printlog("nm/pix = "//nmppx, l_logfile, l_verbose)
    printlog("Central wavelength = "//l_wave, l_logfile, l_verbose)

    # which slit?
    tinfo(mdf//".fits", ttout-)
    nr=tinfo.nrows
    tabpar(mdf//".fits", "NO", 1)
    id1=int(tabpar.value)
    nslit=1
    if ((id1==1 && nr>750) || (id1==51 && nr>350)) {
        nslit=2
    }
    j1=1
    if (id1==751) {
        j1=2
        nslit=2
    }
    # JT: added with the parameter "slit":
    if (nslit==j1) {
        if (slitnum!=0 && slitnum!=j1) {
            printlog("ERROR - GFEXTRACT: exslits="//l_exslits//" requested \
                but the MDF matches the other slit", l_logfile, l_verbose)
            goto error
        }
    } else if (slitnum!=0) { # 2-slit MDF and slit!="*"
        # Here we alter the 2-slit input MDF to correspond to only 1 slit,
        # since we're only extracting one of them! However, we don't set
        # nslit and j1 to slitnum until we have done the 2-slit mode overlap
        # calculations below...
        tselect(mdf//".fits", mdfslit,
            "NO>="//(750*(slitnum-1)+1)//" && NO<="//(750*slitnum))
        imdelete(mdf//".fits", verify-)
        imrename(mdfslit, mdf//".fits", verbose-)
    }

    # get filter information
    fwave1=0.0 ; wmn1=0.0 ; wmn2=0.0
    fwave2=9999.0 ; wmx1=99999.0 ; wmx2=99999.0
    if (l_filter[1] != "" && substr(l_filter[1], 1, 4) != "open")
        match(l_filter[1], l_filterdb, stop-, print+, meta+) | \
            scan(fname, wmn1, wmx1, ffile)

    if (l_filter[2] != "" && substr(l_filter[2], 1, 4) != "open")
        match(l_filter[2], l_filterdb, stop-, print+, meta+) | \
            scan(fname, wmn2, wmx2, ffile)

    if (wmn1 > wmn2) {
        fwave1=wmn1
    } else {
        fwave1=wmn2
    }
    if (wmx1 < wmx2) {
        fwave2=wmx1
    } else {
        fwave2=wmx2
    }

    # determine whether filter or grating limits wavelength coverage
    if (fwave1 > wave1)
        wave1 = fwave1
    if (fwave2 < wave2)
        wave2 = fwave2

    # and red limit of CCD
    if (wave2 > detector_upper_spec_limit[iccd])
        wave2 = detector_upper_spec_limit[iccd]

    # length of spectrum in x
    wmin=wave1
    wmax=wave2
    # in pixels
    sx=nint((wmax-wmin)/nmppx)

    printlog("Filter1: "//l_filter[1], l_logfile, l_verbose)
    printlog("Filter2: "//l_filter[2], l_logfile, l_verbose)
    printlog("Max wavelength coverage: "//wmin//" "//wmax, l_logfile, \
        l_verbose)
    #printlog("Max spectrum length = "//sx, l_logfile, l_verbose)

    # pixel value of central wavelength from left (red) end of spectrum
    wx=sx-nint((l_wave-wmin)/nmppx)
    #printlog("Pixel of central wavelength: "//wx, l_logfile, l_verbose)

    # positions of slits on CCD
    yccd[1]=1
    yccd[2]=1
    xccd[1]=((3072+chipgap)-slitsepp/2)/l_xbin+wavoffset/nmppx
    xccd[2]=xccd[1]+slitsepp/l_xbin
    printlog("Slit separation: "//slitsepp, l_logfile, l_verbose)
    printlog("Positions of slits: "//xccd[1]//" "//xccd[2], l_logfile, \
        l_verbose)

    #if (sx > slitsep && nslit-j1==1)
        #sx=slitsepp

    # find the locations of the two slits ignoring overlap
    for (i=1; i<=2; i+=1) {
        refpix[i]=real(wx)
        x1[i]=xccd[i]-wx+1
        if (x1[i] < 1) {
            refpix[i]=refpix[i]-real(-x1[i]+1)
            x1[i]=1
        }
        x2[i]=xccd[i]-wx+sx
        if (x2[i] > 6218/l_xbin) {
            # deal correctly with binned data
            if (l_xbin==1)
                x2[i]=(6144+(2*chipgap))/l_xbin
            else
                x2[i]=(6144+(2*chipgap_bin))/l_xbin
        }
        nx[i]=x2[i]-x1[i]+1.
        #print(x1[i], x2[i], nx[i])
    }
    # find maximum length that does not overlap in 2-slit mode
    if (l_fl_novlap && nslit-j1==1) {
        #l_perovlap=0.05
        if (x2[1] > x1[2]) {
            novlap=x2[1]-x1[2]+1
            x2[1]=x1[2]-nint(0.5*l_perovlap*novlap)
            x1[2]=x2[1]-novlap+1+nint(0.5*l_perovlap*novlap)
            # JT: think the last line may be wrong but of no
            # practical consequence
            nx[1]=nx[1]-nint((1.+l_perovlap)*novlap)
            nx[2]=nx[2]-nint((1.+l_perovlap)*novlap)
            refpix[2]=refpix[2]-(1.+(0.5*l_perovlap))*novlap+1
        }
    }

    # JT, 2008: calculate remaining pixels from refpix to the end
    # (may be negative if the slit position is beyond the blue end
    # of the spectrum but that's OK); the blue CCD limit was already
    # factored in above when calculating x2, so rempix should never
    # fall off the detector
    for (i=1; i<=2; i+=1) {
        rempix[i]=nx[i]-refpix[i]
    }
    # Extract spectra of the same length from each slit
    refmin=min(refpix[1], refpix[2])
    remmin=min(rempix[1], rempix[2]) # JT: use instead of nxmax for blue limit
    if (l_fl_novlap) 
        nxmax=min(nx[1], nx[2])
    else 
        nxmax=max(nx[1], nx[2])

    # # If there's a problem with the bright ghost-type band in flats at the
    # # bottom of the detector being identified up as fibres by apall, we could
    # # cut it off like this, but for now apall seems to be doing the right
    # # thing with new MDFs so this is commented out rather than making
    # # unnecessary changes without checking all the grating offsets etc. (this
    # # would also need re-checking for any future GMOS-N Hamamatsu CCDs).
    # if (iccd==3) {
    #     # This is where the band tails off for an early B600 test dataset,
    #     # leaving a ~20 pixel buffer for flexure or grating differences.
    #     ymin = 110 / l_ybin
    # } else {
    #     ymin = 1
    # }
    ymin = 1
    ymax = l_naxis2

    if (nslit-j1==1) {
        # Data taken in 2-slit mode

        # JT: print the clean waveband corresponding to each slit, so
        # it is clear to the user what is available/truncated
        for (i=1; i<=2; i+=1) {
            slitl1 = l_wave-nmppx*(nint(xccd[i]+rempix[i])-xccd[i])
            slitl2 = l_wave-nmppx*(nint(xccd[i]-refpix[i]+1.)-xccd[i])
            printf("Slit %d clean waveband: %.1f - %.1f\n", i, slitl1, \
                slitl2) | scan(msg)
        printlog(msg, l_logfile, l_verbose)
        }

        # If we're actually extracting both slits then take the clean
        # region in common:
        if (slitnum==0) {
            # deal correctly with binned data
            # JT: use remmin instead of previous nxmax to solve a bug
            if (l_xbin==1)
                len=min(refmin+remmin, (refmin+((6144+(2*chipgap)) \
                    /l_xbin-xccd[2])))
            else
                len=min(refmin+remmin, (refmin+((6144+(2*chipgap_bin)) \
                    /l_xbin-xccd[2])))

            if (len < 1) {
                printlog("ERROR - GFEXTRACT: no common clean waveband: use \
                    exslits=red/blue or fl_novlap-", l_logfile, l_verbose)
                goto error
            }

            slitreg[1]="["//(nint(xccd[1]-refmin+1.))//":"// \ 
                (nint(xccd[1]-refmin+len))//","//ymin//":"//ymax//"]"
            slitreg[2]="["//(nint(xccd[2]-refmin+1.))//":"// \
                (nint(xccd[2]-refmin+len))//","//ymin//":"//ymax//"]"
        }
        # Otherwise, if only 1 of the 2 slits is being extracted, use its
        # own clean range and then change nslit and j1 to continue as if in
        # 1-slit mode:
        else {
            slitreg[1]="["//(nint(xccd[1]-refpix[1]+1.))//":"// \
                (nint(xccd[1]+rempix[1]))//","//ymin//":"//ymax//"]"
            slitreg[2]="["//(nint(xccd[2]-refpix[2]+1.))//":"// \ 
                (nint(xccd[2]+rempix[2]))//","//ymin//":"//ymax//"]"
            nslit=slitnum
            j1=slitnum
            refmin=refpix[j1]
        }
    } else {
        # 1-slit mode
        for (i=1; i<=2; i+=1) {
            slitreg[i]="["//x1[i]//":"//(x2[i])//","//ymin//":"//ymax//"]"
        }
        if (j1==1) {
            refmin=refpix[1]
        } else { 
            refmin=refpix[2]
        }
    }

    printlog("Slit section(s): "//slitreg[1]//" "//slitreg[2], l_logfile, \
        l_verbose)

    # Parse the slitreg and correct for binning to obtain DETSEC
    for (i=1; i<=2; i+=1) {
        print (slitreg[i]) | scanf ("[%d:%d,%d:%d]", detx1, detx2, dety1, \
            dety2)
        detx1 = ((detx1 - 1) * l_xbin) + 1
        detx2 = detx2 * l_xbin
        dety1 = ((dety1 - 1) * l_ybin) + 1
        dety2 = dety2 * l_ybin
        l_detsec[i] = "["//detx1//":"//detx2//","//dety1//":"//dety2//"]"
    }

    printlog("Detector sections(s): "//l_detsec[1]//" "//l_detsec[2], 
        l_logfile, l_verbose)

    printlog("", l_logfile, l_verbose)

    # Now start working on the image
    #
    # Work on a copy of the input image.  For example, GGAIN works on the
    # input image, we don't want to modify the input image as in some cases
    # it will bite.  Who would realize that running gfextract on a non-reduced
    # but prepared image would change the image gain spoiling it for future
    # reduction.
    
    copy (l_inimage, tmpin//".fits", verbose-)
    
    # Run GGAIN if not done already
    imgets(tmpin//"[0]", "GGAIN", >>& "dev$null")
    if (imgets.value=="0") {
        ggain (tmpin, gaindb=l_gaindb, logfile=l_logfile, 
            key_gain="GAIN", key_ron="RDNOISE", fl_mult+, fl_update+, 
            sci_ext=l_sci_ext, var_ext=l_var_ext, verbose=l_verbose)
        if (ggain.status != 0)
            goto error
    }

    # Sky subtract N&S spectra?
    if (l_fl_gnsskysub) {
        gnsskysub(tmpin, outimages=nsskysub, fl_fixnc=l_fl_fixnc, 
         logfile=l_logfile, verbose=l_verbose)
    # Don't need tmpin any more
        delete (tmpin//'.fits', verify-, >& 'dev$null')
        tmpin=nsskysub
        if (gnsskysub.status != 0)
            goto error
    } 

    # Paste the MEF together
    imgets(tmpin//"[0]", "GMOSAIC", >& "dev$null")
    if (imgets.value == "0" && l_outimage!="") {
        gmosaic (tmpin, outimages=impaste, outpref="", fl_paste=no, 
            fl_vardq=l_fl_vardq, fl_fixpix=l_fl_fixpix, fl_clean=no,
            geointer="linear", gap="default", bpmfile=l_bpmfile,
            statsec="default", obsmode="IMAGE", sci_ext=l_sci_ext,
            var_ext="VAR", dq_ext="DQ", mdf_ext="MDF", key_detsec="DETSEC",
            key_datsec="DATASEC", key_ccdsum="CCDSUM", key_obsmode="OBSMODE",
            logfile=l_logfile, fl_real-, fl_fulldq=l_fl_fulldq,
            dqthresh=l_dqthresh, verbose=l_verbose)
        if (gmosaic.status != 0)
            goto error
        delpaste=yes
    } else {
        impaste=tmpin
        delpaste=no
    }

    # Copy MDF to output image
    wmef(mdf//".fits", l_outimage, extnames="MDF", verbose-, phu=impaste, \
        >& "dev$null")

    if (l_title != "" && l_title != " ") {
        gemhedit(l_outimage//"[0]", "i_title", l_title, "", delete-)
    } else {
        imgets(tmpin//"[0]", "i_title", >>& "dev$null")
        if (imgets.value != "0")
            l_title=imgets.value
    }
    # No need for 'tmpin' anymore
    delete (tmpin, verify-, >& "dev$null")
    delete (tmpin//'.fits', verify-, >& 'dev$null')
    
    # Add APID column to MEF
    tcalc(l_outimage//"[MDF]", "APID", "0", datatype="int", colunits="", 
     colfmt="i4")

    # Loop over slits
    outlog="dev$null"
    if (l_inter) {
        outlog="STDOUT"
    }

    j=0
    for (i=j1; i<=nslit; i+=1) {
        # define temp files used within the loop
        sci=mktemp("tmpsci")
        var=mktemp("tmpvar")
        dq=mktemp("tmpdq")
        tmpmef1=mktemp("tmpmef1")
        tmpmef2=mktemp("tmpmef2")

        printlog("", l_logfile, l_verbose)
        printlog("Extracting slit "//i, l_logfile, l_verbose)
        j=j+1

        # Copy image section since apall won't accept section extensions
        imgsec=outroot//"_"//i
        varsec=outroot//"_var_"//i
        dqsec=outroot//"_dq_"//i
        if (imaccess(imgsec))
            imdelete(imgsec, verify-)
        if (imaccess(varsec))
            imdelete(varsec, verify-)
        if (imaccess(dqsec))
            imdelete(dqsec, verify-)

        imcopy(impaste//"["//l_sci_ext//",1]"//slitreg[i], imgsec, verbose-)
        gemhedit(imgsec, "DISPAXIS", 1, "", delete-)
        if (l_fl_vardq) {
            imcopy(impaste//"["//l_var_ext//",1]"//slitreg[i], varsec, verbose-)
            imcopy(impaste//"["//l_dq_ext//",1]"//slitreg[i], dqsec, verbose-)
            gemhedit(varsec, "DISPAXIS", 1, "", delete-)
            gemhedit(dqsec, "DISPAXIS", 1, "", delete-)
        }

        # Reference
        if (refroot != "") {
            l_reference=refroot//"_"//i
            if (!access("database/ap"//l_reference)) {
                printlog("ERROR - GFEXTRACT: Aperture reference \
                    database/ap"//l_reference//" not found.", l_logfile, \
                    verbose+)
                goto error
            }
        }

        # Make apidtable from MDF
        tselect (l_outimage//"[MDF]", mdfslit, 
            "NO>="//(750*(i-1)+1)//" && NO<="//(750*i)//" && BEAM >= 0")
        tprint (mdfslit, prparam-, prdata+, showrow-, showhdr-, showunits-, 
            columns="NO, BEAM, XINST, YINST, BLOCK", rows="-", option="plain", 
            align+, sp_col="", lgroup=0, > aptable)
        tinfo (aptable, ttout-)
        # Apall
        sr_logfile=specred.logfile
        specred.logfile=l_logfile
        # Make temporary output file
        tmptmpslit = mktemp("tmpslit")
        tmpslit = tmptmpslit//"_"//i//".ms"
        tmpvar = tmptmpslit//"_var_"//i//".ms"
        tmpdq = tmptmpslit//"_dq_"//i//".ms"

        # Adjust old hard defaults based on 5 pix fibre size for Hamamatsu:
        sfact = scale[inst,1] / scale[inst,iccd]
        a_lower = -2.5 * sfact
        a_upper = 2.5 * sfact
        a_width = 5. * sfact
        a_radius = 10. * sfact
        a_minsep = 4.0 * sfact
        a_maxsep = 200. * sfact

        apall (imgsec, tinfo.nrows, output=tmpslit, apertures="", \
            format="multispec", references=l_reference, profiles="", \
            inter=l_inter, find+, recenter=l_recenter, resize-, edit+, \
            trace=l_trace, fittrace+, extract+, extras-, review+, \
            line=l_line, nsum=l_nsum, lower=a_lower, upper=a_upper, \
            apidtable=aptable, width=a_width, radius=a_radius, \
            threshold=l_thresh, minsep=a_minsep, maxsep=a_maxsep, \
            order="increasing", aprecenter="", npeaks=INDEF, shift+, \
            llimit=INDEF, ulimit=INDEF, ylevel=0.1, peak+, bkg-, r_grow=0., \
            avglimits-, t_nsum=l_t_nsum, t_step=l_t_nsum, t_nlost=4, \
            t_function=l_function, t_order=l_order, t_sample="*", \
            t_naverage=1, t_niterate=3, t_low_reject=3., t_high_reject=3., \
            background="none", weights=l_weights, pfit="fit1d", clean-, \
            saturation=INDEF, gain=1., lsigma=4., usigma=4., nsubaps=1)

                # Repeat an identical extraction for the VAR/DQ:
        if (l_fl_vardq) {
            apall (varsec, tinfo.nrows, output=tmpvar, apertures="", \
                format="multispec", references=imgsec, profiles="", \
                inter-, find-, recenter-, resize-, edit-, \
                trace-, fittrace-, extract+, extras-, review-, \
                line=l_line, nsum=l_nsum, lower=a_lower, upper=a_upper, \
                apidtable=aptable, width=a_width, radius=a_radius, \
                threshold=l_thresh, minsep=a_minsep, maxsep=a_maxsep, \
                order="increasing", aprecenter="", npeaks=INDEF, shift+, \
                llimit=INDEF, ulimit=INDEF, ylevel=INDEF, peak+, bkg-, \
                r_grow=0., avglimits-, t_nsum=l_t_nsum, t_step=l_t_nsum, \
                t_nlost=4, t_function=l_function, t_order=l_order, \
                t_sample="*", t_naverage=1, t_niterate=3, t_low_reject=3., \
                t_high_reject=3., background="none", weights=l_weights, \
                pfit="fit1d", clean-, saturation=INDEF, gain=1., lsigma=4., \
                usigma=4., nsubaps=1)

            # Since apall has no option to do a bitwise AND rather than a
            # sum, we have to extract each DQ plane separately and recombine
            # them, to avoid ambiguity.

            # To get the IRAF database files named consistently, we need to
            # feed apall the original dqsec filename but containing each bit
            # plane in turn, so temporarily rename the original dqsec here.
            imrename (dqsec, tmpdqsec, verbose-)

            # Iterate over the range of bits actually used (this is more-
            # or-less cut & pasted from gmosaic but would be fiddly to spin
            # off into a separate script when the main step is different):
            imstat(tmpdqsec, fields="max", format=no) | scan (dqmax)
            if (dqmax < 65536) dqtype = "ushort"
            else dqtype = "uint"
            # Process DQ=0 to get a blank array if that's all there is:
            if (dqmax == 0) dqmin = 0
            else dqmin = 1
            for (n=dqmin; n <= dqmax; n*=2) {
                # Extract bit plane n:
                imexpr("a & b", dqsec, tmpdqsec, n,  outtype="real", verbose-)
                # Process the plane if it has non-zero values:
                imstat(dqsec, fields="npix", lower=1, upper=INDEF, \
                    format=no) | scan (ndq)
                if (ndq > 0 || n==0) {
                    apall (dqsec, tinfo.nrows, output=tmpplane, \
                        apertures="", format="multispec", references=imgsec, \
                        profiles="", inter-, find-, recenter-, resize-, \
                        edit-, trace-, fittrace-, extract+, extras-, review-, \
                        line=l_line, nsum=l_nsum, lower=a_lower, \
                        upper=a_upper, apidtable=aptable, width=a_width, \
                        radius=a_radius, threshold=0., minsep=a_minsep, \
                        maxsep=a_maxsep, order="increasing", aprecenter="", \
                        npeaks=INDEF, shift+, llimit=INDEF, ulimit=INDEF, \
                        ylevel=INDEF, peak+, bkg-, r_grow=0., avglimits-, \
                        t_nsum=l_t_nsum, t_step=l_t_nsum, t_nlost=4, \
                        t_function=l_function, t_order=l_order, t_sample="*", \
                        t_naverage=1, t_niterate=3, t_low_reject=3., \
                        t_high_reject=3., background="none", weights="none", \
                        pfit="fit1d", clean-, saturation=INDEF, gain=1., \
                        lsigma=4., usigma=4., nsubaps=1)
                    imdelete (dqsec, verify-, >& "dev$null")
                    # Round up affected pixels to integer DQ:
                    # (NB. outtype=auto gives 32 bits)
                    imexpr("abs(a)>c ? b : 0", tmpplane2, tmpplane//"[*,*,1]",
                        n, n*l_dqthresh, outtype=dqtype, verbose-)
                    imdelete (tmpplane, verify-, >& "dev$null")
                    # Add transformed plane to accumulated DQ:
                    if (imaccess(tmpaccum)) {
                        imexpr("a | b", tmpplane, tmpaccum, tmpplane2,
                            outtype=dqtype, verbose-)
                        imdelete (tmpplane2, verify-, >& "dev$null")
                        imdelete (tmpaccum, verify-, >& "dev$null")
                        imrename (tmpplane, tmpaccum, verbose-)
                    } else {
                        imrename (tmpplane2, tmpaccum, verbose-)
                    }
                } else {
                    imdelete (dqsec, verify-, >& "dev$null")
                }
                if (n==0) n=1  # otherwise loop gets stuck at 0

            } # end (loop over DQ planes)

            imrename (tmpaccum, tmpdq, verbose-)  # final result
            imrename (tmpdqsec, dqsec, verbose-)  # put input back

        } # end (if l_fl_vardq)

        specred.logfile=sr_logfile

        # check if output image exists, an interactive user may have quit or
        # apall may have crashed...
        if (!imaccess(tmpslit//".fits")) {
            printlog("ERROR - GFEXTRACT: apall output does not exist.",
                l_logfile, verbose+)
            goto error
        }

        # kluge the header
        gemhedit(tmpslit//".fits", "CDELT*", "", "", delete+)
        
        gemhedit(tmpslit//".fits", "CD1_1", 1.0, "", delete-)
        gemhedit(tmpslit//".fits", "CD2_2", 1.0, "", delete-)
        gemhedit(tmpslit//".fits", "CRPIX1", 1, "", delete-)
        gemhedit(tmpslit//".fits", "CRPIX2", 1, "", delete-)
        gemhedit(tmpslit//".fits", "CRVAL1", 1., "", delete-)
        gemhedit(tmpslit//".fits", "CRVAL2", 1., "", delete-)

        if (l_fl_vardq) {
            gemhedit(tmpvar//".fits", "CDELT*", "", "", delete+)

            gemhedit(tmpvar//".fits", "CD1_1", 1.0, "", delete-)
            gemhedit(tmpvar//".fits", "CD2_2", 1.0, "", delete-)
            gemhedit(tmpvar//".fits", "CRPIX1", 1, "", delete-)
            gemhedit(tmpvar//".fits", "CRPIX2", 1, "", delete-)
            gemhedit(tmpvar//".fits", "CRVAL1", 1., "", delete-)
            gemhedit(tmpvar//".fits", "CRVAL2", 1., "", delete-)

            gemhedit(tmpdq//".fits", "CDELT*", "", "", delete+)
        
            gemhedit(tmpdq//".fits", "CD1_1", 1.0, "", delete-)
            gemhedit(tmpdq//".fits", "CD2_2", 1.0, "", delete-)
            gemhedit(tmpdq//".fits", "CRPIX1", 1, "", delete-)
            gemhedit(tmpdq//".fits", "CRPIX2", 1, "", delete-)
            gemhedit(tmpdq//".fits", "CRVAL1", 1., "", delete-)
            gemhedit(tmpdq//".fits", "CRVAL2", 1., "", delete-)
        }

        printlog("Extraction done", l_logfile, l_verbose)

        # Separate the multispec file & optionally re-interpolate chip gaps
        # (and any other missing data) with fit1d if DQ is available:
        if (l_fl_fixgaps && l_fl_vardq) {
            printlog("Re-interpolate gaps", l_logfile, l_verbose)
            imcopy(l_outimage//"[0]", tmpmef1//".fits", verbose-)
            gemhedit(l_outimage//"[0]", "NSCIEXT", 1, "", delete-)
            imcopy(tmpslit//".fits[*,*,1]",
                tmpmef1//"["//l_sci_ext//",1,append]", verbose-)
            imcopy(tmpdq, tmpmef1//"["//l_dq_ext//",1,append]", verbose-)
            growdq(tmpmef1, tmpmef2, radius=l_grow, bitmask=16, logfile="",
                verbose+, >& "dev$null")
            imdelete(tmpmef1, verify-, >& "dev$null")
            gemfix(tmpmef2, tmpmef1, method="fit1d", grow=1.5, bitmask=16,
                axis=1, order=0, low_reject=3., high_reject=3., niter=5,
                fl_inter-, logfile="", verbose-)
            imdelete(tmpmef2, verify-, >& "dev$null")
            imcopy(tmpmef1//"["//l_sci_ext//",1]", sci, verbose-)
            imcopy(tmpvar//".fits[*,*,1]", var, verbose-)
            imcopy(tmpmef1//"["//l_dq_ext//",1]", dq, verbose-)
            imdelete(tmpmef1, verify-, >& "dev$null")
        } else {
            imcopy(tmpslit//".fits[*,*,1]", sci, verbose-)
            if (l_fl_vardq) {
                imcopy(tmpvar//".fits[*,*,1]", var, verbose-)
                imcopy(tmpdq, dq, verbose-)
            }
        }

        # Apply fiber throughput correction 
        if (l_response != "") {
            if (imaccess(l_response//"["//l_sci_ext//","//j//"]")) {
                # Compare the QE correction states of the two files
                if (in_qestate == resp_qestate) {
                    printlog("Applying fiber response", l_logfile, l_verbose)

                    imarith(sci, "/", l_response//"["//l_sci_ext//","//j//"]",\
                        sci, verbose-, noact-)
                } else {
                    printlog ("GFEXTRACT - ERROR: QE correction state "//\
                        "of the input file "//l_input//" ("//in_qestate//\
                        ") and the response file "//l_response//" ("//\
                        resp_qestate//") mismatch.", \
                        l_logfile, l_verbose)
                    goto error
                }
            } else {
                printlog("GFEXTRACT - ERROR: "//l_response//"["//l_sci_ext//",\
                    "//j//"] not found", l_logfile, verbose+)
                goto error
            }
        }
        
        # Update MEF with APID column
        tinfo(mdfslit, ttout-)
        apid=0
        for (k=1; k<=tinfo.nrows; k+=1) {
            tabpar(mdfslit, "NO", k)
            apno=int(tabpar.value)
            tcalc(l_outimage//"[MDF]", "APID", "if NO == "//apno//" then \
                "//k//" else APID")
        }

        # Write to output
        printlog("Writing to output file", l_logfile, l_verbose)
        imgets(l_outimage//"[0]", "NEXTEND")
        nextnd=int(imgets.value)
        if (l_fl_vardq) {
            fxinsert(sci//".fits,"//var//".fits,"//dq//".fits",
                l_outimage//"["//nextnd//"]", "0", verbose-, >& "dev$null")
            gemhedit(l_outimage//"[0]", "NEXTEND", (nextnd+3), "", delete-)
        } else {
            fxinsert(sci//".fits", l_outimage//"["//nextnd//"]", "0", \
                verbose-, >& "dev$null")
            gemhedit(l_outimage//"[0]", "NEXTEND", (nextnd+1), "", delete-)
        }
        gemhedit(l_outimage//"["//(nextnd+1)//"]", "EXTNAME", l_sci_ext, "", 
            delete-)
        gemhedit(l_outimage//"["//(nextnd+1)//"]", "EXTVER", j, "", delete-)
        gemhedit(l_outimage//"["//(nextnd+1)//"]", "i_title", l_title, "", 
            delete-)
        gemhedit(l_outimage//"["//(nextnd+1)//"]", "IFUSLIT", i, "", delete-)
        gemhedit(l_outimage//"["//(nextnd+1)//"]", "REFPIX1", refmin, "", 
            delete-)
        gemhedit(l_outimage//"["//(nextnd+1)//"]", "DETSEC", l_detsec[j], \
            "Detector section(s)",delete-)

        if (l_fl_vardq) {
            gemhedit(l_outimage//"["//(nextnd+2)//"]", "EXTNAME", l_var_ext, \
                "", delete-)
            gemhedit(l_outimage//"["//(nextnd+2)//"]", "EXTVER", j, "", 
                delete-)
            gemhedit(l_outimage//"["//(nextnd+2)//"]", "i_title", 
                "Variance: "//l_title, "", delete-)
            gemhedit(l_outimage//"["//(nextnd+2)//"]", "IFUSLIT", i, "", 
                delete-)
            gemhedit(l_outimage//"["//(nextnd+2)//"]", "REFPIX1", refmin, "", 
                delete-)
            gemhedit(l_outimage//"["//(nextnd+3)//"]", "EXTNAME", l_dq_ext, "",
                delete-)
            gemhedit(l_outimage//"["//(nextnd+3)//"]", "EXTVER", j, "", 
                delete-)
            gemhedit(l_outimage//"["//(nextnd+3)//"]", "i_title", 
                "Data Quality: "//l_title, "", delete-)
            gemhedit(l_outimage//"["//(nextnd+3)//"]", "IFUSLIT", i, "", 
                delete-)
            gemhedit(l_outimage//"["//(nextnd+3)//"]", "REFPIX1", refmin, "", 
                delete-)
            gemhedit(l_outimage//"["//(nextnd+3)//"]", "DETSEC", l_detsec[j], \
                "Detector section(s)",delete-)
        }
        imdelete(imgsec//","//sci//","//var//","//dq, verify-, >& "dev$null")
        if (l_fl_vardq)
            imdelete(varsec//","//dqsec, verify-, >& "dev$null")
        imdelete(tmpslit, verify-, >& "dev$null")
#        imdelete("tmpslit*.ms.fits", verify-, >& "dev$null")
        delete(aptable//","//mdfslit, verify-, >& "dev$null")
    }

    # final header update
    gemhedit(l_outimage//"[0]", "NSCIEXT", j, "", delete-)
    gemdate ()
    gemhedit (l_outimage//"[0]", "GFEXTRAC", gemdate.outdate, 
        "UT Time stamp for GFEXTRACT", delete-)
    gemhedit (l_outimage//"[0]", "GEM-TLM", gemdate.outdate, 
        "UT Last modification with GEMINI", delete-)

    # gsappwave
    if (l_fl_gsappwave) {
        printlog("", l_logfile, l_verbose)
        gsappwave(l_outimage, logfile=l_logfile, gratingdb=l_gratingdb, 
            filterdb=l_filterdb, key_dispaxis="DISPAXIS", dispaxis=1, 
            verbose=l_verbose)
    }

    # clean up
    goto clean

error:
    status=1
    goto clean

clean:
    delete("tmpslit*.fits,"//mdf//".fits"//","//aptable, verify-, \
        >& "dev$null")
    imdelete(sci//","//nsskysub, verify-, >& "dev$null")
    if (l_fl_vardq)
        imdelete(var//","//dq, verify-, >& "dev$null")
    delete(aptable//","//mdfslit, verify-, >& "dev$null")
    delete (tmpin, verify-, >& "dev$null")
    delete (tmpin//'.fits', verify-, >& 'dev$null')
    if (delpaste) 
        imdelete(impaste, verify-, >& "dev$null")
    # close log file
    printlog(" ", l_logfile, l_verbose)
    if (status==0)
        printlog("GFEXTRACT exit status: good", l_logfile, l_verbose)
    else
        printlog("GFEXTRACT exit status: error", l_logfile, l_verbose)
    printlog("----------------------------------------------------------\
        ----------------------", l_logfile, l_verbose )

end

