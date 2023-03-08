# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gfreduce (inimages)

# Reduce GMOS IFU images
#
# Version   Sept 20, 2002 BM v1.4 release
#           Dec 4, 2003   BM generalize for GMOS-S IFU
#           Feb 29, 2004  BM added thresh param and gnsskysub option
#           Mar 19, 2004  BM added fl_novlap and perovlap params
#           Mar 24, 2004  BM fix bug in call to gscalibrate
#           Mar 25, 2004  BM Now exit gracefully if problem with input images
#           Apr 07, 2004  BM Fix print of INDEF, auto-select N&S MDF and
#                            set l_fl_skysub if N&S mask used
#           Oct 01, 2008  JT Automatically set MDF for new data

string  inimages        {prompt="Input images"}
string  outimages       {"",prompt="Output images"}
string  outpref         {"default",prompt="Prefix for output images"}
string  slits           {"header",enum="red|blue|both|header",prompt="Slit(s) used"}
string  exslits         {"*",enum="red|blue|*",prompt="Slit(s) to extract"}
bool    fl_nodshuffle   {no,prompt="Nod & Shuffle output slit mask used"}
bool    fl_inter        {yes,prompt="Interactive use?"}
bool    fl_vardq        {no,prompt="Use variance and DQ planes?"}
bool    fl_addmdf       {yes,prompt="Add mask definition file (run GPREPARE)"}
bool    fl_over         {yes,prompt="Apply overscan correction?"}
bool    fl_trim         {yes,prompt="Trim off the overscan section?"}
bool    fl_bias         {yes,prompt="Apply zero level (bias) correction?"}
bool    fl_scatsub      {no,prompt="Remove scattered light with gfscatsub?"}
bool    fl_qecorr       {no, prompt="QE correct the input images?"}
bool    fl_gscrrej      {yes,prompt="Clean cosmic rays (old GSCRREJ method)?"}
bool    fl_crspec       {no,prompt="Clean cosmic rays (new LA Cosmic method)"}
bool    fl_gnsskysub    {no,prompt="Run GNSSKYSUB?"}
bool    fl_extract      {yes,prompt="Extract spectra?"}
bool    fl_gsappwave    {yes,prompt="Run gsappwave?"}
bool    fl_wavtran      {yes,prompt="Transform cube to apply wavelength solution?"}
bool    fl_skysub       {yes,prompt="Subtract sky?"}
bool    fl_fluxcal      {yes,prompt="Apply flux calibration?"}
bool    fl_fulldq       {no, prompt="Decompose DQ during gmosaic transformation; apply correct handling?"}            # OLDP-3
real    dqthresh        {0.1, min=0.0, max=0.5, prompt="Threshold applied to DQ when fl_fulldq=yes"}
string  rawpath         {"",prompt="Path for input raw images"}
string  key_mdf         {"",prompt="Header keyword for MDF"}
string  mdffile         {"default",prompt="MDF file to use if keyword not found"}
string  mdfdir          {"gmos$data/",prompt="MDF database directory"}
string  key_biassec     {"BIASSEC",prompt="Header keyword for bias section"}
string  key_datasec     {"DATASEC",prompt="Header keyword for data section"}
string  bpmfile         {"gmos$data/chipgaps.dat",prompt="Bad pixel mask for column interpolation"}
real    grow            {1.0, min=0, prompt="Gap growth radius in pixels"}
string  bias            {"",prompt="Bias image name"}
string  reference       {"",prompt="Reference file"}
string  sc_xorder       {"1",prompt="X order for surface fit for gfscatsub"}
string  sc_yorder       {"3",prompt="Y order for surface fit for gfscatsub"}
bool    sc_cross        {no,prompt="Include cross terms in gfscatsub"}
char    qe_refim        {"", prompt="QE wavelength reference image (spectroscopy only)"}
bool    fl_keep_qeim    {no, prompt="Keep QE correction (spectroscopy only)?"}
char    qe_corrpref     {"qecorr", prompt="Prefix for QE correction files (spectroscopy only)"}
char    qe_corrimages   {"", prompt="Name for QE correction data (spectroscopy only)"}
char    qe_data         {"gmosQEfactors.dat", prompt="Data file that contains QE information."}
char    qe_datadir      {"gmos$data/", prompt="Directory containg QE data file."}
string  response        {"",prompt="Fiber response file"}
string  wavtraname      {"",prompt="Names of wavelength calibrations"}
string  sfunction       {"",prompt="Sensitivity function from GSSTANDARD"}
string  extinction      {"",prompt="Extinction file"}
bool    fl_fixnc        {no,prompt="Auto-correct for nod count mismatch?"}
bool    fl_fixgaps      {no,prompt="Re-interpolate chip gaps after extraction?"}
bool    fl_novlap       {yes,prompt="Avoid spectral overlap?"}
real    perovlap        {10.,prompt="Percentage by which to shrink overlapping spectra"}
string  nbiascontam     {"default", prompt="Number of columns removed from overscan region"}
string  biasrows    {"default", prompt="Rows to use for overscan region"}
char    order           {"default", prompt="Order of overscan fitting function"}
real    low_reject      {3., prompt="Low sigma rejection factor in overscan fit"}
real    high_reject     {3., prompt="High sigma rejection factor in overscan fit"}
int     niterate        {2, prompt="Number of rejection iterations in overscan fit"}
int     cr_xorder       {9, prompt="Order of GEMCRSPEC object fit (0=no fit)"}
real    cr_sigclip      {4.5, prompt="Detection limit for GEMCRSPEC (sigma)"}
real    cr_sigfrac      {0.5, prompt="Fractional GEMCRSPEC neighbour detection limit"}
real    cr_objlim       {1.0, prompt="Contrast limit between CR and underlying object"}
int     cr_niter        {4, prompt="Maximum number of GEMCRSPEC iterations"}
int     line            {INDEF,prompt="Line for finding peaks"}
int     nsum            {10,prompt="Number of columns to use for finding apertures"}
bool    trace           {yes,prompt="Trace spectra"}
bool    recenter        {yes,prompt="Recenter apertures?"}
real    thresh          {200.,prompt="Detection threshold for profile centering"}
string  function        {"chebyshev",enum="chebyshev|spline1|spline3|legendre",prompt="Function for trace"}
int     t_order         {5,min=1,prompt="Order of trace fit"}
int     t_nsum          {10,min=1,prompt="Number of dispersion lines to be summed during trace"}
string  weights         {"variance",enum="variance|none",prompt="Weighting during extraction"}
string  gratingdb       {"gmos$data/GMOSgratings.dat",prompt="Gratings database file"}
string  filterdb        {"gmos$data/GMOSfilters.dat",prompt="Filters database file"}
real    xoffset         {INDEF,prompt="X offset in wavelength [nm]"}
string  expr            {"default",prompt="Expression for selection of sky spectra"}
bool    sepslits        {no,prompt="Subtract sky from each slit separately?"}
real    w1              {INDEF,prompt="Starting wavelength"}
real    w2              {INDEF,prompt="Ending wavelength"}
real    dw              {INDEF,prompt="Wavelength interval per pixel"}
int     nw              {INDEF,prompt="Number of output pixels"}
string  observatory     {"default",prompt="Observatory name"}
string  sci_ext         {"SCI",prompt="Name of science extension"}
string  var_ext         {"VAR",prompt="Name of variance extension"}
string  dq_ext          {"DQ",prompt="Name of data quality extension"}
string  logfile         {"",prompt="Logfile name"}
bool    verbose         {yes,prompt="Verbose?"}
int     status          {0,prompt="Exit status (0=good)"}
struct  *scanfile       {"",prompt="Internal use only"}

begin

    # Local variables for input parameters
    string  l_inimages, l_outimages, l_prefix, l_logfile, l_rawpath
    string  l_key_mdf, l_mdffile, l_mdfdir, l_slits, l_exslits
    string  l_bias, l_key_biassec, l_key_datasec
    string  l_response, l_gratingdb, l_filterdb
    string  l_reference, l_bpmfile, l_function, l_weights
    string  l_wavtraname, l_expr, l_nbiascontam
    string  l_sci_ext, l_var_ext, l_dq_ext, l_sc_xorder, l_sc_yorder
    string  l_obs, l_sfunc, l_extinct, l_biasrows, l_order
    real    l_xoffset, l_thresh, l_perovlap, l_low_reject, l_high_reject
    real    l_dqthresh, l_w1, l_w2, l_dw, l_grow
    real    l_cr_sigclip, l_cr_sigfrac, l_cr_objlim
    int     l_t_order, l_line, l_nsum, l_t_nsum, l_nw
    int     l_niterate, atpos, l_cr_xorder, l_cr_niter
    bool    l_verbose, l_fl_inter, l_fl_vardq, l_sepslits
    bool    l_fl_addmdf, l_fl_over, l_fl_bias, l_fl_trim
    bool    l_fl_extract, l_recenter,l_trace, l_fl_gsappwave, l_fl_fulldq
    bool    l_fl_wavtran, l_fl_skysub, l_fl_fluxcal, l_fl_ext, l_fl_gscrrej
    bool    l_fl_crspec, l_fl_gnsskysub, l_fl_fixnc, l_fl_fixgaps, l_fl_novlap
    bool    l_fl_nodshuffle, l_sc_cross

    # Other local variables
    string  inlist, infile[500], outlist, temp1, temp2
    string  img, outimg, suf, root, maskname, mdfprefix, refroot, gapfile
    string  ifu_mdf_stub, slit_name, dettype
    int     i, j, slitnum, nbad, nin, nout, nx, nslit, j1, inst, iccd
    bool    useprefix, usenodsh, havegapfile
    struct  sdate

    # QE local variables
    char    l_qe_refim, l_qe_corrpref, l_qecorr_data
    char    l_qecorr_datadir, l_qe_corrimages
    bool    l_fl_qecorr, l_fl_keep_qeim, l_fl_scatsub

    # Initialize exit status
    status = 0

    # cache some parameter files
    cache ("imgets", "gimverify", "fparse", "gemlogname","gemdate")

    # Initialize local variables
    l_inimages=inimages ; l_outimages=outimages ; l_prefix=outpref
    l_logfile=logfile ; l_verbose=verbose ; l_fl_inter=fl_inter
    l_rawpath=rawpath ; l_key_biassec=key_biassec
    l_fl_addmdf=fl_addmdf ; l_key_mdf=key_mdf ; l_slits=slits
    l_exslits=exslits; l_mdffile=mdffile ; l_mdfdir=mdfdir
    l_fl_over=fl_over ; l_fl_bias=fl_bias ; l_order=order
    l_low_reject=low_reject; l_high_reject=high_reject; l_niterate=niterate
    l_cr_xorder=cr_xorder ; l_cr_sigclip=cr_sigclip ; l_cr_sigfrac=cr_sigfrac
    l_cr_objlim=cr_objlim ; l_cr_niter=cr_niter
    l_bias=bias; l_fl_trim=fl_trim ; l_key_datasec=key_datasec
    l_response=response ; l_reference=reference ; l_bpmfile=bpmfile
    l_grow=grow ; l_sci_ext=sci_ext ; l_var_ext=var_ext ; l_dq_ext=dq_ext
    l_line=line ; l_nbiascontam=nbiascontam; l_biasrows=biasrows
    l_fl_vardq=fl_vardq
    l_gratingdb=gratingdb
    l_filterdb=filterdb
    l_trace=trace
    l_nsum=nsum
    l_t_nsum=t_nsum
    l_t_order=t_order
    l_function=function
    l_weights=weights
    l_fl_gsappwave=fl_gsappwave
    l_recenter=recenter ; l_thresh=thresh
    l_fl_extract=fl_extract
    l_fl_wavtran=fl_wavtran ; l_wavtraname=wavtraname
    l_fl_skysub=fl_skysub ; l_expr=expr ; l_sepslits = sepslits
    l_fl_fluxcal=fl_fluxcal ; l_sfunc=sfunction ; l_extinct=extinction
    l_obs=observatory
    l_fl_gscrrej=fl_gscrrej ; l_fl_crspec=fl_crspec
    l_xoffset=xoffset
    l_w1=w1; l_w2=w2; l_dw=dw; l_nw=nw
    l_fl_gnsskysub=fl_gnsskysub ; l_fl_fixnc=fl_fixnc ; l_fl_fixgaps=fl_fixgaps
    l_fl_novlap=fl_novlap ; l_perovlap=perovlap
    l_fl_nodshuffle=fl_nodshuffle ; l_fl_scatsub=fl_scatsub
    l_sc_xorder=sc_xorder ; l_sc_yorder=sc_yorder ; l_sc_cross=sc_cross
    # GQECORR Parameters:
    l_fl_qecorr=fl_qecorr
    l_qe_refim = qe_refim
    l_fl_keep_qeim = fl_keep_qeim
    l_qe_corrpref = qe_corrpref
    l_qecorr_data = qe_data
    l_qecorr_datadir = qe_datadir
    l_qe_corrimages = qe_corrimages
    # DQ handling
    l_fl_fulldq = fl_fulldq
    l_dqthresh = dqthresh

    # Test the logfile:
    gemlogname (logpar=l_logfile, package="gmos")
    if (gemlogname.status != 0) {
        goto error
    }
    l_logfile = gemlogname.logname

    printlog ("-------------------------------------------------------------\
        ------------------", l_logfile, l_verbose)
    date | scan(sdate)
    printlog ("GFREDUCE -- "//sdate, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    # Moved these initiations here due to it failng when fl_qecorr=yes when
    # QE correction is switched off - MS
    temp1 = mktemp("tmpin")
    temp2 = mktemp("tmpout")

    #The usual kind of checks to make sure everything we need is specified.

    nbad = 0

    if (l_inimages=="" || l_inimages==" ") {
        printlog ("ERROR - GFREDUCE: Input spectra is an empty string",
            l_logfile, yes)
        nbad = nbad+1
    }

    atpos = strstr("@", l_inimages)
    if (atpos >= 1) {
        inlist = substr(l_inimages,atpos + 1, strlen(l_inimages))
        if (!access(inlist)) {
            printlog ("ERROR - GFREDUCE: Input list "//inlist//" not found",
                l_logfile, verbose+)
            nbad = nbad+1
        }
    }

    useprefix = yes
    if (l_outimages!="" && l_outimages!=" ") {
        useprefix = no
    } else if (l_prefix=="default") {
        printlog ("Default prefixes will be used.", l_logfile, l_verbose)
    } else if (l_prefix == "" || l_prefix == " ") {
        printlog ("ERROR - GFREDUCE: Neither outimages or outpref given.",
            l_logfile, yes)
        nbad = nbad+1
    }

    atpos = strstr("@", l_outimages)
    if (atpos >= 1) {
        outlist = substr(l_outimages,atpos + 1, strlen(l_outimages))
        if (!access(outlist)) {
            printlog ("ERROR - GFREDUCE: Output list "//outlist//" not found",
                l_logfile, yes)
            nbad = nbad+1
        }
    }

    #    if (l_fl_gnsskysub && l_fl_skysub) {
    #        printlog("ERROR - GFREDUCE: cannot use both fl_gnsskysub and \
    #            fl_skysub options.",l_logfile,yes)
    #        nbad=nbad+1
    #    }

    sections (l_inimages, > temp1)
    count(temp1) | scan(nin)
    if (!useprefix) {
        sections (l_outimages, > temp2)
        count(temp2) | scan(nout)
        if (nin != nout) {
            printlog ("ERROR - GFREDUCE: Different number of input and \
                output spectra", l_logfile, yes)
            nbad = nbad+1
        }
    } else if (l_prefix != "default") {
        files (l_prefix//"//@"//temp1, sort-, > temp2)
    }

    i = 0
    scanfile = temp1
    while (fscan(scanfile,img) !=EOF) {
        gimverify (l_rawpath//img)
        if (gimverify.status>0) {
            printlog ("ERROR - GFREDUCE: "//img//" not found or not a MEF",
                l_logfile, yes)
            nbad = nbad+1
        } else {
            i = i+1
            # name w/o suffix
            fparse (gimverify.outname, verbose-)
            infile[i] = fparse.root//".fits"
        }
    } #end of while loop over images.
    scanfile = ""

    # Check existence of output images
    if (access(temp2)) {
        scanfile = temp2
        while (fscan(scanfile,img) != EOF) {
            if (imaccess(img//"[0]")) {
                printlog ("ERROR - GFREDUCE: "//img//" already exists",
                    l_logfile, yes)
                nbad = nbad+1
            }
        }
        scanfile = ""
    }

    # If subtracting scattered light, require a reference file or a gaps file
    # derived from it previously (with fl_scatsub-).
    if (l_fl_scatsub) {

        fparse(l_reference, verbose-)
        refroot = fparse.root
        gapfile = refroot//"_gaps"

        if (refroot == "" || !access(gapfile)) {
            havegapfile = no
            gimverify(l_reference)
            if (refroot == "" || gimverify.status != 0) {
                printlog ("ERROR - GFREDUCE: missing reference image required \
                           by fl_scatsub+",
                    l_logfile, yes)
                nbad = nbad+1
            }
        } else {
            havegapfile = yes
	}
    }

    #If anything was wrong then exit.
    if (nbad > 0) {
        printlog ("ERROR - GFREDUCE: "//nbad//" error(s) found with input \
            parameters. Exiting.", l_logfile, yes)
        goto error
    }

    # Determine which instrument
    imgets (l_rawpath//infile[1]//"[0]", "INSTRUME", >>& "dev$null")
    if (imgets.value == "0") {
        printlog ("ERROR - GFREDUCE: Instrument keyword not found.",
            l_logfile, verbose+)
        goto error
    }
    inst = 1 # Default is GMOS-N, support for old data
    if (imgets.value == "GMOS-S") {
        inst = 2
    }

    # Determine which detector type is being used
    iccd = 1
    keypar (l_rawpath//infile[1]//"[0]", "DETTYPE", silent+)
    if (keypar.found) {
        dettype = keypar.value
        if (dettype == "SDSU II CCD") {
            # CCDs GMOS-N and GMOS-S on 2011-07-24
            iccd = 1
        } else if (dettype == "SDSU II e2v DD CCD42-90") {
            # New e2v CCDs GMOS-N
            iccd = 2
        } else if (dettype == "S10892" || dettype == "S10892-N") {
            # GMOS-N and GMOS-S Hamamatsu CCDs
            iccd = 3
        } else {
            printlog ("ERROR - GFREDUCE: DETTYPE not recognised for "//\
                img, l_logfile, verbose+)
            goto error
        }
    }

    # Determine if using nod&shuffle
    usenodsh = no
    imgets (l_rawpath//infile[1]//"[0]", "NODPIX", >& "dev$null")
    if (imgets.value != "0") {
        usenodsh = yes
    }

    if (usenodsh || l_fl_nodshuffle) {
        if (l_fl_skysub) {
            printlog ("WARNING - GFREDUCE: Nod&Shuffle mode on.  Resetting \
                fl_skysub to 'no'", l_logfile, verbose+)
        }
        l_fl_skysub = no
    } else {
        if (l_fl_gnsskysub) {
            printlog ("WARNING - GFREDUCE: Nod&Shuffle mode off.  Resetting \
                fl_gnsskysub to 'no'", l_logfile, verbose+)
        }
        l_fl_gnsskysub = no
    }

    # Get mask name
    imgets (l_rawpath//infile[1]//"[0]", "MASKNAME", >>& "dev$null")
    maskname = imgets.value
    if (maskname == "0") {
        printlog ("ERROR - GFREDUCE: Maskname keyword not found.",
            l_logfile, verbose+)
        goto error
    }

    # Default MDF files
    if (l_mdffile=="default") {
        # Use detector type to determine stub of MDF name
        ifu_mdf_stub = "_mdf"
        if (iccd == 3) {
            if (inst == 2) {
                ifu_mdf_stub = ifu_mdf_stub//"_HAM"
            }
        }

        # Set prefix based on instrument
        if (inst==1) {
            mdfprefix="gnifu_"
        } else {
            mdfprefix="gsifu_"
        }

        # Default to 1 slit
        nslit=1

        # Set the maskname
        if (l_slits=="header") {
            if (maskname == "IFU") {
                printlog ("ERROR - GFREDUCE: Must specify slits=red/blue/both \
                    mdffile for old data", l_logfile, l_verbose)
                goto error
            } else {
                if (maskname == "IFU-2" || \
                    strstr(mdfprefix//"slits"//ifu_mdf_stub, maskname) > 0) {
                    slit_name = "slits"
                    nslit=2
                } else if (maskname == "IFU-B" || \
                    strstr(mdfprefix//"slitb"//ifu_mdf_stub, maskname) > 0) {
                    slit_name = "slitb"
                } else if (maskname == "IFU-R" || \
                    strstr(mdfprefix//"slitr"//ifu_mdf_stub, maskname) > 0) {
                    slit_name = "slitr"
                } else if (maskname == "IFU-NS-2" || \
                    strstr(mdfprefix//"ns_slits"//ifu_mdf_stub, maskname) > 0) {
                    slit_name = "ns_slits"
                    nslit=2
                } else if (maskname == "IFU-NS-B" || \
                    strstr(mdfprefix//"ns_slitb"//ifu_mdf_stub, maskname) > 0) {
                    slit_name = "ns_slitb"
                } else if (maskname == "IFU-NS-R" || \
                    strstr(mdfprefix//"ns_slitr"//ifu_mdf_stub, maskname) > 0) {
                    slit_name = "ns_slitr"
                }
            } # End if (old/new data with non-unique/meaningful MASKNAME)
        } else {
            if (l_slits=="red") {
                slit_name = "slitr"
            } else if (l_slits=="blue") {
                slit_name = "slitb"
            } else {
                slit_name = "slits"
                nslit=2
            }
            if (inst == 2 && (usenodsh || l_fl_nodshuffle)) {
                slit_name = "ns_"//slit_name
            } # End if (nod & shuffle)
        } # End if (MDF from header or GMOS-N or GMOS-S)

        # Set the full name
        l_mdffile = mdfprefix//slit_name//ifu_mdf_stub//".fits"
    } else if (l_slits=="both") {
        nslit=2
    } else {
        nslit=1
    } # End if (use default MDF rather than specified file)

    # For 1-of-2 slits extraction, adjust nslit accordingly, like
    # in gfextract (to allow for 2-slit IFU data with an abnormal
    # wavelength/filter setup)
    j1=1
    if (nslit==2) {
        if (l_exslits=="red") {
            nslit=1
            j1=1
        } else if (l_exslits=="blue") {
            nslit=2
            j1=2
        }
    }

    # Default sky selection region
    if (l_expr=="default") {
        if (inst==1) {
            l_expr = "XINST > 10."
        } else {
            l_expr = "XINST < 10."
        }
    }

    # Default observatory
    if (l_obs=="default") {
        if (inst==1) {
            l_obs = "Gemini-North"
        } else {
            l_obs = "Gemini-South"
        }
    }

    #If we are here then everything should be OK.
    #Write all the relevant info to the logfile:
    #
    printlog ("", l_logfile, l_verbose)
    printlog ("inimages   = "//l_inimages, l_logfile, l_verbose)
    printlog ("outimages  = "//l_outimages, l_logfile, l_verbose)
    printlog ("outpref    = "//l_prefix, l_logfile, l_verbose)
    printlog ("slits      = "//l_slits, l_logfile, l_verbose)
    if (l_fl_nodshuffle || usenodsh) {
        printlog ("nod&shuffle= yes", l_logfile, l_verbose)
    } else {
        printlog ("nod&shuffle= no", l_logfile, l_verbose)
    }
    printlog ("fl_scatsub    = "//l_fl_scatsub, l_logfile, l_verbose)
    printlog ("fl_qecorr     = "//l_fl_qecorr, l_logfile, l_verbose)
    printlog ("mdffile       = "//l_mdffile, l_logfile, l_verbose)
    printlog ("mdfdir        = "//l_mdfdir, l_logfile, l_verbose)
    printlog ("bias          = "//l_bias, l_logfile, l_verbose)
    printlog ("reference     = "//l_reference, l_logfile, l_verbose)
    printlog ("qe_refim      = "//l_qe_refim, l_logfile, l_verbose)
    printlog ("qe_corrimages = "//l_qe_corrimages, l_logfile, l_verbose)
    printlog ("response      = "//l_response, l_logfile, l_verbose)
    printlog ("wavtraname    = "//l_wavtraname, l_logfile, l_verbose)
    printlog ("sfunction     = "//l_sfunc, l_logfile, l_verbose)
    printlog ("extinction    = "//l_extinct, l_logfile, l_verbose)
    printlog ("expr          = "//l_expr, l_logfile, l_verbose)
    printlog ("gratingdb     = "//l_gratingdb, l_logfile, l_verbose)
    printlog ("filterdb      = "//l_filterdb, l_logfile, l_verbose)
    if (l_xoffset != INDEF) {
        printlog ("xoffset    = "//str(l_xoffset), l_logfile, \
            verbose=l_verbose)
    } else {
        printlog ("xoffset    = INDEF", l_logfile, verbose=l_verbose)
    }

    printlog (" ", l_logfile, l_verbose)

    # gprepare
    if (l_fl_addmdf) {
        gprepare ("@"//temp1, rawpath=l_rawpath, outimages="", outpre="g",
            fl_addmdf=l_fl_addmdf, key_mdf=l_key_mdf, mdffile=l_mdffile,
            mdfdir=l_mdfdir, verbose=l_verbose, logfile=l_logfile,
            sci_ext=l_sci_ext)
        if (gprepare.status != 0)
            goto error
        for (i=1; i<=nin; i+=1) {
            infile[i] = "g"//infile[i]
        }
    }

    # gireduce
    if (l_fl_over || l_fl_trim || l_fl_bias) {
        delete (temp1, verify-, >& "dev$null")
        for (i=1; i<=nin; i+=1) {
            print (infile[i], >> temp1)
        }
        gireduce ("@"//temp1, outimage="", outpre="r", fl_over=l_fl_over,
            fl_trim=fl_trim, fl_bias=fl_bias, fl_dark-, fl_qecorr-,
            fl_flat-, order=l_order, low_reject=l_low_reject,
            high_reject=l_high_reject,
            niterate=l_niterate, key_biassec=l_key_biassec,
            key_datasec=l_key_datasec, bias=l_bias, verbose=l_verbose,
            fl_vardq=l_fl_vardq, fl_inter=l_fl_inter, logfile=l_logfile,
            sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext,
            nbiascontam=l_nbiascontam, biasrows=l_biasrows, fl_mult=yes)
        if (gireduce.status != 0)
            goto error
        for (i=1; i<=nin; i+=1) {
            infile[i] = "r"//infile[i]
        }
    }

    # gfscatsub (scattered light)
    if (l_fl_scatsub) {

        # Create a file listing gap locations, based on the apextract limits,
        # if not done previously (and perhaps tweaked manually since).
        if (havegapfile) {
            printlog ("Re-using existing gaps file: "//gapfile,
                l_logfile, l_verbose)
        } else if (nin > 0) {
            printlog ("Running gffindblocks:", l_logfile, l_verbose)
            gffindblocks (infile[1], l_reference, gapfile)
            if (gffindblocks.status != 0)
                goto error
        }
        # Using the gap locations from above, model & subtract the scattered
        # light from each input image:
        for (i=1; i<=nin; i+=1) {
            printlog ("Running gfscatsub", l_logfile, l_verbose)
            gfscatsub (infile[i], mask=gapfile, outimage="b"//infile[i],
                xorder=l_sc_xorder, yorder=l_sc_yorder, cross=l_sc_cross)
            if (gfscatsub.status != 0)
                goto error
            infile[i] = "b"//infile[i]
        }
    }

    # gemcrspec (new GMOS CR cleaning method):
    if (l_fl_crspec) {
        delete (temp1, verify-, >& "dev$null")
        for (i=1; i<=nin; i+=1) {
            print (infile[i], >> temp1)
            infile[i] = "x"//infile[i]
        }
        gemcrspec ("@"//temp1, "x@"//temp1, xorder=l_cr_xorder,
            yorder=-1, sigclip=l_cr_sigclip, sigfrac=l_cr_sigfrac,
            objlim=l_cr_objlim, niter=l_cr_niter, fl_vardq=l_fl_vardq,
            sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext,
            # These aren't set in gireduce either; keep changes minimal for now
            # key_ron=l_key_ron, key_gain=l_key_gain, ron=l_ron, gain=l_gain,
            logfile=l_logfile, verbose=l_verbose)
        if (gemcrspec.status != 0)
            goto error
    }

    # gscrrej (original GMOS CR cleaning method):
    if (l_fl_gscrrej) {
        delete (temp1, verify-, >& "dev$null")
        for (i=1; i<=nin; i+=1) {
            gscrrej (infile[i], "x"//infile[i], verbose=l_verbose,
                fl_inter=l_fl_inter, logfile=l_logfile)
            if (gscrrej.status != 0)
                goto error
            infile[i] = "x"//infile[i]
        }
    }

    # If applying QE correction, re-run gireduce to do it (because the latter
    # includes countless checks that we don't want to repeat here). This is
    # best done after gfscatsub so as not to introduce discontinuities in the
    # background, hence why we don't do it in the first gireduce call.
    if (l_fl_qecorr) {
        delete (temp1, verify-, >& "dev$null")
        for (i=1; i<=nin; i+=1) {
            print (infile[i], >> temp1)
        }
        gireduce ("@"//temp1, outimage="", outpre="q", fl_over-,
            fl_trim-, fl_bias-, fl_dark-, fl_qecorr+, fl_flat-, fl_addmdf-,
            qe_refim=l_qe_refim, fl_keep_qeim=l_fl_keep_qeim,
            qe_corrpref=l_qe_corrpref, qe_corrimage=l_qe_corrimages,
            qe_data=l_qecorr_data, qe_datadir=l_qecorr_datadir,
            key_biassec=l_key_biassec, key_datasec=l_key_datasec,
            verbose=l_verbose, fl_vardq=l_fl_vardq, fl_inter=l_fl_inter,
            logfile=l_logfile, sci_ext=l_sci_ext, var_ext=l_var_ext,
            dq_ext=l_dq_ext, fl_mult-)
        if (gireduce.status != 0)
            goto error
        for (i=1; i<=nin; i+=1) {
            infile[i] = "q"//infile[i]
        }
    }

    # Extract spectra
    if (l_fl_extract) {
        delete (temp1, verify-, >& "dev$null")
        for (i=1; i<=nin; i+=1) {
            gfextract (infile[i], outimage="", outpref="e",
                reference=l_reference, response=l_response, exslits=l_exslits,
                bpmfile=l_bpmfile, grow=l_grow, line=l_line, nsum=l_nsum,
                trace=l_trace, recenter=l_recenter, function=l_function,
                order=l_t_order, t_nsum=l_nsum, weights=l_weights,
                thresh=l_thresh, gratingdb=l_gratingdb, filterdb=l_filterdb,
                xoffset=l_xoffset, fl_gsappwave=l_fl_gsappwave,
                fl_vardq=l_fl_vardq, fl_inter=l_fl_inter, verbose=l_verbose,
                logfile=l_logfile, sci_ext=l_sci_ext, var_ext=l_var_ext,
                dq_ext=l_dq_ext, fl_gnsskysub=l_fl_gnsskysub,
                fl_fixnc=l_fl_fixnc, fl_fixgaps=l_fl_fixgaps,
                fl_novlap=l_fl_novlap, perovlap=l_perovlap, \
                fl_fulldq=l_fl_fulldq, dqthresh=l_dqthresh)
            if (gfextract.status != 0)
                goto error
            infile[i] = "e"//infile[i]
        }
    }

    # Wavelength calibration
    if (l_fl_wavtran) {
        delete (temp1, verify-, >& "dev$null")
        for (i=1; i<=nin; i+=1) {
            print (infile[i], >> temp1)
        }
        gftransform ("@"//temp1, outimage="", outpre="t",
            wavtraname=l_wavtraname, database="database", fl_vardq=l_fl_vardq,
            fl_flux+, logfile=l_logfile, sci_ext=l_sci_ext, var_ext=l_var_ext,
            dq_ext=l_dq_ext, w1=l_w1, w2=l_w2, dw=l_dw, nw=l_nw,
            verbose=l_verbose)
        if (gftransform.status != 0)
            goto error
        for (i=1; i<=nin; i+=1) {
            infile[i] = "t"//infile[i]
        }
    }

    # Sky subtraction
    if (l_fl_skysub) {
        delete (temp1, verify-, >& "dev$null")
        for (i=1; i<=nin; i+=1) {
            print (infile[i], >> temp1)
        }
        gfskysub ("@"//temp1, outimages="", outpre="s", expr=l_expr,
            sepslits=l_sepslits, logfile=l_logfile, verbose=l_verbose,
            fl_inter=l_fl_inter, sci_ext=l_sci_ext, var_ext=l_var_ext,
            dq_ext=l_dq_ext) 
        if (gfskysub.status != 0)
            goto error
        for (i=1; i<=nin; i+=1) {
            infile[i] = "s"//infile[i]
        }
    }

    # Flux calibration
    if (l_fl_fluxcal) {
        delete (temp1, verify-, >& "dev$null")
        for (i=1; i<=nin; i+=1) {
            print (infile[i], >> temp1)
        }
        l_fl_ext = no
        if (access(l_extinct))
            l_fl_ext = yes
        gscalibrate ("@"//temp1, output="", outpref="c", sfunction=l_sfunc,
            fl_vardq=l_fl_vardq, fl_flux=yes, fl_ext=l_fl_ext,
            extinction=l_extinct, logfile=l_logfile, verbose=l_verbose,
            sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext,
            observatory=l_obs)
        if (gscalibrate.status != 0)
            goto error
        for (i=1; i<=nin; i+=1) {
            infile[i] = "c"//infile[i]
        }
    }

    # Rename output files
    if (access(temp2)) {
        printlog ("", l_logfile, l_verbose)
        printlog ("Renaming output files", l_logfile, l_verbose)
        i = 0
        scanfile = temp2
        while (fscan(scanfile,img) !=EOF) {
            i = i+1
            imrename (infile[i], img, verbose=l_verbose)

            gemdate ()
            gemhedit (img//"[0]", "GFREDUCE", gemdate.outdate,
                "UT Time stamp for GFREDUCE", delete-)
            gemhedit (img//"[0]", "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI", delete-)

            # also need to rename aperture and other database files
            if (l_fl_extract) {
                fparse (img, verbose-)
                outimg = fparse.root
                j = stridx("e",infile[i])
                fparse (infile[i], verbose-)
                root = fparse.root
                root = substr(root,j,(strlen(root)))
                for (slitnum=j1; slitnum<=nslit; slitnum+=1) {
                    sed ("-e","s/"//root//"/"//outimg//"/g", \
                        "database/ap"//root//"_"//slitnum, > \
                        "database/ap"//outimg//"_"//slitnum)
                    delete ("database/ap"//root//"_"//slitnum, verify-, \
                        >& "dev$null")
                    if (fl_vardq) {
                        sed ("-e","s/"//root//"/"//outimg//"/g", \
                            "database/ap"//root//"_var_"//slitnum, > \
                            "database/ap"//outimg//"_var_"//slitnum)
                        delete ("database/ap"//root//"_var_"//slitnum, \
                            verify-, >& "dev$null")
                        sed ("-e","s/"//root//"/"//outimg//"/g", \
                            "database/ap"//root//"_dq_"//slitnum, > \
                            "database/ap"//outimg//"_dq_"//slitnum)
                        delete ("database/ap"//root//"_dq_"//slitnum, \
                            verify-, >& "dev$null")
                    }
                }
            }
        }
        scanfile = ""
    }

    ### Note: if 'outpref="default"' there is no easy way to add the
    ###       the GFREDUCE time stamp to the output images.  The output
    ###       images will contain only the time stamps from the GMOS tasks
    ###       used to reduce the data.

    goto clean

error:
    status = 1

clean:
    delete (temp1//","//temp2, verify-, >& "dev$null")
    printlog ("",l_logfile,l_verbose)
    scanfile = ""
    if (status==0)
        printlog ("GFREDUCE exit status: good", l_logfile, l_verbose)
    else
        printlog ("GFREDUCE exit status: error", l_logfile, l_verbose)
    printlog ("------------------------------------------------------------\
        -------------------", l_logfile, l_verbose)

end
