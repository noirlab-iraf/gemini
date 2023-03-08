# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gsflat (inflats,specflat)

# Create a GMOS normalized spectral flatfield (LONGSLIT or MOS) from GCAL flats
#
# Version    Feb 28, 2002  ML,BM,IJ  v1.3 release
#            Aug 15, 2002  BM  add x/yoffset parameters for new gscut
#            Sept 20, 2002    v1.4 release
#            Oct 14, 2002  IJ don't use "i" inside a script
#            Aug 20, 2003  MB add yadd parameter for new gscut
#                             allow detector by detector fitting
#            Sep 11, 2003  KL IRAF2.12 - new/modified parameters
#                             hedit - addonly
#                             imcombine - headers,bpmasks,expmasks,outlimits
#                                 rejmask->rejmasks, plfile->nrejmasks
#            Sep 14, 2004  BM,IJ edge finding of slits, fl_vardq support,
#                                cleaning
#            Mar 10, 2009  JH Add double flat mode for nod-and-shuffle
string  inflats     {prompt="Input flatfields"}
string  specflat    {prompt="Output normalized flat (MEF)"}
bool    fl_slitcorr {no,prompt="Correct output for Illumination/Slit-Function"}
string  slitfunc    {"",prompt="Slit Function (MEF output of gsslitfunc)"}
bool    fl_keep     {no,prompt="Keep imcombined flat?"}
string  combflat    {"",prompt="Filename for imcombined flat"}

bool    fl_over     {yes,prompt="Subtract overscan level"}
bool    fl_trim     {yes,prompt="Trim off overscan region"}
bool    fl_bias     {yes,prompt="Subtract bias image"}
bool    fl_dark     {no,prompt="Subtract (scaled) dark image"}
bool    fl_qecorr   {no, prompt="QE correct the input images?"}
bool    fl_fixpix   {yes,prompt="Interpolate across chip gaps"}
bool    fl_oversize {yes,prompt="Use 1.05x slit length to accommodate distortion?"}
bool    fl_vardq    {no,prompt="Create variance and data quality frames"}
bool    fl_fulldq   {no, prompt="Decompose DQ during gmosaic transformation; apply correct handling?"}            # OLDP-3
real    dqthresh    {0.1, min=0.01, max=0.5, prompt="Threshold applied to DQ when fl_fulldq=yes"}
string  bias        {"",prompt="Bias image"}
string  dark        {"",prompt="Dark image"}
string  key_exptime {"EXPTIME",prompt="Exposure time header keyword"}
string  key_biassec {"BIASSEC",prompt="Header keyword for overscan strip image section"}
string  key_datasec {"DATASEC",prompt="Header keyword for data section (excludes the overscan)"}
string  rawpath     {"",prompt="GPREPARE: Path for input raw images"}
string  sci_ext     {"SCI",prompt="Name of science extension"}
string  var_ext     {"VAR",prompt="Name of variance extension"}
string  dq_ext      {"DQ",prompt="Name of data quality extension"}
string  key_mdf     {"MASKNAME",prompt="Header keyword for the MDF"}
string  mdffile     {"",prompt="MDF to use if keyword not found"}
string  mdfdir      {"gmos$data/", prompt="MDF database directory"}
string  bpm         {"",prompt="Name of bad pixel mask file or image"}
string  gaindb      {"default",prompt="Database with gain data"}
string  gratingdb   {"gmos$data/GMOSgratings.dat",prompt="Gratings database file"}
string  filterdb    {"gmos$data/GMOSfilters.dat",prompt="Filters database file"}
string  bpmfile     {"gmos$data/chipgaps.dat",prompt="Info on location of chip gaps"}
char    refimage        {"",prompt="Reference image for slit positions"}
char    qe_refim        {"", prompt="QE wavelength reference image."}
bool    fl_keep_qeim    {yes, prompt="Keep QE correction?"}
char    qe_corrpref     {"qecorr", prompt="Prefix for QE correction files."}
char    qe_corrimages   {"", prompt="Name for QE correction data."}
char    qe_data         {"gmosQEfactors.dat", prompt="Data file that contains QE information."}
char    qe_datadir      {"gmos$data/", prompt="Directory containg QE data file."}
string  sat         {"default",prompt="Saturation level in raw images"}
real    xoffset     {INDEF,prompt="X offset in wavelength [nm]"}
real    yoffset     {INDEF,prompt="Y offset in unbinned pixels"}
real    yadd        {0,prompt="Additional pixels to add to each end of MOS slitlet lengths"}
real    wave_limit  {INDEF,prompt="Upper wavelength limit of cut spectra (nm). Only for use with fl_detec=no."}
bool    fl_usegrad  {no,prompt="Use gradient method to find MOS slits"}
bool    fl_emis     {no,prompt="mask emission lines from lamp (affected pixels set to 1. in output)"}
string  nbiascontam {"default", prompt="Number of columns removed from overscan region"}
string  biasrows    {"default", prompt="Rows to use for overscan region"}
real    minval     {INDEF, max=1.0, prompt="Minimum pixel value in normalized flat"}
#Response/Profile fitting
bool    fl_inter    {no,prompt="Fit response interactively?"}
bool    fl_answer   {yes,prompt="Continue interactive fitting?"}
bool    fl_detec    {no, prompt="Fit response detector by detector rather than slit by slit?"}
bool    fl_seprows  {yes, prompt="Fit and normalize each row separately?"}
string  function    {"spline3",min="spline3|legendre|chebyshev|spline1", prompt="Fitting function for response"}
string  order       {"15",prompt="Order of fitting function, minimum value=1"}
real    low_reject  {3,prompt="Low rejection in sigma of response fit"}
real    high_reject {3,prompt="High rejection in sigma of response fit"}
int     niterate    {2,prompt="Number of rejection iterations in response fit"}

#Gemcombine
string  combine     {"average",prompt="Combination operation"}
string  reject      {"avsigclip",prompt="Rejection algorithm"}
string  masktype    {"goodvalue", enum="none|goodvalue", prompt="Mask type"}
real    maskvalue   {0., prompt="Mask value"}
string  scale       {"mean",prompt="Image scaling"}
string  zero        {"none",prompt="Image zeropoint offset"}
string  weight      {"none",prompt="Image weights"}
string  statsec     {"",prompt="Statistics section"}
real    lthreshold  {INDEF,prompt="Lower threshold"}
real    hthreshold  {INDEF,prompt="Upper threshold"}
int     nlow        {1,min=0,prompt="minmax: Number of low pixels to reject"}
int     nhigh       {1,min=0,prompt="minmax: Number of high pixels to reject"}
int     nkeep       {0,min=0,prompt="Minimum to keep or maximum to reject"}
bool    mclip       {yes,prompt="Use median in sigma clipping algorithms?"}
real    lsigma      {3.,prompt="Lower sigma clipping factor"}
real    hsigma      {3.,prompt="Upper sigma clipping factor"}
string  key_ron     {"RDNOISE",prompt="Keyword for readout noise in e-"}
string  key_gain    {"GAIN",prompt="Keyword for gain in electrons/ADU"}
real    ron         {3.5,min=0.,prompt="Readout noise rms in electrons"}
real    gain        {2.2,min=0.00001,prompt="Gain in e-/ADU"}
string  snoise      {"0.0",prompt="ccdclip: Sensitivity noise (electrons)"}
real    sigscale    {0.1,prompt="Tolerance for sigma clipping scaling correction"}
real    pclip       {-0.5,prompt="pclip: Percentile clipping parameter"}
real    grow        {0.0,prompt="Radius (pixels) for neighbor rejection"}

#Colbias
bool    ovs_flinter {no,prompt="Interactive overscan fitting?"}
bool    ovs_med     {no,prompt="Use median instead of average in column bias?"}
string  ovs_func    {"chebyshev",min="spline3|legendre|chebyshev|spline1", prompt="Overscan fitting function"}
char    ovs_order   {"default",prompt="Order of overscan fitting function"}
real    ovs_lowr    {3.,prompt="Low sigma rejection factor"}
real    ovs_highr   {3.,prompt="High sigma rejection factor"}
int     ovs_niter   {2,prompt="Number of rejection iterations"}

# Nod and shuffle
bool    fl_double     {no,prompt="Make double flats for nod-and-shuffle science"}
int     nshuffle      {0,prompt="Number of shuffle pixels (unbinned)"}

#General
string  logfile     {"",prompt="Logfile name"}
bool    verbose     {yes,prompt="Verbose"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

begin

    #local variables of parameters
    string  l_inflats, l_specflat, l_logfile, l_slitfunc
    string  l_bias, l_dark, l_bpm, l_key_biassec, l_key_datasec, l_ovs_func
    string  l_combflat, l_bpmfile, l_refimage
    string  l_key_mdf, l_mdffile, l_combine, l_reject, l_scale, l_zero
    string  l_weight, l_statsec, l_key_exptime
    string  l_key_ron, l_key_gain
    string  l_snoise, l_masktype, l_sat, l_ovs_order
    string  l_rfunction, l_sci_ext, l_var_ext, l_dq_ext
    string  l_mdfdir, l_gaindb, l_rawpath, l_gratingdb, l_filterdb
    string  l_order, l_biasrows, l_nbiascontam, pathtest, testfile
    real    l_lthreshold, l_hthreshold, l_lsigma, l_hsigma,l_ron,l_gain
    real    l_sigscale, l_maskvalue
    real    l_pclip, l_grow, l_ovs_lowr, l_ovs_highr
    real    l_rlowrej, l_rhighrej
    real    l_xoffset, l_yoffset, l_yadd, l_upper_wave_limit
    int     l_nlow, l_nhigh, l_nkeep, l_rorder, l_rdorder[12]
    int     l_ovs_niter, l_niterate, l_nshuffle, atposition
    bool    l_slitcorr, l_mclip, l_rinteractive, l_fl_vardq, l_verbose
    bool    l_redkeep
    bool    l_ovs_flinter,l_ovs_med,l_fl_over,l_fl_bias,l_fl_dark,l_fl_trim
    bool    l_fl_fixpix, l_fl_oversize
    bool    l_fl_answer
    bool    l_fl_detec, l_fl_seprows
    bool    l_fl_emis, l_fl_usegrad
    bool    l_fl_gmosaic
    bool    l_fl_double

    #other parameters used within this task
    file    temp1, temp2
    bool    reducing, gpflag, giflag, gsflag, gsfl1, mosflag, mosfl1, florder
    bool    specmode, ishamamatsu
    string  filt1, filt2, filtn1, filtn2, grat1, gratn, cutfile, l_gradimage
    string  grtilt, grtilt1
    string  mdf, inlist, filelist, scilist, img, specsec, imgraw, slitsec
    string  response, suf, combsig
    string  infile[20], headtest
    string  scisec, sciflat, rowavg, flatfit, inscoo, combflatsciext, obsmode
    string  gprep, gired, specsecfile, tmpflat, tmpratsh, tmpspecflat
    real    firstexp, l_sgain, l_sron, gaineff, roneff, expdiff, l_dqthresh, xr
    int     colpos1, colpos2, compos, seclen, xref, x1, y1, x2, y2, nextens
    int     next1, nsciext
    int     n_i, n_j, n_k, l, n, idum, nbad, len, nslits, nim, n_ccd, n_row
    bool    fl_delgrad, l_fl_gsappwave, l_fl_fulldq
    struct  sdate
    char    l_sample, rphend
    int     x11, x12, x21, x22, Xbin, Xmax, ymin, ymax, ycen
    int     norder, gapsize, gapoff, nbcols, fgcol, lgcol, lastcol
    char    keyfound, instrument
    int     junk, xbin, ybin

    # QE correction variables
    char    l_qe_refim, l_qe_corrpref, qecorrapp, l_qecorr_data
    char    l_qecorr_datadir, l_qe_corrimages, l_key_qecorrim
    bool    l_fl_qecorr, l_fl_keep_qeim, flatqecorr_state, inqecorr_state
    bool    l_fl_flat_orig, l_giflat_orig
    int     prev_qecorr

    #Make parameter assignments
    l_inflats=inflats ; l_specflat=specflat
    l_combflat=combflat ; l_redkeep=fl_keep
    l_fl_over=fl_over ; l_fl_bias=fl_bias ; l_fl_dark=fl_dark
    l_fl_fixpix=fl_fixpix ; l_fl_oversize=fl_oversize
    l_bias=bias ; l_dark=dark ; l_bpm=bpm ; l_bpmfile=bpmfile
    l_key_mdf=key_mdf ; l_mdffile=mdffile
    l_key_biassec=key_biassec ; l_key_datasec=key_datasec
    l_ovs_flinter=ovs_flinter ; l_ovs_med=ovs_med
    l_ovs_func=ovs_func ; l_ovs_order=ovs_order
    l_slitcorr=fl_slitcorr ; l_slitfunc=slitfunc
    l_combine=combine ; l_reject=reject
    l_scale=scale ; l_zero=zero ; l_weight=weight
    l_statsec=statsec ; l_key_exptime=key_exptime
    l_lthreshold=lthreshold ; l_hthreshold=hthreshold
    l_nlow=nlow ; l_nhigh=nhigh ; l_nkeep=nkeep
    l_lsigma=lsigma ; l_hsigma=hsigma
    l_verbose=verbose ; l_mclip=mclip
    junk = fscan (key_gain, l_key_gain)
    junk = fscan (key_ron, l_key_ron)
    l_gain=gain ; l_ron=ron ; l_snoise=snoise
    l_sigscale=sigscale ; l_pclip=pclip ; l_grow=grow
    l_sci_ext = sci_ext ; l_dq_ext = dq_ext ; l_var_ext=var_ext
    l_ovs_lowr=ovs_lowr ; l_ovs_highr=ovs_highr ; l_ovs_niter=ovs_niter
    l_fl_vardq=fl_vardq ; l_biasrows=biasrows
    l_rinteractive=fl_inter ; l_rfunction=function;
    l_nbiascontam = nbiascontam ; l_order=order
    l_rlowrej=low_reject ; l_rhighrej=high_reject; l_niterate=niterate
    l_logfile=logfile; l_sat=sat
    l_fl_trim=fl_trim
    l_mdfdir=mdfdir ; l_gaindb=gaindb ; l_rawpath=rawpath
    l_gratingdb=gratingdb ; l_filterdb=filterdb
    l_xoffset=xoffset ; l_yoffset=yoffset; l_yadd=yadd
    l_upper_wave_limit = wave_limit
    l_fl_emis=fl_emis
    l_fl_detec=fl_detec
    l_fl_seprows=fl_seprows
    l_refimage=refimage ; l_fl_usegrad=fl_usegrad
    l_masktype = masktype
    l_maskvalue = maskvalue ; l_fl_double = fl_double
    l_nshuffle = nshuffle

    # DQ handling in gmosaic
    l_fl_fulldq = fl_fulldq
    l_dqthresh = dqthresh

    # Set QE parameters
    l_key_qecorrim = "QECORRIM"

    # Read QE parameters
    l_fl_qecorr = fl_qecorr
    l_qe_refim = qe_refim
    l_fl_keep_qeim = fl_keep_qeim
    l_qe_corrpref = qe_corrpref
    l_qecorr_data = qe_data
    l_qecorr_datadir = qe_datadir
    l_qe_corrimages = qe_corrimages

    status = 0
    if (l_fl_detec) {
        l_fl_gmosaic = no
        l_fl_fixpix = no
        l_fl_gsappwave = no
    } else {
        l_fl_gmosaic = yes
        l_fl_gsappwave = yes
    }

    #Emission Line Masking not enabled yet
    if (l_fl_emis) {
        printlog ("WARNING - GSFLAT: Emission line masking is not yet \
            implemented", l_logfile, verbose+)
        printlog ("                  Setting fl_emis = no", l_logfile, \
            verbose+)
        l_fl_emis = no
    }

    # Keep some parameters from changing by outside world
    cache ("imgets", "gemhedit", "gmos", "tinfo", "fparse", "gimverify", \
        "gemdate")

    #Define temporary files
    specsecfile = mktemp("tmpsecfile")
    scilist = mktemp("tmpscilist")
    tmpflat = mktemp("tmpcombflat")
    scisec = mktemp("tmpscisec")
    if (!l_redkeep)
        l_combflat = mktemp("tmpcombflat")
    if (l_fl_vardq)
        combsig = mktemp("tmpsig")
    response = mktemp("tmpresponse")
    temp1 = mktemp("tmpfilelist")
    temp2 = mktemp("tmpfilelist")
    mdf = mktemp("tmpmdffile")//".fits"

    #Check logfile ...
    if (l_logfile == "STDOUT") {
        l_logfile = ""
        l_verbose = yes
    } else if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GSFLAT: both gsflat.logfile and gmos.logfile \
                are empty.", l_logfile, l_verbose)
            printlog ("                  Using default file gmos.log.",
                l_logfile, l_verbose)
        }
    }
    date | scan(sdate)
    printlog ("",l_logfile,l_verbose)
    printlog ("------------------------------------------------------------\
        --------------------", l_logfile, l_verbose)
    printlog ("GSFLAT -- "//sdate,l_logfile,l_verbose)
    printlog ("",l_logfile,l_verbose)
    printlog ("Input images or list              = "//l_inflats,
        l_logfile,l_verbose)
    printlog ("Output spectral flat              = "//l_specflat,
        l_logfile,l_verbose)
    printlog ("Correct for slit function         = "//l_slitcorr,
        l_logfile,l_verbose)
    if (l_slitcorr)
        printlog ("Slitfunc image from GSSLITFUNCTION= "//l_slitfunc,
            l_logfile,l_verbose)
    printlog ("Keep combined flat combflat       = "//l_redkeep,
        l_logfile,l_verbose)
    if (l_redkeep)
        printlog ("Combined flat (not normalized)    = "//l_combflat,
            l_logfile,l_verbose)
    printlog("",l_logfile,l_verbose)
    printlog("Fitting parameters for Spectral Flatfield: ",l_logfile,l_verbose)
    printlog("  interactive = "//l_rinteractive,l_logfile,l_verbose)
    printlog("  function    = "//l_rfunction,l_logfile,l_verbose)
    printlog("  order       = "//l_order,l_logfile,l_verbose)
    printlog("  low_reject  = "//l_rlowrej,l_logfile,l_verbose)
    printlog("  high_reject = "//l_rhighrej,l_logfile,l_verbose)
    printlog("  niterate    = "//l_niterate,l_logfile,l_verbose)
    printlog("  Fit detectors separately fl_detec = "//l_fl_detec,
        l_logfile,l_verbose)
    printlog("",l_logfile,l_verbose)
    printlog("Overscan Subtraction fl_over   = "//l_fl_over, l_logfile, \
        l_verbose)
    printlog("Trim image           fl_trim   = "//l_fl_trim, l_logfile, \
        l_verbose)
    printlog("Bias Subtraction     fl_bias   = "//l_fl_bias, l_logfile, \
        l_verbose)
    printlog("Dark Subtraction     fl_dark   = "//l_fl_dark, l_logfile, \
        l_verbose)
    printlog("QE Correction        fl_qecorr = "//l_fl_qecorr, l_logfile, \
        l_verbose)
    printlog("VAR & DQ planes      fl_vardq  = "//l_fl_vardq, l_logfile, \
        l_verbose)
    printlog("Fixpix chip gaps     fl_fixpix = "//l_fl_fixpix, l_logfile, \
        l_verbose)
    printlog("Oversize slit length fl_oversize = "//l_fl_oversize, l_logfile, \
        l_verbose)
    printlog("Fit detec. by detec. fl_detec  = "//l_fl_detec, l_logfile, \
        l_verbose)
    printlog("Fit rows separately  fl_seprows = "//l_fl_seprows, \
        l_logfile,l_verbose)
    printlog("Use gradient method  fl_usegrad= "//l_fl_usegrad, l_logfile, \
        l_verbose)
    if (l_upper_wave_limit != INDEF) {
        printlog("User defined upper wavelength limit wave_limit = "//\
            l_upper_wave_limit, l_logfile,l_verbose)
    }
    printlog("Mask emission lines  fl_emis   = "//l_fl_emis, l_logfile, \
        l_verbose)
    printlog("Reference image      refimage  = "//l_refimage, l_logfile, \
        l_verbose)
    printlog("N&S double flatfield fl_double = "//l_fl_double, l_logfile, \
        l_verbose)
    printlog("",l_logfile,l_verbose)
    printlog("bias        = "//l_bias,l_logfile,l_verbose)
    printlog("dark        = "//l_dark,l_logfile,l_verbose)
    printlog("qe_refim    = "//l_qe_refim,l_logfile,l_verbose)
    printlog("qe_corrim   = "//l_qe_corrimages,l_logfile,l_verbose)
    printlog("qe_corrpref = "//l_qe_corrpref,l_logfile,l_verbose)
    printlog("sci_ext     = "//l_sci_ext,l_logfile,l_verbose)
    printlog("var_ext     = "//l_var_ext,l_logfile,l_verbose)
    printlog("dq_ext      = "//l_dq_ext,l_logfile,l_verbose)
    printlog("bpm         = "//l_bpm,l_logfile,l_verbose)
    printlog("key_mdf     = "//l_key_mdf,l_logfile,l_verbose)
    printlog("mdffile     = "//l_mdffile,l_logfile,l_verbose)
    printlog("mdfdir      = "//l_mdfdir,l_logfile,l_verbose)
    printlog("Chip gaps   = "//l_bpmfile,l_logfile,l_verbose)
    printlog("rawpath     = "//l_rawpath,l_logfile,l_verbose)
    printlog("",l_logfile,l_verbose)

    nbad = 0

    #The checks on things begin ...

    # Check that the rawpath has a trailing slash and is a valid entry
    if (l_rawpath != "") {
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
            nbad += 1
            goto error
        }
    }

    if (l_inflats == "" || l_inflats==" ") {
        printlog ("ERROR - GSFLAT: Input Flats not specified", l_logfile, \
            verbose+)
        nbad = nbad+1
    }

    # check existence of list files
    atposition = stridx("@",l_inflats)
    if (atposition > 0) {
        testfile = substr(l_inflats,atposition+1,strlen(l_inflats))
        if (!access(testfile)) {
            printlog ("ERROR - GSFLAT: The input list "//testfile//\
                " does not exist.", logfile=l_logfile, verbose+)
            goto clean
        }
    }

    if (l_specflat == "" || l_specflat==" ") {
        printlog ("ERROR - GSFLAT: output specflat not specified",
            l_logfile, verbose+)
        nbad = nbad+1
    }
    if (l_redkeep && (l_combflat == "" || l_combflat==" ")) {
        printlog ("ERROR - GSFLAT: output combflat not specified",
            l_logfile, verbose+)
        nbad = nbad+1
    }

    # Check that the output file does not already exist. If so, exit.
    gimverify(l_specflat)
    l_specflat=gimverify.outname
    if (gimverify.status != 1) {
        printlog("ERROR - GSFLAT: Output file "//l_specflat//" already \
            exists.", l_logfile, verbose+)
        nbad=nbad+1
    }

    # Check that the output combined flatfield does not already exist. If so,
    # exit.
    gimverify (l_combflat) ; l_combflat=gimverify.outname
    if (gimverify.status != 1) {
        printlog ("ERROR - GSFLAT: Output file "//l_combflat//" already \
            exists.", l_logfile,verbose+)
        nbad = nbad+1
    }

    if (l_sci_ext =="" || l_sci_ext ==" ") {
        printlog ("ERROR - GSFLAT: extension name sci_ext is missing",
            l_logfile,verbose+)
        nbad = nbad+1
    }

    #If dq/var propogation is requested, make sure the names are given
    if(l_fl_vardq) {
        if (l_dq_ext=="" || l_dq_ext ==" ") {
            printlog ("ERROR - GSFLAT: extension name dq_ext is missing",
                l_logfile, verbose+)
            nbad = nbad+1
        } else if (l_var_ext=="" || l_var_ext ==" ") {
            printlog ("ERROR - GSFLAT: extension name var_ext is missing",
                l_logfile, verbose+)
            nbad = nbad+1
        }
    }

    #check no commas in sci_ext, var_ext and dq_ext
    if (stridx(",",l_sci_ext)>0 || stridx(",",l_var_ext)>0 ||
        stridx(",",l_dq_ext)>0 ) {

        printlog ("ERROR - GSFLAT: sci_ext, var_ext or dq_ext contains \
            commas, give root name only", l_logfile, verbose+)
        nbad = nbad+1
    }

    #Now check for slit-function file if slitcorr=yes is specified in the
    #parameters
    if (l_slitcorr) {
        gimverify (l_slitfunc) ; l_slitfunc = gimverify.outname
        if (gimverify.status != 0) {
            printlog ("ERROR - GSFLAT: Input slit-function file : "//\
                l_slitfunc//" either does not exist", l_logfile, verbose+)
            printlog ("                or is not a valid MEF.", l_logfile, \
                verbose+)
            nbad = nbad+1
        } else {
            imgets (l_slitfunc//"[0]", "GSSLITFU")
            if (imgets.value == "0" || imgets.value == " ") {
                printlog ("ERROR - GSFLAT: Input file "//l_slitfunc//" was \
                    not produced by GSSLITFUNCTION.", l_logfile, verbose+)
                nbad = nbad+1
            }
        }
    }
    # slit correction not yet implemented for detector by detector fitting
    if (l_slitcorr && l_fl_detec) {
        printlog ("ERROR - GSFLAT: Slit Function correction not implemented \
            for detector by detector fitting.", l_logfile, verbose+)
        nbad = nbad+1
    }

    #check on existence of grating and filter databases, and chipgaps
    if (!access(l_gratingdb)) {
        printlog ("ERROR - GSFLAT: gratings database file : "//l_gratingdb//\
            " does not exist.", l_logfile, verbose+)
        nbad = nbad+1
    }
    if (!access(l_filterdb)) {
        printlog ("ERROR - GSFLAT: filters database file : "//l_filterdb//\
            " does not exist.", l_logfile, verbose+)
        nbad = nbad+1
    }

    if ((l_gaindb != "") && (l_gaindb != "default")) {
        if (access(l_gaindb) == no) {
            printlog ("ERROR - GSFLAT: gain database file : "//l_gaindb//\
                " does not exist.", l_logfile, verbose+)
            nbad = nbad+1
        }
    }

    if (l_bpmfile != "" && !access(l_bpmfile)) {
        printlog ("ERROR - GSFLAT: Chip gap definition file : "//l_bpmfile//\
            " does not exist.", l_logfile, verbose+)
        nbad = nbad+1
    }

    # Check if user has specified an upper wavelength limit whith fl_detec+
    if (l_upper_wave_limit != INDEF && fl_detec) {
        printlog ("WARNING - GSFLAT: wave_limit has been set to "//\
            l_upper_wave_limit//" but fl_detec=yes.\n"//\
            "                  wave_limit will be ignored.", \
            l_logfile, verbose+)
    }

    #check that the order for fitting the normalization is rational for
    #requested reduction technique (slit vs. detector at a time)
    norder = 0
    for (n_i=1; n_i<=12; n_i+=1)
        l_rdorder[n_i] = 0

    print (l_order) | tokens("STDIN", newl-) | match(",","STDIN",stop+) | \
        count("STDIN") | scan(norder)
    print (l_order) | tokens("STDIN", newl-) | match(",","STDIN",stop+) | \
        fields("STDIN",1,lines="1") | scan(l_rorder)
    if (l_rorder < 1) {
        printlog ("ERROR - GSFLAT: order of fitting function must be greater \
            than zero", l_logfile, verbose+)
        nbad = nbad+1
    }
    if (norder!=1 && norder!=3 && norder!=6 && norder!=12) {
        printlog ("ERROR - GSFLAT: order of fitting function must contain \
            either 1, 3, 6 or 12 elements.", l_logfile, verbose+)
        nbad = nbad+1
    } else if((norder==3 || norder==6) && !l_fl_detec && (l_rorder > 0)) {
        printlog ("WARNING - GSFLAT: Fitting normalization slit by slit; \
            will use only first", l_logfile, verbose+)
        printlog ("                  element, ["//str(l_rorder)//"], of \
            fitting order array, ["//l_order//"]", l_logfile,verbose+)
    } else if (norder==1 && l_fl_detec) {
        for (n_i=1; n_i<=12; n_i+=1)
            l_rdorder[n_i] = l_rorder
    } else if ((norder==3 || norder==6 || norder==12) && (l_rorder > 0)) {
        florder = yes
        for (n_i=1; n_i<=norder; n_i+=1) {
            print (l_order) | tokens("STDIN", newl-) | \
                match(",","STDIN",stop+) | fields("STDIN",1,lines=str(n_i)) | \
                scan(l_rdorder[n_i])
            if (l_rdorder[n_i] < 1)
                florder=no
        }
        if (florder==no) {
            printlog ("ERROR - GSFLAT: order of fitting function must be \
                greater than zero", l_logfile, verbose+)
            nbad = nbad+1
        }
    }

    #Assuming we're ok so far, check that all images in the list exist, and
    # are of the correct filetype

    if (nbad == 0) {

        #Create list of input flat images
        files (l_inflats, > temp1)
        scanfile = temp1

        nim = 0
        while (fscan(scanfile,img) != EOF) {
            nim = nim+1
            gimverify (img) ; img = gimverify.outname//".fits"
            imgraw = gimverify.outname//".fits"
            if (gimverify.status>=1 && l_rawpath!="") {
                gimverify (l_rawpath//img)
                imgraw = gimverify.outname//".fits"
            }
            if (gimverify.status != 0) {
                if (gimverify.status == 1) {
                    printlog ("ERROR - GSFLAT: Input image : "//img//" does \
                        not exist.", l_logfile, verbose+)
                    nbad = nbad+1
                } else if (gimverify.status == 2) {
                    printlog ("ERROR - GSFLAT: Input image : "//img//" is in \
                        OIF format, not MEF.", l_logfile, verbose+)
                    nbad = nbad+1
                } else if (gimverify.status == 3) {
                    printlog ("ERROR - GSFLAT: Input image : "//img//" is a \
                        GEIS file, not MEF.", l_logfile, verbose+)
                    nbad = nbad+1
                } else if (gimverify.status == 4) {
                    printlog ("ERROR - GSFLAT: Input image : "//img//" is a \
                        simple FITS file, not MEF.", l_logfile, verbose+)
                    nbad = nbad+1
                }
            }

        }
        # Exit if problems found
        if (nbad > 0) {
            printlog ("ERROR - GSFLAT: "//nbad//" fatal errors found.",
                l_logfile, verbose+)
            goto error
        }

        scanfile = temp1

        nim = 0
        while (fscan(scanfile,img) != EOF) {
            nim = nim+1
            gimverify (img) ; img = gimverify.outname//".fits"
            imgraw = gimverify.outname//".fits"
            if (gimverify.status>=1 && l_rawpath!="") {
                gimverify(l_rawpath//img)
                imgraw=gimverify.outname//".fits"
            }

            # Check for existing GSREDUCEd image
            if (imaccess("gs"//img)) {
                imgets ("gs"//img//"[0]","GSREDUCE")
                if (imgets.value != "" && imgets.value != " " && \
                    imgets.value != "0") {
                    printlog ("GSFLAT: Using image gs"//img//" previously \
                        processed with GSREDUCE", l_logfile, l_verbose)
                    img = "gs"//img
                    infile[nim] = img
                    if (l_fl_vardq && !imaccess(img//"["//l_var_ext//",1]")) {
                        printlog ("WARNING - GSFLAT: Image "//img,\
                            l_logfile, verbose+)
                        printlog ("                  does not have a \
                            variance plane.", l_logfile, verbose+)
                        printlog ("                  turning off \
                            var/dq propagation.", l_logfile, verbose+)
                        l_fl_vardq = no
                    }
                } else {
                    printlog ("ERROR - GSFLAT: Found output image gs"//img//\
                        " which is not GSREDUCEd", l_logfile, verbose+)
                    nbad = nbad+1
                }
            } else {
                infile[nim] = img
                img = imgraw
            }

            # Based on first img, check if Hamamatsu and assign ranges from
            # which l_sample is constructed, using info from chipgaps.dat
            if (nim == 1) {
                imgets (img//"[0]", "DETTYPE", >& "dev$null")
                if (imgets.value == "S10892" || imgets.value == "S10892-N") { # Hamamatsu data
                    ishamamatsu = yes
                    imgets (img//"[0]", "INSTRUME", >& "dev$null")
                    instrument = imgets.value
                    if (l_bpmfile != "") {
                        if (instrument == "GMOS-S") {
                            fields (l_bpmfile,"1,2",lines="5") | scan(x11,x12)
                            fields (l_bpmfile,"1,2",lines="6") | scan(x21,x22)
                        } else {  # GMOS-N
                            fields (l_bpmfile,"1,2",lines="8") | scan(x11,x12)
                            fields (l_bpmfile,"1,2",lines="9") | scan(x21,x22)
                        }
                    }
                } else {
                    ishamamatsu = no
                    if (l_bpmfile != "") {
                        fields (l_bpmfile,"1,2",lines="2") | scan(x11,x12)
                        fields (l_bpmfile,"1,2",lines="3") | scan(x21,x22)
                    }
                }
            }

            nextens = 0
            hselect (img//"[0]","NEXTEND",yes) | scan(nextens)
            if (nextens==0) {
                fxhead (img, format_file="", long_header-, count_lines-) | \
                    match ("STDIN", "IMAGE", stop-) | count | scan(nextens)
            }

            if (nim == 1)
                next1 = nextens
            if (nextens != next1) {
                printlog ("ERROR - GSFLAT: Input images do not all have the \
                    same number of extensions", l_logfile, verbose+)
                nbad = nbad+1
            }

            imgets (img//"[0]", "MASKTYP", >& "dev$null")
            if (imgets.value == "1" )
                specmode = yes
            else
                specmode = no
            if (!specmode) {
                printlog ("ERROR - GSFLAT: Input file "//img//" not of type \
                    LONGSLIT or MOS.", l_logfile, verbose+)
                nbad = nbad+1
            }

            # -------------
            # check if GSREDUCE has been run if reductions are not
            # requested here

            imgets (img//"[0]", "GSREDUCE", >& "dev$null")
            if (imgets.value == " " || imgets.value == "" || \
                imgets.value == "0") {

                gsflag = no
                if (nim==1)
                    gsfl1=no
                else if (gsflag!=gsfl1) {
                    printlog ("ERROR - GSFLAT: Cannot mix raw/GPREPARE'd \
                        images with GSREDUCE'd images", l_logfile, verbose+)
                    nbad = nbad+1
                }
            } else {
                gsflag = yes
                if (nim == 1)
                    gsfl1 = yes
                else if (gsflag!=gsfl1) {
                    printlog ("ERROR - GSFLAT: Cannot mix raw/GPREPARE'd \
                        images with GSREDUCE'd images", l_logfile, verbose+)
                    nbad = nbad+1
                }
            }

            if ((!gsflag) && \
                (!l_fl_over && !l_fl_bias && !l_fl_dark && !l_fl_trim \
                    && !l_fl_qecorr)) {
                printlog ("ERROR - GSFLAT: Flat Image - "//img//" has not \
                    been GSREDUCED", l_logfile, verbose+)
                printlog ("                and no reduction flags are \
                    requested here.", l_logfile, verbose+)
                nbad = nbad+1
            }

            if (gsflag)
                if (l_fl_over || l_fl_bias || l_fl_dark || l_fl_trim || \
                    l_fl_qecorr) {
                    printlog ("WARNING - GSFLAT: Image "//img//" has been \
                        processed with GSREDUCE", l_logfile, verbose+)
                    printlog ("                  Ignoring the processing \
                        flags", l_logfile,verbose+)
                    l_fl_over=no
                    l_fl_bias=no
                    l_fl_dark=no
                    l_fl_trim=no
                }

            imgets (img//"[0]", "GMOSAIC", >& "dev$null")
            if (imgets.value != " " && imgets.value != "" && \
                imgets.value != "0") {
                printlog ("WARNING - GSFLAT: Image "//img//" has been \
                    processed with GMOSAIC", l_logfile, verbose+)
                mosflag = yes
                if (l_fl_over || l_fl_bias || l_fl_dark || l_fl_trim || \
                    l_fl_qecorr) {
                    printlog ("                  Ignoring the processing \
                        flags", l_logfile,verbose+)
                    l_fl_over=no
                    l_fl_bias=no
                    l_fl_dark=no
                    l_fl_trim=no
                } else {
                    printlog ("                  No mosaicing done",
                        l_logfile,verbose+)
                    l_fl_gmosaic = no
                }
                if (nim==1)
                    mosfl1 = yes
                else if (mosflag!=mosfl1) {
                    printlog ("ERROR - GSFLAT: Cannot mix raw/GPREPARE'd \
                        images with GMOSAIC'd images", l_logfile, verbose+)
                    nbad = nbad+1
                }
                if (l_fl_detec) {
                    printlog ("ERROR - GSFLAT: fl_detec set for fitting \
                        detector by detector, but input image is GMOSAIC'd",
                        l_logfile, verbose+)
                    nbad = nbad+1
                }
            } else {
                mosflag = no
                if (nim==1)
                    mosfl1 = no
                else if (mosflag!=mosfl1) {
                    printlog ("ERROR - GSFLAT: Cannot mix raw/GPREPARE'd \
                        images with GMOSAIC'd images", l_logfile, verbose+)
                    nbad = nbad+1
                }
            }

            # -------------
            # check exposure times
            imgets (img//"[0]", l_key_exptime, >& "dev$null")
            if (imgets.value == "0" || imgets.value == " ") {
                printlog ("ERROR - GSFLAT: Image header parameter "//\
                    l_key_exptime//" not found", l_logfile, verbose+)
                nbad = nbad+1
            }
            if (nim == 1)
                firstexp = real(imgets.value)
            else {
                expdiff = (real(imgets.value)-firstexp)/firstexp
                if (expdiff > 0.1)
                    printlog ("WARNING - GSFLAT: "//img//" - exposure time \
                        different by more than 10%", l_logfile, verbose+)
            }


            #Check that images in list all have the same FILTER1 and
            #GRATING ID's and GRATING tilts and the same central wavelength.
            hselect (img//"[0]", "FILTER1,FILTER2,GRATING,GRTILT", yes) | \
                scan (filtn1,filtn2, gratn, grtilt)

            if (substr(filtn1,1,4) == "open") {
                filtn1 = "open"
            }

            if (substr(filtn2,1,4) == "open") {
                filtn2 = "open"
            }

            if (nim == 1) {
                filt1 = filtn1
                filt2 = filtn2
                grat1 = gratn
                grtilt1 = grtilt
            }
            if (filtn1 != filt1) {
                printlog ("ERROR - GSFLAT: Input images do not all have the \
                    same filter1 ID.", l_logfile, verbose+)
                nbad = nbad+1
            }
            if (filtn2 != filt2) {
                printlog ("ERROR - GSFLAT: Input images do not all have the \
                    same filter2 ID.", l_logfile, verbose+)
                nbad = nbad+1
            }
            if (gratn != grat1){
                printlog ("ERROR - GSFLAT: Input images do not all have the \
                    same grating ID.", l_logfile, verbose+)
                nbad = nbad+1
            }
            if (grtilt != grtilt1){
                printlog ("ERROR - GSFLAT: Input images do not all have the \
                    same grating tilt.", l_logfile, verbose+)
                nbad = nbad+1
            }

        } #End of while loop over images
        scanfile = ""

    } #End of if-loop for "nbad=0 so far"


    # Exit if problems found
    if (nbad > 0) {
        printlog ("ERROR - GSFLAT: "//nbad//" fatal errors found.",
            l_logfile, verbose+)
        goto error
    }

    #END OF BASIC CHECKS ...
    #Make sure the input names has .fits on them

    delete (temp1, verify-, >& "dev$null")
    for (n_i=1; n_i<=nim; n_i+=1)
        files (infile[n_i], >> temp1)

    scanfile = temp1

    #NOW, call GSREDUCE to overscan, bias-subtract, dark-subtract and QE
    # correct if requested (most likely)

    reducing = no
    if (l_fl_over || l_fl_bias || l_fl_dark || l_fl_trim || l_fl_qecorr) {
        reducing = yes
        printlog (" ", l_logfile, l_verbose)
        printlog ("GSFLAT: Calling GSREDUCE.", l_logfile, l_verbose)
        printlog (" ", l_logfile, l_verbose)
        gsreduce ("@"//temp1, outimages="", outpref="gs", logfile=l_logfile,
            verbose=l_verbose, fl_gmosaic=l_fl_gmosaic, fl_fixpix=l_fl_fixpix,
            fl_over=l_fl_over, fl_bias=l_fl_bias, fl_dark=l_fl_dark, fl_flat-,
            fl_gsappwave=l_fl_gsappwave, fl_cut-, bias=l_bias, dark=l_dark, \
            flatim="", fl_title=no, fl_fulldq=l_fl_fulldq, dqthresh=l_dqthresh,\
            bpmfile=l_bpmfile, key_exptime=l_key_exptime,
            key_biassec=l_key_biassec, key_datasec=l_key_datasec,
            fl_vardq=l_fl_vardq, sci_ext=l_sci_ext, var_ext=l_var_ext,
            dq_ext=l_dq_ext, key_mdf=l_key_mdf, mdffile=l_mdffile,
            mdfdir=l_mdfdir, bpm=l_bpm, key_ron=l_key_ron, key_gain=l_key_gain,
            ron=l_ron, gain=l_gain, sat=l_sat, ovs_flinter=l_ovs_flinter,
            ovs_med=l_ovs_med, ovs_func=l_ovs_func, ovs_order=l_ovs_order,
            ovs_lowr=l_ovs_lowr, ovs_highr=l_ovs_highr, ovs_niter=l_ovs_niter,
            fl_trim=l_fl_trim, gaindb=l_gaindb, rawpath=l_rawpath,
            gratingdb=l_gratingdb, filterdb=l_filterdb, gradimage="",
            refimage="", nbiascontam=l_nbiascontam, biasrows=l_biasrows,
            fl_qecorr=l_fl_qecorr, qe_refim=l_qe_refim, \
            fl_keep_qeim=l_fl_keep_qeim, qe_corrpref=l_qe_corrpref,
            qe_corrimages=l_qe_corrimages, qe_data=l_qecorr_data, \
            qe_datadir=l_qecorr_datadir)

        if (gsreduce.status != 0) {
            printlog ("ERROR - there was a problem running GSREDUCE.  \
                Stopping now.", l_logfile, verbose+)
            nbad = nbad+1
            goto error
        }

        printlog ("GSFLAT: Returned from GSREDUCE.", l_logfile, l_verbose)
        for (n_i=1; n_i<=nim; n_i+=1) {
            # This removes any possible diectories or environmental varibales
            # from the filename
            fparse (infile[n_i])
            files("gs"//fparse.root//fparse.extension, >> temp2)
        }
        scanfile = temp2
    }

    #Start Process of combining flats
    #i loop is over the # of images

    n_i = 0
    while (fscan(scanfile,img) !=EOF) {

        n_i = n_i+1

        printlog ("GSFLAT: Image #"//n_i//" - Working on : "//img,
            l_logfile, l_verbose)

        #check gain and readnoise values - input parameters vs. header
        #keywords/values.

        if ( l_key_gain != "" ) {
            keyfound = ""
            if ( l_fl_detec )
                hselect (img//"["//l_sci_ext//",1]", l_key_gain, yes) | \
                    scan (keyfound)
            else
                hselect (img//"[0]", l_key_gain, yes) | scan (keyfound)

            if (keyfound == "") {
                printlog ("WARNING - GSFLAT: keyword "//l_key_gain//" not \
                    found in "//img, l_logfile, l_verbose)
                printlog ("                  Using GAIN = "//str(l_gain),
                    l_logfile, verbose+)
                l_sgain = l_gain
            } else
                l_sgain = real(keyfound)
        } else
            l_sgain = l_gain

        if ( l_key_ron != "" ) {
            keyfound = ""
            if ( l_fl_detec )
                hselect (img//"["//l_sci_ext//",1]", l_key_ron, yes) | \
                    scan (keyfound)
            else
                hselect (img//"[0]", l_key_ron, yes) | scan (keyfound)

            if (keyfound == "") {
                printlog ("WARNING - GSFLAT: keyword "//l_key_ron//" not \
                    found in "//img, l_logfile, l_verbose)
                printlog ("                  Using RON = "//str(l_ron),
                    l_logfile, verbose+)
                l_sron = l_ron
            } else
                l_sron = real(keyfound)
        } else
            l_sron = l_ron

        #Make list for imcombine of science extensions
        print (img, >> scilist)

    } #end of while loop

    #Copy MDF file of the first image to disk for later use
    #Note that the MDFs of the input images should all be identical
    #First check that the MDF exists
    fields (scilist, 1, lines="1") | scan(img)
    tinfo (img//"[MDF]", ttout-, >& "dev$null")
    if (tinfo.tbltype!="fits") {
        printlog ("ERROR - GSFLAT: MDF file does not exist", l_logfile, \
            verbose+)
        printlog ("                If GPREPARE previously called, \
            please check if fl_addmdf=yes.", l_logfile, verbose+)
        nbad = nbad+1
        goto error
    }
    tcopy (img//"[MDF]", mdf, verbose=no, >& "dev$null")

    #Now, combine images
    #special cases for low numbers of images

    printlog (" ", l_logfile, l_verbose)

    if (nim == 1) {
        printlog ("WARNING - GSFLAT: only one image.", l_logfile, verbose+)
        copy (img, l_combflat//".fits", verbose-)
        if (l_fl_vardq)
            combsig = infile[1]//"["//l_var_ext//"]"

    } else {
        if (nim <= 5){
            if (nim == 2){
                printlog ("WARNING - GSFLAT: only combining two images, \
                    turning off rejection", l_logfile, verbose+)
                printlog ("                  and setting combine type to \
                    average.", l_logfile, verbose+)
                l_reject = "none"
                l_combine = "average"
            } else
                printlog ("GSFLAT: combining five or less images.",
                    l_logfile, verbose+)
        }

        #gemcombine takes care of fl_vardq setting
        gemcombine ("@"//scilist, l_combflat, combine=l_combine,
            reject=l_reject, offsets="none", masktype=l_masktype,
            maskvalue=l_maskvalue, bpmfile="", scale=l_scale,
            zero=l_zero, nrejfile="", weight=l_weight, lsigma=l_lsigma,
            hsigma=l_hsigma, statsec=l_statsec, expname=l_key_exptime,
            lthreshold=l_lthreshold, hthreshold=l_hthreshold, nlow=l_nlow,
            nhigh=l_nhigh, nkeep=l_nkeep, mclip=l_mclip, key_ron=l_key_ron,
            key_gain=l_key_gain, ron=l_sron, gain=l_sgain, snoise=l_snoise,
            sigscale=l_sigscale, pclip=l_pclip, grow=l_grow, logfile=l_logfile,
            fl_vardq=l_fl_vardq)

        if (l_fl_vardq)
            combsig = l_combflat//"["//l_var_ext//"]"

        printlog ("GSFLAT: GEMCOMBINE complete.  Output in file : "//\
            l_combflat, l_logfile, l_verbose)
    } # End of more than one image

    delete (scilist, verify-, >& "dev$null")

    printlog (" ", l_logfile, l_verbose)

    # Get detector binning
    imgets (l_combflat//"["//l_sci_ext//",1]","CCDSUM")
    print (imgets.value) | scan (xbin,ybin)

    # Determine the nominal chip gap size, last unbinned equivalent column
    # after mosaicking and any bright columns to exclude at the far ends of
    # the CCDs when fitting. The sizes are from gmosaic, where they are also
    # hard wired, and are adjusted oddly for binning since that's what gmosaic
    # does & therefore what's needed to match chipgaps.dat (which is in
    # unbinned units, like DETSEC).
    nbcols = 0
    if (ishamamatsu) {
        if (instrument == "GMOS-S") {
            if (xbin == 1) {
                gapsize = 61
                lastcol = 6266
            } else {
                gapsize = 60
                lastcol = 6264
            }
        } else { # GMOS-N    ## KL need correct values.
            if (xbin == 1) {
                gapsize = 67
                lastcol = 6266
            } else {
                gapsize = 66
                lastcol = 6264
            }
	    }
        # Since chipgaps.dat does not account for bright columns at the far
        # ends of the Hamamatsu CCDs (only those next to the gaps), take half
        # of the avg "extra" columns around each gap, rounded down, as the
        # number to exclude at the ends (can be overridden interactively).
        if (l_bpmfile != "")
            nbcols = (x12-x11-gapsize+x22-x21-gapsize+2)/4
    } else {
        # These values don't get used as long as sample is set to "*" for EEV
        # data with fl_detec+ but are here if needed in future.
        if (xbin == 1) {
            gapsize = 37
            lastcol = 6218
        } else {
            gapsize = 36
            lastcol = 6216
	    }
    }

    # First and last good columns in the full mosaic, without bright cols:
    fgcol = nbcols + 1
    lgcol = lastcol - nbcols

    ## Branch point for detector by detector or slit by slit fitting
    if (l_fl_detec) { ## fitting detector by detector

        ##l_combflat is a MEF file with 3, 6 or 12 science extensions
        ##At this point, we aren't going to do anything with potential
        ##VAR/DQ planes
        printlog (" ", l_logfile, l_verbose)
        printlog ("GSFLAT: Begin fitting response functions for detectors",
            l_logfile, l_verbose)
        printlog (" ", l_logfile, l_verbose)
        copy (l_combflat//".fits", tmpflat//".fits")
        imgets (tmpflat//".fits[0]", "NSCIEXT")
        nsciext = int(imgets.value)
        imgets (tmpflat//".fits[0]", "OBSMODE")
        print(str(imgets.value)) | scanf("%s", obsmode)
        if (obsmode != "MOS") obsmode="LONGSLIT"

        if (nsciext==6 && (l_rdorder[4] < 1 || l_rdorder[5] < 1\
             || l_rdorder[6] < 1)) {
            printlog ("WARNING- GSFLAT: order of fitting function must be \
                greater than zero", l_logfile, verbose+)
            printlog ("WARNING- GSFLAT: setting order of fitting function \
                to "//str(l_rdorder[1])//" for all amps.", l_logfile, verbose+)
            l_rdorder[2] = l_rdorder[1]
            l_rdorder[3] = l_rdorder[1]
            l_rdorder[4] = l_rdorder[1]
            l_rdorder[5] = l_rdorder[1]
            l_rdorder[6] = l_rdorder[1]
        }

        if (nsciext==12 && (l_rdorder[4] < 1 || l_rdorder[5] < 1 \
             || l_rdorder[6] < 1 || l_rdorder[7] < 1 || l_rdorder[8] < 1 \
             || l_rdorder[9] < 1 || l_rdorder[10] < 1 || l_rdorder[11] < 1 \
             || l_rdorder[12] < 1)) {
            printlog ("WARNING- GSFLAT: order of fitting function must be \
                greater than zero", l_logfile, verbose+)
            printlog ("WARNING- GSFLAT: setting order of fitting function \
                to "//str(l_rdorder[1])//" for all amps.", l_logfile, verbose+)
            l_rdorder[2] = l_rdorder[1]
            l_rdorder[3] = l_rdorder[1]
            l_rdorder[4] = l_rdorder[1]
            l_rdorder[5] = l_rdorder[1]
            l_rdorder[6] = l_rdorder[1]
            l_rdorder[7] = l_rdorder[1]
            l_rdorder[8] = l_rdorder[1]
            l_rdorder[9] = l_rdorder[1]
            l_rdorder[10] = l_rdorder[1]
            l_rdorder[11] = l_rdorder[1]
            l_rdorder[12] = l_rdorder[1]
        }

        for (n_ccd=1; n_ccd<=nsciext; n_ccd+=1) {
            printlog ("GSFLAT: fitting response functions for CCD"//str(n_ccd),
                l_logfile, l_verbose)

            # Get the range of columns corresponding to this amp/CCD in the
            # full, unbinned detector mosaic and compare that with chipgaps.dat
            # to exclude any artifically bright columns bordering the chip gaps
            # (for Hamamatsu detectors) from the fitting sample range. This is
            # disabled for EEV chips (largely to avoid combing through many
            # small regression test differences unnecessarily) but the logic is
            # in place to remove that condition if useful.
            if (ishamamatsu && l_bpmfile != "") {
                imgets (tmpflat//"["//l_sci_ext//", "//str(n_ccd)//"]", \
                    "DETSEC")
                print (imgets.value) | scanf("[%d:%d,%d:%d]", x1, x2, y1, y2)
                # Add chip gap columns to the contiguous DETSEC values to
                # convert to columns in gmosaic output (ignoring sub-pixel
                # adjustments):
                gapoff = int((x1-1)/2048.) * gapsize
                x1 += gapoff
                x2 += gapoff
                xref = x1
                # If either of the limits falls in a masked-out "gap" region
                # (the bright region to one side of the gap), adjust it to the
                # column before or after the appropriate end of the gap. Do
                # likewise for bright columns at the far ends of the array.
                if (x1 < fgcol)
                    x1 = fgcol
                else if (x1 >= x11 && x1 <= x12)
                    x1 = x12 + 1
                else if (x1 >= x21 && x1 <= x22)
                    x1 = x22 + 1
                if (x2 > lgcol)
                    x2 = lgcol
                else if (x2 >= x11 && x2 <= x12)
                    x2 = x11 - 1
                else if (x2 >= x21 && x2 <= x22)
                    x2 = x21 - 1
                # Convert the range of columns to the possibly-binned pixel
                # co-ordinates for the current amp/CCD, rounding to the nearest
                # integer binned pixel:
                xr = real(x1-xref) / xbin + 1
                x1 = int(xr)
                if (frac(xr) > 0.) x1 += 1  # round up unless already int pix
                x2 = int(real(x2-xref) / xbin + 1)  # round down
                # Construct the sample range specification used by the fit:
                l_sample = str(x1)//":"//str(x2)
            } else {
                l_sample = "*"
            }
            printlog ("        sample in fitting "//l_sample, \
                l_logfile, l_verbose)

            # Fit the continuum row-by-row or all together & normalize:
            if (l_fl_seprows) {
                fit1d (tmpflat//"["//l_sci_ext//","//str(n_ccd)//"]",
                    l_combflat//"["//l_sci_ext//","//str(n_ccd)//",overwrite]",
                    type="ratio", axis=1, interactive=l_rinteractive,
                    sample=l_sample, naverage=1, function=l_rfunction,
                    order=l_rdorder[n_ccd], low_reject=l_rlowrej,
                    high_reject=l_rhighrej, niterate=l_niterate, grow=1)

            } else {
                # Divide all the rows by a continuum fit to the average,
                # ignoring any incomplete illumination at the edge rows
                # by rejecting the faintest pixels. Normalize everything to
                # the middle long slit section so the flux is preserved there.
                flatfit = mktemp("tmpfit")
                rowavg = mktemp("tmprowavg")
                if (obsmode=="LONGSLIT") {
                    hselect (tmpflat//"["//l_sci_ext//","//str(n_ccd)//"]",
                        "i_naxis2", yes) | scan(ymax)
                    ycen = nint(ymax/ybin)
                    ymax -= nint(90/ybin) # bridges approx. 45 pix wide
                    ymax = nint(0.85*ymax/6.) # 85% of half length of one sec
                    ymin = ycen - ymax
                    ymax = ycen + ymax
                    slitsec="[*,"//ymin//":"//ymax//"]"
                }
                else {
                    slitsec="[*,*]"
                    print("WARNING: fl_detec+ with fl_seprows- is likely to "//
                        "produce a meaningless")
                    print("         normalization for masks other than the "//
                        "long slits\n")
                }
                imcombine(tmpflat//"["//l_sci_ext//","//str(n_ccd)//"]"//\
                    slitsec,
                    rowavg, headers="", bpmasks="",rejmasks="",nrejmasks="",
                    expmasks="", sigmas="", logfile="dev$null",
                    combine="average", reject="avsigclip", project+,
                    outtype="real", outlimits="", offsets="none",
                    masktype="none", maskvalue=0.,blank=1.,scale="",
                    zero="none", weight="none", statsec="",expname="",
                    lthresh=INDEF, hthresh=INDEF, nlow=nint(10./ybin),
                    nhigh=1, nkeep=1, mclip=yes, lsigma=5., hsigma=5., grow=0.)
                fit1d (rowavg, flatfit, type="fit", axis=1,
                    interactive=l_rinteractive, sample=l_sample,
                    naverage=1, function=l_rfunction, order=l_rdorder[n_ccd],
                    low_reject=l_rlowrej, high_reject=l_rhighrej,
                    niterate=l_niterate, grow=1.)
                imarith(tmpflat//"["//l_sci_ext//","//str(n_ccd)//"]", "/",
                    flatfit,
                    l_combflat//"["//l_sci_ext//","//str(n_ccd)//",overwrite]",
                    title="",verbose-,noact-)
                imdelete (rowavg, verify-, >& "dev$null")
                imdelete (flatfit, verify-, >& "dev$null")
            }
            if (l_fl_vardq) {
                imexpr ("a*b**2/c**2",l_combflat//"["//l_var_ext//"," \
                    //str(n_ccd)//",overwrite]",
                    tmpflat//"["//l_var_ext//","//str(n_ccd)//"]",
                    l_combflat//"["//l_sci_ext//","//str(n_ccd)//"]",
                    tmpflat//"["//l_sci_ext//","//str(n_ccd)//"]",verbose-)
                imcopy (tmpflat//"["//l_dq_ext//","//str(n_ccd)//"]",
                    l_combflat//"["//l_dq_ext//","//str(n_ccd)//",overwrite]",
                    verbose-)
            }
        }

        # Set unilluminated pixels to 1 so they won't make science data noisy.
        if (minval != INDEF)
            imreplace(l_combflat//"["//l_sci_ext//","//str(n_ccd)//"]",
                1.0, lower=INDEF, upper=minval, radius=0.0)

        suf = substr (l_specflat, strlen(l_specflat)-3, strlen(l_specflat))
        if (suf!="fits")
            l_specflat = l_specflat//".fits"
        rename (l_combflat//".fits", l_specflat, field="all")
        rename (tmpflat//".fits", l_combflat//".fits", field="all")

    } else { ## fitting slit by slit

        ##l_combflat is a MEF file with a single mosaiced science extension.

        #VAR/DQ made outside the loop
        #The combined flat is l_combflat.
        #For LONGSLIT mode the MDF will have 3 lines (2 bridges, three spectra
        #sections) We want to keep these in the same SCI extensions however
        #even though we will do some of the calibrations separately for each
        #section.  GSCUT will now do this rather than this task, but some
        #things we still need to take care of here.

        tinfo (mdf, tbltype="fits", subtype="binary", ttout=no)
        nslits = tinfo.nrows
        printlog ("GSFLAT : Found "//str(nslits)//" Slits in MDF.",
            l_logfile, l_verbose)

        # Determine whether MOS or long slit, to select the appropriate
        # region for normalization with fl_seprows-.
        imgets (l_combflat//"[0]", "OBSMODE")
        print(str(imgets.value)) | scanf("%s", obsmode)
        if (obsmode != "MOS") obsmode="LONGSLIT"

        # make dummy images for inserting response sections for both science
        # and variance extensions (if needed)

        printlog (" ", l_logfile, l_verbose)
        combflatsciext = l_combflat//"["//l_sci_ext//",1]"
        imarith (combflatsciext,"/",combflatsciext, response, verbose-)

        #Call gscut to get the images sections for each spectrum in the MDF
        fl_delgrad = no
        l_gradimage = ""
        if (l_fl_usegrad)
            l_gradimage = l_combflat

        printlog ("GSFLAT: Calling GSCUT to determine image sections for \
            spectra.", l_logfile, l_verbose)
        scisec = mktemp("tmpscisec")
        gscut (l_combflat, outimag=scisec, fl_update+, secfile=specsecfile,
            logfile=l_logfile, fl_vard=l_fl_vardq, fl_oversize=l_fl_oversize,
            gratingdb=l_gratingdb, filterdb=l_filterdb, verbose=l_verbose,
            xoffset=l_xoffset, yoffset=l_yoffset, yadd=l_yadd,
            w2=l_upper_wave_limit, refimage=l_refimage, gradimage=l_gradimage)

        if (!access(specsecfile) ) {
            printlog ("ERROR - GSFLAT: GSCUT failed in some way. ",
                l_logfile, verbose+)
            nbad = nbad+1
            goto error
        }
        printlog("GSFLAT: Returned from GSCUT.", l_logfile, l_verbose)

        # Set the sample to avoid the chip gaps, if possible
        if (l_bpmfile != "") {
            # Changed translit command from:
            # translit ("STDIN", '"', "", delete+) to
            # translit ("STDIN", '\"', delete+) for formatting puposes.
            # You don't ned to replace the double quotes with an empty string
            # if you are deleting them. The backslash will never be in CCDSUM
            # and is only being used as an escape character here - MS
            hselect (l_combflat//"["//l_sci_ext//",1]", "CCDSUM", yes) | \
                translit ("STDIN", '\"', delete+) | scan(Xbin)
            hselect (l_combflat//"["//l_sci_ext//",1]", "i_naxis1", yes) | \
                scan(Xmax)
            x11 = int((x11-1)/Xbin+0.5)
            x12 = int((x12-1)/Xbin+0.5)
            x21 = int((x21-1)/Xbin+0.5)
            x22 = int((x22-1)/Xbin+0.5)
            l_sample = "1:"//str(x11-1)//","//str(x12+1)//":" \
                //str(x21-1)//","//str(x22+1)//":"//Xmax
        } else {
            l_sample = "*"
        }
        printlog (" ", l_logfile, l_verbose)
        printlog ("GSFLAT: Begin fitting response functions for slit(s)",
            l_logfile, l_verbose)
        printlog ("        Sample in fitting "//l_sample,
            l_logfile, l_verbose)
        printlog (" ", l_logfile, l_verbose)

        #Now we process each slitlet
        count (specsecfile) | scan (nslits)

        scanfile = specsecfile
        n_i = 1
        while(fscan(scanfile,specsec) !=EOF) {
            # Temporary FITS files re-used within this loop must have
            # a different name at each iteration (FITS Kernel cache is
            # based on file names)

            flatfit = mktemp("tmpfit")
            sciflat = mktemp("tmpflat")
            rowavg = mktemp("tmprowavg")
            inscoo = mktemp("tmpinscoo")

            # Get x1,y1 from MDF, MDFROW from gscut image
            # Parse the specsec we already have instead of reading SECY1
            # from the MDF as previously (was giving the wrong offset for
            # the long slit, where gscut doesn't record it properly!)
            print(specsec) | scanf("[%d:%d,%d:%d]", x1, x2, y1, y2)
            print (x1," ",y1, > inscoo)

            if (l_verbose)
                printlog ("GSFLAT: Fitting response for slit#"//str(n_i)//\
                    " - SpecSec = "//specsec, l_logfile, l_verbose)

            if (l_fl_seprows) {
                fit1d (scisec//"["//l_sci_ext//","//str(n_i)//"]", sciflat,
                    type="ratio", axis=1, interactive=l_rinteractive,
                    sample=l_sample, naverage=1, function=l_rfunction,
                    order=l_rorder, low_reject=l_rlowrej,
                    high_reject=l_rhighrej, niterate=l_niterate, grow=0.)
            } else {
                # Divide all the rows by a continuum fit to the average,
                # ignoring any incomplete illumination at the edge rows
                # by rejecting the faintest pixels. For the long slit,
                # normalize everything to the middle slit section so the
                # flux is preserved there.
                if (obsmode=="LONGSLIT") {
                    hselect (scisec//"["//l_sci_ext//","//str(n_i)//"]",
                      "i_naxis2", yes) | scan(ymax)
                    ycen = nint(ymax/ybin)
                    ymax -= nint(90/ybin) # bridges approx. 45 pix wide
                    ymax = nint(0.85*ymax/6.) # 85% of half length of one sec
                    ymin = ycen - ymax
                    ymax = ycen + ymax
                    slitsec="[*,"//ymin//":"//ymax//"]"
                }
                else slitsec="[*,*]"
                imcombine(scisec//"["//l_sci_ext//","//str(n_i)//"]"//slitsec,
                    rowavg, headers="", bpmasks="",rejmasks="",nrejmasks="",
                    expmasks="",sigmas="", logfile="dev$null",
                    combine="average",reject="minmax", project+,
                    outtype="real",outlimits="",offsets="none",
                    masktype="none",maskvalue=0.,blank=1.,scale="",zero="none",
                    weight="none",statsec="",expname="",lthresh=INDEF,
                    hthresh=INDEF,nlow=nint(10./ybin),nhigh=1,nkeep=1,
                    mclip=yes,lsigma=5., hsigma=5.,grow=0.)
                fit1d (rowavg, flatfit, type="fit", axis=1,
                    interactive=l_rinteractive, sample=l_sample,
                    naverage=1, function=l_rfunction, order=l_rorder,
                    low_reject=l_rlowrej, high_reject=l_rhighrej,
                    niterate=l_niterate, grow=0.)
                imarith(scisec//"["//l_sci_ext//","//str(n_i)//"]", "/",
                    flatfit,sciflat,title="",verbose-,noact-)
            }
            if (l_slitcorr)
                imarith (sciflat,"*",
                    l_slitfunc//"["//l_sci_ext//"]"//specsec,
                    sciflat, verbose-)
            # Set unilluminated pixels to 1 so they won't make data noisy.
            if (minval != INDEF)
                imreplace(sciflat, 1.0, lower=INDEF, upper=minval, radius=0.0)

            # insert sections into larger images
            iminsert (response, sciflat, response, "replace", coordfile=inscoo,
                offset1=0, offset2=0, xcol="c1", ycol="c2")
            imdelete (sciflat, verify-, >& "dev$null")
            imdelete (rowavg, verify-, >& "dev$null")
            imdelete (flatfit, verify-, >& "dev$null")
            delete (inscoo, verify-)
            n_i = n_i+1

            # Give the user a chance to fit the rest non-interactively
            if (l_rinteractive && nslits>1) {
                l_fl_answer = fl_answer
                l_rinteractive = l_fl_answer
            }


        } #end of while loop over slits
        scanfile = ""
        printlog (" ", l_logfile, l_verbose)
        printlog ("GSFLAT: Finished Loop over slits. Now packing output \
            MEF...", l_logfile, l_verbose)

        #Copy the PHU from the first of the "reduced" flatfields
        suf = substr (l_specflat, strlen(l_specflat)-3, strlen(l_specflat))
        if (suf!="fits")
            l_specflat = l_specflat//".fits"

        wmef (input=mdf, output=l_specflat, extnames="MDF",
            phu=l_combflat//".fits", verbose-, >& "dev$null")
        if (wmef.status != 0) {
            printlog ("ERROR - GSFLAT: problem writing output MEF.",
                l_logfile, verbose+)
            nbad = nbad+1
            goto error
        }

        #Now the science extension/flatfield
        fxinsert (response//".fits", l_specflat//"[1]", groups="", verbose=yes,
            >& "dev$null")

        gemhedit (l_specflat//"[2]", "EXTNAME", l_sci_ext, "Extension name", \
            delete-)
        gemhedit (l_specflat//"[2]", "EXTVER", 1, "Extension number", \
            delete-)

        #The VAR/DQ planes if requested ...
        if (l_fl_vardq) {
            imexpr ("a*b**2/c**2",l_specflat//"["//l_var_ext//",1,append]",
                l_combflat//"["//l_var_ext//",1]",
                l_specflat//"["//l_sci_ext//",1]",
                l_combflat//"["//l_sci_ext//",1]", verbose-)
            imcopy (l_combflat//"["//l_dq_ext//",1]",
                l_specflat//"["//l_dq_ext//",1,append]", verbose-)
        }
    } ## end of slit by slit fitting section

    # Delete unnecessary DATAMAX & DATAMIN header keywords fit1d has
    #been adding these keywords and setting them to 0.00000E0 for unknown
    #reasons at this time, will need to eventually look further into this.
    hselect(l_specflat//"[0]","NSCIEXT",yes) | scan(nsciext)
    for (n_ccd=1; n_ccd<=nsciext; n_ccd+=1) {
        keyfound = ""
        hselect (l_specflat//"["//l_sci_ext//","//str(n_ccd)//"]", \
            "DATAMAX", yes) | scan (keyfound)
        if (keyfound != "") {
            gemhedit (l_specflat//"["//l_sci_ext//","//str(n_ccd)//"]", \
                "DATAMAX", "", "", delete+)
            gemhedit (l_specflat//"["//l_sci_ext//","//str(n_ccd)//"]", \
                "DATAMIN", "", "", delete+)
        }
    }

    if (l_fl_vardq) {
        for (n_ccd=1; n_ccd<=nsciext; n_ccd+=1) {
            keyfound = ""
            hselect (l_specflat//"["//l_var_ext//","//str(n_ccd)//"]", \
                "DATAMAX", yes) | scan (keyfound)
            if (keyfound != "") {
                gemhedit (l_specflat//"["//l_var_ext//","//str(n_ccd)//"]", \
                    "DATAMAX", "", "", delete+)
                gemhedit (l_specflat//"["//l_var_ext//","//str(n_ccd)//"]", \
                    "DATAMIN", "", "", delete+)
            }
            keyfound = ""
            hselect (l_specflat//"["//l_dq_ext//","//str(n_ccd)//"]", \
                "DATAMAX", yes) | scan (keyfound)
            if (keyfound != "") {
                gemhedit (l_specflat//"["//l_dq_ext//","//str(n_ccd)//"]", \
                    "DATAMAX", "", "", delete+)
                gemhedit (l_specflat//"["//l_dq_ext//","//str(n_ccd)//"]", \
                    "DATAMIN", "", "", delete+)
            }
        }
    }

    # Doubled flats for nod and shuffle
    if(l_fl_double) {
        if (l_nshuffle == 0) {
            printlog ("GSFLAT: Warning - Must provide a shuffle distance", \
                l_logfile, verbose+)
            printlog ("         to shift Nod-and-shuffle flat.", l_logfile, \
                verbose+)
            printlog ("         nshuffle can usually be found as NODPIX", \
                l_logfile, verbose+)
            printlog ("         header keyword in nod-and-shuffle science \
                data", l_logfile, verbose+)
        }

        printlog("GSFLAT: Making double flats for nod-and-shuffle science", \
            l_logfile, l_verbose)
        hselect(l_combflat//"[0]","NSCIEXT",yes) | scan(nsciext)
        imgets (l_combflat//"["//l_sci_ext//",1]","CCDSUM")
        print (imgets.value) | scan (xbin,ybin)
        for(n_i=1;n_i<=nsciext;n_i+=1) {
            tmpratsh=mktemp("tmpratsh")
            tmpspecflat=mktemp("tmpspecflat")
            imshift(l_combflat//"["//l_sci_ext//","//n_i//"]", tmpratsh, \
                xshift=0, yshift=(-l_nshuffle/ybin), interp_type="nearest", \
                boundary_typ="nearest", constant=0.)
            imshift(l_specflat//"["//l_sci_ext//","//n_i//"]", tmpspecflat, \
                xshift=0, yshift=(-l_nshuffle/ybin), interp_type="nearest", \
                boundary_typ="nearest", constant=0.)
            imexpr("c>d ? a : b",l_specflat//"["//l_sci_ext//","//n_i//", \
                overwrite]", l_specflat//"["//l_sci_ext//","//n_i//"]", \
                tmpspecflat, l_combflat//"["//l_sci_ext//","//n_i//"]", \
                tmpratsh, verbose-)
            imdelete(tmpratsh//","//tmpspecflat,verify-)
        }
    }

    #Update the PHU
    gemdate ()
    gemhedit (l_specflat//"[0]", "NCOMBINE", nim,
        "Number of images in IMCOMBINE", delete-)
    gemhedit (l_specflat//"[0]", "GSFLAT", gemdate.outdate,
        "UT Time stamp for GSFLAT", delete-)
    gemhedit (l_specflat//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Time of last modification with GEMINI", delete-)
    if( l_slitcorr)
        gemhedit (l_specflat//"[0]", "GSLITCOR", gemdate.outdate,
            "UT Time stamp for Slit Function Correction", delete-)
    if (l_fl_detec)
        gemhedit (l_specflat//"[0]", "GSFMODE", "DETECTOR",
            "GSFLAT fitting mode", delete-)
    else
        gemhedit (l_specflat//"[0]", "GSFMODE", "SLIT", \
            "GSFLAT fitting mode", delete-)
    gemhedit (l_specflat//"[0]", "GSFFCT", l_rfunction,
        "GSFLAT fitting function", delete-)
    gemhedit (l_specflat//"[0]", "GSFORDER", l_order, "GSFLAT fitting order", \
        delete-)

    gemhedit (l_specflat//"[0]", "GAINORIG", l_sgain, "Original Gain value", \
        delete-)
    gemhedit (l_specflat//"[0]", "RONORIG", l_sron, \
        "Original Read-noise value", delete-)
    if(l_fl_double) {
        gemhedit(l_specflat//"[0]","GSFDOUB", l_fl_double, \
            "Double flatfield", delete-)
        gemhedit(l_specflat//"[0]","GSFNSHUF", l_nshuffle, \
            "Number of shuffle pixels", delete-)
    }

    gaineff = real(nim)*l_sgain
    roneff = l_sron/sqrt(real(nim))
    if (l_key_gain == "")
        l_key_gain="GAIN"
    if (l_key_ron == "")
        l_key_ron="RON"
    gemhedit (l_specflat//"[0]", l_key_gain, gaineff, "Rescaled Gain value", \
        delete-)
    gemhedit (l_specflat//"[0]", l_key_ron, roneff, \
        "Rescaled Read-noise value", delete-)

    # clean up
    goto clean

error:
    status = 1
    printlog (" ", l_logfile, verbose+)
    printlog ("ERROR - GSFLAT: Program execution failed with "//\
        str(nbad)//" errors.", l_logfile, verbose+)

clean:
    scanfile = ""
    if (status==0) {
        printlog (" ", l_logfile, verbose+)
        date | scan(sdate)
        printlog ("GSFLAT done "//sdate, l_logfile, verbose+)
    }
    if (!l_redkeep)
        delete (l_combflat//".fits", verify-, >& "dev$null")
    imdelete (scisec, verify-, >& "dev$null")
    imdelete (response, verify-, >& "dev$null")
    delete (mdf//","//temp1//","//temp2//","//specsecfile, verify-, \
        >& "dev$null")

    # close log file
    printlog ("-------------------------------------------------------------\
        -------------------", l_logfile, l_verbose)

end
