# Copyright(c) 2001-2006 Association of Universities for Research in Astronomy, Inc.

procedure nsstack (inimages) 

# Version Sept 20, 2002 MT,JJ v1.4 release
#         Aug 19, 2003  KL  IRAF2.12 - new parameter, addonly, in hedit
#                           Fixed even though all hedit calls are commented out
#         Oct 29, 20003 KL  moved from niri to gnirs package

char    inimages    {prompt = "Input images or spectra to combine"}
real    tolerance   {0.5, prompt = "Maximum offset for grouping (arcsec)"}
char    combtype    {"average", prompt = "Combination operation", enum = "average|median"}
char    rejtype     {"sigclip", prompt = "Rejection algorithm", enum = "none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip"}
char    masktype    {"goodvalue", enum="none|goodvalue", prompt="Mask type"}
real    maskvalue   {0., prompt="Mask value"}
char    statsec     {"[*,*]", prompt = "Image section to be used for statistics"}
char    stacksuffix {"_stack", prompt = "Suffix for combined stacks"}
int     dispaxis    {1, prompt = "Dispersion axis (if not in header)"}
real    pixscale    {1., prompt = "Pixel scale in arcsec/pixel (if not in header)"}
char    scale       {"none", prompt = "Image scaling"}
char    zero        {"none", prompt = "Image zeropoint offset"}
char    weight      {"none", prompt = "Image weights"}
real    lthreshold  {INDEF, prompt = "Lower threshold"}
real    hthreshold  {INDEF, prompt = "Upper threshold"}
int     nlow        {1, min = 0, prompt = "minmax: Number of low pixels to reject"}
int     nhigh       {1, min = 0, prompt = "minmax: Number of high pixels to reject"}
int     nkeep       {0, prompt = "Minimum to keep or maximum to reject"}
bool    mclip       {yes, prompt = "Use median in sigma clipping algorithms?"}
real    lsigma      {5., prompt = "Lower sigma clipping factor"}
real    hsigma      {5., prompt = "Upper sigma clipping factor"}
real    ron         {0.0, min = 0., prompt = "Readout noise rms in electrons"}
real    gain        {1.0, min = 0.00001, prompt = "Gain in e-/ADU"}
char    snoise      {"0.0", prompt = "ccdclip: Sensitivity noise (electrons)"}
real    sigscale    {0.1, prompt = "Tolerance for sigma clipping scaling correction"}
real    pclip       {-0.5, prompt = "pclip: Percentile clipping parameter"}
real    grow        {0.0, prompt = "Radius (pixels) for neighbor rejection"}
char    nrejfile    {"", prompt = "Names of rejected pixel count images."}
bool    fl_vardq    {yes, prompt = "Create variance and data quality frames in output images?"}
char    logfile     {"", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose output?"}			
bool    debug       {no, prompt = "Very verbose output?"}			
bool    force       {no, prompt = "Force use with earlier IRAF versions?"}
int     status      {0, prompt = "Exit status (0=good)"}
struct  *scanfile   {"", prompt = "Internal use only"}


begin
        
    char    l_inimages = ""
    real    l_tolerance
    char    l_combtype = ""
    char    l_rejtype = ""
    char    l_masktype = ""
    real    l_maskvalue
    char    l_statsec = ""
    char    l_stacksuffix = ""
    int     l_dispaxis
    real    l_pixscale
    char    l_scale
    char    l_zero	
    char    l_weight	
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
    char    l_logfile = ""
    bool    l_verbose
    bool    l_debug	
    bool    l_force

    # Those headers are not used here, but they are used
    # in NSCHELPER.  Keep the declaration here as documentation.
    #char	l_key_ron = ""
    #char	l_key_gain = ""
    #char	l_key_exptime = ""
    #char	l_key_xoff = ""
    #char	l_key_yoff = ""
    #char	l_sci_ext = ""
    #char	l_var_ext = ""
    #char	l_dq_ext = ""


    int     junk


    junk = fscan (  inimages, l_inimages)
    l_tolerance =   tolerance
    junk = fscan (  combtype, l_combtype)
    junk = fscan (  rejtype, l_rejtype)
    junk = fscan (  masktype, l_masktype)
    l_maskvalue = maskvalue
    l_statsec =     statsec		# may contain spaces
    junk = fscan (  stacksuffix, l_stacksuffix)
    l_dispaxis =    dispaxis
    l_pixscale =    pixscale
    l_scale =       scale
    l_zero =        zero	
    l_weight =      weight	
    l_lthreshold =  lthreshold
    l_hthreshold =  hthreshold
    l_nlow =        nlow	
    l_nhigh	=       nhigh	
    l_nkeep	=       nkeep	
    l_mclip =       mclip
    l_lsigma =      lsigma	
    l_hsigma =      hsigma	
    l_ron =         ron	
    l_gain =        gain
    junk = fscan (  snoise, l_snoise)
    l_sigscale =    sigscale
    l_pclip	=       pclip	
    l_grow =        grow
    junk = fscan (  nrejfile, l_nrejfile)
    l_fl_vardq =    fl_vardq
    junk = fscan (  logfile, l_logfile)
    l_verbose =     verbose	
    l_debug =       debug
    l_force =       force

    status = 1
    cache ("nschelper")

    nschelper (inimages, tolerance = l_tolerance, outimages = "", \
        output_suffix = l_stacksuffix, fl_shift-, fl_cross-, \
        bpm = "", dispaxis = l_dispaxis, pixscale = l_pixscale, \
        fl_keepshift-, fl_shiftint-, interptype = "linear", \
        boundary = "nearest", constant = 0, combtype = l_combtype, \
        rejtype = l_rejtype, masktype=l_masktype, maskvalue=l_maskvalue, \
        statsec = l_statsec, scale = l_scale, \
        zero = l_zero, weight = l_weight, lthreshold = l_lthreshold, \
        hthreshold = l_hthreshold, nlow = l_nlow, nhigh = l_nhigh, \
        nkeep = l_nkeep, mclip = l_mclip, lsigma = l_lsigma, \
        hsigma = l_hsigma, ron = l_ron, gain = l_gain, \
        snoise = l_snoise, sigscale = l_sigscale, pclip = l_pclip, \
        grow = l_grow, nrejfile = l_nrejfile, fl_vardq = l_fl_vardq, \
        fl_inter-, force = l_force, verbose = l_verbose, \
        debug = l_debug, logfile = l_logfile, logname = "NSSTACK")

    status = nschelper.status

end
