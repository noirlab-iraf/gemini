# Copyright(c) 2002-2006 Association of Universities for Research in Astronomy, Inc.

procedure nscombine (inimages) 

# Version Sept 20, 2002 JJ v1.4 release
#         Aug 19, 2003  KL IRAF2.12 - new parameter, addonly, in hedit
#                          Fixed even though the hedit calls are commented out.
#         Oct 29, 2003  KL moved from niri to gnirs

char    inimages    {prompt = "Input images or spectra to shift and combine"}			# OLDP-1-primary-combine-suffix=_add
real    tolerance   {0.5, min = 0., prompt = "Maximum offset for grouping the frames (arcsec)"}	# OLDP-2
char    output      {"", prompt = "File name for output shifted and combined spectra"}		# OLDP-1-output
char    output_suffix   {"_comb", prompt = "Suffix for output shifted and combined image"}		# OLDP-3
char    bpm         {"", prompt = "Name of bad pixel mask."}                # OLDP-1
int     dispaxis    {1, prompt = "Dispersion axis (if not in header)"}      # OLDP-3
real    pixscale    {1., prompt = "Pixel scale in arcsec/pixel (if not in header)"}			# OLDP-3
bool    fl_cross    {no, prompt = "Update WCS from cross-correlation?"}     # OLDP-3
bool    fl_keepshift    {no, prompt = "Keep shifted images?"}               # OLDP-3
bool    fl_shiftint {yes, prompt = "Shift frames by integer pixels?"}       # OLDP-2
char    interptype  {"linear", prompt = "Interpolation type for shifting",  enum = "nearest|linear|poly3|poly5|spline3|sinc|drizzle"}						# OLDP-2
char    boundary    {"nearest", prompt = "Boundary type for shifting", enum = "constant|nearest|reflect|wrap"}									# OLDP-2
real    constant    {0., prompt = "Constant value for boundary extension when shifting"}		# OLDP-2
char    combtype    {"average", prompt = "Combination operation", enum = "average|median"}		# OLDP-2
char    rejtype     {"sigclip", prompt = "Rejection algorithm when combining", enum = "none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip"}					# OLDP-2
char    masktype    {"goodvalue", enum="none|goodvalue", prompt="Mask type"}    # OLDP-3
real    maskvalue   {0., prompt="Mask value"}                               # OLDP-3
char    statsec     {"[*,*]", prompt = "Image section to be used for statistics"}			# OLDP-2
char    scale       {"none", prompt = "Image scaling"}						# OLDP-2
char    zero        {"none", prompt = "Image zeropoint offset"}             # OLDP-2
char    weight      {"none", prompt = "Image weights"}						# OLDP-2
real    lthreshold  {INDEF, prompt = "Lower threshold"}						# OLDP-2
real    hthreshold  {INDEF, prompt = "Upper threshold"}						# OLDP-2
int     nlow        {1, min = 0, prompt = "minmax: Number of low pixels to reject"}			# OLDP-2
int     nhigh       {1, min = 0, prompt = "minmax: Number of high pixels to reject"}		# OLDP-2
int     nkeep       {0, prompt = "Minimum to keep or maximum to reject"}    # OLDP-2
bool    mclip       {yes, prompt = "Use median in sigma clipping algorithms?"}			# OLDP-2
real    lsigma      {5., prompt = "Lower sigma clipping factor"}            # OLDP-2
real    hsigma      {5., prompt = "Upper sigma clipping factor"}            # OLDP-2
real    ron         {0.0, min = 0., prompt = "Readout noise rms in electrons"}			# OLDP-2
real    gain        {1.0, min = 0.00001, prompt = "Gain in e-/ADU"}         # OLDP-2
char    snoise      {"0.0", prompt = "ccdclip: Sensitivity noise (electrons)"}			# OLDP-2
real    sigscale    {0.1, prompt = "Tolerance for sigma clipping scaling correction"}		# OLDP-2
real    pclip       {-0.5, prompt = "pclip: Percentile clipping parameter"} # OLDP-2
real    grow        {0.0, prompt = "Radius (pixels) for neighbor rejection"}			# OLDP-2
char    nrejfile    {"", prompt = "Name of rejected pixel count image."}    # OLDP-3
bool    fl_vardq    {yes, prompt = "Create output variance and data quality frames?"}		# OLDP-2
bool    fl_inter    {no, prompt = "Measure corss-correlation peak interactively?"}			# OLDP-4
char    logfile     {"", prompt = "Logfile"}                                # OLDP-1
bool    verbose     {yes, prompt = "Verbose output?"}                       # OLDP-2
bool    debug       {no, prompt = "Very verbose output?"}                   # OLDP-2
bool    force       {no, prompt = "Force use with earlier IRAF versions?"}  # OLDP-3
int     status      {0, prompt = "Exit status (0=good)"}                    # OLDP-4
struct  *scanfile   {"", prompt = "Internal use only"} 						# OLDP-4

begin
        
    char    l_inimages = ""
    real    l_tolerance
    char    l_output = ""
    char    l_output_suffix = ""
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
    bool    l_verbose
    bool    l_debug
    bool    l_force


    int     junk


    junk = fscan (  inimages, l_inimages)
    l_tolerance =   tolerance
    junk = fscan (  output, l_output)
    junk = fscan (  output_suffix, l_output_suffix)
    junk = fscan (  bpm, l_bpm)
    l_dispaxis =    dispaxis
    l_pixscale =    pixscale
    l_fl_cross =    fl_cross
    l_fl_keepshift =    fl_keepshift
    l_fl_shiftint =     fl_shiftint
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
    l_verbose =     verbose
    l_debug =       debug
    l_force =       force


    status = 1
    cache ("nschelper")

    nschelper (inimages, tolerance = l_tolerance, outimages = l_output, \
        output_suffix = l_output_suffix, fl_shift+, \
        bpm = l_bpm, dispaxis = l_dispaxis, pixscale = l_pixscale, \
        fl_cross = l_fl_cross, \
        fl_keepshift = l_fl_keepshift, fl_shiftint = l_fl_shiftint, \
        interptype = l_interptype, boundary = l_boundary, \
        constant = l_constant, combtype = l_combtype, \
        rejtype = l_rejtype, masktype=l_masktype, maskvalue=l_maskvalue, \
        statsec = l_statsec, scale = l_scale, \
        zero = l_zero, weight = l_weight, lthreshold = l_lthreshold, \
        hthreshold = l_hthreshold, nlow = l_nlow, nhigh = l_nhigh, \
        nkeep = l_nkeep, mclip = l_mclip, lsigma = l_lsigma, \
        hsigma = l_hsigma, ron = l_ron, gain = l_gain, \
        snoise = l_snoise, sigscale = l_sigscale, pclip = l_pclip, \
        grow = l_grow, nrejfile = l_nrejfile, fl_vardq = l_fl_vardq, \
        fl_inter = l_fl_inter, force = l_force, \
        verbose = l_verbose, debug = l_debug, logfile = l_logfile, \
        logname = "NSCOMBINE")

    status = nschelper.status

end
