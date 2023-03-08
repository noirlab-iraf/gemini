# Copyright(c) 2009-2015 Association of Universities for Research in Astronomy, Inc.

procedure gnsdark(inimages,outdark)

# Derives average GMOS dark images from raw or gprepared images
# Updates GAIN and RDNOISE in the final (combined) image (SCI exts)
# Assumes that DARK frames do not contain an MDF
#
# Version  April 22, 2008  JH 

string  inimages    {prompt="Input GMOS dark images or list"}
string  outdark     {prompt="Output dark (zero level) image"}
string  logfile     {"",prompt="Logfile"}
string  rawpath     {"",prompt="GPREPARE: Path for raw input images"}
bool    fl_over     {yes,prompt="Subtract overscan level"}
bool    fl_trim     {yes,prompt="Trim overscan section"}
bool    fl_bias     {no,prompt="  Bias-subtract images?"}
string  bias        {"",  prompt="Bias calibration image"}
string  key_biassec {"BIASSEC",prompt="Header keyword for overscan strip image section"}
string  key_datasec {"DATASEC",prompt="Header keyword for data section (excludes the overscan)"}
string  key_ron     {"RDNOISE",prompt="Header keyword for readout noise"}
string  key_gain    {"GAIN",prompt="Header keyword for gain (e-/ADU"}
real    ron         {3.5,prompt="Readout noise value to use if keyword not found"}
real    gain        {2.2,prompt="Gain value to use if keyword not found"}

# Gprepare
char    gaindb      {"default",prompt="Database with gain data"}
string  sci_ext     {"SCI",prompt="Name of science extension"}
string  var_ext     {"VAR",prompt="Name of variance extension"}
string  dq_ext      {"DQ",prompt="Name of data quality extension"}

string  bpm         {"",prompt="Bad Pixel Mask filename"}
string  sat         {"default",prompt="Saturation level in raw images"}
string  nbiascontam {"default", prompt="Number of columns removed from overscan region"}
string  biasrows    {"default", prompt="Rows to use for overscan region"}
# Colbias
bool    fl_inter    {no,prompt="Interactive overscan fitting"}
bool    median      {no,prompt="Use median instead of average in column bias"}
string  function    {"chebyshev",prompt="Overscan fitting function.",enum="legendre|chebyshev|spline1|spline3"}
char    order       {"default", prompt="Order of overscan fitting function."}
real    low_reject  {3.,prompt="Low sigma rejection factor."}
real    high_reject {3.,prompt="High sigma rejection factor."}
int     niterate    {3,prompt="Number of rejection iterations."}

# Gemcombine
string  combine     {"average",prompt="Type of combination operation",enum="average|median"}
string  reject      {"avsigclip",prompt="Type of rejection algorithm",enum="none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip"}
real    lthreshold  {INDEF,prompt="Lower threshold for rejection before scaling"}
real    hthreshold  {INDEF,prompt="Upper threshold for rejection before scaling"}
string  masktype    {"goodvalue", enum="none|goodvalue", prompt="Mask type"}
real    maskvalue   {0., prompt="Mask value"}
string  scale       {"none",prompt="Image scaling",enum="none|mode|median|mean|exposure|@<file>|!<keyword>"}
string  zero        {"none",prompt="Image zero point offset",enum="none|mode|median|mean|@<file>|!<keyword>"}
string  weight      {"none",prompt="Image weights",enum="none|mode|median|mean|exposure|@<file>|!<keyword>"}
string  statsec     {"[*,*]",prompt="Image region for computing statistics"}
string  key_exptime {"EXPTIME",prompt="Header keyword for exposure time"}
int     nlow        {0,min=0,prompt="minmax: Number of low pixels to reject"}
int     nhigh       {1,min=0,prompt="minmax: Number of high pixels to reject"}
int     nkeep       {1,min=0,prompt="Minimum to keep or maximum to reject"}
bool    mclip       {yes,prompt="Use median in sigma clipping algorithms?"}
real    lsigma      {3.,prompt="Lower sigma clipping factor"}
real    hsigma      {3.,prompt="Upper sigma clipping factor"}
string  snoise      {"0.0",prompt="ccdclip: Sensitivity noise (electrons)"}
real    sigscale    {0.1,prompt="Tolerance for sigma clipping scaling correction"}
real    pclip       {-0.5,prompt="pclip: Percentile clipping parameter"}
real    grow        {0.0,prompt="Radius (pixels) for neighbor rejection"}
bool    fl_vardq    {no,prompt="Create variance and data quality frames?"}

bool    verbose     {yes,prompt="Verbose output?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

begin

    string  l_inimages, l_outdark, l_logfile, l_key_biassec, l_key_datasec
    string  l_tsec[6], l_sci_ext, l_bias, s_empty
    string  l_var_ext, l_dq_ext, l_bpm, l_key_ron, l_key_gain, l_function
    string  l_combine, l_masktype
    string  l_reject, l_scale, l_zero, l_weight,l_expname, l_snoise
    string  l_statsec, l_gaindb, l_rawpath, l_biasrows, l_sat
    string  l_nbiascontam, tmp_prefix, l_order

    bool    l_fl_over, l_fl_trim, l_fl_vardq, l_fl_inter, l_median, l_mclip
    bool    l_verbose, l_fl_bias
    bool    l_fl_keepgprep, l_fl_ovstrim

    int     l_niterate, l_nlow, l_nhigh, l_nkeep

    real    l_ron, l_gain, l_low_reject, l_high_reject
    real    l_lthreshold, l_hthreshold, l_maskvalue
    real    l_lsigma, l_hsigma, l_sigscale, l_pclip, l_grow

    struct  l_struct

    # Runtime variables:
    string  tmplist, inlist, suf, img, test
    string  imgroot, giin, giout, combin
    bool    l_fl_giin
    int     ninp, ngin, len, i, j, k, n, nbad

    # Temporary files
    tmplist = mktemp("tmplist")
    giin = mktemp("tmpgiin")
    giout = mktemp("tmpgiout")
    combin = mktemp("tmpcombin")
    tmp_prefix = mktemp("tmp_prefix")

    # Set values:

    l_inimages=inimages
    l_outdark=outdark
    l_bias=bias
    l_logfile=logfile
    l_key_biassec=key_biassec
    l_key_datasec=key_datasec
    l_sci_ext=sci_ext
    l_var_ext=var_ext
    l_dq_ext=dq_ext
    l_bpm=bpm
    l_key_ron=key_ron
    l_key_gain=key_gain
    l_function=function
    l_nbiascontam = nbiascontam
    l_biasrows = biasrows
    l_combine=combine
    l_masktype=masktype
    l_maskvalue=maskvalue
    l_reject=reject
    l_scale=scale
    l_zero=zero 
    l_gaindb=gaindb
    l_rawpath=rawpath
    l_weight=weight
    l_expname=key_exptime
    l_snoise=snoise
    l_statsec=statsec

    l_fl_bias=fl_bias
    l_fl_over=fl_over
    l_fl_trim=fl_trim
    l_fl_inter=fl_inter
    l_median=median
    l_mclip=mclip
    l_verbose=verbose;

    l_sat=sat
    l_order=order
    l_niterate=niterate
    l_nlow=nlow
    l_nhigh=nhigh
    l_nkeep=nkeep

    l_ron=ron
    l_gain=gain
    l_low_reject=low_reject
    l_high_reject=high_reject
    l_lthreshold=lthreshold
    l_hthreshold=hthreshold
    l_lsigma=lsigma
    l_hsigma=hsigma
    l_sigscale=sigscale
    l_pclip=pclip
    l_grow=grow
    l_fl_vardq=fl_vardq

    status=0
    ninp = 1
    nbad=0

    cache ("imgets", "gimverify", "keypar", "fparse", "gemdate")

    # Tests the logfile:
    if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GNSDARK: both gnsdark.logfile and \
                gmos.logfile are empty.", logfile=l_logfile, verbose=l_verbose)
            printlog ("                Using default file gmos.log.",
                logfile=l_logfile, verbose=yes)
        }
    }

    # Start logging
    date | scan(l_struct)
    printlog ("-----------------------------------------------------------\
        -----------------", logfile=l_logfile, verbose=l_verbose)
    printlog ("GNSDARK -- "//l_struct, logfile=l_logfile, verbose=l_verbose)
    printlog ("", logfile=l_logfile, verbose=l_verbose)

    # Tests the SCI extension name:

    if (l_sci_ext == "" || stridx(" ",l_sci_ext)>0) {
        printlog ("ERROR - GNSDARK: Science extension name parameter sci_ext \
            is missing", logfile=l_logfile, verbose=l_verbose)
        goto error
    }

    # Tests the VAR and DQ extension names, if required:

    if (l_fl_vardq) {
        if (l_var_ext == "" || stridx(" ",l_var_ext)>0 || \
            l_dq_ext == "" || stridx(" ",l_dq_ext)>0) {
            
            printlog ("ERROR - GNSDARK: Variance and/or Data Quality \
                extension name",logfile=l_logfile, verbose=l_verbose)
            printlog ("              parameters (var_ext/dq_ext) are missing",
                logfile=l_logfile, verbose=yes)
            goto error
        }
    }

    # Load up input name list

    inlist = l_inimages

    # Empty string is an error

    if (inlist == "" || stridx(" ",inlist)>0) {
        printlog("ERROR - GNSDARK: Input file not specified.",
            logfile=l_logfile, verbose=l_verbose)
        goto error
    }

    # Test for @. Changed the rest as done for gprepare (avoids the use of sed)
    # KL - initialize 'len' before using it!
    len = strlen(inlist)
    if ((substr(inlist,1,1) == "@") \
        && !access(substr(inlist,2,strlen(inlist))) ) {
        
        printlog ("ERROR - GNSDARK: The input list "//substr(inlist,2,len)//\
            " does not exist.", logfile=l_logfile, verbose=l_verbose)
        goto error
    }

    files(inlist,sort-, > tmplist)

    # Logs the relevant parameters:
    if (substr(inlist,1,1) != "@")
        printlog ("Input images     : "//l_inimages,
            logfile=l_logfile, verbose=l_verbose)
    else {
        printlog ("Input images     :", logfile=l_logfile, verbose=l_verbose)
        scanfile = tmplist
        while (fscan(scanfile,img)!=EOF) {
            printlog ("  "//img, logfile=l_logfile, verbose=l_verbose)
        }
        scanfile=""
    }
    printlog ("Output dark image: "//l_outdark,
        logfile=l_logfile, verbose=l_verbose)

    printlog (" ", logfile=l_logfile, verbose=l_verbose)
    printlog (" ", logfile=l_logfile, verbose=l_verbose)

    scanfile = tmplist

    # Tests the output file
    if (l_outdark == "" || stridx(" ",l_outdark)>0) {
        printlog ("ERROR - GNSDARK: The output filename outdark is not \
            defined", logfile=l_logfile,verbose=l_verbose)
        goto error
    }
    # Must be .fits
    len = strlen(l_outdark)
    suf = substr(l_outdark,len-4,len)
    if (suf !=".fits" ) {
        l_outdark = l_outdark//".fits"
    }

    # Test if output already exists:
    if (imaccess (l_outdark)) {
        printlog ("ERROR - GNSDARK: Output file "//l_outdark//" already \
            exists.", logfile=l_logfile, verbose=l_verbose)
        goto error
    }

    # Test the input images:
    ninp = 0

    # variable 'img' is used to check if exists/mef 
    # variable 'imgroot' is used as input name for input lists
    # they are not the same if prev. gprepared images or if rawpath is used.
    while (fscan(scanfile, img) != EOF) {
        # In most cases, we need to call gireduce. The only time l_fl_giin=no 
        # is if using previously gprepared image AND user doesn't want 
        # over/trim
        l_fl_giin = yes

        gimverify ("g"//img)
        if (gimverify.status==0) { #use previously gprepared image
            img = gimverify.outname//".fits"
            imgroot = img
            printlog ("WARNING - GNSDARK: using the previously \
                gprepared "//img, l_logfile, l_verbose) 
            if (!l_fl_over && !l_fl_trim && !l_fl_bias) 
                l_fl_giin = no
        } else {
            if (l_rawpath != "") { # check rawpath image
                gimverify (l_rawpath//img) # not input to gireduce
            } else { # check raw image in current dir
                gimverify(img)
            }
            if(gimverify.status==1) {
                printlog ("ERROR - GNSDARK: File "//img//" not found.", 
                    l_logfile, verbose+)
                nbad+=1
                goto error
            } else if(gimverify.status>1) {
                printlog ("ERROR - GNSDARK: File "//img//" not a MEF FITS \
                    image.", l_logfile, verbose+)
            nbad+=1
            goto error
            }

            # need to use gimverify again to write image
            gimverify(img)
            img = gimverify.outname//".fits" # here is input to gireduce
            imgroot = img
        }

        #Write into the proper list file
        ngin=0
        if (l_fl_giin) {
            print (imgroot, >> giin)
            print (tmp_prefix//imgroot, >> giout)
            if (l_fl_over || l_fl_trim) {
                print (tmp_prefix//imgroot, >> combin)
            } else {
                # gireduce will only call gprepare
                print ("g"//imgroot, >> combin)
            }
           ngin+=1
        } else {
            print(imgroot, >> combin)
        }
    }
    scanfile=""
    print ("")

    # Call gireduce if needed 
    # Don't create poisson/photon variance with gireduce
    if (ngin > 0) {

        # Input bias image -
        s_empty = ""
        print (l_bias) | scan (s_empty)
        l_bias = s_empty
        gimverify (l_bias)
        l_bias = gimverify.outname
        if (gimverify.status != 0 && l_fl_bias) {
            printlog ("ERROR - GNSDARK: Bias image ("//l_bias//") does not \
                exist.", l_logfile, verbose+)
            goto error 
        }

        gireduce ("@"//giin, outimages="@"//giout, sci_ext=l_sci_ext,
            var_ext=l_var_ext, dq_ext=l_dq_ext, fl_vardq=no,
            fl_over=l_fl_over, key_biassec=l_key_biassec, median=l_median,
            fl_inter=l_fl_inter, function=l_function, order=l_order,
            low_reject=l_low_reject, high_reject=l_high_reject,
            niterate=l_niterate, fl_trim=l_fl_trim, key_datasec=l_key_datasec,
            fl_bias=l_fl_bias, bias=l_bias, fl_dark-, fl_flat-, 
            rawpath=l_rawpath, outpref="g", key_ron=l_key_ron, 
            key_gain=l_key_gain, gain=1.0, gaindb=l_gaindb,
            fl_addmdf-, sat=l_sat, logfile=l_logfile, verbose=l_verbose,
            bpm=l_bpm, nbiascontam=l_nbiascontam, fl_mult=no, 
            biasrows=l_biasrows)
        if (gireduce.status != 0) {
            goto error
        }
    }

    # Combine images if more than 1 # only variance is from combining images
    count(combin) | scan(ngin)
    if (ngin>1) {
        gemcombine ("@"//combin, l_outdark, title="Dark", logfile=l_logfile,
            combine=l_combine, offsets="none", reject=l_reject, 
            masktype=l_masktype, maskvalue=l_maskvalue, scale=l_scale,
            zero=l_zero, weight=l_weight, statsec=l_statsec, expname=l_expname,
            lthreshold=l_lthreshold, hthreshold=l_hthreshold, nlow=l_nlow,
            nhigh=l_nhigh, nkeep=l_nkeep, mclip=l_mclip, lsigma=l_lsigma,
            hsigma=l_hsigma, key_ron=l_key_ron, key_gain=l_key_gain, ron=l_ron,
            gain=l_gain, snoise=l_snoise, sigscale=l_sigscale, pclip=l_pclip,
            grow=l_grow, bpmfile=l_bpm, nrejfile="", sci_ext=l_sci_ext,
            var_ext=l_var_ext, dq_ext=l_dq_ext, fl_vardq=l_fl_vardq,
            fl_dqprop=yes, verbose=l_verbose)
        if (gemcombine.status != 0) {
            goto error
        }
    } else {
        printlog ("GNSDARK: Only one input image", l_logfile, l_verbose)
        type (combin) | scan(img)
        copy (img, l_outdark, verbose-)
        if (l_fl_vardq) {
            printlog("GNSDARK: Cannot create variance plane with one \
                dark image", l_logfile, l_verbose)
            printlog("       Dark frame will not have variance plane.", \
                l_logfile, l_verbose)
        }
    }

    # update PHU of output image
    gemdate ()
    gemhedit (l_outdark//"[0]","GNSDARK", gemdate.outdate, 
        "UT Time stamp for GNSDARK", delete-)
    gemhedit (l_outdark//"[0]","GEM-TLM", gemdate.outdate, 
        "UT Last modification with GEMINI", delete-)

    # Clean up:
    goto clean

    # Exit with error subroutine
error:
    status=1
    goto clean

    #wrap up
clean:
    if (access(giout)) {
        imdelete ("@"//giout,verify-, >>& "dev$null")
    }
    delete (tmplist,verify-, >>& "dev$null")
    delete (giin//","//giout//","//combin,verify-, >>& "dev$null")
    #close log file
    if (status == 0) 
        test="no"
    else 
        test="one or more"

    printlog ("", logfile=l_logfile, verbose=l_verbose)
    printlog ("GNSDARK finished with "//test//" errors.",
        logfile=l_logfile, verbose=l_verbose)
    printlog ("------------------------------------------------------------\
        ----------------", logfile=l_logfile, verbose=l_verbose)

end
