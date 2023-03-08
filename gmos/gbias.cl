# Copyright(c) 2002-2015 Association of Universities for Research in Astronomy, Inc.

procedure gbias(inimages,outbias)

# Derives average GMOS bias images from raw or gprepared images
# Updates GAIN and RDNOISE in the final (combined) image (SCI exts)
# Assumes that BIAS frames do not contain an MDF
#
# Version  Feb 28, 2002  CW,BM,IJ  v1.3 release
#          Aug 10, 2002  IJ   offsets in gemcombine
#          Aug 19, 2002  IJ   make it work on one input image
#          Sep 20, 2002  IJ   v1.4 release
#          Oct 17, 2003  KL   now outputs correct error when input list file
#                             does not exist.
#          Apr 16, 2008 JH Fixed file naming logic. Now works for no
#                       overscan subtraction/trim.

string  inimages    {prompt="Input GMOS bias images or list"}            # OLDP-1-input-primary-combine-suffix=_bias
string  outbias     {prompt="Output bias (zero level) image"}            # OLDP-1-output
string  logfile     {"",prompt="Logfile"}                                   # OLDP-1
string  rawpath     {"",prompt="GPREPARE: Path for raw input images"}       # OLDP-4
bool    fl_over     {yes,prompt="Subtract overscan level?"}                  # OLDP-2
bool    fl_trim     {yes,prompt="Trim overscan section?"}                    # OLDP-2
string  key_biassec {"BIASSEC",prompt="Header keyword for overscan strip image section"}  # OLDP-3
string  key_datasec {"DATASEC",prompt="Header keyword for data section (excludes the overscan)"}  # OLDP-3
string  key_ron     {"RDNOISE",prompt="Header keyword for readout noise"}   # OLDP-3
string  key_gain    {"GAIN",prompt="Header keyword for gain (e-/ADU)"}       # OLDP-3
real    ron         {3.5,prompt="Readout noise value to use if keyword not found"}  # OLDP-3
real    gain        {2.2,prompt="Gain value to use if keyword not found"} # OLDP-3

# Gprepare
char    gaindb      {"default",prompt="Database containing gain data"} # OLDP-2
string  sci_ext     {"SCI",prompt="Name of science extension"}              # OLDP-3
string  var_ext     {"VAR",prompt="Name of variance extension"}             # OLDP-3
string  dq_ext      {"DQ",prompt="Name of data quality extension"}          # OLDP-3
string  bpm         {"",prompt="Bad Pixel Mask filename"}                   # OLDP-2-input
string  sat         {"default",prompt="Saturation level in raw images"}   # OLDP-3
string  nbiascontam {"default", prompt="Number of columns removed from overscan region"}
string  biasrows    {"default", prompt="Rows to use for overscan region"}
# Colbias
bool    fl_inter    {no,prompt="Interactive overscan fitting"}              # OLDP-3
bool    median      {no,prompt="Use median instead of average in column bias"}  # OLDP-3
string  function    {"chebyshev",prompt="Overscan fitting function.",enum="legendre|chebyshev|spline1|spline3"}  # OLDP-3
char    order       {"default",prompt="Order of overscan fitting function."}        # OLDP-3
real    low_reject  {3.,prompt="Low sigma rejection factor"}               # OLDP-3
real    high_reject {3.,prompt="High sigma rejection factor"}              # OLDP-3
int     niterate    {3,prompt="Number of rejection iterations"}          # OLDP-3

# Gemcombine
string  combine     {"average",prompt="Type of combination operation",enum="average|median"}  # OLDP-3
string  reject      {"avsigclip",prompt="Type of rejection algorithm",enum="none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip"}  # OLDP-3
real    lthreshold  {INDEF,prompt="Lower threshold for rejection before scaling"}  # OLDP-3
real    hthreshold  {INDEF,prompt="Upper threshold for rejection before scaling"}  # OLDP-3
string  masktype    {"goodvalue", enum="none|goodvalue", prompt="Mask type"}   # OLDP-3
real    maskvalue   {0., prompt="Mask value"}           # OLDP-3
string  scale       {"none",prompt="Image scaling",enum="none|mode|median|mean|exposure|@<file>|!<keyword>"}  # OLDP-3
string  zero        {"none",prompt="Image zero point offset",enum="none|mode|median|mean|@<file>|!<keyword>"}  # OLDP-3
string  weight      {"none",prompt="Image weights",enum="none|mode|median|mean|exposure|@<file>|!<keyword>"}  # OLDP-3
string  statsec     {"[*,*]",prompt="Image region for computing statistics"}    # OLDP-3

string  key_exptime {"EXPTIME",prompt="Header keyword for exposure time"}   # OLDP-3
int     nlow        {0,min=0,prompt="minmax: Number of low pixels to reject"}   # OLDP-3
int     nhigh       {1,min=0,prompt="minmax: Number of high pixels to reject"}  # OLDP-3
int     nkeep       {1,min=0,prompt="Minimum to keep or maximum to reject"} # OLDP-3

bool    mclip       {yes,prompt="Use median in sigma clipping algorithms?"} # OLDP-3
real    lsigma      {3.,prompt="Lower sigma clipping factor"}               # OLDP-3
real    hsigma      {3.,prompt="Upper sigma clipping factor"}               # OLDP-3
string  snoise      {"0.0",prompt="ccdclip: Sensitivity noise (electrons)"} # OLDP-3
real    sigscale    {0.1,prompt="Tolerance for sigma clipping scaling correction"}  # OLDP-3
real    pclip       {-0.5,prompt="pclip: Percentile clipping parameter"}    # OLDP-3
real    grow        {0.0,prompt="Radius (pixels) for neighbor rejection"} # OLDP-3
bool    fl_vardq    {no,prompt="Create variance and data quality frames?"}  # OLDP-2

bool    verbose     {yes,prompt="Verbose output?"}                          # OLDP-4
int     status      {0,prompt="Exit status (0=good)"}                       # OLDP-4
struct  *scanfile   {"",prompt="Internal use only"}                         # OLDP-4

begin

    string  l_inimages, l_outbias, l_logfile, l_key_biassec, l_key_datasec
    string  l_tsec[12], l_sci_ext
    string  l_var_ext, l_dq_ext, l_bpm, l_key_ron, l_key_gain, l_function
    string  l_combine, l_masktype, l_biasrows
    string  l_reject, l_scale, l_zero, l_weight,l_expname, l_snoise
    string  l_statsec, l_gaindb, l_rawpath, l_sat, l_nbiascontam, l_order

    bool    l_fl_over, l_fl_trim, l_fl_vardq, l_fl_inter, l_median, l_mclip
    bool    l_verbose
    bool    l_fl_keepgprep, l_fl_ovstrim

    int     l_niterate, l_nlow, l_nhigh, l_nkeep, n_files

    real    l_ron, l_gain, l_low_reject, l_high_reject
    real    l_lthreshold, l_hthreshold, l_maskvalue
    real    l_lsigma, l_hsigma, l_sigscale, l_pclip, l_grow

    struct  l_struct

    # Runtime variables:
    string  tmplist, inlist, suf, img, tmpfile, img_root, tmpoutgiredlist
    string  giin,giout,combin, pathtest, rphend, testfile, tmpgiredlist
    bool    l_fl_giin
    int     ninp, ngin, len, i, j, k, n, nbad, atposition, n_lines, l_maxfiles

    # Temporary files
    tmplist = mktemp("tmplist")
    giin = mktemp("tmpgiin")
    giout = mktemp("tmpgiout")
    combin = mktemp("tmpcombin")
    tmpgiredlist = mktemp("tmpgiredlist")
    tmpoutgiredlist = mktemp("tmpoutgiredlist")

    # Set values:

    l_inimages=inimages
    l_outbias=outbias
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

    l_fl_over=fl_over
    l_fl_trim=fl_trim
    l_fl_inter=fl_inter
    l_median=median
    l_mclip=mclip
    l_verbose=verbose

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

    tmpfile = mktemp ("tmpfile")

    # The maximum number of files gireduce can handle
    l_maxfiles = gireduce.maxfiles

    status=0
    ninp = 1
    nbad=0
    ngin = 0

    cache ("imgets", "gimverify", "keypar", "fparse", "gemdate")

    # Tests the logfile:
    if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GBIAS: both gbias.logfile and gmos.logfile \
                are empty.", logfile=l_logfile, verbose+)
            printlog ("                Using default file gmos.log.",
                logfile=l_logfile, verbose+)
        }
    }

    # Start logging
    date | scan(l_struct)
    printlog ("-----------------------------------------------------------\
        -----------------",logfile=l_logfile, verbose=l_verbose)
    printlog ("GBIAS -- "//l_struct, logfile=l_logfile, verbose=l_verbose)
    printlog ("", logfile=l_logfile, verbose=l_verbose)

    # Tests the SCI extension name:
    if (l_sci_ext == "" || stridx(" ",l_sci_ext)>0) {
        printlog ("ERROR - GBIAS: Science extension name parameter sci_ext \
            is missing", logfile=l_logfile, verbose+)
        goto error
    }

    # Tests the VAR and DQ extension names, if required:

    if (l_fl_vardq) {
        if (l_var_ext == "" || stridx(" ",l_var_ext)>0 || \
            l_dq_ext == "" || stridx(" ",l_dq_ext)>0) {

            printlog ("ERROR - GBIAS: Variance and/or Data Quality extension \
                name",logfile=l_logfile, verbose+)
            printlog ("              parameters (var_ext/dq_ext) are missing",
                logfile=l_logfile, verbose+)
            goto error
        }
    }

    # Load up input name list

    inlist = l_inimages

    # Empty string is an error
    if (inlist == "" || stridx(" ",inlist)>0) {
        printlog("ERROR - GBIAS: Input file not specified.",
            logfile=l_logfile, verbose+)
        goto error
    }

    # Test for @. Changed the rest as done for gprepare (avoids the use of sed)
    # KL - initialize 'len' before using it!
    len = strlen(inlist)
    atposition = stridx("@",inlist)
    if (atposition > 0) {
        testfile = substr(inlist,atposition+1,len)
        if (!access(testfile)) {
            printlog ("ERROR - GBIAS: The input list "//testfile//\
                " does not exist.", logfile=l_logfile, verbose+)
            goto error
        }
    }

    files(inlist,sort-, > tmplist)
    count(tmplist) | scan(n_files)

    if (n_files == 0) {
        printlog ("ERROR - GBIAS: No input files supplied", \
            logfile=l_logfile, verbose=yes)
        goto error
    }

    # Logs the relevant parameters:
    printlog ("Input images:", logfile=l_logfile, verbose=l_verbose)
    scanfile = tmplist
    while (fscan(scanfile,img)!=EOF) {
        printlog ("  "//img, logfile=l_logfile, verbose=l_verbose)
    }
    scanfile=""

    printlog ("\nOutput bias image: "//l_outbias,
        logfile=l_logfile, verbose=l_verbose)

    printlog (" ",logfile=l_logfile, verbose=l_verbose)

    scanfile = tmplist

    # Tests the output file
    if (l_outbias == "" || stridx(" ",l_outbias)>0) {
        printlog ("ERROR - GBIAS: The output filename outbias is not defined",
            logfile=l_logfile,verbose+)
        goto error
    }
    # Must be .fits
    len = strlen(l_outbias)
    suf = substr(l_outbias,len-4,len)
    if (suf !=".fits" ) {
        l_outbias = l_outbias//".fits"
    }

    # Test if output already exists:
    if (imaccess (l_outbias)) {
        printlog ("ERROR - GBIAS: Output file "//l_outbias//" already exists.",
            logfile=l_logfile, verbose+)
        goto error
    }

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
            printlog ("ERROR - GBIAS: Cannot access rawpath: "//l_rawpath, \
                l_logfile, verbose+)
            goto error
        }
    }

    # Test the input images:
    ninp = 0

    # variable 'img' is used to check if exists/mef
    # variable 'img_root' is used as input name for input lists
    # they are not the same if prev. gprepared images or if rawpath is used.
    while (fscan(scanfile, img) != EOF) {
        # In most cases, we need to call gireduce. The only time l_fl_giin=no
        # is if using previously gprepared image AND user doesn't want
        # over/trim
        l_fl_giin = yes

        fparse (img)
        img_root = fparse.root//".fits"

        gimverify ("g"//img_root)
        if (gimverify.status==0) { #use previously gprepared image
            printlog ("WARNING - GBIAS: using the previously \
                gprepared "//img_root, l_logfile, l_verbose)
            if (!l_fl_over && !l_fl_trim) {
                l_fl_giin = no
            }

        } else {
            if (l_rawpath != "") { # check rawpath image
                img = l_rawpath//img # not input to gireduce
            }

            gimverify(img)

            if(gimverify.status==1) {
                printlog ("ERROR - GBIAS: File "//img//" not found.",
                    l_logfile, verbose+)
                nbad+=1
                goto error
            } else if(gimverify.status>1) {
                printlog ("ERROR - GBIAS: File "//img//" not a MEF FITS \
                    image.", l_logfile, verbose+)
                nbad+=1
                goto error
            }
        }

        img = gimverify.outname//".fits"

        #Write into the proper list file
        ngin=0
        if (l_fl_giin) {
            print (img, >> giin)
            print (tmpfile//img_root, >> giout)
            if (l_fl_over || l_fl_trim) {
                print (tmpfile//img_root, >> combin)
            } else {
                # gireduce will only call gprepare
                print ("g"//img_root, >> combin)
            }
            ngin+=1
        } else {
            print(img, >> combin)
        }
    }
    scanfile=""
    print ("")

    # Call gireduce if needed
    # Don't create poisson/photon variance with gireduce
    if (ngin > 0) {
        count (giin) | scan (n_lines)
        while (n_lines > 0) {
            tail (giin, nlines=n_lines) |
                head ("STDIN", nlines=l_maxfiles, > tmpgiredlist)
            tail (giout, nlines=n_lines) |
                head ("STDIN", nlines=l_maxfiles, > tmpoutgiredlist)

            gireduce ("@"//tmpgiredlist, outimages="@"//tmpoutgiredlist, sci_ext=l_sci_ext,
                var_ext=l_var_ext, dq_ext=l_dq_ext, fl_vardq=no,
                fl_over=l_fl_over, key_biassec=l_key_biassec, median=l_median,
                fl_inter=l_fl_inter, function=l_function, order=l_order,
                low_reject=l_low_reject, high_reject=l_high_reject,
                niterate=l_niterate, fl_trim=l_fl_trim,
                key_datasec=l_key_datasec,
                fl_bias-, fl_dark-, fl_flat-, rawpath="", outpref="g",
                key_ron=l_key_ron, key_gain=l_key_gain, gain=1.0,
                gaindb=l_gaindb,
                fl_addmdf-, sat=l_sat, logfile=l_logfile, verbose=l_verbose,
                bpm=l_bpm, nbiascontam=l_nbiascontam, fl_mult=no,
                biasrows=l_biasrows)

            delete (tmpgiredlist//","//tmpoutgiredlist, verify-, >& "dev$null")
            n_lines -= l_maxfiles

            if (gireduce.status != 0) {
                goto error
            }
        }
    }
    # Combine images if more than 1 # only variance is from combining images
    count(combin) | scan(ngin)
    if (ngin>1) {
        gemcombine ("@"//combin, l_outbias, title="Bias", logfile=l_logfile,
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
        printlog ("GBIAS: Only one input image", l_logfile, l_verbose)
        type (combin) | scan(img)
        copy (img, l_outbias, verbose-)
        if (l_fl_vardq) {
            printlog("GBIAS: Cannot create variance plane with one \
                bias image", l_logfile, l_verbose)
            printlog("       Bias frame will not have variance plane.", \
                l_logfile, l_verbose)
        }
    }

    # update PHU of output image
    gemdate ()
    gemhedit (l_outbias//"[0]","GBIAS", gemdate.outdate,
        "UT Time stamp for GBIAS", delete-)
    gemhedit (l_outbias//"[0]","GEM-TLM", gemdate.outdate,
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

    printlog ("", logfile=l_logfile, verbose=l_verbose)
    if (status == 0) {
        printlog ("GBIAS: Exit status: GOOD.", \
            logfile=l_logfile, verbose=l_verbose)
    } else  {
        printlog ("GBIAS: Exit status: ERROR.", \
            logfile=l_logfile, verbose=l_verbose)
    }

    printlog ("------------------------------------------------------------\
        ----------------", logfile=l_logfile, verbose=l_verbose)

end
