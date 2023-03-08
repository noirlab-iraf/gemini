# Copyright(c) 2001-2011 Association of Universities for Research in Astronomy, Inc.

procedure nsressky (inimages) 

# script to remove residual sky from GNIRS/NIRI longslit spectra
# 
# Version Sept 20, 2002 MT,JJ v1.4 release
#         Aug 19, 2003  KL  IRAF2.12 - new parameter, addonly, in hedit
#         Oct 29, 2003  KL  moved from niri to gnirs package

char    inimages    {prompt = "Input spectra (2D)"}
char    outspectra  {"", prompt = "Output spectra"}
char    outprefix   {"s", prompt = "Prefix for output images"}
char    sample      {"*", prompt = "Sample to fit"}
char    function    {"chebyshev", prompt = "Function to fit"}
int     order       {2, prompt = "Order for fit"}
int     dispaxis    {1, min=1, max=2, prompt="Dispersion axis if not defined in the header"}
bool    fl_inter    {no, prompt = "Fit interactively?"}
char    logfile     {"", prompt = "Logfile name"}
bool    verbose     {yes, prompt = "Verbose?"}
int     status      {0, prompt = "Exit status (0=good)"}
struct  *scanfile1  {"", prompt = "For internal use only"}
struct  *scanfile2  {"", prompt = "For internal use only"}

begin

    char    l_inimages = ""
    char    l_outspectra = ""
    char    l_outprefix = ""
    char    l_sample = ""
    char    l_function = ""
    int     l_order
    int     l_dispaxis
    bool    l_fl_inter
    char    l_logfile = ""
    bool    l_verbose

    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_dispaxis = ""

    char    filelist, inlist, outlist 
    char    inimg, outimg, tmpfit, sciinlist, scioutlist, insci, outsci
    char    img, suf, badhdr, scilist, keyfound
    int     n, nin, nout, nbad, axis
    bool    bad, mef, debug
    struct  sdate
    int     junk

    debug = no
    status = 1

    junk = fscan (  inimages, l_inimages)
    junk = fscan (  outspectra, l_outspectra)
    junk = fscan (  outprefix, l_outprefix)
    junk = fscan (  sample, l_sample)
    junk = fscan (  function, l_function)
    l_order =       order
    l_dispaxis =    dispaxis
    l_fl_inter =    fl_inter
    junk = fscan (  logfile, l_logfile)
    l_verbose =     verbose

    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"

    cache ("keypar", "gemextn", "gemdate") 

    filelist = mktemp ("tmpfilelist") 
    inlist = mktemp ("tmpinlist") 
    outlist = mktemp ("tmpoutlist") 


    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile)
        if ("" == l_logfile) {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSRESSKY: Both nsressky.logfile and \
                gnirs.logfile are", l_logfile, verbose+) 
            printlog ("                     undefined.  Using " \
                // l_logfile, l_logfile, verbose+)
        }
    }

    date | scan (sdate) 
    printlog ("----------------------------------------------------------\
        ---------------------", l_logfile, l_verbose) 
    printlog ("NSRESSKY -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 
    printlog ("input = " // l_inimages, l_logfile, l_verbose) 
    printlog ("output = " // l_outspectra, l_logfile, l_verbose) 
    printlog ("sample = " // l_sample, l_logfile, l_verbose) 
    printlog ("function = " // l_function, l_logfile, verbose = l_verbose) 
    printlog ("order = " // l_order, l_logfile, verbose = l_verbose) 
    printlog (" ", l_logfile, l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NSRESSKY: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }


    gemextn (l_inimages, check="exists,mef", process="none", index="", \
        extname="", extversion="", ikparams="", omit="kernel,exten", \
        replace="", outfile=inlist, logfile="", glogpars="",
        verbose=l_verbose)
    if (0 < gemextn.fail_count) {
        printlog ("ERROR - NSRESSKY: Bad syntax in inimages.", \
            l_logfile, verbose+) 
        goto clean
    }
    if (0 == gemextn.count) {
        printlog ("ERROR - NSRESSKY: No input images.", \
            l_logfile, verbose+) 
        goto clean
    }
    nin = gemextn.count

    nbad = 0
    scanfile1 = inlist
    while (fscan (scanfile1, inimg) != EOF) {
        keyfound = ""
        hselect (inimg // "[0]", "*PREPAR*", yes) | scan (keyfound)
        if (keyfound == "") {
            printlog ("ERROR - NSRESSKY: Image " // inimg \
                // " not PREPAREd.", l_logfile, verbose+) 
            nbad += 1
        }
    }
    if (nbad > 0) goto clean


    gemextn (l_outspectra, check="absent", process="none", index="", \
        extname="", extversion="", ikparams="", omit="kernel,exten", \
        replace="", outfile=outlist, logfile="", glogpars="", 
        verbose=l_verbose)
    if (gemextn.fail_count != 0) {
        printlog ("ERROR - NSRESSKY: Existing or incorrectly formatted \
            output files", l_logfile, verbose+) 
        goto clean
    }
    if (gemextn.count == 0) {
        gemextn ("%^%" // l_outprefix // "%" // "@" // inlist, \
            check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", omit="kernel,exten", \
            replace="", outfile=outlist, logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - NSRESSKY: No or incorrectly formatted \
                output files", l_logfile, verbose+) 
            goto clean
        }
    }
    nout = gemextn.count

    if (nin != nout) {
        printlog ("ERROR - NSSRESSKY: different numbers of input and \
            output images", l_logfile, verbose+) 
        goto clean
    }


    joinlines (inlist, outlist, output=filelist, delim=" ", verbose-) 
    scanfile1 = filelist
    while (fscan (scanfile1, inimg, outimg) != EOF) {

        sciinlist = mktemp ("tmpsciin")
        scioutlist = mktemp ("tmpsciout")
        scilist = mktemp ("tmpsci")

        gemextn (inimg, check="exists,mef", process="expand", index="", \
            extname=l_sci_ext, extversion="1-", ikparams="", omit="", \
            replace="", outfile=sciinlist, logfile="", glogpars="",
            verbose=l_verbose)

        printlog (inimg // " -> " // outimg, l_logfile, l_verbose) 

        gemextn (inimg, check="", process="expand", index="1-", \
            extname="", extversion="", ikparams="", omit="", \
            replace="", outfile="dev$null", logfile="dev$null", glogpars="",
            verbose-)
        fxcopy (inimg, outimg // ".fits", groups=("0-"//gemextn.count),
            new_file+, verbose-)

        if (0 == gemextn.count) {
            printlog ("WARNING - NSRESSKY: No science data in " \
                // inimg, l_logfile, verbose+) 
        } else {
            gemextn (outimg, check="", process="append", index="", \
                extname=l_sci_ext, extversion="1-" // gemextn.count, \
                ikparams="overwrite", omit="", replace="", \
                outfile=scioutlist, logfile="", glogpars="", 
                verbose=l_verbose)

            joinlines (sciinlist, scioutlist, output=scilist, delim=" ", \
                verbose-) 
            scanfile2 = scilist
            while (fscan (scanfile2, insci, outsci) != EOF) {

                if (debug) print (insci // " => " // outsci)
                tmpfit = mktemp ("tmpfit") 

                keypar (insci, l_key_dispaxis)
                if (keypar.found) axis = int (keypar.value)
                else              axis = l_dispaxis
                axis = 3 - axis # fit along spatial direction

                if (debug) print ("fit1d dispaxis=" // axis)
                fit1d (insci, tmpfit, "fit", axis=axis, \
                    function=l_function, \
                    order=l_order, low_reject=3., high_reject=3., \
                    niter=3, interactive=l_fl_inter, sample=l_sample)

                if (debug) print ("imarith")
                imarith (insci, "-", tmpfit, outsci, verbose-) 

                # Here the fit^2 was added to the variance, but I don't
                # see why, so have dropped it.

                imdelete (tmpfit, verify-, >& "dev$null") 
            }

            gemdate ()
            gemhedit (outimg // "[0]", "NSRESSKY", gemdate.outdate, \
                "UT Time stamp for NSRESSKY") 
            gemhedit (outimg // "[0]", "GEM-TLM", gemdate.outdate, 
                "UT Last modification with GEMINI", delete-) 

        }

        delete (sciinlist, verify-, >& "dev$null")
        delete (scioutlist, verify-, >& "dev$null")
        delete (scilist, verify-, >& "dev$null")

    }

    status = 0

clean:
    scanfile1 = ""
    scanfile2 = ""
    delete (filelist // "," // inlist // "," // outlist, verify-, 
        >& "dev$null") 

    printlog (" ", l_logfile, l_verbose) 
    printlog ("NSRESSKY done", l_logfile, l_verbose) 
    printlog ("--------------------------------------------------------------\
        ------------------", l_logfile, l_verbose) 
        
end
