# Copyright(c) 2005-2011 Association of Universities for Research in Astronomy, Inc.

# acooke 29 jan 2005

# Process a flat frame to give a "pinhole" for nssdits


procedure nfflt2pin (inimage)

char    inimage     {prompt = "Flat frame to use as basis for pinhole"}
char    outimage    {"", prompt = "Output file"}
char    outprefix   {"p", prompt = "Prefix to use if outimage not given"}
char    threshold   {"[*,*]", prompt = "Threshold value or section"}
real    scale       {0.05, prompt = "Fraction of threshold level to use"}
int     dispaxis    {INDEF, prompt="Default dispersion axis (if not in header)"}
char    logfile     {"", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose output?"}
int     status      {0, prompt = "O: Exit status (0 = good)"}

begin
    char    l_inimage = ""
    char    l_outimage = ""
    char    l_outprefix = ""
    char    l_threshold = ""
    real    l_scale
    int     l_dispaxis
    char    l_logfile = ""
    bool    l_verbose

    char    l_sci_ext = ""
    char    l_key_dispaxis = ""

    int     junk, nsig
    bool    debug, havesecn
    char    badhdr, cjunk, sci
    char    tmprot, tmppeak
    real    thresh, midpt, sigma
    struct  sdate, line

    junk = fscan (inimage, l_inimage)
    junk = fscan (outimage, l_outimage)
    junk = fscan (outprefix, l_outprefix)
    l_threshold = threshold
    l_scale = scale
    l_dispaxis = dispaxis
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    status = 1
    debug = no
    nsig = 5 # number of sigma to reject in stats

    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"

    tmprot = mktemp ("tmprot")
    tmppeak = mktemp ("tmppeak")

    cache ("gemextn", "gemdate")


    # Start logging

    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile) 
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NFFLT2PIN: Both nfflt2pin.logfile and \
                gnirs.logfile are empty.", l_logfile, verbose+) 
            printlog ("                   Using default file " \
                // l_logfile // ".", l_logfile, verbose+) 
        }
    }
    date | scan (sdate) 
    printlog ("-----------------------------------------------------------\
        -------------------", l_logfile, verbose = l_verbose) 
    printlog ("NFFLT2PIN -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 


    # Check input data

    if (debug) print ("checking input")
    gemextn (l_inimage, proc="none", check="exists,mef", index="", \
        extname="", extver="", ikparams="", omit="path", replace="", \
        outfile="STDOUT", logfile="", glogpars="", verbose-) \
        | scan (l_inimage)
    if (0 < gemextn.fail_count || no == (1 == gemextn.count)) {
        printlog ("ERROR - NFFLT2PIN: Missing / duplicate input image.", \
            l_logfile, verbose+) 
        goto clean
    }
    gemextn (l_inimage, proc="expand", check="exists", index="", \
        extname=l_sci_ext, extver="-", ikparams="", omit="path", \
        replace="", outfile="STDOUT", logfile="", glogpars="", \
        verbose-) | scan (sci)
    if (0 < gemextn.fail_count || no == (1 == gemextn.count)) {
        printlog ("ERROR - NFFLT2PIN: Missing / duplicate input data.", \
            l_logfile, verbose+) 
        goto clean
    }


    # Check/generate output names

    if (debug) print ("checking output")
    gemextn (l_outimage, proc="none", check="absent", index="", \
        extname="", extver="", ikparams="", omit="", replace="", \
        outfile="STDOUT", logfile="", glogpars="", verbose-) \
        | scan (line)
    if (0 < gemextn.fail_count) {
        printlog ("ERROR - NFFLT2PIN: Output data already exist.", \
            l_logfile, verbose+) 
        goto clean
    }
    if (0 != gemextn.count) {
        if (no == (1 == gemextn.count)) {
            printlog ("ERROR - NFFLT2PIN: Incorrect number of output \
                files.", l_logfile, verbose+) 
            goto clean
        } else {
            l_outimage = line
        }
    } else {
        gemextn (l_outprefix // l_inimage,
            proc="none", check="absent", index="", extname="", \
            extver="", ikparams="", omit="", replace="", \
            outfile="STDOUT", logfile="", glogpars="", verbose-) \
            | scan (l_outimage)
        if (0 < gemextn.fail_count || no == (1 == gemextn.count)) {
            printlog ("ERROR - NFFLT2PIN: Bad output file.", \
                l_logfile, verbose+) 
            goto clean
        }
    }

    if (debug) print ("inimage:  " // l_inimage)
    if (debug) print ("sci:      " // sci)
    if (debug) print ("outimage: " // l_outimage)


    # Parse threshold

    thresh = INDEF
    print (l_threshold) | scan (thresh)
    havesecn = isindef (thresh)

    if (havesecn) {
        if (fscan (l_threshold, cjunk) == 0) {
            printlog ("ERROR - NFFLT2PIN: Missing threshold.", \
                l_logfile, verbose+)
            goto clean
        }
        printlog ("NFFLT2PIN: Measuring threshold from " \
            // sci // l_threshold, l_logfile, l_verbose)
        midpt = INDEF
        sigma = INDEF
        imstat (sci // l_threshold, fields="midpt,stddev", \
            lower=INDEF, upper=INDEF, nclip=0, lsigma=3, usigma=3, \
            binwidth=0.1, format-, cache-) | scan (midpt, sigma)
        if (isindef (midpt) || isindef (sigma)) {
            printlog ("ERROR - NFFLT2PIN: Bad statistics from " \
                // sci // l_threshold, l_logfile, verbose+)
            goto clean
        }
        if (debug) {
            print ("midpt: " // midpt)
            print ("sigma: " // sigma)
        }
        imstat (sci // l_threshold, fields="midpt,stddev", \
            lower=(midpt-nsig*sigma), upper=(midpt+nsig*sigma), nclip=0, \
            lsigma=3, usigma=3, binwidth=0.1, format-, cache-) \
            | scan (midpt, sigma)
        if (isindef (midpt) || isindef (sigma)) {
            printlog ("ERROR - NFFLT2PIN: Bad statistics from " \
                // sci // l_threshold, l_logfile, verbose+)
            goto clean
        }
        printf ("NFFLT2PIN: Level in threshold section is %5.2f \
            +/- %5.2f\n", midpt, sigma) | scan (line)
        printlog (line, l_logfile, l_verbose)

        if (l_scale > 0) {
            thresh = midpt * l_scale
        } else {
            thresh = midpt + sigma * abs (l_scale)
        }
    }

    printf ("NFFLT2PIN: Threshold %5.2f\n", thresh) | scan (line)
    printlog (line, l_logfile, l_verbose)


    # Rotate data if necessary

    hselect (l_inimage // "[0]", l_key_dispaxis, yes) | scan (dispaxis)
    if (isindef (dispaxis)) {
        printlog ("ERROR - NFFLT2PIN: No dispersion axis.", \
            l_logfile, verbose+)
        goto clean
    }
    if (1 == dispaxis) {
        if (debug) print ("rotating counter-clock")
        imtranspose (sci // "[*,-*]", tmprot, len_blk=5112)
        sci = tmprot
    }


    # Generate the new data

    peakhelper (sci, tmppeak, threshold=thresh)

    if (1 == dispaxis) {
        if (debug) print ("rotating clock")
        imdelete (tmprot, verify-, >& "dev$null")
        tmprot = mktemp ("tmprot")
        imtranspose (tmppeak // "[-*,*]", tmprot, len_blk=5112)
        imdelete (tmppeak, verify-, >& "dev$null")
        tmppeak = tmprot
    }

    copy (l_inimage // ".fits", l_outimage // ".fits", verbose-)

    gemextn (l_outimage, proc="expand", check="exists", index="", \
        extname=l_sci_ext, extver="-", ikparams="overwrite", omit="path", \
        replace="", outfile="STDOUT", logfile="", glogpars="", \
        verbose-) | scan (sci)
    imcopy (tmppeak, sci, verbose-)
    
    gemdate ()
    gemhedit (l_outimage//"[0]", "NFFLT2PIN", gemdate.outdate,
        "UT Time stamp for NFFLT2PIN", delete-)
    gemhedit (l_outimage//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI IRAF", delete-)


    status = 0 # Success

clean:

    imdelete (tmprot, verify-, >& "dev$null")
    imdelete (tmppeak, verify-, >& "dev$null")

    if (0 == status) {
        printlog (" ", l_logfile, l_verbose) 
        printlog ("NFFLT2PIN exit status: good.", l_logfile, l_verbose) 
    }

end
