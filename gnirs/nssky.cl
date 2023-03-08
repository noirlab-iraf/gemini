# Copyright(c) 2004-2013 Association of Universities for Research in Astronomy, Inc.

# This code brought to you by acooke, 2 Dec 2003.

# Generate sky frames for each input file, avoiding duplication where
# possible.

# Previously this was part of nsreduce

# Returns a file that has 5 columns:
# - a file from inimages
# - the corresponding sky file
# - the corresponding bpm file (may not exist if not requested)
# - the name of the file containing a list of all files used to make the sky
# - a flag that is true if this sky is repeated from earlier in the list

# Normally the calling routine will need to delete:
# - all the sky files
# - all the bpm files
# - all the list files
# - the returned file
# (the last column is useful in helping avoid duplicate deletions)


procedure nssky (inimages)

char    inimages    {prompt = "Input image(s)"}
char    skyimages   {"", prompt = "Sky image(s) from other nod positions"}

real    distance    {3., prompt = "Radius (arcsec) from ref for target list"}
real    age         {INDEF, prompt = "Maximum time difference (s) (both lists)"}

char    combtype    {"median", enum = "average|median", prompt = "Type of combine operation for sky"}
char    rejtype     {"avsigclip", enum = "none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip", prompt = "Type of rejection for combining sky"}
char    masktype    {"goodvalue", enum = "none|goodvalue", prompt="Mask type"}
real    maskvalue   {0., prompt="Mask value"}
char    scale       {"none", prompt = "Image scaling for combining sky (see imcombine.scale)"}
char    zero        {"median", prompt = "Image zero-point offset for combining sky (see imcombine.zero)"}
char    weight      {"none", prompt = "Image weights for combining sky (see imcombine.weight)"}
char    statsec     {"[*,*]", prompt = "Statistics section"}
real    lthreshold  {INDEF, prompt = "Lower threshold"}
real    hthreshold  {INDEF, prompt = "Upper threshold"}
int     nlow        {1, min = 0, prompt = "minmax: Number of low pixels to reject"}
int     nhigh       {1, min = 0, prompt = "minmax: Number of high pixels to reject"}
int     nkeep       {0, prompt = "Minimum to keep or maximum to reject"}
bool    mclip       {yes, prompt = "Use median in sigma clipping algorithms?"}
real    lsigma      {3., min = 0., prompt = "Lower sigma clipping factor"}
real    hsigma      {3., min = 0., prompt = "Upper sigma clipping factor"}
char    snoise      {"0.0", prompt = "ccdclip: Sensitivity noise (electrons)"}
real    sigscale    {0.1, min = 0., prompt = "Tolerance for sigma clipping scaling correction"}
real    pclip       {-0.5, prompt = "Percentile clipping parameter"}
real    grow        {0.0, min = 0., prompt = "Radius (pixels) for neighbor rejection"}

bool    fl_vardq    {yes, prompt = "Create variance and data quality frames?"}

char    index       {"", prompt = "Output file (if empty, will be generated)"}

char    logfile     {"", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose output?"}
bool    debug       {no, prompt = "Very verbose output?"}
bool    force       {no, prompt = "Force use with earlier IRAF versions?"}

int     status      {0, prompt = "Exit status (0=good)"}

struct  *scanin1    {prompt="Internal use"}
struct  *scanin2    {prompt="Internal use"}

begin
    char    l_inimages = ""
    char    l_skyimages = ""
    real    l_distance, l_age
    char    l_combtype = ""
    char    l_rejtype = ""
    char    l_masktype = ""
    real    l_maskvalue
    char    l_scale = ""
    char    l_zero = ""
    char    l_weight = ""
    char    l_statsec = ""
    real    l_lthreshold, l_hthreshold
    int     l_nlow, l_nhigh, l_nkeep
    bool    l_mclip
    real    l_lsigma, l_hsigma
    char    l_snoise = ""
    real    l_sigscale, l_pclip, l_grow
    bool    l_fl_vardq
    char    l_index = ""
    char    l_logfile = ""
    bool    l_verbose
    bool    l_debug
    bool    l_force

    char    l_key_exptime = ""
    char    l_key_xoff = ""
    char    l_key_yoff = ""
    char    l_key_date = ""
    char    l_key_time = ""
    char    l_key_ron = ""
    char    l_key_gain = ""
    char    l_key_dispaxis = ""
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""

    int     junk, nsky, nbad, nobj, dispaxis
    bool    intdbg, done, already, first
    char    tmpin, tmpskyin, tmpsel, tmpsky, tmpbpm, tmptime, tmpjoin, tmplist
    char    ref, sky, bpm, sel, ref2, img
    char    differ, line, badhdr
    real    t, deltat, agemin, agemax, prevt
    struct  sdate, sline


    junk = fscan (  inimages, l_inimages)
    junk = fscan (  skyimages, l_skyimages)
    l_distance =    distance
    l_age =         age
    junk = fscan (  combtype, l_combtype)
    junk = fscan (  rejtype, l_rejtype)
    junk = fscan (  masktype, l_masktype)
    l_maskvalue =   maskvalue
    junk = fscan (  scale, l_scale)
    junk = fscan (  zero, l_zero)
    junk = fscan (  weight, l_weight)
    junk = fscan (  statsec, l_statsec)
    l_lthreshold =  lthreshold
    l_hthreshold =  hthreshold
    l_nlow =        nlow
    l_nhigh =       nhigh
    l_nkeep =       nkeep
    l_mclip =       mclip
    l_lsigma =      lsigma
    l_hsigma =      hsigma
    junk = fscan (  snoise, l_snoise)
    l_sigscale =    sigscale
    l_pclip =       pclip
    l_grow =        grow
    l_fl_vardq =    fl_vardq
    junk = fscan (  index, l_index)
    junk = fscan (  logfile, l_logfile)
    l_debug =       debug
    l_verbose =     verbose
    l_force =       force

    badhdr = ""
    junk = fscan (  nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (  nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (  nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (  nsheaders.key_ron, l_key_ron)
    if ("" == l_key_ron) badhdr = badhdr + " key_ron"
    junk = fscan (  nsheaders.key_gain, l_key_gain)
    if ("" == l_key_gain) badhdr = badhdr + " key_gain"
    junk = fscan (  nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (  nsheaders.key_xoff, l_key_xoff)
    if ("" == l_key_xoff) badhdr = badhdr + " key_xoff"
    junk = fscan (  nsheaders.key_yoff, l_key_yoff)
    if ("" == l_key_yoff) badhdr = badhdr + " key_yoff"
    junk = fscan (  nsheaders.key_date, l_key_date)
    if ("" == l_key_date) badhdr = badhdr + " key_date"
    junk = fscan (  nsheaders.key_time, l_key_time)
    if ("" == l_key_time) badhdr = badhdr + " key_time"
    junk = fscan (  nsheaders.key_exptime, l_key_exptime)
    if ("" == l_key_exptime) badhdr = badhdr + " key_exptime"

    status = 1
    intdbg = no

    tmpin = mktemp ("tmpin")
    tmpskyin = mktemp ("tmpskyin")
    tmptime = mktemp ("tmptime")
    tmplist = mktemp ("tmplist")

    cache ("gemcombine", "gemextn")


    # Start logging

    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile) 
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSSKY: Both nssky.logfile and \
                gnirs.logfile are empty.", l_logfile, verbose+) 
            printlog ("                   Using default file " \
                // l_logfile // ".", l_logfile, verbose+) 
        }
    }
    date | scan (sdate) 
    printlog ("---------------------------------------------------------\
        ---------------------", l_logfile, verbose = l_verbose) 
    printlog ("NSSKY -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NSSKY: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }


    # Validate input

    if (intdbg) print ("checking input")
    gemextn (l_inimages, check = "exists,mef", process = "none", \
        index = "", extname = "", extversion = "", ikparams = "", \
        omit = "", replace = "", outfile = tmpin, logfile="", glogpars="",
        verbose=l_verbose)
    if (gemextn.fail_count != 0 || gemextn.count == 0) {
        printlog ("ERROR - NSSKY: Bad input images.", l_logfile, 
            verbose+) 
        goto clean
    }
    nobj = gemextn.count

    if (intdbg) print ("checking sky images")
    gemextn (l_skyimages, check = "exists,mef", process = "none", \
        index = "", extname = "", extversion = "", ikparams = "", \
        omit = "", replace = "", outfile = "dev$null", \
        logfile="dev$null", glogpars="", verbose-)
    if (gemextn.fail_count == 0 && gemextn.count == 0) {
        printlog ("WARNING - NSSKY: Will take sky from input images", \
            l_logfile, verbose+)
        tmpskyin = tmpin
        nsky = 0
    } else {
        gemextn (l_skyimages, check = "exists,mef", process = "none", \
            index = "", extname = "", extversion = "", ikparams = "", \
            omit = "", replace = "", outfile = tmpskyin, \
            logfile="", glogpars="", verbose=l_verbose)
        if (gemextn.fail_count != 0) {
            printlog ("ERROR - NSSKY: Bad sky data.", l_logfile, verbose+)
            goto clean
        }
        nsky = gemextn.count
    }

    # Need an output file

    if (intdbg) print ("checking output")
    if ("" == l_index) {
        l_index = mktemp ("index")
        index = l_index
    }
    if (access (l_index)) {
        printlog ("WARNING - NSSKY: Appending to an existing index \
            file.", l_logfile, verbose+) 
    } else {
        touch (l_index)
    }


    # Deal with the "fancy" processing if l_age is INDEF

    if (isindef (l_age)) {

        if (intdbg) print ("fancy processing")

        if (0 == nsky) {

            # Try generating a suitable value for age, if required
            sections ("@"//tmpin//"//.fits\[0]", option="fullname",
                > tmplist)
            hselect ("@"//tmplist, "$I,UT", yes) | sort (> tmptime)
            delete (tmplist, verify-, >& "dev$null")

            scanin1 = tmptime
            first = yes
            while (fscan (scanin1, img, t) != EOF) {
                if (first) {
                    agemin = 24.0
                    agemax = 0.0
                    first = no
                    printf ("NSSKY: Obs time for %-25s: %h\n", img, t) | \
                        scan (sline)
                } else {
                    deltat = t - prevt
                    if (deltat < 0) deltat = 24. + deltat # midnight crossing
                    printf ("NSSKY: Obs time for %-25s: %h %h\n", \
                        img, t, deltat) | scan (sline)
                    # ignore looping through midnight
                    if (deltat < 12) {
                        agemin = min (agemin, deltat)
                        agemax = max (agemax, deltat)
                    }
                }
                printlog (sline, l_logfile, l_verbose)
                prevt = t
            }
            scanin1 = ""
            delete (tmptime, verify-, >& "dev$null")

            printf ("NSSKY: Min/max: %h/%h\n", agemin, agemax) \
                | scan (sline)
            printlog (sline, l_logfile, l_verbose)

            # just slightly bigger than the maximum separation, to include
            # everything (epsilon both as factor and increment to avoid
            # numerical rounding at large and small scales, i think)
            l_age = (1.0 + epsilon) * agemax * 3600 + epsilon

            # we run the risk of including more than one neighbour
            # "from one side" if we can fit two minimum jumps inside
            # the interval. Round that down to be safe
            if ((1.0 - epsilon) * agemin * 2 * 3600 - epsilon < l_age) {
                printlog ("ERROR - NSSKY: Sorry, but the observing times \
                    are too variably spaced", l_logfile, verbose+) 
                printlog ("               for a single age parameter \
                    value to be inferred.", l_logfile, verbose+) 
                goto clean
            }

            printf ("NSSKY: Using observations within %5.1fs as sky.\n", \
                l_age) | scan (sline)
            printlog (sline, l_logfile, verbose+)

        } else {

            if (nsky == 1 || nsky == nobj) {

                head (tmpin) | scan (sline)

                if (1 == nsky) {
                    printlog ("NSSKY: Using single sky for all: " \
                        // sline, l_logfile, l_verbose)
                } else {
                    printlog ("NSSKY: Using explicit sky list.", \
                        l_logfile, l_verbose)
                }

                tmpjoin = mktemp ("tmpjoin")
                joinlines (tmpin // "," // tmpskyin, "", output=tmpjoin, \
                    delim=" ", missing=sline, maxchars=1000, \
                    shortest-, verbose-)
                delete (tmpin, verify-, >& "dev$null")
                tmpin = tmpjoin
                if (intdbg) print ("1-1 sky:")
                if (intdbg) type (tmpjoin)

            } else {

                printlog ("ERROR - NSSKY: Number of sky frames does \
                    not match the ", l_logfile, verbose+)
                printlog ("               number of object frames.", \
                    l_logfile, verbose+)
                goto clean

            }
        }
    }


    # Generate the sky for each input file in turn

    nbad = 0
    tmpbpm = "none"
    printlog ("NSSKY: Grouping images (please wait).", l_logfile, verbose+)
    scanin1 = tmpin
    while (fscan (scanin1, ref, tmpsky) != EOF) {
        tmpsel = mktemp ("tmpsel")
        done = no
        if (intdbg) print ("processing: " // ref)

        if (isindef (l_age)) {

            print (tmpsky, >> tmpsel)

        } else {

            # Infer spatial direction
            # Currently not used / under discussion - do we want distance
            # along the slit, or absolute separation?

            keypar (ref // "[0]", l_key_dispaxis, silent+)
            if (keypar.found) {
                dispaxis = int (keypar.value)
            } else {
                printlog ("WARNING - NSSKY: No " // l_key_dispaxis \
                    // " in " // ref // ".", l_logfile, verbose+)
                printlog ("                 Assuming 2", \
                    l_logfile, verbose+)
                dispaxis = 2
            }

            # List candidates
            if (intdbg) print ("grouping")
            gemoffsetlist (infiles = "@" // tmpskyin, reffile = ref, \
                distance = l_distance, age = l_age, fl_younger = yes, \
                #direction = (3-dispaxis), \
                direction = 3, \
                fl_nearer = no, fl_noref = yes, wcs_source = "phu", \
                targetlist = tmpsel, \
                offsetlist = "dev$null", key_xoff = l_key_xoff, \
                key_yoff = l_key_yoff, key_date=l_key_date, 
                key_time = l_key_time, \
                logfile = l_logfile, verbose = l_debug, force = l_force)

            if (0 != gemoffsetlist.status) {
                printlog ("ERROR - NSSKY: Failure grouping objects.", \
                    l_logfile, verbose+)
                goto clean
            }

            # May have no candidates
            if (no == access (tmpsel)) touch (tmpsel)
            gemextn ("@" // tmpsel, proc="none", check="", index="", \
                extname="", extver="", ikparam="", replace="", omit="", \
                outfile="dev$null", logfile="dev$null", glogpars="", \
                verbose-)
            nsky = gemextn.count
            if (intdbg) type (tmpsel)
            if (intdbg) print ("nsky: " // nsky)

            if (0 == nsky) {

                printlog ("WARNING - NSSKY: No sky found for " // ref,
                    l_logfile, verbose+)
                delete (tmpsel, verify-, >& "dev$null")
                nbad = nbad + 1
                next

            } else if (1 == nsky) {

                if (intdbg) print ("single sky image")
                type (tmpsel) | scan (tmpsky)
                tmpbpm = "none"

            } else {

                # Check whether list already exists
                if (intdbg)
                print ("checking for duplicates in index " // l_index)
                if (intdbg) type (index)
                scanin2 = l_index
                while (EOF != fscan (scanin2, ref2, sky, bpm, sel, already)) {
                    if (no == already) { # Avoid checking duplicates
                        differ = ""
                        diff (sel, tmpsel) | scan (differ)
                        if ("" == differ) {
                            if (intdbg) print (ref // " matches " // ref2)
                            delete (tmpsel, verify-, >& "dev$null")
                            tmpsel = sel
                            tmpsky = sky
                            tmpbpm = bpm
                            done = yes
                            break
                        }
                    }
                }

                # If not, combine
                if (no == done) {
                    if (intdbg) print ("combining")
                    tmpsky = mktemp ("tmpsky")
                    tmpbpm = mktemp ("tmpbpm")
                    gemcombine("@" // tmpsel, tmpsky, logfile = logfile, \
                        title="sky image produced by gemini.gnirs.nssky", \
                        combine=l_combtype, reject=l_rejtype, \
                        masktype=l_masktype, maskvalue=l_maskvalue, \
                        scale=l_scale, zero=l_zero, weight=l_weight, \
                        statsec=l_statsec, expname=l_key_exptime, \
                        lthreshold=l_lthreshold, hthreshold=l_hthreshold, \
                        nlow=l_nlow, nhigh=l_nhigh, nkeep=l_nkeep, \
                        mclip=l_mclip, lsigma=l_lsigma, \
                        hsigma=l_hsigma, key_ron=l_key_ron, \
                        key_gain=l_key_gain, snoise=l_snoise, \
                        sigscale=l_sigscale, pclip=l_pclip, \
                        grow=l_grow, bpmfile=tmpbpm, nrejfile="", \
                        sci_ext=l_sci_ext, var_ext=l_var_ext, \
                        dq_ext=l_dq_ext, fl_vardq=l_fl_vardq, fl_dqprop-, \
                        verbose=l_debug)
                    # don't hide possible error messages - worried about old
                    # iraf versions etc 
                    #verbose = l_verbose, >& "dev$null")
                    if (0 != gemcombine.status) {
                        printlog("ERROR - NSSKY: Error in GEMCOMBINE.", \
                            l_logfile, verbose+)
                        goto clean
                    }
                    if (intdbg) 
                        print ("access " // tmpsky // ": " // access (tmpsky))
                }
            }
        }

        # Record info to avoid repeating this processing
        line = ref // " " // tmpsky // " " // tmpbpm // " " // tmpsel \
            // " " // done
        print (line, >> l_index)
        if (intdbg) print line
    }

    if (nbad > 0) {
        printlog("ERROR - NSSKY - One or more observations has no sky.", \
            l_logfile, verbose+)
        goto clean
    }

    status = 0 # Success

    if (intdbg) print ("final index")
    if (intdbg) type (l_index)

clean:

    delete (files = tmpin, verify-, go_ahead+, >& "dev$null")
    delete (files = tmpskyin, verify-, go_ahead+, >& "dev$null")
    delete (files = tmptime, verify-, go_ahead+, >& "dev$null")

    scanin1 = ""
    scanin2 = ""

    if (0 == status) {
        printlog (" ", l_logfile, l_verbose) 
        printlog ("NSSKY exit status: good.", l_logfile, l_verbose) 
    }

end
