# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.

######################################################
procedure gemscombine (inimages,outimage)

######################################################
# Notes:

# Combine all spectra of input images aperture by aperture (Assumes one aperture
# per extension!)

# Example call:
#

# Version Oct 14, 2005 Bryan Miller (BM) - Add handling of BKG extensions
#         Feb 22, 2006 BM - Change to imcombine and handles dispaxis=3
#
#         Jan 31, 2007 Gwen Rudie (GR) - Generalize to 1D, 2D, or 3D
#         Mar 15, 2007 GR - Added project and section parameters
#         May 19, 2007 GR - Generalize section and make more automatic
#         May 22, 2007 GR - Fix deletion of tmpsci*
#         May 05, 2008 GR - Output DQ plane only indicates pixel bad in all
#                           input images
#         Further updates in revision control

# TODO: Allow multiple apetures per extension?

######################################################
# Parameters:

string inimages     {prompt="Input images"}
string outimage     {prompt="Output image"}
string combine      {"average",enum="average|median|sum",prompt="Combine algorithm"}
string reject       {"avsigclip",enum="none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip",prompt="Rejection algorithm"}
string scale        {"none",prompt="Image scaling"}
string zero         {"none",prompt="Image zeropoint offset"}
string weight       {"none",prompt="Image weights"}
string sample       {"",prompt="Wavelength sample region for statistics"}
string section      {"default",prompt="Image section, e.g. [*,*]"}
real   lthreshold   {INDEF,prompt="Lower threshold"}
real   hthreshold   {INDEF,prompt="Upper threshold"}
int    nlow         {0,prompt="minmax: Number of low pixels to reject"}
int    nhigh        {1,prompt="minmax: Number of high pixels to reject"}
real   lsigma       {3.,prompt="Lower sigma clipping factor"}
real   hsigma       {3.,prompt="Upper sigma clipping factor"}
real   radius       {0.,prompt="Growing radius for bad pixels"}
bool   project      {no,prompt="Project highest dimension of input images?"}
bool   fl_vardq     {no,prompt="Propagate VAR/DQ planes?"}
string logfile      {"",prompt="Logfile name"}
bool   verbose      {yes,prompt="Verbose?"}
int    status       {0,prompt="Exit status (0=good)"}
struct *scanfile    {"",prompt="Internal use only"}

######################################################
# Task:

begin

    ####
    # Variable declaration

    # Task messaging variables
    string task_name, err_msg, warn_msg, task_msg, mdelim, msg_types[4], tmpstr
    string pad_types[4], err_pad, warn_pad, task_pad, mpad_delim, log_wrap
    string pkg_name
    int    n_msg_types, lwrap_len
    struct l_date_struct

    # Declare variables not related to input parameters
    string tmpinlist, inlist, infile[500], extlist, dqlist, img, suf, varlist
    string tmpsci, tmpdq, tmpcombdq, interpdq, tmpvar, tmpcombvar
    string tmpbkg, bkglist, tmpcombbkg, tmpcomb
    string l_key_task, l_key_nsciext, l_key_nextend, l_key_gain
    string l_sci_ext, l_var_ext, l_dq_ext
    string tmpextlist, tmpvarlist, tmpdqlist, tmpbkglist, tmptmpsci
    string tmptmpvar, tmptmpdq, tmptmpbkg, tmptmpcomb, tmptmpcombvar
    string tmptmpcombdq, tmptmpcombbkg, tmpinterpdq
    string l_key_exptime, l_key_rdnoise, secstr
    string inphu, outphu, sciextn, currimgextn, pstr
    string test_list, tlist_file, testfile, tmptodel
    string currphu, tmpcursciextn, dqextn, tmpdqextn, varextn, tmpvarextn
    string bkgextn, tmpbkgextn, l_bkg_ext, l_statsec, exit_status[2]
    string tmptmpcursciextn, tmptmpbkgextn, tmptmpvarextn, tmptmpdqextn
    string tmpnrejmask
    real   wavmin, wavmax, wmin, wmax, dwav, crpix, crval1, cd11, crval3, cd33
    real   gain, rdnoise, gaineff, rdneff, rncomb, exptime, mean, ttime, l1, l2
    int    npix, npixy, x1, x2, ndim, junk, curr_nsci, nsci, i, j, nimages
    int    nextend, atpos
    bool   debug, created_output

    # Local variables for input parameters
    string l_inimages
    string l_outimage
    string l_combine
    string l_reject
    string l_scale
    string l_zero
    string l_weight
    string l_sample
    string l_section
    real   l_lthreshold
    real   l_hthreshold
    int    l_nlow
    int    l_nhigh
    real   l_lsigma
    real   l_hsigma
    real   l_radius
    bool   l_project
    bool   l_fl_vardq
    string l_logfile
    bool   l_verbose

    # Initialize local variables
    junk = fscan (inimages, l_inimages)
    junk = fscan (outimage, l_outimage)
    junk = fscan (combine, l_combine)
    junk = fscan (reject, l_reject)
    junk = fscan (scale, l_scale)
    junk = fscan (zero, l_zero)
    junk = fscan (weight, l_weight)
    l_sample = sample
#    junk = fscan (sample, l_sample)
    junk = fscan (section, l_section)
    l_lthreshold = lthreshold
    l_hthreshold = hthreshold
    l_nlow = nlow
    l_nhigh = nhigh
    l_lsigma = lsigma
    l_hsigma = hsigma
    l_radius = radius
    l_project = project
    l_fl_vardq = fl_vardq
    l_logfile = logfile
#    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    ####
    # Script messages

    # Set task and package names
    task_name = strupr("gemscombine")
    pkg_name = "gemtools"

    # Set up message types
    err_msg = "ERROR - "//task_name//": "
    warn_msg = "WARNING - "//task_name//": "
    task_msg = task_name//": "
    mdelim = "____"//strlwr(task_name)//" -- "

    # Store the message types
    n_msg_types = 4
    msg_types[1] = err_msg
    msg_types[2] = warn_msg
    msg_types[3] = task_msg
    msg_types[4] = mdelim

    # Set up the whitespace padding
    for (i = 1; i <= n_msg_types; i += 1) {
        tmpstr = ""
        for (j = 1; j <= strlen(msg_types[i]); j += 1) {
            tmpstr += " "
        }
        pad_types[i] = tmpstr
    }

    # Sore the whitespace padding for a given message
    err_pad = pad_types[1]
    warn_pad = pad_types[2]
    task_pad = pad_types[3]
    mpad_delim = pad_types[4]

    # Add newline to mdelim
    mdelim = "\n"//mdelim

    # Set up the log file wrapper
    # lwrap_len is the number of "-"
    lwrap_len = 60
    log_wrap = ""
    for (i = 1; i <= lwrap_len; i += 1) {
       log_wrap += "-"
    }

    ####
    # Default values

    # Assume task will fail
    status = 1

    created_output = no

    debug = no
    if (debug) {
        l_verbose = yes
    }

    l_key_task = substr(task_name, 1, 8)
    l_key_nsciext = "NSCIEXT"
    l_key_nextend = "NEXTEND"
    l_key_gain = "GAIN"
    l_sci_ext = "SCI"
    l_var_ext = "VAR"
    l_dq_ext = "DQ"
    l_bkg_ext = "BKG"
    l_key_exptime = "EXPTIME"
    l_key_rdnoise = "RDNOISE"

    ####
    # Temporary files

    inlist = mktemp("tmpinlist")
    tmpinlist = mktemp ("tmpinfiles")

    tmptodel = mktemp ("tmptodel")

    tmpextlist = mktemp ("tmpextlist")
    tmpvarlist = mktemp ("tmpvarlist")
    tmpdqlist = mktemp ("tmpdqlist")
    tmpbkglist = mktemp ("tmpbkglist")
    tmptmpsci = mktemp ("tmpsci")
    tmptmpvar = mktemp ("tmpvar")
    tmptmpdq = mktemp ("tmpdq")
    tmptmpbkg = mktemp ("tmpbkg")
    tmptmpcomb = mktemp ("tmpcomb")
    tmptmpcombvar = mktemp ("tmpcombvar")
    tmptmpcombdq = mktemp ("tmpcombdq")
    tmptmpcombbkg = mktemp ("tmpcombkg")
    tmpinterpdq = mktemp ("tmpinterpdq")

    ########

    # Start task: User input and log file checking.

    ####
    # Check log file
    if (l_logfile == "" || stridx(" ",l_logfile) > 0) {
        ## Cannot easily use pkg_name to determine the logfile to use
        l_logfile = gemtools.logfile
        if (l_logfile == "" || stridx(" ",l_logfile) > 0) {
            l_logfile = pkg_name//".log"

            printlog ("\n"//warn_msg//"both "//task_name//".logfile and "//\
                pkg_name//".logfile are empty.\n"//\
                warn_pad//"Using default file "//pkg_name//".log.",
                l_logfile, l_verbose)
        }
    } # End of logfile check

    ####
    # Start log / show time
    date | scan (l_date_struct)
    printlog ("\n"//log_wrap, l_logfile, l_verbose)
    printlog (task_name//" Started -- "//l_date_struct//"\n", \
        l_logfile, l_verbose)

    ####
    # Record all input parameters
    printlog (task_msg//"Input parameters...\n", l_logfile, l_verbose)
    printlog ("    inimages   = "//l_inimages, l_logfile, l_verbose)
    printlog ("    outimage   = "//l_outimage, l_logfile, l_verbose)
    printlog ("    combine    = "//l_combine, l_logfile, l_verbose)
    printlog ("    reject     = "//l_reject, l_logfile, l_verbose)
    printlog ("    scale      = "//l_scale, l_logfile, l_verbose)
    printlog ("    zero       = "//l_zero, l_logfile, l_verbose)
    printlog ("    weight     = "//l_weight, l_logfile, l_verbose)
    printlog ("    sample     = "//l_sample, l_logfile, l_verbose)
    printlog ("    section    = "//l_section, l_logfile, l_verbose)
    if (l_lthreshold == INDEF) {
        pstr = "INDEF"
    } else {
        pstr = l_lthreshold
    }
    printlog ("    lthreshold = "//pstr, l_logfile, l_verbose)
    if (l_lthreshold == INDEF) {
        pstr = "INDEF"
    } else {
        pstr = l_lthreshold
    }
    printlog ("    hthreshold = "//pstr, l_logfile, l_verbose)
    printlog ("    nlow       = "//l_nlow, l_logfile, l_verbose)
    printlog ("    nhigh      = "//l_nhigh, l_logfile, l_verbose)
    printlog ("    lsigma     = "//l_lsigma, l_logfile, l_verbose)
    printlog ("    hsigma     = "//l_hsigma, l_logfile, l_verbose)
    printlog ("    radius     = "//l_radius, l_logfile, l_verbose)
    printlog ("    project    = "//l_project, l_logfile, l_verbose)
    printlog ("    fl_vardq   = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("    logfile    = "//l_logfile, l_logfile, l_verbose)
    printlog ("    verbose    = "//l_verbose, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    ####
    # User input checks
    # Inputs
    test_list = l_inimages
    atpos = stridx("@", test_list)
    if (atpos > 0) {
        tlist_file = substr(test_list, atpos+1, strlen(test_list))
        if (!access(tlist_file)) {
            printlog (err_msg//"Cannot access list file: "//tlist_file, \
                l_logfile, verbose+)
            goto error
        }
    }

    sections (l_inimages, option="fullname", >> tmpinlist)
    scanfile = tmpinlist

    nsci = -1
    i = 0
    while(fscan(scanfile, testfile) != EOF) {
        i += 1
        gimverify (testfile)
        if (gimverify.status != 0) {
            printlog (err_msg//"input image \""//testfile//\
                "\" not a valid MEF or does not exist", l_logfile, verbose+)
            goto error
        }
        testfile = gimverify.outname//".fits"

        inphu = testfile//"[0]"

        # Check task not already ran on this file
        keypar (inphu, l_key_task, silent+)
        if (keypar.found) {
            printlog (err_msg//"image "//testfile//" has already been "//\
                "processed with "//task_name, l_logfile, verbose+)
            goto error
        }

        # Check NSCIEXT exists
        keypar (inphu, l_key_nsciext, silent+)
        if (!keypar.found) {
            printlog (err_msg//"Cannot find "//l_key_nsciext//\
                " keyword in image "//testfile, l_logfile, verbose+)
            goto error
        }

        # Require all images to have the same number of science extensions
        ##TODO possible bug for MOS Data?
        curr_nsci = int(keypar.value)
        if (nsci == -1) {
            nsci = curr_nsci
        } else if (curr_nsci != nsci) {
            printlog (err_msg//"Mismatch of "//l_key_nsciext//\
                " value for image \""//testfile//"\"", l_logfile, verbose+)
            goto error
        }

        ##TODO Test for var and dq inputs if fl_vardq+; reset variable if they
        ##TODO don't exist
        infile[i] = testfile
    } # End of loop checking inputs
    scanfile = ""

    # Check images have actually been supplied and count them for later use
    nimages = i

    if (nimages == 0) {
        printlog (err_msg//"No input images supplied", l_logfile, verbose+)
        goto error
    } else {
        printlog (task_msg//"Number of inputs = "//nimages, \
            l_logfile, l_verbose)
    }

    # Output checks

    # Check output files do not exist
    gimverify (l_outimage)
    if (gimverify.status != 1) {
        printlog (err_msg//"output image \""//l_outimage//\
            "\" already exists", l_logfile, verbose+)
        goto error
    }
    l_outimage = gimverify.outname//".fits"

    # Checks complete

    ####
    # Hard work part

    # Create output image
    # Assumes the header of the first image is correct
    imcopy (infile[1]//"[0]", l_outimage, verbose-)
    created_output = yes
    outphu = l_outimage//"[0]"
    nextend = 0

    # Copy out the MDF
    # This assumes the same MDF for all images
    tinfo (infile[1]//"[MDF]", ttout-)
    if (tinfo.tbltype == "fits") {
        tcopy (infile[1]//"[MDF]", l_outimage, verbose-)
        nextend += 1
    }

    if (l_lthreshold == INDEF) {
        l_lthreshold = -9999.
    }

    # Combine spectra one extension at a time
    for (i = 1; i <= nsci; i += 1) {
        extlist = tmpextlist//"_"//str(i)
        varlist = tmpvarlist//"_"//str(i)
        dqlist = tmpdqlist//"_"//str(i)
        bkglist = tmpbkglist//"_"//str(i)
        tmpsci = tmptmpsci//"_"//str(i)
        tmpvar = tmptmpvar//"_"//str(i)
        tmpdq = tmptmpdq//"_"//str(i)
        tmpbkg = tmptmpbkg//"_"//str(i)
        tmpcomb = tmptmpcomb//"_"//str(i)//".fits"
        tmpcombvar = tmptmpcombvar//"_"//str(i)//".fits"
        tmpcombdq = tmptmpcombdq//"_"//str(i)//".fits"
        tmpcombbkg = tmptmpcombbkg//"_"//str(i)//".fits"
        interpdq = tmpinterpdq//"_"//str(i)//".fits"
        tmpnrejmask = mktemp("tmpnrejmask")//".pl"

        printlog ("", l_logfile, l_verbose)
        printlog (task_msg//"Working on extension # "//i, \
                  l_logfile, l_verbose)

        # Find min/max wavelength and mean dispersion
        wavmin = 999999.
        wavmax = 0.
        dwav = 0.
        npixy = 1
        cd11 = 0.0
        cd33 = 0.0
        crval3 = 0.0

        # Choose a dispersion axis without the header parameter dispaxis or
        # ctype
        sciextn = "["//l_sci_ext//","//i//"]"
        for (j = 1; j <= nimages; j += 1) {
            currimgextn = infile[j]//sciextn

            hselect (currimgextn, \
                "CRPIX1, i_naxis1, CRVAL1, CD1_1, CD3_3, CRVAL3", "yes") | \
                scan (crpix, npix, crval1, cd11, cd33, crval3)

            keypar (currimgextn, "i_naxis2", silent+)
            if (keypar.found) {
               npixy = int (keypar.value)
            }

            if (crval3 > 3000.) {
                wmin = crval3 + (1. - crpix) * cd33
                wmax = crval3 + (real(npix) - crpix) * cd33
                if (wmin < wavmin)
                    wavmin = wmin
                if (wmax > wavmax)
                    wavmax = wmax
                dwav = dwav + cd33
                ndim = 3
            } else {
                wmin = crval1 + (1. - crpix) * cd11
                wmax = crval1 + (real(npix) - crpix) * cd11
                if (wmin < wavmin)
                    wavmin = wmin
                if (wmax > wavmax)
                    wavmax = wmax
                dwav = dwav + cd11
                if (npixy > 1) {
                    ndim = 2
                } else {
                    ndim = 1
                }
            }
        }
        # End of loop over input images for current extension to determine
        # dispersion axis

        if (l_section == "default") {
            # TODO Can this just be an empty string?
            if (ndim == 3) {
                l_section = "[*,*,*]"
            } else if (ndim == 2) {
                l_section = "[*,*]"
            } else {
                l_section = "[*]"
            }
        }

        dwav = dwav / real(nimages)
        ttime = 0.0
        for (j = 1; j <= nimages; j += 1) {
            currphu = infile[j]//"[0]"

            currimgextn = infile[j]//sciextn
            tmpcursciextn = tmpsci//"_"//str(j)//".fits"
            tmptmpcursciextn = "tmp"//tmpsci//"_"//str(j)//".fits"

            dqextn = infile[j]//"["//l_dq_ext//","//i//"]"
            tmpdqextn = tmpdq//"_"//j//".fits"
            tmptmpdqextn = "tmp"//tmpdq//"_"//j//".fits"

            varextn = infile[j]//"["//l_var_ext//","//i//"]"
            tmpvarextn = tmpvar//"_"//j//".fits"
            tmptmpvarextn = "tmp"//tmpvar//"_"//j//".fits"

            bkgextn = infile[j]//"["//l_bkg_ext//","//i//"]"
            tmpbkgextn = tmpbkg//"_"//j//".fits"
            tmptmpbkgextn = "tmp"//tmpbkg//"_"//j//".fits"

            # Interpolate each SCI extension
            scombine (currimgextn, tmpcursciextn, first-, \
                w1 = wavmin, w2 = wavmax, dw = dwav, logfile="", \
                combine="sum", reject="none", blank=(l_lthreshold - 1.))

            # Copy exposure time to extension
            keypar (currphu, l_key_exptime, silent+)
            exptime = real(keypar.value)
            gemhedit (tmpcursciextn, l_key_exptime, exptime, "", delete-)
            ttime += exptime

            # TODO These booleans should be uncoupled
            if (l_fl_vardq && imaccess(dqextn)) {
                # DQ extension
                # TODO expand the DQ extension into the separate BITs then
                # TODO scombine then combine back into one extension
                scombine (dqextn, tmptmpdqextn, first-, w1=wavmin, w2=wavmax, \
                    dw=dwav, logfile="", combine="sum", reject="none")

                # TODO Use all parameters
                # TDOO set a changable threshold
                imexpr ("a > 0.01 ? 1 : a", tmpdqextn, tmptmpdqextn, \
                    outtype="ushort", verbose-)

                if (l_radius != 0.) {
                    imreplace (tmpdqextn, value=1, lower=1, upper=INDEF, \
                        radius=l_radius)
                }

                imdelete (tmptmpdqextn, verify-, >> "dev$null")
                print (tmpdqextn//l_section, >> dqlist)
                print (tmpdqextn, >> tmptodel)

                # TODO Use all parameters
                imexpr ("b == 1 ? c : a", tmptmpcursciextn, tmpcursciextn, \
                    tmpdqextn, (l_lthreshold - 1.), outtype="ref", refim="a", \
                    verbose-)
                imdelete (tmpcursciextn, verify-, >& "dev$null")
                tmpcursciextn = tmptmpcursciextn

                # Variance extension
                if (imaccess (varextn)) {
                    scombine (varextn, tmptmpvarextn, first-, \
                        w1=wavmin, w2=wavmax, dw=dwav, logfile="", \
                        combine="sum", reject="none", blank=(l_lthreshold-1.))

                    # TODO Use all parameters
                    imexpr ("b == 1 ? c : a", tmpvarextn, tmptmpvarextn, \
                        tmpdqextn, (l_lthreshold - 1.), outtype="ref", \
                        refim="a", verbose-)
                    imdelete (tmptmpvarextn, verify-, >& "dev$null")
                    print (tmpvarextn//l_section, >> varlist)
                    print (tmpvarextn, >> tmptodel)
                }
            }
            print (tmpcursciextn//l_section, >> extlist)
            print (tmpcursciextn, >> tmptodel)

            # interpolate each BKG extension
            if (imaccess(bkgextn)) {
                scombine (bkgextn, tmpbkgextn, first-,
                    w1 = wavmin, w2 = wavmax, dw = dwav, logfile="", \
                    combine="sum", reject="none", blank=(l_lthreshold - 1.))
                gemhedit (tmpbkgextn, l_key_exptime, exptime, "", delete-)
                if (imaccess (tmpdqextn)) {
                    imexpr ("b == 1 ? c : a", tmptmpbkgextn, tmpbkgextn, \
                        tmpdqextn, (l_lthreshold - 1.), outtype="ref", \
                        refim="a", verbose-)
                    imdelete (tmpbkgextn, verify-, >& "dev$null")
                    tmpbkgextn = tmptmpbkgextn
                }
                print (tmpbkgextn//l_section, >> bkglist)
                print (tmpbkgextn, >> tmptodel)
            }
        }
        # End of loop over input images for current extension performing
        # interpolation

        # Calculate statsec in pixels for imcombine
        if (l_sample != "") {
            print (l_sample) | fscan ("%f:%f", l1, l2)

            x1 = nint ((l1 - wavmin) / dwav) + 1
            x2 = nint ((l2 - wavmin) / dwav) + 1
            secstr = x1//":"//x2

            if (ndim == 3) {
                secstr = "*,*,"//secstr
            } else if (ndim == 2 && !l_project) {
                secstr = secstr//",*"
            }

        } else {
            # TODO Do we need this? If not move the l_statsec call to the loop
            # TODO above
            if (ndim == 3) {
                secstr = "*,*,*"
            } else if (ndim == 2) {
                secstr = "*,*"
            } else {
                secstr = "*"
            }
        }
        l_statsec = "["//secstr//"]"
        printlog (task_msg//"STATSEC in pixel space: "//l_statsec, \
            l_logfile, l_verbose)

        # Science extensions

        # TODO Could this be GEMCOMBINE call?
        imcombine ("@"//extlist, tmpcomb, combine=l_combine, \
            nrejmasks=tmpnrejmask, \
            project=l_project, reject=l_reject, nlow=l_nlow, nhigh=l_nhigh, \
            lsigma=l_lsigma, hsigma=l_hsigma, logfile=l_logfile, \
            scale=l_scale, zero=l_zero, weight=l_weight, statsec=l_statsec, \
            lthreshold=l_lthreshold, hthreshold=l_hthreshold, \
            blank=(l_lthreshold - 1.), gain=l_key_gain, rdnoise=l_key_rdnoise)

        # Final exposure time
        keypar (tmpcomb, l_key_exptime, silent+)
        exptime = real(keypar.value)

        # Final interpolation of remaining bad pixels
        imexpr ("a == b ? 1 : 0", interpdq, tmpnrejmask, nimages, \
            outtype="ushort", verbose-)
        imdelete (tmpnrejmask, verify-, >>& "dev$null")

        # TODO Including fixpix too?
        fixpix (tmpcomb, interpdq, cinterp=1)
        delete (interpdq, verify-, >>& "dev$null")

        # VAR extension
        if (l_fl_vardq && access (varlist)) {
            imcombine ("@"//varlist, tmpcombvar, combine=l_combine, \
                nrejmasks=tmpnrejmask, \
                project=l_project, reject="none", nlow=l_nlow, nhigh=l_nhigh, \
                lsigma=l_lsigma, hsigma=l_hsigma, logfile=l_logfile, \
                scale="none", zero="none", weight="none", statsec=l_statsec, \
                lthreshold=l_lthreshold, hthreshold=l_hthreshold, \
                blank=(l_lthreshold - 1.), gain=l_key_gain, \
                rdnoise=l_key_rdnoise)

            # final interpolation of remaining bad pixels
            imexpr ("a == b ? 1 : 0", interpdq, tmpnrejmask, nimages, \
                outtype="ushort", verbose-)
            imdelete (tmpnrejmask, verify-, >>& "dev$null")

            fixpix (tmpcombvar, interpdq, cinterp=1)
            delete (interpdq, verify-, >>& "dev$null")
        }

        # DQ extensions
        if (l_fl_vardq && access(dqlist)) {
            imcombine ("@"//dqlist, tmpcombdq, combine="average", \
                reject="none", logfile=l_logfile, project=l_project, \
                outtype="ushort")
        }

        # BKG extensions
        if (access (bkglist)) {
            imcombine ("@"//bkglist, tmpcombbkg, combine=l_combine, \
                nrejmasks=tmpnrejmask, \
                reject=l_reject, nlow=l_nlow, nhigh=l_nhigh, lsigma=l_lsigma, \
                hsigma=l_hsigma, logfile=l_logfile, scale=l_scale, \
                zero=l_zero, weight=l_weight, statsec=l_statsec, \
                lthreshold=l_lthreshold, hthreshold=l_hthreshold, \
                blank=(l_lthreshold - 1.), gain=l_key_gain, \
                rdnoise=l_key_rdnoise)

            imexpr ("a == b ? 1 : 0", interpdq, tmpnrejmask, nimages, \
                outtype="ushort", verbose-)
            imdelete (tmpnrejmask, verify-, >>& "dev$null")

            # final interpolation of remaining bad pixels
            fixpix (tmpcombbkg, interpdq, cinterp=1)
            delete (interpdq, verify-, >>& "dev$null")
        }

        # Add to output file
        imcopy (tmpcomb, \
            l_outimage//"["//l_sci_ext//","//i//",append]", verbose-)
        nextend += 1

        # Exposure time updating
        if (l_scale == "exposure" || l_weight == "exposure") {
            gemhedit (outphu, l_key_exptime, exptime, "", delete-)
        } else {
            gemhedit (outphu, l_key_exptime, (ttime / real(nimages)), "", \
                delete-)
            gemhedit (l_outimage//"["//l_sci_ext//","//i//"]", l_key_exptime, \
                (ttime / real(nimages)), "", delete-)
        }

        if (l_fl_vardq) {
            imcopy (tmpcombvar, \
                    l_outimage//"["//l_var_ext//","//i//",append]", verbose-)
            imcopy (tmpcombdq, \
                    l_outimage//"["//l_dq_ext//","//i//",append]", verbose-)
            nextend += 2
        }

        if (imaccess (tmpcombbkg)) {
            imcopy (tmpcombbkg, \
                l_outimage//"["//l_bkg_ext//","//i//",append]", verbose-)
            nextend += 1
        }

        # Update gain and read noise
        keypar (l_outimage//"["//l_sci_ext//","//i//"]", "NCOMBINE", silent+)
        rncomb = real(keypar.value)

        keypar (l_outimage//"["//l_sci_ext//","//i//"]", l_key_gain, silent+)
        gain = real(keypar.value)

        keypar (l_outimage//"["//l_sci_ext//","//i//"]", l_key_rdnoise, silent+)
        rdnoise = real(keypar.value)

        gaineff = gain
        rdneff = sqrt(rncomb) * rdnoise

        if (l_combine == "average") {
            gaineff = gain * rncomb
        } else if (l_combine == "median") {
            gaineff = 2. * gain * rncomb / 3.
            rdneff = sqrt(2. * rncomb / 3.) * rdnoise
        }

        if (gaineff != gain) {
            gemhedit (l_outimage//"["//l_sci_ext//","//i//"]", l_key_gain, \
                gaineff, "", delete-)
            if (imaccess (tmpcombbkg)) {
                gemhedit (l_outimage//"["//l_bkg_ext//","//i//"]", l_key_gain, \
                    gaineff, "", delete-)
            }

            # Correct the VAR plane for the gain, then S/N is SCI/sqrt (VAR)
            if (l_fl_vardq) {
                imarith (l_outimage//"["//l_var_ext//","//i//"]", "/", \
                    gaineff, \
                    l_outimage//"["//l_var_ext//","//i//", overwrite+]", \
                    verbose-)
            }
        }

        gemhedit (l_outimage//"["//l_sci_ext//","//i//"]", l_key_rdnoise, \
            rdneff, "", delete-)
        if (imaccess (tmpcombbkg)) {
            gemhedit (l_outimage//"["//l_bkg_ext//","//i//"]", l_key_rdnoise, \
                rdneff, "", delete-)
        }

        # clean
        imdelete ("@"//tmptodel, verify-, >& "dev$null")
        delete (tmptodel, verify-, >& "dev$null")
        delete (extlist//","//dqlist//","//bkglist//","//varlist, \
            verify-, >>& "dev$null")
        delete (tmpcombdq//","//tmpcomb//","//interpdq//","//\
            tmpcombvar//","//tmpcombbkg, \
            verify-, >>& "dev$null")
    } # End of loop over sci extensions

    gemhedit (outphu, l_key_nextend, nextend, "", delete-)
    gemhedit (outphu, l_key_nsciext, nsci, "", delete-)

    gemdate()
    gemhedit (outphu, l_key_task, gemdate.outdate, \
        "UT Time stamp for "//task_name, delete-)
    gemhedit (outphu, "GEM-TLM", gemdate.outdate, \
        "Last modification with GEMINI", delete-)

    status = 0
    goto clean

error:

    if (created_output) {
        imdelete (l_outimage, verify-, >& "dev$null")
    }

clean:

    scanfile = ""

    delete (tmpinlist, verify-, >>& "dev$null")

    # Finish time
    date | scan (l_date_struct)
    printlog ("\n"//task_name//" Finished -- "//l_date_struct//"\n", \
        l_logfile, l_verbose)

    exit_status[1] = "GOOD"
    exit_status[2] = "ERROR"

    printlog (task_name//" exit status:  "//exit_status[status + 1]//\
        ".", l_logfile, l_verbose)

    printlog (log_wrap//"\n", l_logfile, l_verbose)

end
