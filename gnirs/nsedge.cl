# Copyright(c) 2005-2011 Association of Universities for Research in Astronomy, Inc.

procedure nsedge (inimages, bpmfile)

char    inimages        {prompt = "Images to use for detection"}
char    bpmfile         {prompt = "Output bad pixel mask"}
char    flatbpm         {"", prompt = "Input bad pixel mask (from NSFLAT)"}
real    threshold       {0.5, prompt = "Fraction of maximum at which to clip"}
int     grow            {1, prompt = "Extension radius"}
char    logfile         {"", prompt = "Logfile"}
bool    verbose         {yes, prompt = "Verbose output?"}

int     status          {0, prompt = "Exit status (0 = good)"}
struct* scanin1         {"", prompt = "Internal use only"}

begin
        char    l_inimages = ""
        char    l_bpmfile = ""
        char    l_flatbpm = ""
        real    l_threshold
        int     l_grow
        char    l_logfile = ""
        bool    l_verbose

        char    l_sci_ext = ""
        char    l_dq_ext = ""
        char    l_var_ext = ""
        char    l_key_dispaxis = ""
        char    l_key_cut_section = ""

        int     junk, dispaxis, nin, version, nversion, nx, ny
        bool    debug, havebpm
        struct  sline
        char    badhdr, img, angle, reffile, ref, section
        char    tmpin, tmpraw, tmpsci, tmpgrad, tmpmask1, tmpmask2, tmpcomb
        char    tmpmask3

        junk = fscan (  inimages, l_inimages)
        junk = fscan (  bpmfile, l_bpmfile)
        junk = fscan (  flatbpm, l_flatbpm)
        l_threshold =   threshold
        l_grow =        grow
        junk = fscan (  logfile, l_logfile)
        l_verbose =     verbose

        status = 1
        debug = no
        badhdr = ""
        nin = 0

        tmpin = mktemp ("tmpin")
        tmpraw = mktemp ("tmpraw")
        tmpsci = mktemp ("tmpsci")
        tmpgrad = mktemp ("tmpgrad")
        tmpmask1 = mktemp ("tmpmask1")
        tmpmask2 = mktemp ("tmpmask2")
        tmpmask3 = mktemp ("tmpmask3")
        tmpcomb = mktemp ("tmpcomb")

        junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
        if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
        junk = fscan (nsheaders.key_cut_section, l_key_cut_section)
        if ("" == l_key_cut_section) badhdr = badhdr + " key_cut_section"
        junk = fscan (nsheaders.sci_ext, l_sci_ext)
        if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
        junk = fscan (nsheaders.dq_ext, l_dq_ext)
        if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
        junk = fscan (nsheaders.var_ext, l_var_ext)
        if ("" == l_var_ext) badhdr = badhdr + " var_ext"

        cache ("gemextn", "minmax", "gemcombine")

        if ("" == l_logfile) {
            junk = fscan (gnirs.logfile, l_logfile) 
            if (l_logfile == "") {
                l_logfile = "gnirs.log"
                printlog ("WARNING - NSEDGE: Both nsedge.logfile and \
                    gnirs.logfile are empty.", l_logfile, verbose+) 
                printlog ("                  Using default file " \
                    // l_logfile // ".", l_logfile, verbose+) 
            }
        }

        date | scan (sline) 
        printlog ("---------------------------------------------------------\
            ---------------------", l_logfile, verbose = l_verbose) 
        printlog ("NSEDGE -- " // sline, l_logfile, l_verbose) 
        printlog (" ", l_logfile, l_verbose) 


        if ("" != badhdr) {
            printlog ("ERROR - NSEDGE: Parameter(s) missing from \
                nsheaders: " // badhdr, l_logfile, verbose+) 
            goto clean
        }


        # check output available
        if ("" == l_bpmfile) {
            printlog ("ERROR - NSEDGE: No bpmfile specified.", \
                l_logfile, verbose+)
            goto clean
        } else {
            gemextn (l_bpmfile, check="absent", process="none", index="", \
                extname="", extversion="", ikparams="", omit="", \
                replace="", outfile="STDOUT", logfile=l_logfile, glogpars="", \
                verbose=l_verbose) | scan (sline)
            if (gemextn.fail_count > 0 || gemextn.count != 1) {
                printlog ("ERROR - NSEDGE: Output file already exists.", \
                    l_logfile, verbose+)
                goto clean
            } else {
                l_bpmfile = sline
            }
        }

        # check/expand input
        gemextn (l_inimages, check="exists", process="none", index="", \
            extname="", extversion="", ikparams="", omit="", \
            replace="", outfile=tmpin, logfile=l_logfile, glogpars="", \
            verbose=l_verbose)
        if (gemextn.fail_count > 0 || gemextn.count == 0) {
            printlog ("ERROR - NSEDGE: Problems with input images.", \
                l_logfile, verbose+)
            goto clean
        } else {
             nin = gemextn.count
        }


        # take first input as a reference for constructing final bpm
        scanin1 = tmpin
        junk = fscan (scanin1, sline)
        reffile = sline


        # check bpm available
        if (no == ("" == l_flatbpm)) {
            gemextn (l_flatbpm, check="exists,image", process="none", \
                index="", extname="", extversion="", ikparams="", omit="", \
                replace="", outfile="STDOUT", logfile=l_logfile, glogpars="", \
                verbose=l_verbose) | scan (sline)
             if (gemextn.fail_count > 0 || gemextn.count != 1) {
                printlog ("ERROR - NSEDGE: Flat BPM file missing.", \
                    l_logfile, verbose+)
                goto clean
            } else {
                l_flatbpm = sline
                havebpm = yes
            }
        } else {
            havebpm = no
        }


        # combine input
        if (nin > 1) {
            gemcombine ("@" // tmpin, output=tmpraw, title="", \
                combine="median", reject="none", offsets="none", \
                masktype="goodvalue", maskvalue=0., \
                scale="none", zero="none", weight="none", statsec="[*,*]", \
                expname="EXPTIME", lthreshold=INDEF, hthreshold=INDEF, \
                nlow=1, nhigh=1, nkeep=1, mclip+, lsigma=3, hsigma=3, \
                key_ron="RDNOISE", key_gain="GAIN", ron=0, gain=1, \
                snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0, bpmfile="", \
                nrejfile="", sci_ext=l_sci_ext, var_ext=l_var_ext, \
                dq_ext=l_dq_ext, \
                fl_vardq-, logfile=l_logfile, fl_dqprop-, verbose-)
            if (no == (0 == gemextn.status)) goto clean
        } else {
            scanin1 = tmpin
            junk = fscan (scanin1, sline)
            tmpraw = sline
        }


        # need dispersion axis
        keypar (tmpraw // "[0]", l_key_dispaxis, silent+)
        if (no == keypar.found) {
            printlog ("ERROR - NSEDGE: No " // l_key_dispaxis \
                // " for dispersion axis.", l_logfile, verbose+)
            goto clean
        } else {
            dispaxis = keypar.value
        }
        if (1 == dispaxis) angle = "180"
        else               angle = "0"


        # process each science extension in the combined image
        gemextn (tmpraw, check="exists,image", process="expand", \
            index="", extname=l_sci_ext, extversion="1-", ikparams="", \
            omit="", replace="", outfile=tmpsci, logfile=l_logfile, \
            glogpars="", verbose=l_verbose)
        if (gemextn.fail_count > 0 || gemextn.count == 0) {
            printlog ("ERROR - NSEDGE: Problems with combined image.", \
                l_logfile, verbose+)
            goto clean
        } else {
            nversion = gemextn.count
        }


        # construct an output file
        if (1 == nversion) {
            ref = reffile // "[" // l_sci_ext // ",inherit]"
        } else {
            ref = reffile // "[" // l_sci_ext // ",1,inherit]"
        }
        hselect (ref, "naxis1", yes) | scan (nx)
        hselect (ref, "naxis2", yes) | scan (ny)
        hselect (ref, "ORIGXSIZ", yes) | scan (nx)
        hselect (ref, "ORIGYSIZ", yes) | scan (ny)
        mkimage (l_bpmfile, "make", 0, ndim = 2, dims = nx // " " // ny, \
            pixtype = "short") 


        # add the old bpm file
        if (havebpm) {
            addmask (l_bpmfile // "," // l_flatbpm, tmpcomb, "im1 || im2", \
                flags="")
            imdelete (l_bpmfile, verify-, >& "dev$null")
            imrename (tmpcomb, l_bpmfile, >& "dev$null")
            tmpcomb = mktemp ("tmpcomb")
        }


        scanin1 = tmpsci
        version = 0
        while (no == (EOF == fscan (scanin1, img))) {

            version = version + 1
            if (debug) print (version // " " // img)

            # calculate gradient
            imdelete (tmpgrad, verify-, >& "dev$null") 
            tmpgrad = mktemp ("tmpgrad")
            if (debug) display (img, 1)
            gradient (img, tmpgrad, angle, boundary="nearest", constant=0)
            if (debug) display (tmpgrad, 2)

            # find thresholds
            minmax (tmpgrad, force+, update-, verbose-)

            # generate edge mask
            delete (tmpmask1, verify-, >& "dev$null") 
            tmpmask1 = mktemp ("tmpmask1")
            imexpr ("a > " // l_threshold * minmax.maxval \
                // " || a < " // l_threshold * minmax.minval \
                // " ? 1 : 0", tmpmask1, tmpgrad, dims="auto", intype="auto", \
                outtype="int", refim="auto", bwidth=0, btype="nearest", \
                bpixval=0, rangecheck+, verbose-, exprdb="none")
            if (debug) display (tmpmask1, 3)

            # grow by radius
            delete (tmpmask2, verify-, >& "dev$null") 
            tmpmask2 = mktemp ("tmpmask2")
            if (l_grow > 0) {
                crgrow (tmpmask1, tmpmask2, radius=l_grow, inval=INDEF, \
                    outval=INDEF, >& "dev$null")  # prints warning to screen
            } else {
                tmpmask2 = tmpmask1
            }
            if (debug) display (tmpmask2, 4)

            # copy to destination
            if (1 == nversion) {
                ref = reffile // "[" // l_sci_ext // ",inherit]"
            } else {
                ref = reffile // "[" // l_sci_ext // "," // version \
                    // ",inherit]"
            }
            if (debug) print (ref)
            keypar (ref, l_key_cut_section, silent+)
            if (keypar.found) {
                section = keypar.value
                if (debug) print (section)
            } else {
                keypar (ref, "NSCUT", silent+)
                if (nversion > 1 || keypar.found) {
                    printlog ("ERROR - NSEDGE: Cut data, but no cut section " \
                        // l_key_cut_section // ".", l_logfile, verbose+)
                    goto clean
                }
                section = "[*,*]"
            }
            delete (tmpmask3, verify-, >& "dev$null") 
            tmpmask3 = mktemp ("tmpmask3")
            mkimage (tmpmask3, "make", 0, ndim = 2, dims = nx // " " // ny, \
                pixtype = "short") 
            imcopy (tmpmask2, tmpmask3 // section, verbose-)

            addmask (l_bpmfile // "," // tmpmask3, tmpcomb, \
                "im1 || im2", flags="")
            imdelete (l_bpmfile, verify-, >& "dev$null")
            imrename (tmpcomb, l_bpmfile, >& "dev$null")
            tmpcomb = mktemp ("tmpcomb")

        }

        status = 0

clean:

        scanin1 = ""

        delete (tmpin, verify-, >& "dev$null")
        if (nin > 1) imdelete (tmpraw, verify-, >& "dev$null")
        delete (tmpsci, verify-, >& "dev$null")
        imdelete (tmpmask1, verify-, >& "dev$null")
        imdelete (tmpmask2, verify-, >& "dev$null")
        imdelete (tmpmask3, verify-, >& "dev$null")
        imdelete (tmpgrad, verify-, >& "dev$null")

        printlog (" ", l_logfile, l_verbose) 
        if (0 == status) {
            printlog ("NSEDGE exit status:  good.", l_logfile, l_verbose) 
        } else {
            printlog ("NSEDGE exit status:  error.", l_logfile, l_verbose) 
        }
        printlog ("---------------------------------------------------------\
            ----------------------", l_logfile, l_verbose) 

end
