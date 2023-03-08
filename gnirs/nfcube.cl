# Copyright(c) 2004-2011 Association of Universities for Research in Astronomy, Inc.

# Initial version: acooke 13 Nov 2004

# Rearrange IFU data into a cube


procedure nfcube (inimages)

char    inimages    {prompt = "Data to convert to 3D cubes"}
char    outimages   {"", prompt = "Output files"}
char    outprefix   {"c", prompt = "Prefix to use if outimages not given"}
char    process     {"none", enum = "none|rotate|resample", prompt = "Processing for cube"}
int     dispaxis    {INDEF, prompt="Default dispersion axis (if not in header)"}

char    logfile     {"", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose output?"}

int     status      {0, prompt = "O: Exit status (0 = good)"}
struct  *scanin1    {"", prompt = "Internal use only"}

begin
    char    l_inimages = ""
    char    l_outimages = ""
    char    l_outprefix = ""
    char    l_process = ""
    int     l_dispaxis
    char    l_logfile = ""
    bool    l_verbose

    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_pixscale = ""
    char    l_key_dispaxis = ""

    int     junk, nfiles, nx, ny, nz, iextn, i, axis
    bool    debug, havemdf, rotate, resample, havedq
    struct  sdate, sline, sline2
    char    tmpout, tmplist, inimg, outimg, tmpcube
    char    tmpslices, tmpslice, tmpsci, tmpvar, tmpdq, tmproot, tmptrans
    char    tmpscale
    char    files, extns, badhdr, phu, mdf, radec
    char    extn[3]
    real    crpix1, crpix2, crpix3, crval1, crval2, crval3
    real    cd11, cd12, cd21, cd22, cd33
    real    ltm22, ltv2, ltm33, ltv3
    real    fit, sigma, pixscale, apgap
    char    imgsecn, interp
    int     ny1, ny2, pccount, pctotal, trnx, trny, blank
    real    yloas, yhias, yloas2, yhias2, shift, scale, ylo, yhi
    real    x1, x2, y1, y2, dqscale

    junk = fscan (inimages, l_inimages)
    junk = fscan (outimages, l_outimages)
    junk = fscan (outprefix, l_outprefix)
    junk = fscan ( process, l_process)
    l_dispaxis = dispaxis
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    status = 1
    debug = no
    havemdf = no
    pctotal = 10
    dqscale = 1000

    l_process = strlwr (l_process)
    resample = "resample" == l_process
    rotate = resample || "rotate" == l_process

    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"

    extn[1] = l_sci_ext
    extn[2] = l_var_ext
    extn[3] = l_dq_ext

    tmpout = mktemp ("tmpout")
    tmproot = mktemp ("tmproot")
    tmplist = mktemp ("tmplist")
    tmpslice = mktemp ("tmpslice")
    tmpslices = mktemp ("tmpslices")
    tmpsci = mktemp ("tmpsci")
    tmpdq = mktemp ("tmpdq")
    tmpvar = mktemp ("tmpvar")
    tmpcube = mktemp ("tmpcube")


    cache ("gemextn", "wmef", "gemdate")


    # Start logging

    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile) 
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NFCUBE: Both nfcube.logfile and \
                gnirs.logfile are empty.", l_logfile, verbose+) 
            printlog ("                   Using default file " \
                // l_logfile // ".", l_logfile, verbose+) 
        }
    }
    date | scan (sdate) 
    printlog ("--------------------------------------------------------------\
        ----------------", l_logfile, verbose = l_verbose) 
    printlog ("NFCUBE -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 


    # Check input data

    if (debug) print ("checking input")
    gemextn (l_inimages, proc="none", check="", index="", extname="", 
        extver="", ikparams="", omit="path", replace="", outfile=tmproot, 
        logfile="", glogpars="", verbose-)
    if (0 < gemextn.fail_count || 0 == gemextn.count) {
        printlog ("ERROR - NFCUBE: Missing input data.", \
            l_logfile, verbose+) 
        goto clean
    }
    nfiles = gemextn.count


    # Check/generate output names

    if (debug) print ("checking output")
    gemextn (l_outimages, proc="none", check="absent", index="", \
        extname="", extver="", ikparams="", omit="", replace="", \
        outfile=tmpout, logfile="", glogpars="", verbose-)
    if (0 < gemextn.fail_count) {
        printlog ("ERROR - NFCUBE: Output data already exist.", \
            l_logfile, verbose+) 
        goto clean
    }
    if (0 != gemextn.count) {
        if (nfiles != gemextn.count) {
            printlog ("ERROR - NFCUBE: Incorrect number of output \
                files.", l_logfile, verbose+) 
            goto clean
        }
    } else {
        delete (tmpout, verify-, >& "dev$null")
        gemextn ("%^%" // l_outprefix // "%" // "@" // tmproot,
            proc="none", check="absent", index="", extname="", \
            extver="", ikparams="", omit="kernel", replace="", \
            outfile=tmpout, logfile="", glogpars="", verbose-)
        if (0 < gemextn.fail_count || nfiles != gemextn.count) {
            printlog ("ERROR - NFCUBE: Bad output files.", \
                l_logfile, verbose+) 
            goto clean
        }
    }


    # Combine in + out into single list

    joinlines (tmproot, tmpout, output=tmplist, delim=" ", missing="Missing", \
        maxchar=161, shortest=yes, verbose-)

    # Generate the cubes

    scanin1 = tmplist

    while (fscan (scanin1, inimg, outimg) != EOF) {

        if (debug) print (inimg)

        delete (tmpslices, verify-, >& "dev$null")
        tmpslices = mktemp ("tmpslices")
        imdelete (tmpsci, verify-, >& "dev$null")
        tmpsci = mktemp ("tmpsci")

        gemextn (inimg, check="exists", proc="expand", index="", \
            extname=l_sci_ext, extversion="1-", ikparams="", \
            omit="", replace="", outfile=tmpslices, logfile=l_logfile, \
            glogpars="", verbose=l_verbose)

        if (0 == gemextn.count || 0 < gemextn.fail_count) {
            printlog ("ERROR - NFCUBE: Problem with " // inimg, \
                l_logfile, verbose+)
            goto clean
        }

        if (debug) type (tmpslices)
        imstack ("@" // tmpslices, tmpsci, title="*", pixtype="*")
        #!!! no !!! # Invert y order
        #imcopy (tmpsci // "[*,-*,*]", tmpsci, verbose-, >& "dev$null")
        files = tmpsci
        extns = l_sci_ext

        delete (tmpslices, verify-, >& "dev$null")
        tmpslices = mktemp ("tmpslices")
        imdelete (tmpvar, verify-, >& "dev$null")
        tmpvar = mktemp ("tmpvar")

        gemextn (inimg, check="exists", proc="expand", index="", \
            extname=l_var_ext, extversion="1-", ikparams="", \
            omit="", replace="", outfile=tmpslices, logfile=l_logfile, \
            glogpars="", verbose=l_verbose)

        if (0 < gemextn.count && 0 == gemextn.fail_count) {
            imstack ("@" // tmpslices, tmpvar, title="*", pixtype="*")
            #!!! no !!! # Invert y order
            #imcopy (tmpvar // "[*,-*,*]", tmpvar, verbose-, >& "dev$null")
            files = files // "," // tmpvar
            extns = extns // "," // l_var_ext
        } else if (0 < gemextn.fail_count) {
            printlog ("ERROR - NFCUBE: Problem with " // inimg, \
                l_logfile, verbose+)
            goto clean
        }

        delete (tmpslices, verify-, >& "dev$null")
        tmpslices = mktemp ("tmpslices")
        imdelete (tmpdq, verify-, >& "dev$null")
        tmpdq = mktemp ("tmpdq")

        gemextn (inimg, check="exists", proc="expand", index="", \
            extname=l_dq_ext, extversion="1-", ikparams="", \
            omit="", replace="", outfile=tmpslices, logfile=l_logfile, \
            glogpars="", verbose=l_verbose)

        if (0 < gemextn.count && 0 == gemextn.fail_count) {
            imstack ("@" // tmpslices, tmpdq, title="*", pixtype="*")
            #!!! no !!! # Invert y order
            #imcopy (tmpdq // "[*,-*,*]", tmpdq, verbose-, >& "dev$null")
            files = files // "," // tmpdq
            extns = extns // "," // l_dq_ext
        } else if (0 < gemextn.fail_count) {
            printlog ("ERROR - NFCUBE: Problem with " // inimg, \
                l_logfile, verbose+)
            goto clean
        }

        gemextn (inimg, check="exists", proc="expand", index="-", \
            extname="MDF", extversion="", ikparams="", \
            omit="name", replace="", outfile="STDOUT", logfile=l_logfile, \
            glogpars="", verbose-) | scan (mdf)

        if (1 == gemextn.count && 0 == gemextn.fail_count) {
            havemdf = yes
        } else if (0 < gemextn.fail_count) {
            printlog ("ERROR - NFCUBE: Problem with " // inimg, \
                l_logfile, verbose+)
            goto clean
        } else {
            printlog ("WARNING - NFCUBE: No MDF in " // inimg, \
                l_logfile, verbose+)
        }

        if (debug) {
            wmef (input=files, output=outimg, extnames=extns, \
                phu=inimg//"[0]", verbose+)
        } else {
            wmef (input=files, output=outimg, extnames=extns, \
                phu=inimg//"[0]", verbose-, >& "dev$null")
        }
        if (0 != wmef.status) {
            printlog ("ERROR - NFCUBE: Problem constructing MEF.", \
                l_logfile, verbose+)
            goto clean
        }

        if (havemdf) {
            fxinsert (mdf, outimg // ".fits[0]", "", verbose-)
        }


        # Read info to construct WCS

        phu = inimg // "[0]"

        pixscale = INDEF
        hselect (phu, l_key_pixscale, yes) | scan (pixscale)
        if (isindef (pixscale)) {
            hselect (phu, "CD1_1", yes) | scan (pixscale)
            if (isindef (pixscale)) { 
                printlog ("ERROR - NFCUBE: No pixel scale for " // inimg, \
                    l_logfile, verbose+)
                goto clean
            } else {
                printlog ("WARNING - NFCUBE: No pixel scale for " \
                    // inimg, l_logfile, verbose+)
                printlog ("                  so using CD1_1: " \
                    // pixscale, l_logfile, verbose+)
            }
        }

        axis = l_dispaxis
        hselect (phu, l_key_dispaxis, yes) | scan (axis)
        if (isindef (axis)) {
            printlog ("ERROR - NFCUBE: No dispersion axis for " // inimg, \
                l_logfile, verbose+)
            goto clean
        }

        gemextn (inimg, proc="expand", check="", index="", \
            extname=l_sci_ext, extver="1-", ikparams="", \
            omit="", replace="", outfile="STDOUT", logfile="dev$null", \
            glogpars="", verbose-) | scan (sline)

        # Note axis labelling (1-3, x-y) in local variables is for
        # ROTATED cube

        hselect (sline, "CD"//axis//"_"//axis, yes) | scan (cd33)
        hselect (sline, "CRVAL"//axis, yes) | scan (crval3)
        hselect (sline, "CRPIX"//axis, yes) | scan (crpix3)

        gemextn (outimg, proc="expand", check="", index="", \
            extname=l_sci_ext, extver="1-", ikparams="", \
            omit="", replace="", outfile="STDOUT", logfile="dev$null", \
            glogpars="", verbose-) | scan (sline)

        if (2 == axis) {
            hselect (sline, "i_naxis1", yes) | scan (nx)
            hselect (sline, "i_naxis2", yes) | scan (nz)
            hselect (sline, "i_naxis3", yes) | scan (ny)
        } else {
            hselect (sline, "i_naxis1", yes) | scan (nz)
            hselect (sline, "i_naxis2", yes) | scan (nx)
            hselect (sline, "i_naxis3", yes) | scan (ny)
        }

        # Instead of original phu, use central origin

        # Need y spacing from MDF

        gemextn (inimg, proc="expand", check="", index="", \
            extname="MDF", extver="", ikparams="", \
            omit="", replace="%\[%.fits[%", outfile="STDOUT", \
            logfile="dev$null", glogpars="", verbose-) | scan (sline)
        if (0 == gemextn.count || no == (0 == gemextn.fail_count)) {
            printlog ("WARNING - NFCUBE: No MDF, so assuming square \
                pixels for " // inimg, l_logfile, verbose+)
            apgap = pixscale
        } else {
            if (2 == axis) radec = "DEC"
            else           radec = "RA"
            if (debug) {
                print ("fitting to aperture spacing in " // sline)
                tlinear (sline, outtable="STDOUT", xcol="row", \
                    ycol=radec, wcol="", scol="", rows="-", outcoly="", \
                    outcolr="")
            }
            fit = INDEF; sigma = INDEF
            tlinear (sline, outtable="STDOUT", xcol="row", ycol=radec, \
                wcol="", scol="", rows="-", outcoly="", outcolr="") \ 
                | fields ("STDIN", "5", lines="5", quit-, print-) \
                | scan (fit)
            tlinear (sline, outtable="STDOUT", xcol="row", ycol=radec, \
                wcol="", scol="", rows="-", outcoly="", outcolr="") \ 
                | fields ("STDIN", "3", lines="10", quit-, print-) \
                | scan (sigma)
            if (debug) print (fit)
            if (debug) print (sigma)
            if (isindef (fit) || isindef (sigma)) {
                printlog ("WARNING - NFCUBE: Fit to MDF failed, so \
                    square pixels for " // inimg, l_logfile, verbose+)
                apgap = pixscale
            } else if (abs (fit / sigma) < 10.0) {
                printlog ("WARNING - NFCUBE: Fit to MDF poor, so \
                    square pixels for " // inimg, l_logfile, verbose+)
                apgap = pixscale
            } else {
                # HMS to DMS for RA ?
                if (1 == axis) fit = fit * 15 
                apgap = -1.0 * fit * 3600.0
                printf ("NFCUBE: Aperture sepn of %5.3f arcsec for %s\n", \
                    apgap, outimg) | scan (sline)
                printlog (sline, l_logfile, l_verbose)
            }
        }

        cd11 = pixscale
        cd22 = apgap
        cd12 = 0.
        cd21 = 0.
        crval1 = 0.
        crpix1 = (nx + 1) / 2.0
        crval2 = 0.
        crpix2 = (ny + 1) / 2.0


        if (rotate) {

            for (iextn = 1; iextn <= 3; iextn = iextn + 1) {

                gemextn (outimg, proc="expand", check="", index="", \
                    extname=extn[iextn], extver="1", ikparams="", \
                    omit="", replace="", outfile="STDOUT", \
                    logfile="dev$null", glogpars="", verbose-) \
                    | scan (sline)
                if (debug) print (iextn // " " // sline)

                if (1 == gemextn.count) {

                    printlog ("NFCUBE - Please wait, rotating " \
                        // sline, l_logfile, l_verbose)

                    havedq = l_dq_ext == extn[iextn]
                    if (havedq) {
                        interp = "linear"
                        blank = dqscale
                    } else {
                        interp = "spline3"
                        blank = 0
                    }

                    delete (tmpslices, verify-, >& "dev$null")
                    tmpslices = mktemp ("tmpslices")

                    if (resample) {

                        # thinking just about gnirs for the moment
                        # we need to resample along the y axis.
                        # the number of pixels you need if pixels
                        # are all the same size is
                        # ny' = ny * apgap / pixscale
                        # but that won't be an integer, so make it 
                        # ny' = int (ny * apgap / pixscale) + 1

                        # now what about y values themselves?
                        # currently the centre of the extreme pixels
                        # in arcsecs are
                        # yloas = apgap / 2.0
                        # yhias = ny * apgap - apgap / 2.0
                        # and, with the new pixels
                        # yloas' = pixscale / 2.0
                        # yhias' = ny' * pixscale - pixscale / 2.0

                        # these aren't the same for two reasons -
                        # pixel size (with coords at pixel centres)
                        # and the extra extension to get an integer
                        # number of pixels.

                        # what we want as args to transorm, i think,
                        # are the coords of the points identified by
                        # yloas' and thias' on the sky, in terms of
                        # pixels in the current image.

                        # this means a shift in origin and a scale
                        # ylo = (yloas' - shift) * scale
                        # yhi = (yhias' - shift) * scale

                        # now we already know that
                        # 1 = (yloas - shift) * scale
                        # ny = (yhias - shift) * scale
                        # so (ratio)
                        # ny * (yloas - shift) = yhias - shift
                        # ny * yloas - yhias = shift * (ny - 1)
                        # shift = (ny * yloas - yhias) / (ny - 1)
                        # scale = 1 / (yloas - shift)
                        # yloas - shift =
                        #  (yloas * (ny - 1) - ny * yloas + yhias) 
                        #   / (ny - 1)
                        #  (yhias - yloas) / (ny - 1)
                        # scale = (ny - 1) / (yhias - yloas)

                        # (note that ny is actually nz...)

                        ny1 = ny
                        ny2 = int (ny1 * apgap / pixscale) + 1
                        yloas = apgap / 2.0
                        yhias = ny1 * apgap - apgap / 2.0
                        yloas2 = pixscale / 2.0
                        yhias2 = ny2 * pixscale - pixscale / 2.0
                        shift = (ny1 * yloas - yhias) / real (ny1 - 1)
                        scale = real (ny1 - 1) / (yhias - yloas)

                        ylo = (yloas2 - shift) * scale
                        yhi = (yhias2 - shift) * scale

                        if (debug) {
                            print ("resampling:")
                            print ("gap/scale: " // apgap // "/" \
                                // pixscale)
                            print ("ny/ny': " // ny1 // "/" // ny2)
                            print ("ylo, yhi (arcsec): " // yloas \
                                // ", " // yhias)
                            print ("ylo', yhi' (arcsec): " // yloas2 \
                                // ", " // yhias2)
                            print ("shift/scale: " // shift // "/" \
                                // scale)
                            print ("ylo, yhi (pixels): " // ylo \
                                // ", " // yhi)
                        }

                        x1 = INDEF; x2 = INDEF; trnx = INDEF
                        y1 = ylo; y2 = yhi; trny = ny2
                    }

                    pccount = 1

                    for (i = 1; i <= nz; i = i + 1) {

                        if (l_verbose) {
                            if (i / real (nz) > \
                                real (pccount) / real (pctotal) || i == nz) {
                                
                                pccount = pccount + 1
                                printf ("NFCUBE - Rotation %3.0f%% \
                                    complete\n", \
                                    (real (100*i) / real (nz))) \
                                    | scan (sline2)
                                printlog (sline2, l_logfile, l_verbose)
                            }
                        }

                        tmpslice = mktemp ("tmpslice")
                        print (tmpslice, >> tmpslices)
                        if (2 == axis) {
                            # flip y to agree with acq image
                            imgsecn = sline // "[*," // i // ",-*]"
                        } else {
                            # may need -ve sign somewhere....
                            imgsecn = sline // "[" // i // ",*,*]"
                        }

                        # This flip introduces WAXMAP01 and messes up
                        # with the LTM3_3 and LTV3.  Later, when the slices
                        # are stacked, and the cube is rotated in the
                        # process, the LT?3 values are not automatically 
                        # switched with the LT?2 values. This will need to
                        # be corrected below (once we have the final cube).

                        imcopy (imgsecn, tmpslice, >& "dev$null")

                        if (resample) {
                            gemhedit (tmpslice, "CTYPE1", "LINEAR", "", 
                                delete-)
                            gemhedit (tmpslice, "CTYPE2", "LINEAR", "", 
                                delete-)
                            gemhedit (tmpslice, "CRPIX1", "1", "", delete-)
                            gemhedit (tmpslice, "CRPIX2", "1", "", delete-) 
                            gemhedit (tmpslice, "CRVAL1", "1", "", delete-)
                            gemhedit (tmpslice, "CRVAL2", "1", "", delete-)
                            gemhedit (tmpslice, "CD1_1", "1", "", delete-)
                            gemhedit (tmpslice, "CD1_2", "0", "", delete-)
                            gemhedit (tmpslice, "CD2_1", "0", "", delete-)
                            gemhedit (tmpslice, "CD2_2", "1", "", delete-)

                            if (havedq) {
                                tmpscale = mktemp ("tmpscale")
                                imarith (tmpslice, "*", dqscale, \
                                    tmpscale, title="", divzero=0, \
                                    hparams="", pixtype="", calctype="", \
                                    verbose-, noact-)
                                imdelete (tmpslice, verify-, >& "dev$null")
                                imrename (tmpscale, tmpslice)
                            }
                                
                            tmptrans = mktemp ("tmptrans")
                            transform (tmpslice, tmptrans, fitnames="", \
                                minput="", moutput="", database="", \
                                interptype=interp, x1=x1, x2=x2, \
                                dx=INDEF, nx=trnx, xlog-, \
                                y1=y1, y2=y2, dy=INDEF, ny=trny, ylog-, \
                                flux+, blank=blank, logfiles="")
                            imdelete (tmpslice, verify-, >& "dev$null")
                            imrename (tmptrans, tmpslice)

                            if (havedq) {
                                tmpscale = mktemp ("tmpscale")
                                imexpr ("a > 0 ? 1 : 0", tmpscale, \
                                    tmpslice, dims="auto", \
                                    intype="auto", outtype="int", \
                                    refim="auto", bwidth=0, \
                                    btype="nearest", bpixval=0, \
                                    rangecheck+, verbose-, \
                                    exprdb="none")
                                imdelete (tmpslice, verify-, >& "dev$null")
                                imrename (tmpscale, tmpslice)
                            }
                        }
                    }


                    imdelete (tmpcube, verify-, >& "dev$null")
                    tmpcube = mktemp ("tmpcube")
                    imstack ("@" // tmpslices, tmpcube, title="*", \
                        pixtype="*")
                    imdelete ("@" // tmpslices, verify-, >& "dev$null")
                    # DISPAXIS is now 3 (edited below)

                    gemextn (outimg, proc="expand", check="", index="", \
                        extname=extn[iextn], extver="1", \
                        ikparams="overwrite", \
                        omit="", replace="", outfile="STDOUT", \
                        logfile=l_logfile, glogpars="", verbose-) \
                        | scan (sline)
                    imcopy (tmpcube, sline, >& "dev$null")

                    # Update world and logical coordinate systems
                    gemhedit (sline, "DISPAXIS", 3, "", delete-)

                    gemhedit (sline, "WCSDIM", 3, "It's a cube", delete-)
                    gemhedit (sline, "WAT1_001", "wtype=linear", "", delete-)
                    gemhedit (sline, "WAT2_001", "wtype=linear", "", delete-)
                    gemhedit (sline, "WAT3_001", 
                        "wtype=linear label=Wavelength units=angstrom",
                        comment="", delete-)

                    gemhedit (sline, "CTYPE1", "LINEAR", "All axes are linear",
                        delete-)
                    gemhedit (sline, "CTYPE2", "LINEAR", "All axes are linear",
                        delete-)
                    gemhedit (sline, "CTYPE3", "LINEAR", "All axes are linear",
                        delete-)

                    gemhedit (sline, "CRPIX1", crpix1,
                        "Reference pixel for x coord", delete-)
                    gemhedit (sline, "CRPIX2", crpix2, 
                        "Reference pixel for y coord", delete-)
                    gemhedit (sline, "CRPIX3", crpix3, 
                        "Reference pixel for x coord", delete-)

                    gemhedit (sline, "CRVAL1", crval1, 
                        "Reference value for x coord", delete-)
                    gemhedit (sline, "CRVAL2", crval2, 
                        "Reference value for y coord", delete-)
                    gemhedit (sline, "CRVAL3", crval3, 
                        "Reference value for z coord", delete-)

                    gemhedit (sline, "CD1_1", cd11, 
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD1_2", cd12,
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD1_3", 0., 
                        "Coordinate transform matrix", delete-)

                    gemhedit (sline, "CD2_1", cd21, 
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD2_2", cd22, 
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD2_3", 0., 
                        "Coordinate transform matrix", delete-)

                    gemhedit (sline, "CD3_1", 0., 
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD3_2", 0., 
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD3_3", cd33, 
                        "Coordinate transform matrix", delete-)

                    # Remove WAXMAP01 and fix LTM and LTV values
                    # (Problems caused by the axis flip to match the 
                    #  acquisition image and subsequent stacking of
                    #  rotated slices.)

                    gemhedit (sline, "WAXMAP01", "", "", delete=yes)

                    ltm22 = 1. ; ltm33 = 1. ; ltv2 = 0. ; ltv3 = 0.
                    hselect (sline, "LTM3_3", yes) | scan (ltm22)
                    hselect (sline, "LTV3", yes) | scan (ltv2)
                    hselect (sline, "LTM2_2", yes) | scan (ltm33)
                    hselect (sline, "LTV2", yes) | scan (ltv3)
                    gemhedit (sline, "LTM2_2", ltm22, "", delete-)
                    gemhedit (sline, "LTV2", ltv2, "", delete-)
                    gemhedit (sline, "LTM3_3", ltm33, "", delete-)
                    gemhedit (sline, "LTV3", ltv3, "", delete-)

                    # Delete WCS keywords for axis 4, if they exist
                    gemhedit (sline, "CD4_4", "", "", delete=yes)
                    gemhedit (sline, "CTYPE4", "", "", delete=yes)
                    gemhedit (sline, "LTM4_4", "", "", delete=yes)
                    gemhedit (sline, "WAT4_001", "", "", delete=yes)

                }
            }

        } else {

            # Generate WCS for unrotated data
            # This is a verbose mess because of known limitations 
            # with mkheader

            for (iextn = 1; iextn <= 3; iextn = iextn + 1) {

                gemextn (outimg, proc="expand", check="", index="", \
                    extname=extn[iextn], extver="1", ikparams="", \
                    omit="", replace="", outfile="STDOUT", \
                    logfile="dev$null", glogpars="", verbose-) \
                    | scan (sline)
                if (debug) print (iextn // " " // sline)

                if (1 == gemextn.count) {

                    gemhedit (sline, "WCSDIM", 3, "It's a cube", delete-)
                    gemhedit (sline, "WAT1_001", "wtype=linear", "", delete-)
                    gemhedit (sline, "WAT2_001", 
                        "wtype=linear label=Wavelength units=angstrom",
                        comment="", delete-)
                    gemhedit (sline, "WAT3_001", "wtype=linear", "", delete-)
                        
                    gemhedit (sline, "CTYPE1", "LINEAR", 
                        "All axes are linear", delete-)
                    gemhedit (sline, "CTYPE2", "LINEAR", 
                        "All axes are linear", delete-)
                    gemhedit (sline, "CTYPE3", "LINEAR", 
                        "All axes are linear", delete-)

                    gemhedit (sline, "CRPIX1", crpix1, 
                        "Reference pixel for x coord", delete-)
                    gemhedit (sline, "CRPIX2", crpix3, 
                        "Reference pixel for y coord", delete-)
                    gemhedit (sline, "CRPIX3", crpix2, 
                        "Reference pixel for z coord", delete-)

                    gemhedit (sline, "CRVAL1", crval1, 
                        "Reference value for x coord", delete-)
                    gemhedit (sline, "CRVAL2", crval3, 
                        "Reference value for y coord", delete-)
                    gemhedit (sline, "CRVAL3", crval2, 
                        "Reference value for z coord", delete-)

                    gemhedit (sline, "CD1_1", cd11, 
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD1_2", 0., 
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD1_3", cd12, 
                        "Coordinate transform matrix", delete-)

                    gemhedit (sline, "CD2_1", 0., 
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD2_2", cd33, 
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD2_3", 0., 
                        "Coordinate transform matrix", delete-)

                    gemhedit (sline, "CD3_1", cd21, 
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD3_2", 0., 
                        "Coordinate transform matrix", delete-)
                    gemhedit (sline, "CD3_3", cd22, 
                        "Coordinate transform matrix", delete-)

                    # Delete WCS keywords for axis 4, if they exist
                    gemhedit (sline, "CD4_4", "", "", delete=yes)
                    gemhedit (sline, "CTYPE4", "", "", delete=yes)
                    gemhedit (sline, "LTM4_4", "", "", delete=yes)
                    gemhedit (sline, "WAT4_001", "", "", delete=yes)

                }
            }
        }
        gemdate ()
        gemhedit (outimg//"[0]", "NFCUBE", gemdate.outdate,
            "UT Time stamp for NFCUBE", delete-)
        gemhedit (outimg//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI IRAF", delete-)
    }


    status = 0 # Success

clean:

    scanin1 = ""

    delete (tmpout, verify-, >& "dev$null")
    delete (tmplist, verify-, >& "dev$null")
    imdelete (tmpsci, verify-, >& "dev$null")
    imdelete (tmpvar, verify-, >& "dev$null")
    imdelete (tmpdq, verify-, >& "dev$null")
    imdelete (tmpcube, verify-, >& "dev$null")
    delete (tmpslices, verify-, >& "dev$null")
    delete (tmproot, verify-, >& "dev$null")

    printlog (" ", l_logfile, l_verbose) 
    if (0 == status)
        printlog ("NFCUBE exit status: good.", l_logfile, l_verbose) 
    else
        printlog ("NFCUBE exit status: FAILED.", l_logfile, verbose+)

end
