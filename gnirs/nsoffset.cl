# Copyright(c) 2004-2012 Association of Universities for Research in Astronomy, Inc.

# Initial version: acooke 27 aug 2004

# measure the shift of one image relative to another
# *plus* cool (trez cool, mon cherry) extra features like accepting 
# mdfs and updating wcs values

procedure nsoffset (refimage, inimages)

char    refimage    {prompt = "Reference image"}
char    inimages    {prompt = "Images to shift"}
char    outimages   {"", prompt = "Output files"}
char    outprefix   {"o", prompt = "Prefix to use if outimages not given"}
real    shift       {0., prompt = "Measured/required shift"}
int     axis        {INDEF, prompt = "Shift axis"}
real    pixscale    {INDEF, prompt = "Pixel scale (if not in header)"}
bool    fl_measure  {yes, prompt = "Measure shift?"}
bool    fl_apply    {no, prompt = "Apply shift to input images?"}
bool    fl_all      {yes, prompt = "Measure all input images?"}
bool    fl_project  {no, prompt = "Project along dispersion direction?"}
bool    fl_inter    {no, prompt = "Fit peak interactively?"}
char    logfile     {"", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose output?"}

int     status      {0, prompt = "Exit status (0 = good)"}
struct  *scanin1    {"", prompt = "Internal use only"}
struct  *scanin2    {"", prompt = "Internal use only"}

begin
    char    l_refimage = ""
    char    l_inimages = ""
    char    l_outimages = ""
    char    l_outprefix = ""
    real    l_shift
    int     l_axis
    real    l_pixscale
    bool    l_fl_measure
    bool    l_fl_apply
    bool    l_fl_all
    bool    l_fl_project
    bool    l_fl_inter
    char    l_logfile = ""
    bool    l_verbose

    char    l_sci_ext = ""
    char    l_dq_ext = ""
    char    l_key_dispaxis = ""
    char    l_key_pixscale = ""

    int     junk, dispaxis, nfiles, nsecns, isecn, slitid
    int     spatialaxis, fitwidth, dispwidth
    int     nx, ny, fnx, fny
    int     xc, yc, midpix, i, imax
    bool    debug, first, havemdf
    bool    imgmodified
    struct  sline
    char    badhdr, refsci, imgin, imgout, imgsci, mdfin, refphu
    char    slittype, target, defkey, sjunk, mdfout, imgdq, refdq
    char    spec, graphics, word
    char    tmpref, tmprefsci, tmproot, tmpout, tmplist, tmpin
    char    tmpmdfimg, tmparith, tmptapref, tmptaptag, tmpconv
    char    tmpcursor, tmplog1, tmplog2, tmptab, tmpdq, tmpblank, tmpmask
    char    tmpprjtag, tmpprjref, template
    char    cjunk
    real    scale, x1, x2, y1, y2, corner
    real    cmin, cmax
    real    crpix, border, blur, rjunk

    junk = fscan (  refimage, l_refimage)
    junk = fscan (  inimages, l_inimages)
    junk = fscan (  outimages, l_outimages)
    junk = fscan (  outprefix, l_outprefix)
    l_shift =       0.1 * nint (10 * shift) # restrict to 0.1 resoln
    l_axis =        axis
    l_pixscale =    pixscale
    l_fl_measure =  fl_measure
    l_fl_apply =    fl_apply
    l_fl_all =      fl_all
    l_fl_project =  fl_project
    l_fl_inter =    fl_inter
    junk = fscan (  logfile, l_logfile)
    l_verbose =     verbose

    status = 1
    scale = INDEF
    debug = no
    badhdr = ""
    nx = INDEF
    ny = INDEF
    if (l_fl_measure)
        l_shift = INDEF

    # these chosen for ifu mdf, but should work with anything larger
    blur = 2
    border = -5

    # for interactive fitting/plotting
    fitwidth = 3
    dispwidth = 5 * fitwidth

    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"

    cache ("gemextn", "minmax", "gemdate")

    tmpref =    mktemp ("tmpref")
    tmprefsci = mktemp ("tmprefsci")
    tmproot =   mktemp ("tmproot")
    tmpout =    mktemp ("tmpout")
    tmplist =   mktemp ("tmplist")
    tmpin =     mktemp ("tmpin")
    tmpdq =     mktemp ("tmpdq")
    tmpblank =  mktemp ("tmpblank")
    tmpmask =   mktemp ("tmpmask")
    tmpmdfimg = mktemp ("tmpmdfimg")
    tmparith =  mktemp ("tmparith")
    tmptapref = mktemp ("tmptapref")
    tmptaptag = mktemp ("tmptaptag")
    tmpconv =   mktemp ("tmpconv")
    tmpcursor = mktemp ("tmpcursor")
    tmplog1 =   mktemp ("tmplog1")
    tmplog2 =   mktemp ("tmplog2")
    tmptab =    mktemp ("tmptab")
    tmpprjtag = mktemp ("tmpprjtag")
    tmpprjref = mktemp ("tmpprjref")

    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile) 
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSOFFSET: Both nsoffset.logfile and \
                gnirs.logfile are empty.", l_logfile, verbose+) 
            printlog ("                   Using default file " \
                // l_logfile // ".", l_logfile, verbose+) 
        }
    }

    date | scan (sline) 
    printlog ("---------------------------------------------------------\
        ---------------------", l_logfile, verbose = l_verbose) 
    printlog ("NSOFFSET -- " // sline, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NSOFFSET: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }

    if (debug) print ("checking logic")
    if (no == l_fl_measure && no == l_fl_apply) {
        printlog ("ERROR - NSOFFSET: Neither measuring nor applying \
            a shift.", l_logfile, verbose+) 
        goto clean
    }
    if (no == l_fl_measure && l_fl_all) {
        printlog ("ERROR - NSOFFSET: Cannot measure all shifts if \
            measuring none", l_logfile, verbose+) 
        printlog ("                  (fl_all+ but fl_measure-).", \
            l_logfile, verbose+) 
        goto clean
    }

    if (debug) print ("checking reference")
    gemextn (l_refimage, proc="none", check="exists", index="", \
        extname="", extver="", ikparams="", omit="", replace="", \
        outfile=tmpref, logfile="", glogpars="", verbose-)
    if (0 < gemextn.fail_count || 1 != gemextn.count) {
        printlog ("ERROR - NSOFFSET: Error in reference image.", \
            l_logfile, verbose+) 
        goto clean
    }
    scanin1 = tmpref; junk = fscan (scanin1, l_refimage)
    refphu = l_refimage // "[0]"

    if (debug) print ("getting axis")
    spatialaxis = l_axis
    if (debug) print ("axis param: " // spatialaxis)
    dispaxis = INDEF
    hselect (refphu, l_key_dispaxis, yes) | scan (dispaxis)
    if (no == isindef (dispaxis))
        spatialaxis = 3 - dispaxis
    if (isindef (spatialaxis)) {
        printlog ("ERROR - NSOFFSET: No axis parameter value and \
            no " // l_key_dispaxis, l_logfile, verbose+) 
        printlog ("                  in " // refphu, \
            l_logfile, verbose+) 
        goto clean
    }
    if (debug) print ("axis: " // spatialaxis)

    if (l_fl_measure) {

        if (debug) print ("getting first sci extension from ref")
        gemextn (l_refimage, proc="expand", check="exists", index="", \
            extname=l_sci_ext, extver="1-", ikparams="", omit="", \
            replace="", outfile=tmprefsci, logfile="", glogpars="", verbose-)
        if (0 < gemextn.fail_count || 0 == gemextn.count) {
            printlog ("ERROR - NSOFFSET: Missing data in reference.", \
                l_logfile, verbose+) 
            goto clean
        }
        scanin1 = tmprefsci; junk = fscan (scanin1, refsci)
        template = refsci

        if (debug) print ("getting size of reference image")
        hselect (refsci, "i_naxis1", yes) | scan (nx)
        hselect (refsci, "i_naxis2", yes) | scan (ny)
        if (isindef (nx) || isindef (ny)) {
            printlog ("ERROR - NSOFFSET: Cannot read size of " // refsci \
                // ".", l_logfile, verbose+) 
            goto clean
        }
        if (debug) print (nx // "x" // ny)

        if (debug) print ("zeroing bad pixels")
        delete (tmpdq, verify-, >& "dev$null")
        gemextn (l_refimage, proc="expand", check="exists", index="", \
            extname=l_dq_ext, extver="1-", ikparams="", omit="", \
            replace="", outfile=tmpdq, logfile="dev$null", glogpars="",
            verbose-)
        if (0 == gemextn.fail_count && 0 != gemextn.count) {
            scanin2 = tmpdq
            junk = fscan (scanin2, refdq)
            imdelete (tmpmask, verify-, >& "dev$null")
            imexpr ("a>0?0:b", tmpmask, refdq, refsci, \
                dims="auto", intype="auto", outtype="auto", \
                refim="auto", bwidth=0, btype="nearest", \
                bpixval=0, rangecheck+, verbose=debug)
            imdelete (refsci, verify-, >& "dev$null")
            refsci = tmpmask
            if (debug) {
                print ("updated refsci from " // refdq)
                imhead (refsci)
            }
        }

        if (debug) print ("clipping at zero")
        imdelete (tmparith, verify-, >& "dev$null")
        imexpr ("max(0,a)", tmparith, refsci, \
            dims="auto", intype="auto", outtype="auto", \
            refim="auto", bwidth=0, btype="nearest", \
            bpixval=0, rangecheck+, verbose=debug)
        imdelete (tmpmask // "," // refsci, verify-, >& "dev$null")
        tmpmask = tmparith; refsci = tmpmask
        if (debug) imhead (refsci)

        if (l_fl_project) {

            if (debug) print ("projecting reference")
            improject (refsci, tmpprjref, projaxis=(3-spatialaxis), \
                average+, highcut=0, lowcut=0, pixtype="real", verbose-)
            imdelete (refsci, verify-, >& "dev$null")
            refsci = tmpprjref
        }

        if (debug) print ("tapering reference")
        taper (refsci, tmptapref, width="10 %", subtract="none", \
            function="cosbell", verbose-, >& "dev$null")
        if (debug) {
            if (l_fl_project)
                splot (tmptapref)
            else
                display (tmptapref, 2)
        }
    }

    if (debug) print ("checking input")
    gemextn (l_inimages, proc="none", check="exists", index="", \
        extname="", extver="", ikparams="", omit="exten,kernel,index", \
        replace="", outfile=tmproot, logfile="", glogpars="", verbose-)
    if (0 < gemextn.fail_count || 0 == gemextn.count) {
        printlog ("ERROR - NSOFFSET: Missing input data.", \
            l_logfile, verbose+) 
        goto clean
    }
    nfiles = gemextn.count

    if (l_fl_apply) {

        if (debug) print ("checking output")
        gemextn (l_outimages, proc="none", check="absent", index="", \
            extname="", extver="", ikparams="", omit="", replace="", \
            outfile=tmpout, logfile="", glogpars="", verbose-)
        if (0 < gemextn.fail_count) {
            printlog ("ERROR - NSOFFSET: Output data already exist.", \
                l_logfile, verbose+) 
            goto clean
        }
        if (0 != gemextn.count) {
            if (nfiles != gemextn.count) {
                printlog ("ERROR - NSOFFSET: Incorrect number of output \
                    files.", l_logfile, verbose+) 
                status = 2
                goto clean
            }
        } else {
            delete (tmpout, verify-, >& "dev$null")
            gemextn ("%^%" // l_outprefix // "%" // "@" // tmproot,
                proc="none", check="absent", index="", extname="", \
                extver="", ikparams="", omit="", replace="", \
                outfile=tmpout, logfile="", glogpars="", verbose-)
            if (0 < gemextn.fail_count || nfiles != gemextn.count) {
                printlog ("ERROR - NSOFFSET: Bad output files.", \
                    l_logfile, verbose+) 
                status = 2
                goto clean
            }
        }
    } else {
        touch (tmpout)
    }

    joinlines (tmproot // "," // tmpout, "", output=tmplist, delim=" ", \
        missing="Missing", maxchar=1000, shortest-, verbose-)

    first = yes
    scanin1 = tmplist

    while (fscan (scanin1, imgin, imgout) != EOF) {

        imgmodified = no
        if (debug) print ("testing " // imgin)

        delete (tmpin, verify-, >& "dev$null")
        gemextn (imgin, proc="expand", check="exists", index="", \
            extname=l_sci_ext, extver="1-", ikparams="", omit="", \
            replace="", outfile=tmpin, logfile="dev$null", glogpars="", 
            verbose-)
        if (0 == gemextn.fail_count && 0 < gemextn.count) {
            scanin2 = tmpin
            junk = fscan (scanin2, imgsci)
            havemdf = no

            if (debug) print ("have science data: " // imgsci)

        } else {
            delete (tmpin, verify-, >& "dev$null")
            gemextn (imgin, proc="expand", check="table", index="0-", \
                extname="", extver="", ikparams="", omit="", \
                replace="%\[%.fits[%", outfile=tmpin, logfile="dev$null",
                glogpars="", verbose-)
            # ignore errors since many extns not tables
            if (1 == gemextn.count) {
                scanin2 = tmpin
                junk = fscan (scanin2, mdfin)
                havemdf = yes

                if (debug) print ("have mdf: " // mdfin)

            } else {
                printlog ("ERROR - NSOFFSET: " // imgin // " missing \
                    data/MDF.", l_logfile, verbose+)
                goto clean
            }
        }

        if (l_fl_measure && (first || l_fl_all)) {

            if (debug) print ("need to measure offset")

            # first, generate clean image / template

            if (havemdf) {

                if (debug) print ("construct image from mdf")

                if (isindef (scale)) {
                    if (debug) print ("finding pixscale")
                    hselect (refphu, l_key_pixscale, yes) | scan (scale)
                    if (isindef (scale))
                        scale = l_pixscale
                    if (isindef (scale)) {
                        printlog ("ERROR - NSOFFSET: No header for pixel \
                            scale (" // l_key_pixscale // ")", \
                            l_logfile, verbose+) 
                        printlog ("                  in " // refphu, \
                            l_logfile, verbose+) 
                        goto clean
                    }
                    if (debug) print ("pixscale: " // scale)
                }

                if (debug) print ("base image of zeroes")
                imdelete (tmpmdfimg, verify-, >& "dev$null")
                imarith (template, "*", "0", tmpmdfimg)
                imdelete (template, verify-, >& "dev$null")
                gemhedit (tmpmdfimg, l_key_pixscale, scale, "", delete-)
                imgmodified = yes
                if (debug) hselect (tmpmdfimg, l_key_pixscale, yes)

                if (debug) print ("getting table size")
                tinfo (mdfin, ttout-, >& "dev$null")
                if (tinfo.tbltype != "fits") {
                    printlog ("ERROR - NSOFFSET: Bad MDF in " // mdfin \
                        // ".", l_logfile, verbose+)
                    goto clean
                } else {
                    nsecns = tinfo.nrows
                }
                if (debug) print (nsecns // " sections")

                for (isecn = 1; isecn <= nsecns; isecn += 1) {

                    nsmdfhelper (mdfin, isecn, tmpmdfimg, \
                        area="spectrum", pixscale=scale, \
                        dispaxis=(3-spatialaxis), logfile=l_logfile, \
                        logname="NSOFFSET", verbose=l_verbose)
                    if (no == (0 == nsmdfhelper.status)) {
                        printlog ("ERROR - NSOFSFET: Could not read \
                            MDF.", l_logfile, verbose+)
                        goto clean
                    }

                    corner = nsmdfhelper.corner
                    x1 = nsmdfhelper.ixlo
                    y1 = nsmdfhelper.iylo
                    x2 = nsmdfhelper.ixhi
                    y2 = nsmdfhelper.iyhi

                    if (debug) print ("["//x1//":"//x2//","//y1//":"//y2//"]")

                    imdelete (tmparith, verify-, >& "dev$null")
                    if (0 != corner) {
                        imexpr ("I>="//x1//"+J*"//corner/(y2-y1)\
                            //"&&I<="//x2//"-("//y2//"-J)*"//corner\
                            //"/"//(y2-y1)\
                            //"&&J>="//y1//"&&J<="//y2//"?1:a", \
                            tmparith, tmpmdfimg, 
                            dims="auto", intype="auto", outtype="real", \
                            refim="auto", bwidth=0, btype="nearest", \
                            bpixval=0, rangecheck+, verbose=debug)
                    } else {
                        imexpr ("I>="//x1//"&&I<="//x2//"&&J>="//y1\
                            //"&&J<="//y2//"?1:a", tmparith, tmpmdfimg, \
                            dims="auto", intype="auto", outtype="real", \
                            refim="auto", bwidth=0, btype="nearest", \
                            bpixval=0, rangecheck+, verbose=debug)
                    }
                    imdelete (tmpmdfimg, verify-, >& "dev$null")
                    imrename (tmparith, tmpmdfimg, verbose-)
                }
                if (debug) print ("mdf generated: " // tmpmdfimg)

                if (debug) display (tmpmdfimg, 1)
                if (debug) print ("smoothing mdf image")
                imdelete (tmparith, verify-, >& "dev$null")
                gauss (tmpmdfimg, tmparith, blur, ratio=1, theta=0, \
                    nsigma=4, bilinear+, boundary="constant", constant=0)
                imdelete (tmpmdfimg, verify-, >& "dev$null")
                target = tmparith

            } else {

                if (debug) print ("zeroing bad pixels")
                delete (tmpdq, verify-, >& "dev$null")
                gemextn (imgin, proc="expand", check="exists", index="", \
                    extname=l_dq_ext, extver="1-", ikparams="", omit="", \
                    replace="", outfile=tmpdq, logfile="dev$null", 
                    glogpars="", verbose-)
                if (0 != gemextn.fail_count || 0 == gemextn.count) {
                    target = imgsci
                } else {
                    scanin2 = tmpdq
                    junk = fscan (scanin2, imgdq)
                    imdelete (tmpblank, verify-, >& "dev$null")
                    imexpr ("a>0?0:b", tmpblank, imgdq, imgsci, \
                        dims="auto", intype="auto", outtype="auto", \
                        refim="auto", bwidth=0, btype="nearest", \
                        bpixval=0, rangecheck+, verbose=debug)
                    target = tmpblank
                }

                if (debug) print ("clipping at zero")
                imdelete (tmparith, verify-, >& "dev$null")
                imexpr ("max(0,a)", tmparith, target, \
                    dims="auto", intype="auto", outtype="auto", \
                    refim="auto", bwidth=0, btype="nearest", \
                    bpixval=0, rangecheck+, verbose=debug)
                imdelete (tmpblank // "," // target, verify-, >& "dev$null")
                tmpblank = tmparith; target = tmpblank
            }

            # ifu and longslit data can be handled quickly in 1d, but
            # xd data require full 2d cross-correlation

            if (l_fl_project) {

                # next, project along spectral direction, since only want
                # spatial shift

                imdelete (tmpprjtag, verify-, >& "dev$null")
                improject (target, tmpprjtag, projaxis=(3-spatialaxis), \
                    average+, highcut=0, lowcut=0, pixtype="real", \
                    verbose-)
                imdelete (target, verify-, >& "dev$null")
                target = tmpprjtag
            }

            # taper and cross-correlate

            if (debug) print ("tapering target")
            imdelete (tmptaptag, verify-, >& "dev$null")
            taper (target, tmptaptag, width="10 %", subtract="none", \
                function="cosbell", verbose-, >& "dev$null")
            if (debug) {
                if (l_fl_project)
                    splot (tmptaptag)
                else
                    display (tmptaptag, 3)
            }

            if (debug) {
                print ("cross-correlating")
                imhead (tmptapref)
                imhead (tmptaptag)
            }
            imdelete (tmpconv, verify-, >& "dev$null")
            printlog ("NSOFFSET: Calculating cross correlation \
                (please wait).", l_logfile, l_verbose)
            crosscor (tmptapref, tmptaptag, tmpconv, inreal1+, inimag1-, \
                inreal2+, inimag2-, outreal+, outimag-, coord_shift-, \
                center+, chop+, pad+, inmemory+, len_blk=256, \
                verbose=debug)
            if (debug) {
                if (l_fl_project)
                    splot (tmpconv)
                else
                    display (tmpconv, 4)
            }

            # find peak (strongest pixel)

            if (l_fl_project) {
                spec = tmpconv
            } else {
                hselect (tmpconv, "CRPIX" // (3 - spatialaxis), yes) \
                    | scan (midpix)
                if (debug) print (midpix // ", " // spatialaxis)
                if (1 == spatialaxis) {
                    spec = tmpconv // "[*," // midpix // "]"
                } else {
                    spec = tmpconv // "[" // midpix // ",*]"
                }
            }
            gemhedit (tmpconv, "DISPAXIS", "2", "Stop splot asking", delete-)
            imgmodified = yes

            minmax (spec, force+, update-, verbose-)
            printf ("%s\n", minmax.maxpix) | scanf ("[%d]", junk)
            l_shift = junk
            if (debug) print ("shift " // l_shift)

            # see help crosscor for explanation
            if (l_fl_project) {
                hselect (tmpconv, "CRPIX1", yes) | scan (crpix)
            } else {
                hselect (tmpconv, "CRPIX" // spatialaxis, yes) \
                    | scan (crpix)
            }
            if (1 == spatialaxis)
                l_shift = l_shift - crpix
            else
                l_shift = l_shift - crpix

            printlog ("NSOFFSET: Integer offset for " // imgin // " is " \
                // nint (l_shift) // ".", l_logfile, l_verbose)

            # refine peak to give sub-pixel value

            if (debug) print ("refining peak position")
            delete (tmpcursor, verify-, >& "dev$null")

            printf ("%d 0 0 k\n", (l_shift - fitwidth), >> tmpcursor)
            printf ("%d 0 0 k\n", (l_shift + fitwidth), >> tmpcursor)
            if (debug) type (tmpcursor)

            # if interactive, do twice (identical apart from 
            # redirection of output and source of commands)

            if (l_fl_inter)
                imax = 2
            else
                imax = 1

            for (i = 1; i <= imax; i = i + 1) {

                if (1 == i) {

                    # If the WCS keywords are set for system=image (as is the
                    # case for NIRI data) the calls to splot below won't work 
                    # correctly. Delete the appropriate WCS keywords from the
                    # header of the image to be plotted
                    gemhedit (spec, "CD1_1", "", "", delete+)
                    gemhedit (spec, "CD1_2", "", "", delete+)
                    gemhedit (spec, "CD2_2", "", "", delete+)
                    gemhedit (spec, "CD2_1", "", "", delete+)
                    gemhedit (spec, "CTYPE1", "", "", delete+)
                    gemhedit (spec, "CTYPE2", "", "", delete+)

                    splot (images=spec, line=1, band=1, \
                        star_name="cross-correlation", mag=INDEF, \
                        teff=INDEF, next_image="", \
                        new_image="", \
                        overwrite=no, spec2="", constant=0, \
                        wavelength=INDEF, linelist="", wstart=INDEF, \
                        wend=INDEF, dw=INDEF, boxsize=INDEF, units="", \
                        options="auto wreset", xmin=(l_shift-dispwidth), \
                        xmax=(l_shift+dispwidth), \
                        ymin=INDEF, ymax=INDEF, save_file=tmplog1, \
                        graphics="stdgraph", cursor=tmpcursor, \
                        nerrsample=0, \
                        sigma0=INDEF, invgain=INDEF, function="spline3", \
                        order=1, low_reject=2, high_reject=4, \
                        niterate=10, grow=1, markrej+, \
                        caldir=")_.caldir", fnuzero=3.68e-20, \
                        >G "dev$null", > "dev$null")

                    if (debug) type (tmplog1)

                    scanin2 = tmplog1
                    word = ""
                    while (fscan (scanin2, word) != EOF) {
                        if (debug) print ("read: " // word)
                        if ("center" == word) {
                            if (fscan (scanin2, l_shift) == EOF) {
                                printlog ("ERROR - NSOFFSET: Bad splot log.", \
                                    l_logfile, verbose+)
                                goto clean
                            } else {
                                if (debug)
                                    print ("read new centre: " // l_shift)
                            }
                        }
                       word = ""
                    }

                } else {

                    printlog ("NSOFFSET: Use k key to fit profile.", \
                        l_logfile, l_verbose)

                    splot (images=spec, line=1, band=1, \
                        star_name="cross-correlation", mag=INDEF, \
                        teff=INDEF, next_image="", \
                        new_image="", \
                        overwrite=no, spec2="", constant=0, \
                        wavelength=INDEF, linelist="", wstart=INDEF, \
                        wend=INDEF, dw=INDEF, boxsize=INDEF, units="", \
                        options="auto wreset", xmin=(l_shift-dispwidth), \
                        xmax=(l_shift+dispwidth), \
                        ymin=INDEF, ymax=INDEF, save_file=tmplog2, \
                        graphics="stdgraph", cursor="", \
                        nerrsample=0, \
                        sigma0=INDEF, invgain=INDEF, function="spline3", \
                        order=1, low_reject=2, high_reject=4, \
                        niterate=10, grow=1, markrej+, \
                        caldir=")_.caldir", fnuzero=3.68e-20)

                    if (debug) type (tmplog2)

                    # Assume that the last measurement is the one that will be
                    # used. Therefore, the shift corresponds to the value in
                    # the first column on the last row
                    tail (tmplog2, nlines=1) | fields ("STDIN", "1", \
                        lines="1") | scan (l_shift)

                    if (debug) print ("read new centre: " // l_shift)
                }
            }

            # restrict to 0.1 pix resoln
            l_shift = 0.1 * nint (10 * l_shift)
            printlog ("NSOFFSET: Final offset for " // imgin // " is " \
                // l_shift // ".", l_logfile, l_verbose)
        }

        if (l_fl_apply) {

            if (debug) print ("copying input to output")
            if (debug) print (imgin // " -> " // imgout)
            copy (imgin // ".fits", imgout // ".fits", verbose=debug)

            if (havemdf) {

                if (debug) print ("finding table in output")
                delete (tmptab, verify-, >& "dev$null")
                gemextn (imgout, proc="expand", check="table", \
                    index="0-", extname="", extver="", ikparams="", \
                    omit="", replace="%\[%.fits[%", outfile=tmptab, \
                    logfile="dev$null", glogpars="", verbose-)
                scanin2 = tmptab
                junk = fscan (scanin2, mdfout)

                if (debug) print ("updating table")
                if (1 == spatialaxis) {
                    tcalc (mdfout, outcol="x_ccd", \
                        equals="x_ccd-" // l_shift, datatype="real")
                    printlog ("NSOFFSET: Shifting " // mdfout \
                        // " in X by " // l_shift // " pixels.", \
                        l_logfile, l_verbose) 
                } else {
                    tcalc (mdfout, outcol="y_ccd", \
                        equals="y_ccd-" // l_shift, datatype="real")
                    printlog ("NSOFFSET: Shifting " // mdfout \
                        // " in Y by " // l_shift // " pixels.", \
                        l_logfile, l_verbose)
                }
            } else {

                if (debug) print ("copying reference wcs")
                gemwcscopy (imgout // "[0]", refphu, verbose-, \
                    logfile=l_logfile)

                if (debug) print ("shifting wcs")
                hselect (imgout // "[0]", "CRPIX" // spatialaxis, yes) \
                    | scan (crpix)
                gemhedit (imgout // "[0]", "CRPIX" // spatialaxis, \
                    (crpix+l_shift), "", delete-)
                imgmodified = yes
                if (1 == spatialaxis) {
                    printlog ("NSOFFSET: Shifting " // imgout \
                        // " WCS in X by " // l_shift // " pixels.", \
                        l_logfile, l_verbose) 
                } else {
                    printlog ("NSOFFSET: Shifting " // imgout \
                        // " WCS in Y by " // l_shift // " pixels.", \
                        l_logfile, l_verbose)
                }
            }
        }

        if (imgmodified == yes) {
            gemdate ()
            gemhedit (imgout//"[0]", "NSOFFSET", gemdate.outdate,
                "UT Time stamp for NSOFFSET", delete-)
            gemhedit (imgout//"[0]", "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI IRAF", delete-)
        }

        first = no
    }

    status = 0

clean:

    shift = l_shift

    if (0 != status && 2 != status && access (tmpout)) {
        delete ("@" // tmpout, verify-, >& "dev$null")
        imdelete ("@" // tmpout, verify-, >& "dev$null")
    }

    delete (tmpref, verify-, >& "dev$null")
    delete (tmprefsci, verify-, >& "dev$null")
    delete (tmproot, verify-, >& "dev$null")
    delete (tmpout, verify-, >& "dev$null")
    delete (tmplist, verify-, >& "dev$null")
    delete (tmpin, verify-, >& "dev$null")
    delete (tmpdq, verify-, >& "dev$null")
    imdelete (tmpblank, verify-, >& "dev$null")
    imdelete (tmpmask, verify-, >& "dev$null")
    imdelete (tmpmdfimg, verify-, >& "dev$null")
    imdelete (tmparith, verify-, >& "dev$null")
    imdelete (tmptapref, verify-, >& "dev$null")
    imdelete (tmptaptag, verify-, >& "dev$null")
    imdelete (tmpconv, verify-, >& "dev$null")
    delete (tmpcursor, verify-, >& "dev$null")
    delete (tmplog1, verify-, >& "dev$null")
    delete (tmplog2, verify-, >& "dev$null")
    delete (tmptab, verify-, >& "dev$null")
    imdelete (tmpprjtag, verify-, >& "dev$null")
    imdelete (tmpprjref, verify-, >& "dev$null")

    scanin1 = ""
    scanin2 = ""

    printlog (" ", l_logfile, l_verbose) 
    if (0 == status) {
        printlog ("NSOFFSET exit status: good.", l_logfile, l_verbose) 
    } else {
        printlog ("NSOFFSET exit status: error.", l_logfile, l_verbose) 
    }
    printlog ("---------------------------------------------------------\
        ----------------------", l_logfile, l_verbose) 

end
