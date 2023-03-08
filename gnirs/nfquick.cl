# Copyright(c) 2004-2011 Association of Universities for Research in Astronomy, Inc.

# acooke 28 Mar 2004
# brodgers 27 Dec 2004: corrected p and q offsets to match reality (with new TCS),
#                       moved y center down 1 slice (to slice 10)


procedure nfquick (image)

char    image       {prompt = "Acquisition image (int range)"}
char    prefix      {"", prompt = "File prefix (eg S20041227S)"}
char    datadir     {"/net/reggie/staging/perm/", prompt = "Data directory"}
char    outimages   {"", prompt = "Output files"}
char    outprefix   {"a", prompt = "Prefix to use if outimages not given"}

char    mdf         {"gnirs$data/gnirs-ifu-short-32-mdf2.fits", prompt = "MDF that describes mapping"}
real    shiftx      {0.0, prompt = "Shift in X"}
real    shifty      {0.0, prompt = "Shift in Y"}
real    xscale      {0.15, prompt = "Size of pixels in X (arcec)"}
real    yscale      {0.15, prompt = "Size of pixels in Y (arcec)"}
int     dispaxis    {2, prompt = "Dispersion axis"}
bool    fl_inter    {yes, prompt = "Get shift interactively"}
bool    fl_compress {yes, prompt = "Compress to single row?"}
bool    fl_display  {yes, prompt = "Display each image as it is generated"}
int     raw_frame   {1, prompt = "Frame to display raw data"}
int     img_frame   {2, prompt = "Frame to display reconstructed image"}
bool    fl_shift    {yes, prompt = "Get shift to centre star in each image"}

char    logfile     {"", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose output?"}

int     status      {0, prompt = "O: Exit status (0 = good)"}
struct  *scanin1    {"", prompt = "Internal use only"}
struct  *scanin2    {"", prompt = "Internal use only"}

begin
    char    l_image = ""
    char    l_prefix = ""
    char    l_datadir = ""
    char    l_outimages = ""
    char    l_outprefix = ""
    char    l_mdf = ""
    real    l_shiftx, l_shifty, l_xscale, l_yscale
    int     l_dispaxis
    char    l_logfile = ""
    bool    l_fl_inter
    bool    l_fl_compress
    bool    l_fl_display
    int     l_raw_frame
    int     l_img_frame
    bool    l_fl_shift
    bool    l_verbose

    char    l_key_pixscale = ""

    int     junk, nfiles, nx, ny, mx, my
    bool    debug, havex, havey, havexy, first, firstimg
    struct  sdate, sline
    char    tmpin, tmpout, tmplist, inimg, outimg, src, dst, tmpexpand
    char    tmproot, tmpwhole, tmpthin, tmp1, tmp2, tmp3, tmp4, tmplog
    char    eximg, tmptable, tmpimg, badhdr, tmpmdf, tmpstack
    real    x1, y1, dx, dy, x2, y2, dy2, xmid, dmin, poff, qoff
    int     ix1, iy1, idx, idy, ix2, iy2, idx2, idy2
    int     xmin, xmax, ymin, ymax

    l_image = image
    junk = fscan (prefix, l_prefix)
    junk = fscan (datadir, l_datadir)
    junk = fscan (outimages, l_outimages)
    junk = fscan (outprefix, l_outprefix)
    junk = fscan (mdf, l_mdf)
    l_shiftx = shiftx
    l_shifty = shifty
    l_xscale = xscale
    l_yscale = yscale
    l_dispaxis = dispaxis
    junk = fscan (logfile, l_logfile)
    l_fl_inter = fl_inter
    l_fl_compress = fl_compress
    l_fl_display = fl_display
    l_raw_frame = raw_frame
    l_img_frame = img_frame
    l_fl_shift = fl_shift
    l_verbose = verbose

    status = 1
    debug = no
    badhdr = ""

    junk = fscan (nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"

    tmpexpand = mktemp ("tmpexpand")
    tmpin = mktemp ("tmpin")
    tmpout = mktemp ("tmpout")
    tmplist = mktemp ("tmplist")
    tmproot = mktemp ("tmproot")
    tmplog = mktemp ("tmplog")
    tmptable = mktemp ("tmptable")
    tmpmdf = mktemp ("tmpmdf")
    tmpimg = mktemp ("tmpimg")


    cache ("gemextn", "keypar", "nsmdfhelper", "gemdate")


    # Start logging

    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile) 
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NFQUICK: Both gnswedit.logfile and \
                gnirs.logfile are empty.", l_logfile, verbose+) 
            printlog ("                   Using default file " \
                // l_logfile // ".", l_logfile, verbose+) 
        }
    }
    date | scan (sdate) 
    printlog ("-------------------------------------------------------\
        -----------------------", l_logfile, verbose = l_verbose) 
    printlog ("NFQUICK -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NFQUICK: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }


    # Check input data

    if (debug) print ("checking input")
    if (no == ("" == l_prefix)) {
        gemlist (root=l_prefix, range=l_image, >& tmpexpand)
        src = "@" // tmpexpand
    } else {
        src = l_image
    }

    gemextn (src, proc="none", check="", index="", \
        extname="", extver="", ikparams="", omit="path", replace="", \
        outfile=tmproot, logfile=l_logfile, glogpars="", verbose-)
    if (0 < gemextn.fail_count || 0 == gemextn.count) {
        printlog ("ERROR - NFQUICK: Missing input data.", \
            l_logfile, verbose+) 
        goto clean
    }
    nfiles = gemextn.count
    if (debug) print (nfiles)
    while (strlen (l_datadir) > 0 && \
        "/" == substr (l_datadir, strlen (l_datadir), \
        strlen (l_datadir))) {
        
        l_datadir = substr (l_datadir, 1, strlen (l_datadir) - 1)
    }
    if (no == ("" == l_datadir)) l_datadir = l_datadir // "/"
    if (debug) print (l_datadir)
    if (debug) type (tmproot)
    if (debug) {
        gemextn (l_datadir // "@" // tmproot, proc="expand", \
            check="", index="-", extname="", extver="", ikparams="", \
            omit="", replace="", outfile="STDOUT", logfile=l_logfile, \
            glogpars="", verbose+)
    }
    gemextn (l_datadir // "@" // tmproot, proc="expand", \
        check="exist,mef", index="1", \
        extname="", extver="", ikparams="", omit="", replace="", \
        outfile=tmpin, logfile=l_logfile, glogpars="", verbose-)
    if (0 < gemextn.fail_count || nfiles != gemextn.count) {
        printlog ("ERROR - NFQUICK: Bad input data.", l_logfile, verbose+) 
        goto clean
    }
    if (debug) type (tmpin)


    # Check/generate output names

    if (debug) print ("checking output")
    gemextn (l_outimages, proc="none", check="absent", index="", \
        extname="", extver="", ikparams="", omit="", replace="", \
        outfile=tmpout, logfile="", glogpars="", verbose-)
    if (0 < gemextn.fail_count) {
        printlog ("ERROR - NFQUICK: Output data already exist.", \
            l_logfile, verbose+) 
        goto clean
    }
    if (0 != gemextn.count) {
        if (nfiles != gemextn.count) {
            printlog ("ERROR - NFQUICK: Incorrect number of output \
                files.", l_logfile, verbose+) 
            goto clean
        }
    } else {
        delete (tmpout, verify-, >& "dev$null")
        gemextn (l_outprefix // "@" // tmproot,
            proc="none", check="absent", index="", extname="", \
            extver="", ikparams="", omit="kernel", replace="", \
            outfile=tmpout, logfile="", glogpars="", verbose-)
        if (0 < gemextn.fail_count || nfiles != gemextn.count) {
            printlog ("ERROR - NFQUICK: Bad output files.", \
                l_logfile, verbose+) 
            goto clean
        }
    }


    # Combine in + out into single list

    joinlines (tmpin, tmpout, output=tmplist, delim=" ", missing="Missing", \
        maxchar=161, shortest=yes, verbose-)

    # Get output image size from first input, extn 1

    eximg = ""
    head (tmpin) | scan (eximg)
    if (debug) print (eximg)

    nx = 0; ny = 0
    hselect (eximg, "i_naxis1", yes) | scan (nx)
    hselect (eximg, "i_naxis2", yes) | scan (ny)
    if (0 == nx || 0 == ny) {
        printlog ("ERROR - NFQUICK: Cannot get size of first image \
            extension " // eximg // ".", l_logfile, verbose+) 
        goto clean
    }


    # Get centring interactively, if required

    if (l_fl_inter) {

        havex = no; havey = no
        while (no == havex || no == havey) {

            if (no == havey)
                printlog ("NFQUICK: press k to identify the slit \
                    in the y direction", l_logfile, verbose+) 
            if (no == havex)
                printlog ("NFQUICK: press j to identify a central gap \
                    in the x direction, and then q to quit", \
                    l_logfile, verbose+)

            scanin1 = tmplist
            junk = fscan (scanin1, inimg)
            display (inimg, l_raw_frame)
            imexam (inimg, frame=l_raw_frame, image="", logfile=tmplog,
                keeplog+)

            if (debug) type (tmplog)

            scanin2 = tmplog
            while (fscan (scanin2, sline) != EOF) {
                if (1 == stridx (sline, "Column")) {
                    havey = yes
                    junk = fscan (sline, tmp1, tmp2, tmp3, tmp4)
                    if (strlen (tmp3) <= 8) sline = tmp4
                    else sline = substr (tmp3, 8, strlen (tmp3))
                    junk = fscan (sline, l_shifty)
                }
                if (1 == stridx (sline, "Line")) {
                    havex = yes
                    junk = fscan (sline, tmp1, tmp2, tmp3, tmp4)
                    if (strlen (tmp3) <= 8) sline = tmp4
                    else sline = substr (tmp3, 8, strlen (tmp3))
                    junk = fscan (sline, l_shiftx)
                }
            }

            if (debug) print ("trying to delete log")
            delete (tmplog, verify-, >& "dev$null")
        }

        if (2 == l_dispaxis) l_shifty = ny / 2 - l_shifty
        else                 l_shiftx = nx / 2 - l_shiftx

    }


    # Generate the map table (without shift)

    imcopy (eximg, tmpimg)
    gemhedit (tmpimg, l_key_pixscale, l_xscale, "", delete-)
    nsmdfhelper (l_mdf, 1, tmpimg, area="ifu-map", mapfile=tmptable,
        pixscale=INDEF, dispaxis=2, logname="NFQUICK")
    if (no == (0 == nsmdfhelper.status)) {
        printlog ("ERROR - NFQUICK: Cannot get map table.", \
            l_logfile, verbose+) 
        goto clean
    }


    # Generate the images

    scanin1 = tmplist
    firstimg = yes

    while (fscan (scanin1, inimg, outimg) != EOF) {

        if (debug) print (inimg)
        if (debug)
            print ("mkimage " // outimg // "[" // nx // "," // ny // "]")

        mkimage (outimg, option = "make", value = 0.0, ndim = 2, \
            dims = nx // " " // ny, pixtype = "real", slope = 0.0, \
            sigma = 0.0, seed = 0)


        # Get nearest appropriate aperture shift

        if (l_fl_inter && firstimg) {

            if (2 == l_dispaxis) {
                tdump (l_mdf, cdfile="", pfile="", datafile=tmpmdf, \
                    columns="x_ccd", rows="-", pwidth=-1)
            } else {
                tdump (l_mdf, cdfile="", pfile="", datafile=tmpmdf, \
                    columns="y_ccd", rows="-", pwidth=-1)
            }

            scanin2 = tmpmdf
            havex = no
            first = yes
            if (2 == l_dispaxis) dmin = l_shiftx
            else                 dmin = l_shifty

            while (! havex && fscan (scanin2, x1) != EOF) {

                if (first) {
                    first = no
                } else {
                    xmid = 0.5 * (x1 + x2)
                    if (2 == l_dispaxis) dx = xmid - l_shiftx
                    else                 dx = xmid - l_shifty
                    if (abs (dx) < abs (dmin)) dmin = dx
                }
                x2 = x1
            }

            if (2 == l_dispaxis) l_shiftx = dmin
            else                 l_shifty = dmin

            printlog ("NFQUICK: Offsets: x: " // l_shiftx // "; y: " \
                // l_shifty, l_logfile, l_verbose) 
            firstimg = no
        }


        scanin2 = tmptable
        if (debug) type (tmptable)

        first = yes
        while (fscan (scanin2, ix1, iy1, idx, idy, ix2, iy2) != EOF) {

            src = "[" // int (ix1-l_shiftx) // ":" \
                // int (ix1+idx-l_shiftx) \
                // "," \
                // int (iy1-l_shifty) // ":" \
                // int (iy1+idy-l_shifty) // "]"
            src = inimg // src
            if (debug) print (src)

            if (l_fl_compress) {
                if (2 == l_dispaxis) {idy2 = 0; idx2 = idx}
                else                 {idy2 = idy; idx2 = 0}
                tmpwhole = mktemp ("tmpwhole")
                tmpthin = mktemp ("tmpthin")
                tmpstack = mktemp ("tmpstack")
                imcopy (src, tmpwhole, >& "dev$null")
                if (debug) imhead (tmpwhole)
                improject (tmpwhole, tmpthin, projaxis=l_dispaxis, \
                    average+, highcut=0.0, lowcut=0.0, pixtype="real", 
                    verbose-, >& "dev$null")
                if (1 == l_dispaxis) {
                    imstack (tmpthin, tmpstack, title="", pixtype="*")
                    imdelete (tmpthin, verify-, >& "dev$null")
                    tmpthin = mktemp ("tmpthin")
                    imtranspose (tmpstack, tmpthin)
                    if (debug) imhead (tmpthin)
                    src = tmpthin // "[1:1,*]"
                } else {
                    src = tmpthin
                }
            } else {
                idy2 = idy
                idx2 = idx
            }

            dst = "[" // ix2 // ":" // ix2 + idx2 // "," \
                // iy2 // ":" // iy2 + idy2 // "]"

            if (first) {
                xmin = ix2
                xmax = ix2 + idx2
                ymin = iy2
                ymax = iy2 + idy2
                first = no
            } else {
                xmin = min (xmin, ix2)
                xmax = max (xmax, ix2 + idx2)
                ymin = min (ymin, iy2)
                ymax = max (ymax, iy2 + idy2)
            }

            if (debug) print (src // " -> " // outimg // dst)

            if (debug) imstat (src)
            imcopy (src, outimg // dst, >& "dev$null")
            if (debug) imstat (outimg // dst)

            imdelete (tmpwhole, verify-, >& "dev$null")
            imdelete (tmpthin, verify-, >& "dev$null")
            imdelete (tmpstack, verify-, >& "dev$null")

        }


        # Shrink output image

        imcopy (outimg // "[" // xmin // ":" // xmax // "," \
            // ymin // ":" // ymax // "]", outimg, >& "dev$null")

        # Give shift info for star centre

        if (l_fl_shift) {

            havexy = no
            while (no == havexy) {

                printlog ("NFQUICK: press a or x to identify star \
                    and then q to quit", l_logfile, verbose+) 

                display (outimg, l_img_frame)
                imexam (outimg, frame=l_img_frame, image="", \
                    logfile=tmplog, keeplog+)

                if (debug) type (tmplog)

                scanin2 = tmplog
                x1 = -1.0
                y1 = -1.0
                while (! havexy && fscan (scanin2, x1, y1) != EOF) {
                    if (x1 >= 0.0 && y1 >= 0.0) {
                        havexy = yes
                    } else {
                        x1 = -1.0
                        y1 = -1.0
                    }
                }

                if (debug) print ("trying to delete log")
                delete (tmplog, verify-, >& "dev$null")

            }

            # this is all horribly GNIRS specific

            mx = xmax - xmin + 1
            #reduce mid-y by 1 to move farther from bad slice 
            #(bad slice=13, mid-slice=11 but use slice 10)
            my = ymax - ymin 
            x2 = 0.5 * real (mx)
            y2 = 0.5 * real (my)
            dx = x1 - x2
            dy = y1 - y2
            #instrument offsets: 
            poff = dy * l_yscale
            #q offsets increase to the left (-x)
            qoff = dx * l_xscale * (-1.0)
            printlog (" ", l_logfile, verbose+)
            print ("****************************************\
                ***************************************")
            printf ("    GNIRS IFU instrument offsets (arcsec)  p=%-6.3f \
                q=%-6.3f\n", poff, qoff) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            print ("****************************************\
                ***************************************")
            printlog (" ", l_logfile, verbose+)
            printlog ("    Raw data source  " //inimg, l_logfile, verbose+)
            printf ("%35s  (%-6.3f, %-6.3f)\n", \
                "Image size (pixels)", (xmax-xmin+1.0), \
                (ymax-ymin+1.0)) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            printf ("%35s  (%-6.3f, %-6.3f)\n", \
                "Standard location (pixels)", x2, y2) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            printf ("%35s  (%-6.3f, %-6.3f)\n", \
                "Measured peak (pixels)", x1, y1) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            printf ("%35s  (%-6.3f, %-6.3f)\n", \
                "Offset (pixels)", dx, dy) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            printf ("%35s  (%-6.3f, %-6.3f)\n", \
                "Offset (arcsec)", dx * l_xscale, \
                dy * l_yscale) | scan (sline)
            printlog (sline, l_logfile, verbose+)

        } else if (l_fl_display) {
            display (outimg, l_img_frame)
        }

        gemdate ()
        gemhedit (outimg, "NFQUICK", gemdate.outdate,
            "UT Time stamp for NFQUICK", delete-)
        gemhedit (outimg, "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI IRAF", delete-)
    }

    status = 0 # Success

clean:

    scanin1 = ""
    scanin2 = ""

    delete (tmpexpand, verify-, >& "dev$null")
    delete (tmpin, verify-, >& "dev$null")
    delete (tmpout, verify-, >& "dev$null")
    delete (tmplist, verify-, >& "dev$null")
    delete (tmptable, verify-, >& "dev$null")
    delete (tmpmdf, verify-, >& "dev$null")
    delete (tmproot, verify-, >& "dev$null")
    imdelete (tmpimg, verify-, >& "dev$null")

    if (0 == status) {
        printlog (" ", l_logfile, l_verbose) 
        printlog ("NFQUICK exit status:  good.", l_logfile, l_verbose) 
    }

end
