# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

# NFACQUIRE - Collapse NIFS flip-mirror frames into 2D images of the field
#             for acquisition purposes.
#
# Original author: Tracy Beck
#
# Warning: very NIFS-specific
    
procedure nfacquire (image)

int     image       {prompt = "Acquisition image (int range)"}
char    rootname    {"default",prompt="Rootname of image"}
char    datadir     {"adata$", prompt = "Data directory"}
int     sky         {0,  prompt ="Acquisition sky image (int range)"}
char    outimage    {"", prompt = "Output files"}
char    outprefix   {"a", prompt = "Prefix to use if outimage not given"}

char    mdf         {"nifs$data/nifs-mdf.fits", prompt = "MDF that describes mapping"}
real    shiftx      {0.0, prompt = "Shift in X"}
real    shifty      {0.0, prompt = "Shift in Y"}
real    xpos        {0.0, prompt = "Centering position in X"}
real    ypos        {0.0, prompt = "Centering position in Y"}
real    yscale      {0.043, prompt = "Size of pixels in Y (arcec)"}
int     dispaxis    {1, prompt = "Dispersion axis"}
bool    fl_inter    {yes, prompt = "Get shift interactively"}
bool    fl_display  {yes, prompt = "Display each image as it is generated"}
int     raw_frame   {1, prompt = "Frame to display raw data"}
int     img_frame   {2, prompt = "Frame to display reconstructed image"}
bool    fl_shift    {yes, prompt = "Get shift to centre star in each image"}
char    logfile     {"", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose output?"}
int     status      {0, prompt = "Exit status (0 = good)"}
struct  *scanin1    {"", prompt = "Internal use only"}

begin

    int     l_image
    char    l_rootname = ""
    char    l_datadir = ""
    int     l_sky
    char    l_outimage = ""
    char    l_outprefix = ""
    char    l_mdf = ""
    real    l_shiftx, l_shifty, l_yscale, l_xpos, l_ypos
    int     l_dispaxis
    bool    l_fl_inter
    bool    l_fl_display
    int     l_raw_frame
    int     l_img_frame
    bool    l_fl_shift
    char    l_logfile = ""
    bool    l_verbose
    char    l_key_pixscale

    char    lastchar
    char    rootn, inimg, skyimg, outimg
    char    imstring, skystring
    char    dia, mes, zone, ano
    char    ddir, prefix, tmpout2, tmpout4, tmp_sky, tmpfile
    int     junk, nx, ny, mx, my
    bool    debug, havex, havey, havexy, first, firstimg, offonly
    struct  sline
    char    tmpin, tmpout, tmpout1, tmpout3, src, dst
    char    tmpwhole, tmpthin, tmp1, tmp2, tmp3, tmp4, tmplog
    char    eximg, tmpimg, badhdr, tmpmdf, tmpstack, phu
    char    instrument, tmp_inimg, tmp_outimg, guidestate
    real    x1, y1, dx, dy, x2, y2, xmid, dmin, poff, qoff, guidevalue, r
    int     ix1, iy1, idx, idy, ix2, iy2, idx2, idy2, row
    int     i, j, xtest1, ytest1
    int     xmin, xmax, ymin, ymax, mdfrow
    real    x_ccd, y_ccd, slitsize_x, slitsize_y, specord, xdcorner, ra, dec
    real    slicewidth, sliceheight, xtest, ytest
    real    rjunk
    char    cursinput

    cache ("gemextn", "gemdate")

    status = 1
    offonly = no
    debug = no
    badhdr = ""
    row = 28
    
    l_image     = image
    
    junk = fscan (rootname, l_rootname)
    junk = fscan (datadir, l_datadir)
    l_sky       = sky
    junk = fscan (outimage, l_outimage)
    junk = fscan (outprefix, l_outprefix)
    junk = fscan (mdf, l_mdf)
    l_shiftx    = shiftx
    l_shifty    = shifty
    l_yscale    = yscale
    l_xpos      = xpos
    l_ypos      = ypos
    l_dispaxis  = dispaxis
    l_fl_inter  = fl_inter
    l_fl_display = fl_display
    l_raw_frame = raw_frame
    l_img_frame = img_frame
    l_fl_shift  = fl_shift
    junk = fscan (logfile, l_logfile)
    l_verbose   = verbose

    nsheaders("nifs")

    junk = fscan ( nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"

    xtest1=465
    ytest1=80

    # Test the logfile
    if (l_logfile == "") {
        l_logfile = nifs.logfile
        if (l_logfile == "") {
            l_logfile = "nifs.log"
            printlog ("WARNING - NFIMAGE: Both nfacquire.logfile and \
                nifs.logfile are empty.", l_logfile, l_verbose)
            printlog ("                   Using default file nifs.log",
                l_logfile, l_verbose)
        }
    }

    tmpfile = mktemp ("tmpfile")
    tmpout2 = mktemp ("tmpout2")
    tmp_sky = mktemp ("tmp_sky")
    tmpout3 = mktemp ("tmpout3")
    tmpout4 = mktemp ("tmpout4")
    tmpin = mktemp ("tmpin")
    tmpout = mktemp ("tmpout")
    tmpout1 = mktemp ("tmpout1")
    tmplog = mktemp ("tmplog")
    tmpmdf = mktemp ("tmpmdf")
    tmpimg = mktemp ("tmpimg")
    tmp_inimg = mktemp ("tmp_inimg")
    tmp_outimg = mktemp ("tmp_outimg")

    printlog ("------------------------------------------------------------------------", l_logfile, verbose+) 
    printlog ("TASK - NFACQUIRE:  Derive offsets for NIFS flip-mirror images", l_logfile, verbose+) 
    printlog ("------------------------------------------------------------------------", l_logfile, verbose+) 

    # Ensure a trailing '/' on the directory name
    lastchar = substr (l_datadir, strlen(l_datadir), strlen(l_datadir))
    if ((l_datadir != "") && (lastchar != "$") && (lastchar != "/"))
        l_datadir = l_datadir//"/"
    
    ddir = l_datadir

    # Build image names based on date or rootname
    print (l_rootname) | lcase | scan(rootn)
    printf ("%04d\n", l_image) | scan(imstring)
    printf ("%04d\n", l_sky) | scan(skystring)

    if (rootn == "default") {
        date ("-u +'%m %d %Y %Z'") | scan (mes, dia, ano, zone)
        if ((zone != "UTC") && (zone != "GMT")) {
            print ("ERROR: Internal Error")
            print ("ERROR: Date not given as universal time ("//zone//")")
            goto clean
        }
        
        # Build image names
        prefix = "N"//ano//mes//dia//"S"
        inimg = ddir//"N"//ano//mes//dia//"S"//imstring//".fits"
    } else {
        print(rootn) | ucase | scan(rootn)
        prefix = rootn
        inimg = ddir//rootn//""//imstring//".fits"
    }

    if (l_outimage == "") {
        if (l_outprefix == "") {
            if (rootn == "default")
                outimg = "aN"//ano//mes//dia//"S"//imstring//".fits"
            else
                outimg = "a"//rootn//""//imstring//".fits"
        } else {
            if (rootn == "default")
                outimg = l_outprefix//"N"//ano//mes//dia//"S"//imstring//".fits"    
            else
                outimg = l_outprefix//rootn//""//imstring//".fits"
        }
    } else
        outimg = l_outimage

    skyimg = ""
    if (l_sky != 0) {
        if (rootn == "default") {
            skyimg = ddir//"N"//ano//mes//dia//"S"//skystring//".fits"
        } else {
            skyimg = ddir//rootn//""//skystring//".fits"
        }
    }

    #####
    # Check input data

    src = inimg

    gemextn (src, proc="none", check="", index="", extname="", extver="", 
        ikparams="", omit="", replace="", outfile=tmpin, 
        logfile=l_logfile, glogpars="", verbose-)
    if ((0 < gemextn.fail_count) || (0 == gemextn.count)) {
        printlog ("ERROR - NFACQUIRE: Missing input data.",
            l_logfile, verbose+) 
        goto clean
    }

    src = skyimg

    gemextn (src, proc="none", check="exists", index="", extname="", extver="",
        ikparams="", omit="", replace="", outfile=tmpin, logfile=l_logfile,
        glogpars="", verbose-)
    if ((0 < gemextn.fail_count) || (0 == gemextn.count)) {
        printlog ("NFIMAGE: No or Missing input sky data. Not doing sky \
            subtraction", l_logfile, verbose+) 
        skyimg = ""
    }

    phu = inimg // "[0]"

    #####
    # Check output data

    gemextn (outimg, proc="none", check="absent", index="", extname="",
        extver="", ikparams="", omit="", replace="", outfile=tmpout, 
        logfile="", glogpars="", verbose-)
    if (0 < gemextn.fail_count) {
        printlog ("WARNING - NFACQUIRE: Output data already exist.",
            l_logfile, verbose+) 
        printlog ("WARNING - NFACQUIRE: Displaying frame to generate new \
            offsets.", l_logfile, verbose+) 
        l_fl_inter = no
        offonly = yes
        printlog ("----------------------------------------------------------\
            --------------", l_logfile, verbose+) 
    }

   #####
   # Determine guide state for offset value flip.

            hselect (images=phu, fields="AOFOLD", expr=yes) | \
                scan (guidestate)

    if (guidestate == 'IN') {
       guidevalue=-1.0
        printlog ("NACQUIRE:  NIFS Guiding with Altair AO", l_logfile, verbose+) 
    } else {
       guidevalue=1.0
        printlog ("NFACQUIRE:  NIFS Guiding without AO", l_logfile, verbose+) 
    }

    if (skyimg == "") {

        printlog ("NFACQUIRE: ERROR!  Flip mirror acquisition with NO sky defined!", l_logfile, verbose+) 
        printlog ("NFACQUIRE: ERROR! Acquisition filter is too broad to not subtract sky", l_logfile, verbose+) 
        printlog ("NFACQUIRE: Exiting with ERROR.  Take a flip mirror sky frame to continue. ", l_logfile, verbose+) 
        status=1

    } else { # no sky image

    if (offonly == yes) {

        nx = 0
        ny = 0
        hselect (outimg//"[1]", "i_naxis1", yes) | scan (nx)
        hselect (outimg//"[1]", "i_naxis2", yes) | scan (ny)
        printlog ("--------------------------------------------------------\
            ----------------", l_logfile, verbose+) 

        imcopy (outimg//"[1]", tmp_outimg , >& "dev$null")

    } else {   #offonly=no

        imcopy (inimg//"[1]", tmpout1, >& "dev$null")

    #####
    # Pad the NIFS frames in the y-direction

        nfpad (tmpout1//".fits", tmpout2, exttype="SCI", logfile=l_logfile, ver+)

    if (skyimg != "") {
        printlog ("NIFIMAGE:  Input sky file ="//skyimg, l_logfile, verbose+)
        imcopy (skyimg//"[1]", tmp_sky, >& "dev$null")
        nfpad (tmp_sky, tmpout3, exttype="SCI", logfile=l_logfile, ver+)       
    }

    printlog ("NIFIMAGE:  Input file ="//inimg, l_logfile, verbose+) 

        #####

        imcopy (tmpout2, tmp_inimg, >& "dev$null")

        #####
        # Get output image size from first input, extn 1

        eximg = ""
        head (tmpin) | scan (eximg)

        nx = 0
        ny = 0
        hselect (tmp_inimg, "i_naxis1", yes) | scan (nx)
        hselect (tmp_inimg, "i_naxis2", yes) | scan (ny)
        if ((0 == nx) || (0 == ny)) {
            printlog ("ERROR - NFACQUIRE: Cannot get size of first image \
                extension " // eximg // ".", l_logfile, verbose+) 
            goto clean
        }
        
        #####
        #Get MDF for for the bottom full slice:

        tprint (mdf, prparam-, prdata+, showrow-, showhdr-, showunits-,
            col="x_ccd,y_ccd,slitsize_x,slitsize_y,specorder,corner,RA,DEC",
            rows=row, pwidth=160) | scan (x_ccd, y_ccd, slitsize_x, 
            slitsize_y, specord, xdcorner, ra, dec)

        slicewidth = slitsize_x / l_yscale
        sliceheight = slitsize_y / l_yscale

        # Get centring interactively, if required

        hselect (images=phu, fields="INSTRUME", expr=yes) | scan (instrument)

        #below is a copy of the imexam call for NIFS

        if (l_fl_inter) {

            havex = no
            havey = no
            
            while (no == havex || no == havey) {

                printlog ("--------------------------------------------------\
                    ----------------------", l_logfile, verbose+)

                printlog ("NFACQUIRE: INSTRUCTIONS:", l_logfile, verbose+)
                printlog (" ", l_logfile, verbose+)
                printlog ("NFACQUIRE: Go to the BOTTOM LEFT side of the NIFS \
                    image", l_logfile, verbose+) 
                printlog (" where the image slicer pattern is.", l_logfile, verbose+) 
                printlog (" ", l_logfile, verbose+)
                printlog ("NFACQUIRE: for interactive extraction of the \
                    slices,", l_logfile, verbose+) 
                printlog (" follow the below directions for the lowest \
                    *full", l_logfile, verbose+) 
                printlog (" slice* on the NIFS detector - the bottom slice", l_logfile, verbose+)
                printlog (" is only half there.", l_logfile, verbose+) 
                printlog (" ", l_logfile, verbose+)
                printlog ("NFACQUIRE: The slice that you should define using j \
                    and k ", l_logfile, verbose+) 
                printlog (" keys in imexam is identified by the BLUE BOX.", l_logfile, verbose+) 
                if (no == havex)
                printlog (" ", l_logfile, verbose+)
                    printlog ("NFACQUIRE: press j on the upper slice in the blue box \
                        to ", l_logfile, verbose+) 
                printlog (" identify the slice in the x direction.", l_logfile, verbose+) 
                if (no == havey)
                printlog (" ", l_logfile, verbose+)
                    printlog ("NFACQUIRE: press k in the central gap between the two slices in \
                        the ", l_logfile, verbose+)
                printlog (" blue box.  This will identify the slice pattern \
                        in the", l_logfile, verbose+) 
                printlog (" y direction.", l_logfile, verbose+) 
                printlog (" ", l_logfile, verbose+)
                    printlog ("NFACQUIRE: press q to quit", l_logfile, verbose+)

                printlog ("--------------------------------------------------\
                    ----------------------", l_logfile, verbose+) 

                display (tmp_inimg, l_raw_frame, >& "dev$null")

        print (str(xtest1), " ", str(ytest1), "B", >> tmpfile)
    r=0.50
        tvmark.lengths = str(80)//" "//r
        tvmark (l_raw_frame, tmpfile, mark="rectangle", color=206, int-, label-)
                imexam (tmp_inimg, frame=l_raw_frame, image="", logfile=tmplog,
                    keeplog+, >& "dev$null")

                scanin1 = tmplog
                while (fscan (scanin1, sline) != EOF) {
                    if (1 == stridx (sline, "Column")) {
                        havey = yes
                        junk = fscan (sline, tmp1, tmp2, tmp3, tmp4)
                        if (strlen (tmp3) <= 8)
                            sline = tmp4
                        else
                            sline = substr (tmp3, 8, strlen (tmp3))
                        junk = fscan (sline, l_shifty)
                    }
                    if (1 == stridx (sline, "Line")) {
                        havex = yes
                        junk = fscan (sline, tmp1, tmp2, tmp3, tmp4)
                        if (strlen (tmp3) <= 8)
                            sline = tmp4
                        else
                            sline = substr (tmp3, 8, strlen (tmp3))
                        junk = fscan (sline, l_shiftx)
                    }
                }


                delete (tmplog, ver-, >& "dev$null")
                ytest = abs(l_shifty-y_ccd)
                if (ytest > 70.0) {
                    havey = no
                    printlog ("NFACQUIRE: The requested y offset is greater \
                        than the size of a slice.", l_logfile, verbose+)
                    printlog ("           To try again, please press <space> \
                        and make sure you are centered", l_logfile, verbose+)
                    printlog ("           on the bottom-most full slice.",
                        l_logfile, verbose+)
                    printlog ("       OR  Press 'q' in the image display to \
                        abort NFACQUIRE", l_logfile, verbose+)
                    junk = fscan (imcur, rjunk, rjunk, junk, cursinput)
                    if (cursinput == "q") {
                        printlog ("NFACQUIRE: Aborting per user request.",
                            l_logfile, verbose+)
                        status = 1
                        goto clean
                    }
                }     
            }

            printlog ("------------------------------------------------------\
                ------------------", l_logfile, verbose+) 

            xtest = l_shiftx
            ytest = l_shifty

            if (2 == l_dispaxis) l_shifty = y_ccd - l_shifty
            else                 l_shiftx = (l_shiftx - x_ccd)

            if (2 == l_dispaxis) mdfrow = 1
            else                 mdfrow = 28

        }

        # Generate the images - do sky subtraction, if requested

      imdelete (tmpout2, ver-, >& "dev$null")
      tmpout2 = mktemp ("tmpout2")

if (skyimg != "") {
        printlog("NFACQUIRE: Doing the sky subtraction", l_logfile, l_verbose)
        imarith (tmp_inimg, "-", tmpout3, tmpout2, >& "dev$null")
        imcopy (tmpout2//".fits", tmpout4, >& "dev$null")
} else {
         imcopy (tmp_inimg//".fits", tmpout4, >& "dev$null")
}

        firstimg = yes

        mkimage (tmp_outimg, option = "make", value = 0.0, ndim = 2,
            dims = "145 142", pixtype = "real", slope = 0.0, sigma = 0.0,
            seed = 0)

        # Get nearest appropriate aperture shift

        if (l_fl_inter && firstimg) {

            tdump (l_mdf, cdfile="", pfile="", datafile=tmpmdf, columns="y_ccd",
                rows="-", pwidth=-1)

            scanin1 = tmpmdf
            havex = no
            first = yes
            if (2 == l_dispaxis) dmin = l_shiftx
            else                 dmin = l_shifty

            while (! havex && fscan (scanin1, x1) != EOF) {

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


            l_shiftx = xtest - x_ccd
            l_shifty = ytest - y_ccd

            printlog ("NFACQUIRE: Offsets from slice center in MDF: x: "//\
                 l_shiftx// "; y: "//l_shifty, l_logfile, l_verbose) 
            firstimg = no
        }

        first = yes
        i = 1
        j = 1
        while (i <= 29) {
            tprint (mdf, prparam-, prdata+, showrow-, showhdr-, showunits-,
                col="x_ccd,y_ccd,slitsize_x,slitsize_y,specorder,corner,RA,DEC",
                rows=i, pwidth=160) | scan (x_ccd, y_ccd, slitsize_x, 
                slitsize_y, specord, xdcorner, ra, dec)
            i = i+1

            ix1 = x_ccd
            iy1 = y_ccd
            idx = slicewidth
            idy = sliceheight

            if ((iy1+l_shifty) < 1)
                iy1 = abs(l_shifty) + 2

            src = "[" // int (ix1+l_shiftx) // ":" \
                // int (ix1+idx+l_shiftx) \
                // "," \
                // int (iy1+l_shifty) // ":" \
                // int (iy1+idy+l_shifty) // "]"

            src = tmpout4 // src

            idy2 = idy
            idx2 = 0
            tmpwhole = mktemp ("tmpwhole")
            tmpthin = mktemp ("tmpthin")
            tmpstack = mktemp ("tmpstack")
            tmpimg = mktemp ("tmpimg")
            imcopy (src, tmpwhole, >& "dev$null")
            improject (tmpwhole, tmpthin, projaxis=l_dispaxis, average+,
                highcut=0.0, lowcut=0.0, pixtype="real", verbose-,
                >& "dev$null")
            if (1 == l_dispaxis) {
                imstack (tmpthin, tmpstack, title="", pixtype="*")
                imdelete (tmpthin, verify-, >& "dev$null")
                tmpthin = mktemp ("tmpthin")
                imtranspose (tmpstack, tmpthin)
                blkrep (tmpthin, tmpthin, 5, 2)
                imcopy (tmpthin, tmpimg, >& "dev$null")
            }

            idx2 = 4
            idy2 = (2*idy)+1

            ix2 = j
            j = j+idx2
            iy2 = 1
            dst = "[" // ix2 // ":" // ix2 + idx2 // "," \
                // iy2 // ":"  // iy2 + idy2 // "]"

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

            imcopy (tmpimg, tmp_outimg // dst, >& "dev$null")
            imdelete (tmpwhole, ver-, >& "dev$null")
            imdelete (tmpthin, ver-, >& "dev$null")
            imdelete (tmpstack, ver-, >& "dev$null")
            imdelete (tmpimg, ver-, >& "dev$null")
            j = j+1
        }

        # flip the p-axis in the output cube because it is backwards
 
            nx=xmax-xmin+1.0
            ny=ymax-ymin+1.0
        
         imdelete (tmpout, ver-, >& "dev$null")
         imdelete (tmpout2, ver-, >& "dev$null")
      tmpout2 = mktemp ("tmpout2")
      imcopy(tmp_outimg//"[-*,*]",tmpout2, ver-, >& "dev$null")
          imdelete (tmp_outimg, ver-, >& "dev$null")
      tmp_outimg = mktemp ("tmp_outimg")
      imcopy(tmpout2,tmp_outimg, ver-, >& "dev$null")
      imdelete (tmpout2, ver-, >& "dev$null")

   } # offonlyj

        # Give shift info for star centre

        if (l_fl_shift) {

            printlog ("-----------------------------------------------------\
                -------------------", l_logfile, verbose+) 

            havexy = no
            while (no == havexy) {

                printlog ("NFACQUIRE: press a or x to identify star and then \
                    q to quit", l_logfile, verbose+) 
                display (tmp_outimg, l_img_frame, >& "dev$null")
                imexam (tmp_outimg, frame=l_img_frame, image="", logfile=tmplog,
                    keeplog+, >& "dev$null")

                scanin1 = tmplog
                x1 = -1.0
                y1 = -1.0
                while (! havexy && (fscan (scanin1, x1, y1) != EOF)) {
                    if (x1 >= 0.0 && y1 >= 0.0) {
                        havexy = yes
                    } else {
                        x1 = -1.0
                        y1 = -1.0
                    }
                }

                delete (tmplog, ver-, >& "dev$null")

            }

            # this is all horribly NIFS specific 

            mx = nx
            my = ny 
            
            if (l_xpos == 0.0)
                x2 = 0.5 * real (mx)
            else
                x2 = l_xpos

            if (l_ypos == 0.0)
                y2 = 0.5 * real (my)
            else
                y2 = l_ypos

            dx = x1 - x2
            dy = y1 - y2

            #instrument offsets: 
            poff = dx * 0.0218 * guidevalue
            #q offsets increase to the left (-x)
            qoff = dy * 0.0218
            printlog (" ", l_logfile, verbose+)
            print ("****************************************\
                ******************************")
            printf ("    IFU instrument offsets (arcsec)  p=%-6.3f \
                q=%-6.3f\n", poff, qoff) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            print ("**************************************\
                ********************************")
            printlog (" ", l_logfile, verbose+)
            printlog ("    Raw data source  " //inimg, l_logfile, verbose+)
            printf ("%35s  (%-6.2f, %-6.2f)\n", \
                "Expanded Image size (pixels)", (nx), \
                (ny)) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            printf ("%35s  (%-6.3f, %-6.3f)\n", \
                "Offset-to location (pixels)", x2, y2) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            printf ("%35s  (%-6.3f, %-6.3f)\n", \
                "Measured peak (pixels)", x1, y1) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            printf ("%35s  (%-6.3f, %-6.3f)\n", \
                "Offset (pixels)", dx, dy) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            printf ("%35s  (%-6.3f, %-6.3f)\n", \
                "Offset (arcsec)",dx * 0.0218 * guidevalue, dy * 0.0218) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            print ("****************************************\
                ********************************")

        } else if (l_fl_display) {
            display (tmp_outimg, l_img_frame)
        }

    if (offonly == no) {
        printlog ("---------------------------------------------------------\
            ---------------", l_logfile, verbose+) 
        printlog ("NFACQUIRE - Writing output file", l_logfile, verbose+) 
        wmef (tmp_outimg, outimg, extname="SCI", phu=phu, verbose-)
    }

    gemdate ()
    gemhedit (outimg//"[0]", "NFACQUIR", gemdate.outdate,
        "UT Time stamp for NFACQUIRE", delete-)
    gemhedit (outimg//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)
    
    status = 0    # Success

    } # no sky frame

clean:

    scanin1 = ""

    delete (tmpin, ver-, >& "dev$null")
    delete (tmpout, ver-, >& "dev$null")
    delete (tmpfile, ver-, >& "dev$null")
    imdelete (tmpout1, ver-, >& "dev$null")
    imdelete (tmpout2, ver-, >& "dev$null")
    imdelete (tmpout3, ver-, >& "dev$null")
    imdelete (tmpout4, ver-, >& "dev$null")
    imdelete (tmp_sky, ver-, >& "dev$null")
    delete (tmpmdf, ver-, >& "dev$null")
    imdelete (tmpimg, ver-, >& "dev$null")
    imdelete (tmp_inimg, ver-, >& "dev$null")
    imdelete (tmp_outimg, ver-, >& "dev$null")

    printlog ("-------------------------------------------------------------\
        -----------", l_logfile, verbose+) 

    if (0 == status) {
        printlog (" ", l_logfile, l_verbose) 
        printlog ("NFACQUIRE exit status:  good.", l_logfile, l_verbose) 
    }

end
