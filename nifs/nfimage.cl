# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

# NFIMAGE - Collapse NIFS dispersed science spectral frames into 2D images
#           of the field (for quick-look of unprocessed data).
#
# Original author: Tracy Beck
#
# Warning: very NIFS-specific

procedure nfimage (image)

int     image       {prompt = "Acquisition image (int range)"}
char    prefix    {"default", prompt="Prefix of image"}
char    rawpath     {"adata$", prompt = "Data directory"}
int     sky         {0,  prompt ="Acquisition sky image (int range)"}
char    outimage    {"", prompt = "Output file"}
char    outprefix   {"i", prompt = "Prefix to use if outimage not given"}

char    mdf         {"nifs$data/nifs-mdf.fits", prompt = "MDF that describes mapping"}
real    shiftx      {0.0, prompt = "Shift in X"}
real    shifty      {0.0, prompt = "Shift in Y"}
real    xpos        {0.0, prompt = "Centering position in X"}
real    ypos        {0.0, prompt = "Centering position in Y"}
real    yscale      {0.043, prompt = "Size of pixels in Y (arcec)"}
int     dispaxis    {1, prompt = "Dispersion axis"}
bool    fl_inter    {no, prompt = "Get shift interactively"}
bool    fl_display  {yes, prompt = "Display each image as it is generated"}
int     raw_frame   {1, prompt = "Frame to display raw data"}
int     img_frame   {2, prompt = "Frame to display reconstructed image"}
bool    fl_shift    {yes, prompt = "Get shift to centre star in each image"}
bool    fl_imexam   {no, prompt = "Run imexam on final image"}
char    logfile     {"", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose output?"}
int     status      {0, prompt = "Exit status (0 = good)"}
struct  *scanin1    {"", prompt = "Internal use only"}

begin

    int     l_image
    char    l_prefix = ""
    char    l_rawpath = ""
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
	bool    l_fl_imexam
    char    l_logfile = ""
    bool    l_verbose
    char    l_key_pixscale

    char    lastchar
    char    rootn, inimg, skyimg, outimg
    char    dia, mes, zone, ano
    char    ddir, tmpout2
    int     junk, nx, ny, mx, my, xtest1, ytest1
    bool    debug, havex, havey, havexy, first, firstimg, offonly
    struct  sline
    char    tmpin, tmpout, src, dst, tmpfile
    char    tmpwhole, tmpthin, tmp1, tmp2, tmp3, tmp4, tmplog
    char    eximg, tmpimg, badhdr, tmpmdf, tmpstack, phu
    char    instrument, tmp_inimg, tmp_outimg, guidestate,grating
    real    x1, y1, dx, dy, x2, y2, poff, qoff, guidevalue, r
    int     ix1, iy1, idx, idy, ix2, iy2, idx2, idy2, row
    int     i, j
    int     xmin, xmax, ymin, ymax, mdfrow
    real    x_ccd, y_ccd, slitsize_x, slitsize_y, specord, xdcorner, ra, dec
    real    slicewidth, sliceheight, xtest, ytest
    real    rjunk
    char    cursinput, imstring, skystring

    cache ("gemextn", "gemdate")

    status = 1
    offonly = no
    debug = no
    badhdr = ""
    row = 28

    xtest1=1024
    ytest1=79
    
    l_image = image

    junk = fscan (prefix, l_prefix)
    junk = fscan (rawpath, l_rawpath)
    l_sky = sky
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
	l_fl_imexam = fl_imexam
    junk = fscan (logfile, l_logfile)
    l_verbose   = verbose

    # Test the logfile
    if (l_logfile == "") {
        l_logfile = nifs.logfile
        if (l_logfile == "") {
            l_logfile = "nifs.log"
            printlog ("WARNING - NFIMAGE: Both nfimage.logfile and \
                nifs.logfile are empty.", l_logfile, l_verbose)
            printlog ("                   Using default file nifs.log",
                l_logfile, l_verbose)
        }
    }


    tmpout2 = mktemp ("tmpout2")
    tmpin = mktemp ("tmpin")
    tmpfile = mktemp ("tmpfile")
    tmpout = mktemp ("tmpout")
    tmplog = mktemp ("tmplog")
    tmpmdf = mktemp ("tmpmdf")
    tmpimg = mktemp ("tmpimg")
    tmp_inimg = mktemp ("tmp_inimg")
    tmp_outimg = mktemp ("tmp_outimg")

    printlog ("------------------------------------------------------------------------", l_logfile, verbose+) 
    printlog ("TASK - NFIMAGE:  Derive quick images for dispersed NIFS data", l_logfile, verbose+) 
    printlog ("------------------------------------------------------------------------", l_logfile, verbose+) 

    # Ensure a trailing '/' on the directory name
    lastchar = substr (l_rawpath, strlen(l_rawpath), strlen(l_rawpath))
    if ((l_rawpath != "") && (lastchar != "$") && (lastchar != "/"))
        l_rawpath = l_rawpath//"/"
    ddir = l_rawpath
    
    # Build image names based on date or prefix
	printf("%04d\n",l_image) | scan(imstring)
    print (l_prefix) | lcase | scan(rootn)
    if (rootn == "default") {
        date ("-u +'%m %d %Y %Z'") | scan (mes, dia, ano, zone)
        if ((zone != "UTC") && (zone != "GMT")) {
            print ("ERROR: Internal Error")
            print ("ERROR: Date not given as universal time ("//zone//")")
            goto clean
        }

        # Build input image names
        inimg = ddir//"N"//ano//mes//dia//"S"//imstring//".fits"
    } else {
        print(rootn) | ucase | scan(rootn)
        inimg = ddir//rootn//imstring//".fits"
    }


    if (l_outimage == "") {
        if (l_outprefix == "") {
            if (rootn == "default")
                outimg = "aN"//ano//mes//dia//"S"//imstring//".fits"
            else
                outimg = "a"//rootn//imstring//".fits"
        } else {
            if (rootn == "default")
                outimg = l_outprefix//"N"//ano//mes//dia//"S"//imstring//".fits"
            else
                outimg = l_outprefix//rootn//imstring//".fits"
        }
    } else
        outimg = l_outimage
    
    # Build input sky image names
    printf("%04d\n",l_sky) | scan(imstring)
	skyimg = ""
    if (l_sky != 0) {
        if (rootn == "default") {
            skyimg = ddir//"N"//ano//mes//dia//"S"//imstring//".fits"
        } else {
            skyimg = ddir//rootn//imstring//".fits"
        }
    }

    #####
    # Check input data

    src = inimg

    gemextn (src, proc="none", check="exists", index="", extname="", extver="",
        ikparams="", omit="", replace="", outfile=tmpin, logfile=l_logfile,
        glogpars="", verbose-)
    if ((0 < gemextn.fail_count) || (0 == gemextn.count)) {
        printlog ("ERROR - NFIMAGE: Missing input data.", l_logfile, verbose+) 
        goto clean
    }
    
    printlog ("NIFIMAGE:  Input file ="//inimg, l_logfile, verbose+) 

    src = skyimg

    gemextn (src, proc="none", check="exists", index="", extname="", extver="",
        ikparams="", omit="", replace="", outfile=tmpin, logfile=l_logfile,
        glogpars="", verbose-)
    if ((0 < gemextn.fail_count) || (0 == gemextn.count)) {
        printlog ("NFIMAGE: No or Missing input sky data. Not doing sky \
            subtraction", l_logfile, verbose+) 
        skyimg = ""
    }
    

    #####

    phu = inimg // "[0]"

    #####
    # Check output data
        printlog ("---------------------------------------------------------\
            ---------------", l_logfile, verbose+) 

    gemextn (outimg, proc="none", check="absent", index="", extname="",
        extver="", ikparams="", omit="", replace="", outfile=tmpout,
        logfile="", glogpars="", verbose-)
    if (0 < gemextn.fail_count) {
        printlog ("---------------------------------------------------------\
            ---------------", l_logfile, verbose+) 
        printlog ("WARNING - NFIMAGE: Reconstructed IFU image already exists.",
            l_logfile, verbose+) 
        printlog ("WARNING - NFIMAGE: Displaying frame to generate new \
            offsets.", l_logfile, verbose+) 
        l_fl_inter = no
        offonly = yes
        printlog ("---------------------------------------------------------\
            ---------------", l_logfile, verbose+) 
    }

   #####
   # subtract sky, if given

   if (skyimg != "") {
        imarith (inimg//"[1]", "-", skyimg//"[1]", tmpout, >& "dev$null")
        printlog ("NIFIMAGE:  Input sky file ="//skyimg, l_logfile, verbose+) 
   } else {
        imcopy (inimg//"[1]", tmpout, >& "dev$null")
   }

    #####

   #####
   # Determine guide state for offset value flip.

            hselect (images=phu, fields="AOFOLD", expr=yes) | \
                scan (guidestate)

    if (guidestate == 'IN') {
       guidevalue=-1.0
        printlog ("NIFIMAGE:  NIFS Guiding with Altair AO", l_logfile, verbose+) 
    } else {
       guidevalue=1.0
        printlog ("NIFIMAGE:  NIFS Guiding without AO", l_logfile, verbose+) 
    }

   #####
   # Set the offset image if the reconstructed IFU output file already exists

    if (offonly == yes) {

        nx = 0
        ny = 0
        hselect (outimg//"[1]", "i_naxis1", yes) | scan (nx)
        hselect (outimg//"[1]", "i_naxis2", yes) | scan (ny)
        printlog ("--------------------------------------------------------\
            ----------------", l_logfile, verbose+) 

        imcopy (outimg//"[1]", tmp_outimg , >& "dev$null")

    } else {   # offonly == no

   #####
   # Reconstruct the IFU image from raw input data - copied from GNIRS
   # "nfquick".  Is there a better way to do this?  (probably!)

        nfpad (tmpout//".fits", tmpout2, exttype="SCI", logfile=l_logfile, ver+)
        delete (tmpout//".fits", ver-, >& "dev$null")

        imcopy (tmpout2//".fits", tmp_inimg, >& "dev$null")
        delete (tmpout2//".fits", ver-, >& "dev$null")

        # Get output image size from first input, extn 1

        eximg = ""
        head (tmpin) | scan (eximg)

        nx = 0
        ny = 0
        hselect (tmp_inimg, "i_naxis1", yes) | scan (nx)
        hselect (tmp_inimg, "i_naxis2", yes) | scan (ny)
        if ((0 == nx) || (0 == ny)) {
            printlog ("ERROR - NFIMAGE: Cannot get size of first image \
                extension " // eximg // ".", l_logfile, verbose+) 
            goto clean
        }

        #####
        #Get MDF shift for for the bottom full slice:

        tprint (mdf, prparam-, prdata+, pwidth=160, plength=0, showrow-, 
            orig_row+, showhdr-, showunits-,
            col="x_ccd,y_ccd,slitsize_x,slitsize_y,specorder,corner,RA,DEC",
            rows=row, option="plain", align+, sp_col="", lgroup=0) | \
            scan (x_ccd, y_ccd, slitsize_x, slitsize_y, specord, xdcorner, 
            ra, dec)

        slicewidth = slitsize_x / l_yscale
        sliceheight = slitsize_y / l_yscale

        if (debug)
            print (slicewidth, sliceheight)

        if (l_shifty != 0.0) {
             l_fl_inter=no
        }

        if (! l_fl_inter && l_shifty==0.0) {
            hselect (images=phu, fields="GRATING", expr=yes) | \
                scan (grating)

#  Note:  The below numbers need updated periodically because of shifts.
            if (grating == "K_G5605") {
	        l_shifty=-42.6
            }
            if (grating == "H_G5604") {
                l_shifty=-46.5
            }
            if (grating == "J_G5603") {
                l_shifty=-38.5
            }
            if (grating == "Z_G5602") {  
                l_shifty=-41.5
            }
  
            printlog ("NFIMAGE: Grating="//grating//", shifty ="//l_shifty,\
            l_logfile, verbose+) 
        }

        printlog ("--------------------------------------------------------\
            ----------------", l_logfile, verbose+) 

        if (l_fl_inter && shifty == 0.0) {
            # Get centring interactively, if required

            hselect (images=phu, fields="INSTRUME", expr=yes) | \
                scan (instrument)

            #below is a copy of the imexam call for NIFS

            havex = no
            havey = no
            while (no == havex || no == havey) {

                printlog ("-----------------------------------------------\
                    -------------------------", l_logfile, verbose+) 
                printlog ("NFIMAGE: for interactive extraction of the \
                    slices, follow the directions below for the \
                    *bottom-most* full slice on the NIFS detector.", 
                    l_logfile, verbose+) 
                if (no == havey) {
                    printlog ("NFIMAGE INSTRUCTIONS: Go to the \
                     bottom center of the NIFS image, there will be a \
                     blue box drawn defining the lowest gap between the \
                     IFU slices.", l_logfile, verbose+)
                    printlog ("NFIMAGE INSTRUCTIONS: Within the box, \
                       go to the gap between the slices and press k to \
                        identify the gap in the y direction and then q \
                        to quit", l_logfile, verbose+)
                    printlog ("NFIMAGE INSTRUCTIONS: For very short exposures,\
                      or in the Z and J-bands where the background is low and \
                      there aren't many arc lines, \
                      a rough guess of the gap position may be required.",l_logfile, verbose+)
                    printlog ("NFIMAGE INSTRUCTIONS: Approximate positions \
                      of the gap (pixels in y):", l_logfile, verbose+)
                    printlog ("Z-band = y~78 pix", l_logfile, verbose+)
                    printlog ("J-band = y~81 pix", l_logfile, verbose+)
                    printlog ("H-band = y~74 pix", l_logfile, verbose+)
                    printlog ("K-band = y~77 pix", l_logfile, verbose+)
                 }
                printlog ("----------------------------------------------\
                    --------------------------", l_logfile, verbose+) 

                display (tmp_inimg, l_raw_frame, >& "dev$null")

          print (str(xtest1), " ", str(ytest1), "B", >> tmpfile)

	r=0.10
        tvmark.lengths = str(400)//" "//r
        tvmark (l_raw_frame, tmpfile, mark="rectangle", color=206, int-, label-)

                imexam (tmp_inimg, frame=l_raw_frame, image="", \
                    logfile=tmplog, keeplog+, >& "dev$null")

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
                }
                havex = yes
                l_shiftx = 0

                delete (tmplog, ver-, >& "dev$null")
                ytest = abs(l_shifty - y_ccd)

                if (ytest > 70.0) {
                    havey=no
                    printlog ("NFIMAGE: The requested y offset is \
                        greater than the size of a slice.", l_logfile, verbose+)
                    printlog ("         To try again, please press <space> and \
                        make sure you are centered", l_logfile, verbose+)
                    printlog ("         on the bottom-most full slice",
                        l_logfile, verbose+)
                    printlog ("     OR  Press 'q' in the image display to \
                        abort NFIMAGE", l_logfile, verbose+) 
                    junk = fscan (imcur, rjunk, rjunk, junk, cursinput)
                    if (cursinput == "q") {
                        printlog ("NFIMAGE: Aborting per user request.",
                            l_logfile, verbose+)
                        status = 1
                        goto clean
                    }
                }     
            }

            printlog ("-------------------------------------------------\
                -----------------------", l_logfile, verbose+) 

            l_shifty = (l_shifty - y_ccd)

            xtest = l_shiftx
            mdfrow = 28

        }

        #####
        # Generate the reconstructed image

        firstimg = yes

        mkimage (tmp_outimg, option = "make", value = 0.0, ndim = 2,
            dims = "145 142", pixtype = "real", slope = 0.0,
            sigma = 0.0, seed = 0)


        # Get nearest appropriate aperture shift

        if (firstimg) {

            tdump (l_mdf, cdfile="", pfile="", datafile=tmpmdf,
                columns="y_ccd", rows="-", pwidth=-1)

            printlog ("NFIMAGE: Offset from slice center in MDF: y: "//l_shifty,
                l_logfile, l_verbose) 
        }

        first = yes
        i = 1
        j = 1
        while (i <= 29) {
            tprint (mdf, prparam-, prdata+, pwidth=160, plength=0, showrow-,
                showhdr-, showunits-,
                col="x_ccd,y_ccd,slitsize_x,slitsize_y,specorder,corner,RA,DEC",
                rows=i, option="plain", align+, sp_col="", lgroup=0) | \
                scan (x_ccd, y_ccd, slitsize_x, slitsize_y, specord, xdcorner,
                ra, dec)
            i = i+1

            ix1 = 1
            iy1 = y_ccd
            idx = (nx-1)
            idy = sliceheight

            if ((iy1+l_shifty) < 1)
                iy1 = (abs(l_shifty)+2)

            src = "[" // int (ix1+l_shiftx) // ":" \
                // int (ix1+idx+l_shiftx) \
                // "," \
                // int (iy1+l_shifty) // ":" \
                // int (iy1+idy+l_shifty) // "]"
            printlog (src, l_logfile, l_verbose)

            src = tmp_inimg // src

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
            idy2 = (2*idy) + 1

            ix2 = j
            j = j + idx2
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
        
        #####
        # flip the p-axis in the output cube because it is backwards
 
            nx=xmax-xmin+1.0
            ny=ymax-ymin+1.0

        delete (tmpout, ver-, >& "dev$null")
        imcopy (tmp_outimg//"[-*,*]", tmpout2, ver-, >& "dev$null")
        delete (tmp_outimg//".fits", ver-, >& "dev$null")
        imcopy (tmpout2, tmp_outimg, ver-, >& "dev$null")
        delete (tmpout2//".fits", ver-, >& "dev$null")

      } #end else for offonly - needed to merge offset sections -tlb 01/24/06


        ##### 
        # Get shift info for target centering

        if (l_fl_shift) {

            printlog ("-----------------------------------------------------\
                -------------------", l_logfile, verbose+) 

            havexy = no
            while (no == havexy) {

                printlog ("NFIMAGE: press a or x to identify star and then \
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
            poff = dx * 0.0218*guidevalue
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
            printf ("%35s  (%-6.2f, %-6.2f)\n", "Expanded Image size (pixels)",
                nx, ny) | scan (sline)
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
                "Offset (arcsec)", dx * 0.0218*guidevalue, dy * 0.0218) | \
                scan (sline)
            printlog (sline, l_logfile, verbose+)
            print ("****************************************\
                ********************************")

        } else if (l_fl_display) {
            display (tmp_outimg, l_img_frame)
		    if (l_fl_imexam) {
			    imexam()
		    }
        }


    if (offonly == no) {
        printlog ("----------------------------------------------------------\
            --------------", l_logfile, verbose+) 
        printlog ("NFIMAGE - Writing output file: "//outimg,
            l_logfile, verbose+) 
        wmef (tmp_outimg, outimg, extname="SCI", phu=phu, verbose-)
    }
    
    gemdate ()
    gemhedit (outimg//"[0]", "NFIMAGE", gemdate.outdate,
        "UT Time stamp for NFIMAGE", delete-)
    gemhedit (outimg//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)
    
    status = 0	# Success

clean:

    scanin1 = ""

    delete (tmpfile, ver-, >& "dev$null")
    delete (tmpin, ver-, >& "dev$null")
    delete (tmpout, ver-, >& "dev$null")
    delete (tmpout2, ver-, >& "dev$null")
    delete (tmpmdf, ver-, >& "dev$null")
    imdelete (tmpimg//".fits", ver-, >& "dev$null")
    imdelete (tmp_inimg//".fits", ver-, >& "dev$null")
    imdelete (tmp_outimg//".fits", ver-, >& "dev$null")

    printlog ("-------------------------------------------------------------\
        -----------", l_logfile, verbose+) 

    if (0 == status) {
        printlog (" ", l_logfile, l_verbose) 
        printlog ("NFIMAGE exit status:  good.", l_logfile, l_verbose) 
    }

end
