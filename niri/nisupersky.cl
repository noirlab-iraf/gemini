# Copyright(c) 2006-2009 Association of Universities for Research in Astronomy, Inc.

procedure nisupersky(inimage)

# Identify objects in a combined, stacked image, and extract
# object masks for use with NISKY.  It is assumed that the
# combined image was created using IMCOADD, so that the header
# contains the names of the input images.  It is also assumed
# that the offsets/WCS in the headers are correct.
#
# Possible upgrades:
#  - deal with any input image with a valid WCS, and a list of images
#    for which masks are to be made
#  - Providing input image/output mask names explicitly instead of 
#    relying on header list would be good.
#  - We need all the good pixels, not just those cropped by imcoadd
#  - currently the input combined image is assumed to be the same size
#    as the output image.  This will need to be fixed if imcoadd is
#    modified.
#  - more extensive error checking/handling. 
# 
# Version  May 27, 2005  JJ, first version


char    inimage     {prompt="Input combined NIRI image"}                # OLDP-1-input-primary
real    threshold   {3.5,min=1.5,prompt="Threshold in sigma for object detection"} # OLDP-2
int     ngrow       {3,prompt="Number of iterations to grow objects into the wings"} # OLDP-2
real    agrow       {3., prompt="Area limit for growing objects into the wings"} # OLDP-2
int     minpix      {6,prompt="Minimum number of pixels to be identified as an object"} # OLDP-3
char    sci_ext     {"SCI",prompt="Name or number of science extension"}      # OLDP-3
char    dq_ext      {"DQ",prompt="Name or number of data quality extension"}  # OLDP-3
char    logfile     {"", prompt="Logfile"}                              # OLDP-1
bool    verbose     {yes, prompt="Verbose"}                             # OLDP-4
int     status      {0, prompt="Exit status (0=good)"}                  # OLDP-4
struct  *scanfile   {prompt="Internal use"}                             # OLDP-4

begin

    char    l_inimage = ""
    char    l_logfile = ""
	char    l_sci_ext = ""
	char    l_dq_ext = ""
    bool    l_verbose
    int     l_ngrow
    int     l_minpix
    real    l_threshold
    real    l_agrow
    int     maxfiles = 200
    char    l_temp = ""
    char    tmpout, tmpmask, tmpshift, mask, tmpcoordin, istr
    int     junk
    int     i
    int     xsize, ysize
    int     x1, y1
    real    xoff, yoff
    real    raref, decref, xrefpix, yrefpix
    struct  l_struct

    status = 0
    tmpmask = mktemp("tmpmask")
    mask=""

    cache("niri", "imgets", "gemdate")

    # set the local variables
    
    junk = fscan (inimage, l_inimage)
    junk = fscan (logfile, l_logfile)
    junk = fscan (sci_ext, l_sci_ext)
    junk = fscan (dq_ext, l_dq_ext)    
	l_threshold = threshold
	l_ngrow = ngrow
	l_agrow = agrow
    l_minpix = minpix
    l_verbose = verbose
    
    #------------------------------------------------------------------------
    # Check for package log file or user-defined log file
    if (l_logfile == "") {
        junk = fscan (niri.logfile, l_logfile)
        if (l_logfile == "") {
            l_logfile="niri.log"
            printlog ("WARNING - NISUPERSKY: Both nisupersky.logfile and \
                niri.logfile are empty.", l_logfile, verbose+)
            printlog ("                     Using default file niri.log.",
                l_logfile, verbose+)
        }
    }
    # Open log file
    date | scan (l_struct)
    printlog ("-------------------------------------------------------------\
        ---------------", l_logfile, l_verbose)
    printlog ("NISUPERSKY -- "//l_struct, l_logfile, l_verbose)
    printlog(" ",l_logfile, l_verbose)

    #-----------------------------------------------------------------------
    # Check for existence of input image
    if (!imaccess(l_inimage)) {
        printlog ("ERROR - NISUPERSKY: Input image does not exist", l_logfile, verbose+)
        status = 1
        goto clean
    }

    printlog ("Processing masks for "//l_inimage, l_logfile, l_verbose)

    #--------------------------------------------------------------------------
    # Create the masks

    # First create the master mask from the combined image
    objmasks(l_inimage//"["//l_sci_ext//"]",tmpmask//".pl",omtype="boolean",
        hsigma=l_threshold,minpix=l_minpix,ngrow=l_ngrow,agrow=l_agrow,
        >& "dev$null")
    
    # Determine reference RA, dec for master mask
    tmpcoordin = mktemp("tmpcoordin")
    raref = 0.0; decref = 0.0
    hselect(tmpmask,"CRVAL1,CRVAL2", yes) | scan(raref,decref)
    if (2 != nscan ()) {
      printlog("NISUPERSKY - ERROR: No WCS info for combined image.",
          l_logfile,verbose+)
      status=1
      goto clean
    } 
    printf("%12.8f %12.8f\n", raref, decref, > tmpcoordin)
    xrefpix = 0.0; yrefpix = 0.0
    hselect(tmpmask,"CRPIX1,CRPIX2", yes) | scan(xrefpix,yrefpix)
    if (2 != nscan ()) {
      printlog("NISUPERSKY - ERROR: No WCS info for combined image.",
          l_logfile,verbose+)
      status=1
      goto clean
    } 
    
#    # Set first NISKY mask name (this is needed due to a bug in NISKY)
#    hselect(tmpmask,"NISKYMSK",yes) | scan(mask)

    # Main loop: Make the individual masks
    for(i=1; i<maxfiles; i+=1) {
      tmpout = mktemp("tmpout")
      tmpshift = mktemp("tmpshift")
      l_temp=""
      printf ("%02d\n", i) | scan (istr)
      hselect(l_inimage//"[0]","IMAGE"//istr,yes) | scan(l_temp)    
      if(l_temp!="") {
        hselect(l_temp//"[0]","NISKYMSK",yes) | scan(mask)
        printlog("   "//mask,l_logfile,l_verbose)

        # shift the mask and copy into full sized frame
        imarith(tmpmask,"*",0.,tmpshift//".pl")

        # determine image size
        xsize=0; ysize=0
        hselect (tmpmask, "naxis1", yes) | scan (xsize)
        hselect (tmpmask, "naxis2", yes) | scan (ysize)
 
        # determine offsets in pixels
        wcsctran(tmpcoordin,"STDOUT",l_temp//"[0]","world","logical",
            columns="1 2", units="", formats="", min_sigdigit=7,
            verbose-) | scan(xoff,yoff)
        if (2 != nscan ()) {
          printlog("NISUPERSKY - ERROR: No WCS info for image "//l_temp,
              l_logfile,verbose+)
          status=1
          goto clean
        } 
        xoff=xrefpix-xoff
        yoff=yrefpix-yoff
        x1=nint(xoff) ; y1=nint(yoff)

        # chop out the appropriate piece and paste into the new mask
        imcopy(tmpmask//"["//(max(x1,1))//":"//(min(x1+xsize,xsize))//","\
            //(max(y1,1))//":"//(min(y1+ysize,ysize))//"]",
            tmpshift//".pl["//(max(1-x1,1))//":"//(min(xsize-x1,xsize))//","\
            //(max(1-y1,1))//":"//(min(ysize-y1,ysize))//"]", >& "dev$null")

        if(access(mask)) {
        # combine with old mask
          addmasks(mask//","//tmpshift,tmpout//".pl","im1 || im2")

        } else if(imaccess(l_temp//"["//l_dq_ext//"]")) {
        # generate new mask, including DQ plane of original if possible
          addmasks(l_temp//"["//l_dq_ext//"],"//tmpshift,tmpout//".pl",
              "im1 || im2")
        } else {
          imcopy(tmpshift,tmpout//".pl", >& "dev$null")
        }       

        if(imaccess(mask)) {
          printlog("WARNING - NISUPERSKY: overwriting mask "//mask,
              l_logfile,l_verbose)
          imdelete(mask,ver-, >& "dev$null")
        }
        imrename(tmpout,mask)
     
        imdelete(tmpout,ver-, >& "dev$null")
        imdelete(tmpshift,ver-, >& "dev$null")

        gemdate ()
        gemhedit (mask, "NISUPSKY", gemdate.outdate, "UT Time stamp for NISUPERSKY")
        gemhedit (mask, "GEM-TLM", gemdate.outdate, "UT Last modification with GEMINI")

      } # end if
    } # end main for loop



    #---------------------------------------------------------------------------
    # Clean up
clean:
    printlog (" ", l_logfile, l_verbose)
    if (status == 0)
        printlog("NISUPERSKY exit status:  good.", l_logfile, l_verbose)
    else
        printlog("NISUPERSKY exit status:  failed.", l_logfile, l_verbose)
    printlog("----------------------------------------------------------------------------",l_logfile,l_verbose)

    scanfile = ""
    imdelete (tmpmask, ver-, >& "dev$null")
    delete (tmpcoordin, ver-, >& "dev$null")
end


