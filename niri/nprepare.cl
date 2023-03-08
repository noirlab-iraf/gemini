# Copyright(c) 2001-2011 Association of Universities for Research in Astronomy, Inc.

procedure nprepare(inimages)

# Procedure to take raw NIRI data, with a single MEF extension,
# and create the preliminary VAR and DQ frames.
#
# Data is "fixed" to account for low-noise reads, digital averages,
# and coadds as follows:
#    for data before Nov. 2001, output = input / (n reads)
#    noise = (value from nprepare.dat) * sqrt(n coadds)
#    exptime = input * (n coadds)
#    gain = 11.5 e-
#    saturation = x * (n coadds) where x=200,000 e- for bias=-0.6,
#                                        280,000 e- for bias=-0.9
#    non-linear for > 0.7*saturation
#    (these are nominal values; the data file must contain entries with
#    the names readnoise, gain, shallowbias, shallowwell, deepbias, deepwell
#    which are used in the appropriate places in the calculations above).
#
# The variance frame is generated as:
#	var = (read noise/gain)^2 + max(data,0.0)/gain
#
# The preliminary DQ frame is constructed by using the bad pixel
# mask to set the 1 bit and saturation level to set the 4 bit.
#
# Also sets FILTER keyword by alphabetically sorting the contents
# of header keywords FILTER1/2/3 for NIRI (ignoring "open").
# Gemini filter IDs are removed before parsing.
# The SPECSECn header keywords are set based on a lookup table
# of slits in the focal plane mask and CAMERA setting.
#
# The zero point of the WCS is fixed for subarray readouts by subtracting
# half the difference in size between the subarray and the full array
# from CRPIX1 and CRPIX2.
#
# Version Jun 19, 2001  JJ, beta-release
#           5 Jul 2001  fixed input list + rawpath parsing (sed's back), JJ
#          15 Oct 2001  JJ, v1.2 release
#          16 Nov 2001  JJ, version to handle fixed NIRI LNRS scaling;
#                       disallowed header updates when already NPREPARED
#                       added EXTVER to headers of output images
#          19 Nov 2001  JJ, combined two loops, changed header access
#                       to output to improve speed
#          20 Nov 2001  JJ, combined all header output to one mkheader call
#          12 Feb 2002  JJ, header output file not being cleared, fixed.
#                       Also fixed bpm=none bug.
#          25 Feb 2002  MT,TB,IJ fix bug, if(oldtimer) but had lnrs=1, then 
#                       need to copy image again. 
#         Feb 28, 2002  JJ, v1.3 release
#          24 May 2002  JJ, check size of BPM
#          12 Jun 2002  JJ, remove multiple identical input files from lists
#          10 Jul 2002  JJ, added trimming of specsecs depending on ROI
#                           NOTE that only centered ROIs are supported!
#          Sep 10, 2002  IJ parameter encoding
#          Sept 20, 2002 JJ  v1.4 release
#          Aug 18, 2003  KL IRAF2.12 - new parameter, addonly, in hedit
#          Jul 6, 2004   JJ added medium and low readnoise modes to nprepare.dat
#                           left high readnoise as "readnoise" for back compatibility
#          Nov 23, 2004  JJ added fix to WCS zero point for subarrays prior to Feb. 2005

char  inimages     {prompt="Input NIRI image(s)"}               # OLDP-1-input-primary-single-prefix=n
char  rawpath      {"",prompt="Path for input raw images"}      # OLDP-4
char  outimages    {"",prompt="Output image(s)"}                # OLDP-1-output
char  outprefix    {"n",prompt="Prefix for output image(s)"}    # OLDP-4
char  bpm          {"",prompt="Bad pixel mask file"}            # OLDP-1-input
char  logfile      {"",prompt="Logfile"}                        # OLDP-1
char  sci_ext      {"SCI",prompt="Name or number of science extension"}     # OLDP-3
char  var_ext      {"VAR",prompt="Name or number of variance extension"}    # OLDP-3
char  dq_ext       {"DQ",prompt="Name or number of data quality extension"} # OLDP-3
bool  fl_vardq     {yes,prompt="Create variance and data quality frames?"}  # OLDP-3
char  key_ron      {"RDNOISE",prompt="New header keyword for read noise (e-)"}  # OLDP-3
char  key_gain     {"GAIN",prompt="New header keyword for gain (e-/ADU)"}   # OLDP-3
char  key_sat      {"SATURATI",prompt="New header keyword for saturation (ADU)"}  # OLDP-3
char  key_nonlinear {"NONLINEA",prompt="New header keyword for non-linear regime (ADU)"}  # OLDP-3
char  key_filter   {"FILTER",prompt="New header keyword for parsed filter"} # OLDP-3
char  specsec      {"[1:1024,1:1024]",prompt="Default region of image to use for spectra when database cannot be found"}  # OLDP-3
char  database     {"niri$nprepare.dat",prompt="NPREPARE database to use"}  # OLDP-3
bool  verbose      {yes,prompt="Verbose"}            # OLDP-4
int   status       {0,prompt="Exit status (0=good)"} # OLDP-4
struct *scanfile   {"",prompt="Internal use only"}   # OLDP-4

begin

    char   l_inimages, l_outimages, l_filter, l_filter1, l_filter2, l_filter3
    char   l_rawpath, l_prefix, l_logfile, l_temp, tmpfile, tmpfile2
    char   in[1000], inpath[1000], out[1000]
    char   tmpsci, tmpvar, tmpdq, tmpinimage, keyfound
    int    i, nimages, noutimages, maxfiles, nbad
    char   l_key_ron, l_key_gain, l_key_filter, l_sci_ext, l_var_ext, l_dq_ext
    char   l_key_sat, l_key_nonlinear, l_bpm, l_pupil
    char   l_varexpression, l_tempexpression, l_database, l_fpmask, l_camera
    char   l_specsec, l_specsec1, l_specsec2, l_specsec3, tmphead
    real   l_ron, l_ronref, l_ronmedref, l_ronlowref, l_gain
    int    l_lnrs, l_coadds, l_ndavgs, l_sat, l_linlimit
    real   l_vdduc, l_vdet, l_biasvolt, l_exptime
    real   l_deepwell, l_shallowwell, l_deepbias, l_shallowbias
    real   l_cd11, l_cd12, l_cd21, l_cd22, l_pixscale, l_linear
    bool   l_verbose, l_fl_vardq, alreadyfixed[1000], bad
    bool   oldtimer, fixwcs
    int    bpmxsize, bpmysize, xsize, ysize, sx1, sx2, sy1, sy2
    real   xcen, ycen
    struct l_struct

    status=0
    maxfiles=1000
    # cache imgets - used throughout the script
    cache("imgets", "gemdate")

    # set the local variables
    l_inimages=inimages ; l_outimages=outimages ; l_rawpath=rawpath
    l_verbose=verbose ; l_prefix=outprefix ; l_logfile=logfile
    l_key_ron=key_ron ; l_key_gain=key_gain ; l_key_filter=key_filter
    l_key_sat=key_sat ; l_key_nonlinear=key_nonlinear ; l_fl_vardq=fl_vardq
    l_sci_ext=sci_ext ; l_var_ext=var_ext ; l_dq_ext=dq_ext ; l_bpm=bpm
    l_database=database ; l_specsec=specsec

    # open temp files
    tmpfile = mktemp("tmp1")
    tmpfile2 = mktemp("tmp2")
    tmphead = mktemp("tmphead")
    
    # assign other tmp files variables dummy names.
    # (Just need to give the variable a value. Real names will be created later)
    tmpsci = "dummy"
    tmpvar = "dummy"
    tmpdq = "dummy"
    tmpinimage = "dummy"

    #----------------------------------------------------------------------

    # Test that name of logfile makes sense
    cache("niri")
    print(l_logfile) | scan(l_logfile)
    if (l_logfile == "" || l_logfile == " ") {
        l_logfile = niri.logfile
        print (l_logfile) | scan (l_logfile)
        if ((l_logfile=="") || (l_logfile==" ")) {
            l_logfile = "niri.log"
            printlog ("WARNING - NPREPARE: Both nprepare.logfile and \
                niri.logfile are empty.",logfile=l_logfile, verbose+)
            printlog("                  Using default file niri.log.",
                logfile=l_logfile, verbose+)
        }
    }
    
    # Open log file
    date | scan(l_struct)
    printlog ("--------------------------------------------------------------\
        --------------", l_logfile, l_verbose)
    printlog ("NPREPARE -- "//l_struct, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    #----------------------------------------------------------------------
    # Set a few array-dependent characteristics (must be in nprepare.dat)

    l_deepwell=0. ; l_shallowwell=0. ; l_deepbias=0. ; l_shallowbias=0.0
    l_ronref=0. ; l_ronmedref=0. ; l_ronlowref=0. ; l_gain=0. ; l_linear=0.

    if (!access(l_database)) {
        printlog("ERROR - NPREPARE: Database file "//l_database//" not found.",
	    l_logfile,verbose+)
        status=1
        goto clean
    } else {
        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
	    match("readnoise",stop-) | scan(l_temp,l_ronref)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
	    match("medreadnoise",stop-) | scan(l_temp,l_ronmedref)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
	    match("lowreadnoise",stop-) | scan(l_temp,l_ronlowref)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
	    match("gain",stop-) | scan(l_temp,l_gain)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
	    match("shallowwell",stop-) | scan(l_temp,l_shallowwell)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
	    match("deepwell",stop-) | scan(l_temp,l_deepwell)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
	    match("shallowbias",stop-) | scan(l_temp,l_shallowbias)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
	    match("deepbias",stop-) | scan(l_temp,l_deepbias)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
	    match("linearlimit",stop-) | scan(l_temp,l_linear)
    }

    if (l_ronref==0. || l_ronmedref==0. || l_ronlowref==0. || l_gain==0. \
          || l_deepwell==0. || l_shallowwell==0. \
	  || l_deepbias==0. || l_shallowbias==0. || l_linear==0.) {
        printlog("ERROR - NPREPARE: Array characteristic entry not found in",
	    l_logfile,verbose+) 
        printlog("                  data file "//l_database,l_logfile,verbose+)
        status=1
        goto clean
    }

    #----------------------------------------------------------------------
    # Check the rawpath name for a final /
    if (substr(l_rawpath,(strlen(l_rawpath)),(strlen(l_rawpath))) != "/")
        l_rawpath=l_rawpath//"/"
    if (l_rawpath=="/" || l_rawpath==" ")
        l_rawpath=""

    #-----------------------------------------------------------------------
    # Load up arrays of input name lists
    # This version handles both *s, / and commas in l_inimages

    # check that list file exists
    if (substr(l_inimages,1,1)=="@") {
        l_temp=substr(l_inimages,2,strlen(l_inimages))
        if (!access(l_temp) && !access(l_rawpath//l_temp)) {
	    printlog("ERROR - NPREPARE:  Input file "//l_temp//" not found.",
	        l_logfile,verbose+)
	    status=1
	    goto clean
        }
    }

    # parse wildcard and comma-separated lists
    if (substr(l_inimages,1,1)=="@") 
        scanfile=substr(l_inimages,2,strlen(l_inimages))
    else {
        print(l_inimages, >tmpfile)
        sed("-e",'s/\,/\ /g',tmpfile) | words("STDIN", >tmpfile2)
        scanfile=tmpfile2
        delete(tmpfile,verify-, >>& "dev$null")
    }  
    while (fscan(scanfile,l_temp) != EOF) {
        files(l_rawpath//l_temp,sort-) | unique("STDIN", >> tmpfile)
    }
    delete(tmpfile2,verify-, >>& "dev$null")
    scanfile=tmpfile

    nimages=0
    nbad=0
    while (fscan(scanfile,l_temp) != EOF) {
        bad=no
        # remove rawpath
        l_temp=substr(l_temp,strlen(l_rawpath)+1,strlen(l_temp))

        # remove .fits if present
        if (substr(l_temp,strlen(l_temp)-4,strlen(l_temp)) == ".fits")
	    l_temp=substr(l_temp,1,(strlen(l_temp)-5))

        # check to see if the file exists (or if it's not MEF FITS)
        gimverify(l_rawpath//l_temp)
        if (gimverify.status==1) {
	    printlog("ERROR - NPREPARE: Input image "//l_rawpath//l_temp//" not found.",
	        l_logfile,verbose+)
	    nbad+=1
	    bad=yes
        } else if (gimverify.status>1) {
	    printlog("ERROR - NPREPARE: Input image not multiple-extension FITS.",
	        l_logfile,verbose+)
	    nbad+=1
	    bad=yes
        } else {
	    nimages=nimages+1
	    if (nimages > maxfiles) {
	        printlog("ERROR - NPREPARE: Maximum number of input images exceeded ("//str(maxfiles)//")",
		    l_logfile,verbose+)
	        status=1
	        goto clean
	    }
	    # check to see if already fixed with nprepare
	    if (!bad) {
	        keyfound=""
		hselect(l_rawpath//l_temp//"[0]","*PREPAR*",yes) | \
		    scan(keyfound)
		if (keyfound != "") {
                    alreadyfixed[nimages]=yes
                    printlog("WARNING - NPREPARE: Image "//l_rawpath//l_temp//" already fixed using *PREPARE.",
		        l_logfile,verbose+) 
                    printlog("                    Data scaling, correction for number of reads, and header",
		        l_logfile,verbose+)
                    printlog("                    updating not performed.",
		        l_logfile,verbose+)
	        } else
	            alreadyfixed[nimages]=no
	        in[nimages]=l_temp
	    } # end if(!bad)
	    # trim off path from in[i], if present
	    i=strlen(in[nimages])-1
            inpath[nimages] = ""
	    while (i>1) {
	        if (substr(in[nimages],i,i) == "/") {
                    inpath[nimages]=substr(in[nimages],1,i)
                    in[nimages]=substr(in[nimages],i+1,strlen(in[nimages]))
                    i=0
	        }
	        i=i-1
	    } # end while(i>1) path-chopping loop
        } # end else
    } # end while(fscan) loop

    # Exit if problems found
    if (nbad > 0) {
        printlog("ERROR - NPREPARE: "//nbad//" image(s) either do not exist or are not MEF files.",
	    l_logfile,verbose+)
        status=1
        goto clean
    }

    printlog("Processing "//nimages//" files",l_logfile,l_verbose)
    scanfile="" ; delete(tmpfile,ver-, >>& "dev$null")

    # Now for the output images
    # outimages could contain legal * if it is of a form like %st%stX%*.imh

    noutimages=0
    nbad=0
    print(l_outimages) | scan(l_outimages)
    if (l_outimages!="" && l_outimages!=" ") {
        if (substr(l_outimages,1,1)=="@")
	    scanfile=substr(l_outimages,2,strlen(l_outimages))
        else { 
	    files(l_outimages,sort-) | unique("STDIN", > tmpfile)
	    scanfile=tmpfile
        }

        while (fscan(scanfile,l_temp) != EOF) {
            # remove .fits if present
	    if (substr(l_temp,strlen(l_temp)-4,strlen(l_temp)) == ".fits")
	        l_temp=substr(l_temp,1,(strlen(l_temp)-5))
	        noutimages=noutimages+1
	        if (noutimages > maxfiles) {
	            printlog("ERROR - NPREPARE: Maximum number of output images exceeded ("//str(maxfiles)//")",
		        l_logfile,verbose+ )
	            status=1
	            goto clean
	        }
	        out[noutimages]=l_temp 
	        if (imaccess(out[noutimages])) {
	            printlog("ERROR - NPREPARE: Output image "//out[noutimages]//" exists",
                        l_logfile,verbose+)
	            nbad+=1
	        }
        } # end while
    }
    scanfile="" ; delete(tmpfile,ver-, >>& "dev$null")

    # Exit if problems found
    if (nbad > 0) {
        printlog("ERROR - NPREPARE: "//nbad//" image(s) already exist.",
	    l_logfile,verbose+)
        status=1
        goto clean
    }

    # if there are too many or too few output images, and any defined
    # at all at this stage - exit with error
    if (nimages!=noutimages && l_outimages!="") {
        printlog("ERROR - NPREPARE: Number of input and output images unequal.",
	    l_logfile,verbose+)
        status=1
        goto clean
    }

    # If prefix is to be used instead
    if (l_outimages=="" || l_outimages==" ") {
        print(l_prefix) | scan(l_prefix)
        if (l_prefix=="" || l_prefix==" ") {
	    printlog("ERROR - NPREPARE: Neither output image name nor output prefix is defined.",
	        l_logfile,verbose+)
	    status=1
	    goto clean
        }
        i=1
        nbad=0
        while (i<=nimages) {
	    out[i]=l_prefix//in[i]
	    if (imaccess(out[i])) {
	        printlog("ERROR - NPREPARE: Output image "//out[i]//" already exists.",
                    l_logfile,verbose+)
	        nbad+=1
	    }
	    i=i+1
        }
        if (nbad > 0) {
	    printlog("ERROR - NPREPARE: "//nbad//" image(s) already exist.",
	        l_logfile,verbose+)
	    status=1
	    goto clean
        }
    }

    #-------------------------------------------------------------------------
    # Check for existence of bad pixel mask
    # still need to add check for .pl or .fits, and convert .pl to fits
    # if required by addmasks

    if (l_fl_vardq && !imaccess(l_bpm) && l_bpm!="" && stridx(" ",l_bpm)<=0) {
        printlog("WARNING - NPREPARE: Bad pixel mask "//l_bpm//" not found.",
	    l_logfile,verbose+)
        printlog("                    Only saturated pixels will be flagged in the DQ frame.",
	    l_logfile,verbose+)
        l_bpm=""
    } else if (l_fl_vardq && (l_bpm=="" || stridx(" ",l_bpm)>0)) {
        printlog("WARNING - NPREPARE: Bad pixel mask is either an empty string or contains",
	    l_logfile,verbose+)
        printlog("                    spaces.  Only saturated pixels will be flagged in the",
	    l_logfile,verbose+)
        printlog("                    DQ frame.",l_logfile,verbose+)
        l_bpm=""
    }
    if (l_bpm=="") l_bpm="none"

    # Check to make sure BPM and input images are same size
    if (l_bpm!="none") {
        imgets(l_bpm,"i_naxis1", >& "dev$null")
        if (imgets.value!="0")
	    bpmxsize=int(imgets.value)
        imgets(l_bpm,"i_naxis2", >& "dev$null")
        if (imgets.value!="0")
	    bpmysize=int(imgets.value)
        imgets(l_rawpath//inpath[1]//in[1]//"[1]","i_naxis1")
        if (int(imgets.value)==bpmxsize) {
	    imgets(l_rawpath//inpath[1]//in[1]//"[1]","i_naxis2")
	    if (int(imgets.value)!=bpmysize) {
	        printlog("WARNING - NPREPARE: Input images and BPM are not the same size.",
		    l_logfile,verbose+)
	        printlog("                    Not using BPM to generage data quality plane.",
		    l_logfile,verbose+)
	        l_bpm="none"
	    }
        } else {
	    printlog("WARNING - NPREPARE: Input images and BPM are not the same size.",
	        l_logfile,verbose+)
	    printlog("                    Not using BPM to generage data quality planes.",
	        l_logfile,verbose+)
	    l_bpm="none"
        }        
    }

    #--------------------------------------------------------------------------
    # start output
    printlog(" ",l_logfile,l_verbose)
    printlog("  n      input file -->      output file",l_logfile,l_verbose)
    printlog("                 filter     focal plane       input BPM   RON  gain     sat",
        l_logfile,l_verbose)
    printlog(" ",l_logfile,l_verbose)

    l_filter1="" ; l_filter2="" ; l_filter3=""
    l_pupil="none" ; l_fpmask=""

    # The loop:  create VAR and DQ if fl_vardq=yes
    #            determine read noise, gain, saturation, exposure, etc.

    i=1
    while (i<=nimages) {
        # Create tmp FITS file names used within this loop
	tmpsci = mktemp("tmpsci")
	tmpvar = mktemp("tmpvar")
	tmpdq = mktemp("tmpdq")
	tmpinimage = mktemp("tmpinimage")

        in[i]=l_rawpath//inpath[i]//in[i]

        # set read noise, gain, and saturation for the data
        imgets(in[i]//"[0]","NDAVGS", >>& "dev$null")
        if (imgets.value=="0") {
	    printlog("ERROR - NPREPARE: Could not read number of digital averages from header",
	        l_logfile,verbose+)
	    printlog("                  of file "//in[i],l_logfile,verbose+)
	    status=1
	    goto clean
        } else
	    l_ndavgs=real(imgets.value)
        imgets(in[i]//"[0]","COADDS", >>& "dev$null")
        if (imgets.value=="0") {
	    printlog("ERROR - NPREPARE: Could not read number of coadds from header",
	        l_logfile,verbose+)
	    printlog("                  of file "//in[i],l_logfile,verbose+)
	    status=1
	    goto clean
        } else
	    l_coadds=real(imgets.value)
 
        imgets(in[i]//"[0]","LNRS", >>& "dev$null")
        if (imgets.value=="0") {
	    printlog("ERROR - NPREPARE: Could not read number of non-destructive reads from header",
	        l_logfile,verbose+)
	    printlog("                  of file "//in[i],l_logfile,verbose+)
	    status=1
	    goto clean
        } else
	    l_lnrs=real(imgets.value)
        if(l_lnrs==1 && l_ndavgs==1) l_ron=l_ronref*sqrt(l_coadds)
        else if(l_lnrs==1 && l_ndavgs==16) l_ron=l_ronmedref*sqrt(l_coadds)
        else if(l_lnrs==16 && l_ndavgs==16) l_ron=l_ronlowref*sqrt(l_coadds)
        else {
          printlog("WARNING - NPREPARE: Unmatched obseration mode.  Using",
            l_logfile,verbose+)
          printlog("                    medium readnoise mode.",l_logfile,verbose+)
            l_ron=l_ronmedref*sqrt(l_coadds)
        }

        imgets(in[i]//"[0]","A_VDDUC", >>& "dev$null")
        if (imgets.value=="0") {
	    printlog("ERROR - NPREPARE: Could not read A_VDDUC from header.",
	        l_logfile,verbose+)
	    status=1
	    goto clean
        } else
	    l_vdduc=real(imgets.value)

        imgets(in[i]//"[0]","A_VDET", >>& "dev$null")
        if (imgets.value=="0") {
	    printlog("ERROR - NPREPARE: Could not read A_VDET from header.",
	        l_logfile,verbose+)
	    status=1
	    goto clean
        } else
	    l_vdet=real(imgets.value)

        l_biasvolt=l_vdduc-l_vdet
        if (abs(l_biasvolt-l_shallowbias) < 0.05)
	    l_sat=int(l_shallowwell*l_coadds/l_gain)
        else if(abs(l_biasvolt-l_deepbias) < 0.05)
	    l_sat=int(l_deepwell*l_coadds/l_gain)
        else {
	    printlog("ERROR - NPREPARE: Cannot determine saturation level from bias voltage.",
	        l_logfile,verbose+)
	    status=1
	    goto clean
        }
        l_linlimit=int(l_sat*l_linear)

        #-----------------
        # Fix the data if number of non-destructive reads > 1 for data before 
	# Nov 2001
        if ((l_lnrs != 1) && (!alreadyfixed[i])) {
	    oldtimer=no
	    imgets(in[i]//"[0]","DATE-OBS", >>& "dev$null")
	    if (imgets.value=="0" || imgets.value=="" || imgets.value==" ") {
                printlog("WARNING - NPREPARE: Cannot read observation date.  No LNRS scaling.",
		    l_logfile,verbose+)
                printlog("                    (only important for data taken before Nov. 2001)",
		    l_logfile,verbose+)
                oldtimer=no
	    } else {
                if (real(substr(imgets.value,1,4)) < 2001) oldtimer=yes
                if (real(substr(imgets.value,1,4)) > 2001) oldtimer=no
                if ((real(substr(imgets.value,1,4)) == 2001) && \
                    (real(substr(imgets.value,6,7)) < 11)) oldtimer=yes
	    }
	    if (oldtimer) {
                imexpr("a/"//l_lnrs,tmpinimage,in[i]//"[1]",outtype="real",
		    >>& "dev$null")
                printlog("Data taken prior to Nov. 2001:",l_logfile,l_verbose)
                printlog("Dividing "//in[i]//" by number of non-destructive reads: "//l_lnrs,
		    l_logfile,l_verbose)
	    } else
                imcopy(in[i]//"[1]",tmpinimage,verbose-)
        } else
	   imcopy(in[i]//"[1]",tmpinimage,verbose-)

        #-------------
        # fix the WCS zero point for subarray readout for data taken before 2005
        fixwcs = no
        oldtimer = no
        imgets (tmpinimage, "i_naxis2")
        if ((int(imgets.value) != 1024) && (alreadyfixed[i] == no)) {
            imgets (in[i]//"[0]","DATE-OBS", >>& "dev$null")
	    if (imgets.value=="0" || imgets.value=="" || imgets.value==" ") {
                printlog("WARNING - NPREPARE: Cannot read observation date.  \
                    No WCS fixing.", l_logfile, verbose+)
                printlog("          (only important for subarray data taken \
                    before Feb. 2005)", l_logfile, verbose+)
                oldtimer = no
	    } else {
                if (real(substr(imgets.value,1,4)) < 2005) oldtimer = yes
                if (real(substr(imgets.value,1,4)) > 2005) oldtimer = no
                if ((real(substr(imgets.value,1,4)) == 2005) && \
                    (real(substr(imgets.value,6,7)) > 1)) oldtimer = no
	    }
            if (oldtimer) {
                imgets (tmpinimage, "CRPIX1")
                if (imgets.value != "") xcen = real(imgets.value)
                imgets (tmpinimage, "CRPIX2")
                if (imgets.value!="") ycen = real(imgets.value)
                imgets (tmpinimage, "i_naxis1")
                xcen = xcen - (1024-int(imgets.value)) / 2
                imgets (tmpinimage, "i_naxis2")
                ycen = ycen - (1024-int(imgets.value)) / 2
                gemhedit (tmpinimage, "CRPIX1", xcen, comment="", delete-)
                gemhedit (tmpinimage, "CRPIX2", ycen, comment="", delete-)
                fixwcs = yes
                printlog ("Subarray data taken prior to Feb. 2005:",
                    l_logfile,l_verbose)
                printlog ("Shifting WCS zero point in "//in[i]//" : ",
                    l_logfile,l_verbose)
            }
        }

        #-------------
        if (l_fl_vardq) {
	    # create the variance frame, if it doesn't already exist
	    # The variance frame is generated as:
	    #	var = (read noise/gain)^2 + max(data,0.0)/gain

            if (!imaccess(in[i]//"["//l_var_ext//"]")) {
	        imgets(in[i]//"[0]",l_key_gain,>>& "dev$null")
	        if (imgets.value != "0") {
                    l_gain=real(imgets.value)
                    printlog("WARNING - NPREPARE: Gain already set in image header.  Using a gain",
		        l_logfile,verbose+)
                    printlog("                    of "//l_gain//" electrons per ADU.",
		        l_logfile,verbose+)
	        }   

	        imgets(in[i]//"[0]",l_key_ron,>>& "dev$null")
	        if (imgets.value != "0") {
                    l_ron=real(imgets.value)
                    printlog("WARNING - NPREPARE: Read noise already set in image header.",
		        l_logfile,verbose+)
                    printlog("                    Using a read noise of "//l_ron//" electrons.",
		        l_logfile,verbose+)
	        }

	        imgets(in[i]//"[0]",l_key_sat,>>& "dev$null")
	        if (imgets.value != "0") {
                    l_sat=real(imgets.value)
                    printlog("WARNING - NPREPARE: Saturation level already set in image header.",
		        l_logfile,verbose+)
                    printlog("                    Using a saturation level of "//l_sat//" ADU.",
		        l_logfile,verbose+)
	        }
	        l_varexpression="((max(a,0.0))/"//l_gain//"+("//l_ron//"/"//l_gain//")**2)"
	        imexpr(l_varexpression,tmpvar,tmpinimage,outtype="real",verbose-)
            } else {
	        printlog("WARNING - NPREPARE: Variance frame already exists for "//in[i]//".",
		    l_logfile,verbose+) 
	        printlog("                    New variance frame not created.",
		    l_logfile,verbose+)
	        imcopy(in[i]//"["//l_var_ext//"]",tmpvar,verbose-)
            }

	    #-------------
	    # create the DQ frame, if it doesn't already exist
	    # The preliminary DQ frame is constructed by using the bad pixel
	    # mask to set bad pixels to 1, pixels in the non-linear regime to 2,
	    # and saturated pixels to 4.

            if (!imaccess(in[i]//"["//l_dq_ext//"]")) {
	        l_tempexpression="(a>"//l_sat//") ? 4 : ((a>"//l_linlimit//") ? 2 : 0)"
	        imexpr(l_tempexpression,tmpsci,tmpinimage,outtype="short",
		    verbose-)

                # If there's no BPM, then just keep the saturated pixels
	        if (l_bpm=="none")
                    addmasks(tmpsci//".fits",tmpdq//".fits","im1")
	        else
                    addmasks(tmpsci//".fits,"//l_bpm,tmpdq//".fits","im1 || im2")
	        imdelete(tmpsci,ver-)
            } else {
		printlog("WARNING - NPREPARE: Data quality frame already exists for "//in[i]//".",
		    l_logfile,verbose+) 
		printlog("                    New DQ frame not created.",
		    l_logfile,verbose+)
		imcopy(in[i]//"["//l_dq_ext//"]",tmpdq//".fits",verbose-)
            }

	    #-------------
	    # Pack up the results and clean up

            wmef(tmpinimage//","//tmpvar//".fits,"//tmpdq//".fits",out[i],
	        extnames=l_sci_ext//","//l_var_ext//","//l_dq_ext,
		phu=in[i]//".fits[0]",verbose-, >>& "dev$null")
            if (wmef.status !=0) {
	        printlog("ERROR - NPREPARE: Could not write final MEF file (WMEF).",
		    l_logfile,verbose+)
	        status=1
	        goto clean
            }
            gemhedit(out[i]//"["//l_sci_ext//"]","EXTVER",1,"")
            gemhedit(out[i]//"["//l_var_ext//"]","EXTVER",1,"")
            gemhedit(out[i]//"["//l_dq_ext//"]","EXTVER",1,"")
            imdel(tmpvar//","//tmpdq,verify-)

	    # end if(l_fl_vardq) 
        } else {
            wmef(tmpinimage,out[i],extnames=l_sci_ext, phu=in[i]//".fits[0]",
	        verbose-, >>& "dev$null")
            if (wmef.status !=0) {
	        printlog("ERROR - NPREPARE: Could not write final MEF file (WMEF).",
		    l_logfile,verbose+)
	        status=1
	        goto clean
            }
            gemhedit(out[i]//"["//l_sci_ext//"]","EXTVER",1,"")
        }
        imdelete(tmpinimage,verify-)

	#-----------------
	# Fix up the headers
        gemhedit (out[i]//"[0]", "WMEF", "", comment="", delete=yes)
        gemhedit(out[i]//"[0]", "GEM-TLM", "", comment="", delete=yes)

        # Fix the WCS zero point in the PHU, if necessary
        if (fixwcs) {
            imgets (out[i]//"["//l_sci_ext//"]", "CRPIX1")
            xcen = real(imgets.value)
            imgets (out[i]//"["//l_sci_ext//"]", "CRPIX2")
            ycen = real(imgets.value)
            gemhedit (out[i]//"[0]", "CRPIX1", xcen, comment="", delete-)
            gemhedit (out[i]//"[0]", "CRPIX2", ycen, comment="", delete-)
        }

        # Add comment to NEXTEND
        imgets(out[i]//"[0]","NEXTEND", >>& "dev$null")
        if (imgets.value != "0") {
            gemhedit(out[i]//"[0]", "NEXTEND", "", comment="", delete=yes)
            printf("%-8s= %20.0f / %-s\n","NEXTEND",real(imgets.value),
	        "Number of extensions", >> tmphead)
        }

        if (!alreadyfixed[i]) {
            printf("%-8s= %20.2f / %-s\n",l_key_ron,l_ron,
	        "Estimated read noise (electrons)", >> tmphead)
	    printf("%-8s= %20.2f / %-s\n",l_key_gain,l_gain,
	        "Gain (electrons/ADU)", >> tmphead)
	    printf("%-8s= %20.3f / %-s\n","BIASVOLT",l_biasvolt,
	        "Array bias voltage (V)", >> tmphead)
	    printf("%-8s= %20.0f / %-s\n",l_key_sat,l_sat,
	        "Saturation level in ADU", >> tmphead)
	    printf("%-8s= %20.0f / %-s\n",l_key_nonlinear,l_linlimit,
	        "Non-linear regime in ADU", >> tmphead)

            if ((l_coadds != 1) && (!alreadyfixed[i])) {
	        imgets(out[i]//"[0]","EXPTIME", >>& "dev$null")
	        l_exptime = real(imgets.value)
	        printf("%-8s= %20.5f / %-s\n","COADDEXP",l_exptime,
		    "Exposure time (s) for each frame", >> tmphead)
	        l_exptime = l_exptime * l_coadds
                # Can't seem to change this hedit to nhedit without causing a 
                # glibc error in nresidual, so leaving it as hedit - EH
                hedit(out[i]//"[0]","EXPTIME",add-,addonly-,delete+,verify-,
                    show-,update+) 
	        printf("%-8s= %20.5f / %-s\n","EXPTIME",l_exptime,
		    "Exposure time (s) for sum of all coadds", >> tmphead)
	        printlog("Scaling exposure time for "//out[i]//". Coadded exposure time is "//l_exptime//" s.",
		    l_logfile,l_verbose)
            }

            # Parse the individual filter entries into a single header keyword
            imgets(out[i]//"[0]",l_key_filter,>>& "dev$null")
            if (imgets.value == "0") {
	        imgets(out[i]//"[0]","FILTER1",>>& "dev$null")
	        l_filter1=imgets.value
		imgets(out[i]//"[0]","FILTER2",>>& "dev$null")
		l_filter2=imgets.value
		imgets(out[i]//"[0]","FILTER3",>>& "dev$null")
		l_filter3=imgets.value
		if (l_filter1=="0" || l_filter1=="OPEN" || l_filter1=="open" || l_filter1=="Open" || l_filter1=="INVALID")
		    l_filter1=""
		if (l_filter2=="0" || l_filter2=="OPEN" || l_filter2=="open" || l_filter2=="Open" || l_filter2=="INVALID")
		    l_filter2=""
	        if (l_filter3=="0" || l_filter3=="OPEN" || l_filter3=="open" || l_filter3=="Open" || l_filter3=="INVALID")
		    l_filter3=""
	        if (substr(l_filter3,1,3)=="pup") {
	            l_pupil=l_filter3
	            if (substr(l_pupil,strlen(l_pupil)-5,strlen(l_pupil)-4) == "_G")
                        l_pupil=substr(l_pupil,1,(strlen(l_pupil)-6))
	            l_filter3=""
	        }
	        scanfile=tmpfile2

                # strip the Gemini filter numbers before making the combined 
		# filter keyword
	        if (substr(l_filter1,strlen(l_filter1)-5,strlen(l_filter1)-4) == "_G")
                    l_filter1=substr(l_filter1,1,(strlen(l_filter1)-6))
	        if (substr(l_filter2,strlen(l_filter2)-5,strlen(l_filter2)-4) == "_G")
                    l_filter2=substr(l_filter2,1,(strlen(l_filter2)-6))
	        if (substr(l_filter3,strlen(l_filter3)-5,strlen(l_filter3)-4) == "_G")
                    l_filter3=substr(l_filter3,1,(strlen(l_filter3)-6))
	        l_filter=""
	        if ((l_filter1=="blank")||(l_filter2=="blank")||(l_filter3=="blank")) {
	            print("blank", > tmpfile2)
	            l_pupil="none"
	        } else {
	            printf("%s\n%s\n%s\n",l_filter1,l_filter2,l_filter3) | \
	                sort("STDIN",col=0,ignore+,num-,rev-, > tmpfile2)
	        }
	        while (fscan(scanfile,l_temp) != EOF) l_filter=l_filter+l_temp
	        printf("%-8s= \'%-18s\' / %-s\n",l_key_filter,l_filter,
		    "Filter name combined from all 3 wheels", >> tmphead)
	        printf("%-8s= \'%-18s\' / %-s\n","PUPILMSK",l_pupil,
		    "Name of pupil mask", >> tmphead)
		scanfile=""
		delete(tmpfile,ver-, >>& "dev$null")
		delete(tmpfile2,ver-, >>& "dev$null")

		# end if(imgets.value)
            } else {
	        printlog("WARNING - NPREPARE: filter keyword "//l_key_filter//" already exists.",
		    l_logfile,verbose+)
	        printlog("                    NIRI filters not parsed.",
		    l_logfile,verbose+)
	        l_filter=""
            } # end else

	    #-------------
	    # compute pixscale and add to header
	    imgets(out[i]//"[0]","CD1_1",>>& "dev$null")
	    l_cd11=real(imgets.value)
	    imgets(out[i]//"[0]","CD1_2",>>& "dev$null")
	    l_cd12=real(imgets.value)
	    imgets(out[i]//"[0]","CD2_1",>>& "dev$null")
	    l_cd21=real(imgets.value)
	    imgets(out[i]//"[0]","CD2_2",>>& "dev$null")
	    l_cd22=real(imgets.value)
	    l_pixscale=3600.*(sqrt(l_cd11*l_cd11+l_cd12*l_cd12)+sqrt(l_cd21*l_cd21+l_cd22*l_cd22))/2.
	    l_pixscale=real(int(l_pixscale*10000.)/10000.)

            if (l_pixscale != 0.)
	        printf("%-8s= %20.5f / %-s\n","PIXSCALE",l_pixscale,
		    "Pixel scale in arcsec/pixel", >> tmphead)

	    #-------------
	    # Parse the CAMERA and FPMASK keywords to set the appropriate 
	    # SPECSECs
	    l_specsec1="none"
	    l_specsec2="none"
	    l_specsec3="none"

            imgets(out[i]//"[0]","SPECSEC1",>>& "dev$null")
            if (imgets.value!="0") {
	        printlog("WARNING - NPREPARE: Spectral region keyword SPECSEC1 already exists.",
		    l_logfile,verbose+)
	        printlog("                    FPMASK not parsed and SPECSECs not modified.",
		    l_logfile,verbose+)
            } else {
	        imgets(out[i]//"[0]","FPMASK",>>& "dev$null")
	        l_fpmask=imgets.value
	        if (substr(imgets.value,strlen(imgets.value)-5,strlen(imgets.value)-4) == "_G")
	            imgets.value=substr(imgets.value,1,(strlen(imgets.value)-6))
	        if ((imgets.value!="Open") && (imgets.value!="open") && \
                    (imgets.value!="OPEN") && (imgets.value!="f6-cam") && \
                    (imgets.value!="f14-cam") && (imgets.value!="f32-cam") && \
                    (imgets.value!="INVALID")) {
                    # Set DISPAXIS
	            printf("%-8s= %20.0f / %-s\n","DISPAXIS",1,
		        "Dispersion axis", >> tmphead)
	            imgets(out[i]//"[0]","CAMERA",>>& "dev$null")
	            l_camera=imgets.value
	            l_camera=l_camera+l_fpmask

	            if (!access(l_database)) {
                        printlog("WARNING - NPREPARE: Database file "//l_database//" not found.  Using",
			    l_logfile,verbose+)
                        printlog("                    default value of "//l_specsec//" for SPECSEC1.",
			    l_logfile,verbose+)
                        printf("%-8s= \'%-18s\' / %-s\n","SPECSEC1",l_specsec,
			    "Region of image containing spectrum", >> tmphead)
	            } else {
                        fields(l_database,"1-4",lines="1-",quit_if_miss-,
			    print_file_n-) | match(l_camera,stop-) | \
                            scan(l_camera,l_specsec1,l_specsec2,l_specsec3)
                        print(l_specsec1) | scan(l_specsec1)
                        print(l_specsec2) | scan(l_specsec2)
                        print(l_specsec3) | scan(l_specsec3)

                        # Adjust the SPECSECs for smaller ROI images
                        imgets(out[i]//"["//l_sci_ext//"]","i_naxis1",
			    >& "dev$null")
                        if (imgets.value!="0")
        	            xsize=int(imgets.value)
                        imgets(out[i]//"["//l_sci_ext//"]","i_naxis2",
			    >& "dev$null")
                        if (imgets.value!="0")
        	            ysize=int(imgets.value)
                        if (l_specsec1!="none") {
                            sx1=int(substr(l_specsec1,2,
			        stridx(":",l_specsec1)-1))
                            sx2=int(substr(l_specsec1,stridx(":",l_specsec1)+1,
			        stridx(",",l_specsec1)-1))
                            l_temp=substr(l_specsec1,stridx(",",l_specsec1)+1,
			        stridx("]",l_specsec1)-1)
                            sy1=int(substr(l_temp,1,stridx(":",l_temp)-1))
                            sy2=int(substr(l_temp,stridx(":",l_temp)+1,
			        strlen(l_temp)))
                            sx1=sx1-int(512.-xsize/2.)
        		    sx2=sx2-int(512.-xsize/2.)
        		    if (sx1<1) sx1=1
        		    if (sx2>xsize) sx2=xsize
        		    sy1=sy1-int(512.-ysize/2.)
        		    sy2=sy2-int(512.-ysize/2.)
        		    if (sy1<1) sy1=1
        		    if (sy2>xsize) sy2=ysize
        		    l_specsec1="["//sx1//":"//sx2//","//sy1//":"//sy2//"]"
        		}
        		if (l_specsec2!="none") {
        		    sx1=int(substr(l_specsec2,2,
			        stridx(":",l_specsec2)-1))
        		    sx2=int(substr(l_specsec2,stridx(":",l_specsec2)+1,
			        stridx(",",l_specsec2)-1))
        		    l_temp=substr(l_specsec2,stridx(",",l_specsec2)+1,
			        stridx("]",l_specsec2)-1)
        		    sy1=int(substr(l_temp,1,stridx(":",l_temp)-1))
        		    sy2=int(substr(l_temp,stridx(":",l_temp)+1,
			        strlen(l_temp)))
        		    sx1=sx1-int(512.-xsize/2.)
        		    sx2=sx2-int(512.-xsize/2.)
        		    if (sx1<1) sx1=1
        		    if (sx2>xsize) sx2=xsize
        		    sy1=sy1-int(512.-ysize/2.)
        		    sy2=sy2-int(512.-ysize/2.)
        		    if (sy1<1) sy1=1
        		    if (sy2>xsize) sy2=ysize
        		    l_specsec2="["//sx1//":"//sx2//","//sy1//":"//sy2//"]"
        		}
        		if (l_specsec3!="none") {
        		    sx1=int(substr(l_specsec3,2,
			        stridx(":",l_specsec3)-1))
        		    sx2=int(substr(l_specsec3,stridx(":",l_specsec3)+1,
			        stridx(",",l_specsec3)-1))
        		    l_temp=substr(l_specsec3,stridx(",",l_specsec3)+1,
			        stridx("]",l_specsec3)-1)
        		    sy1=int(substr(l_temp,1,stridx(":",l_temp)-1))
        		    sy2=int(substr(l_temp,stridx(":",l_temp)+1,
			        strlen(l_temp)))
        		    sx1=sx1-int(512.-xsize/2.)
        		    sx2=sx2-int(512.-xsize/2.)
        		    if(sx1<1) sx1=1
        		    if(sx2>xsize) sx2=xsize
        		    sy1=sy1-int(512.-ysize/2.)
        		    sy2=sy2-int(512.-ysize/2.)
        		    if(sy1<1) sy1=1
        		    if(sy2>xsize) sy2=ysize
        		    l_specsec3="["//sx1//":"//sx2//","//sy1//":"//sy2//"]"
        		}

                        if ((l_specsec1 !="none")&&(l_specsec1 !=""))
                            printf("%-8s= \'%-18s\' / %-s\n","SPECSEC1",
			        l_specsec1,
				"Region of image containing spectrum #1",
				>> tmphead)
                        if ((l_specsec2 !="none")&&(l_specsec2 !=""))
                            printf("%-8s= \'%-18s\' / %-s\n","SPECSEC2",
			        l_specsec2,
				"Region of image containing spectrum #2",
				>> tmphead)
                        if ((l_specsec3 !="none")&&(l_specsec3 !=""))
                            printf("%-8s= \'%-18s\' / %-s\n","SPECSEC3",
			        l_specsec3,
				"Region of image containing spectrum #3",
				>> tmphead)
                        if ((l_specsec1 == "none") || (l_specsec1 == "")) {
                            printlog("WARNING - NPREPARE: FPMASK does not match a database slit name, or",
			        l_logfile,verbose+)
        		    printlog("                    the combination of CAMERA and FPMASK is not in",
			        l_logfile,verbose+)
        		    printlog("                    the database, so the default "//l_specsec,
			        l_logfile,verbose+)
        		    printlog("                    will be entered in the header.",
			        l_logfile,verbose+)
        		    printf("%-8s= \'%-18s\' / %-s\n","SPECSEC1",
			        l_specsec,"Region of image containing spectrum",
				>> tmphead)

                        }
	            }
	        }
            } # end else

            # Time stamps
	    gemdate ()
	    printf("%-8s= \'%-18s\' / %-s\n","GEM-TLM",gemdate.outdate,
	        "UT Last modification with GEMINI", >> tmphead)
	    printf("%-8s= \'%-18s\' / %-s\n","PREPARE",gemdate.outdate,
	        "UT Time stamp for NPREPARE", >> tmphead)
	    printf("%-8s= \'%-18s\' / %-s\n","BPMFILE",l_bpm,
	        "Input bad pixel mask file", >> tmphead)

            # put all the new stuff in the header
            mkheader(out[i]//"[0]",tmphead,append+,verbose-)
            delete(tmphead,verify-, >>& "dev$null")

            # print output to logfile
            imgets(out[i]//"[0]","FPMASK", >>& "dev$null")
            if (substr(imgets.value,strlen(imgets.value)-5,strlen(imgets.value)-4) == "_G")
	        l_fpmask=substr(imgets.value,1,(strlen(imgets.value)-6))
            else
	        l_fpmask=imgets.value

            printf("%3.0d %15s --> %16s \n",i,in[i],out[i]) | scan(l_struct)
            printlog(l_struct,l_logfile,l_verbose)
            printf("      %17s %15s %15s %5.1f %5.1f %7.0d \n",l_filter,
	        l_fpmask,l_bpm,l_ron,l_gain,l_sat) | scan(l_struct)
            printlog(l_struct,l_logfile,l_verbose)

        } # end if(!alreadyfixed[i])
	
	# Delete tmp files created in this loop
	imdelete(tmpsci//","//tmpvar//","//tmpdq,ver-, >>& "dev$null")
	
        i=i+1
    } # end loop


clean:
    #---------------------------------------------------------------------------
    # Clean up
    if (status==0) {
	printlog(" ",l_logfile,l_verbose)
	printlog("NPREPARE exit status:  good.",l_logfile,l_verbose)
    }
    printlog("----------------------------------------------------------------------------",
        l_logfile,l_verbose)

    scanfile="" 
    delete(tmpfile//","//tmpfile2//","//tmphead,ver-, >>& "dev$null")
    imdelete(tmpsci//","//tmpvar//","//tmpdq,ver-, >>& "dev$null")
    imdelete(tmpinimage,verify-, >>& "dev$null")

end
