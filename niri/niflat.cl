# Copyright(c) 2001-2013 Association of Universities for Research in Astronomy, Inc.

procedure niflat(lampson)

# Derive NIRI bad pixel mask and flat field image
#
# Version  June 19, 2001  JJ, beta-release
#    21 Jun 2001  fixed single-image-without-vardq bug, JJ
#    22 Jun 2001  set saturated pixels, fixed fl_rmstars bugs, JJ
#    27 Jun 2001  fixed minor delete(filelist) bug, JJ
#    29 Jun 2001  added checks for bad imstat output (INDEF), JJ
#     3 Jul 2001  added check for no input files with wildcard, JJ
#    11 Oct 2001  removed BPM from VAR extension header, JJ
#    12 Oct 2001  allow logging from imcomb, some speed improvement 
#                 when fl_vardq-, BM
#     7 Nov 2001  remove logging from imcomb, fix sky flat file list bug, JJ
#    14 Nov 2001  rewrite to write directly to MEF, improve speed,
#                 add default name for flatfile if not specified, JJ
#    19 Nov 2001  replaced gemhedits with one mkheader,
#                 checked for output flat in filelist written to header, JJ
#    03 Jan 2002  fixed call to nisky
#    28 Feb 2002  fixed mask combining bug that trashed the flat bad pix, JJ
#                 also added alreadywarned for various warnings, to avoid repeats
#    Feb 28,2002  JJ v1.3 release
#    24 May 2002  JJ, add output of BPM as .pl file with full header
#    29 May 2002  JJ, minor bug fix on when bpmfile is printed out
#    07 Jun 2002  JJ, fix bug for repeated files in input lists
#                 and bug with output bpm filename
#    10 Jun 2002  JJ, added check for existence of bpm
#    06 Aug 2002  JJ, tiny bug fix in output header format
#    Sept 20, 2002 JJ v1.4 release
#    Aug 18 2003  KL, IRAF2.12 - new/modified parameters
#                       hedit: addonly
#                       imstat: nclip, lsigma, hsigma, cache
#                       imcombine: headers,bpmasks,expmasks,outlimits
#                                  rejmask->rejmasks, plfile->nrejmasks
#    Sep 16, 2004 JJ, changed parameters for nisky
#
char   lampson      {"",prompt="Input images (sky flats lamps on)"}	# OLDP-1-input-primary-combine-suffix=_flat
char   flatfile     {"",prompt="Output flat field image"}		# OLDP-1-output
char   lampsoff     {"",prompt="Input images (lamps off)"}		# OLDP-2-input
char   darks        {"",prompt="Input images (short darks)"}		# OLDP-2-input
char   flattitle    {"default",prompt="Title for output flat image"}	# OLDP-3
char   bpmtitle     {"default",prompt="Title for output bad pixel mask file"}	# OLDP-3
char   bpmfile      {"default",prompt="Name of output bad pixel mask PL file"}  # OLDP-3
char   logfile      {"",prompt="Logfile name"}					# OLDP-1    
real   thresh_flo   {0.80,prompt="Lower bad-pixel threshold (fraction of peak) for flats"}	# OLDP-3
real   thresh_fup   {1.25,prompt="Upper bad-pixel threshold (fraction of peak) for flats"}	# OLDP-3
real   thresh_dlo   {-20.,prompt="Lower bad-pixel threshold (ADU) for darks"}	# OLDP-3
real   thresh_dup   {100.,prompt="Upper bad-pixel threshold (ADU) for darks"}	# OLDP-3
bool   fl_inter     {no,prompt="Set bad pixel cut levels interactively?"}	# OLDP-2
bool   fl_fixbad    {yes,prompt="Fix bad pixels in the output flat field image?"}	# OLDP-2
real   fixvalue     {1.0,prompt="Replace bad pixels in output flat with this value"}	# OLDP-3
char   normstat     {"mean",enum="mean|midpt",prompt="Statistic to use for normalization."}		# OLDP-2
char   combtype     {"default",enum="default|average|median",prompt="Type of combine operation"}	# OLDP-2
char   rejtype      {"avsigclip",enum="none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip",prompt="Type of rejection"}	# OLDP-2  
char   scale        {"none",prompt="Image scaling"}			# OLDP-3
char   zero         {"none",prompt="Image zeropoint offset"}		# OLDP-3
char   weight       {"none",prompt="Image weights"}			# OLDP-3
char   statsec      {"[100:924,100:924]",prompt="Statistics section"}	# OLDP-2
char   key_exptime  {"EXPTIME",prompt="Exposure time header keyword"}	# OLDP-3
real   lthreshold   {INDEF,prompt="Lower threshold"}			# OLDP-3
real   hthreshold   {INDEF,prompt="Upper threshold"}			# OLDP-3
int    nlow         {1,min=0,prompt="minmax: Number of low pixels to reject"}	# OLDP-3
int    nhigh        {1,min=0,prompt="minmax: Number of high pixels to reject"}	# OLDP-3
int    nkeep        {0,min=0,prompt="Minimum to keep or maximum to reject"}  	# OLDP-3
bool   mclip        {yes,prompt="Use median in sigma clipping algorithms?"}	# OLDP-3
real   lsigma       {3.,min=0.,prompt="Lower sigma clipping factor"}		# OLDP-3
real   hsigma       {3.,min=0.,prompt="Upper sigma clipping factor"}		# OLDP-3
char   snoise       {"0.0",prompt="ccdclip: Sensitivity noise (electrons)"}	# OLDP-3
real   sigscale     {0.1,min=0.,prompt="Tolerance for sigma clipping scaling correction"}	# OLDP-3
real   pclip        {-0.5,prompt="pclip: Percentile clipping parameter"}		# OLDP-3
real   grow         {0.0,min=0.,prompt="Radius (pixels) for neighbor rejection"}	# OLDP-3
char   key_ron      {"RDNOISE",prompt="Keyword for readout noise in e-"}		# OLDP-3
char   key_gain     {"GAIN",prompt="Keyword for gain in electrons/ADU"}			# OLDP-3
char   key_sat      {"SATURATI",prompt="Keyword for saturation in ADU"}			# OLDP-3
char   key_nonlinear {"NONLINEA",prompt="Header keyword for non-linear regime (ADU)"}	# OLDP-3
char   sci_ext      {"SCI",prompt="Name or number of science extension"}	# OLDP-3
char   var_ext      {"VAR",prompt="Name or number of variance extension"}	# OLDP-3
char   dq_ext       {"DQ",prompt="Name or number of data quality extension"}	# OLDP-3
bool   fl_rmstars   {no,prompt="Remove stars in flats using NISKY?"}		# OLDP-2
bool   fl_keepmasks {no,prompt="Keep object masks for each input image? (NISKY)"}		# OLDP-2
int    ngrow        {3,prompt="Number of iterations to grow objects into the wings (NISKY)"} # OLDP-2
real   agrow        {3., prompt="Area limit for growing objects into the wings (NISKY)"} # OLDP-2
int    minpix       {6,prompt="Minimum number of pixels to be identified as an object (NISKY)"} # OLDP-3
bool   fl_vardq     {yes,prompt="Create output variance and data quality frames?"}		# OLDP-2
bool   verbose      {yes,prompt="Verbose output?"}			# OLDP-4
int    status       {0,prompt="Exit status (0=good)"}			# OLDP-4
struct *scanfile    {"",prompt="Internal use only"} 			# OLDP-4

begin
        
    char   l_lampson, l_lampsoff, l_darks, l_flatfile, l_bpmfile
    char   l_bpmtitle, l_flattitle, l_statsec, l_key_exptime, l_normstat
    char   l_combtype, l_rejtype, l_scale, l_zero, l_weight
    char   l_key_ron, l_key_gain, l_key_sat
    char   l_snoise, l_logfile, l_title, l_key_nonlinear
    real   l_ron, l_gain, l_sat, l_lthreshold, l_hthreshold, l_lsigma, l_hsigma
    real   l_grow, l_sigscale, l_pclip, dark1, dark2, flat1, flat2
    real   l_lower, l_upper, l_sigthresh, l_nonlinear, l_fixvalue
    int    l_nlow, l_nhigh, l_nkeep, nfiles
    bool   l_fl_dark, l_fl_lampsoff, l_fl_fixbad
    bool   l_fl_rmstars, l_fl_keepmasks, l_fl_vardq
    int    l_minpix, l_ngrow
    real   l_agrow
    bool   l_verbose, l_mclip, l_fl_inter, find, alreadywarned, warning
    char   l_sci_ext, l_var_ext, l_dq_ext
    char   l_sgain, l_sron, lthreshold_char, hthreshold_char
    real   l_thresh_flo, l_thresh_fup, l_thresh_dlo, l_thresh_dup
    char   filelist, scilist, dqlist, infiles[3]
    char   img, suf, tmphead, tmpdq, combout[3]
    char   combvar, combpl, keyfound
    real   norm, stddev, flatnorm, temp
    real   etime[3], rtime
    int    i, k, nbad, len, dum
    struct sdate
    
    bool    intdebug

    intdebug = no

    #------------------------------------------------------------------------
    # Set local parameters
    l_lampson=lampson ; l_lampsoff=lampsoff ; l_darks=darks
    l_flatfile=flatfile ; l_logfile=logfile
    l_bpmtitle=bpmtitle ; l_flattitle=flattitle
    l_combtype=combtype ; l_rejtype=rejtype ; l_normstat=normstat
    l_scale=scale ; l_zero=zero ; l_weight=weight
    l_statsec=statsec ; l_key_exptime=key_exptime ; l_key_sat=key_sat
    l_lthreshold=lthreshold ; l_hthreshold=hthreshold
    l_nlow=nlow ; l_nhigh=nhigh ; l_nkeep=nkeep
    l_lsigma=lsigma ; l_hsigma=hsigma
    l_verbose=verbose ; l_mclip=mclip
    l_key_gain=key_gain ; l_key_ron=key_ron
    l_snoise=snoise ; l_sigscale=sigscale ; l_pclip=pclip ; l_grow=grow
    l_sci_ext = sci_ext ; l_dq_ext = dq_ext ; l_var_ext=var_ext
    l_thresh_flo=thresh_flo ; l_thresh_fup=thresh_fup
    l_thresh_dlo=thresh_dlo ; l_thresh_dup=thresh_dup
    l_fl_inter=fl_inter ; l_fl_rmstars=fl_rmstars
    l_fl_keepmasks=fl_keepmasks
    l_fl_vardq=fl_vardq
    l_key_nonlinear=key_nonlinear ; l_fl_fixbad=fl_fixbad 
    l_fixvalue=fixvalue ; l_bpmfile=bpmfile
    l_ngrow=ngrow ; l_agrow=agrow ; l_minpix=minpix
    alreadywarned=no
    warning=no
    status=0

    # Keep parameters from changing by outside world
    cache ("imgets", "gemdate")

    #------------------------------------------------------------------------
    # Check input images and logic

    if ((l_lampson=="") || (l_lampson==" ")) {
        printlog("ERROR - NIFLAT: No input lamps-on (or sky) images specified.",
            l_logfile,verbose+)
        status=1
        goto clean
    }

    if ((l_lampsoff=="") || (l_lampsoff==" "))
        l_fl_lampsoff=no
    else
        l_fl_lampsoff=yes

    if ((l_darks=="") || (l_darks==" "))
        l_fl_dark=no
    else
        l_fl_dark=yes

    # Check if fl_vardq=no and fl_dark=yes
    if (l_fl_dark && !l_fl_vardq) {
        printlog("WARNING - NIFLAT: Darks not needed if no bad pixel mask is to",
	    l_logfile,verbose+)
        printlog("                  be written (fl_vardq=no).",
	    l_logfile,verbose+)
        l_fl_dark=no
    }

    #------------------------------------------------------------------------
    # Set up log file
    if ((l_logfile=="") || (l_logfile==" ")) {   
	l_logfile=niri.logfile
	if ((l_logfile=="") || (l_logfile==" ")) {  
            l_logfile="niri.log"
            printlog("WARNING - NIFLAT: Both niflat.logfile and niri.logfile are",
	        l_logfile,verbose+) 
            printlog("                  undefined.  Using niri.log.",
	        l_logfile,verbose+)
	}
    }

    date | scan(sdate)
    printlog("----------------------------------------------------------------------------",
        l_logfile,l_verbose)
    printlog("NIFLAT -- "//sdate,l_logfile,l_verbose)
    printlog(" ",l_logfile,l_verbose)

    #------------------------------------------------------------------------
    # Define temporary files
    filelist=mktemp("tmpfilelist")
    infiles[1]=l_lampson 
    if (l_fl_lampsoff) {
        infiles[2]=l_lampsoff 
        combout[2]=mktemp("tmpoff")
    } else
        infiles[2]=""

    if (l_fl_dark) {  
        infiles[3]=l_darks
        combout[3]=mktemp("tmpdark")
    } else
        infiles[3]=""
    scilist=mktemp("tmpscilist")
    combout[1]=mktemp("tmpon")
    combpl=mktemp("tmppl")
    tmpdq=mktemp("tmpdq")
    tmphead=mktemp("tmphead")


    #------------------------------------------------------------------------
    # Verify that the SCI extension name is not empty, otherwise exit gracefully
    if (l_sci_ext=="" || l_sci_ext==" ") {
        printlog("ERROR - NIFLAT: Extension name sci_ext is missing",
	    l_logfile,verbose+)
        status=1
        goto clean
    } 

    # Set titles
    if (l_flattitle=="default" || l_flattitle=="" || l_flattitle==" ") 
        l_flattitle="FLAT FIELD IMAGE from gemini.niri.niflat"
    if (l_bpmtitle=="default" || l_bpmtitle=="" || l_bpmtitle==" ") 
        l_bpmtitle="BAD PIXEL MASK from gemini.niri.niflat"

    # Check if key_exptime is set
    if (l_key_exptime=="" || l_key_exptime==" ") {
        printlog("WARNING - NIFLAT: key_exptime not set, cannot check exposure times.",
            l_logfile,verbose+)
    }

    # Check threshold values if fl_inter=no
    if (!l_fl_inter) {
        if ((l_thresh_flo > l_thresh_fup) || (l_thresh_dlo > l_thresh_dup)) {
	    printlog("ERROR - NIFLAT: Lower threshold must be lower than the upper threshold.",
	        l_logfile,verbose+)
            status=1
            goto clean
        }
    }

    # Check if fl_fixbad=yes and fl_vardq=no
    if (l_fl_fixbad && !l_fl_vardq) {
        printlog("WARNING - NIFLAT: Cannot fix bad pixels without bad pixel mask, but",
	    l_logfile,verbose+)
        printlog("                  fl_vardq=no.  Setting fl_fixbad=no.",
	    l_logfile,verbose+)
        l_fl_fixbad=no
    }

    #--------------------------------------------------------------------------
    # check that files exist and are MEF
    nbad=0
    for (i = 1; i<=3; i+=1) {
        if ((i==1) || ((i==2) && (l_fl_lampsoff)) || ((i==3) && (l_fl_dark))) {
	    if (infiles[i] == "" || infiles[i]==" ") {
	        printlog("ERROR - NIFLAT: input files not specified",
		    l_logfile,verbose+)
	        status=1
	        goto clean
	    } 

            # check existence of list files
	    if (substr(infiles[i],1,1) == "@") {
	        len=strlen(infiles[i])
	        if (!access(substr(infiles[i],2,len))) {
                    printlog("ERROR - NIFLAT: "//substr(infiles[i],2,len)//" does not exist.",
		        l_logfile,verbose+)
                    status=1
                    goto clean
	        }
	    }
	    files(infiles[i],sort-) | unique("STDIN", > filelist)

            # Check that all images in the list exist
	    scanfile=filelist
	    nfiles=0
	    while (fscan(scanfile, img) != EOF) {
	        gimverify(img)
	        if (gimverify.status==1) {
                    printlog("ERROR - NIFLAT: File "//img//" not found.",
		        l_logfile,verbose+)
                    nbad+=1
	        } else if (gimverify.status>1) {
                    printlog("ERROR - NIFLAT: File "//img//" not a MEF FITS image.",
		        l_logfile,verbose+)
                    nbad+=1
	        } else {
		    keyfound=""
		    hselect(img//"[0]","*PREPAR*",yes) | scan(keyfound)
		    if (keyfound == "") {
                        printlog("ERROR - NIFLAT: Image "//img//" not *PREPAREd.",
			    l_logfile,verbose+)
                        nbad+=1
                    }
	        }
	        nfiles+=1
	    } #end of while loop

	    # check for empty file list
	    if (nfiles==0) {
	        printlog("ERROR - NIFLAT:  No input images meet wildcard criteria.",
		    l_logfile,verbose+)
	        status=1
	        goto clean
	    }

	    scanfile=""
	    delete(filelist,verify-, >>& "dev$null")
        } #end if(i==1) && etc.
    }  #end of for loop

    # Exit if problems found
    if (nbad > 0) {
        printlog("ERROR - NIFLAT: "//nbad//" image(s) either do not exist, are not MEF files, or",
	    l_logfile,verbose+)
        printlog("                 have not been run through *PREPARE.",
	    l_logfile,verbose+)
	status=1
	goto clean
    }

    #---------------------------------------------------------------------------
    # Input info to logfile
    printlog("lampson = "//l_lampson,l_logfile,l_verbose)
    if (l_fl_lampsoff) printlog("lampsoff = "//l_lampsoff,l_logfile,l_verbose)
    if (l_fl_dark) printlog("darks = "//l_darks,l_logfile,l_verbose)
    if (l_flatfile!="" && l_flatfile!=" ")
        printlog("flatfile = "//l_flatfile,l_logfile,l_verbose)
    if (!l_fl_inter) {
        printlog("thresh_flo = "//l_thresh_flo,l_logfile,l_verbose)
        printlog("thresh_fup = "//l_thresh_fup,l_logfile,l_verbose)
        if (l_fl_dark)
	    printlog("thresh_dlo = "//l_thresh_dlo,l_logfile,l_verbose)
        if (l_fl_dark)
	    printlog("thresh_dup = "//l_thresh_dup,l_logfile,l_verbose)
    }
    printlog("normstat = "//l_normstat,l_logfile,l_verbose)
    if (l_fl_rmstars) {
        printlog("Removing stars from input images.",l_logfile,l_verbose)
        printlog("  threshold for object identification = 3 sigma",l_logfile,l_verbose)
        printlog("  minimum pixels for object identification = "//l_minpix,l_logfile,l_verbose)
        printlog("  iterations for growing objects into the wings = "//l_ngrow,l_logfile,l_verbose)
        printlog("  maximum area for growing objects into the wings = "//l_agrow,l_logfile,l_verbose)
    }

    #---------------------------------------------------------------------------
    # Combine, lamps on, lamps off, then darks
    for (i=1; i<=3; i+=1) {
	# Create tmp FITS file names re-used within this loop
	combvar=mktemp("tmpvar")

        if (infiles[i] != "" && infiles[i] != " ") {
	    files(infiles[i],sort-) | unique("STDIN", > filelist)
	    k=0
	    scanfile=filelist
	    # strip out science and dq planes
	    while (fscan(scanfile, img) != EOF) {
	        k=k+1
	        suf = substr(img,strlen(img)-4,strlen(img))
	        # Add the .fits suffix to the input file:
	        if (suf!=".fits" && imaccess(img//".fits"))
                    img = img//".fits"
	        if ((i==1) && (k==1)) {
                    if ((l_flatfile=="") || (l_flatfile==" ")) {
                        suf = substr(img,strlen(img)-4,strlen(img))
                        if (suf==".fits")
        	            l_flatfile=substr(img,1,strlen(img)-5)//"_flat" 
                        else
        	            l_flatfile=img//"_flat"
                        printlog("flatfile = "//l_flatfile,l_logfile,l_verbose)
                        if (l_fl_vardq) {
        	            if (l_bpmfile=="default" || l_bpmfile=="" || l_bpmfile==" ") {
        	                l_bpmfile=substr(l_flatfile,1,
				    strlen(l_flatfile)-4)+"bpm.pl"
        	                printlog("bpmfile = "//l_bpmfile,
				    l_logfile,l_verbose)
        	            }
                        }
                    } else {
                        suf = substr(l_flatfile,strlen(l_flatfile)-4,
			    strlen(l_flatfile))
                        if (suf==".fits") {
        	            l_flatfile=substr(l_flatfile,1,
			        strlen(l_flatfile)-5)
                        }
                        if (l_bpmfile=="default" || l_bpmfile=="" || l_bpmfile==" ") {
        	            l_bpmfile=l_flatfile+"_bpm.pl"
        	            printlog("bpmfile = "//l_bpmfile,l_logfile,l_verbose)
                        }
                    }
                    if (imaccess(l_flatfile)) {
                        printlog("ERROR - NIFLAT: Output file "//l_flatfile//" already exists.",
			    l_logfile,verbose+)
                        status=2
                        goto clean
                    }
                    if (fl_vardq && imaccess(l_bpmfile)) {
                        printlog("ERROR - NIFLAT: Output file "//l_bpmfile//" already exists.",
			    l_logfile,verbose+)
                        status=2
                        goto clean
                    }
	            imcopy(img//"[0]",l_flatfile,verbose-)
	        }

                # check gain, readnoise, and saturation values
	        imgets(img//"[0]",l_key_gain, >>& "dev$null")
	        if (imgets.value == "0" || imgets.value=="") {
                    printlog("ERROR - NIFLAT: Could not get gain from header of image "//img,
		        l_logfile,verbose+)
                    status=1
                    goto clean
	        } else if ((k==1) && (i==1)) {
                    l_gain=real(imgets.value)
                    l_sgain=imgets.value
	        } else if ((i!=3)&&(!alreadywarned)) {
                    if (abs(real(imgets.value)-l_gain) > 0.5) {
                        printlog("WARNING - NIFLAT: gain values are different.  Continuing, but the gain in",
			    l_logfile,verbose+)
                        printlog("                  the output header will be wrong.",
			    l_logfile,verbose+)
                        warning=yes
                    }
	        }

	        imgets(img//"[0]",l_key_ron, >>& "dev$null")
	        if (imgets.value == "0" || imgets.value=="") {
                    printlog("ERROR - NIFLAT: Could not get read noise from header of image "//img,l_logfile,verbose+)
                    status=1
                    goto clean
	        } else if ((k==1) && (i==1)) {
                    l_ron=real(imgets.value)
                    l_sron=imgets.value
	        } else if ((i!=3)&&(!alreadywarned)) {
                    if (abs(real(imgets.value)-l_ron) > 0.5) {
                        printlog("WARNING - NIFLAT: read noise values are different.  Continuing, but the",
			    l_logfile,verbose+)
                        printlog("                  read noise in the output header will be wrong.",
			    l_logfile,verbose+)
                        warning=yes
                    }
	        }

	        imgets(img//"[0]",l_key_sat, >>& "dev$null")
	        if (imgets.value == "0" || imgets.value=="") {
                    printlog("ERROR - NIFLAT: Could not get saturation level from header of image "//img,
		        l_logfile,verbose+)
                    status=1
                    goto clean
	        } else if ((k==1) && (i==1)) {
                    l_sat=real(imgets.value)
	        } else if ((i!=3)&&(!alreadywarned)) {
                    if (abs(real(imgets.value)-l_sat) > 1.) {
                        printlog("WARNING - NIFLAT: saturation values are different.  Continuing, but the",
			    l_logfile,verbose+)
                        printlog("                  saturation in the output header will be wrong.",
			    l_logfile,verbose+)
                        warning=yes
                    }
	        }

	        imgets(img//"[0]",l_key_nonlinear, >>& "dev$null")
	        if (imgets.value == "0" || imgets.value=="") {
                    printlog("ERROR - NIFLAT: Could not get non-linear level from header of image "//img,
		        l_logfile,verbose+)
                    status=1
                    goto clean
	        } else if ((k==1) && (i==1)) {
                    l_nonlinear=real(imgets.value)
	        } else if ((i!=3)&&(!alreadywarned)) {
                    if (abs(real(imgets.value)-l_nonlinear) > 1.) {
                        printlog("WARNING - NIFLAT: nonlinear values are different.  Continuing, but the",
			    l_logfile,verbose+)
                        printlog("                  nonlinear level in the output header will be wrong.",
			    l_logfile,verbose+)
                        warning=yes
                    }      
	        }
	        if (warning) alreadywarned=yes

	        # science extension
	        print(img//"["//l_sci_ext//"]", >> scilist)
	        # DQ extension - set up for imcombine
	        if (imaccess(img//"["//l_dq_ext//"]")) {
                    imcopy(img//"["//l_dq_ext//"]",tmpdq//"_"//k//".pl",verbose-)
                    gemhedit(img//"["//l_sci_ext//"]", "BPM", tmpdq//"_"//k//".pl", 
                        "Bad Pixel Mask")
	        }
	    } # end while(fscan) 
	    scanfile=""

            # combine images----------------------------

	    #special cases for low numbers of images
	    if (k == 1) {
	        printlog("WARNING - NIFLAT: only one input image.",
		    l_logfile,verbose+)
	        imcopy(img//"[0]",combout[i],verbose-)
	        imcopy(img//"["//l_sci_ext//"]",
		    combout[i]//"["//l_sci_ext//",append]",verbose-)
	        if (imaccess(img//"["//l_dq_ext//"]") && imaccess(img//"["//l_var_ext//"]"))  {
                    imcopy(img//"["//l_dq_ext//"]",
		        combout[i]//"["//l_dq_ext//",append]",verbose-)
                    imcopy(img//"["//l_var_ext//"]",
		        combout[i]//"["//l_var_ext//",append]",verbose-)        
	        } else if (l_fl_vardq) {
                    l_fl_vardq=no
                    l_fl_fixbad=no
                    l_fl_dark=no
                    printlog("WARNING - NIFLAT: only one input image without VAR and DQ planes.",
		        l_logfile,verbose+)
                    printlog("                  Proceeding with fl_vardq=no and fl_fixbad=no",
		        l_logfile,verbose+)
	        }
	    } else {
	        if (k <= 4) {
                    if (k == 2) {
                        printlog("WARNING - NIFLAT: only two images to combine, turning off rejection.",
        	            l_logfile,verbose+)
                        l_rejtype="none"
                    } else {
                        printlog("WARNING - NIFLAT: four or less images to combine.",
			    l_logfile,verbose+)
                    }
	        }        
	        if (i==1 && l_fl_rmstars) {
                    printlog(" ",l_logfile,l_verbose)
                    printlog("NIFLAT calling NISKY",l_logfile,l_verbose)
                    nisky("@"//filelist,outimage=combout[i],outtitle="default",
                        combtype=l_combtype, rejtype=l_rejtype,logfile=l_logfile,
                        nlow=0,nhigh=1,lsigma=l_lsigma,hsigma=l_hsigma,
			threshold=3.,key_ron=l_key_ron,
			key_gain=l_key_gain,
			masksuffix="msk",ngrow=l_ngrow,agrow=l_agrow,minpix=l_minpix,
			fl_keepmasks=l_fl_keepmasks,
                        sci_ext=l_sci_ext,var_ext=l_var_ext,dq_ext=l_dq_ext,
                        fl_vardq=l_fl_vardq,fl_dqprop=no,verbose=l_verbose)

                    if (nisky.status !=0) {
                        printlog("ERROR - NIFLAT: Problem in NISKY prevents determination of sky flat.",
			    l_logfile,verbose+)
                        status=1
                        goto clean
                    }
	        } else {
                    if (l_combtype=="default") l_combtype="average"
                    if (l_fl_vardq==no) {
        	        combpl=""
        	        combvar=""
                    }
                    imcopy(l_flatfile,combout[i],verbose-)
                    imcombine("@"//scilist,combout[i]//"["//l_sci_ext//",append]",
	                headers="",bpmasks="",rejmasks="",nrejmasks=combpl,
			expmasks="",sigmas=combvar,logfile="STDOUT",
			combine=l_combtype,reject=l_rejtype,project=no,
			outtype="real",outlimits="",offsets="none",
			masktype="goodvalue",maskvalue=0,blank=0,scale=l_scale,
			zero=l_zero,weight=l_weight,statsec=l_statsec,
			expname=l_key_exptime,lthreshold=l_lthreshold,
			hthreshold=l_hthreshold,nlow=l_nlow,nhigh=l_nhigh,
			nkeep=l_nkeep,mclip=l_mclip,lsigma=l_lsigma,
	                hsigma=l_hsigma,rdnoise=l_sron,gain=l_sgain,
			snoise=l_snoise,sigscale=l_sigscale,pclip=l_pclip,
			grow=l_grow,>>& "dev$null")

                    # check the imcombine produced something
                    if (!imaccess(combout[i])) {
	               status=1
	               printlog("NIFLAT - ERROR: imcombine error",l_logfile,yes)
                       goto clean
                    }

                    if (l_fl_vardq) {
                        # make variance image by squaring combvar
                        imarith(combvar,"*",combvar,
			    combout[i]//"["//l_var_ext//",append]",
			    pixtype="real",verbose-)
                        # bad pixels have dqsum=n
                        imexpr("(a=="//k//") ? 1 : 0",
			    combout[i]//"["//l_dq_ext//",append]",combpl,
			    outtype="short",verbose-)
                    }
	        } # end else imcombine
	    } # end else 

        scanfile=scilist
        while (fscan(scanfile, img) != EOF) {
            gemhedit (img, "BPM", "", "", delete=yes)
        }
        scanfile=""

#            gemhedit ("@"//scilist, "BPM", "", "", delete=yes)
	    delete(filelist,verify-, >>& "dev$null")
	    delete(scilist,verify-, >>& "dev$null")
	    imdelete(tmpdq//"*.pl",verify-, >>& "dev$null")

            gemhedit (combout[i]//"[0]", "BPM", "", "", delete=yes)
            gemhedit (combout[i]//"["//l_sci_ext//"]", "BPM", "", "", delete=yes)
            gemhedit (combout[i]//"["//l_var_ext//"]", "EXTVER", "1",
                "Extension version")
        } # end if(infiles[i]!="" etc)

	imdelete(combvar//","//combpl,verify-, >>& "dev$null")
    } # end of for loop

    #---------------------------------------------------------------------------
    # subtract and normalize flats, combine bad pixel masks

    # darks-------------------------

    if (l_fl_dark) {
        if (l_fl_inter) {
	    imhistogram(combout[3]//"["//l_sci_ext//"]",z1=-100,z2=1000,
	        binwidth=INDEF,nbins=512,autoscale+,top_closed-,
		hist_type="normal",listout-,plot_type="line",logy+,
		device="stdgraph")
	    print("Dark histogram")
	    print("Hit any key to mark lower limit for good pixels")
	    dum=fscan(gcur,l_lower)
	    print("Hit any key to mark upper limit for good pixels")
	    dum=fscan(gcur, l_upper)
	    print(" ")
	    if (l_lower>l_upper) {
	        temp=l_upper
	        l_upper=l_lower
	        l_lower=temp
	    }
        } else {
	    l_lower=l_thresh_dlo
	    l_upper=l_thresh_dup
        }
        imexpr("((b<"//l_lower//") || (b>"//l_upper//") || (a>0)) ? 1 : 0",
	    combout[3]//"["//l_dq_ext//",overwrite]",
	    combout[3]//"["//l_dq_ext//"]",combout[3]//"["//l_sci_ext//"]",
	    outtype="short",verbose-)
        dark1=l_lower
        dark2=l_upper
    } # end if(l_fl_dark)

    # flat--------------------------

    if (l_fl_lampsoff) {
        imstat(combout[1]//"["//l_sci_ext//"]"//l_statsec,fields="midpt,stddev",
	    lower=100.,upper=l_sat,nclip=0,lsigma=INDEF,usigma=INDEF,
	    binwidth=0.1,format-,cache-) | scan(temp,stddev)
        imstat(combout[1]//"["//l_sci_ext//"]"//l_statsec,fields="midpt",
	    lower=(temp-3*stddev),upper=(temp+3*stddev),nclip=0,lsigma=INDEF,
	    usigma=INDEF,binwidth=0.1,format-,cache-) | scan(temp)
        imstat(combout[2]//"["//l_sci_ext//"]"//l_statsec,fields="midpt,stddev",
	    lower=100.,upper=l_sat,nclip=0,lsigma=INDEF,usigma=INDEF,
	    binwidth=0.1,format-,cache-) | scan(norm,stddev)
        imstat(combout[2]//"["//l_sci_ext//"]"//l_statsec,fields="midpt",
	    lower=(norm-3*stddev),upper=(norm+3*stddev),nclip=0,lsigma=INDEF,
	    usigma=INDEF,binwidth=0.1,format-,cache-) | scan(norm)
        if (norm==INDEF || temp==INDEF) {
	    printlog("ERROR - NIFLAT: Statistics failed, possibly due to a bad statsec.",l_logfile,verbose+)
	    status=1
	    goto clean
        }

        if (temp<norm) {
	    printlog("WARNING - NIFLAT:  Flux level in lampson is lower than in lampsoff.",
                l_logfile,verbose+)
	    printlog("                   Proceeding assuming these should be reversed.",
                l_logfile,verbose+)
	    imarith(combout[2]//"["//l_sci_ext//"]","-",
	        combout[1]//"["//l_sci_ext//"]",
		combout[1]//"["//l_sci_ext//",overwrite]",pixtype="real",
		verbose-)
        } else {
	    imarith(combout[1]//"["//l_sci_ext//"]","-",
	        combout[2]//"["//l_sci_ext//"]",
		combout[1]//"["//l_sci_ext//",overwrite]",pixtype="real",
		verbose-)
        }
        if (l_fl_vardq) {
	    imarith(combout[1]//"["//l_var_ext//"]","+",
	        combout[2]//"["//l_var_ext//"]",
		combout[1]//"["//l_var_ext//",overwrite]",pixtype="real",
		verbose-)
	    addmasks(combout[1]//"["//l_dq_ext//"],"//combout[2]//"["//l_dq_ext//"]",
	        combout[1]//"["//l_dq_ext//",overwrite]","im1 || im2",flags=" ")
        }
        #adjust saturation value
        l_sat=l_sat-norm
        l_nonlinear=l_nonlinear-norm
    }

    # just use combout[1] for now, not dsub
    # Get normalization and good pixel limits for flat
    if (l_fl_inter) {
        imhistogram(combout[1]//"["//l_sci_ext//"]",z1=INDEF,z2=INDEF,
	    binwidth=INDEF,nbins=512,autoscale+,top_closed-,hist_type="normal",
	    listout-,plot_type="line",logy+,device="stdgraph")
        print("Flat histogram")
        print("Hit any key to mark lower limit for good pixels")   
        dum=fscan(gcur,l_lower)
        print("Hit any key to mark upper limit for good pixels")
        dum=fscan(gcur, l_upper)
	print(" ")
	if (l_lower>l_upper) {
	    temp=l_upper
	    l_upper=l_lower
	    l_lower=temp
	}
        imstat(combout[1]//"["//l_sci_ext//"]"//l_statsec,fields=l_normstat,
	    lower=l_lower,upper=l_upper,nclip=0,lsigma=INDEF,usigma=INDEF,
	    binwidth=0.1,format-,cache-) | scan(norm)
    } else {
        imstat(combout[1]//"["//l_sci_ext//"]"//l_statsec,fields="mode,stddev",
	    lower=100.,upper=l_sat,nclip=0,lsigma=INDEF,usigma=INDEF,
	    binwidth=0.1,format-,cache-) | scan(norm,stddev)
        imstat(combout[1]//"["//l_sci_ext//"]"//l_statsec,fields="midpt",
	    lower=(norm-stddev*3),upper=(norm+stddev*3),nclip=0,lsigma=INDEF,
	    usigma=INDEF,binwidth=0.1,format-,cache-) | scan(norm)
        if (norm==INDEF) {
	    printlog("ERROR - NIFLAT: Statistics failed, possibly due to a bad statsec.",
	        l_logfile,verbose+)
	    status=1
	    goto clean
        }

        l_lower=norm*l_thresh_flo
        l_upper=norm*l_thresh_fup
        if (l_upper > l_sat) {
	    printlog("WARNING - NIFLAT: Upper threshold for good pixels in the flat field image",
	        l_logfile,verbose+)
	    printlog("                  is larger than the saturation value.  Resetting upper",
	        l_logfile,verbose+)
	    printlog("                  limit to the saturation value.",
	        l_logfile,verbose+)
	    l_upper=l_sat
        }
        imstat(combout[1]//"["//l_sci_ext//"]"//l_statsec,fields=l_normstat,
	    lower=l_lower,upper=l_upper,nclip=0,lsigma=INDEF,usigma=INDEF,
	    binwidth=0.1,format-,cache-) | scan(norm)
        if (norm==INDEF) {
	    printlog("ERROR - NIFLAT: Statistics failed, possibly due to a bad statsec.",
	        l_logfile,verbose+)
	    status=1
	    goto clean
        }
    }
    flat1=l_lower
    flat2=l_upper

    # combine masks------------------
    if (l_fl_vardq) {
        imexpr("(b>"//l_sat//") ? 4 : (b>"//l_nonlinear//") ? 2 : ((b<"//l_lower//") || (b>"//l_upper//") || (a>0)) ? 1 : 0",
	    combout[1]//"["//l_dq_ext//",overwrite]",
	    combout[1]//"["//l_dq_ext//"]",combout[1]//"["//l_sci_ext//"]",
	    outtype="short",verbose-)
        if (l_fl_dark) {
	    addmasks(combout[1]//"["//l_dq_ext//"],"//combout[3]//"["//l_dq_ext//"]",
	        combout[1]//"["//l_dq_ext//",overwrite]","im1 || im2",flags=" ")
        }
    }

    # normalize flat-----------------
    printlog(" ",l_logfile,l_verbose)
    printlog("Normalizing flat by "//norm,l_logfile,l_verbose)
    imarith(combout[1]//"["//l_sci_ext//"]","/",norm,
        combout[1]//"["//l_sci_ext//",overwrite]",divzero=0.,pixtype="real",
	verbose-)
    imcopy(combout[1]//"["//l_sci_ext//"]",
        l_flatfile//"["//l_sci_ext//",append]",verbose-)
    if (l_fl_vardq) {
        imarith(combout[1]//"["//l_var_ext//"]","/",(norm*norm),
	    combout[1]//"["//l_var_ext//",overwrite]",divzero=0.,pixtype="real",
	    verbose-)
        imcopy(combout[1]//"["//l_var_ext//"]",
	    l_flatfile//"["//l_var_ext//",append]",verbose-)
        imcopy(combout[1]//"["//l_dq_ext//"]",
	    l_flatfile//"["//l_dq_ext//",append]",verbose-)
    }

    # fix bad pixels in flat image
    if (l_fl_fixbad) {
        imexpr("(b==1) ? "//l_fixvalue//" : a",
	    l_flatfile//"["//l_sci_ext//",overwrite]",
	    l_flatfile//"["//l_sci_ext//"]",l_flatfile//"["//l_dq_ext//"]",
	    outtype="real",verbose-)
    }

    #--------------------------------------------------------------------------
    # update headers in final image PHU
    
	print('default_pars after="", before="", add+, addonly-, del-, ver-, \
        show-, update+\n', >> tmphead)
    
        gemhedit (l_flatfile//"[0]", "i_title", l_flattitle, "")
        gemhedit (l_flatfile//"["//l_sci_ext//"]", "i_title", l_flattitle, "")
    if (l_fl_vardq) {
        gemhedit(l_flatfile//"["//l_var_ext//"]", "i_title", l_flattitle, "")
        gemhedit(l_flatfile//"["//l_dq_ext//"]", "i_title", l_bpmtitle, "")
        gemhedit(l_flatfile//"["//l_var_ext//"]", "BPM", "", "", delete=yes)
    }
    l_flatfile=l_flatfile//"[0]"
    gemdate ()
    gemhedit(l_flatfile,"GEM-TLM",gemdate.outdate,"UT Last modification with GEMINI")
    printf("%-8s  \'%-18s\'   \'%-s\'\n","NIFLAT",gemdate.outdate,
		"UT Time stamp for NIFLAT", >> tmphead) 
    printf("%-8s  %20.5f   \'%-s\'\n","NIFLTNRM",norm,
        "Normalization constant for NIFLAT", >> tmphead)
    printf("%-8s  \'%-18s\'   \'%-s\'\n","NIFLTNST",l_normstat,
        "Normalization statistic for NIFLAT", >> tmphead)
    printf("%-8s  %20.5f   \'%-s\'\n","NIFLTFLO",flat1,
        "Lower limit for good pixels in flat (ADU)", >> tmphead)
    printf("%-8s  %20.5f   \'%-s\'\n","NIFLTFUP",flat2,
        "Upper limit for good pixels in flat (ADU)", >> tmphead)

    if (l_fl_dark) {
	dark1=(int(dark1*100.))/100.
	dark2=(int(dark2*100.))/100.
	printf("%-8s  %20.5f   \'%-s\'\n","NIFLTDLO",dark1,
	    "Lower limit for good pixels in dark (ADU)", >> tmphead)
	printf("%-8s  %20.5f   \'%-s\'\n","NIFLTDUP",dark2,
	    "Upper limit for good pixels in dark (ADU)", >> tmphead)
    }
    if (l_fl_fixbad)
	printf("%-8s  %20.5f   \'%-s\'\n","NIFLTFIX",l_fixvalue,
	    "Bad pixels replaced with this value", >> tmphead)
    if (l_fl_rmstars) {
	printf("%-8s  \'%-18b\'   \'%-s\'\n","NIFLTRMS",l_fl_rmstars,
	    "Stars removed by NISKY", >> tmphead)
         gemhedit(l_flatfile, "NISKYMSK", "", "", delete=yes)
    } else {
	printf("%-8s  \'%-18s\'   \'%-s\'\n","NIFLTCOM",l_combtype,"Combine method used by NIFLAT", >> tmphead)
	if (l_zero!="none")
	    printf("%-8s  \'%-18s\'   \'%-s\'\n","NIFLTZER",l_zero,
	        "Statistic used to compute additive offsets", >> tmphead)
	if (l_scale!="none")
	    printf("%-8s  \'%-18s\'   \'%-s\'\n","NIFLTSCL",l_scale,
	        "Statistic used to compute scale factors", >> tmphead)
	if (l_weight!="none")
	    printf("%-8s  \'%-18s\'   \'%-s\'\n","NIFLTWEI",l_weight,
	        "Statistic used to compute relative weights", >> tmphead)
	printf("%-8s  \'%-18s\'   \'%-s\'\n","NIFLTSTA",l_statsec,
	    "Statistics section used by NIFLAT", >> tmphead)
	if (isindef(l_lthreshold))  lthreshold_char="+INDEF"
	else  lthreshold_char=str(l_lthreshold)	
    if (isindef(l_hthreshold))  hthreshold_char="+INDEF"
    else  hthreshold_char=str(l_hthreshold)
	if (l_rejtype!="none") {
	    printf("%-8s \'%-18s\'  \'%-s\'\n","NIFLTREJ",l_rejtype,
	        "Rejection algorithm used by NIFLAT", >> tmphead)
            printf("%-8s \'%-18s\'   \'%-s\'\n","NIFLTLTH",l_lthreshold,
	        "Lower threshold before combining", >> tmphead)
            printf("%-8s \'%-18s\'   \'%-s\'\n","NIFLTHTH",l_hthreshold,
	        "Upper threshold before combining", >> tmphead)
            printf("%-8s %20.5f   \'%-s\'\n","NIFLTGRW",l_grow,
	        "Radius for additional pixel rejection", >> tmphead)
	}
	if (l_rejtype=="minmax") {
            printf("%-8s %20.0f   \'%-s\'\n","NIFLTNLO",l_nlow,
	        "Low pixels rejected (minmax)", >> tmphead)
            printf("%-8s %20.0f   \'%-s\'\n","NIFLTNHI",l_nhigh,
	        "High pixels rejected (minmax)", >> tmphead)
	}
	if (l_rejtype=="sigclip" || l_rejtype=="avsigclip" || l_rejtype=="ccdclip" \
              || l_rejtype=="crreject" || l_rejtype=="pclip") {
	    printf("%-8s %20.5f   \'%-s\'\n","NIFLTLSI",l_lsigma,
	        "Lower sigma for rejection", >> tmphead)
            printf("%-8s %20.5f   \'%-s\'\n","NIFLTHSI",l_hsigma,
	        "Upper sigma for rejection", >> tmphead)
            printf("%-8s %20.0f   \'%-s\'\n","NIFLTNKE",l_nkeep,
	        "Min(max) number of pixels to keep(reject)", >> tmphead)
            printf("%-8s \'%-18b\'   \'%-s\'\n","NIFLTMCL",l_mclip,
	        "Use median in clipping algorithms", >> tmphead)
	} 
	if (l_rejtype=="sigclip" || l_rejtype=="avsigclip" || l_rejtype=="ccdclip" \
              || l_rejtype=="crreject")
	    printf("%-8s %20.5f   \'%-s\'\n","NIFLTSSC",l_sigscale,
	        "Tolerance for sigma clip scaling correction", >> tmphead)
	if (l_rejtype=="pclip")
            printf("%-8s %20.5f   \'%-s\'\n","NIFLTPCL",l_pclip,
	        "Percentile clipping factor used by pclip", >> tmphead)
	if (l_rejtype=="ccdclip")
            printf("%-8s \'%-18s\'   \'%-s\'\n","NIFLTSNO",l_snoise,
	        "Sensitivity noise (e) used by ccdclip", >> tmphead)
    } # end else

    delete(filelist,verify-, >>& "dev$null")
    files(infiles[1],sort-) | unique("STDIN", > filelist)     
    i=1
    scanfile=filelist
    while (fscan(scanfile,img) != EOF) {
	if (substr(img,1,strlen(img)-5) != substr(l_flatfile,1,strlen(l_flatfile)-3)) {
	    printf("%-8s \'%-18s\'   \'%-s\'\n","NIFLON"//str(i),img,
	        "Input sky or GCAL flat (open shutter)", >> tmphead)
	    i+=1
	} 
    }
    delete(filelist,verify-, >& "dev$null")
    if (l_fl_lampsoff) {
	files(infiles[2],sort-) | unique("STDIN", > filelist)
	k=1
	scanfile=filelist
	while (fscan(scanfile,img) != EOF) {
	    if (substr(img,1,strlen(img)-5) != substr(l_flatfile,1,strlen(l_flatfile)-3)) {
                printf("%-8s \'%-18s\'   \'%-s\'\n","NIFLOF"//str(k),img,
		    "Input GCAL flat (closed shutter)", >> tmphead)
                k+=1
	    }
	}
    }

    if (intdebug) {
        print ("l_flatfile="//l_flatfile)
        copy (tmphead, tmphead//".debug")
        print ("Running fxhead BEFORE nhedit on "//\
            substr(l_flatfile,1,strlen(l_flatfile)-3))
        fxhead (substr(l_flatfile,1,strlen(l_flatfile)-3))
    }
        
    # put all the new stuff in the header
    # KL (Jan2011):  This call to nhedit should work on both the fitsutil
    #               and the imutil version of nhedit.
    nhedit(l_flatfile, comfile=tmphead)
    if (intdebug) {
        print ("Running fxhead AFTER nhedit on "//\
            substr(l_flatfile,1,strlen(l_flatfile)-3))
        fxhead (substr(l_flatfile,1,strlen(l_flatfile)-3))
    }
    # fix up header values of read noise, gain, and saturation for flat
    l_gain=i*l_gain*norm
    l_ron=l_ron*l_ron
    l_sat=l_sat/norm
    l_nonlinear=l_nonlinear/norm
    if (l_fl_lampsoff)
        l_ron=sqrt((l_ron*i)+(l_ron*k))/norm
    else
        l_ron=sqrt(l_ron*i)
    l_sat=(int(l_sat*100.))/100.
    l_nonlinear=(int(l_nonlinear*100.))/100.
    l_ron=(int(l_ron*10000.))/10000.
    l_gain=(int(l_gain*10.))/10.
    if (intdebug) {
        print ("Running fxhead BEFORE gain/ron/sat/nonlinear nhedit on "//\
            substr(l_flatfile,1,strlen(l_flatfile)-3))
        fxhead (substr(l_flatfile,1,strlen(l_flatfile)-3))
    }
    gemhedit (l_flatfile, l_key_ron, l_ron, "", delete-)
    gemhedit (l_flatfile, l_key_sat, l_sat, "", delete-)
    gemhedit (l_flatfile, l_key_nonlinear, l_nonlinear, "", delete-)
    gemhedit (l_flatfile, l_key_gain, l_gain, "", delete-)
    
    if (intdebug) {
        print ("Running fxhead AFTER gain/ron/sat/nonlinear nhedit on "//\
            substr(l_flatfile,1,strlen(l_flatfile)-3))
        fxhead (substr(l_flatfile,1,strlen(l_flatfile)-3))
    }

    #--------------------------------------------------------------------------
    # Write the output BPM file, with compete header

    if (l_fl_vardq) {
        l_flatfile=substr(l_flatfile,1,strlen(l_flatfile)-3)
        imcopy(l_flatfile//"["//l_dq_ext//",inherit]",l_bpmfile, >& "dev$null")
    }

    #--------------------------------------------------------------------------
    # clean up

clean:
    delete(filelist//","//scilist//","//tmphead,verify-, >>& "dev$null")
    for (i=1; i<=3; i+=1) {
	imdelete(combout[i],verify-, >>& "dev$null")
    }
    imdelete(tmpdq//"*.pl",verify-, >>& "dev$null")
    if (status==1) imdelete(l_flatfile,verify-, >>& "dev$null")
    scanfile=""
    # close log file
    if (status==0) {
	printlog(" ",l_logfile,l_verbose)
	printlog("NIFLAT exit status:  good.",l_logfile,l_verbose)
    }
    printlog("----------------------------------------------------------------------------",
        l_logfile,l_verbose)

end
