# Copyright(c) 2001-2006 Association of Universities for Research in Astronomy, Inc.

procedure oreduce(inimages)

# Basic reductions of OSCIR images. 
# Difference and add the chop-nod postions.
# Optional flat fielding.
#
# Flatfield is selected on the basis key_filter
# 6 different flat fields can be used.
#
# Version: Sept 14, 2002 IJ,BR  Release v1.4
#          Aug 20, 2003  KL  IRAF2.12 - new/modified parameters
#                              imcombine: headers,bpmasks,expmasks,outlimits
#                                    rejmask->rejmasks, plfile->nrejmasks

char  inimages     {prompt="Input OSCIR image(s)"}
char  outimages    {"",prompt="Output image(s)"}
char  outprefix    {"r",prompt="Prefix for output image(s)"}
char  normtype     {"sec",min="sec|raw|frame",
                      prompt="Normalization of output image(s)"}
char  logfile      {"",prompt="Logfile"}
bool  fl_chopnod   {yes,prompt="Combine chop and nod positions"}
char  bnod         {"",prompt="Bad nodsets"}
char  bset         {"",prompt="Bad savesets in the bad nodsets"}
char  frameeff     {"oscir$data/oscirframeeff",prompt="Frame efficiency image"}
bool  fl_flat      {no,prompt="Do flat fielding"}
char  key_filter   {"FILTER",prompt="Keyword for filter id"}
char  flatimage1   {"",prompt="Flat field image no.1"}
char  filter1      {"N_wide",prompt="Filter for flat field no.1"}
char  flatimage2   {"",prompt="Flat field image no.2"}
char  filter2      {"IHW_(17-19)",prompt="Filter for flat field no.2"}
char  flatimage3   {"",prompt="Flat field image no.3"}
char  filter3      {"Q3_(20-22)",prompt="Filter for flat field no.3"}
char  flatimage4   {"",prompt="Flat field image no.4"}
char  filter4      {"S_10.3",prompt="Filter for flat field no.4"}
char  flatimage5   {"",prompt="Flat field image no.5"}
char  filter5      {"S_11.7",prompt="Filter for flat field no.5"}
char  flatimage6   {"",prompt="Flat field image no.6"}
char  filter6      {"S_12.5",prompt="Filter for flat field no.6"}
bool  verbose      {no,prompt="Verbose"}
int   status       {0,prompt="Exit status (0=good)"}
struct* scanfile

begin

    char l_inimages, l_outimages, l_flatimage, l_filter, l_bnod, l_bset
    char l_frameeff
    char l_prefix, l_logfile, l_temp, tmpin, tmpfile
    char in[100], out[100]
    char l_flatimage1, l_flatimage2, l_flatimage3, l_flatimage4, l_flatimage5 
    char l_flatimage6, l_filter1, l_filter2, l_filter3, l_filter4, l_filter5
    char l_filter6, l_keyfilter, l_normtype
    int  i, nimages, noutimages, maxfiles
    bool l_fl_chopnod, l_fl_flat, l_verbose, l_fl_skytemp, l_fl_flattemp
    bool flatok, l_fl_print, l_fl_dollar
    int  i_bnod[100], i_bset[100], n_bnod, n_bset, i_temp
    real n_norm, n_exptime, n_frmtime, n_efffrm
    struct l_struct, n_title

    int  n_saveset, n_nodset, n_frmcoadd, n_chpcoadd, n_i, n_j, n_jj, n_total
    real n_ave1, n_ave2
    char n_dif1, n_dif2, tmpimage, tmpimage2, tmpimage3
    char l_test, l_junk, tmpbnod


    status=0
    unlearn("imexpr")
    maxfiles=100  # This is the maximum number of files that can
        	   # be reduced with one call to oreduce
    tmpfile = mktemp("tmpin")
    tmpin = mktemp("tmpin")
    tmpbnod = mktemp("tmpbnod")

    # cache imgets - used throughout the script
    cache("imgets", "gemdate")

    # set the local variables
    l_inimage=inimage ; l_outimage=outimage 
    l_fl_chopnod=fl_chopnod ; l_bnod=bnod ; l_bset=bset ; l_frameeff=frameeff
    l_fl_flat=fl_flat ; l_verbose=verbose ; l_prefix=outprefix
    l_normtype=normtype ; l_logfile=logfile 
    l_filter1=filter1 ; l_filter2=filter2 ; l_filter3=filter3
    l_filter4=filter4
    l_filter5=filter5 ; l_filter6=filter6 ; l_keyfilter=key_filter
    l_flatimage1=flatimage1 ; l_flatimage2=flatimage2 ; l_flatimage3=flatimage3 
    l_flatimage4=flatimage4 ; l_flatimage5=flatimage5 ; l_flatimage6=flatimage6

    # Check for package log file or user-defined log file
    cache("oscir")
    if ((l_logfile=="") || (l_logfile==" ")) {
        l_logfile=oscir.logfile
        if ((l_logfile=="") || (l_logfile==" ")) {
	    l_logfile="oscir.log"
	    printlog("WARNING - OREDUCE:  Both oreduce.logfile and \
	        oscir.logfile are empty.",logfile=l_logfile, verbose+)
	    printlog("                 Using default file oscir.log.",
		logfile=l_logfile, verbose+)
       }
    }
    # Open log file
    date | scan(l_struct)
    printlog("--------------------------------------------------------------\
        --------------",logfile=l_logfile, verbose=l_verbose)
    printlog("OREDUCE -- "//l_struct, logfile=l_logfile, verbose=l_verbose)
    printlog(" ",logfile=l_logfile, verbose=l_verbose)

    #-----------------------------------------------------------------------
    # Load up arrays of input name lists

    # Make list if * in inimages
    if (stridx("*",l_inimages)>0) {
	files(l_inimages, > tmpin)
	l_inimages="@"//tmpin
    }

    #
    nimages=0
    if (substr(l_inimages,1,1)=="@") 
	scanfile=substr(l_inimages,2,strlen(l_inimages))
    else {
	files(l_inimages,sort-, > tmpfile)
	scanfile=tmpfile
    }

    l_fl_dollar=no
    while (fscan(scanfile,l_temp) != EOF) {
	if (substr(l_temp,strlen(l_temp)-4,strlen(l_temp))==".fits" )
	    l_temp=substr(l_temp,1,strlen(l_temp)-5)

	if (!imaccess(l_temp))
	    printlog("WARNING - OREDUCE: Input image "//l_temp//" not found.",
	        logfile=l_logfile,verbose+)
	else {
	    nimages=nimages+1
	    if (nimages > maxfiles) {
		printlog("ERROR - OREDUCE: Maximum number of input images \
		    exceeded ("//str(maxfiles)//")",logfile=l_logfile,verbose+)
		status=1
		goto clean
	    }
	    in[nimages]=l_temp 

	    # Catch $  and / in input
	    if (stridx("$",in[nimages])!=0)
		l_fl_dollar=yes
	    if (stridx("/",in[nimages])!=0)
		l_fl_dollar=yes
	}
    }
    printlog("Processing "//nimages//" file(s).",
        logfile=l_logfile,verbose=l_verbose)
    scanfile=""
    delete(tmpfile//","//tmpin,ver-, >& "dev$null")
    if (nimages==0) {
	printlog("ERROR - OREDUCE: No existing input images defined",
	    l_logfile,verbose+)
	status=1
	goto clean
    }

    # Now for the output images
    # outimages could contain legal * if it is of a form like %st%stX%*.imh

    noutimages=0
    if (l_outimages!="" ) {
	if (substr(l_outimages,1,1)=="@") 
	    scanfile=substr(l_outimages,2,strlen(l_outimages))
	else if (stridx("*",l_outimages)>0)  {
	    files(l_outimages,sort-) | \
	        match(".hhd",stop+,print-,metach-, > tmpfile)
	    scanfile=tmpfile
	} else {
	    files(l_outimages,sort-, > tmpfile)
	    scanfile=tmpfile
	}

	while (fscan(scanfile,l_temp) != EOF) {
	    noutimages=noutimages+1
	    if (noutimages > maxfiles) {
		printlog("ERROR - OREDUCE: Maximum number of input images \
		    exceeded ("//str(maxfiles)//").",logfile=l_logfile,verbose+)
		status=1
		goto clean
	    }
	    out[noutimages]=l_temp 
	    if (imaccess(out[noutimages])) {
		printlog("ERROR - OREDUCE: Output image "//out[noutimages]//\
		    " exists.",logfile=l_logfile,verbose+)
		status=1
		goto clean
	    }
	}
    }
    scanfile="" 
    delete(tmpfile,ver-, >& "dev$null")

    # if there are too many or too few output images, and any defined
    # at all at this stage - exit with error
    if (nimages!=noutimages && l_outimages!="") {
	printlog("ERROR - OREDUCE: Number of input and output images mismatch.",
	    logfile=l_logfile,verbose+)
	status=1
	goto clean
    }

    # If prefix is to be used instead
    if (l_outimages=="" ) {
	if ((l_prefix=="") || (l_prefix==" ")) {
	    printlog("ERROR - OREDUCE: Neither output image name or output \
	        prefix is defined.",logfile=l_logfile,verbose+)
	    status=1
	    goto clean
	}
	if (l_fl_dollar) {
	    printlog("ERROR - OREDUCE: Cannot use outprefix with path as \
	        part of input image names",l_logfile,verbose+)
	    printlog("                 Set outimages to avoid this error",
		l_logfile,verbose+)
	    status=1
	    goto clean
	}
	i=1
	while (i<=nimages) {
	    out[i]=l_prefix//in[i]
	    if (imaccess(out[i])) {
		printlog("ERROR - OREDUCE: Output image "//out[i]//" exists.",
		    logfile=l_logfile,verbose+)
		status=1
		goto clean
	    }
	    i=i+1
	}
    }

    # Load arrays for bnod and bset if !="" and only one input image
    if (nimages>1 && (l_bnod!="" || l_bset!="")) {
	printlog("WARNING - OREDUCE: bnod and bset cannot be set if more \
	    than one input image",l_logfile,verbose+)
	printlog("                   Settings ignored.",l_logfile,verbose+)
	l_bnod="" ; l_bset=""
    }

    if (l_bnod!="" && l_bnod!=" ") {
	files(l_bnod,sort-, > tmpbnod)
	n_bnod=0
	scanfile=tmpbnod
	while (fscan(scanfile,i_temp)!=EOF) {
	    n_bnod+=1
	    i_bnod[n_bnod]=i_temp
	}
	scanfile=""
    } else {
	i_bnod[1]=0
	n_bnod=1
    }
    delete(tmpbnod,verify-, >& "dev$null")

    if (l_bset!="" && l_bset!=" ") {
	files(l_bset,sort-, > tmpbnod)
	n_bset=0
	scanfile=tmpbnod
	while (fscan(scanfile,i_temp)!=EOF) {
	    n_bset+=1
	    i_bset[n_bset]=i_temp
	}
	scanfile=""
    } else {
	for (n_i=1;n_i<=n_bnod;n_i+=1)
	    i_bset[n_i]=0
	n_bset=n_bnod
    }
    delete(tmpbnod,verify-, >& "dev$null")

    if (n_bset!=n_bnod) {
	printlog("WARNING - OREDUCE: Entries for bset must match entries for \
	    bnod",l_logfile,verbose+)
	printlog("                   Ignoring bset",l_logfile,verbose+)
	l_bset=""
	for (n_i=1;n_i<=n_bnod;n_i+=1) 
	    i_bset[n_i]=0
    }

    if (!imaccess(l_frameeff) && l_frameeff!="" && l_frameeff!=" " && \
        l_normtype=="sec" && l_fl_chopnod) {
	printlog("WARNING - OREDUCE: Frame efficiency image not found",
	    l_logfile,verbose+)
	l_frameeff=""
    }

    if ((l_frameeff=="" || l_frameeff==" ") && l_normtype=="sec" && l_fl_chopnod) {
	printlog("Frame efficiency for normalization to seconds: Using \
	    average value",l_logfile,l_verbose)
	l_frameeff=""
    } else if (l_normtype=="sec" && l_fl_chopnod) {
	printlog("Frame efficiency image for normalization to seconds: "//\
	    l_frameeff,l_logfile,l_verbose)
    }

    #--------------------------------------------------------------------------
    # The math and bookkeeping:  (MAIN LOOP)

    i=1
    while (i<=nimages) {
        # Create tmp FITS file names used within this loop
	n_dif1 = mktemp("tmpdif1")
        n_dif2 = mktemp("tmpdif2")
	tmpimage = mktemp("tmpimag")
	tmpimage2 = mktemp("tmpimag")
	tmpimage3 = mktemp("tmpimag")


   	printlog(" ",l_logfile,l_verbose)
	printlog("Input  image : "//in[i],l_logfile,l_verbose)
	printlog("Output image : "//out[i],l_logfile,l_verbose)

	# use this to figure out if the file is valid FITS file
	l_junk="" ; l_test=""
	imhead(in[i],imlist = "*.imh,*.fits,*.pl,*.qp,*.hhh",long-,
	    userfields+) |& scan(l_junk,l_test)
	if (l_test=="Negative") {
	    printlog("WARNING - OREDUCE: Image "//in[i]//" not a valid FITS \
	        file, not processed",l_logfile,verbose+)
	    goto nextimage
	}

	#----------------
	# check for previous chopnod combination turn it off if necessary
	if (l_fl_chopnod) {
	    l_fl_skytemp=yes
	    imgets(in[i],"i_naxis",>& "dev$null")
	    if (imgets.value == "2") {
		l_fl_skytemp=no
		printlog("WARNING - OREDUCE: Chop-nods already combined for \
		    image",logfile=l_logfile,verbose+)
		printlog("                   Chop-nod combination not \
		    performed.",logfile=l_logfile,verbose+)
	    }
	    # end if (l_fl_chopnod)
	} else
	    l_fl_skytemp=no

	#----------------
	# check for previous flat fielding, and turn it off if necessary
	l_filter="none" ; l_flatimage="none"  # for logging
	if (l_fl_flat) {
	    imgets(in[i],"FLATIMAG",>& "dev$null")
	    if (imgets.value != "0") {
		l_fl_flattemp=no
		l_flatimage="none"
		printlog("WARNING - OREDUCE: Image "//in[i]//" has already \
		    been flat-fielded.",logfile=l_logfile,verbose+)
		printlog("                   by oreduce.  Flat-fielding not \
		    performed.",logfile=l_logfile,verbose+)
		l_filter="none"  # set l_filter to something, for logging
	    } else {
		l_fl_flattemp=yes
		imgets(in[i],l_keyfilter,>& "dev$null")
		print(imgets.value) | scan(l_filter)  # strip off encoder position
		if (i==1)
		    l_temp=imgets.value

		# Find the right flat field by comparing the filter name 
		# from the header
		flatok=no
		if (l_filter == filter1 && flatimage1!="") {
		    l_flatimage=flatimage1 ; flatok=yes
		    goto flatset
		}
		if (l_filter == filter2 && flatimage2!="") {
		    l_flatimage=flatimage2 ; flatok=yes
		    goto flatset
		}
		if (l_filter == filter3 && flatimage3!="") {
		    l_flatimage=flatimage3 ; flatok=yes
		    goto flatset
		}
		if (l_filter == filter4 && flatimage4!="") {
		    l_flatimage=flatimage4 ; flatok=yes
		    goto flatset
		}
		if (l_filter == filter5 && flatimage5!="") {
		    l_flatimage=flatimage5 ; flatok=yes
		    goto flatset
		}
		if (l_filter == filter6 && flatimage6!="") {
		    l_flatimage=flatimage6 ; flatok=yes
		    goto flatset
		}

flatset:
		if (!flatok) {
		    printlog("WARNING - OREDUCE: Image "//in[i]//" is taken \
		        in filter "//l_filter,logfile=l_logfile,verbose+)
		    printlog("                   Flat field not defined for \
		        this filter.",logfile=l_logfile,verbose+)
		    printlog("                   Flat-fielding not performed.",
		        logfile=l_logfile,verbose+)
		    l_fl_flattemp = no ; l_flatimage="none"
		} else {  # START of FLATOK

		    if ((i==1) || (l_temp != l_filter)) {
			if(!imaccess(l_flatimage)) {
			    printlog("ERROR - OREDUCE: Flat field "//\
			        l_flatimage//" not found.",
				logfile=l_logfile,verbose+)
			    status=1
			    goto clean
			}
		    }
		    l_temp=l_filter
		}
	    }
        } else  {	# NO flat fielding
	    l_fl_flattemp=no
	    l_flatimage="none"
	}

	# Check if attempting to flat field raw image - disallowed
	if (!l_fl_skytemp && l_fl_flat) {
	    imgets(in[i],"i_naxis",>& "dev$null")
	    if (imgets.value != "2") {
		l_fl_flattemp=no
		printlog("WARNING - OREDUCE: Attempting to flat field raw \
		    image before combining",l_logfile,verbose+)
		printlog("                   chop/nod positions. Need to set \
		    fl_chop=yes",l_logfile,verbose+)
	    }
	}

	#-------------
	# Check if both flags are no
	if (!l_fl_flattemp && !l_fl_skytemp) {
	    printlog("WARNING - OREDUCE: Image "//in[i]//" no reductions \
	        selected",l_logfile,l_verbose)
	    goto nextimage
	}

	#-------------
	# Combine chops and nods
	if (l_fl_skytemp) {
	    # Check/Get the dimensions
	    imgets(in[i],"i_naxis3")
	    if(imgets.value!="2") {
		printlog("WARNING - OREDUCE: Number of chop positions not 2",
		    l_logfile,verbose+)
		goto nextimage
	    }
	    imgets(in[i],"i_naxis5")
	    if (imgets.value!="2") {
		printlog("WARNING - OREDUCE: Number of nod positions not 2",
		    l_logfile,verbose+)
		goto nextimage
	    }

	    imgets(in[i],"i_naxis4")
	    if (imgets.value=="0") {
		printlog("WARNING - OREDUCE: No savesets",l_logfile,verbose+)
		goto nextimage
	    } else
		n_saveset=int(imgets.value)
		
	    imgets(in[i],"i_naxis6")
	    if (imgets.value=="0") {
		printlog("WARNING - OREDUCE: No nodsets",l_logfile,verbose+)
		goto nextimage
	    } else
		n_nodset=int(imgets.value)

	    printlog("Savesets="//n_saveset//"   Nodsets="//n_nodset,
	        l_logfile,l_verbose)
	    printlog("Bad nodsets="//l_bnod,l_logfile,l_verbose)
	    printlog("Individually defined bad savesets="//l_bset,
	        l_logfile,l_verbose)

	    # Combine all differences from nod position 1
	    for (n_j=1;n_j<=n_nodset;n_j+=1) {
		for (n_i=1;n_i<=n_saveset;n_i+=1) {
		    l_fl_print=yes
		    for (n_jj=1;n_jj<=n_bnod;n_jj+=1) {
			if (n_j==i_bnod[n_jj] && (n_i==i_bset[n_jj] || i_bset[n_jj]==0) ) 
			    l_fl_print=no
		    }

		    if (l_fl_print) {
			print(in[i]//"[*,*,1,"//str(n_i)//",1,"//str(n_j)//"]",
			    >> tmpin)
			print(in[i]//"[*,*,2,"//str(n_i)//",1,"//str(n_j)//"]",
			    >> tmpfile)
		    }

		}
	    }
	    imcombine("@"//tmpin,n_dif1,headers="",bpmasks="",rejmasks="",
		nrejmasks="",expmasks="",sigmas="",logfile="",combine="average",
		reject="none",project=no,outtype="double",outlimits="",
		offsets="none",masktype="none",maskvalue=0.,blank=0.,
		scale="none",zero="none",weight="none",statsec="",expname="",
		lthreshold=INDEF,hthreshold=INDEF)
	    imcombine("@"//tmpfile,tmpimage,headers="",bpmasks="",rejmasks="",
		nrejmasks="",expmasks="",sigmas="",logfile="",combine="average",
		reject="none",project=no,outtype="double",outlimits="",
		offsets="none",masktype="none",maskvalue=0.,blank=0.,
		scale="none",zero="none",weight="none",statsec="",expname="",
		lthreshold=INDEF,hthreshold=INDEF)
	    imarith(n_dif1,"-",tmpimage,n_dif1)
	    imdelete(tmpimage,verify-)
	    delete(tmpin//","//tmpfile,verify-)

	    # Combine all differences from nod position 2
	    for (n_j=1;n_j<=n_nodset;n_j+=1) {
		for (n_i=1;n_i<=n_saveset;n_i+=1) {
		    l_fl_print=yes
		    for (n_jj=1;n_jj<=n_bnod;n_jj+=1) {
			if (n_j==i_bnod[n_jj] && (n_i==i_bset[n_jj] || i_bset[n_jj]==0) )
			    l_fl_print=no
		    }

		    if (l_fl_print) {
			print(in[i]//"[*,*,2,"//str(n_i)//",2,"//str(n_j)//"]",
			    >> tmpin)
			print(in[i]//"[*,*,1,"//str(n_i)//",2,"//str(n_j)//"]",
			    >> tmpfile)
		    }
		}
	    }
	    imcombine("@"//tmpin,n_dif2,headers="",bpmasks="",rejmasks="",
	        nrejmasks="",expmasks="",sigmas="",logfile="",combine="average",
		reject="none",project=no,outtype="double",outlimits="",
		offsets="none",masktype="none",maskvalue=0.,blank=0.,
		scale="none",zero="none",weight="none",statsec="",expname="",
		lthreshold=INDEF,hthreshold=INDEF)
	    imcombine("@"//tmpfile,tmpimage2,headers="",bpmasks="",rejmasks="",
	        nrejmasks="",expmasks="",sigmas="",logfile="",combine="average",
		reject="none",outtype="double",outlimits="",offsets="none",
		masktype="none",maskvalue=0.,blank=0.,scale="none",zero="none",
		weight="none",statsec="",expname="",lthreshold=INDEF,
		hthreshold=INDEF)
	    imarith(n_dif2,"-",tmpimage2,n_dif2)
	    imdelete(tmpimage2,verify-)
	    count(tmpin) | scan(n_total)
	    delete(tmpin//","//tmpfile,verify-)

	    # Derive final chop-nod combination
	    # n_total is the actual number of n_nodset*n_saveset that were 
	    # combined
	    imexpr("(a+b)*"//real(n_total),out[i],n_dif1,n_dif2,verbose-)
	    imdelete(n_dif1//","//n_dif2,verify-)

	    #-------------
	    # Exposure time and normalization
	    imgets(in[i],"FRMCOADD")
	    n_frmcoadd=int(imgets.value)
	    imgets(in[i],"CHPCOADD")
	    n_chpcoadd=int(imgets.value)
	    if (n_frmcoadd*n_chpcoadd<=0) {
		printlog("WARNING - OREDUCE: Header info missing for FRMCOADD \
		    or CHPCOADD",l_logfile,verbose+)
		n_frmcoadd=1
		n_chpcoadd=1
	    }
	    imgets(in[i],"FRMTIME") ; n_frmtime=real(imgets.value)/1000.  # sec
	    #imgets(in[i],"EFF_FRM") ; n_efffrm=real(imgets.value)/100.    # fraction
	    n_efffrm=1018.5/1024.0

	    # n_total is the actual number of n_nodset*n_saveset that were 
	    # combined
	    n_exptime=n_frmcoadd*n_chpcoadd*n_frmtime*n_efffrm
	    if (l_normtype=="sec")
	        n_norm=n_exptime*n_total*2.
	    else if (l_normtype=="frame")
		n_norm=n_total*2.
	    else
		n_norm=1.

	    n_exptime=n_exptime*n_total*2.
	    if (n_norm!=1. && (l_frameeff=="" || l_normtype!="sec")) 
		imarith(out[i],"/",n_norm,out[i])

	    if (n_norm!=1. && l_frameeff!="" && l_frameeff!=" " && l_normtype=="sec") {
		imrename(out[i],tmpimage3,verbose-)
		imexpr("a/(b*"//str(n_frmcoadd*n_chpcoadd*n_frmtime*n_total*2.)//")",
		    out[i],tmpimage3,l_frameeff,verbose-)
		imdelete(tmpimage3,verify-)
	    }

	    chpixtype(out[i],out[i],"real",oldpixtype="all",verbose-)

	    printf("Number of on-source frames obtained at each nod position = %4d\n",
	        (n_nodset*n_saveset)) | scan(l_struct)
	    printlog(l_struct,l_logfile,l_verbose)
	    printf("Number of on-source frames used from each nod position   = %4d\n",
	        n_total) | scan(l_struct)
	    printlog(l_struct,l_logfile,l_verbose)
	    printf("Total average exposure time of co-added image  = %10.4f\n",
	        n_exptime) | scan(l_struct)
	    printlog(l_struct,l_logfile,l_verbose)
	    printf("Normalization value for co-added image         = %10.4f\n",
	        n_norm) | scan(l_struct)
	    printlog(l_struct,l_logfile,l_verbose)
	    printlog("Normalization type for co-added image  : "//l_normtype,
	        l_logfile,l_verbose)

        } # end of fl_skytemp

        #-------------

        # Flat fielding
	if (l_fl_flattemp && l_fl_skytemp) {
	    printlog("Flat field : "//l_flatimage//"  Filter : "//l_filter,
		l_logfile,l_verbose)
	    imarith(out[i]//".fits","/",l_flatimage//".fits",out[i]//".fits",
		pixtype="real",verbose=no,title="")
	} else if (l_fl_flattemp && !l_fl_skytemp) {
	    printlog("Flat field : "//l_flatimage//"  Filter : "//l_filter,
		l_logfile,l_verbose)
	    imarith(in[i]//".fits","/",l_flatimage//".fits",out[i]//".fits",
		pixtype="real",verbose=no,title="")
	}

	# and update the header
	gemdate ()
	# date stamp the modification
	gemhedit(out[i],"GEM-TLM",gemdate.outdate,"UT Last modification with GEMINI")
	gemhedit(out[i],"OREDUCE",gemdate.outdate,"UT Time stamp for oreduce")
	if (l_fl_skytemp) {
	    gemhedit(out[i],"EXPCOADD",n_exptime,
	        "Total exposure time for co-add")
	    gemhedit(out[i],"NORMTYPE",l_normtype,
	        "Normalization type for co-add")
	    gemhedit(out[i],"NORMVAL",n_norm,"Normalization value for co-add")
	}
	if (l_fl_flattemp)
	    gemhedit(out[i],"FLATIMAG",l_flatimage,"Flat image used by oreduce")

nextimage:   # Jump point for "invalid" image
	i=i+1
    } # end the main loop

clean:
    #---------------------------------------------------------------------------
    # Clean up
    if (status==0)
	printlog("OREDUCE exit status:  good.",
	    logfile=l_logfile, verbose=l_verbose)
    printlog("--------------------------------------------------------------\
        --------------", logfile=l_logfile, verbose=l_verbose)
    scanfile=""
    delete(tmpfile//","//tmpin,ver-, >& "dev$null")

end
