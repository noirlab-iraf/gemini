# Copyright(c) 2000-2011 Association of Universities for Research in Astronomy, Inc.

procedure qflat (lampson)

# Derive QUIRC bad pixel mask and flat
#
# Version   Jan 26, 2001 JJ,BM,IJ
#           Aug 20, 2003 KL  IRAF2.12 - new/modified parameters
#                              hedit: addonly
#                              imcombine: headers,bpmasks,expmasks,outlimits
#                                    rejmask->rejmasks, plfile->nrejmasks
#                              imstat: nclip,lsigma,usigma,cache

string lampson      {prompt="Input images (sky flats or lamps on)"}
string lampsoff     {"",prompt="Input images (lamps off)"}
string darks        {"",prompt="Input images (darks)"}
string flatimage    {"",prompt="Output flat field image"}
string flattitle    {"default",prompt="Title for flat image"}
string bpmimage     {"",prompt="Output bad pixel file"}
string bpmtitle     {"default",prompt="Title for bad pixel image"}
string logfile      {"",prompt="Logfile name"}
real   thresh_flo   {1000.,prompt="Lower bad-pixel threshold for flats"}
real   thresh_fup   {40000.,prompt="Upper bad-pixel threshold for flats"}
real   thresh_dlo   {-50.,prompt="Lower bad-pixel threshold for darks"}
real   thresh_dup   {200.,prompt="Upper bad-pixel threshold for darks"}
string normstat     {"mean",enum="mean|midpt",prompt="Statistic to use for normalization."}
string combtype     {"default",enum="default|average|median",prompt="Type of combine operation"}
string rejtype      {"avsigclip",enum="none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip",
                        prompt="Type of rejection"}
string scale        {"none",prompt="Image scaling"}
string zero         {"none",prompt="Image zeropoint offset"}
string weight       {"none",prompt="Image weights"}
string statsec      {"[*,*]",prompt="Statistics section"}
string key_exptime  {"EXPTIME",prompt="Header keyword for exposure time"}
real   lthreshold   {INDEF,prompt="Lower threshold"}
real   hthreshold   {INDEF,prompt="Upper threshold"}
int    nlow         {1,min=0,prompt="minmax: Number of low pixels to reject"}
int    nhigh        {1,min=0,prompt="minmax: Number of high pixels to reject"}
int    nkeep        {0,prompt="Minimum to keep or maximum to reject"}  
bool   mclip        {yes,prompt="Use median in sigma clipping algorithms?"}
real   lsigma       {3.,min=0.,prompt="Lower sigma clipping factor"}
real   hsigma       {3.,min=0.,prompt="Upper sigma clipping factor"}
string key_ron      {"RON",prompt="Keyword for readout noise in e-"}
string key_gain     {"GAIN",prompt="Keyword for gain in electrons/ADU"}
real   ron          {15.0,min=0.,prompt="Default read noise if header keyword absent"}
real   gain         {1.85,min=0.01,prompt="Default gain to use if header keyword absent"}
string snoise       {"0.0",prompt="ccdclip: Sensitivity noise (electrons)"}
real   sigscale     {0.1,min=0.,prompt="Tolerance for sigma clipping scaling correction"}
real   pclip        {-0.5,prompt="pclip: Percentile clipping parameter"}
real   grow         {0.0,min=0.0,prompt="Radius (pixels) for neighbor rejection"}
bool   fl_rmstars   {no,prompt="Remove stars in flats?"}
real   fwhmpsf      {4.,min=0.5,prompt="PSF FWHM in pixels for finding stars"}
bool   fl_gemseeing {no,prompt="Use GEMSEEING to improve QSKY results?"}
bool   fl_keepmasks {no,prompt="Keep QSKY object masks for each input image?"}
bool   verbose      {no,prompt="Verbose output?"}
int    status       {0,prompt="Exit status (0=good)"}
struct  *flist      {"",prompt="Internal use only"} 

begin
        
    string l_lampson,l_lampsoff,l_darks,l_bpmfile, l_statsec, l_expname
    string l_combine, l_reject, l_scale, l_zero, l_weight, l_flatfile
    string l_key_ron, l_key_gain, l_snoise, l_bpmtitle, l_flattitle, l_logfile
    real   l_ron, l_gain, l_lthreshold, l_hthreshold, l_lsigma, l_hsigma
    real   l_grow, l_sigscale, l_pclip
    real   l_lower, l_upper, l_sigthresh
    int    l_nlow, l_nhigh, l_nkeep
    bool   l_fl_flat, l_fl_dark, l_fl_bpm, l_fl_lampsoff
    bool   l_verbose, l_mclip, l_fl_rmstars, l_fl_keepmasks, l_fl_gemseeing
    string l_sgain, l_sron, l_normstat
    real   l_flat_lo, l_flat_up, l_thresh_flo, l_thresh_fup
    real   l_thresh_dlo, l_thresh_dup, l_fwhmpsf
    string filelist,infiles[3]
    string img,dqsum,dqsumold,suf
    string tmpsci,tmpdq,combout[3],combsig[3],combdq[3]
    string dsub,mask[2]
    real   norm,stddev
    real   etime[3],rtime
    int    i,k,n,nbad,len
    struct sdate

    # Query parameters
    l_lampson=lampson ; l_lampsoff=lampsoff ; l_darks=darks
    l_bpmfile=bpmimage ;  l_logfile=logfile ; l_flatfile=flatimage
    l_bpmtitle=bpmtitle ; l_flattitle=flattitle
    l_combine=combtype ; l_reject=rejtype
    l_scale=scale ; l_zero=zero ; l_weight=weight
    l_statsec=statsec ; l_expname=key_exptime
    l_lthreshold=lthreshold ; l_hthreshold=hthreshold
    l_nlow=nlow ; l_nhigh=nhigh ; l_nkeep=nkeep
    l_lsigma=lsigma ; l_hsigma=hsigma
    l_verbose=verbose ; l_mclip=mclip
    l_key_gain=key_gain ; l_key_ron=key_ron
    l_gain=gain ; l_ron=ron ;l_snoise=snoise 
    if (l_snoise=="" || l_snoise==" ")
	l_snoise="0.0"
    l_sigscale=sigscale ; l_pclip=pclip ; l_grow=grow
    l_thresh_flo=thresh_flo ; l_thresh_fup=thresh_fup ; l_thresh_dup=thresh_dup
    l_fl_rmstars=fl_rmstars ; l_normstat=normstat ; l_thresh_dlo=thresh_dlo
    l_fwhmpsf=fwhmpsf ; l_fl_keepmasks=fl_keepmasks
    l_fl_gemseeing=fl_gemseeing

    status=0

    # Keep imgets parameters from changing by outside world
    cache ("imgets", "gemdate")

    # logfile name
    if ((l_logfile=="") || (l_logfile==" ")) {   
	l_logfile=quirc.logfile
	if ((l_logfile=="") || (l_logfile==" ")) {  
            l_logfile="quirc.log"
            printlog("WARNING - QFLAT:  Both qflat.logfile and quirc.logfile \
	        are",l_logfile, yes) 
            printlog("                 empty.  Using quirc.log.",l_logfile,yes)
	}
    }

    # Start log file
    date | scan(sdate)
    printlog("--------------------------------------------------------------\
        ------------",l_logfile,verbose=l_verbose)
    printlog("QFLAT -- "//sdate,l_logfile,verbose=l_verbose)
    printlog(" ",l_logfile,verbose=l_verbose)

    # Set titles
    if (l_flattitle=="default" || l_flattitle=="" || l_flattitle==" ") 
	l_flattitle="FLAT FIELD IMAGE from gemini.quirc.qflat"
    if (l_bpmtitle=="default" || l_bpmtitle=="" || l_bpmtitle==" ") 
	l_bpmtitle="BAD PIXEL MASK from gemini.quirc.qflat"

    # Define temporary files
    filelist=mktemp("tmpfilelist")
    infiles[1]=mktemp("tmponlist")
    infiles[2]=mktemp("tmpofflist")
    infiles[3]=mktemp("tmpdarklist")
    combout[1]=mktemp("tmpon")
    combout[2]=mktemp("tmpoff")
    combout[3]=mktemp("tmpdark")
    combsig[1]=mktemp("tmponsig")
    combsig[2]=mktemp("tmpoffsig")
    combsig[3]=mktemp("tmpdrksig")
    combdq[1]=mktemp("tmpcombdqon")
    combdq[2]=mktemp("tmpcombdqoff")
    combdq[3]=mktemp("tmpcombdqdrk")
    dsub=mktemp("tmpdsub")
    mask[1]=mktemp("tmpmask1")
    mask[2]=mktemp("tmpmask2")

    infiles[1]=l_lampson ; infiles[2]=l_lampsoff ; infiles[3]=l_darks

    #------------------------------------------------------------------------
    # Check output files for empty strings, files that already exist
    if ((l_flatfile=="") || (l_flatfile==" ")) {
	printlog("WARNING - QFLAT:  No output flatfield image was specified.  \
	    Flat will not",l_logfile,yes)
	printlog("                  be written.",l_logfile,yes)
	l_fl_flat=no
    } else {
	l_fl_flat=yes
	if (imaccess(l_flatfile)) {
	    printlog("ERROR - QFLAT: Output file "//l_flatfile//" already \
	        exists.",l_logfile,yes)
	    goto error
	}
    }

    if ((l_bpmfile=="") || (l_bpmfile==" ")) {
	printlog("WARNING - QFLAT:  No output bad pixel mask was specified.  \
	    Mask will not",l_logfile,yes)
	printlog("                  be written.",l_logfile,yes)
	l_fl_bpm=no
    } else {
	l_fl_bpm=yes
	if (imaccess(l_bpmfile)) {
	    printlog("ERROR - QFLAT: Output file "//l_bpmfile//" already \
	        exists.",l_logfile,yes)
	    goto error
	}
    }

    if (!l_fl_flat && !l_fl_bpm) {
	printlog("ERROR - QFLAT:  No output files specified.",l_logfile,yes)
	goto error
    }
    if ((l_lampson=="") || (l_lampson==" ")) {
	printlog("ERROR - QFLAT:  No input lamps-on images specified.",
	    l_logfile,yes)
	goto error
    }

    if ((l_lampsoff=="") || (l_lampsoff==" "))
	l_fl_lampsoff=no
    else
	l_fl_lampsoff=yes

    if ((l_darks=="") || (l_darks==" "))
	l_fl_dark=no
    else
	l_fl_dark=yes

    # check that files exist 
    nbad=0
    for (i = 1; i<=3; i+=1) {
	#check input files
	if (infiles[i] != "" && infiles[i] != " ") {
	    # check existence of list file
            if (substr(infiles[i],1,1) == "@") {
        	len=strlen(infiles[i])
        	if (!access(substr(infiles[i],2,len))) {
                    printlog("ERROR - QFLAT: "//substr(infiles[i],2,len)//\
		        " does not exist.",l_logfile,yes)
                    goto error
        	}
            }

            files(infiles[i],sort-, > filelist)

            # Check that all images in the list exist
            flist=filelist
            while (fscan(flist, img) != EOF) {
        	if (!imaccess(img)) {
                    printlog("WARNING - QFLAT: "//img//" not found",logfile,yes)
                    nbad=nbad+1
        	}
            }
            flist=""
            delete(filelist,verify-, >& "dev$null")
	} 
    } 
    #end of for loop

    # Exit if problems found
    if (nbad > 0) {
	printlog("ERROR - QFLAT: "//nbad//" image(s) do not exist",
	    l_logfile,yes)
	goto error
    }

    # Check if key_exptime is set
    if (l_expname=="" || l_expname==" ")
	printlog("WARNING - QFLAT: key_exptime not set, cannot check exposure \
	    times.",l_logfile,yes)

    # Input info to logfile
    printlog("lampson = "//l_lampson,l_logfile,l_verbose)
    if (l_fl_lampsoff)
        printlog("lampsoff = "//l_lampsoff,l_logfile,l_verbose)
    if (l_fl_dark)
    	printlog("darks = "//l_darks,l_logfile,l_verbose)
    if (l_fl_bpm)
    	printlog("bpmfile = "//l_bpmfile,l_logfile,l_verbose)
    if (l_fl_flat)
    	printlog("flatfile = "//l_flatfile,l_logfile,l_verbose)
    printlog("thresh_flo = "//l_thresh_flo,l_logfile,verbose=l_verbose)
    printlog("thresh_fup = "//l_thresh_fup,l_logfile,verbose=l_verbose)
    if (l_fl_dark)
    	printlog("thresh_dlo = "//l_thresh_dlo,l_logfile,verb=l_verbose)
    if( l_fl_dark)
    	printlog("thresh_dup = "//l_thresh_dup,l_logfile,verb=l_verbose)
    printlog("normstat = "//l_normstat,l_logfile,l_verbose)
    if (l_fl_rmstars) {
	printlog("Removing stars from input images.",l_logfile,l_verbose)
	printlog("fwhmpsf = "//l_fwhmpsf,l_logfile,l_verbose)
    }

    #--------------------------------------------------------------------------
    # images ok, so combine, lamps on, lamps off, then darks
    n=0
    for (i = 1; i<=3; i+=1) {       
	if (infiles[i] != "" && infiles[i] != " ") {
            files(infiles[i],sort-, > filelist)

            l_sgain=l_key_gain
            l_sron=l_key_ron
            k=0
            flist=filelist
            while (fscan(flist, img) != EOF) {
        	n=n+1
        	k=k+1
        	# check gain and readnoise values 
        	if (l_reject=="ccdclip") {
                    if (l_key_gain != "" || l_key_gain != " ") {
                	if (l_sgain != str(l_gain)) {
                            imgets(img,l_key_gain, >& "dev$null")
                            if (imgets.value == "0" || imgets.value == "") {
                        	printlog("WARNING - QFLAT: keyword "//
				    l_key_gain//"not found in %s. ",
				    l_logfile,yes)
                        	printlog("Using gain = "//l_gain,l_logfile,yes)
                        	l_sgain=str(l_gain)
                            }
                	}
                    } else
                	l_sgain=str(l_gain)
			
                    if (l_key_ron != "" || l_key_ron != " ") {
                	if (l_sron != str(l_ron)) {
                            imgets(img,l_key_ron, >& "dev$null")
                            if (imgets.value == "0" || imgets.value == "") {
                        	printlog("WARNING - QFLAT: keyword"//\
				    l_key_ron//"not found in "//img,
				    l_logfile,yes) 
                        	printlog("Using ron = "//l_ron,l_logfile,yes)
                        	l_sron=str(l_ron)
                            }
                	}
                    } else
                	l_sron=str(l_ron)
        	}
            }
            flist=""

	    #combine images
            #special cases for low numbers of images
            if (k == 1) {
        	printlog("WARNING - QFLAT: only one image.",l_logfile,yes)
        	combout[i]=infiles[i]
            } else {
        	if (k <= 5){
                    if (k == 2){
                	printlog("WARNING - QFLAT: only two images to combine, \
			    turning off rejection.",l_logfile,yes)
                	l_reject="none"
                    } else
                	printlog("WARNING - QFLAT: five or less images to \
			    combine.",l_logfile,yes)
        	}
        	if (i==1 && l_fl_rmstars==yes) {
                    qsky("@"//filelist,combout[i],outtitle="default",
                        combtype=l_combine,rejtype=l_reject,logfile=l_logfile,
                        nlow=0,nhigh=1,lsigma=l_lsigma,hsigma=l_hsigma,
			threshold=4.5,fwhmpsf=l_fwhmpsf,datamax=l_thresh_fup,
			key_ron=l_key_ron,key_gain=l_key_gain,ron=l_ron,
			gain=l_gain,key_filter="FILTER",key_airmass="AIRMASS",
			masksuffix="msk",maskfactor=1.,
			fl_gemseeing=l_fl_gemseeing,fl_keepmasks=l_fl_keepmasks,
                        verbose=l_verbose)
                    combdq[i]=combout[i]//"msk"
        	} else {
        	    if (l_combine=="default")
		        l_combine="average"

        	    imcombine("@"//filelist,combout[i],headers="",bpmasks="",
		        rejmasks="",nrejmasks=combdq[i],expmasks="",
			sigmas=combsig[i],logfile="",combine=l_combine,
			reject=l_reject,project=no,outtype="real",outlimits="",
			offsets="none",masktype="goodvalue",maskvalue=0.,
			blank=0.,scale=l_scale,zero=l_zero,weight=l_weight,
		        statsec=l_statsec,expname=l_expname,
			lthreshold=l_lthreshold,hthreshold=l_hthreshold,
			nlow=l_nlow,nhigh=l_nhigh,nkeep=l_nkeep,mclip=l_mclip,
			lsigma=l_lsigma,hsigma=l_hsigma,rdnoise=l_sron,
			gain=l_sgain,snoise=l_snoise,sigscale=l_sigscale,
			pclip=l_pclip,grow=l_grow, >& "dev$null")

		    # bad pixels have dqsum[i]=k
		    imcalc(combdq[i]//".pl",
		        combdq[i]//".pl","if (im1 == "//k//") then 1 else 0",
			verbose-,pixtype="int")
        	}

        	# exposure time
        	imgets(combout[i],l_expname, >& "dev$null")
        	etime[i]=real(imgets.value)

        	#edit header
        	gemhedit (combout[i], "BPM", "", "", delete=yes)
            }
            delete(filelist,verify-)
	}
    }  #end of for loop

    #-------------------------------------------------------------------------- 
    # Subtract and normalize flats, combine bad pixel masks

    # subtract lampsoff from lampson
    if (l_fl_lampsoff) {
        if (abs(etime[1]-etime[2]) > 0.1) {
	    printlog("WARNING - QFLAT:  Exposure times for lamps-on and \
	        lamps-off images are",l_logfile,yes)
	    printlog("                  different. Thermal emission will not \
	        be subtracted",l_logfile,yes)
	    printlog("                  correctly. ",l_logfile,yes)
       }
       imcalc(combout[1]//","//combout[2],dsub,"im1-im2",verbose-,
           pixtype="real")
    } else
	imcopy(combout[1],dsub,verbose-)

    # Make flat
    imstat(dsub//statsec,fields="mode,stddev",lower=l_thresh_flo,
	upper=l_thresh_fup,nclip=0,lsigma=INDEF,usigma=INDEF,binwidth=0.1,
     	  format-,cache-) | scan(norm,stddev)
    l_flat_lo=norm-5*stddev
    l_flat_up=norm+5*stddev
    imstat(dsub//statsec,fields=l_normstat,lower=l_flat_lo,upper=l_flat_up,
        nclip=0,lsigma=INDEF,usigma=INDEF,binwidth=0.1,format-,cache-) | \
	scan(norm)
    l_upper=l_thresh_fup
    l_lower=l_thresh_flo
    if (l_lower < 0.0)
	l_lower=1.0

    # Normalize and put the right info in the headers
    if (l_fl_flat) {
	imarith(dsub,"/",norm,l_flatfile,divzero=0.,verbose-)
    gemdate ()
	gemhedit(l_flatfile,"GEM-TLM",gemdate.outdate,"UT Last modification with GEMINI")
	gemhedit(l_flatfile,"QFLAT",gemdate.outdate,"UT Time stamp for qflat")
	gemhedit(l_flatfile,"QFLTNORM",norm,"Normalization for qflat")
	if (l_fl_rmstars) {
	    gemhedit(l_flatfile,"QFLTRMST",l_fl_rmstars,"Stars removed by qsky")
	    gemhedit (l_flatfile, "QSKYMASK", "", "", delete=yes)
	} else {
	    gemhedit(l_flatfile,"QFLTCOMB",l_combine,
	        "Calculation method used by qflat")
	    if (l_zero!="none")
		gemhedit(l_flatfile,"QFLTZERO",l_zero,
		    "Statistic used to compute additive offsets") 
	    if (l_scale!="none")
		gemhedit(l_flatfile,"QFLTSCAL",l_scale,
        	    "Statistic used to compute scale factors")
	    if (l_weight!="none")
		gemhedit(l_flatfile,"QFLTWEIG",l_weight,
        	    "Statistic used to compute relative weights") 
	    gemhedit(l_flatfile,"QFLTSTAT",l_statsec,
		"Statistics section used by qflat")
	    if (l_reject!="none") {
		gemhedit(l_flatfile,"QFLTREJE",l_reject,
        	    "Rejection algorithm used by qflat")
		gemhedit(l_flatfile,"QFLTLTHR",l_lthreshold,
        	    "Lower threshold before combining")
		gemhedit(l_flatfile,"QFLTHTHR",l_hthreshold,
        	    "Upper threshold before combining")
		gemhedit(l_flatfile,"QFLTGROW",l_grow,
        	    "Radius for additional pixel rejection")
	    }
	    if (l_reject=="minmax") {
		gemhedit(l_flatfile,"QFLTNLOW",l_nlow,
		    "Low pixels rejected (minmax)")
		gemhedit(l_flatfile,"QFLTNHIG",l_nhigh,
		    "High pixels rejected (minmax)")
	    }
	    if (l_reject=="sigclip" || l_reject=="avsigclip" || \
	        l_reject=="ccdclip" || l_reject=="crreject" || \
		l_reject=="pclip") {
		gemhedit(l_flatfile,"QFLTLSIG",l_lsigma,
		    "Lower sigma for rejection")
		gemhedit(l_flatfile,"QFLTHSIG",l_hsigma,
		    "Upper sigma for rejection")
		gemhedit(l_flatfile,"QFLTNKEE",l_nkeep,
        	    "Min(max) number of pixels to keep(reject)")
		gemhedit(l_flatfile,"QFLTMCLI",l_mclip,
        	    "Use median in clipping algorithms")
	    } 
	    if (l_reject=="sigclip" || l_reject=="avsigclip" || \
	        l_reject=="ccdclip" || l_reject=="crreject")
		gemhedit(l_flatfile,"QFLTSIGS",l_sigscale,
        	    "Tolerance for sigma clipping scaling correction")
	    if (l_reject=="pclip")
		gemhedit(l_flatfile,"QFLTPCLI",l_pclip,
		    "Percentile clipping factor used by pclip")
	    if (l_reject=="ccdclip")
		gemhedit(l_flatfile,"QFLTSNOI",l_snoise,
        	    "Sensitivity noise (e) used by ccdclip")
	} # end else

	if (l_flattitle != "" && l_flattitle != " ")
	    gemhedit (l_flatfile, "i_title", l_flattitle, "", delete-)

	delete(filelist,verify-, >& "dev$null")
	files(infiles[1],sort-, > filelist)     
	i=1
	flist=filelist
	while (fscan(flist,img) != EOF) {
	    gemhedit(l_flatfile,"QFLTON"//str(i),img,
	        "Input sky or dome flat with lamps on")
	    i+=1
	}
	delete(filelist,verify-, >& "dev$null")
	if (l_fl_lampsoff) {
	    files(infiles[2],sort-, > filelist)
	    i=1
	    flist=filelist
	    while (fscan(flist,img) != EOF) {
		gemhedit(l_flatfile,"QFLTOF"//str(i),img,
		    "Input dome flat with lamps off")
		i+=1
	    }
	}
    } # end if(l_fl_flat)

    #--------------------------------------------------------------------------
    # Make bad pixel mask
    if (l_fl_bpm) {
	imcalc(dsub,mask[1],
	    "if (im1 < "//l_lower//") || (im1 > "//(l_upper)//") then 1 else 0",
	    pixtype="int",verbose-)

        if (l_fl_lampsoff)
	    imcalc(combdq[1]//".pl,"//combdq[2]//".pl",combdq[1]//".pl",
	        "if((im1==1) || (im2==1)) then 1 else 0",pixtype="int",verbose-)
        if (l_fl_dark) {
	    imcalc(combdq[1]//".pl,"//combdq[3]//".pl",combdq[1]//".pl",
	        "if((im1==1) || (im2==1)) then 1 else 0",pixtype="int",verbose-)
	    imcalc(combout[3],mask[2],
	        "if((im1<"//l_thresh_dlo//") || (im1>"//l_thresh_dup//")) then 1 else 0",
	        pixtype="int",verbose-)
	    imcalc(mask[1]//","//mask[2]//","//combdq[1]//".pl",
	        l_bpmfile//".pl",
		"if((im1==1) || (im2==1) || (im3==1)) then 1 else 0",
	        pixtype="int",verbose-)
        } else
	    imcalc(mask[1]//","//combdq[1]//".pl",l_bpmfile//".pl",
	        "if((im1==1) || (im2==1)) then 1 else 0",pixtype="int",verbose-)
        l_bpmfile=l_bpmfile+".pl"
        gemdate ()
        gemhedit(l_bpmfile,"GEM-TLM",gemdate.outdate,"UT Last modification with GEMINI")
        gemhedit(l_bpmfile,"QFLAT",gemdate.outdate,"UT Time stamp for qflat")
        gemhedit(l_bpmfile,"QFLTTFLO",l_thresh_flo,
            "Lower flat cutoff for bad pixels")
        gemhedit(l_bpmfile,"QFLTTFUP",l_thresh_fup,
            "Upper flat cutoff for bad pixels")
        if (l_fl_dark) {
	    gemhedit(l_bpmfile,"QFLTTDLO",l_thresh_dlo,
                "Lower dark cutoff for bad pixels")
	    gemhedit(l_bpmfile,"QFLTTDUP",l_thresh_dup,
                "Upper dark cutoff for bad pixels")
	}
	if (l_bpmtitle != "" && l_bpmtitle != " ")
	    gemhedit (l_bpmfile, "i_title", l_bpmtitle, "", delete-)
	if (l_fl_rmstars)
	    gemhedit (l_flatfile, "QSKYMASK", "", "", delete=yes)
	delete(filelist,verify-, >& "dev$null")
	files(infiles[1],sort-, > filelist)     
	i=1
	flist=filelist
	while (fscan(flist,img) != EOF) {
	    gemhedit(l_bpmfile,"QFLTON"//str(i),img,
	        "Input sky or dome flat with lamps on")
	    i+=1
	}
	delete(filelist,verify-, >& "dev$null")
	if (l_fl_lampsoff) {
	    files(infiles[2],sort-, > filelist)
	    i=1
	    flist=filelist
	    while (fscan(flist,img) != EOF) {
		gemhedit(l_bpmfile,"QFLTOF"//str(i),img,
		    "Input dome flat with lamps off")
		i+=1
	    }
	}
	if (l_fl_dark) {
	    delete(filelist,verify-, >& "dev$null")
	    files(infiles[3],sort-, > filelist)     
	    i=1
	    flist=filelist
	    while (fscan(flist,img) != EOF) {
		gemhedit(l_bpmfile,"QFLTDK"//str(i),img,"Input dark image")
		i+=1
	    }
	}
    }

    # clean up
    goto clean

error:
    status=1
    goto clean

clean:
    delete(filelist,verify-, >& "dev$null")
    for (i=1; i<=3; i+=1) {
	imdelete(combout[i]//","//combsig[i]//","//combdq[i],verify-,
	    >& "dev$null")
    }
    imdelete(dsub//","//mask[1]//","//mask[2]//",*msk.pl",verify-,
        >& "dev$null")
    flist=""
    # close log file
    printlog(" ",l_logfile,l_verbose)
    if (status==0) printlog("QFLAT exit status:  good.",l_logfile,l_verbose)
    printlog("--------------------------------------------------------------\
        ------------",l_logfile,l_verbose)

end
