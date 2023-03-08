# Copyright(c) 2002-2015 Association of Universities for Research in Astronomy, Inc

procedure gbpm(flatlong,flatshort,biaslist,bpmfile)

# Takes GMOS raw flats and bias (OBSMODE=IMAGE), 3 or 6 SCI extensions 
# and calculates the Bad Pixel Mask flagging CCD dead/hot/non-linear pixels 
# Sets bad pixels to 1 bit
#
# Version  Feb 28, 2002  JJ v1.3 release, not tested on spec flats
#          Sept 20, 2002 IJ v1.4 release
#          Aug 25, 2003  KL IRAF2.12 - new parameters
#                             hedit: addonly
#                             imstat: nclip, lsigma, usigma, cache
#          Dec 9, 2003   IJ Reduction changed to agree with current use, 
#                           defaults changed, use giflat to get normalized flat
#                           fields, small bugfixes 

char    flatlong    {prompt="Input list of long exposure flats"}
char    flatshort   {prompt="Input list of short exposure flats"}
char    biaslist    {prompt="Input list of bias images"}
char    bpmfile     {prompt="Output Bad Pixel Mask filename"}
char    rawpath     {"",prompt="Path for input raw images"}
char    gp_outpref  {"g",prompt="GPREPARE: Prefix for output images"}
char    title       {"BPM",prompt="Title for Bad Pixel Mask"} 
char    sci_ext     {"SCI",prompt="Name of science extension"}
char    dq_ext      {"DQ",prompt="Name of data quality extension"}
char    gaindb      {"default",prompt="Database with gain data"}
char    key_ron     {"RDNOISE",prompt="Header keyword for readout noise"}
char    key_gain    {"GAIN",prompt="Header keyword for gain [e-/ADU]"}
real    ron         {3.5,prompt="Readout noise in electrons"}
real    gain        {2.2,prompt="Gain in e-/ADU"}
char    key_biassec {"BIASSEC",prompt="Header keyword for overscan image section"}
char    key_datasec {"DATASEC",prompt="Header keyword for data section (excludes the overscan)"}
bool    fl_trim     {yes,prompt="Trim overscan section"}
bool    fl_over     {no,prompt="Reduce images using overscan subtraction"}
bool    fl_inter    {no,prompt="Set cut levels interactively"}
real    bhigh_cut   {20.,prompt="Hot pixel (high bias) threshold [ADU] above the median bias level"}
real    fdead_cut   {0.85,prompt="Low-sensitivity pixel threshold from long flats"}
real    flow_cut    {0.95,prompt="Lower non-linear bad-pixel threshold"}
real    fhigh_cut   {1.05,prompt="Upper non-linear bad-pixel threshold"}
char    sat         {"default", prompt="Saturation level in raw images [ADU/pixel]"}
bool    inter       {no,prompt="Interactive overscan fitting"}
bool    median      {no,prompt="Use median instead of average in column bias"}
char    function    {"chebyshev",prompt="Overscan fitting function.",enum="legendre|chebyshev|spline1|spline3"}
int     order       {1,prompt="Order of overscan fitting function"}
real    low_reject  {3.,prompt="Low sigma rejection factor"}
real    high_reject {3.,prompt="High sigma rejection factor"}
int     niterate    {3,prompt="Number of rejection iterations"}
bool    fl_qecorr   {no, prompt="QE correct the input images?"}
char    qe_data     {"gmosQEfactors.dat", prompt="Data file that contains QE information."}
char    qe_datadir  {"gmos$data/", prompt="Directory containg QE data file."}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"} 
struct  *scanfile   {"",prompt="Internal use only"}

begin

    # Local variables:
    char    l_flatshort, l_flatlong, l_biaslist, l_bpmfile, l_title, l_logfile
    char    l_sci_ext, l_key_ron, l_key_gain, l_key_datasec, l_func 
    char    l_key_biassec, l_dq_ext, l_rawpath, l_gp_outpref, l_gaindb, l_sat
    bool    l_fl_inter, l_verbose, l_inter, l_med, l_fl_over
    bool    l_fl_keepgprep, l_fl_trim, l_fl_short, l_fl_bias
    int     l_order, l_niter
    real    satvalue
    real    l_ron, l_gain, l_flow_cut, l_fhigh_cut, l_bhigh_cut, l_fdead_cut
    real    l_lowr, l_highr
    struct  l_struct
    char    inlist[3], flist[3], tmplist, tmplist2, img, mskexp, tmpsci
    char    tmpratio, combout[3], mask[3,50], tmpbpm[50], fxexp, obsmode
    char    obsmode_old, temp, keyfound, l_qe_data, l_qe_datadir
    int     i, j, nsci[3], len, nbad, nfiles[3]
    real    upp, low, stddev, hotpix, deadpix, dum, test, l_median
    bool    l_fl_raw, l_fl_newbpm, l_fl_qecorr

    # Make temporary files
    tmplist = mktemp("tmplist")

    # Assign dummy value to tmp file variables
    # (just need to be given any value; the real value will be given later)
    for (i=1; i<=3; i+=1) {
        flist[i]="dummy"
    }

    # Read QE correction parameters
    l_fl_qecorr = fl_qecorr
    l_qe_data = qe_data
    l_qe_datadir = qe_datadir

    # Set local variables
    l_flatshort=flatshort; l_flatlong=flatlong; l_biaslist=biaslist
    l_bpmfile=bpmfile
    l_title=title; l_logfile=logfile; l_sci_ext=sci_ext; l_key_ron=key_ron
    l_key_gain=key_gain; l_key_datasec=key_datasec; l_func=function
    l_key_biassec=key_biassec; l_fl_inter=fl_inter; l_verbose=verbose 
    l_inter=inter; l_med=median; l_dq_ext=dq_ext
    l_fl_trim=fl_trim ; l_fl_over=fl_over
    l_order=order; l_niter=niterate
    l_ron=ron; l_gain=gain; l_flow_cut=flow_cut; l_fhigh_cut=fhigh_cut
    l_bhigh_cut=bhigh_cut; l_fdead_cut=fdead_cut; l_sat=sat; l_lowr=low_reject
    l_highr=high_reject
    l_rawpath=rawpath ; l_gp_outpref=gp_outpref ; l_gaindb=gaindb

    status = 0

    cache("imgets", "gimverify", "gemhedit", "gemdate")
    
    l_fl_short=yes
    l_fl_bias=yes
    l_fl_newbpm=no

    # check input
    print(l_flatshort) | scan(l_flatshort)
    if (l_flatshort=="")
        l_fl_short=no
    print(l_biaslist) | scan(l_biaslist)
    if (l_biaslist=="")
        l_fl_bias=no

    #-----------------------------------------------------------------------
    # Check the logfile:
    if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GBPM: both gbpm.logfile and gmos.logfile \
                are empty.", l_logfile,verbose+)
            printlog ("                Using default file gmos.log.",
                l_logfile,verbose+)
        }
    }

    date | scan(l_struct)
    printlog ("-----------------------------------------------------------\
        -----------------", l_logfile, l_verbose)
    printlog ("GBPM -- "//l_struct, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)
    printlog ("Long flats list = "//l_flatlong, l_logfile, l_verbose)
    if (l_fl_short)
        printlog ("Short flats list = "//l_flatshort, l_logfile, l_verbose)
    if (l_fl_bias)
        printlog ("Bias list = "//l_biaslist, l_logfile, l_verbose)
    printlog ("Output BPM = "//l_bpmfile, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    #------------------------------------------------------------------------
    # Tests the SCI extension name:
    if (l_sci_ext=="" || stridx(" ",l_sci_ext)>0) {
        printlog ("ERROR - GBPM: Science extension name sci_ext is missing",
            l_logfile,verbose+)
        goto crash
    }
    if (l_dq_ext=="" || stridx(" ",l_dq_ext)>0) {
        printlog ("ERROR - GBPM: Data quality extension name dq_ext is \
            missing", l_logfile,verbose+)
        goto crash
    }

    # Tests the RDNOISE and GAIN keywords
    if (l_key_ron == "" || stridx(" ",l_key_ron)>0) {
        printlog ("ERROR - GBPM: Readout noise keyword parameter key_ron \
            is missing.", l_logfile, verbose+)
        goto crash
    }
    if (l_key_gain == "" || stridx(" ",l_key_gain)>0) {
        printlog ("ERROR - GBPM: Gain keyword parameter key_gain is missing.",
            l_logfile, verbose+)
        goto crash
    } 

    # Tests the output BPM file:
    print(l_bpmfile) | scan(l_bpmfile)
    if (l_bpmfile == "" || stridx(" ",l_bpmfile)>0) {
        printlog ("ERROR - GBPM: The output filename bpmfile is not defined",
            l_logfile, verbose+)
        goto crash
    }
    if (imaccess(l_bpmfile)) {
        printlog ("ERROR - GBPM: Output file "//l_bpmfile//" already exists.",
            l_logfile, verbose+)
        goto crash
    }

    len = strlen(l_bpmfile)
    if (substr(l_bpmfile,len-2,len)==".pl")
        l_bpmfile = substr(l_bpmfile,1,len-2)
    if (substr(l_bpmfile,len-4,len)!=".fits")
        l_bpmfile = l_bpmfile//".fits"

    # Define the lists:
    if(l_fl_short)
        inlist[3] = l_flatshort
    inlist[2] = l_flatlong
    if(l_fl_bias)
        inlist[1] = l_biaslist
    combout[1] = mktemp("tmpcmb1")
    combout[2] = mktemp("tmpcmb2")
    combout[3] = mktemp("tmpcmb3")

    #-------------------------------------------------------------------------
    # Check input files

    print(l_flatlong) | scan(l_flatlong)
    if (l_flatlong=="" || l_flatlong==" ") {
        printlog ("ERROR - GBPM: Input file list for long flats is missing",
            l_logfile, verbose+)
        goto crash
    }
    if (!l_fl_short) {
        printlog ("ERROR - GBPM: Input file list for short flats is missing",
            l_logfile, verbose+)
        goto crash
    }
    if (!l_fl_bias) {
        printlog ("ERROR - GBPM: Input file list for bias images is missing",
            l_logfile, verbose+)
        goto crash
    }

    nbad=0
    for (i=1; i<=3; i+=1) {
        if ((i==2) || ((i==1) && (l_fl_bias)) || ((i==3) && (l_fl_short))) {
            if(substr(inlist[i],1,1) == "@") {
                len = strlen(inlist[i])
                if (!access(substr(inlist[i],2,len))) {
                    printlog ("ERROR - GBPM: "//substr(inlist[i],2,len)//\
                        " does not exist.", l_logfile, verbose+)
                    goto crash
                }
            }
            tmplist2 = mktemp("tmplist2")
            flist[i] = tmplist2
            files (inlist[i],sort-, > tmplist)
            scanfile = tmplist
            nfiles[i] = 0
            while (fscan(scanfile,img) != EOF) {
                l_fl_raw = no
                gimverify (img)
                img = gimverify.outname
                if (gimverify.status>=1 && l_rawpath!="") {
                    gimverify (l_rawpath//img)
                    l_fl_raw=yes
                }
                if(gimverify.status==1) {
                    printlog("ERROR - GBPM: File "//img//" not found.",
                        l_logfile,verbose+)
                    nbad+=1
                } else if (gimverify.status>1) {
                    printlog ("ERROR - GBPM: File "//img//" not a MEF FITS \
                        image.", l_logfile, verbose+)
                    nbad+=1
                } else {
                    keyfound = ""
                    if (l_fl_raw)
                        hselect (l_rawpath//img//"[0]", "*PREPAR*", yes) | \
                            scan (keyfound)
                    else
                        hselect (img//"[0]", "*PREPAR*", yes) | scan (keyfound)

                    if (l_fl_raw && keyfound != "") {
                        printlog("ERROR - GBPM: Input image "//img//" is \
                            in raw data directory "//l_rawpath, l_logfile, \
                            verbose+)
                        printlog("              and is gprepared", l_logfile, \
                            verbose+)
                        goto crash
                    }
                    if (keyfound == "") {
                        # Use already gprepared image if it exists
                        if (access(l_gp_outpref//img//".fits")) {
                            hselect (l_gp_outpref//img//"[0]", "*PREPAR*", \
                                yes) | scan (keyfound)
                            if (keyfound == "") {
                                printlog ("ERROR - GBPM: gprepare outpref \
                                    parameter points to existing un-gprepared \
                                    file", l_logfile, verbose+)
                                goto crash
                            } else {
                                printlog("Using already gprepared image "//\
                                    l_gp_outpref//img, l_logfile, l_verbose)
                            }
                        } else {
                            gprepare (img, rawpath=l_rawpath, outimages="",
                                outpref=l_gp_outpref, fl_addmdf=no,
                                logfile=l_logfile, sci_ext=l_sci_ext,
                                gaindb=l_gaindb, key_ron=l_key_ron,
                                key_gain=l_key_gain, verbose=l_verbose)

                            if (gprepare.status!=0) {
                                printlog("ERROR - GBPM: gprepare returned \
                                    with error", l_logfile, verbose+)
                                goto crash
                            }
                        }
                        img = l_gp_outpref//img
                    }
                    nfiles[i] = nfiles[i]+1
                    print (img, >> tmplist2)
                    # get first header from long flats for PHU of output 
                    # image, open dummy image
                    if (nfiles[i]==1 && i==2 && nbad==0) {
                        imcopy (img//"[0]", l_bpmfile, verbose-)
                        gemhedit (l_bpmfile, "BPMMASK", "", "", delete+)
                        l_fl_newbpm=yes
                    }
                }
            } # end while
            scanfile = ""
            
            if (nfiles[i]==0) {
                printlog ("ERROR - GBPM: No input images meed wildcard \
                    criteria.", l_logfile, verbose+)
                goto crash
            }
            delete (tmplist, verify-, >& "dev$null") 
        } # end if(i==2) etc.
    } # end for loop

    # Exit if problems found
    if (nbad > 0) {
        printlog ("ERROR - GBPM: "//nbad//" image(s) either do not exist, \
            are not MEF files, or", l_logfile, verbose+)
        printlog ("              have not been run through GPREPARE.",
            l_logfile, verbose+)
        goto crash
    }

    #--------------------------------------------------------------------------
    # Call gbias to trim and combine the bias - have to add a test to 
    # skip this part if a reduced bias is supplied
    if (l_fl_bias) {
        gbias ("@"//flist[1], combout[1], logfile=l_logfile, fl_over=l_fl_over,
            key_biassec=l_key_biassec, key_datasec=l_key_datasec,
            fl_trim=l_fl_trim, key_ron=l_key_ron, key_gain=l_key_gain,
            ron=l_ron, gain=l_gain, fl_vardq=no, sci_ext=l_sci_ext,
            var_ext="VAR", dq_ext=l_dq_ext, bpm="", sat=l_sat, 
            fl_inter=l_inter, median=l_med, function=l_func, order=l_order, 
            low_reject=l_lowr, high_reject=l_highr, niterate=l_niter, 
            combine="average", reject="minmax", masktype="goodvalue", 
            maskvalue=0., lthreshold=INDEF, hthreshold=INDEF, scale="none",
            zero="none", weight="none", statsec="[*,*]", key_exptime="EXPTIME",
            nlow=0, nhigh=1, nkeep=1, verbose=l_verbose, rawpath="")
        if (gbias.status!=0) {
            printlog ("ERROR - GBPM: Error in GBIAS prevents completion \
                of GBPM.", l_logfile, verbose+)
            goto crash
        }
    }

    # At this point all the images have already been gprepared, use OBSMODE
    obsmode_old = ""
    for (i=2; i<=3; i+=1) {
        if ((i==2) || ((i==3) && (l_fl_short))) {
            scanfile = flist[i]
            while (fscan (scanfile,img) != EOF) {
                obsmode = "UNKNOWN"
                hselect (img//"[0]", "OBSMODE", yes) | scan (obsmode)
                if (obsmode != "IMAGE") {
                    printlog ("ERROR - GBPM: Only imaging flats can be used \
                        for GBPM at this point.", l_logfile, verbose+)
                    goto crash
                }
                if (obsmode_old != "" && obsmode != obsmode_old) {
                    printlog ("ERROR - GBPM: The flatfield images are of \
                        mixed OBSMODES.", l_logfile, verbose+)
                    goto crash
                }
                obsmode_old = obsmode
            }
            scanfile = ""
            delete (tmplist, verify-, >& "dev$null") 

            # Call giflat or gsflat to trim and combine the flats
            giflat ("@"//flist[i], combout[i], normsec="default", fl_scale=yes,
                sctype="mean", statsec="default", key_ron=l_key_ron,
                key_gain=l_key_gain, fl_vardq=no, sci_ext=l_sci_ext,
                var_ext="VAR", dq_ext=l_dq_ext, bpm="", sat=l_sat,
                fl_inter=l_inter, median=l_med, function=l_func, 
                key_biassec=l_key_biassec, key_datasec=l_key_datasec,
                fl_trim=l_fl_trim, fl_bias=yes, bias=combout[1],
                logfile=l_logfile, fl_over=l_fl_over, order=l_order,
                low_reject=l_lowr, high_reject=l_highr, niterate=l_niter,
                combine="average", reject="avsigclip", lthreshold=INDEF,
                hthreshold=INDEF, nkeep=1, lsigma=3., hsigma=3., 
                verbose=l_verbose, rawpath="", gp_outpref=l_gp_outpref,
                gi_outpref="r", fl_qecorr=l_fl_qecorr, \
                qe_data=l_qe_data, qe_datadir=l_qe_datadir)

            if (giflat.status!=0) {
                printlog("ERROR - GBPM: Error in GIFLAT prevents completion \
                    of GBPM.", l_logfile, verbose+)
                goto crash
            }
        }
    } # end for(i=2,3)

    # All data must have the same number of SCI extensions
    for (i=1; i<=3;i+=1) {
        if ((i==2) || ((i==1) && (l_fl_bias)) || ((i==3) && (l_fl_short))) {
            imgets (combout[i]//"[0]", "NSCIEXT", >>& "dev$null")
            nsci[i] = int(imgets.value)
        } else 
            nsci[i]=-1
    }

    if ((nsci[1]!=-1 && nsci[1]!=nsci[2]) \
        || (nsci[3]!=-1 && nsci[3]!=nsci[2])) {
        
        printlog ("ERROR - GBPM: The number of SCI extensions in input \
            images is different.", l_logfile, verbose+)
        goto crash
    }

    #--------------------------------------------------------------
    # Build the BPM :

    # Get upper limit from bias - these would be the hot pixels:
    if (l_fl_bias) {
        mskexp = "((a > b) ? 1 : 0)"
        for (j=1; j<=nsci[2]; j+=1) {
            tmpsci = combout[1]//"["//l_sci_ext//","//j//"]"
            imstat (tmpsci, fields="midpt", format-, lower=INDEF, upper=INDEF,
                nclip=0., lsigma=3, usigma=3., binwidth=0.1, cache=no) | \
                scan(l_median)
            if (l_fl_inter) {
                imhistogram (tmpsci, z1=-100.+l_median, z2=100.+l_median,
                    binwidth=1, nbins=512, autoscale=no, top_closed=no,
                    hist_type="normal", listout=no, plot_type="box", logy+,
                    device="stdgraph")
                print (" ")
                print ("Bias histogram - there will be one for each SCI \
                    extension")
                print ("Set cursor position and hit any key to mark UPPER \
                    limit for good pixels")
                dum = fscan(gcur,test)
                hotpix = test
            } else {
                hotpix = l_bhigh_cut+l_median
            }
            
            # Set saturation level
            if (l_sat == "default") { 
                 gsat (combout[1], extension="["//l_sci_ext//","//j//"]", \
                     gaindb=l_gaindb, bias="static", \
                     pixstat="midpt", statsec="default", \
                     gainval="default", logfile=l_logfile, \
                     verbose=no)
                 if (gsat.status != 0) {
                     printlog ("ERROR - GBPM: GSAT returned a "//\
                         "non-zero status. Exiting", \
                         l_logfile, verbose+)
                     status = 1
                     goto crash
                 }
                 satvalue = real(gsat.saturation)
            } else {
                # option for user-defined value
                satvalue = real(l_sat)
            }

            if (hotpix > satvalue)
                hotpix = satvalue
            gemhedit (l_bpmfile//"[0]", "GBPMHOT"//j, hotpix,
                "Limit for hot (high-bias) pixels in ext "//j//" (ADU)", 
                delete-)
            printlog ("Bias for extension "//j//". Limit for hot (high-bias) \
                pixels is "//hotpix, l_logfile, l_verbose)
            mask[1,j] = mktemp("tmpmsk1")
            imexpr (mskexp, mask[1,j], tmpsci, hotpix, dims="auto",
                intype="auto", outtype="short", refim="auto", bwidth=0,
                btype="nearest", bpixval=0., rangecheck+, verbose-,
                exprdb="none", lastout="")
        }
    }

    # Get lower limit for long flats - these would be the dead pixels:
    if (l_fl_short) {
        mskexp = "((a < b) ? 1 : 0)"
        for (j=1; j<=nsci[2];j+=1) {
            tmpsci = combout[2]//"["//l_sci_ext//","//j//"]"
            if (l_fl_inter) {
                imhistogram (tmpsci, z1=0.2, z2=2.0, binwidth=0.01, nbins=512,
                    autoscale=no, top_closed=no, hist_type="normal", 
                    listout=no, plot_type="box", logy+, device="stdgraph")
                print (" ")
                print ("Long flat field (normalized) histogram - there will \
                    be one for each SCI extension")
                print ("Set cursor position and hit any key to mark LOWER \
                    limit for good pixels")
                dum = fscan(gcur,test)
                deadpix = test
            } else {
                deadpix = l_fdead_cut # flat field is already normalized
            }
            if (deadpix < 0.)
                deadpix = 0.
            gemhedit (l_bpmfile//"[0]", "GBPMLOW"//j, deadpix,
                "Limit for low-sensitivity pixels, ext "//j//" (ADU)", delete-)
            printf ("Long flat for extension %s. Limit for low-sensitivity \
                pixels is %9.2f\n", j, deadpix) | scan(l_struct)
            printlog (l_struct,l_logfile,l_verbose)
            mask[2,j] = mktemp("tmpmsk2")
            imexpr (mskexp, mask[2,j], tmpsci, deadpix, dims="auto",
                intype="auto", outtype="short", refim="auto", bwidth=0,
                btype="nearest", bpixval=0., rangecheck+, verbose-, 
                exprdb="none", lastout="")
        }

        # Get upper and lower limits for ratio short/long flats - 
        # these would be the non-linear pixels
        mskexp = "((a < b) || (a > c) ? 1 : 0)"
        for (j=1; j<=nsci[2]; j+=1) {
            tmpratio = mktemp("tmpratio")
            imarith (combout[3]//"["//l_sci_ext//","//j//"]","/",
                combout[2]//"["//l_sci_ext//","//j//"]", tmpratio, verbose-)
            if (l_fl_inter) {
                imhistogram (tmpratio, z1=0.2, z2=2.0, binwidth=0.01, 
                    nbins=512, autoscale=no, top_closed=no, 
                    hist_type="normal", listout=no, plot_type="box", logy+, 
                    device="stdgraph")
                print (" ")
                print ("Flat field ratio (normalized) histogram - there will \
                    be one for each SCI extension")
                print ("Set cursor position and hit any key to mark UPPER \
                    limit for good pixels")
                dum = fscan(gcur,test)
                upp = test
                print ("Set cursor position and hit any key to mark LOWER \
                    limit for good pixels")
                dum = fscan(gcur,test)
                low = test
            } else {
                upp = l_fhigh_cut # flats are normalized
                low = l_flow_cut
            }
            if (upp < low) {
                test = low
                low = upp
                upp = test
            }

            # Set the saturation limit again - according to amplifier and
            # DETTYPE 
            if (l_sat == "default") { 
                 gsat (combout[3], extension="["//l_sci_ext//","//j//"]", \
                    gaindb=l_gaindb, bias="static", \
                    pixstat="midpt", statsec="default", \
                    gainval="default", logfile=l_logfile, \
                    verbose=no)
                if (gsat.status != 0) {
                    printlog ("ERROR - GBPM: GSAT returned a "//\
                        "non-zero status. Exiting", \
                        l_logfile, verbose+)
                    status = 1
                    goto crash
                }
                satvalue = real(gsat.saturation)
            } else {
                # option for user-defined value
                satvalue = real(l_sat)
            }

            if (low < 0.0)
                low = 0.0
            if (upp > satvalue)
                upp = satvalue

            gemhedit (l_bpmfile//"[0]","GBPMNUP"//j, upp,
                "Upper limit for non-linear pixels, ext "//j, delete-)
            gemhedit (l_bpmfile//"[0]", "GBPMNLO"//j, low,
                "Lower limit for non-linear pixels, ext "//j, delete-)
            printlog ("Flat field ratio for extension "//j//". Upper and \
                lower non-linearity", l_logfile, l_verbose) 
            printlog ("thresholds: "//upp//" and "//low, l_logfile, l_verbose)
            mask[3,j] = mktemp("tmpmsk3")
            imexpr (mskexp, mask[3,j], tmpratio, low, upp, dims="auto",
                intype="auto", outtype="short", refim="auto", bwidth=0,
                btype="nearest", bpixval=0., rangecheck+, verbose-, 
                exprdb="none", lastout="")
            imdelete (tmpratio, verify-, >>& "dev$null")
        }
    } # end if(l_fl_short)

    # Delete the combined images:
    imdelete (combout[1]//","//combout[2]//","//combout[3], verify-, 
        >>& "dev$null")

    #------------------------------------------------------------------------
    #Combine the masks:
    if (l_fl_bias && l_fl_short) {
        mskexp = "((a == 1) || (b ==1) || (c == 1) ? 1 : 0)"
        for (j=1; j<=nsci[2]; j+=1) {
            tmpbpm[j] = mktemp("tmpbpm")
            imexpr (mskexp, tmpbpm[j], mask[1,j], mask[2,j], mask[3,j],
                dims="auto", intype="auto", outtype="short", refim="auto",
                bwidth=0, btype="nearest", bpixval=0., rangecheck+, verbose-,
                exprdb="none", lastout="")
            imdelete (mask[1,j]//","//mask[2,j]//","//mask[3,j], verify-,
                >>& "dev$null")
        }
    } else if (l_fl_bias) {
        mskexp = "((a == 1) || (b ==1) ? 1 : 0)"
        for (j=1; j<=nsci[2]; j+=1) {
            tmpbpm[j] = mktemp("tmpbpm")
            imexpr (mskexp, tmpbpm[j], mask[1,j], mask[2,j], dims="auto",
                intype="auto", outtype="short", refim="auto", bwidth=0,
                btype="nearest", bpixval=0., rangecheck+, verbose-,
                exprdb="none", lastout="")
            imdelete (mask[1,j]//","//mask[2,j], verify-, >>& "dev$null")
        }
    } else if (l_fl_short) {
        mskexp = "((a == 1) || (b ==1) ? 1 : 0)"
        for (j=1; j<=nsci[2]; j+=1) {
            tmpbpm[j] = mktemp("tmpbpm")
            imexpr (mskexp, tmpbpm[j], mask[2,j], mask[3,j], dims="auto",
                intype="auto", outtype="short", refim="auto", bwidth=0,
                btype="nearest", bpixval=0., rangecheck+, verbose-,
                exprdb="none", lastout="")
            imdelete (mask[2,j]//","//mask[3,j], verify-, >>& "dev$null")
        }
    } else {
        for (j=1; j<=nsci[2]; j+=1) {
            tmpbpm[j] = mktemp("tmpbpm")
            imcopy (mask[2,j], tmpbpm[j], verbose-)
        }
    }

    # And pack the MEF:
    fxexp = tmpbpm[1]//".fits"
    j = 2
    if (nsci[2]>1) {
        while (j <= nsci[2]) {
            fxexp = fxexp//","//tmpbpm[j]//".fits"
            j+=1
        }
    }

    # Make the output mask
    fxinsert (fxexp, l_bpmfile//"[0]", group="", verbose-, >>& "dev$null")

    # Update the headers and clean up temporaries
    j = 1
    while (j<=nsci[2]) {
        # they all be DQ?
        gemhedit (l_bpmfile//"["//j//"]", "EXTNAME", l_dq_ext, "", delete-)
        gemhedit (l_bpmfile//"["//j//"]", "EXTVER", j, "", delete-)
        flpr
        imdelete (tmpbpm[j],verify-, >& "dev$null")
        j = j+1
    }

    # Update the PHU
    date | scan(l_struct)
    if (l_title != "" && l_title != " ")
        gemhedit (l_bpmfile//"[0]", "i_title", l_title, "Image title", delete-)
    else
        gemhedit (l_bpmfile//"[0]", "i_title",
            "GMOS Bad Pixel Mask - "//l_struct, "Image title", delete-)

    gemhedit (l_bpmfile//"[0]", "NEXTEND", nsci[2], "Number of extensions", \
        delete-)
    
    gemdate ()
    gemhedit (l_bpmfile//"[0]", "GBPM", gemdate.outdate, 
        "UT Time stamp for GBPM", delete-)
    gemhedit (l_bpmfile//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)

    if (obsmode != "")
        gemhedit (l_bpmfile//"[0]", "OBSMODE", obsmode,
            "Observing mode (IMAGE|IFU|MOS|LONGSLIT)", delete-)
    goto clean

    # Exit with error
crash:
    status=1
    if (l_fl_newbpm)
        delete (l_bpmfile, verify-, >>& "dev$null")
    goto clean

    # Clean up and exit
clean:
    delete (flist[1]//","//flist[2]//","//flist[3], verify-, >& "dev$null")
    printlog ("", l_logfile, l_verbose)
    if (status == 0) 
        printlog ("GBPM exit status: good.", l_logfile, l_verbose)
    else
        printlog ("GBPM exit status: error.", l_logfile, l_verbose)
    printlog ("------------------------------------------------------------\
        ----------------", l_logfile, l_verbose)

end
