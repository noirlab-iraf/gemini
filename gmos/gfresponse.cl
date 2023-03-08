# Copyright(c) 2002-2012 Association of Universities for Research in Astronomy, Inc.

procedure gfresponse (inimage,outimage)

# Derive the relative fiber response for the GMOS IFU
# inimage must be gfextracted
#
# Version   Sept 20, 2002 BM v1.4 release
#           Aug 25, 2003  KL IRAF2.12 - new/modified parameters
#                              hedit: addonly
#                              imcombine: headers,bpmasks,expmasks,outlimits
#                                 rejmask->rejmasks, plfile->nrejmasks
#                              imstat: nclip,lsigma,usigma,cache
#           Feb 29, 2004  BM Update for GMOS-S N&S
#           Mar 19, 2004  BM Fix ungraceful exit
#           Mar 25, 2004  BM Fix handling of blue slit
#           Apr 07, 2004  BM Fix problem in scaling of slits, put mktemps
#                            in the right place

string inimage   {prompt="Input file"}
string outimage  {prompt="Ouput fiber response file"}
string title     {"",prompt="Title for output science extension"}
string skyimage  {"",prompt="Sky response for illumination correction"}
bool   fl_inter  {no,prompt="Fit interactively?"}
bool   fl_fit    {no,prompt="Smooth fit to final flat?"}
string function  {"chebyshev",enum="spline1|spline3|legendre|chebyshev",prompt="Function for fit"}
int    order     {7,min=1,prompt="Order of 1D fit"}
string sample    {"*",prompt="Sample points to use in fit"}
string sci_ext   {"SCI",prompt="Name of science extension"}
string var_ext   {"VAR",prompt="Name of variance extension"}
string dq_ext    {"DQ",prompt="Name of data quality extension"}
string logfile   {"",prompt="Logfile name"}
bool   verbose   {yes,prompt="Verbose output?"}
int    status    {0,prompt="Exit status (0=good)"}

begin

    string l_inimage,l_outimage,l_title,l_logfile,l_skyimage,l_function
    string improj,imnorm,flatsky,pflatsky,stack,transpose,projfit,mdf
    string l_sci_ext,l_var_ext,l_dq_ext,l_sample, l_key_qecorrim
    real   pmean,mean[2]
    int    i,nx,nr,nslit,j1,id1,j,l_order
    bool   l_verbose,l_fl_inter,l_fl_fit, in_qestate, sky_qestate
    struct sdate

    # Query parameters
    l_inimage=inimage ; l_outimage=outimage ; l_title=title
    l_skyimage=skyimage; l_logfile=logfile; l_verbose=verbose
    l_fl_inter=fl_inter
    l_function=function ; l_order=order ; l_sample=sample
    l_fl_fit=fl_fit
    l_sci_ext=sci_ext ; l_var_ext=var_ext ; l_dq_ext=dq_ext

    status=0

    # Keep parameters from changing by outside world
    cache("imgets", "fparse","gemlogname","gimverify","gemdate")

    # Define temporary files
    mdf=mktemp("tmpmdf")
    # Define these to "" for now
    improj=""
    imnorm=""
    flatsky=""
    pflatsky=""
    stack=""
    transpose=""
    projfit=""

    # Set QE correction image keyword
    l_key_qecorrim = "QECORRIM"

    # Start logging to file
    if (l_logfile == "STDOUT") {
        l_logfile = ""
        l_verbose=yes
    }
    # Test the logfile:
    gemlogname(logpar=l_logfile,package="gmos")
    if (gemlogname.status != 0)
        goto error
    l_logfile=gemlogname.logname

    date | scan(sdate)
    printlog("----------------------------------------------------------------\
        ----------------", l_logfile,verbose=l_verbose)
    printlog("GFRESPONSE -- "//sdate,l_logfile,verbose=l_verbose)
    printlog(" ",l_logfile,verbose=l_verbose)
    printlog("inimage  = "//l_inimage,l_logfile,verbose=l_verbose)
    printlog("outimage = "//l_outimage,l_logfile,verbose=l_verbose)
    printlog("skyimage = "//l_skyimage,l_logfile,verbose=l_verbose)
    printlog("title    = "//l_title,l_logfile,verbose=l_verbose)
    printlog("order    = "//l_order,l_logfile,verbose=l_verbose)
    printlog("fl_fit   = "//l_fl_fit,l_logfile,verbose=l_verbose)
    printlog(" ",l_logfile,verbose=l_verbose)

    #check that there are input files
    if (l_inimage == "" || l_inimage == " "){
        printlog("ERROR - GFRESPONSE: input files not specified", l_logfile, \
            verbose+)
        goto error
    }

    # check existence of list file
    if (substr(l_inimage,1,1) == "@") {
        printlog("ERROR - GFRESPONSE: lists are currently not supported",
            l_logfile, verbose+)
        goto error
    }

    # check existence of input file
    gimverify(l_inimage)
    if (gimverify.status != 0) {
      printlog("ERROR - GFRESPONSE: "//l_inimage//" does not exist or is not \
          a MEF", l_logfile, verbose+)
      goto error
    }
    l_inimage=gimverify.outname//".fits"

    #check that an output file is given
    if (l_outimage == "" || l_outimage == " "){
        printlog("ERROR - GFRESPONSE: output files not specified", l_logfile, \
            verbose+)
        goto error
    }

    #if not given, put a .fits extension on output file name
    fparse(l_outimage,verbose-)
    if (fparse.extension == "")
        l_outimage=l_outimage//".fits"

    # Check that the output file does not already exist. If so, exit.
    if (imaccess(l_outimage)) {
      printlog("ERROR - GFRESPONSE: Output file "//l_outimage//" already "//\
          "exists.", l_logfile, verbose+)
      goto error
    }

    # Check whether input has been gsreduced

    # If gsreduced, check whether the extraction has taken place
    imgets(l_inimage//"[0]","GFEXTRAC", >& "dev$null")
    if (imgets.value == "0") {
        printlog("ERROR - GFRESPONSE: Input image must be processed by "//\
            "GFEXTRACT.", l_logfile, verbose+)
        goto error
    }

    # Read the QE correction state of the input image
    keypar (l_inimage//"[0]", l_key_qecorrim, silent+)
    if (keypar.found) {
        in_qestate = yes
    } else {
        in_qestate = no
    }

    # If skyimage given, verify it
    if (l_skyimage != "" && l_skyimage != " ") {
        gimverify(l_skyimage)
        if (gimverify.status != 0) {
            printlog("ERROR - GFRESPONSE: "//l_skyimage//" does not exist or \
                is not a MEF", l_logfile, verbose+)
            goto error
        }
        l_skyimage=gimverify.outname//".fits"

        # Check the QE correction state of the sky file
        keypar (l_skyimage//"[0]", l_key_qecorrim, silent+)
        if (keypar.found) {
            sky_qestate = yes
        } else {
            sky_qestate = no
        }
    }

    # which slit?
    tinfo(l_inimage//"[MDF]",ttout-)
    nr=tinfo.nrows
    tabpar(l_inimage//"[MDF]","NO",1)
    id1=int(tabpar.value)
    nslit=1
    if ((id1==1 && nr>750) || (id1==51 && nr>350)) {
        nslit=2
    }
    j1=1
    if (id1==751) {
        j1=2
        nslit=2
    }

    # Normalize extracted spectra, looping over slits
    j=0
    for (i=j1; i<=nslit; i+=1) {
        j=j+1
        if (imaccess(l_inimage//"["//l_sci_ext//","//j//"]")) {
            # Create the tmp img name for the images used only within this loop
            improj=mktemp("tmpimproj")
            imnorm=mktemp("tmpimnorm")
            flatsky=mktemp("tmpflatsky")
            pflatsky=mktemp("tmppflatsky")
            stack=mktemp("tmpstack")
            transpose=mktemp("tmptranspose")
            projfit=mktemp("tmpprojfit")

            if (l_verbose)
                printlog("Normalizing the fiber response for slit "//i,
                    l_logfile,l_verbose)
            # average by projecting to 1D
            imcombine(l_inimage//"["//l_sci_ext//","//j//"]",improj,headers="",
                bpmasks="",rejmasks="",nrejmasks="",expmasks="",sigmas="",
                logfile="dev$null",combine="average",reject="avsigclip",
                project+,outtype="real",outlimits="",offsets="none",
                masktype="none",maskvalue=0.,blank=0.,scale="",zero="none",
                weight="none",statsec="",expname="",lthresh=INDEF,
                hthresh=INDEF,nlow=1,nhigh=1,nkeep=0,mclip=yes,lsigma=5.,
                hsigma=5.,grow=0.)

            # fit the response function for the average spectrum
            fit1d(improj,projfit,"fit",inter=l_fl_inter,sample=l_sample,
                naverage=1,function=l_function, order=l_order,low_reject=3.,
                high_reject=3.,niterate=3,grow=0.)

            # divide by the average fit
            imarith(l_inimage//"["//l_sci_ext//","//j//"]","/",projfit,imnorm,
                title=l_title,verbose-,noact-)

            # correct by sky if available
            if (imaccess(l_skyimage//"["//l_sci_ext//","//j//"]")) {
                # Compare the QE correction states of the two files
                if (in_qestate != sky_qestate) {
                    printlog ("GFEXTRACT - ERROR: QE correction state "//\
                        "of the input file "//l_input//" ("//in_qestate//\
                        ") and the response file "//l_skyimage//" ("//\
                        sky_qestate//") mismatch.", \
                        l_logfile, l_verbose)
                    goto error
                }
                # divide the sky by the lamp flat
                imarith(l_skyimage//"["//l_sci_ext//","//j//"]","/",imnorm,\
                    flatsky,verbose-,noact-)
                imgets(flatsky,"i_naxis1")
                nx=int(imgets.value)

                # average the result in x
                improject(flatsky,pflatsky,1,average+,highcut=0.,lowcut=0.,
                    verbose-)

                # normalize
                imstat(pflatsky,field="mean",lower=INDEF,upper=INDEF,nclip=0,
                    lsigma=INDEF,usigma=INDEF,binwidth=0.1,format-,cache-) | \
                    scan(pmean)
                mean[i]=pmean
                imarith(pflatsky,"/",pmean,pflatsky)
                imdelete(flatsky, verify-, >& "dev$null")

                # Stack to add a second dimension
                imstack(pflatsky,stack,title="*",pixtype="real")

                # Transpose to get the image orientated the correct way
                imtranspose(stack,transpose)

                # Block replecate to the size of the input sky image
                blkrep (transpose, flatsky, nx, 1)

                # multiply the lamp flat by the smoothed flat sky
                imarith(imnorm,"*",flatsky,imnorm,verbose-,noact-)
                imdelete(flatsky//","//pflatsky//","//stack//","//transpose,
                    verify-)
            } else {
                # relative fluxes if no sky
                imstat(improj,field="mean",lower=INDEF,upper=INDEF,nclip=0,
                    lsigma=INDEF,usigma=INDEF,binwidth=0.1,format-,cache-) | \
                    scan(pmean)
                mean[i]=pmean
            }

            # Smooth final flat?
            if (l_fl_fit) {
                fit1d(imnorm,imnorm,"fit",interactive=l_fl_inter,axis=1,
                    sample="*",naverage=1,function=l_function, order=l_order,
                    low_reject=3.,high_reject=3.,niterate=3,grow=0.)
            }

            # write the normalized spectrum to the output file
            if (!imaccess(l_outimage)) {
                wmef(imnorm, l_outimage, extnames=l_sci_ext, phu=l_inimage, \
                    verbose-, >>& "dev$null")
                if (wmef.status > 0)
                    goto error
            } else {
                fxinsert(imnorm//".fits",l_outimage//"[1]","0",verbose-,
                    >& "dev$null")
                gemhedit(l_outimage//"["//j//"]", "EXTNAME", l_sci_ext, "",
                    delete-)
                gemhedit(l_outimage//"[0]", "NEXTEND", j, "", delete-)
            }
            gemhedit(l_outimage//"["//j//"]", "EXTVER", j, "", delete-)

            imdelete(improj//","//imnorm//","//projfit,verify-)
        }
    }

    # If two slits, correct for relative throughput of the slits
    # The reference is the red slit (slit 1)
    if (j1==1 && nslit==2) {
        printlog("",l_logfile,l_verbose)
        printlog("Correcting for the relative throughputs of the slits",
            l_logfile,l_verbose)
        printlog("Slit 2/Slit 1 = "//(mean[2]/mean[1]),l_logfile,l_verbose)
        flpr
        imarith(l_outimage//"["//l_sci_ext//",2]","*",(mean[2]/mean[1]),
            l_outimage//"["//l_sci_ext//",2,overwrite+]",verbose=no,noact=no)
    }

    # slice out the MDF, add to outimage
    tcopy(l_inimage//"[MDF]",mdf//".fits",verbose-)
    fxinsert(mdf//".fits",l_outimage//"["//j//"]",1,verbose-, >>& "dev$null")

    # final header update
    gemdate ()
    gemhedit (l_outimage//"[0]", "GFRESPON", gemdate.outdate,
        "UT Time stamp for GFRESPONSE", delete-)
    gemhedit(l_outimage//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)

    # clean up
    goto clean

error:
    status=1
    goto clean

clean:
    imdelete(improj//","//imnorm//","//projfit,verify-, >& "dev$null")
    imdelete(flatsky//","//mdf,verify-, >& "dev$null")
    # close log file
    printlog(" ",l_logfile,l_verbose)
    if(status==0)
        printlog("GFRESPONSE exit status: good",l_logfile,l_verbose)
    else
        printlog("GFRESPONSE exit status: error",l_logfile,l_verbose)
    printlog("----------------------------------------------------------------\
        ----------------", l_logfile,l_verbose)

end
