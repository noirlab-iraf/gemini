# Copyright(c) 2002-2011 Association of Universities for Research in Astronomy, Inc.

procedure gfapsum (inimages)

# Combine spectra from a GMOS IFU datacube
# 
# Version  Sept 20, 2002 BM   v1.4 release
#          Aug 25, 2003  KL  IRAF2.12 - new parameter, addonly, in hedit
#          Mar 19, 2004  BM  fix bug in ap selection, 
#                            allow use on non-GSTRANSFORMED 2-slit images
#          Apr 07, 2004  BM  default expr based on INSTRUME, relocate mktemp
#          Oct 28, 2005  BM  fix "rang e=" error in line 394

string  inimages    {prompt="Input images"}
string  outimages   {"",prompt="Output images"}
string  outprefix   {"a",prompt="Prefix for output images"}
string  apertures   {"",prompt="Aperture list"}
string  expr        {"default",prompt="Expression for aperture selection"}
string  combine     {"sum",enum="average|median|sum",prompt="scombine algorithm"}
string  reject      {"none",enum="none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip",prompt="Rejection algorithm"}
string  scale       {"none",prompt="Image scaling"}
string  zero        {"none",prompt="Image zeropoint offset"}
string  weight      {"none",prompt="Image weights"}
real    lthreshold  {INDEF,prompt="Lower threshold"}
real    hthreshold  {INDEF,prompt="Upper threshold"}
int     nlow        {1,min=0,prompt="minmax: Number of low pixels to reject"}
int     nhigh       {1,min=0,prompt="minmax: Number of high pixels to reject"}
int     nkeep       {0,min=0,prompt="Minimum to keep or maximum to reject"}  
bool    mclip       {yes,prompt="Use median in sigma clipping algorithms?"}
real    lsigma      {3.,prompt="Lower sigma clipping factor"}
real    hsigma      {3.,prompt="Upper sigma clipping factor"}
string  key_ron     {"RDNOISE",prompt="Keyword for readout noise in e-"}
string  key_gain    {"GAIN",prompt="Keyword for gain in electrons/ADU"}
string  snoise      {"0.0",prompt="Sensitivity noise (fraction), ccdclip and crreject"}
real    sigscale    {0.1,prompt="Tolerance for sigma clipping scaling correction"}
real    pclip       {-0.5,prompt="pclip: Percentile clipping parameter"}
real    grow        {0.0,prompt="Radius (pixels) for neighbor rejection"}
real    blank       {0.0,prompt="Value if there are no pixels"}
string  sci_ext     {"SCI",prompt="Name of science extension"}
string  var_ext     {"VAR",prompt="Name of variance extension"}
string  dq_ext      {"DQ",prompt="Name of data quality extension"}
bool    fl_inter    {yes,prompt="Select spectra interactively?"}
string  logfile     {"",prompt="Logfile name"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile1  {"",prompt="Internal use only"}
struct  *scanfile2  {"",prompt="Internal use only"}

begin

    # Local variables for input parameters
    string  l_inimages, l_outimages, l_prefix, l_logfile
    string  l_sci_ext, l_var_ext, l_dq_ext
    string  l_expr, l_aper
    bool    l_verbose, l_fl_inter, l_mclip
    string  l_combine, l_reject, l_scale, l_zero, l_weight
    string  l_key_ron, l_key_gain, l_snoise, l_gain, l_rdnoise
    real    l_lthreshold, l_hthreshold, l_lsigma, l_hsigma
    real    l_grow, l_sigscale, l_pclip, l_blank
    int     l_nlow, l_nhigh, l_nkeep

    # Other local variables
    string  inlist, outlist, temp1, temp2, temp3, mdf, origmdf
    string  img, outimg, suf, aplist, apsel, range, sdum, tmpout
    string  tmpmsa, tmpmsb, msjoin, tmpjoin, sec1, sec2
    real    etime, ref1, ref2
    int     ii, jj, nbad, nin, nout, nsky, prev, rstart, nrange, dum, nextnd
    int     nsci, len, refpix, nx1, nx2, apid, inst
    bool    useprefix
    struct  sdate, sline

    # Initialize exit status
    status = 0

    # cache some parameter files
    cache ("imgets", "gmos", "gimverify", "tinfo", "tabpar", "gemdate")

    # Initialize local variables
    l_inimages=inimages ; l_outimages=outimages ; l_prefix=outprefix
    l_logfile=logfile ; l_verbose=verbose ; l_fl_inter=fl_inter
    l_sci_ext=sci_ext ; l_var_ext=var_ext ; l_dq_ext=dq_ext 
    l_expr=expr ; l_combine=combine ; l_reject=reject
    l_scale=scale ; l_zero=zero ; l_weight=weight
    l_lthreshold=lthreshold ; l_hthreshold=hthreshold
    l_nlow=nlow ; l_nhigh=nhigh ; l_nkeep=nkeep
    l_lsigma=lsigma ; l_hsigma=hsigma ; l_mclip=mclip
    l_gain=key_gain ;l_snoise=snoise ; l_rdnoise=key_ron
    l_sigscale=sigscale ; l_pclip=pclip ; l_grow=grow ; l_blank=blank
    l_aper=apertures

    # Start logfile
    gemlogname (logpar=l_logfile, package="gmos")
    if (gemlogname.status != 0)
        goto error
    l_logfile = gemlogname.logname

    printlog ("-------------------------------------------------------------\
        ------------------", l_logfile, l_verbose)
    date | scan(sdate)
    printlog ("GFAPSUM -- "//sdate, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    #The usual kind of checks to make sure everything we need is specified. 
    nbad = 0
    if (l_inimages=="" || l_inimages==" ") {
        printlog ("ERROR - GFAPSUM: Input spectra is an empty string",
            l_logfile, yes)
        nbad=nbad+1 
    }

    if (substr(l_inimages,1,1)=="@") {
        inlist = substr (l_inimages, 2, strlen(l_inimages))
        if (!access(inlist)) {
            printlog ("ERROR - GFAPSUM: Input list "//inlist//" not found",
                l_logfile, verbose+)
            nbad=nbad+1
        }
    }

    if ((l_outimages=="" || l_outimages==" ") \
        && (l_prefix=="" || l_prefix==" ")) {
        
        printlog ("ERROR - GFAPSUM: Neither the output spectra nor prefix is \
            specified.", l_logfile, yes)
        nbad = nbad+1
    } else if ((l_outimages!="" && l_outimages!=" "))
        useprefix = no
    else
        useprefix = yes

    if (substr(l_outimages,1,1)=="@") {
        outlist = substr (l_outimages,2,strlen(l_outimages))
        if (!access(outlist)) {
            printlog ("ERROR - GFAPSUM: Output list "//outlist//" not found",
                l_logfile, yes)
            nbad = nbad+1
        }
    }

    # Temporary files
    temp1 = mktemp("tmpin")
    temp2 = mktemp("tmpout") 
    temp3 = mktemp("tmpfilelist") 
    apsel = mktemp("tmpapsel")
    # Temporary files that will be defined later.
    tmpout=""
    mdf=""
    origmdf=""
    tmpjoin=""
    tmpmsa=""
    tmpmsb=""
    msjoin=""

    files (l_inimages,sort-, > temp1)
    count (temp1) | scan (nin)
    if (!useprefix) { 
        files (l_outimages,sort-, > temp2)
        count(temp2) | scan(nout)
        if (nin != nout) {
            printlog("ERROR - GFAPSUM: Different number of input and output \
                spectra", l_logfile, yes)
            nbad = nbad+1
        }
    } else {
        files (l_prefix//"//@"//temp1, sort-, > temp2)
    }

    # check existence of output files
    scanfile2=temp2
    while (fscan(scanfile2,img) != EOF) {
        if (imaccess(img)) {
        printlog ("ERROR - GFSKYSUB: "//img//" already exists", l_logfile, yes)
        nbad = nbad+1
        }
    }
    scanfile2 = ""

    #i=0
    scanfile1=temp1
    while (fscan(scanfile1,img) !=EOF) { 
        gimverify (img)
        if (gimverify.status>0) {
            nbad = nbad+1
        } else {
            #i=i+1
            # name w/o suffix
            #infile[i]=gimverify.outname
        }

    } #end of while loop over images. 
    scanfile1 = ""

    #If anything was wrong then exit. 
    if (nbad > 0) { 
        printlog("ERROR - GFAPSUM: "//nbad//" errors found with input \
            parameters. Exiting.", l_logfile, yes) 
        goto error
    } 

    #If we are here then everything should be OK. 
    #Write all the relevant info to the logfile:
    #
    printlog ("inimages  = "//l_inimages, l_logfile, l_verbose) 
    printlog ("outimages = "//l_outimages, l_logfile, l_verbose) 
    printlog ("outprefix = "//l_prefix ,l_logfile, l_verbose) 
    printlog ("apertures = "//l_aper, l_logfile, l_verbose) 
    printlog ("expr      = "//l_expr, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    # Make list file for main loop
    joinlines (temp1//","//temp2, output=temp3, delim=" ", \
        missing="Missing", maxchar=161, shortest=yes, verbose=no)
    scanfile1 = temp3 
    while (fscan(scanfile1,img,outimg) != EOF) {
        # Create temporary file names, used within this loop, for FITS files.
        tmpout = mktemp("tmpout")
        mdf = mktemp("tmpmdf")
        origmdf = mktemp("tmporigmdf")
        tmpjoin = mktemp("tmpjoin")
        tmpmsa = mktemp("tmpmsa")
        tmpmsb = mktemp("tmpmsb")
        msjoin = mktemp("tmpmsjoin")

        # check .fits
        suf = substr(img,strlen(img)-4,strlen(img))
        if (suf!=".fits")
            img = img//".fits"

        suf = substr(outimg,strlen(outimg)-4,strlen(outimg))
        if (suf!=".fits")
            outimg=outimg//".fits"

        # Determine which instrument
        if (l_expr=="default") {
            imgets (img//"[0]", "INSTRUME", >>& "dev$null")
            if (imgets.value == "0") {
                printlog ("ERROR - GFSKYSUB: Instrument keyword not found.",
                    l_logfile, verbose+)
                goto error
            }
            inst = 1 # Default is GMOS-N, support for old data
            if (imgets.value == "GMOS-S")
                inst=2

            # Default selection region is the large IFU field
            if (inst==1)
                l_expr="XINST < 10."
            else
                l_expr="XINST > 10."
        }

        # Copy MDF
        tcopy (img//"[MDF]", origmdf//".fits", verbose-)

        # If 2 separate slits, join
        imgets (img//"[0]", "NSCIEXT", >>& "dev$null")
        nsci = int(imgets.value)
        if (nsci == 2) {
            printlog ("Joining slits", l_logfile, l_verbose)
            imgets (img//"["//l_sci_ext//",1]", "i_naxis1", >& "dev$null")
            nx1 = int(imgets.value)
            imgets (img//"["//l_sci_ext//",1]", "refpix1", >& "dev$null")
            ref1 = real(imgets.value)
            imgets (img//"["//l_sci_ext//",2]", "i_naxis1", >& "dev$null")
            nx2 = int(imgets.value)
            imgets (img//"["//l_sci_ext//",2]", "refpix1", >& "dev$null")
            ref2 = real(imgets.value)

            len = min(nx1,nx2)
            refpix = min(ref1,ref2)
            sec1 = "["//(nint(ref1-refpix+1.))//":\
                "//(nint(ref1-refpix+len))//",*]"
            sec2 = "["//(nint(ref2-refpix+1.))//":\
                "//(nint(ref2-refpix+len))//",*]"
            #print(sec1, sec2)

            imcopy (img//"["//l_sci_ext//",1]"//sec1, tmpmsa//".ms", verbose-)
            imcopy (img//"["//l_sci_ext//",2]"//sec2, tmpmsb//".ms", verbose-)
            scopy (tmpmsa//".ms,"//tmpmsb//".ms",msjoin//".ms", renum+, merge+,
                clobber+, verbose-, rebin-)
            wmef (msjoin//".ms.fits", tmpjoin, extnames=l_sci_ext, 
                phu=img//"[0]", verbose-, >>& "dev$null")
            gemhedit (tmpjoin//"[1]", "EXTVER", 1, "", delete-)

            # Update APID in MEF
            tinfo (origmdf//".fits", ttout-)
            apid = 0
            for (ii=1; ii<=tinfo.nrows; ii+=1) {
                tabpar (origmdf//".fits", "BEAM", ii)
                if (int(tabpar.value) != -1) {
                    apid = apid+1
                    partab (apid,origmdf//".fits", "APID", ii)
                }
            }
            img = tmpjoin//".fits"
        } 

        # Select the spectra
        if (l_aper == "" || l_aper == " ") {
            tcalc (origmdf//".fits",
                "BEAM","if BEAM != -1 && "//l_expr//" then 1 else 0")
            tselect (origmdf//".fits",mdf//".fits", "BEAM > 0")
            aplist = ""
            tinfo (mdf//".fits", ttout-)
            nsky = 0
            prev = 0
            range = ""
            nrange = 0
            jj = 0
            for (ii=1; ii<=tinfo.nrows; ii+=1) {
                tabpar (mdf//".fits", "BEAM", ii)
                if (tabpar.value=="1") {
                    tabpar (mdf//".fits", "APID", ii)
                    jj = int(tabpar.value)
                    if (prev==0) {
                        prev = jj
                        rstart = jj
                        nrange = 1
                    } else if (jj-prev != 1) {
                        range = str(rstart)
                        if (nrange>1)
                            range = range//"-"//str(prev)
                        if (aplist=="")
                            aplist = range
                        else
                            aplist = aplist//","//range
                        if (ii==tinfo.nrows)
                            aplist = aplist//","//jj
                        prev = jj
                        rstart = jj
                        range = ""
                        nrange = 1
                    } else if (ii==tinfo.nrows) {
                        range = str(rstart)
                        if (nrange>1)
                            range = range//"-"//str(jj)
                        if (aplist=="")
                            aplist = range
                        else
                            aplist = aplist//","//range
                    } else {
                        nrange = nrange+1
                        prev = jj
                    }
                    nsky = nsky+1
                }
            }
            printlog ("Found "//nsky//" spectra", l_logfile, l_verbose)
            #print (aplist)

            if (l_fl_inter) {
                specplot (img//"["//l_sci_ext//",1]", apertures=aplist,
                    autolayout-, logfile=apsel)

                # Make new aperture list 
                aplist = ""
                nsky = 0
                prev = 0
                range = ""
                nrange = 0
                scanfile2 = apsel
                #read header lines
                for (ii=1; ii<=5; ii+=1) {
                    dum = fscan(scanfile2,sline)
                }
                dum = fscan(scanfile2,sline)
                while (sline != "") {
                    print (sline) | translit ("STDIN","()"," ",delete-) | \
                        scan (dum,sdum,ii)
                    if (prev==0) {
                        prev = ii
                        rstart = ii
                        nrange = 1
                    } else if (ii-prev != 1) {
                        range = str(rstart)
                        if (nrange>1)
                            range = range//"-"//str(prev)
                        if (aplist=="")
                            aplist = range
                        else
                            aplist = aplist//","//range
                        prev = ii
                        rstart = ii
                        range = ""
                        nrange = 1
                    } else {
                        nrange = nrange+1
                        prev = ii
                    }
                    nsky = nsky+1
                    dum = fscan(scanfile2,sline)
                }
                range = str(rstart)
                if (nrange>1)
                    range = range//"-"//str(prev)
                if (aplist=="")
                    aplist = range
                else
                    aplist = aplist//","//range
                scanfile2 = ""
                #print(aplist)
            }
        } else
            aplist = l_aper

        printlog ("Final aperture list: "//aplist, l_logfile, l_verbose)

        # Make the combined spectrum
        scombine (img//"["//l_sci_ext//",1]", tmpout, noutput="", logfile="",
            apertures=aplist, group="images", comb=l_combine, reject=l_reject,
            scale=l_scale, zero=l_zero, weight=l_weight,
            lthreshold=l_lthreshold, hthreshold=l_hthreshold, nlow=l_nlow,
            nhigh=l_nhigh, nkeep=l_nkeep, mclip=l_mclip, lsigma=l_lsigma,
            hsigma=l_hsigma, rdnoise=l_rdnoise, gain=l_gain, snoise=l_snoise,
            pclip=l_pclip, grow=l_grow, blank=l_blank)

        wmef (tmpout//".fits", outimg, extnames=l_sci_ext, phu=img//"[0]",
            verbose-, >& "dev$null")
        gemhedit (outimg//"[1]", "EXTVER", 1, "", delete-)
        if (l_combine=="sum") {
            # get exposure time from input image and update combined spectrum
            imgets (img//"[0]", "EXPTIME")
            etime = real(imgets.value)
            gemhedit (outimg//"[1]", "EXPTIME", etime, "", delete-)
        }
        # Add MDF
        fxinsert (mdf//".fits", outimg//"[1]", groups="1", verbose-,
             >& "dev$null")
        gemhedit (outimg//"[0]", "NEXTEND", 2, "", delete-)

        # update phu
        gemdate ()
        gemhedit (outimg//"[0]", "GFAPSUM", gemdate.outdate, 
            "UT Time stamp for GFAPSUM", delete-)
        gemhedit (outimg//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)

        # clean up
        delete (apsel//","//tmpout//".fits", verify-, >& "dev$null")
        delete (mdf//".fits,"//origmdf//".fits", verify-, >& "dev$null")
        imdelete (tmpjoin, verify-, >& "dev$null")
        imdelete (tmpmsa//".ms,"//tmpmsb//".ms,"//msjoin//".ms", verify-,
            >& "dev$null")
    }
    scanfile1 = ""

    goto clean

error:
    status=1

clean:
    delete (temp1//","//temp2//","//temp3, verify-, >& "dev$null") 
    delete (apsel//","//tmpout//".fits", verify-, >& "dev$null")
    delete (mdf//".fits,"//origmdf//".fits", verify-, >& "dev$null")
    imdelete (tmpmsa//".ms,"//tmpmsb//".ms,"//msjoin//".ms,"//tmpjoin,
        verify-, >>& "dev$null")
    printlog ("", l_logfile, l_verbose)
    scanfile1 = ""
    if (status==1)
        printlog ("GFAPSUM exit status: error", l_logfile, l_verbose)
    else
        printlog ("GFAPSUM exit status: good", l_logfile, l_verbose)
    printlog ("------------------------------------------------------------\
        -------------------", l_logfile, l_verbose)

end
