# Copyright(c) 2002-2015 Association of Universities for Research in Astronomy, Inc.

procedure gfskysub (inimages)

# Subtract a combined sky spectrum from a GMOS IFU datacube
# 
# Version   Sept 20, 2002 BM v1.4 release
#           Aug 25, 2003  KL IRAF2.12 - new parameter, addonly, in hedit
#           Apr 07, 2004  BM fix bug in ap selection, 
#                            allow use on non-GSTRANSFORMED 2-slit images,
#                            default expr based on INSTRUME, mktemps moved
#           Dec 17, 2014  JT opt. to subtract 2 slits separately; minor fixes

string inimages     {prompt="Input images"}
string outimages    {"",prompt="Output images"}
string outpref      {"s",prompt="Prefix for output images"}
string apertures    {"",prompt="Aperture list"}
string expr         {"default",prompt="Expression for selection of sky spectra"}
string combine      {"average",enum="average|median",prompt="Combine operation"}
string reject       {"avsigclip",enum="none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip",prompt="Rejection algorithm"}
string scale        {"none",prompt="Image scaling"}
string zero         {"none",prompt="Image zeropoint offset"}
string weight       {"none",prompt="Image weights"}
bool   sepslits     {no,prompt="Subtract each slit separately?"}
real   lthreshold   {INDEF,prompt="Lower threshold"}
real   hthreshold   {INDEF,prompt="Upper threshold"}
int    nlow         {1,min=0,prompt="minmax: Number of low pixels to reject"}
int    nhigh        {1,min=0,prompt="minmax: Number of high pixels to reject"}
int    nkeep        {0,min=0,prompt="Minimum to keep or maximum to reject"}  
bool   mclip        {yes,prompt="Use median in sigma clipping algorithms?"}
real   lsigma       {3.,prompt="Lower sigma clipping factor"}
real   hsigma       {3.,prompt="Upper sigma clipping factor"}
string key_ron      {"RDNOISE",prompt="Keyword for readout noise in e-"}
string key_gain     {"GAIN",prompt="Keyword for gain in electrons/ADU"}
string snoise       {"0.0",prompt="Sensitivity noise (fraction), ccdclip and crreject"}
real   sigscale     {0.1,prompt="Tolerance for sigma clipping scaling correction"}
real   pclip        {-0.5,prompt="pclip: Percentile clipping parameter"}
real   grow         {0.0,prompt="Radius (pixels) for neighbor rejection"}
real   blank        {0.0,prompt="Value if there are no pixels"}
string  sci_ext     {"SCI",prompt="Name of science extension"}
string  var_ext     {"VAR",prompt="Name of variance extension"}
string  dq_ext      {"DQ",prompt="Name of data quality extension"}
string logfile      {"",prompt="Logfile name"}
bool   fl_inter     {yes,prompt="Select sky spectra interactively?"}
bool   verbose      {yes,prompt="Verbose?"}
int    status       {0,prompt="Exit status (0=good)"}
struct *scanfile1   {"",prompt="Internal use only"}
struct *scanfile2   {"",prompt="Internal use only"}

begin

    # Local variables for input parameters
    string l_inimages,l_outimages,l_prefix,l_logfile
    string l_sci_ext, l_var_ext, l_dq_ext
    string l_expr, l_aper
    bool   l_verbose,l_fl_inter, l_mclip, l_sepslits
    string l_combine, l_reject, l_scale, l_zero, l_weight
    string l_key_ron, l_key_gain, l_snoise, l_gain, l_rdnoise
    real   l_lthreshold, l_hthreshold, l_lsigma, l_hsigma
    real   l_grow, l_sigscale, l_pclip, l_blank
    int    l_nlow, l_nhigh, l_nkeep

    # Other local variables
    string inlist,outlist,temp1,temp2,temp3,mdf,mdfbak,tmplist
    string img,outimg,suf,aplist,apsel,range,sdum,sky,slitexpr,slitrange[2]
    string tmpmsa,tmpmsb,msjoin,tmpjoin,sec1,sec2,tmpsub
    real   etime,ref1,ref2
    int    row,apid,nbad,nin,nout,nsky,previd,rstart,nrange,dum,nextnd
    int    len,refpix,nx1,nx2,inst,nslit,orow,firstorow[2],lastorow[2]
    int    nrows,firstrow,firstslit,lastslit,firstap,lastap,naps
    bool   useprefix
    struct sdate,sline

    # Initialize exit status
    status=0

    # cache some parameter files
    cache("imgets","gmos","gimverify","tinfo","tabpar","keypar","gemdate")

    # Initialize local variables
    l_inimages=inimages ; l_outimages=outimages ; l_prefix=outpref
    l_logfile=logfile ; l_verbose=verbose ; l_fl_inter=fl_inter
    l_sci_ext=sci_ext ; l_var_ext=var_ext ; l_dq_ext=dq_ext 
    l_expr=expr
    l_combine=combine ; l_reject=reject
    l_scale=scale ; l_zero=zero ; l_weight=weight ; l_sepslits=sepslits
    l_lthreshold=lthreshold ; l_hthreshold=hthreshold
    l_nlow=nlow ; l_nhigh=nhigh ; l_nkeep=nkeep
    l_lsigma=lsigma ; l_hsigma=hsigma ; l_mclip=mclip
    l_gain=key_gain ;l_snoise=snoise ; l_rdnoise=key_ron
    l_sigscale=sigscale ; l_pclip=pclip ; l_grow=grow ; l_blank=blank
    l_aper=apertures

    # Start logfile
    gemlogname(logpar=l_logfile,package="gmos")
    if (gemlogname.status != 0)
        goto error
    l_logfile=gemlogname.logname

    printlog("---------------------------------------------------------------\
        ----------------",l_logfile,l_verbose)
    date | scan(sdate)
    printlog("GFSKYSUB -- "//sdate,l_logfile,l_verbose)
    printlog(" ",l_logfile,l_verbose)

    #The usual kind of checks to make sure everything we need is specified. 
    nbad=0
    if (l_inimages=="" || l_inimages==" ") {
        printlog("ERROR - GFSKYSUB: Input spectra is an empty string",
            l_logfile, verbose+)
        nbad=nbad+1 
    }

    if (substr(l_inimages,1,1)=="@") {
        inlist=substr(l_inimages,2,strlen(l_inimages))
        if (!access(inlist)) {
            printlog("ERROR - GFSKYSUB: Input list "//inlist//" not found",
                l_logfile,verbose+)
            nbad=nbad+1
        }
    }

    if ((l_outimages=="" || l_outimages==" ") && \
        (l_prefix=="" || l_prefix==" ")) {
        printlog("ERROR - GFSKYSUB: Neither the output spectra nor prefix is \
            specified.", l_logfile, verbose+)
        nbad=nbad+1
    } else if ((l_outimages!="" && l_outimages!=" ")) {
        useprefix=no
    } else
        useprefix=yes

    if (substr(l_outimages,1,1)=="@") {
        outlist=substr(l_outimages,2,strlen(l_outimages))
        if (!access(outlist)) {
            printlog("ERROR - GFSKYSUB: Output list "//outlist//" not found",
                l_logfile, verbose+)
            nbad=nbad+1
        }
    }

    # Temporary files
    temp1=mktemp("tmpin")
    temp2=mktemp("tmpout") 
    temp3=mktemp("tmpfilelist") 
    apsel=mktemp("tmpapsel")
    # Temporary files that will be defined later
    mdf=""
    sky=""
    tmpjoin=""
    tmpmsa=""
    tmpmsb=""
    msjoin=""

    files(l_inimages,sort-, > temp1)
    count(temp1) | scan(nin)
    if (!useprefix) { 
        files(l_outimages,sort-, > temp2)
        count(temp2) | scan(nout)
        if (nin != nout) {
            printlog("ERROR - GFSKYSUB: Different number of input and output \
                spectra", l_logfile, verbose+)
            nbad=nbad+1
        }
    } else {
        files(l_prefix//"//@"//temp1,sort-, > temp2)
    }

    # check existence of output files
    scanfile2=temp2
    while (fscan(scanfile2,img) != EOF) {
        if (imaccess(img)) {
            printlog("ERROR - GFSKYSUB: "//img//" already exists", l_logfile, \
                verbose+)
            nbad=nbad+1
        }
    }
    scanfile2=""

    scanfile1=temp1
    while (fscan(scanfile1,img) !=EOF) { 
        gimverify(img)
        if (gimverify.status>0)
            nbad=nbad+1
    
    } #end of while loop over images. 
    scanfile1=""

    #If anything was wrong then exit. 
    if (nbad > 0) { 
        printlog("ERROR - GFSKYSUB: "//nbad//" errors found with input \
            parameters. Exiting.", l_logfile, verbose+) 
        goto error
    } 

    #If we are here then everything should be OK. 
    #Write all the relevant info to the logfile:
    #
    printlog("inimages = "//l_inimages,l_logfile,l_verbose) 
    printlog("outimages = "//l_outimages,l_logfile,l_verbose) 
    printlog("outpref = "//l_prefix,l_logfile,l_verbose) 
    printlog("apertures = "//l_aper,l_logfile,l_verbose) 
    printlog("expr = "//l_expr,l_logfile,l_verbose) 
    printlog(" ",l_logfile,l_verbose)

    # Make list file for main loop
    joinlines (temp1//","//temp2, output=temp3, delim=" ", \
        missing="Missing", maxchar=161, shortest=yes, verbose=no)

    # Process each image
    scanfile1=temp3 
    while (fscan(scanfile1,img,outimg) !=EOF) {
        # temp files used within the loop
        mdf=mktemp("tmpmdf")
        mdfbak=mktemp("tmpmdfbak")
        sky=mktemp("tmpsky")
        tmpjoin=mktemp("tmpjoin")
        tmpmsa=mktemp("tmpmsa")
        tmpmsb=mktemp("tmpmsb")
        msjoin=mktemp("tmpmsjoin")
        tmpsub=mktemp("tmpsub")
        
        # check .fits
        suf=substr(img,strlen(img)-4,strlen(img))
        if (suf!=".fits")
            img=img//".fits"

        suf=substr(outimg,strlen(outimg)-4,strlen(outimg))
        if (suf!=".fits")
            outimg=outimg//".fits"

        # Determine which instrument
        if (l_expr=="default") {
            imgets(img//"[0]","INSTRUME", >>& "dev$null")
            if (imgets.value == "0") {
                printlog ("ERROR - GFSKYSUB: Instrument keyword not found.",
                    l_logfile, verbose+)
                goto error
            }
            inst=1 # Default is GMOS-N, support for old data
            if (imgets.value == "GMOS-S")
                inst=2

            # Default sky selection region
            if (inst==1)
                l_expr="XINST > 10."
            else
                l_expr="XINST < 10."
        }

        # Check that input has been wavelength calibrated (transformed):
        imgets(img//"[0]","GFTRANSF", >>& "dev$null")
        if (imgets.value == "0") {
            printlog ("ERROR - GFSKYSUB: "//img//" not rectified with \
                GFTRANSFORM",
                l_logfile, verbose+)
            goto error
        }

        # Check the number of input rows to avoid trying to copy too many
        # spectra in cases where the extraction has gone wrong (gfextract
        # isn't accounting properly for omissions in the MDF APID column).
        imgets(img//"["//l_sci_ext//",1]","i_naxis2", >>& "dev$null")
        naps = int(imgets.value)
        if (imgets.value == "0") {
            printlog ("ERROR - GFSKYSUB: "//img//" malformed ", \
                l_logfile, verbose+)
            goto error
        }

        # Copy image to output
        # (Note that the VAR/DQ are currently just propagated unchanged from
        # the input. The VAR should really be recalculated here but the error
        # from not doing so will normally be very small, when many fibres are
        # averaged to get the sky - JT.)
        copy(img,outimg,verbose-)

        # How many rows in the input MDF (number gets reduced later) and
        # starting from which fibre?
        tinfo(outimg//"[MDF]",ttout-)
        nrows = tinfo.nrows
        tabpar(outimg//"[MDF]","NO",1)
        firstrow = int(tabpar.value)

        # What are the original MDF rows corresponding to each slit (or
        # both slits if sepslits=no)?
        if (l_sepslits) {

            # Figure out which slit(s) we're looping over:
            if (firstrow < 751) firstslit=1
            else firstslit=2
            if (nrows < 751) lastslit=firstslit
            else lastslit=2

            for (nslit=firstslit; nslit <= lastslit; nslit+=1) {
	        firstorow[nslit] = 1+(nslit-1)*750
	        lastorow[nslit] = nslit*750
            }

        } else {
            firstslit=1  # treat everything as 1 slit now it's transformed
            lastslit=1
            if (firstrow < 751)
                firstorow[1] = 1
            else
                firstorow[1] = 751
            if (nrows < 751)
                lastorow[1] = firstorow[1] + 749
            else
                lastorow[1] = 1500
        }

        # Loop over slits and use the MDF to determine the range of apertures
        # (image rows) corresponding to each slit:
        for (nslit=firstslit; nslit <= lastslit; nslit+=1) {

            # Extract MDF rows for good fibres in this slit:
            slitexpr = "NO >= "//firstorow[nslit]//" && NO <= "//\
                lastorow[nslit]
            tselect(outimg//"[MDF]", mdf//".fits", "BEAM != -1 && "//slitexpr)

            # Look up corresponding range of image rows (== aperture numbers
            # after gftransform renumbering) from APID column:
            tinfo(mdf//".fits",ttout-)
            nrows = tinfo.nrows
            tabpar(mdf//".fits","APID",1)
            firstap=int(tabpar.value)
            tabpar(mdf//".fits","APID",nrows)
            lastap=int(tabpar.value)

            # Deal with any missing apertures the best we can (result will
            # be imperfect but good for quick-look & without crashing).
            if (lastap > naps) {
                printlog ("WARNING - GFSKYSUB: "//img//": defined \
                    apertures are missing;", l_logfile, verbose+)
                printlog ("                    obj/sky/slit mapping will \
                    be approximate", l_logfile, verbose+)
                printlog ("                    (check extraction & MDF)", \
                    l_logfile, verbose+)
                lastap = naps
	    }

            # Record the renumbered aperture range for this slit (including
            # object spectra) so we can separate and sky subtract them later
            # using scopy:
            slitrange[nslit] = "[*,"//firstap//":"//lastap//"]"

            # Clean up temp MDF:
            imdelete(mdf//".fits", verify-, >>& "dev$null")
        }

        # Create MDF with only good sky fibres:
        # This modified obj/sky BEAM number gets propagated to the output
        # as in the previous gfskysub:
        tcalc (outimg//"[MDF]","BEAM",
            "if BEAM != -1 && "//l_expr//" then 0 else BEAM")
        tselect(outimg//"[MDF]",mdf//".fits","BEAM == 0")
        tinfo(mdf//".fits",ttout-)
        nrows = tinfo.nrows

        # Create a temporary copy of the modified output MDF, since it gets
        # corrupted somewhere below, probably by a FITS kernel bug:
        tcopy(outimg//"[MDF]", mdfbak//".fits", verbose-)

        # Loop over the slits:
        row=1
        for (nslit=firstslit; nslit <= lastslit; nslit+=1) {

            # Automatically select sky fibres?
            if (l_aper == "" || l_aper == " ") {

                nsky=0
                previd=0
                range=""
                nrange=0
                apid=0
                aplist=""

                # print "slit "//nslit//"; starting row "//row

                # Loop over the sky MDF rows corresponding to this slit
                # (NB. not the original row numbers):
                for (; row<=nrows; row+=1) {

                    # Get the aperture number & original MDF row:
                    #
                    # The "aperture" parameters in onedspec tasks seem to match
                    # against the first value of each APNUM header keyword,
                    # which at least for recent data correspond to the APID in
                    # the MDF (after renumbering by gftransform).
                    tabpar(mdf//".fits","APID",row)
                    apid=int(tabpar.value)
                    tabpar(mdf//".fits","NO",row)
                    orow=int(tabpar.value)

                    if (orow <= lastorow[nslit]) {  # apid within this slit

                        if (previd==0) {
                            # Start first scombine fibre range:
                            previd=apid
                            rstart=apid
                            nrange=1
                        } else if (apid-previd != 1) {
                            # Start a new range when non-consecutive:
                            range=str(rstart)
                            if (nrange>1)
                                range=range//"-"//str(previd)
                            if (aplist=="")
                                aplist=range
                            else
                                aplist=aplist//","//range
                            if (row==nrows)
                                aplist=aplist//","//apid
                            previd=apid
                            rstart=apid
                            range=""
                            nrange=1
                        } else if (row==nrows) {
                            # Complete slit's last range when the final
                            # aperture of the final slit is reached:
                            range=str(rstart)
                            if (nrange>1)
                                range=range//"-"//str(apid)
                            if (aplist=="")
                                aplist=range
                            else
                                aplist=aplist//","//range
                        } else {
                            # Increment current range while consecutive:
                            nrange=nrange+1
                            previd=apid
                        }
                        nsky=nsky+1

                    } else {  # past end of slit

                        # Complete slit's last aperture range when a new
                        # slit is reached:
                        range=str(rstart)
                        if (nrange>1)
                            range=range//"-"//str(previd)
                        if (aplist=="")
                            aplist=range
                        else
                            aplist=aplist//","//range

                        # Continue with next slit:
                        break  # (from inner loop over rows)

                    } # end if (apid within this slit)

                } # end loop over MDF rows

                printlog("Found "//nsky//" spectra",l_logfile,l_verbose)
                # print (aplist)  # debug

                # Modify the above list of apertures interactively (slow):
                if (l_fl_inter) {
                    specplot(img//"["//l_sci_ext//",1]",apertures=aplist,
                        autolayout-,logfile=apsel)

                    # Make new aperture list 
                    aplist=""
                    nsky=0
                    previd=0
                    range=""
                    nrange=0
                    scanfile2=apsel
                    #read header lines
                    for (apid=1; apid<=5; apid+=1) {  # not real apids
                        dum=fscan(scanfile2,sline)
                    }
                    dum=fscan(scanfile2,sline)
                    while (sline != "") {
                        print(sline) | translit("STDIN","()"," ",delete-) | \
                            scan(dum,sdum,apid)
                        if (previd==0) {
                            previd=apid
                            rstart=apid
                            nrange=1
                        } else if (apid-previd != 1) {
                            range=str(rstart)
                            if (nrange>1)
                                range=range//"-"//str(previd)
                            if (aplist=="")
                                aplist=range
                            else
                                aplist=aplist//","//range
                            previd=apid
                            rstart=apid
                            range=""
                            nrange=1
                        } else {
                            nrange=nrange+1
                            previd=apid
                        }
                        nsky=nsky+1
                        dum=fscan(scanfile2,sline)
                    }  # end while loop over lines from specplot log
                    range=str(rstart)
                    if (nrange>1)
                        range=range//"-"//str(previd)
                    if (aplist=="")
                        aplist=range
                    else
                        aplist=aplist//","//range
                    scanfile2=""
                    #print(aplist)

                    delete(apsel, verify-, >>& "dev$null")

                } # end (fl_inter+)

            } else {   # user specified sky apertures manually
                aplist=l_aper
            }

            printlog("Final aperture list: "//aplist,l_logfile,l_verbose)

            # Make the combined sky
            scombine(outimg//"["//l_sci_ext//",1]", sky//"_"//nslit,
                noutput="", logfile="", apertures=aplist, group="images",
                comb=l_combine, reject=l_reject, scale=l_scale, zero=l_zero,
                weight=l_weight, lthreshold=l_lthreshold,
                hthreshold=l_hthreshold, nlow=l_nlow, nhigh=l_nhigh,
                nkeep=l_nkeep, mclip=l_mclip, lsigma=l_lsigma,
                hsigma=l_hsigma, rdnoise=l_rdnoise, gain=l_gain,
                snoise=l_snoise, pclip=l_pclip, grow=l_grow, blank=l_blank)

            # Copy the image section for the current slit:
            imcopy(img//"[sci,1]"//slitrange[nslit], tmpsub//"_"//nslit)

            # Subtract the sky
            imarith(tmpsub//"_"//nslit,"-",sky//"_"//nslit//".fits",
                tmpsub//"_"//nslit)

        } # end loop over slits

 	# Join the sky-subtracted images back together:
        tmplist=tmpsub//"_"//firstslit
        for (nslit=firstslit+1; nslit <= lastslit; nslit+=1) {
            tmplist=tmplist//","//tmpsub//"_"//nslit
        }
        imjoin(tmplist, tmpsub, join_dim=2)
        imcopy(tmpsub, outimg//"[SCI,1,overwrite+]",verbose-)
        imdelete(tmpsub//","//tmplist, verify-, >>& "dev$null")

        gemdate ()
        gemhedit (outimg//"[0]", "GFSKYSUB", gemdate.outdate,
            "UT Time stamp for GFSKYSUB", delete-)
        gemhedit (outimg//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)

        # Put sky image in output
        nextnd=1
        for (nslit=firstslit; nslit <= lastslit; nslit+=1) {
            imcopy(sky//"_"//nslit//".fits",outimg//"[SKY,"//nextnd//",append]",
                verbose-)
            gemhedit(outimg//"[SKY,"//nextnd//"]","OBJECT", "MEAN SKY", "",
                delete-)
            nextnd += 1
            imdelete(sky//"_"//nslit//".fits",verify-, >>& "dev$null")
        }
        imgets(outimg//"[0]","NEXTEND")
        nextnd=int(imgets.value)+(lastslit-firstslit+1)
        gemhedit(outimg//"[0]","NEXTEND",nextnd,"",delete-)

        # Add the uncorrupted MDF back to the output image (see above):
        tcopy(mdfbak//".fits", outimg//"[MDF,overwrite]", verbose-)

        # Clean up
        delete(mdf//".fits",verify-, >>& "dev$null")
        delete(mdfbak//".fits",verify-, >>& "dev$null")
        imdelete(tmpjoin,verify-,>>& "dev$null")
        imdelete(tmpmsa//".ms,"//tmpmsb//".ms"//msjoin//".ms",verify-,
            >>& "dev$null")
    }
    scanfile1=""

    goto clean

error:
    status=1

clean:
    delete(temp1//","//temp2//","//temp3//","//mdf//".fits",verify-,
        >& "dev$null") 
    imdelete(tmpmsa//".ms,"//tmpmsb//".ms,"//msjoin//".ms,"//tmpjoin,verify-,
        >>& "dev$null")
    delete(apsel//","//sky//".fits",verify-, >>& "dev$null")
    printlog("",l_logfile,l_verbose)
    scanfile1="" ; scanfile2=""
    if (status==0)
        printlog("GFSKYSUB exit status: good",l_logfile,l_verbose)
    else
        printlog("GFSKYSUB exit status: error",l_logfile,l_verbose)
    printlog("----------------------------------------------------------------\
        ---------------",l_logfile,l_verbose)

end
