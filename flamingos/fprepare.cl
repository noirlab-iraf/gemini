# Copyright(c) 2001-2009 Association of Universities for Research in Astronomy, Inc.

procedure fprepare(inimages)

# Procedure to take raw Flamingos data, convert it to MEF 
# and create the preliminary VAR and DQ frames if requested.
# Script based on niri.nprepare
#
# Data is "fixed" to account for low-noise reads, digital averages,
# and coadds as follows:
#    output = input / (n reads)
#    noise = [? e- / ( sqrt(n reads)*sqrt(n digital avs) ) * sqrt(n coadds)
#    exptime = input * (n coadds)
#    gain = ? e-
#    saturation = x * (n coadds) where x depends on the bias voltage 
#
#    non-linear for > y*saturation  where y needs to be determined
#
#    (these are nominal values; the data file must contain entries with
#    the names readnoise, gain, shallowbias, deepbias, fullwell
#    which are used in the appropriate places in the calculations above).
#
# The variance frame is generated as:
#	var = (read noise/gain)^2 + max(data,0.0)/gain
#
# The preliminary DQ frame is constructed by using the bad pixel
# mask to set the 1 bit and saturation level to set the 4 bit.
#
# Version July 19, 2002, ML  , Release v1.4
#         Aug 20, 2003  KL  IRAF2.12 - new parameter, addonly, in hedit
#
#

char    inimages    {prompt="Input Flamingos image(s)"}
char    rawpath     {"",prompt="Path for input raw images"}
char    outimages   {"",prompt="Output image(s)"}
char    outprefix   {"f",prompt="Prefix for output image(s)"}
char    bpm         {"",prompt="Bad pixel mask file"}
char    sci_ext     {"SCI",prompt="Name or number of science extension"}
char    var_ext     {"VAR",prompt="Name or number of variance extension"}
char    dq_ext      {"DQ",prompt="Name or number of data quality extension"}
bool    fl_vardq    {no,prompt="Create variance and data quality frames?"}
char    key_ron     {"RDNOISE",prompt="New header keyword for read noise (e-)"}
char    key_gain    {"GAIN",prompt="New header keyword for gain (e-/ADU)"}
char    key_sat     {"SATURATI",prompt="New header keyword for saturation (ADU)"}
char    key_nonlinear   {"NONLINEA",prompt="New header keyword for non-linear regime (ADU)"}
char    specsec1    {"",prompt="Region for 1st band of spectra (Longslit)"} 
char    specsec2    {"",prompt="Region for 2nd band of spectra (Longslit)"}
char    database    {"flamingos$fprepare.dat",prompt="Database to use"}
bool    verbose     {yes,prompt="Verbose"}
char    observer    {"",prompt="Observer"}
char    ssa         {"",prompt="SSA"}
char    logfile     {"",prompt="Logfile"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

begin

    char    l_inimages, l_outimages
    char    l_rawpath, l_prefix, l_logfile, l_temp, tmpfile, tmpfile2
    char    in[1000], out[1000]
    char    tmpinimage, tmpphu, tmpphu2
    char    l_observer, l_ssa, suf, l_specsec1, l_specsec2 
    int     i, j, nimages, noutimages, maxfiles, nbad, nimproc
    char    l_key_ron, l_key_gain, l_sci_ext, l_var_ext, l_dq_ext
    char    l_key_sat, l_key_nonlinear, l_bpm, l_pupil, meftmp[1000] 
    char    l_varexpression, l_tempexpression, l_database
    char    l_filter,specmode,l_grism,slitsize,specid,l_specid,fpmask 
    char    spattype,spatwcs
    real    l_ron, l_gain, l_ronrefimg, l_gainimg, l_ronrefspec, l_gainspec
    int     l_lnrs, l_coadds, l_ndavgs, l_sat, l_linlimit
    real    l_biasvolt, l_exptime
    real    l_fullwell, l_deepbias, l_shallowbias
    real    l_pixscale, l_linear
    bool    l_verbose, l_fl_vardq, alreadyfixed[1000], bad, mefconv[1000]
    bool    useprefix,skip[1000]
    struct  l_struct

    status = 0
    maxfiles = 1000
    nimages = 0
    nbad = 0
    bad = no

    # cache imgets - used throughout the script
    cache ("imgets","fparse")

    # set the local variables
    l_inimages=inimages ; l_outimages=outimages ; l_rawpath=osfn(rawpath)
    l_verbose=verbose ; l_prefix=outprefix ; l_logfile=logfile
    l_key_ron=key_ron ; l_key_gain=key_gain 
    l_key_sat=key_sat ; l_key_nonlinear=key_nonlinear ; l_fl_vardq=fl_vardq
    l_sci_ext=sci_ext ; l_var_ext=var_ext ; l_dq_ext=dq_ext ; l_bpm=bpm
    l_database=database ; l_observer=observer ; l_ssa=ssa 
    l_specsec1=specsec1 ; l_specsec2=specsec2  

    #For now, it looks like the number of Non-destructive digital averages 
    #is always 1 for FLAMINGOS
    #Leave it in here just in case this changes later on.
    l_ndavgs=1

    # open temp files
    tmpfile = mktemp("tmp1")
    tmpfile2 = mktemp("tmp2")
    tmpphu = mktemp("tmpphu")
    tmpphu2 = mktemp("tmpphu2")
    #----------------------------------------------------------------------


    # Test that name of logfile makes sense
    cache("flamingos", "gemdate")
    print(l_logfile) | scan(l_logfile)
    if (l_logfile=="" || l_logfile==" ") {
        l_logfile=flamingos.logfile
        print(l_logfile) | scan(l_logfile)
        if ((l_logfile=="") || (l_logfile==" ")) {
            l_logfile="flamingos.log"
            printlog("WARNING - FPREPARE: Both fprepare.logfile and \
                flamingos.logfile are empty.",logfile=l_logfile, verbose+)
            printlog("                    Using default; flamingos.log.",
                logfile=l_logfile, verbose+)
        }
    }
    # Open log file
    date | scan(l_struct)
    printlog("--------------------------------------------------------------\
        --------------",l_logfile,l_verbose)
    printlog("FPREPARE -- "//l_struct,l_logfile,l_verbose)
    printlog(" ",l_logfile,l_verbose)

    #----------------------------------------------------------------------

    l_fullwell=0. ; l_deepbias=0. ; l_shallowbias=0.0
    l_ronrefimg=0. ; l_ronrefspec=0. ; l_gainimg=0. ; l_gainspec=0.
    l_linear=0.

    if (!access(l_database)) {
        printlog("ERROR   - FPREPARE: Database file "//l_database//" not \
            found.",l_logfile,verbose+)
        status=1
        goto clean
    } else {
        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
            match("readnoiseimg",stop-) | scan(l_temp,l_ronrefimg)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
            match("readnoisespec",stop-) | scan(l_temp,l_ronrefspec)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
            match("gainimg",stop-) | scan(l_temp,l_gainimg)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
            match("gainspec",stop-) | scan(l_temp,l_gainspec)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
            match("fullwell",stop-) | scan(l_temp,l_fullwell)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
            match("shallowbias",stop-) | scan(l_temp,l_shallowbias)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
            match("deepbias",stop-) | scan(l_temp,l_deepbias)

        fields(l_database,"1-2",lines="1-",quit_if_miss-,print_file_n-) | \
            match("linearlimit",stop-) | scan(l_temp,l_linear)
    }

    if (l_ronrefimg==0. || l_ronrefspec==0. || l_gainimg==0. || l_gainspec==0. || l_fullwell==0. || l_deepbias==0. || l_shallowbias==0. || l_linear==0.) {
        printlog("ERROR   - FPREPARE: Array characteristic entry not found in",
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
            printlog("ERROR   - FPREPARE:  Input file "//l_temp//" not found.",
                l_logfile,verbose+)
            status=1
            goto clean
        }
    }


    # parse wildcard and comma-separated lists
    if (substr(l_inimages,1,1)=="@") 
        scanfile=substr(l_inimages,2,strlen(l_inimages))
    else {
        files(l_inimages,sort-, >tmpfile)
        scanfile=tmpfile
    }  

    while (fscan(scanfile,l_temp) != EOF) {
        files(l_rawpath//l_temp, >> tmpfile2)
    }

    scanfile=tmpfile2 
    delete(tmpfile,verif-, >& "dev$null")

    while (fscan(scanfile,l_temp) != EOF) {
        nimages+=1
        # remove rawpath
        suf=substr(l_temp,strlen(l_temp)-4,strlen(l_temp))
        if (suf!=".fits" && imaccess(l_temp//".fits")) 
            l_temp=l_temp//".fits"
        in[nimages]=l_temp

        # check to see if the file exists (or if it's not MEF FITS)
        #This is the normal condition to start with as FLAMINGOS data is not 
        #inherently MEF, but it is possible that one called FPREPARE a 2nd 
        #time on converted FLAMINGOS data to add VAR/DQ planes 

        gimverify(in[nimages])
        skip[nimages]=no
        if (gimverify.status==1) {
            printlog("ERROR   - FPREPARE: Input image "//l_rawpath//l_temp//\
                " not found.",l_logfile,verbose+)
            nbad+=1
            bad=yes
        } else if (gimverify.status == 0) { 
            mefconv[nimages]=no
            if (imaccess(in[nimages]//"["//l_var_ext//"]")) {
                printlog("ERROR   - FPREPARE: Image: "//in[nimages]//\
                    " already has attached VAR/DQ planes.",l_logfile,l_verbose)
                printlog("                    There is nothing to do for this\
                    image.",l_logfile,l_verbose)
                skip[nimages]=yes
            }
            if ((!l_fl_vardq) && (!imaccess(in[nimages]//"["//l_var_ext//"]"))) { 
                printlog("ERROR   - FPREPARE: Image : "//in[nimages]//\
                    " is already in MEF and",l_logfile,l_verbose)
                printlog("                    you have not requested to create \
                    VAR/DQ planes. Nothing to do.",l_logfile,l_verbose)
                skip[nimages]=yes
            }
        } else if (gimverify.status==4)
            mefconv[nimages]=yes

        # check to see if already fixed with fprepare

        if (!bad || !mefconv[nimages]) {
            imgets(in[nimages]//"[0]","FPREPARE", >& "dev$null")
            if (imgets.value != "0") {
                alreadyfixed[nimages]=yes
                printlog("WARNING - FPREPARE: Image "//l_rawpath//l_temp//\
                    " already fixed using FPREPARE.",l_logfile,verbose+) 
                printlog("                    Data scaling and correction for \
                    number of reads (etc.)",l_logfile,verbose+)
                printlog("                    will not be performed.",
                    l_logfile,verbose+)
            } else 
                alreadyfixed[nimages]=no
        } # end if(!bad)
    } # end while(fscan) loop

    nimproc=0
    for (i=1;i<=nimages;i+=1) { 
        if (!skip[i]) nimproc+=1
    } 

    printlog(" ",l_logfile,l_verbose)
    printlog("Processing "//nimproc//" files",l_logfile,l_verbose)
    printlog(" ",l_logfile,l_verbose)

    if (nimages==1 && skip[1]) { 
        printlog("ERROR   - FPREPARE: Nothing to do for image: "//in[1]//\
            ", Exiting.",l_logfile,l_verbose)
        nbad+=1
        status=1
        goto clean
    }

    # Exit if problems found
    if (nbad > 0) {
        printlog("ERROR   - FPREPARE: "//nbad//" image(s) either do not \
            exist or conversion to MEF failed.",l_logfile,verbose+)
        status=1
        goto clean
    }

    scanfile="" ; delete(tmpfile//","//tmpfile2,verif-, >& "dev$null")


    # Now for the output images
    # outimages could contain legal * if it is of a form like %st%stX%*.imh

    noutimages=0
    nbad=0
    print(l_outimages) | scan(l_outimages)
    if (l_outimages!="" && l_outimages!=" ") {
        useprefix=no
        if (substr(l_outimages,1,1)=="@")
            scanfile=substr(l_outimages,2,strlen(l_outimages))
        else { 
            files(l_outimages,sort-, >> tmpfile)
            scanfile=tmpfile
        }
    } else
        useprefix=yes

    # If prefix is to be used instead
    if (useprefix) { 
        print(l_prefix) | scan(l_prefix)
        if (l_prefix=="" || l_prefix==" ") {
            printlog("ERROR   - FPREPARE: Neither output image name nor \
                output prefix is defined.",l_logfile,verbose+)
            status=1
            goto clean
        } else { 
            i=1
            nbad=0
            while (i<=nimages) {
                fparse(in[i],verbose-)
                print(l_prefix//fparse.root//".fits", >> tmpfile)
                i+=1
            } 
            scanfile=tmpfile
        } 
    }

    noutimages=0 
    while (fscan(scanfile,l_temp) != EOF) {
        noutimages=noutimages+1
        if (noutimages > maxfiles) {
            printlog("ERROR   - FPREPARE: Maximum number of output images \
                exceeded ("//str(maxfiles)//")",l_logfile,verbose+ )
            status=1
            goto clean
        }
        out[noutimages]=l_temp
        if (imaccess(out[noutimages])) {
            printlog("ERROR   - FPREPARE: Output image "//out[noutimages]//\
                " exists", l_logfile,verbose+)
            nbad+=1
        }
    } # end while
    scanfile="" ; delete(tmpfile,ver-, >& "dev$null")


    # Exit if problems found
    if (nbad > 0) {
        printlog("ERROR   - FPREPARE: "//nbad//" image(s) already exist.",
            l_logfile,verbose+)
        status=1
        goto clean
    }

    # if there are too many or too few output images at this stage - 
    # exit with error
    if (nimages!=noutimages && l_outimages!="") {
        nbad+=1 
        printlog("ERROR   - FPREPARE: Number of input and output images \
            unequal.",l_logfile,verbose+)
        status=1
        goto clean
    }

    #-------------------------------------------------------------------------
    # Check for existence of bad pixel mask

    if (l_fl_vardq && !imaccess(l_bpm) && l_bpm!="" && stridx(" ",l_bpm)<=0) {
        printlog("WARNING - FPREPARE: Bad pixel mask "//l_bpm//" not found.",
            l_logfile,verbose+)
        printlog("                    Only saturated pixels will be flagged \
            in the DQ frame.",l_logfile,verbose+)
        l_bpm=""
    } else if (l_fl_vardq && (l_bpm=="" || stridx(" ",l_bpm)>0)) {
        printlog("WARNING - FPREPARE: Bad pixel mask is either an empty \
            string or contains",l_logfile,verbose+)
        printlog("                    spaces.  Only saturated pixels will be \
            flagged in the",l_logfile,verbose+)
        printlog("                    DQ frame.",l_logfile,verbose+)
        l_bpm=""
    } else if (!l_fl_vardq && (l_bpm != "" || l_bpm != " ")) { 
        printlog("WARNING - FPREPARE: Bad pixel mask specified but fl_vardq \
            is turned off.",l_logfile,verbose+)
        printlog("                    Nothing to be done with the badpixel \
            mask.",l_logfile,verbose+) 
        l_bpm=""
    }

    #--------------------------------------------------------------------------
    # The big loop:  create VAR and DQ if fl_vardq=yes
    #                  determine read noise, gain, saturation, exposure, etc.
    #If already MEF, VAR/DQ planes do not exist, and fl_vardq=yes, then create 
    #them in a new output MEF 

    #-----------------------------------------------------------------------
    printlog(" ",l_logfile,l_verbose)
    printlog("  n      input file -->      output file",l_logfile,l_verbose)
    printlog("                 filter      grism   input BPM       RON   \
        gain     sat",l_logfile,l_verbose)
    printlog(" ",l_logfile,l_verbose)

    i=1
    while (i<=nimages) { 
        # Create tmp FITS file names used within this loop
        tmpinimage = mktemp("tmpinimage")

        if (!skip[i]) { 
            if (alreadyfixed[i])
                in[i]=in[i]//"[0]"
            imgets(in[i],"COADDS", >& "dev$null")
            if (imgets.value=="0") {
                printlog("ERROR   - FPREPARE: Could not read number of \
                    coadds from header",l_logfile,verbose+)
                printlog("                  of file "//in[i],l_logfile,verbose+)
                status=1
                goto clean
            } else
                l_coadds=real(imgets.value)

            imgets(in[i],"NREADS", >& "dev$null")
            if (imgets.value=="0") {
                printlog("ERROR   - FPREPARE: Could not read number of \
                    non-destructive reads from header",l_logfile,verbose+)
                printlog("                  of file "//in[i],l_logfile,verbose+)
                status=1
                goto clean
            } else
                l_lnrs=real(imgets.value)

            imgets(in[i],"BIAS", >& "dev$null")
            if (imgets.value=="0") { 
                printlog("ERROR   - FPREPARE: Could not find BIAS voltage \
                    setting from header",l_logfile,verbose+)
                printlog("                  of file "//in[i],l_logfile,verbose+)
                status=1
                goto clean
            } else
                l_biasvolt=real(imgets.value)

            if (abs(l_biasvolt-l_shallowbias) < 0.05) {
                l_gain=l_gainspec
                l_ron=(l_ronrefspec*sqrt(l_coadds)) / (sqrt(l_lnrs)*sqrt(l_ndavgs))
                l_sat=int(l_fullwell*l_coadds/l_gain)
            } else if (abs(l_biasvolt-l_deepbias) < 0.05) {
                l_gain=l_gainimg
                l_ron=(l_ronrefimg*sqrt(l_coadds)) / (sqrt(l_lnrs)*sqrt(l_ndavgs))
                l_sat=int(l_fullwell*l_coadds/l_gain)
            } else {
                printlog("ERROR   - FPREPARE: Cannot determine saturation \
                    level from bias voltage.",l_logfile,verbose+)
                status=1
                goto clean
            }

            l_linlimit=int(l_sat*l_linear)

            # Fix the data if number of non-destructive reads > 1
            if ((l_lnrs != 1) && (!alreadyfixed[i])) {
                printlog("Dividing "//in[i]//" by number of non-destructive \
                    reads: "//l_lnrs,l_logfile,l_verbose) 
                imexpr("a/"//l_lnrs,tmpinimage,in[i],outtype="real",
                    >& "dev$null")
                in[i]=tmpinimage
            }

            suf=substr(out[i],strlen(out[i])-3,strlen(out[i]))
            if (suf!="fits")
                out[i]=out[i]//".fits"
            if (alreadyfixed[i]) 
                in[i]=substr(in[i],1,strlen(in[i])-3)
            if (mefconv[i]) { 
                imcopy (in[i]//"[0]", out[i], verbose=no)
                imcopy(in[i],out[i]//"[1,append+]",verbose-)
                hedit(out[i]//"[0]","EXTEND","T",add-,addonly-,delete-,verify-,
                    show-,update+)
                gemhedit(out[i]//"[0]","NEXTEND",1,"Number of extensions")
                gemhedit(out[i]//"[0]","NSCIEXT",1,
                    "Number of science extensions")
                gemhedit(out[i]//"[0]","OBSERVAT","Gemini-South",
                    "Name of Observatory")
                gemhedit(out[i]//"[0]","NDAVGS",1,"Number of digital averages")
                gemhedit(out[i]//"[0]","PIXSCALE",0.078,
                    "Pixel scale in arcsec/pixel")
                imgets(out[i]//"[0]","DECKER", >& "dev$null")
                specmode=imgets.value
                imgets(in[i],"GRISM",>& "dev$null")
                l_grism=(imgets.value)
                if (specmode == "imaging" && (l_grism == "JH" || l_grism == "HK"))
                    specmode="slitless"  
                if (specmode == "slit" || specmode == "mos" || specmode == "slitless") {
                    imgets(in[i],"FILTER",>& "dev$null")
                    l_filter=(imgets.value)
                    imgets(in[i],"SLIT",>& "dev$null") 
                    slitsize=(imgets.value)
                    specid=l_grism//"+"//slitsize
                    fpmask=specid
                    if (!access(l_database)) {
                        printlog("WARNING - FPREPARE: Database file "//\
                            l_database//" not found.  Using",l_logfile,verbose+)
                        printlog("                    value of "//l_specsec1//\
                            " for SPECSEC1.",l_logfile,verbose+)
                        printlog("WARNING - FPREPARE: Database file "//\
                            l_database//" not found.  Using",l_logfile,verbose+)
                        printlog("                    value of "//l_specsec2//\
                            " for SPECSEC2.",l_logfile,verbose+)
                    } else { 
                        fields(l_database,"1-4",lines="1-",quit_if_miss-,
                            print_file_n-) | match(specid,stop-) | \
                            scan(l_specid,l_specsec1,l_specsec2) 
                        print(l_specsec1) | scan(l_specsec1)
                        print(l_specsec2) | scan(l_specsec2)
                    }             

                    #We should preserve the WCS spatial info when in 
                    #spectroscopic mode since these keywords will be deleted 
                    #later on by nsappwave when we apply a wavelength solution
                    
                    imgets(in[i],"CTYPE1",>& "dev$null")
                    spattype=(imgets.value)
                    gemhedit(out[i]//"[0]","CTYPE1_S",spattype,
                        "Original WCS (Spatial)")
                    imgets(in[i],"CTYPE2",>& "dev$null")
                    spattype=(imgets.value)
                    gemhedit(out[i]//"[0]","CTYPE2_S",spattype,
                        "Original WCS (Spatial)")
                    imgets(in[i],"CRVAL1",>& "dev$null")
                    spatwcs=(imgets.value)
                    gemhedit(out[i]//"[0]","CRVAL1_S",spatwcs,
                        "Original WCS (Spatial)")
                    imgets(in[i],"CRVAL2",>& "dev$null")
                    spatwcs=(imgets.value)
                    gemhedit(out[i]//"[0]","CRVAL2_S",spatwcs,
                        "Original WCS (Spatial)")
                    imgets(in[i],"CRPIX1",>& "dev$null")
                    spatwcs=(imgets.value)
                    gemhedit(out[i]//"[0]","CRPIX1_S",spatwcs,
                        "Original WCS (Spatial)")
                    imgets(in[i],"CRPIX2",>& "dev$null")
                    spatwcs=(imgets.value)
                    gemhedit(out[i]//"[0]","CRPIX2_S",spatwcs,
                        "Original WCS (Spatial)")
                    imgets(in[i],"CDELT1",>& "dev$null")
                    spatwcs=(imgets.value)
                    gemhedit(out[i]//"[0]","CDELT1_S",spatwcs,
                        "Original WCS (Spatial)")
                    imgets(in[i],"CDELT2",>& "dev$null")
                    spatwcs=(imgets.value)
                    gemhedit(out[i]//"[0]","CDELT2_S",spatwcs,
                        "Original WCS (Spatial)")
                    imgets(in[i],"CROTA1",>& "dev$null")
                    spatwcs=(imgets.value)
                    gemhedit(out[i]//"[0]","CROTA1_S",spatwcs,
                        "Original WCS (Spatial)")
                    hedit(out[i]//"[0]",
                        "COM_A,COM_AA,CTYPE1,CTYPE2,CRVAL1,CRVAL2,CRPIX1,CRPIX2,CDELT1,CDELT2,CROTA1",
                        add-,addonly-,delete+,verify-,show-,update+)

                    #Now finish header updates
                    gemhedit(out[i]//"[0]","SPECSEC1",l_specsec1,
                        "Region for 1st band of spectra") 
                    gemhedit(out[i]//"[0]","SPECSEC2",l_specsec2,
                        "Region for 2nd band of spectra") 
                    gemhedit(out[i]//"[1]","SPECSEC1",l_specsec1,
                        "Region for 1st band of spectra") 
                    gemhedit(out[i]//"[1]","SPECSEC2",l_specsec2,
                        "Region for 2nd band of spectra") 
                    gemhedit(out[i]//"[0]","DISPAXIS",1,
                        "Dispersion axis (1=along lines, 2=along columns)")
                    gemhedit(out[i]//"[1]","DISPAXIS",1,
                        "Dispersion axis (1=along lines, 2=along columns)")
                    gemhedit(out[i]//"[0]","SPECMODE",specmode,
                        "Spectroscopy mode (slit|mos|slitless)")
                    gemhedit(out[i]//"[0]","FPMASK","+"//fpmask,
                        "SpecSetup: Focal-Plane Mask ID (NIRI)")
                    gemhedit(out[i]//"[1]","FPMASK","+"//fpmask,
                        "SpecSetup: Focal-Plane Mask ID (NIRI)")
                    gemhedit(out[i]//"[0]","FLAMMASK",
                        l_filter//"+"//fpmask//"-SETUP","SpecSetup: Flamingos")  
                    gemhedit(out[i]//"[1]","FLAMMASK",
                        l_filter//"+"//fpmask//"-SETUP","SpecSetup: Flamingos")  
                    gemhedit(out[i]//"[0]","CAMERA","",
                        "Dummy keyword to emulate NIRI") 
                }
            } else {  #if already MEF but adding VAR/DQ, create new output MEF file 
                fxcopy(in[i],out[i],groups="",new_file+,verbose-, >& "dev$null")
            } #end if(mefconv[i]) 

            imgets(out[i]//"[0]","EXP_TIME", >& "dev$null")
            l_exptime=real(imgets.value)
            gemhedit(out[i]//"[0]","EXPTIME",l_exptime,"Integration time (sec)")
            gemhedit(out[i]//"[1]","EXPTIME",l_exptime,"Integration time (sec)")

            if (l_observer != "" && l_observer != " ")
                gemhedit(out[i]//"[0]","OBSERVER",l_observer,"Observer")
            if (l_ssa != "" && l_ssa != " ")
                gemhedit(out[i]//"[0]","SSA",l_ssa,"SSA")

            if (l_fl_vardq) {
                # create the variance frame, if it doesn't already exist
                # The variance frame is generated as:
                #  var = (read noise/gain)^2 + max(data,0.0)/gain

                l_varexpression=\
                    "((max(a,0.0))/"//l_gain//"+("//l_ron//"/"//l_gain//")**2)"

                imgets(out[i]//"[0]",l_key_gain,>& "dev$null")
                if (imgets.value != "0") {
                    l_gain=real(imgets.value)
                    printlog("WARNING - FPREPARE: Gain already set in image \
                        header.  Using a gain",l_logfile,verbose+)
                    printlog("                    of "//l_gain//" electrons \
                        per ADU for image "//in[i],l_logfile,verbose+)
                    printlog(" ",l_logfile,l_verbose)
                }   

                imgets(out[i]//"[0]",l_key_ron,>& "dev$null")
                if (imgets.value != "0") {
                    l_ron=real(imgets.value)
                    printlog("WARNING - FPREPARE: Read noise already set in \
                        image header.",l_logfile,verbose+)
                    printlog("                    Using a read noise of "//\
                        l_ron//" electrons.",l_logfile,verbose+)
                    printlog(" ",l_logfile,l_verbose)
                }

                imgets(out[i]//"[0]",l_key_sat,>& "dev$null")
                if (imgets.value != "0") {
                    l_sat=real(imgets.value)
                    printlog("WARNING - FPREPARE: Saturation level already \
                        set in image header.",l_logfile,verbose+)
                    printlog("                    Using a saturation level \
                        of "//l_sat//" ADU.",l_logfile,verbose+)
                    printlog(" ",l_logfile,l_verbose)
                }
                if (alreadyfixed[i])
                    imexpr(l_varexpression,out[i]//"[2,append+,dupnam+]",
                        in[i]//"[1]",outtype="real",verbose-)  
                else 
                    imexpr(l_varexpression,out[i]//"[2,append+,dupnam+]",
                        in[i],outtype="real",verbose-)  

                #-------------
                # create the DQ frame, if it doesn't already exist
                # The preliminary DQ frame is constructed by using the bad pixel
                # mask to set bad pixels to 1, pixels in the non-linear regime 
                # to 2, and saturated pixels to 4.

                l_tempexpression=\
                    "(a>"//l_sat//") ? 4 : ((a>"//l_linlimit//") ? 2 : 0)"
                if (alreadyfixed[i])
                    imexpr(l_tempexpression,out[i]//"[3,append+,dupnam+]",
                        in[i]//"[1]",outtype="short",verbose-)
                else
                    imexpr(l_tempexpression,out[i]//"[3,append+,dupnam+]",
                        in[i],outtype="short",verbose-)

                # If there's no BPM, then just worry about the saturated pixels
                if (l_bpm=="")
                    addmasks(out[i]//"[3]",out[i]//"[3,overwrite+]","im1")
                else
                    addmasks(out[i]//"[3],"//l_bpm,out[i]//"[3,overwrite+]",
                        "im1 || im2")

                #-------------
                # Edit the headers and clean up 

                #hedit(out[i]//"[0]","EXTEND",value="T",add+,addonly-,del-,ver-,show-,updat+)
                hedit(out[i]//"[2]","EXTNAME",value=l_var_ext,add+,addonly-,
                    delete-,verify-,show-,update+)
                hedit(out[i]//"[2]","EXTVER",value=1,add+,addonly-,delete-,
                    verify-,show-,update+)
                hedit(out[i]//"[3]","EXTNAME",value=l_dq_ext,add+,addonly-,
                    delete-,verify-,show-,update+)
                hedit(out[i]//"[3]","EXTVER",value=1,add+,addonly-,delete-,
                    verify-,show-,update+)
                hedit(out[i]//"[0]","NEXTEND",3,add-,addonly-,delete-,verify-,
                    show-,update+)
                gemhedit(out[i]//"[2]","EXPTIME",l_exptime,
                    "Integration time (sec)")
                gemhedit(out[i]//"[3]","EXPTIME",l_exptime,
                    "Integration time (sec)")
            } # end if(l_fl_vardq) 

            imdelete (tmpinimage//".fits",verify-, >& "dev$null")

            if ((l_coadds != 1) && (!alreadyfixed[i])) {
                imgets(out[i]//"[0]","EXPTIME", >& "dev$null")
                l_exptime = real(imgets.value)
                gemhedit(out[i]//"[0]","COADDEXP",l_exptime,
                    "Exposure time (s) for each frame")
                l_exptime = l_exptime * l_coadds
                hedit(out[i]//"[0]","EXPTIME",add-,addonly-,delete+,verify-,
                    show-,update+, >& "dev$null")
                gemhedit(out[i]//"[0]","EXPTIME",l_exptime,
                    "Exposure time (s) for sum of all coadds")
                gemhedit(out[i]//"[1]","EXPTIME",l_exptime,
                    "Exposure time (s) for sum of all coadds")
                if (l_fl_vardq) {
                    gemhedit(out[i]//"[2]","EXPTIME",l_exptime,
                        "Exposure time (s) for sum of all coadds")
                    gemhedit(out[i]//"[3]","EXPTIME",l_exptime,
                        "Exposure time (s) for sum of all coadds")
                }
                printlog("Scaling exposure time for "//out[i]//\
                    ". Coadded exposure time is "//l_exptime//" s.",
                    l_logfile,l_verbose)
            }

            #Now, fixup additional edits of header and clean up 

            # update the header
            hedit(out[i]//"[1]","EXTNAME",value=l_sci_ext,add+,addonly-,delete-,
                verify-,show-,update+)
            hedit(out[i]//"[1]","EXTVER",value=1,add+,addonly-,delete-,verify-,
                show-,update+)

            gemdate ()
            gemhedit(out[i]//"[0]",l_key_ron,l_ron,
                "Estimated read noise (electrons)")
            gemhedit(out[i]//"[0]",l_key_gain,l_gain,
                "Gain (electrons/ADU)")
            gemhedit(out[i]//"[0]",l_key_sat,l_sat,
                "Saturation level in ADU")
            gemhedit(out[i]//"[0]",l_key_nonlinear,l_linlimit,
                "Non-linear regime in ADU")
            gemhedit(out[i]//"[0]","FPREPARE",gemdate.outdate,
                "UT Time stamp for FPREPARE")
            gemhedit(out[i]//"[0]","NPREPARE",gemdate.outdate,
                "UT Time stamp for NPREPARE (be NIRI)")
            gemhedit(out[i]//"[0]","GEM-TLM",gemdate.outdate,
                "UT Last modification with GEMINI")

            # l_bpm cannot be an empty string when put in the header
            # also set to none if no DQ plane was created
            if (l_bpm=="" || !l_fl_vardq) { 
                l_bpm="none"
                gemhedit(out[i]//"[0]","BPMFILE",l_bpm,"Input BPM file")
                l_bpm=""
            }
            imgets(out[i]//"[0]","FILTER", >& "dev$null")
            l_filter=(imgets.value)
            imgets(out[i]//"[0]","GRISM", >& "dev$null")
            l_grism=(imgets.value)
            imgets(out[i]//"[0]",l_key_ron, >& "dev$null")
            l_ron=real(imgets.value)
            imgets(out[i]//"[0]",l_key_gain, >& "dev$null")
            l_gain=real(imgets.value)
            imgets(out[i]//"[0]",l_key_sat, >& "dev$null")
            l_sat=int(imgets.value)

        } # end if (!skip[i]) check 

        if (skip[i]){
            printf("%3.0d %15s --> %17s\n",i,in[i],"Skipped!") | scan(l_struct)
            printlog(l_struct,l_logfile,l_verbose)
        } else { 
            printf("%3.0d %15s --> %17s\n",i,in[i],out[i]) | scan(l_struct)
            printlog(l_struct,l_logfile,l_verbose)
        }
        printf("                  %2s         %5s  %15s %5.1f %5.1f %7.0d \n",
            l_filter,l_grism,l_bpm,l_ron,l_gain,l_sat) | scan(l_struct)
        printlog(l_struct,l_logfile,l_verbose)

        i=i+1

    } # end while loop


clean:
    #--------------------------------------------------------------------------
    # Clean up
    if (status==0) {
        printlog(" ",l_logfile,l_verbose)
        printlog("FPREPARE exit status:  good.",l_logfile,l_verbose)
    } else 
        printlog("FPREPARE exit status:  bad.  Exited with "//nbad//" errors.",
            l_logfile,l_verbose)
    printlog("---------------------------------------------------------------\
        -------------",l_logfile,l_verbose)

    scanfile="" 
    delete(tmpinimage,verif-, >& "dev$null")
    delete(tmpfile//","//tmpfile2,ver-, >& "dev$null")
    delete(tmpphu//","//tmpphu2,ver-, >& "dev$null")

end
