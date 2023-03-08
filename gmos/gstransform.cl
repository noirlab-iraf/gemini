# Copyright(c) 2002-2015 Association of Universities for Research in Astronomy, Inc.

procedure gstransform(inimages)

# Rectify GMOS spectra using available S-distortion correction. 
#
# Version   Feb 28, 2002  ML,BM  v1.3 release
#           Apr 17, 2002  BM     gemhedit bug fix for s-distortion
#           Aug 27, 2002  IJ avoid appending to existing output image, removed infile
#                            local parameter which is not used
#                            missing s-dist/wavelength transformations is an error
#                            simplify snum naming
#           Sept 20, 2002    v1.4 release
#           Aug 26, 2003  KL IRAF2.12 - new parameter, addonly, in hedit
#                              Replaced 'hedit' calls with 'gemhedit' calls:
#                              from the arguments it is clear that they should
#                              be 'gemhedit' calls.  [calls deal with variance
#                              and DQ headers ('if (imacess(tmpvar))' block)]
#           May 16, 2004  BM  Separate VAR and DQ plane processing, fix vardq bugs
#           May 19, 2004  BM  add parameters and call to GSAPPWAVE
string  inimages    {prompt="Input GMOS spectra"}
string  outimages   {"",prompt="Output spectra"}
string  outprefix   {"t",prompt="Prefix for output spectra"}
bool    fl_stran    {no,prompt="Apply S-distortion correction"}
string  sdistname   {"",prompt="Names of S-distortions calibrations"}
bool    fl_wavtran  {yes,prompt="Apply wavelength calibration from arc spectrum"}
string  wavtraname  {"",prompt="Names of wavelength calibrations"} 
string  database    {"database",prompt="Directory for calibration files"}
bool    fl_vardq    {no,prompt="Transform variance and data quality planes"} 
string  interptype  {"linear",enum="nearest|linear|poly3|poly5|spline3",prompt="Interpolation type for transform"} 
real    lambda1     {INDEF, prompt="First output wavelength for transform (Ang)"}
real    lambda2     {INDEF, prompt="Last output wavelength for transform (Ang)"}
real    dx          {INDEF, prompt="Output wavelength to pixel conversion ratio for transform (Ang/pix)"}
real    nx          {INDEF, prompt="Number of output pixels for transform (pix)"}
bool    lambdalog   {no,prompt="Logarithmic wavelength coordinate for transform"}
bool    ylog        {no,prompt="Logarithmic y coordinate for transform"}
bool    fl_flux     {yes,prompt="Conserve flux per pixel in the transform"} 
string  gratingdb   {"gmos$data/GMOSgratings.dat",prompt="Gratings database file"}
string  filterdb    {"gmos$data/GMOSfilters.dat",prompt="Filters database file"}
string  key_dispaxis    {"DISPAXIS",prompt="Keyword for dispersion axis"}
int     dispaxis    {1,min=1,max=2,prompt="Dispersion axis"}
string  sci_ext     {"SCI",prompt="Name of science extension"}
string  var_ext     {"VAR",prompt="Name of variance extension"}
string  dq_ext      {"DQ",prompt="Name of data quality extension"}
string  logfile     {"",prompt="Logfile"}
#bool    fl_inter    {yes,prompt="Interactive?"}
bool    verbose     {yes,prompt="Verbose"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use"}

begin

    #Local variable declarations 
    string  l_inimages, l_outimages, l_prefix, l_sdistname
    string  l_wavtran, l_wavtraname, woutlist
    string  l_grating, l_filter, l_gratingdb, l_filterdb, l_key_dispaxis
    string  l_database, l_interptype, l_sci_ext, l_var_ext, l_dq_ext
    string  l_logfile, l_fitcfunc, size_print
    bool    l_fl_vardq, l_fl_wavtran, l_fl_stran
    real    l_lambda1, l_lambda2, l_dx, l_nx
    bool    l_fl_inter, l_lambdalog, l_ylog, l_flux, l_verbose 
    int     l_spatbin, l_specbin, l_fitcxord, l_fitcyord, l_status
    int     l_dispaxis


    #Other variables used within this task 

    file    temp1, temp2, temp3, wavlist, sdistlist, tmpsci, tmpvar, tmpdq, mdf
    string  inlist,outlist,img,outimg,snum
    string  wcalfile,scalfile,wcal,scal,fitname
    string  obsmode,dbprefix, suf, l_stdout
    string  ffile,fname,gname, coordtmp
    int     nbad,nin,nout,nwavfil,ndistfil,nsciext,nextens,i,j  
    int     nsciin,nsciout,sel,nextnd,mdfrow, pos_wavecal[4]
    int     wavheight, wavlength, imgheight, imglength
    bool    useprefix,go, outmade
    struct  sdate

    # Set the local variables
    l_inimages=inimages ; l_outimages=outimages ; l_prefix=outprefix
    l_fl_stran=fl_stran; l_sdistname=sdistname 
    l_fl_wavtran=fl_wavtran ; l_wavtraname=wavtraname ; l_database=database
    #l_fl_inter=fl_inter
    l_lambda1 = lambda1 ; l_lambda2 = lambda2 ; l_dx = dx ; l_nx = nx
    l_fl_vardq=fl_vardq ; l_interptype=interptype ; l_lambdalog=lambdalog
    l_ylog=ylog ; l_flux=fl_flux
    l_sci_ext=sci_ext ; l_var_ext=var_ext ; l_dq_ext=dq_ext 
    l_logfile=logfile ; l_verbose=verbose
    l_gratingdb=gratingdb ; l_filterdb=filterdb
    l_key_dispaxis=key_dispaxis ; l_dispaxis=dispaxis

    # Initialize exit status
    status = 0
    outmade = no

    l_stdout = ""
    if (l_verbose)
        l_stdout = "STDOUT,"

    # cache some parameter files
    cache ("imgets", "gmos", "gemhedit", "gimverify", "gemdate")

    # Check logfile
    if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GSTRANSFORM: Both gstransform.logfile and \
                gmos.logfile fields are empty", l_logfile, verbose+)
            printlog ("                    Using default file gmos.log",
                l_logfile, verbose+)
        }
    }

    # temporary files
    wavlist = mktemp("tmpwavlist") 
    sdistlist = mktemp("tmpsdistlist")
    temp3 = mktemp("tmpfilelist")

    printlog ("-------------------------------------------------------------\
        ------------------", l_logfile, l_verbose)
    date | scan(sdate)
    printlog ("GSTRANSFORM -- "//sdate, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)


    #The usual kind of checks to make sure everything we need is specified. 
    nbad = 0

    if (l_inimages=="" || l_inimages==" ") {
        printlog ("ERROR - GSTRANSFORM: Input spectra is an empty string",
            l_logfile, verbose+)
        nbad = nbad+1 
    }

    if (strstr("@",l_inimages) != 0) {
        inlist = substr (l_inimages, stridx("@",l_inimages)+1, \
                     strlen(l_inimages))
        if (!access(inlist)) {
            printlog ("ERROR - GSTRANSFORM: Input list "//inlist//" not found",
                l_logfile, verbose+)
            nbad = nbad+1
        }
    }

    if ((l_outimages=="" || l_outimages==" ") && (l_prefix=="" || \
        l_prefix==" ")) {
        printlog ("ERROR - GSTRANSFORM: Neither the output spectra nor \
            prefix is specified.", l_logfile, verbose+)
        nbad = nbad+1
    } else if ((l_outimages!="" && l_outimages!=" "))
        useprefix = no
    else
        useprefix = yes

    if (substr(l_outimages,1,1)=="@") {
        outlist = substr (l_outimages, 2, strlen(l_outimages))
        if (!access(outlist)) {
            printlog ("ERROR - GSTRANSFORM: Output list "//outlist//\
                " not found", l_logfile, verbose+)
            nbad = nbad+1
        }
    }

    temp1 = mktemp("tmpin")
    temp2 = mktemp("tmpout") 
    if (strstr("@",l_inimages) != 0) {
        sections (l_inimages, > temp1)
    } else {
        files (l_inimages, sort-, > temp1)
    }
    count(temp1) | scan(nin)
    if (!useprefix) { 
        files (l_outimages, sort-, > temp2)
        count(temp2) | scan(nout)
        if (nin != nout) {
            printlog ("ERROR - GSTRANSFORM: Different number of input and \
                output spectra", l_logfile, verbose+)
            nbad = nbad+1
        }
    } else
        files (l_prefix//"//@"//temp1, sort-, > temp2)

    # Are we going to do anything?
    if (l_fl_wavtran==no && l_fl_stran==no) { 
        printlog ("ERROR - GSTRANSFORM: Neither S-distortion or wavelength \
            correction have been requested", l_logfile, verbose+)
        nbad = nbad+1
    } 


    # S-distortion calibrations
    if (l_fl_stran) { 
        if (substr(l_sdistname,1,1) == "@") {
            soutlist = substr (l_sdistname, 2, strlen(l_sdistname))
            if (!access(soutlist)) {
                printlog ("ERROR - GSTRANSFORM: Output list "//soutlist//\
                    " not found", l_logfile, verbose+)
                nbad = nbad+1
            }
        }
        if (nbad == 0) {
            files (l_sdistname, sort-, > sdistlist)
            count(sdistlist) | scan(ndistfil)
            if (ndistfil != 1 || ndistfil != nin) { 
                printlog ("ERROR - GSTRANSFORM: Number of S-distortion images \
                    must be 1, or ", l_logfile, verbose+)
                printlog ("                     equal to the number of input \
                    images.", l_logfile, verbose+)
                nbad = nbad+1
            } 
            if (ndistfil == 1 && nbad==0) {
                delete (sdistlist, verify-, >& "dev$null")
                for (i=1; i<=nin; i+=1) { 
                    print (l_sdistname, >> sdistlist) 
                } 
            }
        }
    } else {
        for (i=1; i<=nin; i+=1) { 
            print("none", >> sdistlist) 
        } 
    }

    #Wavelength calibrations
    if (l_fl_wavtran) { 
        if (substr(l_wavtraname,1,1) == "@") {
            woutlist = substr (l_wavtraname, 2, strlen(l_wavtraname))
            if (!access(woutlist)) {
                printlog ("ERROR - GSTRANSFORM: Output list "//woutlist//\
                    " not found", l_logfile, verbose+)
                nbad = nbad+1
            }
        }
        if (nbad == 0) { 
            files (l_wavtraname, sort-, > wavlist)
            count(wavlist) | scan(nwavfil)
            if (nwavfil != 1 && nwavfil != nin) { 
                printlog ("ERROR - GSTRANSFORM: Number of wavelength \
                    calibration images must be 1, or ", l_logfile, verbose+)
                printlog ("                     equal to the number of input \
                    images.", l_logfile, verbose+)
                nbad = nbad+1
            } 
            if (nwavfil == 1 && nbad==0) {
                delete (wavlist, verify-, >& "dev$null")
                for (i=1; i<=nin; i+=1) { 
                    files (l_wavtraname, sort-, >> wavlist) 
                } 
            }
        }
    } else {
        for (i=1; i<=nin; i+=1) { 
            print ("none", >> wavlist) 
        } 
    }

    if (l_sci_ext=="" || l_sci_ext==" ") {
        printlog ("ERROR - GSTRANSFORM: extension name sci_ext is missing",
            l_logfile, verbose+)
        nbad = nbad+1
    }

    #If var propogation is requested, make sure the names are given
    if (l_fl_vardq) {
        if (l_dq_ext=="" || l_dq_ext ==" ") {
            printlog ("ERROR - GSTRANSFORM: extension name dq_ext is missing",
                l_logfile, verbose+)
            nbad = nbad+1
        } else if (l_var_ext=="" || l_var_ext ==" ") {
            printlog ("ERROR - GSTRANSFORM: extension name var_ext is missing",
                l_logfile, verbose+)
            nbad = nbad+1
        }
    }

    # check no commas in sci_ext, var_ext and dq_ext
    if (stridx(",",l_sci_ext)>0 || stridx(",",l_var_ext)>0 || \
        stridx(",",l_dq_ext)>0 ) {
        printlog ("ERROR - GSTRANSFORM: sci_ext, var_ext or dq_ext contains \
            commas, give root name only", l_logfile, verbose+)
        nbad = nbad+1
    }

    scanfile = temp1
    while (fscan(scanfile,img) !=EOF) { 
        gimverify (img)
        if (gimverify.status > 0) {
            printlog ("ERROR - GSTRANSFORM: Input image "//img//" does not \
                exist or is not MEF", l_logfile, verbose+)
            nbad = nbad+1
        } 
    }   #end of while loop over input images. 
    scanfile = ""

    scanfile = temp2
    while (fscan(scanfile,img) !=EOF) { 
        gimverify (img)
        if (gimverify.status!=1) {
            printlog ("ERROR - GSTRANSFORM: Output image "//img//" exists",
                l_logfile, verbose+)
            nbad = nbad+1
        } 
    }   #end of while loop over output images. 
    scanfile = ""

    #If anything was wrong then exit. 
    if (nbad > 0) { 
        printlog ("ERROR - GSTRANSFORM: "//nbad//" errors found with input \
            parameters. Exiting.", l_logfile, verbose+) 
        goto error
    } 


    #If we are here then everything should be OK. 
    #Write all the relevant info to the logfile:
    #
    printlog ("",l_logfile,l_verbose)
    printlog ("inimages   = "//l_inimages, l_logfile, l_verbose) 
    printlog ("outimages  = "//l_outimages, l_logfile, l_verbose) 
    printlog ("outprefix  = "//l_prefix, l_logfile, l_verbose) 
    printlog ("fl_stran   = "//l_fl_stran, l_logfile, l_verbose)
    if (l_fl_stran)
        printlog ("sdistname  = "//l_sdistname, l_logfile, l_verbose)
    printlog ("fl_wavtran = "//l_fl_wavtran ,l_logfile, l_verbose)
    if (l_fl_wavtran)
        printlog ("wavtraname = "//l_wavtraname, l_logfile, l_verbose)
    printlog ("database   = "//l_database, l_logfile, l_verbose)
    #printlog ("fl_inter = "//l_fl_inter, l_logfile, l_verbose)
    printlog ("fl_vardq   = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("interptype = "//l_interptype, l_logfile, l_verbose)
    if (l_lambda1!=INDEF)
        printlog ("lambda1    = "//str(l_lambda1), l_logfile, 
            verbose=l_verbose)
    else
        printlog ("lambda1    = INDEF", l_logfile, verbose=l_verbose)
    if (l_lambda2!=INDEF)
        printlog ("lambda2    = "//str(l_lambda2), l_logfile, \
            verbose=l_verbose)
    else
        printlog ("lambda2    = INDEF", l_logfile, verbose=l_verbose)
    if (l_dx!=INDEF)
        printlog ("dx         = "//str(l_dx), l_logfile, verbose=l_verbose)
    else
        printlog ("dx         = INDEF", l_logfile, verbose=l_verbose)
    if (l_nx!=INDEF)
        printlog ("nx         = "//str(l_nx), l_logfile, verbose=l_verbose)
    else
        printlog ("nx         = INDEF", l_logfile, verbose=l_verbose)
    printlog ("lambdalog  = "//l_lambdalog, l_logfile, l_verbose)
    printlog ("ylog       = "//l_ylog, l_logfile, l_verbose)
    printlog ("fl_flux    = "//l_flux, l_logfile, l_verbose)
    printlog ("gratingdb  = "//l_gratingdb, l_logfile, l_verbose)
    printlog ("filterdb   = "//l_filterdb, l_logfile, l_verbose)
    printlog ("key_dispaxis = "//l_key_dispaxis, l_logfile, l_verbose)
    printlog ("dispaxis   = "//l_dispaxis, l_logfile, l_verbose)
    printlog ("sci_ext    = "//l_sci_ext, l_logfile, l_verbose)
    if (l_fl_vardq) { 
        printlog ("var_ext = "//l_var_ext, l_logfile, l_verbose)
        printlog ("dq_ext  = "//l_dq_ext, l_logfile, l_verbose)
    } 
    printlog ("", l_logfile, l_verbose)

    # Check if the l_database has variable has a trailing slash
    # If so, remove it - it's hardcoded later to add it back in.
    if (substr(l_database,(strlen(l_database)),strlen(l_database)) == "/") {
        l_database = substr(l_database,1,(strlen(l_database)-1))
    }

    #
    # Join lists together. The logic here is a tad confusing but it generates 
    # a structured filelist which includes the input file, outputfile name 
    # (MEF), and the transformation database(s) to use for each image in the 
    # input list - all on one line per image to be transformed. 
    #

    # Make list file for main loop
    joinlines (temp1//","//temp2//","//sdistlist//","//wavlist, output=temp3, \
        delim=" ", missing="Missing", maxchar=161, shortest=yes, verbose=no)
    scanfile = temp3 
    while (fscan(scanfile,img,outimg,scal,wcal) != EOF) { 
        # check .fits
        suf = substr (img, strlen(img)-4, strlen(img))
        if (suf != ".fits")
            img = img//".fits"

        suf = substr (outimg, strlen(outimg)-4, strlen(outimg))
        if (suf!=".fits")
            outimg = outimg//".fits"

        suf = substr(wcal, strlen(wcal)-4, strlen(wcal))
        if (suf == ".fits")
            wcal = substr(wcal, 1, strlen(wcal)-5)

        suf = substr(scal, strlen(scal)-4, strlen(scal))
        if (suf == ".fits")
            scal = substr(scal, 1, strlen(scal)-5)

        # copy MDF
        mdf = mktemp("tmpmdf")
        tcopy (img//"[MDF]", mdf//".fits", verbose-)
        wmef (mdf//".fits", outimg, extnames="MDF", verbose-, phu=img,
            >& "dev$null")
        if (access(outimg)) {
            outmade = yes
        }
        imdelete (mdf, verify-, >& "dev$null")

        # get obsmode
        imgets (img//"[0]", "OBSMODE")
        obsmode = imgets.value
        
        # check that GSAPPWAVE has been run
        imgets (img//"[0]", "GSAPPWAVE", >& "dev$null")
        if (imgets.value == "0") {
            gsappwave (img, gratingdb=l_gratingdb, filterdb=l_filterdb,
                key_dispaxis=l_key_dispaxis, dispaxis=l_dispaxis,
                logfile=l_logfile, verbose=l_verbose)
        }

        # get number of rows in MDF
        imgets (img//"[0]", "NSCIEXT")
        nsciin = int(imgets.value)
        # Loop over extensions
        nsciout = 0
        for (i=1; i<=nsciin; i+=1) {
            sel = 1
            mdfrow = i
            if (obsmode == "MOS") {
                imgets (img//"["//l_sci_ext//","//i//"]", "MDFROW",
                    >& "dev$null")
                if (imgets.value != "0") {
                    mdfrow = int(imgets.value)
                    tabpar (img//"[MDF]", "SELECT", mdfrow)
                    if (tabpar.undef == no)
                        sel = int(tabpar.value)
                }
            } 
            if (sel==1) {
                #snum = "000"+mdfrow
                printf ("%03d\n", mdfrow) | scan(snum)
                if (obsmode=="IFU") {
                    dbprefix = "id"
                    l_fl_stran = no
                } else
                    dbprefix = "fc"
                    
                scalfile = scal//"_"//snum
                wcalfile = wcal//"_"//snum
                if (l_fl_stran && l_fl_wavtran) {
                    fitname = scalfile//","//wcalfile
                    go = access(l_database//"/"//dbprefix//wcalfile)
                    go = access(l_database//"/"//dbprefix//scalfile)
                } else if (l_fl_stran && !l_fl_wavtran) {
                    fitname = scalfile
                    go = access (l_database//"/"//dbprefix//scalfile)
                } else if (!l_fl_stran && l_fl_wavtran) {
                    fitname = wcalfile
                    go = access (l_database//"/"//dbprefix//wcalfile)
                }
                if (!go) {
                    printlog ("ERROR - GSTRANFORM: Wavelength or S-distortion \
                        calibration file not found", l_logfile, verbose+)
                    printlog ("ERROR - GSTRANFORM: for image "//img//\
                        "["//l_sci_ext//","//i//"]", l_logfile, verbose+)
                    printlog ("ERROR - GSTRANFORM: output image "//outimg//\
                        " will be incomplete", l_logfile, verbose+)
                    printlog ("Continuing with remaining extensions and \
                        images", l_logfile,l_verbose)
                    status = 1
                }

                if (go && imaccess(img//"["//l_sci_ext//","//i//"]")) {
                    nsciout = nsciout+1

                    printlog ("Transforming "//img//"["//l_sci_ext//",\
                        "//i//"]", l_logfile, l_verbose)
                    printlog ("MDF row: "//mdfrow, l_logfile, l_verbose)

                    if (obsmode != "IFU") {                    

                        keypar (img//"["//l_sci_ext//","//i//"]", "i_naxis2", \
                            silent+)
                        if (keypar.found) {
                            imgheight = int(keypar.value)
                        } else {
                            printlog ("ERROR - GSTRANSFORM: Cannot access \
                                NAXIS2 keyword for "//img//"["//l_sci_ext//\
                                ","//i//"]. Exiting.", l_logfile, verbose+)
                            goto error
                        }

                        # Parse out the image dimensions of the wavcal file
                        # from the calibrration database file
                        for (j = 1; j <= 4; j += 1) {
                            coordtmp = ""
                            fields (files=l_database//"/"//dbprefix//wcalfile,\
                                fields="1", lines=str(10 + j), \
                                quit_if_miss=no,print_file_n=no) | \
                                scan (coordtmp)
                            if (coordtmp != "") {
                                # 1=x1; 2=x2; 3=y1, 4=y2
                                pos_wavecal[j] = nint(real(coordtmp))
                            } else {
                                printlog ("ERROR - GSTRANSFORM: Cannot parse \
                                    dimensions of wavelength cailbration \
                                    from database file. Exiting.", \
                                    l_logfile, verbose+)
                                goto error
                            }
                        }
                        wavheight = pos_wavecal[4] - (pos_wavecal[3] - 1)
                        wavlength = pos_wavecal[2] - (pos_wavecal[1] - 1)

                        # Read the naxis keywords for the image to be 
                        # transformed
                        keypar (img//"["//l_sci_ext//","//i//"]", "i_naxis1", \
                            silent+)
                        if (keypar.found) {
                            imglength = int(keypar.value)
                        } else {
                            printlog ("ERROR - GSTRANSFORM: Cannot access \
                                NAXIS1 keyword for "//img//"["//l_sci_ext//\
                                ","//i//"]. Exiting.", l_logfile, verbose+)
                            goto error
                        }
                        # Compare the dimensions of the wavelength calibration
                        # and the image to be transformed
                        if (imgheight != wavheight || imglength != wavlength) {
                            if (imgheight > wavheight || \
                                imglength > wavlength) {
                                # Image to transform is bigger than wavecal
                                size_print = "smaller"
                            } else if (imgheight < wavheight || \
                                imglength < wavlength) {
                                # Image to transform is smaller than wavecal
                                size_print = "larger"
                            }
                            printlog ("ERROR - GSTRANSFORM: The \
                                wavelength calibration image "//\
			        wcal//".fits["//l_sci_ext//","//\
                                str(i)//"]"//\
                                "\n                     appears to \
                                have "//size_print//" dimensions than \
                                the image to be"//\
                                "\n                     transformed "//\
                                img//"["//l_sci_ext//","//i//"].", \
                                l_logfile, verbose+)

                            printlog ("                     Please confirm \
                                that the detector sections (DETSEC keyword)"//\
                                "\n                     for the input "//\
                                l_sci_ext//" extension(s) to be transformed \
                                match the"//\
                                "\n                     detector \
                                section(s) for the corresponding "//\
                                l_sci_ext//" extensions"//\
                                "\n                     from the image \
                                used to create the wavelength calibration.", \
                                l_logfile, verbose+)
                            printlog ("                     If they do not \
                                match, please trim them accordingly. In"//\
                                "\n                     case you need to \
                                trim the image used to create the"//\
                                "\n                     wavelength \
                                calibration image, please re-run \
                                gswavlength on"//\
                                "\n                     that image after \
                                trimming.", \
                                l_logfile, verbose+)
                            printlog ("                     NOTE: DETSEC \
                                refers to raw pixels on the CCDs, i.e., \
                                it is"//\
                                "\n                     not binned. Please \
                                convert into binned space before trimming.", \
                                l_logfile, verbose+)
                            printlog ("                     Exiting.", \
                                l_logfile, verbose+)
                            goto error
                        }
                    }
                                        
                    tmpsci = mktemp("tmpsci") 
                    tmpvar = mktemp("tmpvar") 
                    tmpdq = mktemp("tmpdq") 

                    if (obsmode == "IFU") {
                        # science plane
                        gemhedit (img//"["//l_sci_ext//","//i//"]", "REFSPEC1",
                            wcalfile, "", delete-)
                        dispcor (img//"["//l_sci_ext//","//i//"]", tmpsci,
                            linear+, database=l_database, table="", w1=INDEF,
                            w2=INDEF, dw=INDEF, nw=INDEF, log-, flux=l_flux,
                            samedisp+, global-, ignoreaps-, listonly-,
                            verbose=l_verbose, logfile=l_logfile)
                        if (l_fl_vardq && imaccess(img//"["//l_var_ext//",\
                            "//i//"]")) {
                            # variance plane
                            gemhedit (img//"["//l_var_ext//","//i//"]", \
                                "REFSPEC1", wcalfile, "", delete-)
                            dispcor (img//"["//l_var_ext//","//i//"]", tmpvar,
                                linear+, database=l_database, table="",
                                w1=INDEF, w2=INDEF, dw=INDEF, nw=INDEF, log-,
                                flux=l_flux, samedisp+, global-, ignoreaps-,
                                listonly-, verbose=l_verbose, 
                                logfile=l_logfile)
                            # dq plane
                            gemhedit (img//"["//l_dq_ext//","//i//"]", \
                                "REFSPEC1", wcalfile, "", delete-)
                            dispcor (img//"["//l_dq_ext//","//i//"]", tmpdq,
                                linear+, database=l_database, table="",
                                w1=INDEF, w2=INDEF, dw=INDEF, nw=INDEF, log-,
                                flux=l_flux, samedisp+, global-, ignoreaps-,
                                listonly-, verbose=l_verbose, 
                                logfile=l_logfile)
                        }
                    } else { 
                        # science extension
                        transform (img//"["//l_sci_ext//","//i//"]", tmpsci,
                            fitnames=fitname, database=l_database,
                            interptype=l_interptype, x1=l_lambda1, 
                            x2=l_lambda2, dx=l_dx, nx=l_nx, xlog=l_lambdalog, 
                            y1=INDEF, y2=INDEF, dy=INDEF, ny=INDEF, 
                            ylog=l_ylog, flux=l_flux, 
                            logfiles=l_stdout//l_logfile)
                        if (l_fl_vardq && imaccess(img//"["//l_var_ext//",\
                            "//i//"]")) { 
                            # variance plane
                            transform (img//"["//l_var_ext//","//i//"]", 
                                tmpvar, fitnames=fitname, database=l_database,
                                interptype=l_interptype, x1=l_lambda1, 
                                x2=l_lambda2, dx=l_dx, nx=l_nx, 
                                xlog=l_lambdalog, y1=INDEF, y2=INDEF, 
                                dy=INDEF, ny=INDEF, ylog=l_ylog, flux=l_flux, 
                                logfiles=l_stdout//l_logfile) 
                        }
                        if (l_fl_vardq && imaccess(img//"["//l_dq_ext//",\
                            "//i//"]")) { 
                            # dq plane
                            transform (img//"["//l_dq_ext//","//i//"]", tmpdq,
                                fitnames=fitname, database=l_database,
                                interptype=l_interptype, x1=l_lambda1, 
                                x2=l_lambda2, dx=l_dx, nx=l_nx,  
                                xlog=l_lambdalog, y1=INDEF, y2=INDEF,  
                                dy=INDEF, ny=INDEF, ylog=l_ylog, flux=no,  
                                logfiles=l_stdout//l_logfile) 
                            imreplace (tmpdq, 1, lower=1, upper=INDEF, 
                                radius=0.) 
                        } 
                    }
                    ## Put transformed slit in output file 
                    imgets (outimg//"[0]", "NEXTEND")
                    nextnd = int(imgets.value)
                    fxinsert (tmpsci//".fits",outimg//"["//nextnd//"]", "0",
                        verbose-, >& "dev$null")
                    nextnd = nextnd+1
                    gemhedit (outimg//"[0]", "NEXTEND", nextnd, "", delete-)
                    gemhedit (outimg//"["//nextnd//"]", "EXTNAME", l_sci_ext, 
                        "", delete-)
                    gemhedit (outimg//"["//nextnd//"]", "EXTVER", nsciout, "", 
                        delete-)
                    if (l_fl_stran)
                        gemhedit (outimg//"["//nextnd//"]", "SDISTNAM",
                            scalfile, "Name of S-distortion", delete-)
                    if (l_fl_wavtran)
                        gemhedit (outimg//"["//nextnd//"]", "WAVTRAN", 
                            wcalfile, "Name of wavelength transformation", 
                            delete-)
                    # Variance and DQ
                    if (imaccess(tmpvar)) {
                        fxinsert (tmpvar//".fits",
                            outimg//"["//nextnd//"]", "0", verbose-,
                            >& "dev$null")
                            
                        # Variance header
                        nextnd = nextnd + 1
                        gemhedit (outimg//"[0]", "NEXTEND", nextnd, "", 
                            delete-)
                        gemhedit (outimg//"["//nextnd//"]", "EXTNAME", \
                            l_var_ext, "", delete-)
                        gemhedit (outimg//"["//nextnd//"]", "EXTVER", \
                            nsciout, "", delete-)
                        if (l_fl_stran)
                            gemhedit (outimg//"["//nextnd//"]", "SDISTNAM",
                                scalfile, "Name of S-distortion", delete-)
                        if (l_fl_wavtran)
                            gemhedit (outimg//"["//nextnd//"]", "WAVTRAN",
                                wcalfile, "Name of wavelength transformation", 
                                delete-)
                    }
                    if (imaccess(tmpdq)) {
                        fxinsert (tmpdq//".fits",
                            outimg//"["//nextnd//"]", "0", verbose-,
                            >& "dev$null")
                        
                        # DQ header
                        nextnd = nextnd + 1
                        gemhedit (outimg//"[0]", "NEXTEND", nextnd, "", 
                            delete-)
                        gemhedit (outimg//"["//nextnd//"]", "EXTNAME", \
                            l_dq_ext, "", delete-)
                        gemhedit (outimg//"["//nextnd//"]", "EXTVER", \
                            nsciout, "", delete-)
                        if (l_fl_stran)
                            gemhedit (outimg//"["//nextnd//"]", "SDISTNAM",
                                scalfile, "Name of S-distortion", delete-)
                        if (l_fl_wavtran)
                            gemhedit (outimg//"["//nextnd//"]", "WAVTRAN",
                                wcalfile, "Name of wavelength transformation", 
                                delete-)
                    }
                    imdelete (tmpsci//","//tmpvar//","//tmpdq, verify-,
                        >& "dev$null")
                }
            }
        } # end loop over extensions
        gemhedit (outimg//"[0]", "NSCIEXT", nsciout, "", delete-)
               
        gemdate ()
        gemhedit (outimg//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
        gemhedit (outimg//"[0]", "GSTRANSF", gemdate.outdate,
            "UT Time stamp for GSTRANSFORM", delete-)
    } #end of while loop over input image list
    scanfile = ""

    goto clean

error:
    status = 1
    if (outmade) {
        if (imaccess(outimg)) {
            imdelete (outimg, verify=no, >& "dev$null")
        }
    }
    goto clean

clean: 
    #Clean Up 
    printlog ("", l_logfile, l_verbose)
    scanfile = ""
    if (status==0)
        printlog ("GSTRANSFORM exit status: good", l_logfile, l_verbose)
    else
        printlog ("GSTRANSFORM exit status: error", l_logfile, l_verbose)

    printlog ("-------------------------------------------------------------\
        ------------------", l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)
    delete (temp1//","//temp2//","//temp3//","//wavlist//","//sdistlist,
        verify-, >& "dev$null") 

end
