# Copyright(c) 2002-2016 Association of Universities for Research in Astronomy, Inc.

procedure gscalibrate (input)

# Apply the (relative) flux calibration to GMOS spectra. All input
# spectra must be processed with GSREDUCE and GSTRANSFORM. It is
# assumed that the spectra could have been processed with GSSKYSUB,
# but this is not a requirement. The spectra may or may not have
# been processed with GSEXTRACT. The spectra may be 1D or 2D.
#
# Version  Feb 28, 2002   RC  v1.3 release
#          May 23, 2002   BM  update for IFU
#          Sept 20, 2002         v1.4 release
#          Aug 26, 2003   KL  IRAF2.12 - new parameter, addonly, in hedit
#          Oct 28, 2005   BM  fix imexpr calls for DQ/VAR propagation

string  input       {prompt="Input spectra to calibrate"}
string  output      {"",prompt="Output calibrated spectra"}
string  outpref     {"c",prompt="Output prefix"}
string  sfunction   {"sens",prompt="Input image root name for sensitivity function"}
string  sci_ext     {"SCI",prompt="Name of science extension"}
string  var_ext     {"VAR",prompt="Name of variance extension"}
string  dq_ext      {"DQ",prompt="Name of data quality extension"}
string  key_airmass {"AIRMASS",prompt="Airmass header keyword"}
string  key_exptime {"EXPTIME",prompt="Exposure time header keyword"}
bool    fl_vardq    {no,prompt="Propagate VAR/DQ planes"}
bool    fl_ext      {no,prompt="Apply extinction correction to input spectra"}
bool    fl_flux     {yes,prompt="Apply flux calibration to input spectra"}
bool    fl_scale    {yes,prompt="Multiply output with fluxscale"}
real    fluxscale   {1.0E15,prompt="Value of the flux scale (fl_scale=yes)"}
bool    ignoreaps   {yes,prompt="Ignore aperture numbers in flux calibration"}
bool    fl_fnu      {no,prompt="Create spectra having units of FNU"}
string  extinction  {"",prompt="Extinction file"}
string  observatory {"Gemini-North",prompt="Observatory"}
string  logfile     {"",prompt="Logfile name"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}
string  *scanfile   {prompt="For internal use only"}

begin

    # Local variable definitions
    string  l_input, l_output, l_prefix, l_sfunction, l_sci_ext
    string  l_var_ext, l_dq_ext, l_extinction, l_observatory, l_logfile
    string  l_key_airmass, l_key_exptime
    real    l_fluxscale
    bool    l_fl_ext, l_fl_flux, l_fl_fnu, l_ignoreaps,l_fl_vardq
    bool    l_verbose,l_fl_scale

    # Other variables
    string  tmpsci, tmpvar, tmpdq, tmpin, tmpout, mdffile, tmpnoise
    string  joinlst, suf, img, oimg, valexp, inlst, outlst, valpar
    string  inimg[200], outimg[200], obsmode[200], observat[200]
    string  l_key_qecorrim, l_key_qestate
    int     nerror, ii, j, i, ninpimg, noutimg, numext, nima, nextnd
    int     mdfpos[200], nsciext[200], msktype
    real    airmass[200], exptime[200]
    bool    pref, mdf, sfunc_qestate, inqecorr_state
    struct  sdate

    # Query parameters
    l_input = input; l_output=output; l_prefix=outpref
    l_sci_ext = sci_ext; l_var_ext=var_ext; l_dq_ext=dq_ext
    l_sfunction=sfunction; l_key_airmass=key_airmass;
    l_key_exptime=key_exptime; l_fl_ext=fl_ext; l_fl_flux=fl_flux
    l_ignoreaps=ignoreaps; l_fl_fnu=fl_fnu; l_extinction=extinction
    l_observatory=observatory; l_verbose=verbose
    l_logfile=logfile; l_fl_vardq=fl_vardq; l_fluxscale=fluxscale
    l_fl_scale=fl_scale

    status = 0
    sfunc_qestate = no
    inqecorr_state = no
    l_key_qecorrim = "QECORRIM"
    l_key_qestate = "QESTATE"

    # Cache some important tasks
    cache ("imgets", "gimverify", "gextverify", "gemdate")

    # define some temporary files
    tmpin = mktemp("tmpin")
    tmpout = mktemp("tmpout")
    joinlst = mktemp("tmpjoin")

    #
    # Checking procedure.
    #
    # First of all: check the logfile
    #
    # Check that logfile is not an empty string, otherwise,
    #   use default: gmos.log
    #
    if (l_logfile == "STDOUT") {
        l_logfile = ""
        l_verbose = yes
    } else if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GSCALIBRATE: Both gscalibrate.log and \
                gmos.logfile are empty.", l_logfile, verbose+)
            printlog ("                       Using default file gmos.log.",
                l_logfile, verbose+)
        }
    }

    # Setup date
    date | scan(sdate)

    # Logfile: what will be done:
    printlog ("-------------------------------------------------------------\
        ---------------", l_logfile,verbose+)
    printlog ("GSCALIBRATE -- "//sdate, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)
    printlog ("input        = "//l_input, l_logfile, l_verbose)
    printlog ("output       = "//l_output, l_logfile, l_verbose)
    printlog ("outpref      = "//l_prefix, l_logfile, l_verbose)
    printlog ("sfunction    = "//l_sfunction, l_logfile, l_verbose)
    if (l_fluxscale > 1E8) {
        printf ("fluxscale    = %7.3g\n", l_fluxscale) | scan(sdate)
        printlog (sdate, l_logfile, l_verbose)
    } else
        printlog ("fluxscale    = "//l_fluxscale, l_logfile, l_verbose)

    printlog ("fl_ext       = "//l_fl_ext, l_logfile, l_verbose)
    printlog ("fl_flux      = "//l_fl_flux, l_logfile, l_verbose)
    printlog ("fl_vardq     = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("sci_ext      = "//l_sci_ext, l_logfile, l_verbose)
    if (l_fl_vardq) {
        printlog ("var_ext      = "//l_var_ext, l_logfile, l_verbose)
        printlog ("dq_ext       = "//l_dq_ext, l_logfile, l_verbose)
    }
    printlog ("ignoreaps    = "//l_ignoreaps, l_logfile, l_verbose)
    printlog ("fl_fnu       = "//l_fl_fnu, l_logfile, l_verbose)
    printlog ("extinction   = "//l_extinction, l_logfile, l_verbose)
    printlog ("key_airmass  = "//l_key_airmass, l_logfile, l_verbose)
    printlog ("key_exptime  = "//l_key_exptime, l_logfile, l_verbose)
    printlog ("observatory  = "//l_observatory, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    # Check that input and output image is not an empty string. Check that
    # it exists.
    nerror = 0

    # Check if the input file is not empty string
    if ((l_input=="") || (l_input==" ")) {
        printlog ("ERROR - GSCALIBRATE: input image(s) or list is not \
            specified", l_logfile, verbose+)
        nerror = nerror+1
    }

    # Check existence of input list

    if (strstr("@",l_input) != 0) {
        inlst = substr(l_input,stridx("@",l_input)+1,strlen(l_input))
        if (!access(inlst)) {
            printlog ("ERROR - GSCALIBRATE: Input list "//inlst//" not found",
                l_logfile, verbose+)
            nerror = nerror+1
        }
    }

    if ((l_output=="") || (l_output==" "))
        pref = yes
    else if ((l_output!="") || (l_output!=" "))
        pref = no

    if (pref) {
        if ((l_prefix=="") || (l_prefix==" ")) {
            printlog ("ERROR - GSCALIBRATE: outprefix is not specified.",
                l_logfile, verbose+)
            nerror = nerror+1
        } else if (substr(l_output,1,1)=="@") {
            outlst = substr (l_output,2,strlen(l_output))
            if (!access(outlst)) {
                printlog ("ERROR - GSCALIBRATE: Output list "//outlst//\
                    " not found", l_logfile, verbose+)
                nerror = nerror+1
            }
        }
    }

    # Check if input images exists and count the number of input images
    if (strstr("@",l_input) != 0) {
        sections (l_input, > tmpin)
    } else
        files (l_input, sort-, > tmpin)

    # Count the number of input images
    count (tmpin) | scan (ninpimg)
    scanfile = tmpin
    while (fscan(scanfile,img) !=EOF) {
        gimverify (img)
        if (gimverify.status>0) {
            printlog("ERROR - GSCALIBRATE: Input image "//img//" does not \
                exist.", l_logfile, verbose+)
            nerror = nerror+1
        }
    }

    scanfile = ""

    # Check if outimages exist and if the the number of input images are
    # equal to the number of output images
    if (!pref) {
        if (substr(l_output,1,1)=="@") {
            outlst = substr(l_output, 2, strlen(l_output))
            sections ("@"//outlst, > tmpout)
        } else
            files (l_output, sort-, > tmpout)

        count (tmpout) | scan(noutimg)
        scanfile = tmpout
        if (ninpimg != noutimg) {
            printlog ("ERROR - GSCALIBRATE: Number of input and output \
                images does not match.", l_logfile, verbose+)
            nerror = nerror+1
        }
        while (fscan(scanfile,img) !=EOF) {
            if (imaccess(img)) {
                printlog ("ERROR - GSCALIBRATE: Output image "//img//\
                    " exists.", l_logfile, verbose+)
                nerror = nerror+1
            }
        }
    }
    if (pref) {
        sections (l_prefix//"@"//tmpin, > tmpout)
        scanfile = tmpout
        while (fscan(scanfile,img) !=EOF) {
            if (imaccess(img)) {
                printlog ("ERROR - GSCALIBRATE: Output image "//img//" exist.",
                    l_logfile, verbose+)
                nerror = nerror+1
            }
        }
    } # End check outprefix and outimages

    # Check if the sensitivity function is not an empty string
    if ((l_sfunction == "") || (l_sfunction == " ")) {
        printlog ("ERROR - GSCALIBRATE : sensitivity function is not \
            specified", l_logfile, l_verbose)
        nerror = nerror+1
    } else {
        # Read the QE state of the sensitivity function
        keypar (l_sfunction//"[0]", l_key_qestate, silent+)
        if (keypar.found) {
            if (keypar.value == "yes") {
                sfunc_qestate = yes
            } else if (keypar.value == "no") {
                sfunc_qestate = no
            } else {
                printlog ("ERROR - GSCALIBRATE: keyword "//l_key_qestate//\
                    " is neither 'yes' or 'no'.", \
                    l_logfile, verbose+)
                nerror += 1
            }
        } else {
            printlog ("WARNING - GSCALIBRATE: Cannot find the keyword "//\
                l_key_qestate//" in "//l_sfunction//".\n"//\
                "                       Assuming the data "//\
                "used to create the sensitivity function \n"//\
                "                       "//l_sfunction//\
                " was not QE corrected.", \
                l_logfile, verbose+)
            sfunc_qestate = no
        }
    }

    # Check AIRMASS, EXPTIME and OBSERVATORY parameters
    if ((l_key_airmass== "") || (l_key_airmass == " ")) {
        printlog ("ERROR - GSCALIBRATE : Airmass image header parameter name \
            is missing.", l_logfile, verbose+)
        nerror = nerror+1
    }
    if ((l_key_exptime == "") || (l_key_exptime == " ")) {
        printlog ("ERROR - GSCALIBRATE : Exptime image header parameter name \
            is missing.", l_logfile, verbose+)
        nerror = nerror+1
    }
    if (l_observatory == "" || l_observatory == " ") {
        printlog ("ERROR - GSCALIBRATE: parameter observatory is not \
            specified", l_logfile, verbose+)
        nerror = nerror+1
    }

    # Verify that the extension name SCI is not an empty string.
    if (l_sci_ext=="" || l_sci_ext==" ") {
        printlog ("ERROR - GSCALIBRATE: sci_ext is an empty string.",
            l_logfile, verbose+)
        nerror = nerror+1
    }

    # If fl_vardq=yes, verify that the extension names VAR and DQ are not
    # an empty strings
    if (l_fl_vardq) {
        if (l_var_ext=="" || l_var_ext ==" ") {
            printlog ("ERROR - GSCALIBRATE: var_ext is an empty string.",
                l_logfile, verbose+)
            nerror = nerror+1
        }
        if (l_dq_ext=="" || l_dq_ext ==" ") {
            printlog ("ERROR - GSCALIBRATE: dq_ext is an empty string.",
                l_logfile, verbose+)
            nerror = nerror+1
        }
    }

    # End of first checks. # If nerror!=0, go outside
    if (nerror > 0)
        goto outerror

    scanfile = ""
    nerror = 0

    # Check each input image
    scanfile = tmpin
    i = 0
    mdf = no
    while (fscan(scanfile,img) !=EOF) {
        i+=1
        suf = substr(img,strlen(img)-3,strlen(img))
        if (suf!="fits")
            img = img//".fits"

        # Check if the images are MEF
        imgets (img//"[0]", "EXTEND", >& "dev$null")
        if (imgets.value=="F" || imgets.value=="0") {
            printlog ("ERROR - GSCALIBRATE: image "//img//\
                " is not a MEF file.",
                l_logfile, verbose+)
            nerror = nerror+1
        }

        # Check image type
        imgets (img//"[0]","MASKTYP", >& "dev$null")
        msktype = int(imgets.value)
        if (msktype != 1 && msktype != -1) {
            printlog ("ERROR - GSCALIBRATE: "//img//" has MASKTYP other than \
                LONGSLIT , MOS, or IFU mode.", l_logfile, verbose+)
            nerror = nerror+1
        }
        # Check the MDF file.  Check where is it.
        imgets (img//"[0]", "NEXTEND", >& "dev$null")
        numext = int(imgets.value)
        for (k=1;k<=numext;k+=1) {
            keypar (img//"["//k//"]", keyword="EXTNAME", silent=yes)
            if (keypar.value == "MDF" && keypar.found) {
                mdf = yes
                mdfpos[i] = k
            }
        }
        if (!mdf) {
            printlog ("ERROR - GSCALIBRATE: Input file "//img//" does not \
                have an attached MDF table.", l_logfile, verbose+)
            nerror = nerror+1
        }
        # Get obsmode
        imgets (img//"[0]","OBSMODE",>& "dev$null")
        obsmode[i] = imgets.value
        if (imgets.value == "" || imgets.value == " " || imgets.value == "0") {
            printlog ("ERROR - GSCALIBRATE: could not find the OBSMODE \
                keyword parameter in "//img, l_logfile, verbose+)
            nerror = nerror+1
        }
        if (obsmode[i] != "LONGSLIT" && obsmode[i] != "MOS" && \
            obsmode[i] != "IFU") {
                printlog ("ERROR - GSCALIBRATE: "//img//" has OBSMODE other"//\
                    "than LONGSLIT, MOS, or IFU mode.", l_logfile, verbose+)
            nerror = nerror+1
        }
        # Check how many SCI extension we have
        imgets (img//"[0]","NSCIEXT", >& "dev$null")
        nsciext[i] = int(imgets.value)
        if (nsciext[i] == 0) {
            printlog ("ERROR - GSCALIBRATE: Number of SCI extensions unknown \
                for image "//img, l_logfile, verbose+)
            nerror = nerror+1
        } else {
            for (j=1;j<=nsciext[i];j+=1) {
                if (!imaccess(img//"["//l_sci_ext//","//str(j)//"]")) {
                    printlog ("ERROR - GSCALIBRATE: Could not access "//img//\
                        "["//l_sci_ext//","//str(j)//"]", l_logfile, verbose+)
                    nerror = nerror+1
                } else {
                    # Check if the MDFROW keyword exist inside each extension
                    # (MOS only)
                    if (obsmode[i] == "MOS") {
                        imgets (img//"["//l_sci_ext//","//str(j)//"]",
                            "MDFROW", >& "dev$null")
                        if (imgets.value == "" || imgets.value == " " || \
                            imgets.value == "0") {
                                printlog ("ERROR - GSCALIBRATE: could not "//\
                                    "find the MDFROW keyword parameter in "//\
                                    img//"["//l_sci_ext//","//str(j)//"]", \
                                    l_logfile, verbose+)
                            nerror = nerror+1
                        }
                    }
                    # If fl_vardq=yes, check if there are access the images
                    if (l_fl_vardq) {
                        if (!imaccess(img//"["//l_var_ext//","//str(j)//"]")) {
                            printlog ("WARNING - GSCALIBRATE: Could not \
                                access "//img//"["//l_var_ext//","//\
                                str(j)//"]", l_logfile, verbose+)
                            l_fl_vardq=no
                        }
                        if (!imaccess(img//"["//l_dq_ext//","//str(j)//"]")) {
                            printlog ("WARNING - GSCALIBRATE: Could not \
                                access "//img//"["//l_dq_ext//","//j//"]",
                                l_logfile, verbose+)
                            l_fl_vardq=no
                        }
                        if (!l_fl_vardq)
                            printlog ("                       "//\
                                      "fl_vardq disabled", l_logfile, verbose+)
                    }
                }
            }
        }

        # Check AIRMASS and EXPTIME keyword inside the header
        imgets (img//"[0]",l_key_airmass,>&"dev$null")
        if (imgets.value=="" || imgets.value==" " || imgets.value=="0") {
            printlog ("ERROR - GSCALIBRATE: Image header parameter "//\
                l_key_airmass//" not found in "//img, l_logfile, verbose+)
            nerror = nerror+1
        } else {
            airmass[i] = real(imgets.value)
            if (airmass[i]==0.0) {
                printlog ("WARNING - GSCALIBRATE: "//l_key_airmass//\
                    " for "//img//" is 0.0 Using airmass=1.", l_logfile,
                    verbose+)
                airmass[i] = 1.0
            }
        }

        imgets (img//"[0]", l_key_exptime, >&"dev$null")
        if (imgets.value=="" || imgets.value==" " || imgets.value=="0") {
            printlog ("ERROR - GSCALIBRATE: Image header parameter "//\
                l_key_exptime//" not found in "//img, l_logfile, verbose+)
            nerror = nerror+1
        } else
            exptime[i] = real(imgets.value)

        imgets (img//"[0]","OBSERVAT",>&"dev$null")
        if (imgets.value=="" || imgets.value==" " || imgets.value=="0")
            observat[i] = l_observatory
        else
            observat[i] = imgets.value

        # Check if the images where GSTRANSFORMED. If not, out!
        if (obsmode[i] == "IFU") {
            imgets (img//"[0]", "GFTRANSF", >& "dev$null")
            if (imgets.value == "" || imgets.value == " " || \
                imgets.value=="0") {
                    printlog ("ERROR - GSCALIBRATE: Image "//img//\
                        " has not been transformed. Run GFTRANSFORM first.",\
                        l_logfile, verbose+)
                    nerror = nerror+1
            }
        } else {
            imgets (img//"[0]", "GSTRANSF", >& "dev$null")
            if (imgets.value == "" || imgets.value == " " || \
                imgets.value=="0") {
                    printlog ("ERROR - GSCALIBRATE: Image "//img//\
                        " has not been transformed. Run GSTRANSFORM first.", \
                        l_logfile, verbose+)
                    nerror = nerror+1
            }
        }

        # Read the QE state of the input image
        keypar (img//"[0]", l_key_qecorrim, silent+)
        if (keypar.found) {
            # Record the QE state of the inout images (default is no)
            inqecorr_state = yes
        } else {
            inqecorr_state = no
        }

        # Compare the QE state of the input image to the sensitivity function
        if (inqecorr_state != sfunc_qestate) {
            printlog ("ERROR - GSCALIBRATE: QE correction states of "//\
                img//" and "//l_sfunction//" do not match.", \
                l_logfile, verbose+)
            nerror += 1
        }

    } # End the big loop over the images

    if (nerror > 0) goto outerror

    joinlines (tmpin//","//tmpout, output=joinlst, delim=" ", \
        missing="Missing", maxchars=161, shortest-, verbose-)
    delete (tmpin//","//tmpout, verify-, >& "dev$null")

    scanfile = joinlst

    i = 0
    while (fscan(scanfile,img,oimg) != EOF) {
        i+=1
        inimg[i] = img
        outimg[i] = oimg
        suf = substr(inimg[i],strlen(inimg[i])-3,strlen(inimg[i]))
        if (suf != "fits")
            inimg[i] = inimg[i]//".fits"

        suf = substr(outimg[i],strlen(outimg[i])-3,strlen(outimg[i]))
        if (suf != "fits")
            outimg[i]=outimg[i]//".fits"
    }
    delete (joinlst, verify-, >& "dev$null")

    nima = i
    i = 0

    #Setup SPECRED defaults which should not be changed by the outside
    specred.dispaxis = 1
    specred.nsum = 10
    specred.interp = "poly5"
    specred.verbose = l_verbose
    specred.logfile = l_logfile

    for (j=1; j<=nima; j+=1) {
        # create tmp FITS file names used only within this loop.
        mdffile = mktemp("tmpmfd")

        printlog ("", l_logfile, l_verbose)
        #
        # copy MDF to outfile
        #
        tcopy (inimg[j]//"["//mdfpos[j]//"]", mdffile//".fits", verbose-)
        wmef (mdffile//".fits", outimg[j], extnames="MDF", verbose-, \
            phu=inimg[j], >& "dev$null")

        # loop over the number of SCI extensions for each image. If LONGSLIT,
        # then nsciext=1.  Each image is extinction corrected, if needed, and
        # flux calibrated.

        printlog ("Input image : "//inimg[j], l_logfile, l_verbose)
        printlog ("Output image : "//outimg[j], l_logfile, l_verbose)
        printlog ("Number of science extension (slits) : "//nsciext[j],
            l_logfile, l_verbose)
        printlog (" ", l_logfile, l_verbose)

        for (ii=1;ii<=nsciext[j];ii+=1) {
            # Create tmp FITS file name used only within this loop
            tmpsci = mktemp("tmpsci")
            if (l_fl_vardq) {
                tmpvar = mktemp("tmpvar")
                tmpdq = mktemp("tmpdq")
            }
            tmpnoise = mktemp("tmpnoise")

            printlog ("Slit #    "//ii, l_logfile, l_verbose)
            specred.calibrate (inimg[j]//"["//l_sci_ext//","//str(ii)//"]",
                tmpsci, airmass[j], exptime[j], extinct=l_fl_ext,
                flux=l_fl_flux, extinction=l_extinction,
                observatory=observat[j], ignoreaps=l_ignoreaps,
                sensitivity=l_sfunction, fnu=l_fl_fnu, >> l_logfile)

            # Update the airmass keyword in the science and variance extensions
            # specred.standard writes them to the science extension; which gets
            # used by specred.calibrate. If missing from the header AIRMASS
            # gets set to INDEF which is not FITS standard - MS
            # 2015-01-08 - MS: Comment out the call to update the input - IRAF
            # doesn't appear to add the AIRMASS keyword to the input when only
            # calibrating the input; only update it for the output.
#            gemhedit (inimg[j]//"["//l_sci_ext//","//j//"]", l_key_airmass, \
#                airmass[j], "Mean airmass for the observation", delete-, \
#                upfile="")

            # Repeat for output file
            gemhedit (tmpsci, l_key_airmass, \
                airmass[j], "Mean airmass for the observation", delete-, \
                upfile="")

            # Optional multiplication by the value given in l_fluxscale
            # (default 1E15) to avoid floating point errors in using the
            # spectra. Using imarith.
            if (l_fl_scale) {
                imarith (tmpsci,"*", l_fluxscale, tmpsci, verbose-,
                    >& "dev$null")
            } else {
                l_fluxscale=1.
            }
            # Make var plane of the spectra (if fl_vardq=yes)
            #
            # JT: the old calculation headed "copied from nscalibrate" is
            # replaced below with a more direct analogue of how the SCI is
            # processed, since it wasn't working properly for multiple
            # reasons:
            # - Imexpr was combining some extreme numerical values to produce
            #   an output spectrum with all zeroes (at least for IFU data).
            # - Sarith was crashing on IFU data due to an IRAF header buffer
            #   overflow with > 1024 APID keywords (fixed in Ureka).
            # - After fixing the IRAF crash, when multiplying by an (erroneous)
            #   spectrum with all zeros, sarith was reproducing its input
            #   (ie. *1.0) for some reason, such that the final variance ended
            #   up the same as the input variance while the SCI had been
            #   scaled by orders of magnitude into flux units.
            # - Also, it was missing one or two second-order terms like
            #   extinction that go into the SCI calculation.
            if (l_fl_vardq) {
                # Propagate data quality extension if fl_vardq=yes (copy the
                # original DQ to a temp file)
                imcopy (inimg[j]//"["//l_dq_ext//","//str(ii)//"]", tmpdq,
                    verbose-)

                # Convert variance to noise (and scale as for the SCI data):
                imexpr("a*sqrt(b)", tmpnoise, l_fluxscale,
                    inimg[j]//"["//l_var_ext//","//str(ii)//"]", dims="auto",
                    intype="auto", outtype="real", refim="auto", bwidth=0,
                    btype="nearest", bpixval=0., rangecheck=yes, verbose=no,
                    exprdb="none", lastout="", mode="ql")

                # Perform the same operation as on the SCI data:
                specred.calibrate (tmpnoise, tmpnoise, airmass[j], exptime[j],
                    extinct=l_fl_ext, flux=l_fl_flux, extinction=l_extinction,
                    observatory=observat[j], ignoreaps=l_ignoreaps,
                    sensitivity=l_sfunction, fnu=l_fl_fnu, >> l_logfile)

                # Convert back to variance
                imexpr("a*a", tmpvar, tmpnoise, dims="auto", intype="auto",
                    outtype="real", refim="auto", bwidth=0, btype="nearest",
                    bpixval=0., rangecheck=yes, verbose=no, exprdb="none",
                    lastout="", mode="ql")

                imdelete (tmpnoise, verify-, >& "dev$null")
            }
            #
            # Insert flux calibrated spectra into outfile
            #
            imgets (outimg[j]//"[0]", "NEXTEND")
            nextnd = int(imgets.value)
            if (l_fl_vardq) {
                fxinsert (tmpsci//".fits,"//tmpvar//".fits,"//tmpdq//".fits",
                    outimg[j]//"["//nextnd//"]", "0", verbose-, >& "dev$null")
                # update the header
                gemhedit (outimg[j]//"[0]", "NEXTEND", (nextnd+3), "",
                    delete-)
                gemhedit (outimg[j]//"["//(nextnd+1)//"]", "EXTNAME", \
                    l_sci_ext, "", delete-)
                gemhedit (outimg[j]//"["//(nextnd+1)//"]", "EXTVER", ii,
                    "", delete-)
                gemhedit (outimg[j]//"["//(nextnd+1)//"]", "FLUXSCAL",
                    fluxscale, "", delete-)
                gemhedit (outimg[j]//"["//(nextnd+2)//"]", "EXTNAME", \
                    l_var_ext, "", delete-)
                gemhedit (outimg[j]//"["//(nextnd+2)//"]", "EXTVER", ii,
                    "", delete-)
                gemhedit (outimg[j]//"["//(nextnd+3)//"]", "EXTNAME", l_dq_ext,
                    "", delete-)
                gemhedit (outimg[j]//"["//(nextnd+3)//"]", "EXTVER", ii,
                    "", delete-)
                gemhedit (outimg[j]//"["//l_var_ext//","//ii//"]", \
                    l_key_airmass, airmass[j], \
                    "Mean airmass for the observation", \
                    delete-, upfile="")
                gemhedit (outimg[j]//"["//l_dq_ext//","//ii//"]", \
                    l_key_airmass, airmass[j], \
                    "Mean airmass for the observation", \
                    delete-, upfile="")

                imdelete (tmpsci//".fits,"//tmpvar//".fits,"//tmpdq//".fits",
                    verify-, >& "dev$null")
            } else {
                fxinsert (tmpsci//".fits", outimg[j]//"["//nextnd//"]", "0",
                    verbose-, >& "dev$null")
                gemhedit (outimg[j]//"[0]", "NEXTEND", (nextnd+1), "",
                    delete-)
                gemhedit (outimg[j]//"["//(nextnd+1)//"]", "EXTNAME", \
                    l_sci_ext, "", delete-)
                gemhedit (outimg[j]//"["//(nextnd+1)//"]", "EXTVER", ii,
                    "", delete-)
                imdelete (tmpsci//".fits", verify-, >& "dev$null")
            }
        }
        # Update the final header
        gemdate ()
        gemhedit (outimg[j]//"[0]", "GSCALIBR", gemdate.outdate,
            "UT Time stamp for GSCALIBRATE", delete-)
        gemhedit (outimg[j]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
        gemhedit (outimg[j]//"[0]", "NSCIEXT", nsciext[j], "", delete-)
        imdelete (mdffile//".fits", verify-, >& "dev$null")
        printlog ("", l_logfile, l_verbose)
    } # End of big while loop

    status = 0
    goto clean

outerror:
    # Exit with error
    delete (tmpin//","//tmpout, verify-, >& "dev$null")
    status = 1
    goto clean

clean:
    # clean up
    scanfile = ""
    date | scan (sdate)
    if (status > 0) {
        printlog ("ERROR - GSCALIBRATE: "//nerror//" error(s) found. Exiting.",
            l_logfile, verbose+)
        printlog ("", l_logfile, verbose+)
        printlog ("GSCALIBRATE done. Exit status bad. --"//sdate,
            l_logfile,verbose+)
    } else {
        printlog ("GSCALIBRATE done. Exit status good --"//sdate,
            l_logfile, verbose+)
    }
    printlog ("------------------------------------------------------------\
        -------------------", l_logfile, l_verbose)
    printlog ("", l_logfile, verbose+)

end
