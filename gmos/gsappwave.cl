# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gsappwave (inimages)

# Wavelength calibrate GMOS spectra
#
# Version   Feb 28, 2002  BM,IJ   v1.3 release
#           Sept 20, 2002         v1.4 release
#           Oct 14, 2002  IJ      don't use ijk from within script
#           Mar 20, 2003  BM      pixel scale for both instruments
#           May 9, 2003   IJ      change in instrument logic, support for old GMOS-N data
#           Aug 26, 2003  KL      IRAF2.12 - new parameter, addonly, in hedit
#           Mar 26, 2004  BM      2 filter support like gscut/gfextract
#           Sept 30, 2008 JT      Print ref. pixel & range correctly for the IFU

string  inimages    {prompt="Input images"}
string  gratingdb   {"gmos$data/GMOSgratings.dat",prompt="Gratings database file"}
string  filterdb    {"gmos$data/GMOSfilters.dat",prompt="Filters database file"}
string  key_dispaxis    {"DISPAXIS",prompt="Keyword for dispersion axis"}
int     dispaxis    {1,min=1,max=2,prompt="Dispersion axis"}
string  logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose output?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

begin

    string  l_inimages, l_logfile, l_filterdb, l_gratingdb, l_linelist
    string  l_grating, l_filter[2], l_key_dispaxis
    string  obsmode, infiles, inlist, img, mdf
    string  gname, fname, ffile, gscut, gfextract, gmosaic
    string  dettype, sci_extn, var_extn, dq_extn
    real    l_wave1, wave1, wave2, fwave1, fwave2, gwave1, gwave2, wmin, wmax
    real    cwave, gblaze, gR, gcoverage, nmppx, scale[2,3], pi
    real    tilt, gtilt, a, greq
    real    refpix, pixcwave, wavoffset
    real    wmn1, wmn2, wmx1, wmx2
    int     detector_upper_spec_limit[3]
    int     l_xbin, l_ybin, grule, speclen, width, center, grorder
    int     xccd[2], yccd[2], slitsep, ifuslitx[2], ifuslitnum
    int     nexten, l_dispaxis, nxpix, nypix, n_i, inst, nr, id1, iccd
    bool    l_verbose
    struct  sdate, l_wave_range

    # Query parameters
    l_inimages = inimages
    l_logfile = logfile
    l_gratingdb = gratingdb
    l_filterdb = filterdb
    l_key_dispaxis = key_dispaxis
    l_dispaxis = dispaxis
    l_verbose = verbose
    status = 0
    pi = 3.14159265
    
    grorder = 1

    # Pixel scale: scale[inst,iccd]
    scale[1,1] = 0.0727 # GMOS-N EEV CCDs
    scale[1,2] = 0.07288 # GMOS-N e2vDD CCDs
    scale[1,3] = 0.0807 # GMOS-N Hamamatsu CCDs 
    scale[2,1] = 0.073 # GMOS-S EEV CCDs
    scale[2,3] = 0.0800 # GMOS-S Hamamatsu CCDs ##M PIXEL_SCALE

    # Set the upper (red) spectral limit for each detector type
    detector_upper_spec_limit[1] = 1025
    detector_upper_spec_limit[2] = 1050 ##M
    detector_upper_spec_limit[3] = 1080 ##M

    # Keep imgets parameters from changing by outside world
    cache ("imgets", "fparse", "gemdate", "gimverify")

    # Test the logfile:
    if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GSAPPWAVE: Both gsappwave.logfile and "//\
                "gmos.logfile fields are empty", l_logfile, verbose+)
            printlog ("                    Using default file gmos.log", \
                l_logfile, verbose+)
        }
    }

    # Temporary files
    infiles = mktemp("tmpinfiles")
    mdf=mktemp("tmpmdf")

    # Start logging to file
    if (l_logfile == "STDOUT") {
        l_logfile = ""
        l_verbose = yes
    }

    date | scan(sdate)
    printlog ("-------------------------------------------------"//\
        "-------------------------------", l_logfile, verbose=l_verbose)
    printlog ("GSAPPWAVE -- "//sdate, l_logfile, verbose=l_verbose)
    printlog (" ", l_logfile, verbose=l_verbose)
    printlog ("inimages = "//l_inimages, l_logfile, verbose=l_verbose)

    #check that there are input files
    if (l_inimages == "" || l_inimages == " ") {
        printlog ("ERROR - GSAPPWAVE: input files not specified", \
            l_logfile, verbose+)
        goto error
    }

    # check existence of list file
    if (substr(l_inimages,1,1) == "@") {
        inlist = substr (l_inimages,2,strlen(l_inimages))
        if (!access(inlist)) {
            printlog ("ERROR - GSWAPPWAVE: Input list "//inlist//" not found",
                l_logfile, verbose+)
            goto error
        }
    }

    # Loop over inimages
    files(l_inimages,sort-, > infiles)
    scanfile=infiles
    while (fscan(scanfile, img) != EOF) {
        # Check existence of file
        gimverify(img)
        # if (!imaccess(img)) {
        if (gimverify.status!=0) {
            printlog ("ERROR - GSAPPWAVE: "//img//" does not exist or is"//
                "not a MEF file.", l_logfile, verbose+)
            goto error
        }
        img = gimverify.outname # .fits needed if looking at the MDF

        # get header information
        # Which instrument?
        imgets (img//"[0]", "INSTRUME", >& "dev$null")
        if (imgets.value == "0") {
            printlog ("ERROR - GSAPPWAVE: Instrument keyword not found ("//
                img//")", l_logfile, verbose+)
            goto error
        }
        inst = 1 # default is GMOS-N, support for old data
        if (imgets.value == "GMOS-S") {
            inst=2
        }

        # Read the detector type and set iccd accordingly
        keypar (img//"[0]", "DETTYPE", silent+) 
        if (keypar.found) {
            dettype = keypar.value
        } else {
            printlog ("ERROR - GSAPPWAVE: DETTYPE keyword not found.", \
                l_logfile, verbose+)
            goto error
        }
        if (dettype == "S10892" || dettype == "S10892-N") {# Hamamatsu CCDs
            iccd = 3
        } else if (dettype == "SDSU II e2v DD CCD42-90") {# e2vDD CCDs
            iccd = 2
        } else if (dettype == "SDSU II CCD") {# EEV CCDs
            iccd = 1
        } else {
            printlog ("ERROR - GSAPPWAVE: DETTYPE keyword "//dettype//\
                " not known", l_logfile, verbose+)
            goto error
        }

        imgets (img//"[0]","OBSMODE", >>& "dev$null")
        obsmode = imgets.value

        imgets (img//"[0]","GMOSAIC", >>& "dev$null")
        gmosaic = imgets.value

        imgets (img//"[0]","GSCUT", >>& "dev$null")
        gscut = imgets.value

        imgets (img//"[0]","GFEXTRAC", >>& "dev$null")
        gfextract = imgets.value

        # Throw an error if the data haven't been processed
        # sufficiently to get a meaningful result & without crashing
        if (gmosaic=="0") {
          printlog ("ERROR - GSAPPWAVE: run gmosaic before gsappwave",
              l_logfile, verbose+)
          goto error
        }

        imgets (img//"[0]","GRATING")
        l_grating = imgets.value

        imgets (img//"[sci,1]","CCDSUM")
        l_xbin = int(substr(imgets.value,1,1))
        l_ybin = int(substr(imgets.value,3,3))

        imgets (img//"[0]","FILTER1")
        l_filter[1] = imgets.value

        imgets (img//"[0]","FILTER2")
        l_filter[2] = imgets.value

        imgets (img//"[0]","GRWLEN")
        cwave = real(imgets.value)

        imgets (img//"[0]","GRTILT")
        tilt = real(imgets.value)*pi/180.
        
        imgets (img//"[0]", "GRORDER")
        grorder = int(imgets.value)

        imgets (img//"[0]","NSCIEXT")
        nexten = int(imgets.value)

        imgets (img//"[sci,1]","i_naxis1")
        nxpix = int(imgets.value)

        imgets (img//"[sci,1]","i_naxis2")
        nypix = int(imgets.value)

        # get grating information
        match (l_grating, l_gratingdb, stop-, print+, meta+) |\
            scan (gname, grule, gblaze, gR, gcoverage, gwave1, gwave2, \
            wavoffset)

        print (((cwave*grule)/1.e6)) | \
            interp ("gmos$data/gratingeq.dat", "STDIN", int_mode="spline",
            curve_gen-) | scan (greq, gtilt)

        gtilt = gtilt * pi/180.
        a = sin(gtilt+0.872665)/sin(gtilt)
        gR = 206265. * greq/(0.5*81.0*sin(gtilt)) # Assumes 0.5'' slit
        nmppx = (a*scale[inst,iccd]*real(l_xbin)*cwave*81.0*sin(gtilt)) / \
            (206265.*greq)
        slitsep = 175.*1.611444/(scale[inst,iccd]*a) # For the IFU

        wave1 = gwave1
        wave2 = gwave2
        
        if (grorder == 2) {
        	gR = gR*4.
        	nmppx = nmppx/2.
        }

        printlog ("Grating: "//gname, l_logfile, l_verbose)
        printlog ("Grating central wavelength: "//cwave//" [nm]", \
            l_logfile, l_verbose)
        printlog ("Resolution (0.5'' slit): "//gR, l_logfile, l_verbose)
        printlog ("Anamorphic factor: "//a, l_logfile, l_verbose)
        printlog ("Grating tilt (header): "//(tilt*180./pi), \
            l_logfile, l_verbose)
        printlog ("Calculated tilt: "//(gtilt*180./pi), l_logfile, l_verbose)
        printlog ("nm/pix = "//nmppx, l_logfile, l_verbose)

        # get filter information
        fwave1 = 0.0 ; wmn1 = 0.0 ; wmn2 = 0.0
        fwave2 = 9999.0 ; wmx1 = 99999.0 ; wmx2 = 99999.0
        if (l_filter[1] != "" && substr(l_filter[1],1,4) != "open" ) 
            match (l_filter[1], l_filterdb, stop-, print+, meta+) | \
                scan (fname, wmn1, wmx1, ffile)

        if (l_filter[2] != "" && substr(l_filter[2],1,4) != "open")
            match (l_filter[2], l_filterdb, stop-, print+, meta+) | \
                scan (fname, wmn2, wmx2, ffile)

        if (wmn1 > wmn2) {
            fwave1 = wmn1
        } else {
            fwave1 = wmn2
        }
        if (wmx1 < wmx2) {
            fwave2 = wmx1
        } else {
            fwave2 = wmx2
        }

        # determine whether filter or grating limits wavelength coverage
        if (fwave1 > wave1) {
            wave1 = fwave1
        }
        if (fwave2 < wave2) {
            wave2 = fwave2
        }
        # and red limit of CCD
        if (wave2 > real(detector_upper_spec_limit[iccd])) {
            wave2 = real(detector_upper_spec_limit[iccd])
        }

        # Det. length of spectrum in pixels and reference pixel (of cwave)
        if (obsmode == "IFU") {
          if (gfextract=="0") {
            # If gfextract hasn't been run yet, calculate the slit
            # position and wavelength range here, assuming the central
            # pixel of the data has not changed (nothing other than
            # gfextract should have cut out a sub-region).

            ifuslitx[1] = (3109-slitsep/2)/l_xbin+wavoffset/nmppx ##M 
            ifuslitx[2] = ifuslitx[1]+slitsep/l_xbin
            # (these are integers because that's how gfextract does it)

            speclen = nxpix

            # Figure out which IFU slit was used:
            ifuslitnum=0
            imgets (img//"[0]","MASKNAME")
            if (imgets.value == "IFU") {
              # I haven't found any reliable way to determine whether a
              # binary table exists in a file (apparently gfextract barfs
              # on this too), so we won't try to figure out the IFU slit
              # from the MDF here for older data -> comment this out.
              # tcopy(img//".fits[mdf]",mdf//".fits",verbose-, >>& "dev$null")
              # if (access(mdf//".fits")) {
              #  tinfo(mdf//".fits",ttout-)
              #  nr=tinfo.nrows
              #  tabpar(mdf//".fits","NO",1)
              #  id1=int(tabpar.value)
              #  delete (mdf//".fits")
              #  if ((id1==1 && nr<=750) || (id1==51 && nr<=350)) {
              #    ifuslitnum = 1
              #  } else if (id1==751) {
              #    ifuslitnum = 2
              #  }
              # }
            } else if (imgets.value == "IFU-R" || imgets.value == "IFU-NS-R") {
                ifuslitnum = 1
            } else if (imgets.value == "IFU-B" || imgets.value == "IFU-NS-B") {
                ifuslitnum = 2
            }
            if (ifuslitnum==0) {
              # For 2-slit mode, give a warning because we can't write a
              # single meaningful wavelength WCS here and then set the
              # central wavelength to the centre of the detector (ie.
              # the average of the slits)
              printlog ("WARNING - GSAPPWAVE: IFU in 2-slit (or "//
                  "undetermined) configuration;", l_logfile, l_verbose)
              printlog ("  using an average wavelength solution, correct "//
                  "to ~50% of each slit's", l_logfile, l_verbose)
              printlog ("  range; run gfextract first to get a proper "//
                  "solution", l_logfile, l_verbose)
              refpix = 0.5 * (ifuslitx[1] + ifuslitx[2])
            } else {
              refpix = ifuslitx[ifuslitnum]
            }
          } else {
            # If gfextract was run, get the range/refpix from its output
            imgets (img//"[SCI,1]", "REFPIX1")
            refpix = real(imgets.value)
            if (imgets.value == "0") {
                printlog ("ERROR - GSAPPWAVE: REFPIX1 keyword from "//
                  "gfextract not found", l_logfile, verbose+)
                goto error
            }
            imgets (img//"[SCI,1]","i_naxis1")
            speclen = int(imgets.value)
          }
          # Recalculate the wavelength limits
          wave2 = cwave+(refpix-1)*nmppx
          wave1 = wave2-(nxpix-1)*nmppx
          wave1 = max(wave1, fwave1)
          wave2 = min(wave2, fwave2, real(detector_upper_spec_limit[iccd]))
        # For non-IFU modes, calculate reference pixel and length
        # as previously, based on whether the spectra are cut out:
        } else if (gscut == "0") { # 0 means not found in the header
            speclen = nxpix
            refpix = real(nxpix)/2. + wavoffset/nmppx
        } else {
            speclen = nint((wave2-wave1)/nmppx)
            refpix = real(speclen) - (cwave-wavoffset-wave1)/nmppx
        }

        printlog ("Filter1: "//l_filter[1], l_logfile, l_verbose)
        printlog ("Filter2: "//l_filter[2], l_logfile, l_verbose)
        printlog ("Approximate available wavelength coverage "//\
            "based on grating,\n    filters and detector type is: "//\
            wave1//" - "//wave2//" [nm]\n", l_logfile, l_verbose)

        # Loop over extensions
        for (n_i=1; n_i<=nexten; n_i+=1) {
            sci_extn = "[SCI,"//n_i//"]"
            var_extn = "[VAR,"//n_i//"]"
            dq_extn  = "[DQ,"//n_i//"]"

            printlog ("Calibrating: "//sci_extn, l_logfile, l_verbose)
            if (imaccess(img//sci_extn)) {
                gemhedit (img//sci_extn, "WAT1_001", \
                    "wtype=linear label=Wavelength units=Angstroms", "", \ 
                    delete-)
                imgets (img//sci_extn, "REFPIX1", >& "dev$null")
                if (imgets.value!="0") {
                    refpix = real(imgets.value)
                }

                # Recalculate the sepctrum length in pixels and the 
                # wavelength coverage for each sci extension
                # Obtain the NAXIS1 value
                keypar (img//sci_extn, "i_naxis1", silent+)
                nxpix = int(keypar.value)

                # The spectrum length in pixel
                speclen = nxpix

                # Wavelength coverage (reassign wave1 and wave 2 using the 
                # following equations, remember red is on the left of the 
                # detector
                
                # Printlog statements
                printlog ("    Spectrum length in pixels: "//speclen, \
                    l_logfile, l_verbose)
                printlog ("    Approximate location of grating central \n"//\
                    "        wavelength in pixels: "//refpix, \
                    l_logfile, l_verbose)
                printlog ("    Approximate wavelength coverage based on", \
                    l_logfile, l_verbose)
                printf ("        spectrum length and calculated dispersion "//\
                    "is: %4.0f - %4.0f [nm]\n", wave1, wave2) | \
                        scan (l_wave_range)
                printlog (l_wave_range, l_logfile, l_verbose)

                # Delete LTV1/2 keywords as, LTV1 messes up gswavelength when 
                # SECX1 is not 1. LTV1 maps the offset of the cut position
                # to the origin (in physical units) of the image from which
                # it was cut. - MS 
                # Delete the physical offset (LTV) keywords so the WCS is now
                # in the frame of reference of the cut image and now not 
                # mapping to the parent image coordinate system. This leaves 
                # LTM keywords in place, as is. This is fine as they represent 
                # the scale/rotation matrix which is still needed. - MS
                gemhedit (img//sci_extn, "LTV1", "", "", delete+)
                gemhedit (img//sci_extn, "LTV2", "", "", delete+)

                # Update header keywords
                gemhedit (img//sci_extn, "CTYPE1", "LINEAR", "", delete-)
                gemhedit (img//sci_extn, "CRVAL1", (cwave*10.), "", delete-)
                gemhedit (img//sci_extn, "CRVAL2", 1., "", delete-)
                gemhedit (img//sci_extn, "CRPIX1", refpix, "", delete-)
                gemhedit (img//sci_extn, "CRPIX2", 1., "", delete-)
                gemhedit (img//sci_extn, "CD1_1", (-10.*nmppx), "", delete-)
                gemhedit (img//sci_extn, "CD2_2", 1.0, "", delete-)
                if (obsmode != "IFU") {
                    gemhedit (img//sci_extn, "CTYPE2", "LINEAR", "", delete-)
                    gemhedit (img//sci_extn, "WAT0_001", "system=image", "", 
                        delete-)
                    gemhedit (img//sci_extn, "WAT2_001", "wtype=linear", "", 
                        delete-)
                    gemhedit (img//sci_extn, "CD1_2", 0.0, "", delete+)
                    gemhedit (img//sci_extn, "CD2_1", 0.0, "", delete+)
                    gemhedit (img//sci_extn, "CD2_2", 1.0, "", delete-)
                }
                gemhedit (img//sci_extn, l_key_dispaxis, l_dispaxis,
                    "Dispersion axis", delete-)

            }
            if (imaccess(img//var_extn)) {
                gemhedit (img//var_extn, "WAT1_001",
                    "wtype=linear label=Wavelength units=Angstroms", "", 
                    delete-)
                gemhedit (img//var_extn, "CTYPE1", "LINEAR", "", 
                    delete-)
                gemhedit (img//var_extn, "CRVAL1", (cwave*10.), "", 
                    delete-)
                gemhedit (img//var_extn, "CRPIX1", refpix, "", 
                    delete-)
                gemhedit (img//var_extn, "CD1_1", (-10.*nmppx), "", 
                    delete-)
                if (obsmode != "IFU") {
                    gemhedit (img//var_extn, "WAT0_001", 
                        "system=image", "", delete-)
                    gemhedit (img//var_extn, "WAT2_001", 
                        "wtype=linear", "", delete-)
                    gemhedit (img//var_extn, "CD2_1", 0.0, "", 
                        delete-)
                    gemhedit (img//var_extn, "CD2_2", 1.0, "", 
                        delete-)
                }
                gemhedit (img//var_extn, "CD2_2", 1.0, "", delete-)
                gemhedit (img//var_extn, l_key_dispaxis, l_dispaxis,
                    "Dispersion axis", delete-)
            }
            if (imaccess(img//dq_extn)) {
                gemhedit (img//dq_extn, "WAT1_001",
                    "wtype=linear label=Wavelength units=Angstroms", "", 
                    delete-)
                gemhedit (img//dq_extn, "CTYPE1", "LINEAR", "", delete-)
                gemhedit (img//dq_extn, "CRVAL1", (cwave*10.), "", delete-)
                gemhedit (img//dq_extn, "CRPIX1", refpix, "", delete-)
                gemhedit (img//dq_extn, "CD1_1", (-10.*nmppx), "", delete-)
                if (obsmode != "IFU") {
                    gemhedit (img//dq_extn, "WAT0_001", "system=image",
                        "", delete-)
                    gemhedit (img//dq_extn, "WAT2_001", "wtype=linear",
                        "", delete-)
                    gemhedit (img//dq_extn, "CD2_1", 0.0, "", delete-)
                    gemhedit (img//dq_extn, "CD2_2", 1.0, "", delete-)
                }
                gemhedit (img//dq_extn, "CD2_2", 1.0, "",
                    delete-)
                gemhedit (img//dq_extn, l_key_dispaxis, l_dispaxis,
                    "Dispersion axis", delete-)
            }
        }

        # final header update
        gemdate ()
        gemhedit (img//"[0]", "GSAPPWAV", gemdate.outdate,
            "UT Time stamp for GSAPPWAVE", delete-)
        gemhedit (img//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
    }
    scanfile = ""
    # clean up
    goto clean

error:
    status = 1
    goto clean

clean:
    scanfile = ""
    delete (infiles, verify-, >& "dev$null")
    # close log file
    printlog (" ", l_logfile, l_verbose)
    printlog ("GSAPPWAVE done", l_logfile, l_verbose)
    printlog ("------------------------------------------------------------\
        --------------------", l_logfile, l_verbose )

end
