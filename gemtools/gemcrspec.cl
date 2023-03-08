# Copyright(c) 2015-2017 Association of Universities for Research in Astronomy, Inc.

procedure gemcrspec (inimages, outimages)

# gemcrspec is an MEF wrapper for the spectroscopic version of the Laplacian
# Cosmic Ray Identification routine by Pieter G. van Dokkum (see
# http://www.astro.yale.edu/dokkum/lacosmic/) lacos_spec
# (http://www.astro.yale.edu/dokkum/lacosmic/lacos_spec.cl) 
#
# gemcrspec removes cosmic rays from spectroscopic data
#
# The lacos_spec and gemcrspec tasks should be defined as follows, before
# running gemcrspec
# 
#   ecl> task lacos_spec=/your/path/to/lacos_spec.cl
#   ecl> task gemcrspec=/your/path/to/gemcrspec.cl 
#
# Version Sep 26, 2002  Bryan Miller  Original version (gscrspec)
#         Nov  4, 2011  Emma Hogan    Standardized script
#         Apr 26, 2013  James Turner  Fix loop over files; use correct DQ=8
#         Sep 24, 2014  James Turner  Try generalizing a bit from GMOS; rename;
#                                     deal with missing input VAR/DQ
#         Feb 20, 2017  James Turner  Don't delete existing output files

char    inimages      {prompt="Input images"}
char    outimages     {prompt="Output images"}
int     xorder       {9, prompt="Order of object fit (0=no fit)"}
int     yorder       {-1, prompt="Order of sky line fit (0=no fit, -1=autosize)"}
real    sigclip      {4.5, prompt="Detection limit for cosmic rays (sigma)"}
real    sigfrac      {0.5, prompt="Fractional detection limit for neighbouring pix"}
real    objlim       {1.0, prompt="Contrast limit between CR and underlying object"}
int     niter        {4, prompt="Maximum number of iterations"}
bool    fl_vardq     {no, prompt="Propagate variance and data quality frames?"}
char    sci_ext      {"SCI", prompt="Name of science extension"}
char    var_ext      {"VAR", prompt="Name of variance extension"}
char    dq_ext       {"DQ", prompt="Name of data quality extension"}
char    key_ron      {"RDNOISE", prompt="Header keyword for readout noise"}
char    key_gain     {"GAIN", prompt="Header keyword for gain"}
real    ron          {3.5, prompt="Readout noise value in electrons (if not in header)"}
real    gain         {2.2, prompt="Gain value in e-/ADU (if not in header)"}
char    logfile      {"", prompt="Logfile name"}
bool    verbose      {yes, prompt="Verbose?"}
int     status       {0, prompt="Exit status (0=good)"}
struct  *flist       {"", prompt="Internal use only"}

begin

    # Define local variables
    char    l_inimages, l_outimages, l_sci_ext, l_var_ext, l_dq_ext, mdf
    char    tmpclean, tmpmask, l_key_gain, l_key_ron, l_logfile, masksum
    real    l_gain, l_ron, l_sigclip, l_sigfrac, l_objlim
    bool    l_fl_vardq, l_verbose
    int     l_xorder, l_yorder, l_niter
    char    infiles, infile[500], outfiles, outfile[500], img, suf
    int     nsciext, i, j, nbad, nimg, ny, dispaxis
    struct  sdate

    # Initialize exit status
    status = 0

    # Make temporary files
    mdf = mktemp("tmpmdf")
    infiles = mktemp("tmpinfiles")
    outfiles = mktemp("tmpoutfiles")

    # Initialize local variables
    l_inimages = inimages
    l_outimages = outimages
    l_xorder = xorder
    l_yorder = yorder
    l_sigclip = sigclip
    l_sigfrac = sigfrac 
    l_objlim = objlim
    l_niter = niter
    l_fl_vardq = fl_vardq
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_key_ron = key_ron
    l_key_gain = key_gain
    l_ron = ron
    l_gain = gain
    l_logfile = logfile
    l_verbose = verbose

    # Cache some parameter files
    cache("keypar", "gemhedit", "gimverify")

    # Test the logfile:
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
        l_logfile = gemtools.logfile
        if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
            l_logfile = "gemtools.log"
            printlog ("WARNING - GEMCRSPEC: Both gemcrspec.logfile and \
                gemtools.logfile fields are empty", l_logfile, verbose+)
            printlog ("                    Using default file gemtools.log",
                l_logfile, verbose+)
        }
    }

    # Start logging
    date | scan(sdate)
    printlog ("---------------------------------------------------------------\
        -------------", l_logfile, l_verbose)
    printlog ("GEMCRSPEC -- "//sdate, l_logfile, l_verbose)

    # Check that we actually have the LACosmic dependency available:
    if (!deftask("lacos_spec")) {
        printlog ("ERROR - GEMCRSPEC: lacos_spec.cl dependency not found!",
            l_logfile, verbose+)
        printlog ("                   See \
            http://www.astro.yale.edu/dokkum/lacosmic", l_logfile, verbose+)
        goto crash
    }

    # Get input files
    sections(l_inimages, option="fullname", > infiles)
    sections(l_outimages, option="fullname", > outfiles)

    # Check the input images
    nbad = 0
    i = 0
    flist = infiles
    while (fscan(flist, img) != EOF) {
        gimverify(img)
        if (gimverify.status > 0) {
            printlog("ERROR - GEMCRSPEC: Image "//img//" doesn't exist \
                or not MEF", l_logfile, verbose+)
            nbad = nbad + 1
        } else {
            i = i + 1
            # Name w/o suffix
            infile[i] = gimverify.outname//".fits"
            if (i == 1) {
                keypar (infile[i]//"[0]", "NSCIEXT", silent+)
                if (keypar.found) {
                    nsciext = int(keypar.value)
                } else {
                    printlog ("ERROR - GEMCRSPEC: Cannot find the number of \
                        science extensions", l_logfile, verbose+)
                    nbad = nbad + 1
                }
                keypar (infile[i]//"["//l_sci_ext//"]", "DISPAXIS", silent+)
                if (keypar.found) {
                    dispaxis = int(keypar.value)
                } else {
                    printlog ("ERROR - GEMCRSPEC: Missing DISPAXIS header \
                        keyword in "//l_sci_ext//" ext", l_logfile, verbose+)
                    nbad = nbad + 1
                }
                if (dispaxis != 1) {
                    printlog ("ERROR - GEMCRSPEC: Only dispersion along x \
                        axis is currently supported", l_logfile, verbose+)
                    nbad = nbad + 1
	        }
            }
        }
    }
    flist = ""
    nimg = i
    if (nbad > 0) {
        goto crash
    }

    # Check the output images
    i = 0
    flist = outfiles
    while (fscan(flist, img) != EOF) {
        gimverify(img)
        if (gimverify.status != 1) {
            printlog("ERROR - GEMCRSPEC: Output image "//img//" exists", \
                l_logfile, verbose+)
            nbad = nbad + 1
        } else {
            i = i + 1
            outfile[i] = gimverify.outname//".fits"
        }
    } # end loop over output files
    flist = ""
    if (nbad > 0) {
        goto crash
    }

    # Make sure we have the same number of input & output images:
    if (nimg == 0 || i != nimg) {
        printlog("ERROR - GEMCRSPEC: Different number of (or zero) input & \
          output images", l_logfile, verbose+)
        goto crash
    }

    # Loop over the number of input images
    for (j = 1; j <= nimg ; j += 1) {

        # Copy the MDF to the output image
        tcopy (infile[j]//"[MDF]", mdf//".fits", verbose-)
        wmef (mdf//".fits", outfile[j], extnames="MDF", verbose-, \
            phu=infile[j], >& "dev$null")
        delete (mdf//".fits", verify-, >& "dev$null")
        if (wmef.status != 0) {
            printlog("ERROR - GEMCRSPEC: failed to create initial "//outfile[j],
                     l_logfile, verbose+)
            imdelete (outfile[j], verify-, >& "dev$null")
            goto crash
        }

        # Loop over the number of extensions in each input image
        for (i = 1; i <= nsciext; i += 1) {
            tmpclean = mktemp("tmpclean")
            tmpmask = mktemp("tmpmask")

            # Use the gain as defined in the header, if it exists
            keypar (infile[j]//"["//l_sci_ext//","//i//"]", l_key_gain, 
                silent+)
            if (keypar.found) {
                l_gain = real(keypar.value)
                printlog ("Using the gain value from the header = "//l_gain, 
                    l_logfile, l_verbose)
            } else {
                printlog ("Using the default gain value = "//l_gain, 
                    l_logfile, l_verbose)
            }

            # Use the readout noise as defined in the header, if it exists
            keypar (infile[j]//"["//l_sci_ext//","//i//"]", l_key_ron, silent+)
            if (keypar.found) {
                l_ron = real(keypar.value)
                printlog ("Using the readout noise value from the header = "//
                    l_ron, l_logfile, l_verbose)
            } else {
                printlog ("Using the default readout noise value = "//l_ron, 
                    l_logfile, l_verbose)
            }

            # Get height and set order appropriately
            keypar (infile[j]//"["//l_sci_ext//","//i//"]",
                "i_naxis"//dispaxis, silent+)
            if (keypar.found) {
                ny = int(keypar.value)
            }
            if (yorder < 0) {
                if (ny < 50) {
                    l_yorder = 2
                } else if (ny < 80) {
                    l_yorder = 3
                } else {
                    l_yorder = 5
                }
            }
            printlog ("Extension: " // i // " ny = " // ny // " order = \
                " // l_yorder, l_logfile, l_verbose)

            # Run lacos_spec
            lacos_spec (infile[j]//"["//l_sci_ext//","//i//"]", tmpclean, 
                tmpmask, gain=l_gain, readn=l_ron, xorder=l_xorder,
                yorder=l_yorder, sigclip=l_sigclip, sigfrac=l_sigfrac,
                objlim=l_objlim, niter=l_niter, verbose=l_verbose)

            # Combine masks
            if (l_fl_vardq) {
                chpixtype (tmpmask, tmpmask, "ushort", verbose-)
                imreplace (tmpmask, 8, imaginary=0, lower=0.5, upper=INDEF,
                    radius=0.)  # DQ=8 is used for cosmic rays
                flpr  # avoid "pixel storage file truncated" from addmasks
                if (imaccess(infile[j]//"["//l_dq_ext//","//i//"]")) {
                    masksum = mktemp("tmpmasksum")
                    addmasks(infile[j]//"[DQ,"//i//"],"//tmpmask//".fits", \
                        masksum//".fits","im1 || im2", flags=" ")
                    delete(tmpmask//".fits", verify-, >& "dev$null")
                    tmpmask = masksum
                } else {
                    printlog ("WARNING - GEMCRSPEC: No input "//l_dq_ext// \
                        " to propagate; creating from scratch\n", l_logfile,
                        verbose+)
                }
            }

            # Insert the science extensions into the output image
            imcopy (tmpclean//".fits", outfile[j]//"["//l_sci_ext//","//i//
                ",append]", verbose-, >& "dev$null")

            # Insert the var and dq extensions into the output image, if
            # appropriate 
            if (l_fl_vardq) {
	        if (imaccess(infile[j]//"["//l_var_ext//","//i//"]")) {
                    imcopy (infile[j]//"["//l_var_ext//","//i//"]", 
                        outfile[j]//"["//l_var_ext//","//i//",append]",
                        verbose-, >& "dev$null")
                } else {
                    printlog ("WARNING - GEMCRSPEC: No input "//l_var_ext// \
                        "; skipping propagation\n", l_logfile,
                        verbose+)
		}
                imcopy (tmpmask//".fits", outfile[j]//"["//l_dq_ext//","//i//
                    ",append]", verbose-, >& "dev$null")
            }

            # Clean up
            imdelete(tmpclean//","//tmpmask, verify-, >& "dev$null")

        } # end i <= nsciext

        # Update headers
        gemdate ()
        gemhedit (outfile[j]//"[0]", "GEM-TLM", gemdate.outdate, 
            "UT Last modification with GEMINI", delete-)
        gemhedit (outfile[j]//"[0]", "GEMCRSPEC", gemdate.outdate, 
            "UT Time stamp for GEMCRSPEC", delete-)

    } # end loop over input images
    goto clean

crash:
    # Exit with error subroutine
    status = 1
    goto clean

clean:
    # Clean up
    delete(infiles, verify-, >& "dev$null")
    delete(outfiles, verify-, >& "dev$null")

    # Close log file
    if (status == 0) {
        printlog ("", l_logfile, l_verbose)
        printlog ("GEMCRSPEC exit status: Good", l_logfile, l_verbose)
        printlog ("-----------------------------------------------------------\
            -----------------", l_logfile, l_verbose)
    } else {
        printlog ("", l_logfile, l_verbose)
        printlog ("GEMCRSPEC exit status: Error", l_logfile, l_verbose)
        printlog ("-----------------------------------------------------------\
            -----------------", l_logfile, l_verbose)
    }

end
