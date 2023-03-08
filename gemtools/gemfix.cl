# Copyright(c) 2013-2015 Association of Universities for Research in Astronomy, Inc.
#
# Version Apr 26, 2013  James E.H. Turner  (for my IFU data).
#         May 12, 2014  JT  Stop fixpix shifting DQ erroneously WRT SCI.
#

procedure gemfix (inimages, outimages)

# Take a list of input MEF files with DQ extensions and run fixpix on all
# the SCI extensions wherever DQ is not good, creating new output files.

char    inimages     {prompt="Input images"}
char    outimages    {prompt="Output images"}
#char    linterp      {"INDEF", prompt="Mask values for line interpolation"}
#char    cinterp      {"INDEF", prompt="Mask values for column interpolation"}
char    method       {"fixpix", enum="fit1d|fixpix",
                      prompt="Interpolation method"}
real    grow         {1.5,min=0,prompt="Substitution growth radius in pix"}
int     bitmask      {65535, prompt="Mask for DQ bits to ignore=0 or use=1"}
int     axis         {1,min=1,max=2,prompt="Fitting axis for method=fit1d"}
int     order        {0,min=0,prompt="Order for fit1d (0=default)"}
real    low_reject   {3., prompt="Low sigma rejection threshold (fit1d)"}
real    high_reject  {2.3, prompt="High sigma rejection threshold (fit1d)"}
int     niterate     {5, prompt="Number of rejection iterations (fit1d)"}
bool    fl_inter     {no, prompt="Fit interactively (for method=fit1d)?"}
char    logfile      {"", prompt="Logfile name"}
bool    verbose      {yes, prompt="Verbose?"}
int     status       {0, prompt="Exit status (0=good)"}
struct  *flist       {"", prompt="Internal use only"}

begin

    # Define local variables
    char    l_inimages, l_outimages, l_method, l_logfile, l_sci_ext, l_dq_ext
    bool    l_verbose, l_fl_inter
    real    l_grow
    int     l_bitmask, l_axis, l_niterate, l_order
    real    l_low_reject, l_high_reject

    int     i, j, nbad, nimg, nsciext
    char    infiles, infile[500], outfiles, outfile[500], img, tmpclean
    char    sciext, sciout, dqext, l_datasec, scireg, dqreg, tmpsci, tmpdq
    struct  sdate
    bool    trimmed

    # Initialize local variables
    l_inimages = inimages
    l_outimages = outimages
    l_method = method
    l_grow = grow
    l_bitmask = bitmask
    l_axis = axis
    l_order = order
    l_low_reject = low_reject
    l_high_reject = high_reject
    l_niterate = niterate
    l_fl_inter = fl_inter
    l_logfile = logfile
    l_verbose = verbose
    l_sci_ext = "SCI"
    l_dq_ext = "DQ"

    # Initialize exit status
    status = 0

    # Cache some parameter files
    cache("keypar", "gemhedit", "gimverify")

    # Make temporary files
    infiles = mktemp("tmpinfiles")
    outfiles = mktemp("tmpoutfiles")

    # Test the logfile:
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
        l_logfile = gemtools.logfile
        if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
            l_logfile = "gemtools.log"
            printlog ("WARNING - GEMFIX: Both gemfix.logfile and \
                gemtools.logfile fields are empty", l_logfile, verbose+)
            printlog ("                    Using default file gemtools.log",
                l_logfile, verbose+)
        }
    }

    # Start logging
    date | scan(sdate)
    printlog ("---------------------------------------------------------------\
        -------------", l_logfile, l_verbose)
    printlog ("GEMFIX -- "//sdate, l_logfile, l_verbose)

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
            printlog("ERROR - GEMFIX: Image "//img//" doesn't exist"+\
                " or not MEF", l_logfile, verbose+)
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
                    printlog ("ERROR - GEMFIX: Cannot find the number of \
                        science extensions", l_logfile, verbose+)
                    goto crash
                }
            }
        }
    }
    flist = ""
    nimg = i
    if (nbad > 0) goto crash

    # Check the output images
    i = 0
    flist = outfiles
    while (fscan(flist, img) != EOF) {
        gimverify(img)
        if (gimverify.status != 1) {
            printlog("ERROR - GEMFIX: Output image "//img//" exists", \
                l_logfile, verbose+)
            nbad = nbad + 1
        } else {
            i = i + 1
            outfile[i] = gimverify.outname//".fits"
	}
    } # end loop over output files
    flist = ""
    if (nbad > 0) goto crash

    # Make sure we have the same number of input & output images:
    if (nimg == 0 || i != nimg) {
        printlog("ERROR - GSCRSPEC: Different number of (or zero) input & "+
          "output images", l_logfile, verbose+)
        goto crash
    }

    # Loop over the number of input images
    for (j = 1; j <= nimg ; j += 1) {

        tmpclean = mktemp("tmpclean")

        # Copy the input image before replacing values:
        copy (infile[j], tmpclean//".fits", verbose-, >& "dev$null")

        # Loop over the number of extensions in each input image
        for (i = 1; i <= nsciext; i += 1) {

            sciext = tmpclean//"["//l_sci_ext//","//i//"]"
            sciout = tmpclean//"["//l_sci_ext//","//i//",overwrite]"
            dqext = tmpclean//"["//l_dq_ext//","//i//"]"

            # Should check here that the SCI & DQ dimensions match, otherwise
            # the task will crash at one step or another.

            # Don't do this for now -- it shouldn't be necessary and becomes
            # overcomplicated when we end up with sections of sections below.
            #
            # # Only fix the non-overscan region, if applicable (though really
            # # we expect to be operating on trimmed data). DQ should match the
            # # science data, ie. be padded with overscan if necessary. The
            # # fixpix step below can operate on an image subsection in place,
            # # so no need to copy bits of images from different places.
            # keypar (sciext, "TRIMSEC", silent+)
            # if (keypar.found) {
            #     l_datasec = ""
            # } else {
            #     keypar (sciext, "DATASEC", silent+)
            #     if (keypar.found) {
            #         l_datasec = keypar.value
            #     } else {
	    #         l_datasec = ""
            #     }
	    # }
            # scireg = sciext//l_datasec   # replace sciext with this below
            # dqreg = dqext//l_datasec

            # Run fixpix on the SCI extension as long as there's a matching DQ
	    if (imaccess(dqext)) {

	        printlog ("  "//infile[j]//"["//l_sci_ext//","//i//"]", \
                    l_logfile, verbose+)

                # Running fixpix using the original DQ isn't all that
                # effective because the very edges of CRs aren't getting
                # masked out at the LACosmic step, so interpolating from
                # those pixels produces higher counts than it should (the CR
                # flux is still greatly reduced but the area of contamination
                # is not). However, a fixed 5x5 median like LACosmic uses is
                # also liable to include contamination from large CRs. What
                # we do here is grow the DQ mask by convolution just for the
                # purpose of replacing contaminated values, without actually
                # masking out more pixels. This seems to work pretty well.
                # The replacement value doesn't have to be spot on, just
                # reliably close.

                tmpdq = mktemp("tmpdq")
                tmpsci = mktemp("tmpsci")

                # Make a 3rd copy of the DQ with the bitmask applied and
                # where we can grow the DQ regions as needed:
                imexpr("a & b", tmpdq, dqext, l_bitmask, dims="auto",
                  intype="auto", outtype="auto", refim="a", bwidth=0,
                  rangecheck+, verbose-, exprdb="none")

                # Commented out this convolve step because we have to run
                # imreplace after convolving anyway (to deal with flux
                # non-conservation or fractional values) and imreplace
                # already has its own built-in radius option.
                # convolve (dqext, tmpdq, kernel="", xkernel="1 1 1",
                #  ykernel="1 1 1", bilinear+, radsym-, boundary="nearest")

                # Grow the DQ, converting all bad values to 1:
                imreplace (tmpdq, 1.0, imaginary=0.0, lower=0.5, upper=INDEF,
                    radius=l_grow)

                # Make sure DQ has the same WCS as SCI, since fixpix, in its
                # dubious wisdom, applies relative shifts otherwise:
                wcscopy(tmpdq, sciext, verbose-)

                if (l_method=="fixpix") {

                    # First we'll fix a copy of the science image, including
                    # growth pixels, then later we'll just use the values
                    # from it corresponding to the original DQ before growth:
                    imcopy(sciext, tmpsci, verbose-)

                    # proto.fixpix interpolates in the narrowest dimension by
                    # default as long as linterp & cinterp are INDEF:
                    proto.fixpix (tmpsci, tmpdq, linterp="INDEF", \
                      cinterp="INDEF", verbose-, pixels-)

	        } else if (l_method=="fit1d") {

                    # Use a low-to-intermediate-order fit (16 for GMOS 1x1):
		    if (l_order==0) {
                        keypar(sciext, "i_naxis"//l_axis, silent+)
                        l_order = max(2, nint (real(keypar.value) / 128.))
                        printlog ("    order="//l_order, l_logfile, l_verbose)
		    }

                    fit1d(sciext, tmpsci, type="fit", bpm=tmpdq, axis=l_axis,
                      interactive=l_fl_inter, sample="*", naverage=1,
                      function="chebyshev", order=l_order,
                      low_reject=l_low_reject, high_reject=l_high_reject,
                      niterate=l_niterate, grow=0., graphics="stdgraph",
                      cursor="")

		} else {
                    # IRAF should make this impossible, but still...
                    printlog ("ERROR - GEMFIX: invalid method", l_logfile, \
                      verbose+)
                    imdelete(tmpdq, verify-, >& "dev$null")
                    goto crash
		} # end (select which method)

                # Where the input DQ matches the bitmask, replace the science
                # pixel values with the ones from the interpolated image:
                imexpr("(a & b) == 0 ? c : d", sciout, dqext, l_bitmask,
                  sciext, tmpsci, dims="auto", intype="auto", outtype="auto",
                  refim="c", bwidth=0, rangecheck+, verbose-, exprdb="none")

                imdelete(tmpdq, verify-, >& "dev$null")
                imdelete(tmpsci, verify-, >& "dev$null")

            } else {
                printlog ("ERROR - GEMFIX: "//dqext//" not found", \
                    l_logfile, verbose+)
                imdelete(tmpclean, verify-, >& "dev$null")
                goto crash
	    } # end (if there's a DQ ext)

        } # end i <= nsciext

        # Rename the working image to the final output name:
        rename(tmpclean//".fits", outfile[j], field="all", >& "dev$null")

    } # end loop over input images

    goto clean


crash:
    # Exit with error subroutine
    status = 1

clean:
    # Clean up
    delete(infiles, verify-, >& "dev$null")
    delete(outfiles, verify-, >& "dev$null")

    # Close log file
    if (status == 0) {
        printlog ("", l_logfile, l_verbose)
        printlog ("GEMFIX exit status: Good", l_logfile, l_verbose)
        printlog ("-----------------------------------------------------------\
            -----------------", l_logfile, l_verbose)
    } else {
        printlog ("", l_logfile, l_verbose)
        printlog ("GEMFIX exit status: Error", l_logfile, l_verbose)
        printlog ("-----------------------------------------------------------\
            -----------------", l_logfile, l_verbose)
    }

end
