# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
#
# Version May 16, 2014  James E.H. Turner  (for my IFU data).

procedure growdq (inimages, outimages)

# Take a list of input MEF files with DQ extensions and run fixpix on all
# the SCI extensions wherever DQ is not good, creating new output files.

char    inimages     {prompt="Input images"}
char    outimages    {prompt="Output images"}
real    radius       {1.5, min=0, prompt="Adjacent growth radius in pixels"}
int     bitmask      {65535, prompt="Mask for DQ bits to ignore=0 or use=1"}
char    logfile      {"", prompt="Logfile name"}
bool    verbose      {yes, prompt="Verbose?"}
int     status       {0, prompt="Exit status (0=good)"}
struct  *flist       {"", prompt="Internal use only"}

begin

    # Define local variables
    char    l_inimages, l_outimages, l_logfile, l_sci_ext, l_dq_ext
    bool    l_verbose
    real    l_radius
    int     l_bitmask

    int     i, j, n, nbad, nimg, nsciext, dqmin, dqmax, ndq
    char    infiles, infile[500], outfiles, outfile[500], img, tmpout
    char    dqout, dqext, dqtype, tmpdq1, tmpdq2, tmpaccum, kernel
    struct  sdate

    # Initialize local variables
    l_inimages = inimages
    l_outimages = outimages
    l_radius = radius
    l_bitmask = bitmask
    l_logfile = logfile
    l_verbose = verbose
    l_sci_ext = "SCI"
    l_dq_ext = "DQ"

    # Initialize exit status
    status = 0

    # Cache some parameter files
    cache("keypar", "gimverify")

    # Make temporary files
    infiles = mktemp("tmpinfiles")
    outfiles = mktemp("tmpoutfiles")

    # Test the logfile:
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
        l_logfile = gemtools.logfile
        if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
            l_logfile = "gemtools.log"
            printlog ("WARNING - GROWDQ: Both growdq.logfile and \
                gemtools.logfile fields are empty", l_logfile, verbose+)
            printlog ("                    Using default file gemtools.log",
                l_logfile, verbose+)
        }
    }

    # Start logging
    date | scan(sdate)
    printlog ("---------------------------------------------------------------\
        -------------", l_logfile, l_verbose)
    printlog ("GROWDQ -- "//sdate, l_logfile, l_verbose)

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
            printlog("ERROR - GROWDQ: Image "//img//" doesn't exist"+\
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
                    printlog ("ERROR - GROWDQ: Cannot find the number of \
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
            printlog("ERROR - GROWDQ: Output image "//img//" exists", \
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

    # Temporary filenames:
    tmpdq1 = mktemp("tmpdq1")
    tmpdq2 = mktemp("tmpdq2")
    tmpaccum = mktemp("tmpaccum")

    # Loop over the number of input images
    for (j = 1; j <= nimg ; j += 1) {

        tmpout = mktemp("tmpout")

        # Copy the input image before replacing values:
        copy (infile[j], tmpout//".fits", verbose-, >& "dev$null")

        # Loop over the number of extensions in each input image
        for (i = 1; i <= nsciext; i += 1) {

            dqext = tmpout//"["//l_dq_ext//","//i//"]"
            dqout = tmpout//"["//l_dq_ext//","//i//",overwrite]"

            # Process the extension if it actually exists, otherwise complain:
            if (imaccess(dqext)) {

                printlog ("  "//infile[j]//"["//l_dq_ext//","//i//"]", \
                    l_logfile, verbose+)

                # Process each applicable bit plane of the DQ array
                # (do DQ=0 only if everything is blank, so that the output
                # still gets created):
                imstat(dqext, fields="max", format=no) | scan (dqmax)
                if (dqmax == 0) dqmin = 0
                else dqmin = 1

                # Decide on output data type:
                if (dqmax < 65536) dqtype = "ushort"
                else dqtype = "uint"

                # Start with the existing input DQ:
                imcopy(dqext, tmpaccum//".fits", verbose-)

                # Iterate over the range of bits actually used (from 1):
                for (n=dqmin; n <= dqmax; n*=2) {

                    # Extract bit plane n, if included in the bitmask:
                    imexpr("a & b & c", tmpdq1//".fits", dqext, n,
                        l_bitmask, dims="auto", intype="auto",
                        outtype=dqtype, refim="a", bwidth=0, rangecheck+,
                        verbose-, exprdb="none")

                    # Process plane if it has non-zero values:
                    imstat(tmpdq1//".fits", fields="npix", lower=1,
                        upper=INDEF, format=no) | scan (ndq)

                    if (ndq > 0) {

                        # Grow the DQ by the specified radius, replacing
                        # affected pixels with the appropriate bit value:
                        imreplace (tmpdq1, n, imaginary=0.0, lower=0.5,
                            upper=INDEF, radius=l_radius)

                        # Add transformed plane to accumulated DQ:
                        imexpr("a | b", tmpdq2//".fits", tmpaccum//".fits",
                            tmpdq1//".fits", outtype=dqtype, verbose-)
                        imdelete (tmpdq1//".fits", verify-, >& "dev$null")
                        imdelete (tmpaccum//".fits", verify-, >& "dev$null")
                        imrename (tmpdq2//".fits", tmpaccum//".fits", verbose-)

		    } else {  # (no pixels to process in this plane)

                        imdelete(tmpdq1, verify-, >& "dev$null")

                    } # end (process bitplane with non-zero values)

                    if (n==0) n=1  # otherwise *2 loop gets stuck at 0

                } # End n <= dqmax

                # Update the output DQ if one or more planes was non-zero
                # after applying the bitmask:
                imcopy(tmpaccum//".fits", dqout, verbose-)
                imdelete(tmpaccum//".fits", verify-, >& "dev$null")

            } else {
                printlog ("ERROR - GROWDQ: "//dqext//" not found", \
                    l_logfile, verbose+)
                imdelete(tmpout, verify-, >& "dev$null")
                goto crash
	    } # end (if there's a DQ ext)

        } # end i <= nsciext

        # Rename the working image to the final output name:
        rename(tmpout//".fits", outfile[j], field="all", >& "dev$null")

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
        printlog ("GROWDQ exit status: Good", l_logfile, l_verbose)
        printlog ("-----------------------------------------------------------\
            -----------------", l_logfile, l_verbose)
    } else {
        printlog ("", l_logfile, l_verbose)
        printlog ("GROWDQ exit status: Error", l_logfile, l_verbose)
        printlog ("-----------------------------------------------------------\
            -----------------", l_logfile, l_verbose)
    }

end
