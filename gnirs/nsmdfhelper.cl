# Copyright(c) 2001-2013 Association of Universities for Research in Astronomy, Inc.

procedure nsmdfhelper (mdf, row, image)

# Calculate various regions from an MDF file.  Two different MDF
# formats are used here.  The first is based on information in
# http://www.gemini.edu/sciops/instruments/gmos/gmosMOSobserver.html
# and
# http://www.gemini.edu/sciops/instruments/gmos/gmosmaskmakingv104.ps.gz
# PLUS the columns specorder and corner.  The second is historical -
# the MDF used in the first gnirs package release - is based on some
# confused conversations and reverse engineering and is defined only
# by the code below.

char    mdf             {prompt = "MDF file (a fits table)"}
int     row             {prompt = "MDF row to use"}
char    image           {prompt = "File used for header and size information"}

char    area            {"spectrum", enum="spectrum|image|ifu-map", prompt="Calculation to perform"}
char    mapfile         {"", prompt="Output file for IFU map"}
real    pixscale        {INDEF, prompt="Default pixel scale (if not in header)"}
int     dispaxis        {INDEF, prompt="Default dispersion axis (if not in header)"}

int     ixlo            {INDEF, prompt = "Lower left x coord (output)"}
int     iylo            {INDEF, prompt = "Lower left y coord (output)"}
int     ixhi            {INDEF, prompt = "Upper right x coord (output)"}
int     iyhi            {INDEF, prompt = "Upper right y coord (output)"}
int     ixovershoot     {INDEF, prompt = "Overshoot in x when very slanted (output)"}

real    xlo             {INDEF, prompt = "Lower left x coord (output)"}
real    ylo             {INDEF, prompt = "Lower left y coord (output)"}
real    xhi             {INDEF, prompt = "Upper right x coord (output)"}
real    yhi             {INDEF, prompt = "Upper right y coord (output)"}

real    corner          {INDEF, prompt = "Corner region for XD data (output)"}
real    specorder       {INDEF, prompt = "Spectral order from MDF (output)"}
real    slitwidth       {INDEF, prompt = "Slit width from MDF (output)"}

char    logfile         {"", prompt = "Logfile"}
char    logname         {"NSMDFHELPER", prompt = "Name printed in messages"}
bool    verbose         {yes, prompt = "Verbose output?"}
int     status          {0, prompt = "Exit status (0=good)"}

begin
    char    l_mdf = ""
    int     l_row
    char    l_image = ""
    char    l_area = ""
    char    l_mapfile = ""
    real    l_pixscale
    int     l_dispaxis
    char    l_logfile = ""
    char    l_logname = ""
    bool    l_verbose

    char    l_key_pixscale = ""
    char    l_key_dispaxis = ""

    int     junk, axis
    bool    debug, historical, haveheader, warnedya, even
    real    x_ccd, y_ccd, slitpos_mx, slitpos_my, slitsize_mx
    real    slitsize_my, slittilt_m, slitsize_mr, slitsize_mw
    real    slitsize_x, slitsize_y
    real    xdcorner, pixscl, dx, dy, ra, dec, rabar, decbar
    char    slittype, header, sparea, imarea, mparea, badhdr
    int     slitid, specord, undefined_order, nx, ny
    int     nrows, row0, rown
    real    x1, x2, y1, y2
    int     ix1, ix2, iy1, iy2, ixover, idx, idy, width, half

    cache ("gemextn", "tinfo", "tintegrate", "keypar")

    junk = fscan (mdf, l_mdf)
    l_row = row
    junk = fscan (image, l_image)
    junk = fscan (area, l_area)
    junk = fscan (mapfile, l_mapfile)
    l_pixscale = pixscale
    l_dispaxis = dispaxis
    junk = fscan (logfile, l_logfile)
    junk = fscan (logname, l_logname)
    l_verbose = verbose

    badhdr = ""
    junk = fscan (nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"

    status = 1
    debug = no
    undefined_order = -9999 # -1 and larger values are ok
    sparea = "spectrum"
    imarea = "image"
    mparea = "ifu-map"
    warnedya = no

    if (l_logname == "") l_logname = "NSMDFHELPER"

    if (l_logfile == "") {
        junk = fscan (gnirs.logfile, l_logfile) 
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - "//l_logname//": Both nscombine.logfile \
                and gnirs.logfile are empty.", l_logfile, verbose+) 
            printlog ("                     Using default file " \
                // l_logfile, l_logfile, verbose+) 
        }
    }

    if ("" != badhdr) {
        printlog ("ERROR - " // l_logname // ": Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }

    # Expand image to include PHU by inheritance
    haveheader = (no == ("" == l_image))
    if (haveheader) {
        # (assumes image was a named section, like Sxxx[SCI,1])
        gemextn (l_image, check="exists,mef", process="expand", index="", \
            extname="", ikparams="inherit", omit="", replace="", \
            outfile="STDOUT", logfile="dev$null", glogpars="", \
            verbose=l_verbose) | scan (header)
        if (debug)
            print (header)
        if (no == (1 == gemextn.count) || 0 < gemextn.fail_count) {
            # try direct access
            gemextn (l_image, check="exists", process="none", index="", \
                extname="", ikparams="", omit="", replace="", \
                outfile="STDOUT", logfile="dev$null", glogpars="", \
                verbose=l_verbose) | scan (header)
            if (debug)
                print (header)
            if (no == (1 == gemextn.count) || 0 < gemextn.fail_count) {
                printlog ("ERROR - " // l_logname // ": Problem with \
                image " // image, l_logfile, verbose+)
                goto clean
            }
        }
        if (debug && no == isindef (header))
            print (header)
    } else {
        printlog ("ERROR - " // l_logname // ": No reference image.", \
            l_logfile, verbose+)
        goto clean
    }

    # Read image info if available (priority to original size)
    nx = INDEF; ny = INDEF; pixscl = l_pixscale; axis = l_dispaxis
    hselect (header, "ORIGXSIZ", yes) | scan (nx)
    if (isindef (nx)) {
        hselect (header, "i_naxis1", yes) | scan (nx)
    } else {
        printlog (l_logname // ": Using original X size " // nx, \
            l_logfile, l_verbose)
    }
    hselect (header, "ORIGYSIZ", yes) | scan (ny)
    if (isindef (ny)) {
        hselect (header, "i_naxis2", yes) | scan (ny)
    } else {
        printlog (l_logname // ": Using original Y size " // ny, \
            l_logfile, l_verbose)
    }
    hselect (header, l_key_pixscale, yes) | scan (pixscl)
    hselect (header, l_key_dispaxis, yes) | scan (axis)
    if (debug) {
        hselect (header, "$I,i_naxis2,i_naxis1,"//l_key_dispaxis, yes)
        imhead (header, long-)
        print (pixscl)
        print (axis)
        print (nx)
        print (ny)
    }

    if (isindef (nx) || isindef (ny) \
        || isindef (pixscl) || isindef (axis)) {
        printlog ("ERROR - " // l_logname // ": Missing header data in " \
            // header, l_logfile, verbose+)
        goto clean
   }

    if (debug)
        print (nx // ", " // ny // ", " // pixscl)

    # Auto-detect historical format
    junk = 0
    tlcol (mdf, nlist=1) | match ("slitsize_x") | count | scan (junk)
    historical = (0 == junk)

    # Handle IFU map
    if (mparea == l_area) {
        if (historical) {
            printlog ("ERROR - " // l_logname // ": Cannot generate \
                map file", l_logfile, verbose+)
            printlog ("        from historical MDF format.", \
                l_logfile, verbose+)
            goto clean
        }

        # Output file for GNIRS IFU map
        if ("" == l_mapfile) {
            l_mapfile = mktemp ("mapfile")
            mapfile = l_mapfile
        }
        touch (l_mapfile)

        tinfo (mdf, ttout-)
        nrows = tinfo.nrows
        row0 = 1
        rown = nrows

        tintegrate (mdf, "RA", "", >& "dev$null")
        rabar = tintegrate.integral * 15.0 / real (nrows)
        tintegrate (mdf, "DEC", "", >& "dev$null")
        decbar = tintegrate.integral / real (nrows)

    } else {
        row0 = l_row
        rown = l_row
    }

    for (l_row = row0; l_row <= rown; l_row = l_row + 1) {
        if (historical) {
            # Handle historical format
            if (no == (2 == axis)) {
                printlog ("ERROR - " // l_logname // ": The historical \
                    MDF format is only", l_logfile, verbose+)
                printlog ("        supported for GNIRS instruments \
                    (DISPAXIS = 2)", l_logfile, l_verbose)
                goto clean
            }
            if (no == warnedya) {
                printlog (l_logname // ": Using historical MDF format.", \
                    l_logfile, l_verbose)
                warnedya = yes
            }
            specord = undefined_order
            xdcorner = 0

            tprint (mdf, prparam-, prdata+, showrow-, showhdr-, showunits-, \
                col="x_ccd, y_ccd, slittype, slitid, slitpos_mx, slitpos_my, \
                slitsize_mx, slitsize_my, slittilt_m, slitsize_mr, \
                slitsize_mw, specorder, corner", rows=l_row, pwidth=160) | \
                scan (x_ccd, y_ccd, slittype, slitid, slitpos_mx, slitpos_my, \
                slitsize_mx, slitsize_my, slittilt_m, slitsize_mr, \
                slitsize_mw, specord, xdcorner)

            if (undefined_order == specord) {
                printlog ("ERROR - " // l_logname // ": No data found \
                    in " // mdf, l_logfile, verbose+)
                goto clean
            }
            if (debug)
                print ("(" // x_ccd // "," // y_ccd // ") + " \
                    // slitpos_mx // "; " // slitsize_mx \
                    // "; " // xdcorner)

            if (sparea == strlwr (l_area)) {

                # The following calc follows what I understand from 
                # Inger's email, dated 03-11-11
                # (Which was wrong and inconsistent with the email. Duh)

                # shift x centre to compensate for object offset
                # the sign here may be incorrect, depending on convention
                x_ccd = x_ccd + slitpos_mx / pixscl
                # rewritten to avoid changing by +/- 1 pixel with rounding
                # Adding the 0.000001 as per the discussion below by JH. CL is
                # seeing the integers are very close to the number it should be
                # but is recording them as ~1E-13 below so when 'int' is used
                # it's rounding down. PyRAF is fine. - MS
                ix1 = x_ccd - 0.5 * slitsize_mx / pixscl
                dx = slitsize_mx / pixscl
                ix2 = ix1 + dx - 1 + 0.000001
                # y covers whole chip
                iy1 = 1
                iy2 = ny

                # Clip rectangle (do y anyway, in case algorithm 
                # above changes)
                ixover = 0
                if (ix1 < 1) ix1 = 1
                if (ix2 > nx) {
                    ixover = ix2 - nx
                    ix2 = nx
                }
                if (iy1 < 1) iy1 = 1
                if (iy2 > ny) iy2 = ny

                xlo = ix1
                ylo = iy1
                xhi = ix2
                yhi = iy2

                ixlo = ix1
                iylo = iy1
                ixhi = ix2
                iyhi = iy2
                ixovershoot = ixover

                corner = xdcorner
                specorder = specord
                slitwidth = ixhi - ixlo - corner

            } else {
                printlog ("ERROR - " // l_logname // ": Historical MDF \
                    supports only spectrum mode.", l_logfile, verbose+)
                goto clean
            }

        } else {
            # Handle correct format
            specord = undefined_order
            xdcorner = 0

            tprint (mdf, prparam-, prdata+, showrow-, showhdr-, \
                showunits-, col="x_ccd, y_ccd, slitsize_x, \
                slitsize_y, specorder, corner, RA, DEC", rows=l_row, \
                pwidth=160) | scan (x_ccd, y_ccd, slitsize_x, \
                slitsize_y, specord, xdcorner, ra, dec)

            if (undefined_order == specord) {
                printlog ("ERROR - " // l_logname // ": No data found \
                    in " // mdf, l_logfile, verbose+)
                goto clean
            }

            x1 = x_ccd - 0.5 * slitsize_x / pixscl
            dx = slitsize_x / pixscl
            x2 = x1 + dx
            y1 = y_ccd - 0.5 * slitsize_y / pixscl
            dy = slitsize_y / pixscl
            y2 = y1 + dy

            slitwidth = slitsize_x

            # what we want here is something consistent with the 
            # following:
            # - as small a spatial extent as possible
            #   (requested by tracy)
            # - the same number of pixels for all apertures of a
            #   given width in arcsecs
            #   (needed to let the IFU carry forwards wavelength solns)
            # - centred on the aperture

            # for simplicity, keep floating point (exact) and integer
            # (pixels) calculations completely separate.  the exact
            # values are easy (and already calculated above)

            # should the aperture include partially illuminated 
            # pixels?  i think the real requirement is to guarantee
            # that fully illuminated pixels are always included.  the
            # rest doesn't matter.

            # noting that int() rounds down for +ve floats, the width
            # that guarantees (with correct centring) all full pixels
            # is int (real_width), from considering the case where
            # the aperture starts at exactly one pixel and extends
            # to exactly one pixel (when the expression evaluates to 
            # real_width) or any amount larger than that, but less 
            # than a pixel more.

            # we could position this by simply starting at the first
            # completely filled pixel.  but that would give aperture
            # masks biased in one direction.  instead, we should
            # position by centring.

            # it seems like there are two distinct cases, depending
            # on whether the aperture width corresponds to an odd or
            # or even number of pixels.

            # let's use 
            # int width = int (real_width)
            # int half = width / 2 (integer division, rounds down)

            # the mid point of the aperture is x_ccd 
            # (asuming dispaxis=2).

            # pixel coords in iraf are at the centre of pixels.

            # for odd number widths, this is easy.  find the mid pixel
            # and add/subtract width.  mid pixel is nint (x_ccd)
            # where nint (x) = int (x + 0.5)

            # for even number widths, you have to see the diagram
            # on my whiteboard, which shows that the pixel boundary
            # nearest the centre is int(centre)+0.5 and that the
            # extent of the aperture is therefore
            # int (centre) + 1 - half   to   int (centre) + half

            # so, here goes...
            
            # Solaris is viewing reals as their value minus a
            # a small fraction. When converting to int they get rounded
            # down and that small fraction is changing the values of 
            # the final int...ie. 31 = int(32.) because solaris thinks 
            # 32. is 31.9999999999.
            # Better add a small number to all reals that get converted
            # to ints. (JH)

            dx = dx + 0.000001
            dy = dy + 0.000001
            x_ccd = x_ccd + 0.000001
            y_ccd = y_ccd + 0.000001
            if (2 == axis) width = int (dx)
            else           width = int (dy)
            half = width / 2
            even = (half * 2 == width)

            if (even) {
                if (2 == axis) {
                    ix1 = int (x_ccd) + 1 - half
                    ix2 = int (x_ccd) + half
                } else {
                    iy1 = int (y_ccd) + 1 - half
                    iy2 = int (y_ccd) + half
                }
            } else {
                if (2 == axis) {
                    ix1 = int (x_ccd + 0.5) - half
                    ix2 = int (x_ccd + 0.5) + half
                } else {
                    iy1 = int (y_ccd + 0.5) - half
                    iy2 = int (y_ccd + 0.5) + half
                }
            }

            if (sparea == strlwr (l_area)) {

                # Ignore shifting of object relative to slit
                # (ie slitpos... columns)

                # dispersion direction covers whole chip
                if (2 == axis) {
                    y1 = 1
                    y2 = ny
                    iy1 = 1
                    iy2 = ny
                } else {
                    x1 = 1
                    x2 = nx
                    ix1 = 1
                    ix2 = nx
                }

                xlo = x1
                ylo = y1
                xhi = x2
                yhi = y2

                ixover = 0
                if (ix1 < 1) ix1 = 1
                if (ix2 > nx) {
                    ixover = ix2 - nx
                    ix2 = nx
                }
                if (iy1 < 1) iy1 = 1
                if (iy2 > ny) iy2 = ny

                ixlo = ix1
                iylo = iy1
                ixhi = ix2
                iyhi = iy2
                ixovershoot = ixover

                corner = xdcorner
                specorder = specord
                slitwidth = ixhi - ixlo - corner

            } else if (imarea == strlwr (l_area) || \
                mparea == strlwr (l_area)) {

                # Ignore shifting of object relative to slit
                # (ie slitpos... columns)

                xlo = x1
                ylo = y1
                xhi = x2
                yhi = y2
                
                if (2 == axis) width = int (dy)
                else           width = int (dx)
                half = width / 2
                even = (half * 2 == width)

                if (even) {
                    if (2 == axis) {
                        iy1 = int (y_ccd) + 1 - half
                        iy2 = int (y_ccd) + half
                    } else {
                        ix1 = int (x_ccd) + 1 - half
                        ix2 = int (x_ccd) + half
                    }
                } else {
                    if (2 == axis) {
                        iy1 = int (y_ccd + 0.5) - half
                        iy2 = int (y_ccd + 0.5) + half
                    } else {
                        ix1 = int (x_ccd + 0.5) - half
                        ix2 = int (x_ccd + 0.5) + half
                    }
                }

                ixover = 0
                if (ix1 < 1) ix1 = 1
                if (ix2 > nx) {
                    ixover = ix2 - nx
                    ix2 = nx
                }
                if (iy1 < 1) iy1 = 1
                if (iy2 > ny) iy2 = ny

                ixlo = ix1
                iylo = iy1
                ixhi = ix2
                iyhi = iy2
                ixovershoot = ixover

                corner = xdcorner
                specorder = specord
                slitwidth = ixhi - ixlo - corner

            } else {
                printlog ("ERROR - " // l_logname // ": Unexpected area " \
                    // l_area, l_logfile, verbose+)
                goto clean
            }
        }

        if (mparea == l_area) {
            idx = ix2 - ix1 + 1
            idy = iy2 - iy1 + 1

            x1 = (ra * 15.0 - rabar) * 3600.0
            x1 = x1 * cos (decbar * 3.14159 / 180.0)
            x1 = x1 / pixscl
            x1 = x1 - (real (idx) / 2.0) + (real (nx) / 2.0)
            y1 = (dec - decbar) * 3600.0
            y1 = y1 / pixscl
            y1 = y1 - (real (idy) / 2.0) + (real (ny) / 2.0)

            if (debug)
                print (y1 // " " // int (y1))

            # the -0.1 below is a hack that avoids (hopefully)
            # moire-like interference/bug when DEC evenly spaced
            # by 1/n

            print (ix1 // " " // iy1 // " " // idx // " " // idy // " " \
                // int (x1-0.1) // " " // int (y1-0.1), >> l_mapfile)
        }
    }

    # Completed OK
    status = 0

clean:
    if (no == (0 == status))
        delete (l_mapfile, verify-, >& "dev$null")
    
end
