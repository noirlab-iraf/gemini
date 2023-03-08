# Copyright(c) 2005-2017 Association of Universities for Research in Astronomy, Inc.
#
# GMSKIMG - Produce a fake image for use with gmmps
#
# Original author: Rachel Johnson
#
# User inputs image to be transformed, R.A. and Dec of input objects, x,y
# of input objects, the required hemisphere, the R.A. and Dec of the required
# field centre and the (OT) position angle for the required field.
# The input image must have WCS values in the header, namely CRPIX, CRVAl and
# the CD matrix. The routine uses these WCS values to turn the input R.A.,Dec
# into x,y in the input image. It then finds the mapping between input x,y and
# x,y and uses this to transform the image. This means that for a
# transformed image where the input objects are in the correct place the WCS
# of the input image must correspond to the R.A.,Dec of the input objects -
# i.e. must be such that the input x,y calculated for these objects is correct.
# It is intended that this routine is called by gmskcreate.cl

# 10/07 A. Wong: Added support for image transformation using predetermined
#                  R.A./Dec object grid
#                Removed hemisphere, replaced with selection for instrument
#                Added place holder FLAMINGOS-2, will exit if FLAMINGOS-2 is
#                  selected, all FLAMINGOS-2 related coefficients are set to 0.
#                Added instrument dependency for size of CCD for geotran
#                Moved WCS check and header keyword edits from gmskcreate to
#                  here
# 02/08 A. Wong: Changed value for added header keyword DATE-OBS to 1900-01-01
#                  for GSA ingestion purposes.
# Note the array size of the x,y is currently 1000
#
#
# 2017-03-28 mischa:
#                -- Updated detector geometries and pixel scale for GMOS-S Hamamatsu and F2.
#                -- Added a 25 point circular grid for F2 with 3 arcmin radius
#                -- Default 2x2 binning for gmos
#                -- Writing DETTYPE keyword, required by gmscreate;
# 2017-04-21 mischa: Added pseudo-image support for F2
#

procedure gmskimg (inimage, gprgid, instrument, rafield, decfield, pa)

char    inimage     {prompt="Image to be transformed (if MEF, specify extension)"}
char    gprgid      {prompt="Your Gemini program ID (e.g. GN-2007A-Q-4)"}
char    instrument  {enum="gmos-n|gmos-s|flamingos2", prompt="Instrument (gmos-n|gmos-s|flamingos2)"}
real    rafield     {prompt="R.A. value of field center"}
real    decfield    {prompt="Dec value of field center"}
real    pa          {prompt="PA of field if required by OT"}
char    fraunits    {"hours", enum="hours|degrees", prompt="Field centre R.A. units (hours|degrees)"}
char    outimage    {"",prompt="Output transformed image"}
bool    fl_inter    {no, prompt="Interactive mode"}
char    logfile     {"", prompt="Logfile name"}
bool    verbose     {yes, prompt="Verbose"}
bool    fl_debug    {no, prompt="Print debugging information"}
int     status      {0, prompt="Exit status (0=good)"}

struct  *scanfile1   {"", prompt="Internal use only"}
struct  *scanfile2   {"", prompt="Internal use only"}

begin

    char    l_inimage = ""
    char    l_gprgid = ""
    char    l_instrument = ""
    char    l_logfile = ""
    char    l_fraunits = ""
    char    l_outimage = ""
    bool    l_fl_inter
    bool    l_verbose
    bool    l_fl_debug
    real    l_rafield
    real    l_decfield
    real    l_pa

    char    instrument_hdr, paramstr, errmsg, tmpgrid, tmpxy, tmpskxy, tmpcrds
    char    tmpastro, tmpxyxy, tmpdb, tmpstr, tmplog, s1, s2, s3, s4, id
    int     junk, pri, lowrow, lowcol, highrow, highcol, x1, x2, y1, y2
    real    pi, ra0, dec0, x, y, radiff, xi, eta, pxi, peta, theta, pixsc, bin
    real    xgemini, ygemini, ra, dec, mag, ssx, ssy, st, spy
    struct  tmpstruct

    # outpref is applied regardless of whether outimage is set or not

    char    outpref = "GMI"

    # Coefficients used to generate coordinate grid to do image transformation.
    # These coefficients represent the arcsec offsets
    real    xgrid_gmos[25] = 0, 0, 50, 50, 50, 0, -50, -50, -50, 0, 100, \
        100, 100, 0, -100, -100, -100, 0, 150, 150, 150, 0, -150, -150, -150
    real    ygrid_gmos[25] = 0, 50, 50, 0, -50, -50, -50, 0, 50, 100, 100, 0, \
        -100, -100, -100, 0, 100, 150, 150, 0, -150, -150, -150, 0, 150

    real    xgrid_flamingos[25] = 0, 0, 42, 60, 42, 0, -42, -60, -42, 0, 84, \
	    120, 84, 0, -84, -120, -84, 0, 126, 180, 126, 0, -126, -180, -126
    real    ygrid_flamingos[25] = 0, 60, 42, 0, -42, -60, -42, 0, 42, 120, 84, \
	    0, -84, -120, -84, 0, 84, 180, 126, 0, -126, -180, -126, 0, 126

    real    xgrid[25]
    real    ygrid[25]

    real tmpgridx, tmpgridy

    junk = fscan (inimage, l_inimage)
    junk = fscan (gprgid, l_gprgid)
    junk = fscan (instrument, l_instrument)
    l_rafield = rafield
    l_decfield = decfield
    l_pa = pa
    junk = fscan (fraunits, l_fraunits)
    junk = fscan (outimage, l_outimage)
    l_fl_inter = fl_inter
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose
    l_fl_debug = fl_debug

    # Initialize

    status = 0
    pi = 4 * atan2(1,1)

    # Create names for temp files

    tmpgrid = mktemp ("tmpgrid")
    tmpxy = mktemp ("tmpxy")
    tmpskxy = mktemp ("tmpskxy")
    tmpcrds = mktemp ("tmpcrds")
    tmpastro = mktemp ("tmpastro")
    tmpxyxy = mktemp ("tmpxyxy")
    tmpdb = mktemp ("tmpdb")

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.

    paramstr =  "inimage        = "//inimage.p_value//"\n"
    paramstr += "gprgid         = "//gprgid.p_value//"\n"
    paramstr += "instrument     = "//instrument.p_value//"\n"
    paramstr += "rafield        = "//rafield.p_value//"\n"
    paramstr += "decfield       = "//decfield.p_value//"\n"
    paramstr += "pa             = "//pa.p_value//"\n"
    paramstr += "fraunits       = "//fraunits.p_value//"\n"
    paramstr += "outimage       = "//outimage.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value//"\n"
    paramstr += "fl_debug       = "//fl_debug.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.

    gloginit (l_logfile, "gmskimg", "mostools", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Does the Gemini program ID have a valid structure?

    print (l_gprgid) | translit ("STDIN", "-", " ") | scan (s1,s2,s3,s4)

    if ( ((s1 == "GN") || (s1 == "GS")) && \
        (substr(s2,1,1) == "2") && \
        ((substr(s2,5,5) == "A") || (substr(s2,5,5) == "B")) && \
        ((s3 == "Q") || (s3 == "C") || (s3 == "SV") || (s3 == "DD") || \
         (s3 == "ENG") || (s3 == "LP") || (s3 == "FT")) ) {

        glogprint (l_logfile, "gmskimg", "status", type="string",
            str="Gemini program ID "//gprgid, verbose=l_verbose)

    } else {
        errmsg = "The Gemini program ID is not valid."
        status = 121
        glogprint (l_logfile, "gmskimg", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto exit
    }

    if (l_outimage == "")
        l_outimage = outpref//l_inimage

    fparse (l_outimage)
    l_outimage = fparse.directory//fparse.root//fparse.extension

    # Define CCD size for given instrument and grid coordinates

    if (l_instrument == "gmos-n") {
  	    lowrow  = 1
    	lowcol  = 1
        highrow = 4176
        highcol = 6278
        pixsc   = 0.0807
		# New default: 2x2 binning for GMOS
        bin     = 2.
        pixsc   = pixsc * bin
        highrow = highrow / bin
        highcol = highcol / bin		

        for (i=1; i<26; i=i+1) {
            xgrid[i] = xgrid_gmos[i]
            ygrid[i] = ygrid_gmos[i]
        }

    } else if (l_instrument == "gmos-s") {
        lowrow  = 1
        lowcol  = 1
        highrow = 4176
        highcol = 6264
        pixsc   = 0.0800
		# New default: 2x2 binning for GMOS
        bin     = 2.
        pixsc   = pixsc * bin
        highrow = highrow / bin
        highcol = highcol / bin		

        for (i=1; i<26; i=i+1) {
            xgrid[i] = xgrid_gmos[i]
            ygrid[i] = ygrid_gmos[i]
        }

    } else if (l_instrument == "flamingos2") {
        lowrow  = 1
        lowcol  = 1
        highrow = 2048
        highcol = 2048
        pixsc   = 0.1792

        for (i=1; i<26; i=i+1) {
            xgrid[i] = xgrid_flamingos[i]
            ygrid[i] = ygrid_flamingos[i]
        }
    }
    glogprint (l_logfile, "gmskimg", "status", type="string",
        str="Using CCD, pixel scale and grid values for "//l_instrument,
        verbose=l_verbose)

    ######## check that the input image has the WCS set correctly ############
    # The WCS is required by gmskimg to turn R.A.,Dec into original image x,y
    # The required header keywords are:
    #       CRPIX1, CRPIX2, CRVAL1, CRVAL2, CD1_1, CD1_2, CD2_1, CD2_2

    errmsg = ""
    imgets (l_inimage, 'CRPIX1')
    if (imgets.value == "0") {
        errmsg = "CRPIX1 "
        status = 131
    }
    imgets (l_inimage, 'CRPIX2')
    if (imgets.value == "0") {
        errmsg += "CRPIX2 "
        status = 131
    }
    imgets (l_inimage, 'CRVAL1')
    if (imgets.value == "0") {
        errmsg += "CRVAL1 "
        status = 131
    }
    imgets (l_inimage, 'CRVAL2')
    if (imgets.value == "0") {
        errmsg += "CRVAL2 "
        status = 131
    }
    imgets (l_inimage, 'CD1_1')
    if (imgets.value == "0") {
        errmsg += "CD1_1 "
        status = 131
    }
    imgets (l_inimage, 'CD1_2')
    if (imgets.value == "0") {
        errmsg += "CD1_2 "
        status = 131
    }
    imgets (l_inimage, 'CD2_1')
    if (imgets.value == "0") {
        errmsg += "CD2_1 "
        status = 131
    }
    imgets (l_inimage, 'CD2_2')
    if (imgets.value == "0") {
        errmsg += "CD2_2 "
        status = 131
    }
    if (status != 0) {
        glogprint (l_logfile, "gmskimg", "status", type="error",
            errno=status, str="The following WCS parameters were missing \
            in the input image header: "//errmsg, verbose+)
        goto exit
    }

    ###################### end of WCS check ##################################

    glogprint (l_logfile, "gmskimg", "status", type="string",
        str="Generating object grid for transformation", verbose=l_verbose)

    if (l_fraunits == "hours") {

        # Convert to degrees and use DMS for all coordinates to avoid
        # confusion with units

        l_rafield = l_rafield * (360./24.)
    }

    for (i=1; i<26; i=i+1) {

        tmpgridx = l_rafield + \
            (xgrid[i]/3600.)/cos((l_decfield + (ygrid[i]/3600.))*(2.*pi)/360.)
        tmpgridy = l_decfield + (ygrid[i]/3600.)

        print (i, tmpgridx, tmpgridy, 28, 1, 1.0, 5.0, 0, 0, >> tmpgrid)

        if (l_fl_debug) {
            print (i, tmpgridx, tmpgridy) | scan (tmpstruct)
            glogprint (l_logfile, "gmskimg", "engineering", type="string",
                str=tmpstruct, verbose=yes)
        }
    }

    glogprint (l_logfile, "gmskimg", "status", type="fork", fork="forward",
        child="gmskxy", verbose=l_verbose)

    # Determine the x and y pixel coordinates of the R.A. and Dec grid for
    # the pseudo image

    gmskxy (tmpgrid, l_instrument, l_rafield, l_decfield, l_pa,
        iraunits="degrees", fraunits="degrees", outcoords=tmpxy,
        outtab=tmpskxy, logfile=l_logfile, verbose=l_verbose,
        fl_debug=l_fl_debug, status=status)

    glogprint (l_logfile, "gmskimg", "status", type="fork", fork="backward",
        child="gmskxy", verbose=l_verbose)

    scanfile1 = tmpgrid
    scanfile2 = tmpxy

    # j counts the total number of objects used for the transformation
    # (needs to be >6)

    j = 0

    while ((fscan (scanfile1, id, ra, dec, mag, pri, ssx, ssy, st, spy)) != \
        EOF) {

        # Determine the x and y pixel coordinates of the R.A. and Dec grid for
        # inimage

        rd2xy (l_inimage, ra, dec, hour=no) | scan(tmpstruct)

        junk = fscanf (tmpstruct, "%s %s %g %s %s %s %g", tmpstr, tmpstr, x,
            tmpstr, tmpstr, tmpstr, y)

        if (l_fl_debug) {
            printf ("Source #%d : R.A.(degrees)=%.6f , Dec=%.6f, \
                (x,y)=(%.2f,%.2f)\n", i, ra, dec, x, y) | scan (tmpstruct)
            glogprint (l_logfile, "gmskimg", "engineering", type="string",
                str=tmpstruct, verbose=yes)
        }

        # Create the list of x and y pixel coordinates for geomap

        junk = fscan (scanfile2, xgemini, ygemini)

        print (xgemini, ygemini, x, y, >> tmpxyxy)

        # All the x and y pixel coordinates of the R.A. and Dec grid determined
        # for the pseudo image will lie within the bounds of the pseudo image
        # by definition. Check whether the x and y pixel coordinates of the
        # R.A. and Dec grid determined for inimage lies within the bounds of
        # inimage.

        x1 = 1
        keypar (l_inimage, "i_naxis1", silent+)
        if (keypar.found)
            x2 = int (keypar.value)
        else
            glogprint (l_logfile, "gmskimg", "status", type="error",
                errno=status, str="NAXIS1 keyword not found", verbose+)

        y1 = 1
        keypar (l_inimage, "i_naxis2", silent+)
        if (keypar.found)
            y2 = int (keypar.value)
        else
            glogprint (l_logfile, "gmskimg", "status", type="error",
                errno=status, str="NAXIS2 keyword not found", verbose+)

        # Default object grid provides 25 objects, at least 6 must fall in the
        # field of view

        if ((x >= x1) && (x <= x2) && (y >= y1) && (y <= y2)) {
            j += 1
        }
    }
    scanfile1 = ""
    scanfile2 = ""

    if (j < 6) {
        status = 121
        glogprint (l_logfile, "gmskimg", "status", type="string",
            str="Less than 6 objects fell on to the field of view of "//\
            l_instrument//".", verbose=l_verbose)
        glogprint (l_logfile, "gmskimg", "status", type="string",
            str="Your field center falls outside the field of view of your \
            input image, please redefine your field center. ",
            verbose=l_verbose)
        glogprint (l_logfile, "gmskimg", "status", type="string",
            str="Geomap can not determine a proper fit to do image \
            transformation.", verbose=l_verbose)
        glogprint (l_logfile, "gmskimg", "status", type="string",
            str="No pseudo-"//l_instrument//" image will be created",
            verbose=l_verbose)
        glogprint (l_logfile, "gmskimg", "status", type="error", errno=status,
            verbose=yes)
        goto exit
    } else {
        glogprint (l_logfile, "gmskimg", "status", type="string",
            str="Using " //j// " objects to determine transformation using \
            geomap", verbose=l_verbose)
    }

    # Determine the transformation between the x and y pixel coordinates in the
    # pseudo image (determined from the R.A. and Dec grid) and the x and y
    # pixel coordinates in inimage (determined from the same R.A. and Dec grid)

    glogprint (l_logfile, "gmskimg", "visual", type="visual", vistype="empty",
        verbose=l_verbose)
    glogprint (l_logfile, "gmskimg", "status", type="string", str="Finding \
        transformation from x,y ->  x,y to create fake "//l_instrument//" \
        image", verbose=l_verbose)

    tmplog = mktemp ("tmplog")
    geomap (tmpxyxy, tmpdb, lowcol, highcol, lowrow, highrow,
        transforms="imtrans", results="", fitgeometry="general",
        function="polynomial", xxorder=2, yxorder=2,xxterms="half", xyorder=2,
        yyorder=2, yxterms="half", maxiter=0, reject=3., calctype="real",
        verbose=l_verbose, interactive=l_fl_inter, >& tmplog)
    glogprint (l_logfile, "gmskimg", "science", type="file", str=tmplog,
        verbose=l_verbose)
    delete (tmplog, verify-, >& "dev$null")

    # Now transform the image

    glogprint (l_logfile, "gmskimg", "visual", type="visual", vistype="empty",
        verbose=l_verbose)
    glogprint (l_logfile, "gmskimg", "task", type="string",
        str="Transforming image...", verbose=l_verbose)

    tmplog = mktemp ("tmplog")
    geotran (l_inimage, l_outimage, tmpdb, "imtrans", geometry="geometric",
        xin=INDEF, yin=INDEF, xshift=INDEF, yshift=INDEF, xout=INDEF,
        yout=INDEF, xmag=INDEF, ymag=INDEF, xrotation=INDEF, yrotation=INDEF,
        xmin=INDEF, xmax=INDEF, ymin=INDEF, ymax=INDEF, xscale=1.0, yscale=1.0,
        ncols=INDEF, nlines=INDEF, xsample=1.0, ysample=1.0,
        interpolant="linear", boundary="constant", constant=0.0,
        fluxconserve=yes, nxblock=512, nyblock=512, verbose=l_verbose,
        >& tmplog)

    glogprint (l_logfile, "gmskimg", "engineering", type="file", str=tmplog,
        verbose=l_verbose)
    delete (tmplog, verify-, >& "dev$null")

    ## Adding Gemini required keywords to the image header ##
    # R.A. and DEC in header are in decimal degrees. Dec is already in this
    # format. Convert R.A. if necessary. Checks to see that gmskimg completed
    # properly before adding headers to output image

    if ((l_instrument == "gmos-n") || (l_instrument == "gmos-s")) {
        # Mask tracking database requires upper case letters for GMOS
        instrument_hdr = strupr(l_instrument)
    } else if (l_instrument == "flamingos2") {
        # GMMPS requires "F2"
        instrument_hdr = "F2"
    }

    gemhedit (l_outimage, "INSTRUME", instrument_hdr,
        "Instrument used to acquire data", delete-)
    gemhedit (l_outimage, "OBJECT", "pseudo"//l_instrument, "Object name",
        delete-)
    gemhedit (l_outimage, "GEMPRGID", l_gprgid, "Gemini program ID", delete-)
    gemhedit (l_outimage, "RA", l_rafield, "Right Ascension", delete-)
    gemhedit (l_outimage, "DEC", l_decfield, "Declination of Target", delete-)
    gemhedit (l_outimage, "DATE-OBS", "1900-01-01",
        "UT Date of observation (YYYY-MM-DD)", delete-)
    gemhedit (l_outimage, "TIME-OBS", "dummy", "Time of observation", delete-)
    gemhedit (l_outimage, "TIME-OBS", "00:00:00", "Time of observation",
        delete-)
    # GMOS and F2 pseudo-images have different binning:
    # Also NEED DETTYPE keyword for odf2mdf
    if (l_instrument == "gmos-n") {
        gemhedit (l_outimage, "CCDSUM", "2 2", "CCD sum", delete-)
        gemhedit (l_outimage, "DETTYPE", "S10892-N", "Detector type", delete-)
    } else if (l_instrument == "gmos-s") {
        gemhedit (l_outimage, "CCDSUM", "2 2", "CCD sum", delete-)
        gemhedit (l_outimage, "DETTYPE", "S10892", "Detector type", delete-)
    } else if (l_instrument == "flamingos2") {
        gemhedit (l_outimage, "CCDSUM", "1 1", "CCD sum", delete-)
        gemhedit (l_outimage, "DETTYPE", "Hawaii-2 2048", "Detector type", delete-)
    }
    gemhedit (l_outimage, "PIXSCALE", pixsc, "Pixel scale", delete-)

    glogprint (l_logfile, "gmskimg", "status", type="string", str="Required \
        Gemini header keywords added to output transformed "//\
        l_instrument//" image", verbose=l_verbose)

    ############# end output header update ################################

exit:
    scanfile1 = ""
    scanfile2 = ""
    delete (tmpgrid // "," // tmpxy // "," // tmpskxy // ".tab," // tmpcrds //
        "," // tmpastro // "," // tmpxyxy // "," // tmpdb, verify-,
        >& "dev$null")

    if (status == 0)
        glogclose (l_logfile, "gmskimg", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "gmskimg", fl_success-, verbose=l_verbose)

exitnow:
    ;

end
