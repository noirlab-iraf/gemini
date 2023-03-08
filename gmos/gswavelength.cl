# Copyright(c) 2002-2015 Association of Universities for Research in Astronomy, Inc.

procedure gswavelength (inimages)

# Wavelength calibrate GMOS spectra
#
# Version   Feb 28, 2002  BM,IJ   v1.3 release
#           Aug 27, 2002  IJ     check for required gmosaic, report exit status
#                                fix the naming of the calibrations (mdfrow)
#           Aug 28, 2002  IJ     make the mdfrow thing work for IFU/Longslit
#           Sept 20, 2002    v1.4 release
#           Oct 14, 2002  IJ  don't use ijk from within script

string  inimages    {prompt="Input images"}
string  crval       {"CRVAL1",prompt="Approximate wavelength at coordinate reference pixel"}
string  cdelt       {"CD1_1",prompt="Approximate dispersion"}
string  crpix       {"CRPIX1",prompt="Coordinate reference pixel"}
char    key_dispaxis    {"DISPAXIS",prompt="Header keyword for dispersion axis"}
int     dispaxis    {1,min=1,max=2,prompt="Dispersion axis"}
#
char    database    {"database",prompt="Directory for files containing feature data"}
char    coordlist   {"gmos$data/CuAr_GMOS.dat",prompt="User coordinate list, line list"}
string  gratingdb   {"gmos$data/GMOSgratings.dat",prompt="Gratings database file"}
string  filterdb    {"gmos$data/GMOSfilters.dat",prompt="Filters database file"}
char    fl_inter    {"yes",enum="yes|no|YES|NO",prompt="Examine identifications interactively"}
#
char    section     {"default",prompt="Image section for running identify"}
int     nsum        {10,min=1,prompt="Number of lines or columns to sum"}
char    ftype       {"emission",min="emission|absorption",prompt="Feature type"}
real    fwidth      {10.,min=2.,prompt="Feature width in pixels"}
real    gsigma      {1.5,min=0.0,prompt="Gaussian sigma for smoothing"}
real    cradius     {12.,min=2.,prompt="Centering radius in pixels"}
real    threshold   {0.,prompt="Feature threshold for centering"}
real    minsep      {5.,prompt="Minimum pixel separation for features"}
real    match       {-6.,prompt="Coordinate list matching limit, <0 pixels, >0 user"}
char    function    {"chebyshev",min="legendre|chebyshev|spline1|spline3",prompt="Coordinate fitting function"}
int     order       {4,min=1,prompt="Order of coordinate fitting function"}
char    sample      {"*",prompt="Coordinate sample regions"}
int     niterate    {10,min=0,prompt="Rejection iterations"}
real    low_reject  {3.,min=0,prompt="Lower rejection sigma"}
real    high_reject {3.,min=0,prompt="Upper rejection sigma"}
real    grow        {0.,min=0,prompt="Rejection growing radius"}
bool    refit       {yes,prompt="Refit coordinate function when running reidentify"}
int     step        {10,prompt="Steps in lines or columns for reidentification"}
bool    trace       {yes,prompt="Use fit from previous step rather than central aperture"}
int     nlost       {15,min=0,prompt="Maximum number of lost features"}
int     maxfeatures {150,min=3,prompt="Maximum number of features"}
int     ntarget     {30,min=3,prompt="Number of features used for autoidentify"}
int     npattern    {5,min=3,max=10,prompt="Number of features used for pattern matching (autoidentify)"}
bool    fl_addfeat  {yes,prompt="Allow features to be added by reidentify"}
char    aiddebug    {"",prompt="Debug parameter for aidpars"}
char    fl_dbwrite  {"YES",enum="yes|no|YES|NO",prompt="Write results to database?"}
bool    fl_overwrite    {yes,prompt="Overwrite existing database entries?"}
bool    fl_gsappwave    {no,prompt="Run GSAPPWAVE on all images?"}
#
string  fitcfunc    {"chebyshev",enum="chebyshev|legendre",prompt="Function for fitting coordinates"}
int     fitcxord    {4,prompt="Order of fitting function in X-axis"}
int     fitcyord    {4,prompt="Order of fitting function in Y-axis"}
string  logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose output?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

begin

    string  l_inimages
    char    l_section, l_database, l_coordlist, l_units, l_ftype
    char    l_key_specsec, l_specsec, l_filterdb, l_gratingdb
    char    l_function, l_sample
    char    l_key_filter, l_key_fpmask, l_nsappwavedb, l_npreparedb
    char    l_sci_ext, l_logfile
    bool    l_refit, l_trace
    bool    l_fl_trimmed, l_fl_overwrite, l_verbose, l_fl_gsappwave
    bool    l_fl_addfeat
    char    l_crval, l_cdelt, l_crpix
    real    l_fwidth, l_gsigma, l_cradius, l_threshold, l_minsep, l_match
    real    l_low_reject, l_high_reject, l_grow
    int     l_dispaxis, l_nsum, l_order, l_niterate, l_step, l_nlost
    int     l_maxfeat, l_ntarget, l_npattern
    char    l_key_dispaxis, l_aiddebug
    char    l_fl_inter, l_fl_dbwrite
    string  l_fitcfunc
    int     l_fitcxord, l_fitcyord, nbad, n_i, atpos

    string  obsmode, imgroot, snum, imgext, inlist, infiles, img
    string  tmpcoord_compare, tmpdispcor, tmp2dispcor, fitcoord_log
    string  cd1_1_ident
    int     nexten, sel, mdfrow, nxpix, nypix
    real    crpix1_gsapp, crval1_gsapp, cd1_1_gsapp
    bool    flag_int, compare
    struct  sdate

    # Query parameters
    l_inimages = inimages
    l_crval = crval
    l_cdelt = cdelt
    l_crpix = crpix

    l_units = "Angstroms"
    l_dispaxis = dispaxis
    l_key_dispaxis = key_dispaxis
    l_section = section
    if ((l_section=="default" || l_section=="" || l_section==" ") \
        && l_dispaxis==1)
            l_section = "middle line"
    if ((l_section=="default" || l_section=="" || l_section==" ") \
        && l_dispaxis==2)
            l_section = "middle column"

    l_database=database ;l_gratingdb=gratingdb ; l_filterdb=filterdb
    l_coordlist=coordlist ; l_fl_inter=fl_inter ; l_nsum=nsum
    l_ftype=ftype ; l_fwidth=fwidth ; l_cradius=cradius ; l_threshold=threshold
    l_minsep=minsep ; l_match=match ; l_function=function ; l_order=order
    l_sample=sample ; l_niterate=niterate
    l_low_reject=low_reject ; l_high_reject=high_reject ; l_grow=grow
    l_refit=refit ; l_step=step ; l_trace=trace ; l_nlost=nlost
    l_fl_dbwrite=fl_dbwrite ; l_fl_overwrite=fl_overwrite
    l_aiddebug=aiddebug ; l_fl_gsappwave=fl_gsappwave
    l_logfile=logfile
    l_verbose=verbose
    l_fitcfunc=fitcfunc ; l_fitcxord=fitcxord ; l_fitcyord=fitcyord
    l_gsigma=gsigma

    l_maxfeat=maxfeatures ; l_ntarget=ntarget ; l_npattern=npattern
    l_fl_addfeat=fl_addfeat

    flag_int = no
    if (l_fl_inter == "yes" || l_fl_inter == "YES") flag_int = yes

    # Set compare flag
    compare = no

    # Reset status flag
    status = 0

    # Keep imgets parameters from changing by outside world
    cache ("imgets", "fparse", "gemdate")

    # Check logfile
    if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GSWAVELENGTH: Both gswavelength.logfile and \
                gmos.logfile fields are empty", l_logfile, yes)
            printlog ("                    Using default file gmos.log",
                l_logfile, yes)
        }
    }

    # Temporary files
    infiles = mktemp("tmpinfiles")

    # Start logging to file
    if (l_logfile == "STDOUT") {
        l_logfile = ""
        l_verbose = yes
    }

    date | scan(sdate)
    printlog ("--------------------------------------------------------------\
        ------------------", l_logfile, verbose=l_verbose)
    printlog ("GSWAVELENGTH -- "//sdate, l_logfile, verbose=l_verbose)
    printlog (" ", l_logfile, verbose=l_verbose)
    printlog ("inimages = "//l_inimages, l_logfile, verbose=l_verbose)

    #check that there are input files
    if (l_inimages == "" || l_inimages == " ") {
        printlog ("ERROR - GSWAVELENGTH: input files not specified",
            l_logfile, yes)
        goto error
    }

    # check existence of list file
    atpos = strstr("@", l_inimages)
    if (atpos >= 1) {
        inlist = substr (l_inimages, atpos + 1, strlen(l_inimages))
        if (!access(inlist)) {
            printlog ("ERROR - GSWAVELENGTH: Input list "//inlist//\
                " not found",
                l_logfile, verbose+)
            goto error
        }
    }

    # Loop over inimages
    sections (l_inimages, > infiles)
    nbad = 0
    scanfile = infiles
    while (fscan(scanfile, img) != EOF) {
        # check if input file is MEF
        gimverify (img)
        if (gimverify.status > 0) {
            printlog ("ERROR - GSWAVELENGTH: "//img//" does not exist or is \
                not a MEF", l_logfile, yes)
            nbad+=1
        }
        imgroot = gimverify.outname
        img = imgroot//".fits"
        if (gimverify.status==0) {
            keypar (imgroot//"[0]", "GMOSAIC", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GSWAVELENGTH: "//img//" has not been \
                    processed with GMOSAIC", l_logfile, yes)
                nbad+=1
            }
        }

    }
    scanfile = ""

    if (nbad > 0)
        goto error

    scanfile = infiles
    while (fscan(scanfile, img) != EOF) {
        gimverify (img)
        imgroot = gimverify.outname
        img = imgroot//".fits"

        # check that gsappwave has been called
        imgets (img//"[0]", "GSAPPWAV", >& "dev$null")
        if (imgets.value=="0" || l_fl_gsappwave) {
            gsappwave (img, logfile=l_logfile, gratingdb=l_gratingdb,
                filterdb=l_filterdb, key_dispaxis=l_key_dispaxis,
                dispaxis=l_dispaxis, verbose=l_verbose)
            if (gsappwave.status != 0) {
                goto error
            }
        }

        # get header information
        imgets (img//"[0]","NSCIEXT")
        nexten = int(imgets.value)

        imgets (img//"[0]", "OBSMODE")
        obsmode = imgets.value

        #Loop over extensions
        for (n_i=1; n_i<=nexten; n_i+=1) {
            printlog ("Calibrating extension: "//n_i, l_logfile, l_verbose)
            sel = 1
            if (obsmode == "MOS") {
                imgets (img//"[SCI,"//n_i//"]", "MDFROW", >& "dev$null")
                if (imgets.value != "0") {
                    mdfrow = int(imgets.value)
                    tabpar (img//"[MDF]", "SELECT", mdfrow)
                    if (tabpar.undef==no) {
                        sel = int(tabpar.value)
                    }
                    printlog ("MDF row: "//mdfrow, l_logfile, l_verbose)
                }
            } else {
                mdfrow = n_i # Longslit and IFU cases
            }
            if (sel==1) {
                #snum = "000"+mdfrow
                printf ("%03d\n", mdfrow) | scan(snum)
                imgext = imgroot//"_"//snum
                # delete temp image if it exists
                if (imaccess(imgext)) {
                    imdelete (imgext, verify-)
                }
                imcopy (img//"[SCI,"//n_i//"]", imgext, verbose-)


                if (compare) {
                    printlog ("--------- COMPARE 1 ---------", \
                        l_logfile, verbose+)

                    printlog ("[SCI,"//n_i//"]", l_logfile, verbose+)

                    # Read some parameters for the compare statement after the
                    # solution is found only for 2-D spectra...
                    hselect (imgext//"[0]", "i_naxis1,i_naxis2,CRPIX1,"//\
                        "CRVAL1,CD1_1", yes) | scan (nxpix, nypix, \
                        crpix1_gsapp, crval1_gsapp, cd1_1_gsapp)

                    # Create tempory file for the pixel coordinates to evaluate
                    tmpcoord_compare = mktemp("tmpcoord_compare")

                    print (\
                        "1 "//nint(nypix/2)//"\n"//\
                        nxpix//" "//nint(nypix/2)//"\n"//\
                        crpix1_gsapp//" "//nint(nypix/2)//"\n"//\
                        "1 1\n"//\
                        nxpix//" 1\n"//\
                        "1 "//nypix//"\n"//\
                        nxpix//" "//nypix,\
                        >> tmpcoord_compare)

                    printlog ("GSAPPWAVE CRPIX1: "//crpix1_gsapp//\
                        " CRVAL1: "//crval1_gsapp//" CD1_1: "//cd1_1_gsapp, \
                        l_logfile, verbose+)
                    printlog ("    Calculated wavelength values for ends of"//\
                        " spectrum at middle line of image, "//\
                        "\n    using GSAPPWAVE solution:",\
                        l_logfile, verbose+)
                    printlog ("        X Y Wavelength[Ang] Y", \
                        l_logfile, verbose+)
                    printlog ("        1 "//nint(nypix/2)//" "//\
                        real(crval1_gsapp - ((crpix1_gsapp - 1) * \
                        cd1_1_gsapp))//" 1.0", \
                        l_logfile, verbose+)
                    printlog ("        "//nxpix//" "//nint(nypix/2)//" "//\
                        real(crval1_gsapp - ((crpix1_gsapp - nxpix) * \
                        cd1_1_gsapp))//" 1.0", \
                        l_logfile, verbose+)
                    printlog ("--------- COMPARE 1 ---------", \
                        l_logfile, verbose+)
                }

                # smooth spectrum
                if (l_gsigma > 0.0) {
                    gauss (imgext,imgext, l_gsigma, ratio=0.0, theta=0.0,
                        nsigma=4., bilinear+, boundary="nearest", constant=0.0)
                }
                # autoidentify, reidentify, then fitcoords to establish wav.
                # calibration
                autoidentify (imgext, l_crval, l_cdelt, "yes",
                    coordlist=l_coordlist, units=l_units,
                    interactive=l_fl_inter, aidpars="", section=l_section,
                    nsum=l_nsum, ftype=l_ftype, fwidth=l_fwidth,
                    cradius=l_cradius, threshold=l_threshold,
                    minsep=l_minsep, match=l_match, function=l_function,
                    order=l_order, sample=l_sample, niterate=l_niterate,
                    low_reject=l_low_reject, high_reject=l_high_reject,
                    grow=l_grow, dbwrite=l_fl_dbwrite, \
                    overwrite=l_fl_overwrite,\
                    database=l_database, verbose=l_verbose, logfile=l_logfile,
                    plotfile="", graphics="stdgraph", cursor="",
                    aidpars.debug=l_aiddebug, aidpars.crpix=l_crpix,
                    aidpars.cddir="sign", aidpars.ntarget=l_ntarget,
                    aidpars.npattern=l_npattern)

                # Check the file was created
                if (!access(l_database//"/id"//imgext)) {
                    printlog ("ERROR - GSWAVELENGTH: AUTOIDENTIFY did not \
                        find a solution -"//\
                        "\n                      a database file was not \
                        created for "//\
                        "\n                      "//\
                        img//"[SCI,"//n_i//"]", l_logfile, verbose+)

                    imdelete (imgext, verify-, >& "dev$null")
                    goto error
                }

                if (compare) {
                    printlog ("--------- COMPARE 2 ---------", \
                        l_logfile, verbose+)

                    # Use the newly created wavelength database to calculate
                    # new dispersion

                    # Make a tmp file to write the fits output of dispcor too
                    # it is not needed.
                    tmpdispcor = mktemp("tmpdispcor")//".fits"
                    tmp2dispcor = mktemp("tmp2dispcor")

                    dispcor (input=imgext, output=tmpdispcor, \
                        database=l_database, linear+, listonly+, verbose+, \
                        >> tmp2dispcor)

                    # Read the output from dispcor
                    fields (files=tmp2dispcor, fields="13", lines="2",\
                        quit_if_miss=no, print_file_n=no) | \
                        scan (cd1_1_ident)

                    # Value has a comma at the end so strip it off
                    cd1_1_ident = substr(cd1_1_ident, 1, \
                        (strlen(cd1_1_ident) - 1))

                    printlog ("Determined linearised wavelength solution "//\
                        "CD1_1 is: "//(-1.0 * real(cd1_1_ident)), \
                        l_logfile, verbose+)

                    delete (tmpdispcor//", "//tmp2dispcor, verify-, \
                        >& "dev$null")

                    printlog ("--------- COMPARE 2 ---------", \
                        l_logfile, verbose+)
                }

                reidentify (imgext, imgext, coordlist=l_coordlist,
                    interactive=l_fl_inter, section=l_section, newaps=no,
                    refit=l_refit, trace=l_trace, step=l_step, nsum=l_nsum,
                    shift="0.", search=0., nlost=l_nlost, cradius=l_cradius,
                    threshold=l_threshold, addfeatures=l_fl_addfeat,
                    match=l_match, maxfeatures=l_maxfeat, minsep=l_minsep,
                    override=l_fl_overwrite, database=l_database,
                    verbose=l_verbose, logfile=l_logfile, plotfile="",
                    graphics="stdgraph", cursor="")

                # Check the file was created - this may be redundant.
                if (!access(l_database//"/id"//imgext)) {
                    printlog ("ERROR - GSWAVELENGTH: REIDENTIFY did not \
                        find a solution -"//\
                        "\n                      a database file was not \
                        created for "//\
                        "\n                      "//\
                        img//"[SCI,"//n_i//"]", l_logfile, verbose+)

                    imdelete (imgext, verify-, >& "dev$null")
                    goto error
                }

                if (obsmode != "IFU") {

                    if (l_verbose) {
                        fitcoord_log = "STDOUT,"//l_logfile
                    } else {
                        fitcoord_log = l_logfile
                    }

                    fitcoords (imgext, fitname="", interactive=flag_int,
                        combine=no, database=l_database, deletions="",
                        function=l_fitcfunc, xorder=l_fitcxord,
                        yorder=l_fitcyord, logfiles=fitcoord_log,\
                        plotfile="", graphics="stdgraph", cursor="")
                }

                if (compare) {
                    # Compare statement to test wavelength solution against
                    # gsappwave solution
                    printlog ("--------- COMPARE 3 ---------", \
                        l_logfile, verbose+)

                    # Use the newly created fitcoord wavelength database to
                    # evaluate the wavelength of the pixel coordinates in
                    # tmpcoord

                    if (obsmode != "IFU") {

                        printlog ("Wavelengths evaluated by FCEVAL for:", \
                            l_logfile, verbose+)
                        printlog ("        LHS of Central line\n"//\
                            "    RHS of central line\n"//\
                            "    CRPIX1 from GSAPPWAVE of central line"//\
                            "\n"//\
                            "    Bottom left corner\n"//\
                            "    Bottom right corner\n"//\
                            "    Top left corner\n"//\
                            "    Top right corner\n",\
                            l_logfile, verbose+)
                        printlog ("X Y Wavelength[Ang] Y", \
                            l_logfile, verbose+)
                        fceval (input=tmpcoord_compare, output="STDOUT", \
                            fitnames=imgext, database=l_database)
                        fceval (input=tmpcoord_compare, output="STDOUT", \
                            fitnames=imgext, database=l_database, >> l_logfile)
                    }

                    printlog ("--------- COMPARE 3 ---------", \
                        l_logfile, verbose+)
                    delete (tmpcoord_compare, verify-, \
                        >& "dev$null")
                }
                # clean
                imdelete (imgext, verify-, >& "dev$null")
            }
        } # end for-loop over extensions

        # final header update
        gemdate ()
        gemhedit (img//"[0]", "GSWAVELE", gemdate.outdate,
            "UT Time stamp for GSWAVELENGTH")
        gemhedit (img//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last update with GEMINI")
    } # end loop over images
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
    if (status == 0)
        printlog ("GSWAVELENGTH exit status: good", l_logfile, l_verbose)
    else
        printlog ("GSWAVELENGTH exit status: error", l_logfile, l_verbose)
    printlog ("--------------------------------------------------------------\
        ------------------", l_logfile, l_verbose )

end
