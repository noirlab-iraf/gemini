# Copyright(c) 2003-2015 Association of Universities for Research in Astronomy, Inc.

procedure girmfringe(inimages,fringe)

# Scale and subtract a fringe frame from a GMOS gireduced image
#
# Version Jul 30, 2003  BM  created
#         Aug 27, 2003  KL  IRAF2.12 - new parameters
#                             imstat: nclip, lsigma, hsigma, cache
#         Dec 08, 2003 PLG  Changed scale = 1.0 as default
#         Feb 22, 2004 BM   Changed default scale back to 0.0, needed or OLDP
#                           and exposure time scaling is needed in general

string  inimages    {prompt="Input MEF images"}                      # OLDP-1-input-primary-single-prefix=f
string  outimages   {"",prompt="Output MEF images"}                     # OLDP-1-output
string  outpref     {"f",prompt="Prefix for output images"}             # OLDP-4
string  fringe      {prompt="Fringe frame"}                          # OLDP-1-input
bool    fl_statscale    {no,prompt="Scale by statistics rather than exposure time?"} # OLDP-2
string  statsec     {"default",prompt="Extension and region for satistics."} # OLDP-2
real    scale       {0.0,prompt="Override auto-scaling if not 0.0"}     # OLDP-2
bool    fl_propfvardq    {no,prompt="Propagate fringe variance and data quality information?"}
string  logfile     {"",prompt="Log file"}                              # OLDP-1
bool    verbose     {yes,prompt="Verbose output?"}                      # OLDP-4
int     status      {0,prompt="Exit status (0=good)"}                   # OLDP-4
struct  *scanfile   {"",prompt="Internal use only"}                     # OLDP-4

begin

    string  l_inimages, l_outimages, l_fringe, l_statsec
    string  l_logfile, l_outpref
    real    l_scale
    int     l_status
    bool    l_verbose, l_fl_statsc, l_fl_propfvardq

    string  filelist, img, suf, infile[500], outfile[500]
    real    expfringe, expfile, ratio, med, stda, stdb
    int     ii, kk, nbad, nsciext, nf, len, idx
    struct  sdate

    char    pixels_to_use, tmpgtile, stat_file, l_sci_ext, l_var_ext, l_dq_ext

    # Query parameters
    l_inimages=inimages ; l_outimages=outimages ; l_outpref=outpref
    l_fringe=fringe ; l_logfile=logfile ; l_statsec=statsec
    l_scale=scale
    l_verbose=verbose ; l_fl_statsc=fl_statscale
    l_fl_propfvardq = fl_propfvardq

    # Define science extention names
    l_sci_ext = "SCI"
    l_var_ext = "VAR"
    l_dq_ext = "DQ"

    status = 0

    # Keep imgets parameters from changing by outside world
    cache ("imgets","gemlogname", "gtile")

    # Define the name of the logfile
    gemlogname (logpar=l_logfile, package="gmos")
    l_logfile = gemlogname.logname

    date | scan(sdate)
    printlog ("---------------------------------------------------------------\
        -----------------", l_logfile, verbose=l_verbose)
    printlog ("GIRMFRINGE -- "//sdate, l_logfile, verbose=l_verbose)
    printlog (" ", l_logfile, verbose=l_verbose)
    printlog ("inimages = "//l_inimages, l_logfile, verbose=l_verbose)
    printlog ("outimages = "//l_outimages, l_logfile, verbose=l_verbose)
    printlog ("outpref = "//l_outpref, l_logfile, verbose=l_verbose)
    printlog ("fringe = "//l_fringe, l_logfile, verbose=l_verbose)
    printlog ("fl_statscale = "//l_fl_statsc, l_logfile, verbose=l_verbose)
    printlog ("statsec = "//l_statsec, l_logfile, verbose=l_verbose)
    printlog ("scale = "//l_scale, l_logfile, verbose=l_verbose)
    printlog (" ", l_logfile, verbose=l_verbose)

    # Define temporary files
    filelist = mktemp("tmpfilelist")

    #check that a fringe frame is given
    if (l_fringe == "" || l_fringe == " "){
        printlog ("ERROR - GIRMFRINGE: fringe frame not specified",
            l_logfile, verbose+)
        goto error
    }

    # Check that the fringe frame exists
    if (!imaccess(l_fringe)) {
        printlog ("ERROR - GIRMFRINGE: Fringe frame "//l_fringe//" does \
            not exist.", l_logfile, verbose+)
        goto error
    }
    imgets (l_fringe//"[0]","NSCIEXT", >& "dev$null")
    nsciext = int(imgets.value)

    #check that there are input files
    if (l_inimages == "" || l_inimages == " ") {
        printlog ("ERROR - GIRMFRINGE: inimages files not specified",
            l_logfile, verbose+)
        goto error
    }

    #check statistics section
    if (l_fl_statsc && l_scale != 0.0 && (l_statsec=="" || l_statsec==" ")) {
        printlog ("ERROR - GIRMFRINGE: statistics section not given",
            l_logfile, verbose+)
        goto error
    }

    # check existence of input list file
    idx = stridx("@",l_inimages)
    if (idx != 0) {
        len = strlen(l_inimages)
        if (!access(substr(l_inimages,(idx+1),len))) {
            printlog ("ERROR - GIRMFRINGE: "//substr(l_inimages,(idx+1),len)//\
                " does not exist.", l_logfile, verbose+)
            goto error
        }
    }
    files (l_inimages,sort-, > filelist)

    # Check that all images in the list exist and are MEF
    nbad = 0
    kk = 0
    scanfile = filelist
    while (fscan(scanfile, img) != EOF) {
        gimverify (img)
        if (gimverify.status>0) {
            nbad=nbad+1
        } else {
            kk = kk+1
            if (kk > 500) {
                printlog ("ERROR - GIRMFRINGE: Number of images exceeds array \
                    limit.", l_logfile, l_verbose)
                goto error
            }
            # name w/o suffix
            infile[kk] = gimverify.outname//".fits"
            imgets (infile[kk]//"[0]","NSCIEXT", >& "dev$null")
            if (nsciext != int(imgets.value)) {
                printlog ("ERROR - GIRMFRINGE: "//infile[kk]//" has a \
                    different number of extensions than the fringe frame.",
                    l_logfile, verbose+)
                goto error
            }
        }
    }
    nf = kk
    scanfile = ""
    delete (filelist, verify-, >& "dev$null")

    # Exit if problems found
    if (nbad > 0) {
        printlog ("ERROR - GIRMFRINGE: "//nbad//" image(s) either do not \
            exist or are not MEF files.", l_logfile, verbose+)
        goto error
    }

    # Check or create output file names
    if (l_outimages != "" && l_outimages != " "){
        # list file?
        idx = stridx("@",l_outimages)
        if (idx != 0) {
            len = strlen(l_outimages)
            if (!access(substr(l_outimages,(idx+1),len))) {
                printlog ("ERROR - GIRMFRINGE: "//\
                    substr(l_outimages,(idx+1),len)//" does not exist.",
                    l_logfile, verbose+)
                goto error
            }
        }
        files (l_outimages,sort-, > filelist)
        scanfile = filelist
        nbad = 0
        kk = 0
        while(fscan(scanfile, img) != EOF) {
            if (imaccess(img)) {
                nbad = nbad+1
            }
            kk = kk+1
            if (kk > 500) {
                printlog ("ERROR - GIRMFRINGE: Number of output images \
                    exceeds array limit.", l_logfile, l_verbose)
                goto error
            }
            #put suffix on output image
            len = strlen(img)
            suf = substr(img,(len-4),len)
            if (suf != ".fits") {
                outfile[kk] = img//".fits"
            } else {
                outfile[kk] = img
            }
        }
        scanfile = ""
        delete (filelist,verify-, >& "dev$null")
        if (nbad > 0) {
            printlog ("ERROR - GIRMFRINGE: "//nbad//" output images already \
                exist", l_logfile, verbose+)
            goto error
        }
        if (kk != nf) {
            printlog ("ERROR - GIRMFRINGE: number of input and output images \
                not the same.", l_logfile, verbose+)
            goto error
        }
    } else {
        if (l_outpref=="" || l_outpref==" ") {
            printlog ("ERROR - GIRMFRINGE: output prefix not given",
                l_logfile, verbose+)
            goto error
        }
        for (ii=1; ii<=nf; ii+=1) {
            outfile[ii] = l_outpref//infile[ii]
        }
    }


    # Main loop
    if (l_statsec == "default") {
        # Obtain dettype and set
#        hselect (l_fringe//"[0]","DETTYPE","yes") | scan(dettype)
#        if (dettype == "SDSU II CCD") { # Current EEV CCDs
#             l_statsec = "[SCI,2][*,*]"
#        }else if (dettype == "SDSU II e2v DD CCD42-90") {
#             New ev2DD CCDs
#             l_statsec = "[SCI,2][*,*]"
#        }else if (dettype == "S10892") { # Hamamatsu CCDs
#             l_statsec = "[SCI,6][*,*]"
#        }
        pixels_to_use = "[0][*,*]"
    }


    # Get exposure time for fringe frame
    hselect (l_fringe//"[0]","EXPTIME","yes") | scan(expfringe)
    if (l_scale == 0.0 && l_fl_statsc) {
        # Put a call into gtile if l_statsec == "default"
        if (l_statsec == "default") {
            tmpgtile = mktemp("tmpgtile")

            # Tile CCD2 extensions
            gtile (inimages=l_fringe, outimages=tmpgtile, out_ccds="2", \
                fl_stats_only=yes, fl_tile_det=no, sci_ext=l_sci_ext, \
                var_ext="VAR", dq_ext="DQ", mdf_ext="MDF", \
                key_detsec="DETSEC", key_ccdsec="CCDSEC", \
                key_datasec="DATASEC", key_biassec="BIASSEC", \
                key_ccdsum="CCDSUM", fl_verbose=l_verbose, logfile=l_logfile)

            stat_file = tmpgtile//pixels_to_use
        } else {
            stat_file = l_fringe//l_statsec
        }

        imstat(stat_file, fields="midpt,stddev", lower=INDEF,
            upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, binwidth=0.1,
            format-, cache-) | scan (med, stda)
        # Delete the tmpgtile file
        if (l_statsec == "default") {
            delete (tmpgtile//"*.fits", verify-, >& "dev$null")
        }

    }
    ratio = l_scale
    for (ii=1; ii<=nf; ii+=1) {
        #Exposure time
        if (l_scale == 0.0) {
            if (l_fl_statsc) {
                # Put a call into gtile if l_statsec == "default"
                if (l_statsec == "default") {
                    tmpgtile = mktemp("tmpgtile")

                    # Tile CCD2 extensions
                    gtile (inimages=infile[ii], outimages=tmpgtile, \
                        out_ccds="2", \
                        fl_stats_only=yes, fl_tile_det=no, sci_ext=l_sci_ext, \
                        var_ext="VAR", dq_ext="DQ", mdf_ext="MDF", \
                        key_detsec="DETSEC", key_ccdsec="CCDSEC", \
                        key_datasec="DATASEC", key_biassec="BIASSEC", \
                        key_ccdsum="CCDSUM", fl_verbose=l_verbose, \
                        logfile=l_logfile)

                    # Set file to run imstats on
                    stat_file = tmpgtile//pixels_to_use
                } else {
                    # Set file to run imstats on
                    stat_file = infile[ii]//l_statsec
                }

                imstat (stat_file, fields="midpt,stddev",
                    lower=0.01, upper=INDEF, nclip=0, lsigma=INDEF,
                    usigma=INDEF, binwidth=0.1, format-, cache-) | \
                    scan (med, stdb)
                imstat (stat_file, fields="midpt,stddev",
                    lower=(med-3.*stdb), upper=(med+2.5*stdb), nclip=0,
                    lsigma=INDEF, usigma=INDEF, binwidth=0.1, format-,
                    cache-) | scan (med, stdb)

                # Delete the tmpgtile file
                if (l_statsec == "default") {
                    delete (tmpgtile//"*.fits", verify-, >& "dev$null")
                }

                ratio = stdb/stda
            } else {
                hselect (infile[ii]//"[0]", "EXPTIME", "yes") | scan (expfile)
                ratio = expfile/expfringe
            }
        }
        printlog (infile[ii]//" "//outfile[ii]//" ratio = "//ratio,
            l_logfile, l_verbose)

        # TODO This is a hack - this whole script requires a re-write
        # TODO Handle VAR and DQ better
        if ((imaccess(l_fringe//"["//l_var_ext//",1]") && \
                imaccess(l_fringe//"["//l_dq_ext//",1]")) && \
                l_fl_propfvardq) {
            # Use gemexpr to update the VAR and DQ planes too!
            # Don't use 'a=', 'b=' etc. CL doesn't like it.
            # Use int(a[DQ]) as input DQ planes not always integers
            gemexpr ("a-(b*c)", outfile[ii], infile[ii], l_fringe, ratio, \
                var_expr="a[VAR] - (b[VAR] * c * c)", \
                dq_expr="int(a[DQ]) | int(b[DQ])", fl_vardq=l_fl_propfvardq, \
                verbose=l_verbose)
        } else {
            copy (infile[ii], outfile[ii], verbose-)
            for (kk = 1; kk <= nsciext; kk += 1) {
                imexpr ("a - (b * c)", \
                    outfile[ii]//"["//l_sci_ext//","//kk//",overwrite+]", \
                    infile[ii]//"["//l_sci_ext//","//kk//"]", \
                    l_fringe//"["//l_sci_ext//","//kk//"]", \
                    ratio, verbose-)
            }
        }

        gemdate ()
        gemhedit (outfile[ii]//"[0]", "GIRMFRIN", gemdate.outdate,
            "UT Time stamp form GIRMFRINGE", delete-)
        gemhedit (outfile[ii]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
    }
    goto clean

error:
    status = 1
    goto clean

clean:
    delete (filelist, verify-, >& "dev$null")
    # close log file
    printlog (" ", l_logfile, l_verbose)
    printlog ("GIRMFRINGE done", l_logfile, l_verbose)
    printlog ("-----------------------------------------------------------\
        ---------------------", l_logfile, l_verbose )

end
