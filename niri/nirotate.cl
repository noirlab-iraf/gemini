# Copyright(c) 2004-2011 Association of Universities for Research in Astronomy, Inc.

procedure nirotate(inimages)

# Simple script to figure out rotations from NIRI WCS info
# and rotate the images.  Primary application is Altair data
# taken with the cass rotator fixed.  
#
# 28 Oct 2004  first working version, J. Jensen
# 15 Nov 2004  WCS angles using images.immatch.skymap, J. Jensen
#  7 Feb 2005  K. Labrie: complete overhaul to make it releaseable 
#                         + add OLDP encoding

char    inimages    {prompt="Input NIRI images"}        # OLDP-1-input-primary-single-prefix=rot
char    outimages   {"", prompt="Output images"}        # OLDP-1-output
char    outprefix   {"rot", prompt="Output prefix"}      # OLDP-4
char    interpolant {"linear", prompt="Interpolant for ROTATE (nearest,linear,poly3,poly5,spline3,sinc,lsinc,drizzle)"}  # OLDP-2
bool    fl_vardq    {yes, prompt="Propagate variance and data quality frames"} # OLDP-3
char    logfile     {"", prompt="Name of the logfile"}  # OLDP-1
pset    glogpars    {"", prompt="Logging preferences"}  # OLDP-4
bool    verbose     {yes, prompt="Verbose output"}      # OLDP-4
char    sci_ext     {"SCI", prompt="Name of science extension"} # OLDP-3
char    var_ext     {"VAR", prompt="Name of variance extension"} # OLDP-3
char    dq_ext      {"DQ", prompt="Name of data quality extension"} # OLDP-3
int     status      {0, prompt="Exit status (0=good)"}  # OLDP-4
struct  *scanfile   {"", prompt="Internal use only"}    # OLDP-4

begin

    char    l_inimages
    char    l_outimages = ""
    char    l_outprefix = ""
    char    l_logfile = ""
    char    l_interp, l_sci_ext, l_var_ext, l_dq_ext
    bool    l_verbose, l_fl_vardq
    
    char    paramstr, errmsg, msgstr, outputstr
    char    tmpfile, tmpinimg, tmpoutimg
    char    tmpimg, tmpvar, tmpdq
    char    tmpfxlog, tmprotlog
    char    in[100], out[100], filename
    char    refsci, crpixstr1, crpixstr2
    char    field, valuestr
    int     naxis1, naxis2
    int     nimages, noutimages, maximages, nmissing, warn
    int     i, junk
    int     xmax, ymax, deltax, deltay
    real    rotation, xcen, ycen
    bool    fl_found

    l_inimages = inimages
    junk = fscan (outimages, l_outimages)
    junk = fscan (outprefix, l_outprefix)
    junk = fscan (interpolant, l_interp)
    l_fl_vardq = fl_vardq
    junk = fscan (logfile, l_logfile)
    junk = fscan (sci_ext, l_sci_ext)
    junk = fscan (var_ext, l_var_ext)
    junk = fscan (dq_ext, l_dq_ext)
    l_verbose = verbose
    
    # Initialize
    status = 0
    warn = 0
    nimages = 0
    noutimages = 0
    nmissing = 0
    maximages = 100
    
    cache ("gemextn", "gemdate")
        
    # Temporary files
    tmpfile = mktemp("tmpfile")
    tmpinimg = mktemp("tmpinimg")
    tmpoutimg = mktemp("tmpoutimg")

    # Create the list of parameter/value pairs. One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outprefix      = "//outprefix.p_value//"\n"
    paramstr += "interpolant    = "//interpolant.p_value//"\n"
    paramstr += "fl_vardq       = "//fl_vardq.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "glogpars       = "//glogpars.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value//"\n"
    paramstr += "sci_ext        = "//sci_ext.p_value//"\n"
    paramstr += "var_ext        = "//var_ext.p_value//"\n"
    paramstr += "dq_ext         = "//dq_ext.p_value
    
    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "nirotate", "niri", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Load up the array of input file names
    gemextn (l_inimages, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension", replace="",
        outfile=tmpfile, logfile=l_logfile, verbose=l_verbose)
    gemextn ("@"//tmpfile, check="exist,mef", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpinimg, logfile=l_logfile, verbose=l_verbose)
    nimages = gemextn.count
    delete (tmpfile, ver-, >& "dev$null")
    
    if ((gemextn.fail_count > 0) || (nimages == 0) || \
        (nimages > maximages)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" images were not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input images defined."
            status = 121
        } else if (nimages > maximages) {
            errmsg = "Maximum number of input images ["//str(maximages)//"] \
                has been exceeded."
            status = 121
        }
        
        glogprint (l_logfile, "nirotate", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto clean
        
    } else {
        scanfile = tmpinimg
        i = 0
        while (fscan(scanfile, filename) != EOF) {
            i += 1
            in[i] = filename
        }
        scanfile = ""
        if (i != nimages) {
            status = 99
            glogprint (l_logfile, "nirotate", "status", type="error",
                errno=status, str="Error while counting the input images.",
                verbose+)
            goto clean
        }
    }
    # If fl_vardq is yes, ensure that all input images do have 
    # VAR and DQ planes
    if (l_fl_vardq) {
        gemextn ("@"//tmpinimg, process="append", extname="VAR", \
            check="ext=exists", index="", extversion="", ikparams="", omit="",
            replace="", outfile="dev$null", logfile=l_logfile,
            verbose=l_verbose)
        nmissing = gemextn.fail_count
        gemextn ("@"//tmpinimg, process="append", extname="DQ", \
            check="ext=exists", index="", extversion="", ikparams="", omit="",
            replace="", outfile="dev$null", logfile=l_logfile, 
            verbose=l_verbose)
        nmissing += gemextn.fail_count
        if (nmissing > 0) {
            warn = 123
            glogprint (l_logfile, "nirotate", "status", type="warning",
                errno=warn, str="Not all images have variance and data \
                quality planes.  Re-setting fl_vardq to 'no'",
                verbose=l_verbose)
            l_fl_vardq = no
        }
    }
    
    # Load up the array of output file names
    if (l_outimages != "")
        outputstr = l_outimages
    else if (l_outprefix != "") {
        gemextn ("@"//tmpinimg, check="", process="none", index="", extname="",
            extversion="", ikparams="", omit="path", replace="",
            outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
        outputstr = l_outprefix//"@"//tmpoutimg
    } else {
        status = 121
        glogprint (l_logfile, "nirotate", "status", type="error", errno=status,
            str="Neither output image name nor output prefix is defined.",
            verbose+)
        goto clean
    }
    
    gemextn (outputstr, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension", replace="",
        outfile=tmpfile, logfile=l_logfile, verbose=l_verbose)
    delete (tmpoutimg, ver-, >& "dev$null")
    gemextn ("@"//tmpfile, check="absent", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
    noutimages = gemextn.count
    delete (tmpfile, ver-, >& "dev$null")
    
    if ((gemextn.fail_count > 0) || (noutimages == 0) || \
        (noutimages != nimages)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" image(s) already exist(s)."
            status = 102
        } else if (noutimages == 0) {
            errmsg = "No output images defined."
            status = 121
        } else if (noutimages != nimages) {
            errmsg = "Different number of input images ("//nimages//") and \
                output images ("//noutimages//")."
            status = 121
        }
        
        glogprint (l_logfile, "nirotate", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto clean
    } else {
        scanfile=tmpoutimg
        i = 0
        while (fscan(scanfile, filename) != EOF) {
            i += 1
            out[i] = filename//".fits"
        }
        scanfile=""
        if (i != noutimages) {
            status = 99
            glogprint (l_logfile, "nirotate", "status", type="error",
                errno=status, str="Error while counting the output images.",
                verbose+)
            goto clean
        }
    }
    delete (tmpinimg, ver-, >& "dev$null")
    delete (tmpoutimg, ver-, >& "dev$null")
    
    # Now do the work
    for (i=1; i <= nimages; i+=1) {

        if (i == 1) {
            refsci = in[i]//"["//l_sci_ext//"]"
            crpixstr1 = "" ; crpixstr2 = ""
            hselect (refsci, "CRPIX1", yes) | scan (crpixstr1)
            hselect (refsci, "CRPIX2", yes) | scan (crpixstr2)
            if ((crpixstr1 == "") || (crpixstr2 == "")) {
                status = 131
                glogprint (l_logfile, "nirotate", "status", type="error",
                    errno=status, str="Missing WCS data in "//refsci, 
                    verbose+)
                goto clean
            }
            xcen = real(crpixstr1)
            ycen = real(crpixstr2)
            
            hselect (refsci, "i_naxis1", yes) | scan (naxis1)
            hselect (refsci, "i_naxis2", yes) | scan (naxis2)
            xmax = int( real(naxis1) * 1.41 )
            ymax = int( real(naxis2) * 1.41 )
            deltax = int( real(xmax - naxis1) / 2. )
            deltay = int( real(ymax - naxis2) / 2. )
            
            rotation = 0.
            
        } else {
            images.immatch.skymap(in[i]//"["//l_sci_ext//"]", refsci,
                database=tmpfile, transforms="", results="", xmin=INDEF,
                ymin=INDEF, xmax=INDEF, ymax=INDEF, nx=10, ny=10, wcs="world",
                fitgeometry="rotate", function="polynomial", interactive=no,
                verbose=no)
            
            # We are interested in the value of 'xrotation' in database
            # 'tmpfile'.  There does not seem to be CL tasks to access
            # databases.  Going down the list until 'xrotation' is found.

            scanfile = tmpfile
            fl_found = no
            while (fscan(scanfile, field, valuestr) != EOF) {
                if (field == "xrotation") {
                    rotation = real(valuestr)
                    fl_found = yes
                    break
                }
            }
            scanfile = ""
            delete (tmpfile, ver-, >& "dev$null")
            if (fl_found == no) {
                status = 99
                glogprint (l_logfile, "nirotate", "status", type="error",
                    errno=status, str="Unable to find rotation angle for \
                    image "//in[i], verbose+)
                goto clean
            }
        }
        
        msgstr = "   "//in[i]//"   "//rotation  
        glogprint (l_logfile, "nirotate", "science", type="string",
            str=msgstr, verbose=l_verbose)

        # ----- rotate the image - SCI extension ----- #
        tmpfxlog = mktemp("tmpfxlog")
        tmprotlog = mktemp("tmprotlog")
        tmpimg = mktemp("tmpimg")
        
        # Copy PHU to output file
        fxcopy (in[i]//"[0]", out[i], >& tmpfxlog)
        glogprint (l_logfile, "nirotate", "task", type="file", str=tmpfxlog,
            verbose=l_verbose)
        delete (tmpfxlog, ver-, >& "dev$null")
        
        # Rotate and insert the SCI extension to the output file
        rotate (in[i]//"["//l_sci_ext//"]", tmpimg, rotation, xin=xcen,
            yin=ycen, xout=xcen+deltax, yout=ycen+deltay, interpolant=l_interp,
            ncols=xmax, nlines=ymax, boundary="nearest", >>& tmprotlog)
        fxinsert (tmpimg//".fits", out[i]//"[1]", groups="", ver-, >>& tmpfxlog)
        glogprint (l_logfile, "nirotate", "engineering", type="file",
            str=tmprotlog, verbose=l_verbose)
        glogprint (l_logfile, "nirotate", "task", type="file", str=tmpfxlog,
            verbose=l_verbose)
        delete (tmpfxlog//","//tmprotlog, ver-, >& "dev$null")
        imdelete (tmpimg, ver-, >& "dev$null")
        
        
        
        # Take care of the VAR/DQ planes
        if (l_fl_vardq) {
            tmpfxlog = mktemp("tmpfxlog")
            tmprotlog = mktemp("tmprotlog")
            tmpvar = mktemp("tmpvar")
            tmpdq = mktemp("tmpdq")
            
            # Rotate and insert the VAR plane to the output file
            rotate (in[i]//"["//l_var_ext//"]", tmpvar, rotation, xin=xcen,
                yin=ycen, xout=xcen+deltax, yout=ycen+deltay,
                interpolant=l_interp, ncols=xmax, nlines=ymax,
                boundary="nearest", >& tmprotlog)
            fxinsert (tmpvar//".fits", out[i]//"[2]", groups="", ver-,
                >>& tmpfxlog)
            glogprint (l_logfile, "nirotate", "engineering", type="file",
                str=tmprotlog, verbose=l_verbose)
            glogprint (l_logfile, "nirotate", "task", type="file", str=tmpfxlog,
                verbose=l_verbose)
            delete (tmpfxlog//","//tmprotlog, ver-, >& "dev$null")
            imdelete (tmpvar, ver-, >& "dev$null")
            
            # Rotate and insert the DQ plane to the output file
            rotate (in[i]//"["//l_dq_ext//"]", tmpdq, rotation, xin=xcen,
                yin=ycen, xout=xcen+deltax, yout=ycen+deltay,
                interpolant="nearest", ncols=xmax, nlines=ymax,
                boundary="constant", constant=1., >& tmprotlog)
            fxinsert (tmpdq//".fits", out[i]//"[3]", groups="", ver-,
                >>& tmpfxlog)
            glogprint (l_logfile, "nirotate", "engineering", type="file",
                str=tmprotlog, verbose=l_verbose)
            glogprint (l_logfile, "nirotate", "task", type="file", str=tmpfxlog,
                verbose=l_verbose)
            delete (tmpfxlog//","//tmprotlog, ver-, >& "dev$null")
            imdelete (tmpdq, ver-, >& "dev$null")
        }
        
        # update PHU
        gemdate ()
        gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate, "", delete-)
        gemhedit (out[i]//"[0]", "ROTATE", gemdate.outdate, 
            "UT Time stamp for NIROTATE", delete-)
        
        # update SCI/VAR/DQ headers
        gemhedit (out[i]//"[1]", "EXTNAME", l_sci_ext, "Extension name", delete-)
        gemhedit (out[i]//"[1]", "EXTVER", 1, "Extension version", delete-)
        gemhedit (out[i]//"[1]", "INHERIT", 'F', "Inherits global header", delete-)
        gemhedit (out[i]//"[1]", "NIROTATE", rotation,
            "Angle of rotation (degrees)", delete-)
        
        if (l_fl_vardq) {
            gemhedit (out[i]//"[2]", "EXTNAME", l_var_ext, 
                "Extension name", delete-)
            gemhedit (out[i]//"[2]", "EXTVER", 1, "Extension version", delete-)
            gemhedit (out[i]//"[2]", "INHERIT", 'F', 
                "Inherits global header", delete-)
            gemhedit (out[i]//"[2]", "NIROTATE", rotation,
                "Angle of rotation (degrees)", delete-)
                
            gemhedit (out[i]//"[3]", "EXTNAME", l_dq_ext, 
                "Extension name", delete-)
            gemhedit (out[i]//"[3]", "EXTVER", 1, "Extension version", delete-)
            gemhedit (out[i]//"[3]", "INHERIT", 'F', 
                "Inherits global header", delete-)
            gemhedit (out[i]//"[3]", "NIROTATE", rotation,
                "Angle of rotation (degrees)", delete-)
        }        
    }  

clean:
    scanfile = ""
    delete (tmpfile, ver-, >& "dev$null")
    delete (tmpinimg, ver-, >& "dev$null")
    delete (tmpoutimg, ver-, >& "dev$null")
    
    if (status == 0) {
        if (warn != 0) {
            glogprint (l_logfile, "nirotate", "status", type="warning",
                errno=warn, str="An important warning was issued.  \
                Please review logs.", verbose+)
        }
        glogclose (l_logfile, "nirotate", fl_success+, verbose=l_verbose)
    } else
        glogclose (l_logfile, "nirotate", fl_success-, verbose=l_verbose)

exitnow:
    ;

end

