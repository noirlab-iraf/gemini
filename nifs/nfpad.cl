# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

# NFPAD - Subtract off reference pixels and pad NIFS data arrays
#
# Original author: Tracy Beck

procedure nfpad (inimage, outimage) 

# The following steps are taken:
# - NIFS data files are padded by the appropriate pixels to take into account 
#   the MDF spectral overfilling
# - NIFS data are subtracted for the reference pixels.
#
# NOTE: "nfpad" is written to be a simple script only to be called from other 
#       tasks (i.e., nfprepare and nfacquire).  It does not have the full 
#       functionality of some of the other IRAF scripts because it doesn't 
#       use "@" lists or comma separated input images.

char    inimage     {prompt = "Input NIFS image"}        
char    outimage    {prompt = "Output NIFS image"}
char    exttype     {"SCI", prompt = "SCI, VAR or DQ?"}
char    logfile     {"", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose?"}
int     status      {0, prompt = "Exit status (0 = good)"}

begin

    char    l_inimage = ""
    char    l_outimage = ""
    char    l_logfile = ""
    char    l_exttype = ""
    bool    l_verbose
    
    char    tmpsci3, tmpvar3, tmpdq3, tmpinimage3
    int     mefvalue, numrow, i 
    int     n_ext, nx, ny

    cache ("gemextn", "gemdate")

    status = 0
    
    l_inimage = inimage
    l_outimage = outimage
    l_logfile = logfile
    l_exttype = exttype
    l_verbose = verbose

    if (l_logfile == "") {
        l_logfile = "nifs.log"
    }

    printlog ("--------------------------------------------------------------\
        --------------", l_logfile, verbose=l_verbose)    
    printlog ("NFPAD: Padding NIFS data in the y-dimension for extended slices",
        l_logfile, verbose=l_verbose)

    n_ext = 0
    mefvalue = 0

    gemextn (l_inimage, check="exists", process="none", index="", extname="",
        extversion="", ikparams="", omit="", replace="", outfile="STDOUT",
        logfile="", verbose-, >& "dev$null")
    if (gemextn.fail_count == 1) {
        printlog ("ERROR: NFPAD - Input file does not exist", l_logfile,
            verbose+)
        status = 1
        goto clean
    }
    gemextn (l_inimage, check="mef", process="expand", index="", extname="",
        extversion="", ikparams="", omit="", replace="", outfile="dev$null",
        logfile="", verbose-, >& "dev$null")
    if (gemextn.fail_count != 0) {
        hselect (l_inimage, "i_naxis1", yes) | scan (nx)
        hselect (l_inimage, "i_naxis2", yes) | scan (ny)
        if ((nx == 2048) && (ny == 2048)) {
            printlog ("NFPAD: Input file is a non-MEF NIFS frame", 
                l_logfile, verbose=l_verbose)
            mefvalue=0
        } else {
            printlog ("ERROR: NFPAD Exiting, MEF file is not in a recognizable \
                format", l_logfile, verbose+)
            status = 1
            goto clean
        }
    } else {
        printlog ("NFPAD: Input file is a MEF", l_logfile, ver+)
        mefvalue = 1

        if (gemextn.count == 2) {
            printlog ("NFPAD: Input file is a MEF with a PHU and one image \
                extension", l_logfile, verbose+)
            n_ext = 1
        } else if (gemextn.count == 4) {
            printlog ("NFPAD: Input file is a MEF with a PHU and three image \
                extensions", l_logfile, verbose+)
            n_ext = 3
        } else if (gemextn.count == 5) {
            printlog ("NFPAD: Input file is a MEF with a PHU, MDF, and three \
                image extensions", l_logfile, verbose+)
            n_ext = 4
        }
    }
    
    gemextn (l_outimage, check="exists", process="none", index="", extname="",
        extversion="", ikparams="", omit="", replace="", outfile="STDOUT",
        logfile="", verbose-, >& "dev$null")
    if (gemextn.fail_count == 0) {
        printlog ("ERROR: NFPAD  Output file already exists",
            l_logfile, verbose+)
        status = 1
        goto clean
    }


    #--------------------------------------------------------------------------
    # start output

    # Create tmp FITS file names used within this loop

    tmpsci3 = mktemp("tmpsci3")
    tmpvar3 = mktemp("tmpvar3")
    tmpdq3 = mktemp("tmpdq3")
    tmpinimage3 = mktemp("tmpinimage3")

    #make padding array

    mkpattern (tmpinimage3, output="", pattern="constant", option="replace",
        v1=0., v2=1., size=1, title="", pixtype="real", ndim=2, ncols=2040, 
        nlines=2080, header="")
    imcopy (tmpinimage3, tmpsci3, ver-)
    imcopy (tmpinimage3, tmpvar3, ver-)
    imcopy (tmpinimage3, tmpdq3, ver-)

    #copy NIFS image extensions into padding array


    if (l_exttype == "") {
    numrow=25
        if (mefvalue == 0) {
            imcopy (l_inimage//"[5:2044,5:2044]",
                tmpinimage3//"[1:2040,26:2065]", ver-)
            for (i=1; i<=numrow; i+=1) {
               imcopy (l_inimage//"[5:2044,5:5]",
                   tmpinimage3//"[1:2040,"//i//":"//i//"]", ver-)
            }
        }
    }
   
    if (l_exttype == "SCI") {
        if (mefvalue == 0) {
            imcopy (l_inimage//"[5:2044,5:2044]",
                tmpinimage3//"[1:2040,26:2065]", ver-)
    numrow=25
            for (i=1; i<=numrow; i+=1) {
               imcopy (l_inimage//"[5:2044,5:5]",
                   tmpinimage3//"[1:2040,"//i//":"//i//"]", ver-)
            }
        }

        if ((mefvalue == 1) && (n_ext == 1)) {
            imcopy (l_inimage//"[1][5:2044,5:2044]",
                tmpinimage3//"[1:2040,26:2065]", ver-)
            for (i=1; i<=numrow; i+=1) {
               imcopy (l_inimage//"[5:2044,5:5]",
                   tmpinimage3//"[1:2040,"//i//":"//i//"]", ver-)
            }
        } 
        if ((mefvalue == 1) && (n_ext == 3)) {
            imcopy (l_inimage//"[1][5:2044,5:2044]",
                tmpsci3//"[1:2040,26:2065]", ver-)
            imcopy (l_inimage//"[2][5:2044,5:2044]",
                tmpvar3//"[1:2040,26:2065]", ver-)
            imcopy (l_inimage//"[3][5:2044,5:2044]",
                tmpdq3//"[1:2040,26:2065]", ver-)
            imreplace (tmpdq3//"[1:2040,1:20]", value=1.0, imaginary=0.,
                lower=INDEF, upper=INDEF, radius=0.)
            imreplace (tmpdq3//"[1:2040,2061:2080]", value=1.0, imaginary=0.,
                lower=INDEF, upper=INDEF, radius=0.)
        } 
        if ((mefvalue == 1) && (n_ext == 4)) {
            imcopy (l_inimage//"[2][5:2044,5:2044]",
                tmpinimage3//"[1:2040,26:2065]", ver-)
            imcopy (l_inimage//"[3][5:2044,5:2044]",
                tmpvar3//"[1:2040,26:2065]", ver-)
            imcopy (l_inimage//"[4][5:2044,5:2044]",
                tmpdq3//"[1:2040,26:2065]", ver-)
            imreplace (tmpdq3//"[1:2040,1:20]", value=1.0, imaginary=0., 
                lower=INDEF, upper=INDEF, radius=0.)
            imreplace (tmpdq3//"[1:2040,2061:2080]", value=1.0, imaginary=0.,
                lower=INDEF, upper=INDEF, radius=0.)
        }
    }

    if (l_exttype == "VAR") {
        if (mefvalue == 0) {
            imcopy (l_inimage//"[5:2044,5:2044]",
                tmpinimage3//"[1:2040,26:2065]", ver-)
        } else {
            printlog ("ERROR: NFPAD - Unsupported file structure for VARIANCE \
                plane", l_logfile, verbose+)
            status = 1
            goto clean
        }
    }

    if (l_exttype == "DQ") {
        if (mefvalue == 0) {
            imcopy (l_inimage//"[5:2044,5:2044]",
                tmpinimage3//"[1:2040,26:2065]", ver-)
            imreplace (tmpinimage3//"[1:2040,1:20]", value=1.0, imaginary=0.,
                lower=INDEF, upper=INDEF, radius=0.)
            imreplace (tmpinimage3//"[1:2040,2061:2080]", value=1.0, 
                imaginary=0., lower=INDEF, upper=INDEF, radius=0.)
        } else {
            printlog ("ERROR: NFPAD - Unsupported file structure for DQ plane",
                l_logfile, verbose+)
            status = 1
            goto clean
        }
    }

    #write output file

    if (mefvalue == 0){
        imcopy (tmpinimage3, l_outimage, verbose-)
    }
    if ((mefvalue == 1) && (n_ext == 1)) {
        wmef (tmpinimage3, output=l_outimage, extname="SCI,VAR,DQ", 
            phu=l_inimage, verbose-, >& "dev$null")
        if (wmef.status != 0) {
            printlog ("ERROR: NFPAD - Could not write final MEF file for \
                unknown reason.", l_logfile, verbose+)
            status = 1
            goto clean
        }
    }
    if ((mefvalue == 1) && (n_ext != 1)) {
        wmef (tmpsci3//","//tmpvar3//","//tmpdq3, output=l_outimage, 
            extnames="SCI,VAR,DQ", phu=l_inimage, verbose-, >& "dev$null")
        if (wmef.status != 0) {
            printlog ("ERROR: NFPAD - Could not write final MEF file for \
                unknown reason.", l_logfile, verbose+)
            status = 1
            goto clean
        }           
    }
    
    gemdate()
    gemhedit (l_outimage//"[0]", "NFPAD", gemdate.outdate, 
        "UT Time stamp for NFPAD", delete-)
    gemhedit (l_outimage//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)

clean:
    imdelete (tmpinimage3, verify-, >& "dev$null")
    imdelete (tmpsci3, verify-, >& "dev$null")
    imdelete (tmpvar3, verify-, >& "dev$null")
    imdelete (tmpdq3, verify-, >& "dev$null")

    #-------------------------------------------------------------------------
    # Clean up
    if (status==0) {
        printlog ("NFPAD exit status:  good.", l_logfile, verbose=l_verbose)
    } else {
        printlog ("NFPAD exited with errors.", l_logfile, verbose=l_verbose)
    }
    printlog ("--------------------------------------------------------------\
        --------------", l_logfile, verbose=l_verbose)

end

