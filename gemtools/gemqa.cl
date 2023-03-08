# Copyright(c) 2000-2011 Association of Universities for Research in Astronomy, Inc.

procedure gemqa (inimages, finalint, finaliq, piconfig, pitarget)

# Add data quality assessemnt keywords to reduced images
# 
# Version  Oct 12, 2001 IJ, v1.2 release
#          Sept 20,2002 IJ, v1.4 release

char    inimages    {prompt="Images to be updated"}  # OLDP-1-input-primary-update
char    finalint    {min="pass|fail|unknown",prompt="Final integration time (pass|fail)"}  # OLDP-1
char    finaliq     {min="pass|fail|unknown",prompt="Final image quality (pass|fail)"}  # OLDP-1
char    piconfig    {min="pass|fail|unknown",prompt="Used correct instrument configuration (pass|fail)"}  # OLDP-1
char    pitarget    {min="pass|fail|unknown",prompt="Observed correct target (pass|fail)"}  # OLDP-1
char    sci_ext     {"0",prompt="If MEF, name or number of extension to update"}  # OLDP-3
char    logfile     {"gemqa.log",prompt="Name of log file"}  # OLDP-1
bool    verbose     {no,prompt="Verbose actions"}  # OLDP-4
int     status      {0,prompt="Exit status (0=good)"}  # OLDP-4
struct  *scanfile   {"",prompt="Internal use only"}  # OLDP-4

begin

    # declare local variables
    char    l_inimages, l_finalint, l_finaliq, l_piconfig, l_pitarget
    char    l_sci_ext, l_logfile
    char    tmpimglist, tmpextnlist, tmpfile, img, imgextn
    char    projID, newprojID, prefix, suffix, datalabel, obsID
    char    keyfound
    bool    l_verbose
    struct  l_struct

    status = 0
    scanfile = ""
    cache ("imgets", "gemextn", "gemdate")

    # Set local variables to input parameters
    l_inimages=inimages 
    l_finalint=finalint ; l_finaliq=finaliq ; l_piconfig=piconfig
    l_pitarget=pitarget
    l_sci_ext=sci_ext ; l_logfile=logfile
    l_verbose=verbose

    # Check for package log file or user-defined log file
    if ((l_logfile=="") || (l_logfile==" ")) {
        l_logfile = "gemqa.log"
        printlog ("WARNING - GEMQA: No logfile defined. Using default file \
            gemqa.log", l_logfile, verbose+)
    }

    # Open log file
    date | scan(l_struct)
    printlog ("--------------------------------------------------------------\
        --------------", l_logfile, l_verbose)
    printlog ("GEMQA -- "//l_struct, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    
    # Make temporary files
    tmpimglist = mktemp ("tmpimglist")
    tmpextnlist = mktemp ("tmpextnlist")
    tmpfile = mktemp ("tmpfile")
    
    # Get input file names, make sure the images and the extensions
    # do exist.
    gemextn (l_inimages, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension", replace="",
        outfile=tmpimglist, logfile=l_logfile, glogpars="", verbose=l_verbose)
    sections ("@"//tmpimglist//"//\["//l_sci_ext//"]", > tmpfile)
    delete (tmpimglist, verify-, >& "dev$null")
    gemextn ("@"//tmpfile, check="exist,mef", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpextnlist, logfile=l_logfile, glogpars="", verbose=l_verbose)
    delete (tmpfile, verify-, >& "dev$null")
    
    if (gemextn.fail_count > 0) {
        printlog ("ERROR -- GEMQA: "//gemextn.fail_count//" images/extensions \
            were not found.", l_logfile, verbose+)
        status = 1
        goto clean
    } else if (gemextn.count == 0) {
        printlog ("ERROR -- GEMQA: No input images defined.", l_logfile,
            verbose+)
        status = 1
        goto clean
    }
    
    # now that we've check that the extensions exist, prune the list
    # of the extension designation.  This is necessary to allow the 
    # update of GEM-TLM for example which should always be in the PHU
    # (extension 0) regarless of the value of 'l_sci_ext'
    gemextn ("@"//tmpextnlist, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="index,kernel", replace="",
        outfile=tmpimglist, logfile=l_logfile, glogpars="", verbose=l_verbose)
    if (gemextn.fail_count > 0) { # this should never happen
        printlog ("ERROR -- GEMQA: Internal error.", l_logfile, verbose+)
        status = 1
        goto clean
    }

    printlog ("Updating images:", l_logfile, l_verbose)
    if (l_verbose) type(tmpimglist)
    type(tmpimglist, >> l_logfile)

    printlog ("FINALINT = "//l_finalint, l_logfile, l_verbose)
    printlog ("FINALIQ  = "//l_finaliq, l_logfile, l_verbose)
    printlog ("PICONFIG = "//l_piconfig, l_logfile, l_verbose)
    printlog ("PITARGET = "//l_pitarget, l_logfile, l_verbose)

    #-------------------------------------------------------------------------
    # Update the headers
    scanfile = tmpimglist
    while (fscan(scanfile,img)!=EOF) {
        imgextn = img//"["//l_sci_ext//"]"

        # Make DATALAB a unique identifier

        # Parse the image name. get prefix and suffix
        hselect (imgextn, "DATALAB", yes) | scan(datalabel)
        prefix = ""
        suffix = ""
        if ( strstr("N",datalabel) != 0 ) {
            if ( stridx("N",img) > 1 )
                prefix = "-"//substr(img,1,stridx("N",img)-1)
        } else if ( strstr("S",datalabel) != 0 ) {
            if ( stridx("S",img) > 1 )
                prefix = "-"//substr(img,1,stridx("S",img)-1)
        } else {
            printlog ("ERROR - GEMQA: project ID, "//datalabel//", is invalid",
                l_logfile, verbose+)
            status = 1
            goto clean
        }
        if ( stridx("_",img) != 0 ) {
            suffix = "-"//substr(img, stridx("_",img)+1, strlen(img) )
            print(suffix) | translit("STDIN", "_", "-", delete-, collapse-) | \
                scan(suffix)
        }

        # Parse the datalabel to extract the 'obsID-XXX' string; rebuild string
        hselect (imgextn, "OBSID", yes) | scan(obsID)
        projID = substr(datalabel, 1, strlen(obsID)+4)
        datalabel = projID//prefix//suffix
        hedit (imgextn,"DATALAB", datalabel, add+, addonly-, delete-, verify-,
            show-, update+)

        # Edit the other keywords
        gemhedit (imgextn, "FINALINT", l_finalint,
            "Requested total integration time")
        gemhedit (imgextn, "FINALIQ", l_finaliq,
            "Requested final image quality")
        gemhedit (imgextn, "PICONFIG", l_piconfig,
            "Requested instrument configuration")
        gemhedit (imgextn, "PITARGET", l_pitarget,
            "Requested target")
        
        # Edit GEM-TLM if this is a processed image (raw images should
        # not have a GEM-TLM
        
        keyfound = ""
        hselect (img//"[0]", "GEM-TLM", yes) | scan (keyfound)
        if (keyfound != "") { # keyword present => processed image => update
            gemdate ()
            gemhedit (img//"[0]", "GEM-TLM", gemdate.outdate, 
                "UT Last modification with GEMINI", delete-)
        }

    }
    scanfile = ""

    #------------------------------------------------------------------------

    # Clean up
clean:
    delete (tmpextnlist, verify-, >& "dev$null")
    delete (tmpimglist, verify-, >& "dev$null")
    delete (tmpfile, verify-, >& "dev$null")
    
    if (status==0)
        printlog ("GEMQA exit status: good", l_logfile, l_verbose)
    printlog ("--------------------------------------------------------------\
        --------------", l_logfile, l_verbose)

end

