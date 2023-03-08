# Copyright(c) 2002-2012 Association of Universities for Research in Astronomy, Inc.

procedure gnsskysub(inimages)

# Generate a sky subtracted nod & shuffled image.
#
# Version Sept 20, 2002  RA,IJ  v1.4 release
#         Feb 25, 2003   MB: generalized the loop over sciext for 6amp mode
#         May 09, 2003   IJ, removed all the old code
#         Jun 03, 2003   IJ  check sci,1 for binning
#         Aug 20, 2003   IJ  fool proof imshift call

char    inimages    {prompt="Input GMOS images or list"}
char    outimages   {"",prompt="Output images or list"}
char    outpref     {"n",prompt="Prefix for output images"}
bool    fl_fixnc    {no,prompt="Auto-correct for nod count mismatch?"}
char    sci_ext     {"SCI", prompt="Name of science extension"}
char    var_ext     {"VAR", prompt="Name of variance extension"}
char    dq_ext      {"DQ", prompt="Name of data quality extension"}
char    mdf_ext     {"MDF", prompt="Mask definition file extension name"}
bool    fl_vardq    {no, prompt="Create variance and data quality frames"}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

begin

    char    l_inimages, l_outimages, l_logfile, l_dir, l_outpref, l_mdf_ext
    bool    l_verbose, l_fl_fixnc, inadjusted, l_fl_vardq
    struct  l_struct
    int     l_test, i, j, k, ninp, nout, nbad, num_extn, l_pixtype
    char    img, filelist, newfilelist, temp1, temp2, temp3, temp4
    char    l_temp, extn[3], inextn, appendextn, l_sci_ext, l_var_ext, l_dq_ext
    int     shuffle, binning, anodcnt, nsci, bnodcnt, dummy, maxfiles
    char    inimg[200], outimg[200], gouttype

    # Initialize
    status = 0
    maxfiles = 200

    # Localize global variables

    l_inimages=inimages; l_outimages=outimages; l_logfile=logfile
    l_verbose=verbose; l_outpref=outpref
    l_fl_fixnc=fl_fixnc
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_mdf_ext = mdf_ext
    l_fl_vardq = fl_vardq

    # Keep task parameters from changing from the outside
    cache ("imgets", "gemdate")

    # Create any needed temporary file names
    filelist = mktemp("tmpfile2")
    newfilelist = mktemp("tmpfile3")

    # Test the logfile:
    if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GNSSKYSUB: Both gnsskysub.logfile and \
                gmos.logfile fields are empty", l_logfile,verbose+)
            printlog ("                    Using default file gmos.log",
                l_logfile,verbose+)
        }
    }

    # Start logging
    date | scan(l_struct)
    printlog ("-----------------------------------------------------------\
        -----------------", l_logfile, l_verbose)
    printlog ("GNSSKYSUB -- "//l_struct, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)
    printlog ("Input list    = "//l_inimages, l_logfile, l_verbose)
    printlog ("Output list   = "//l_outimages, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    # Load up input name list: @list, * and ?, comma separated
    if (l_inimages == "" || stridx(" ",l_inimages)>0) {
        printlog ("ERROR - GNSSKYSUB: Input file not specified", l_logfile, \
            verbose+)
        goto crash
    }

    # Test for @filelist
    if ((substr(l_inimages,1,1) == "@") \
        && !access(substr(l_inimages,2,strlen(l_inimages))) ) {

        printlog ("ERROR - GNSSKYSUB: Input list "//\
            substr(l_inimages,2,strlen(l_inimages))//" does not exist",
            l_logfile,verbose+)
        goto crash
    }

    # parse wildcard and comma-separated lists
    if (substr(l_inimages,1,1)=="@") {
        scanfile = substr(l_inimages,2,strlen(l_inimages))
        while (fscan(scanfile,l_temp) != EOF) {
            files (l_temp, >> filelist)
        }
    } else {
        if (stridx(",",l_inimages)==0)
            files (l_inimages, > filelist)
        else {
            l_test = 9999
            while (l_test!=0) {
                l_test = stridx(",",l_inimages)
                if (l_test>0)
                    files (substr(l_inimages,1,l_test-1), >> filelist)
                else
                    files (l_inimages, >> filelist)
                l_inimages = substr(l_inimages,l_test+1,strlen(l_inimages))
            }
        }
    }
    scanfile = ""
    scanfile = filelist

    # Define the input name list
    ninp = 0
    nbad = 0
    while (fscan(scanfile, img) != EOF) {

        # split off directory path
        fparse (img, verbose-)
        img = fparse.root
        l_dir = fparse.directory

        # Must be there, but cannot be GEIS or OIF
        gimverify (l_dir//img)
        if (gimverify.status==0) {
            ninp = ninp+1
            if (ninp > maxfiles) {
                printlog ("ERROR - GNSSKYSUB: Maximum number of input images \
                    exceeded", l_logfile, verbose+)
                goto crash
            }
            inimg[ninp] = gimverify.outname//".fits" # now has path if relevant
        } else if (gimverify.status==1) {
            printlog ("ERROR - GNSSKYSUB: Input image "//l_dir//img//" does \
                not exist", l_logfile, verbose+)
            nbad+=1
        } else {
            printlog ("ERROR - GNSSKYSUB Input image "//l_dir//img//" is not \
                MEF", l_logfile, verbose+)
            nbad+=1
        }
    } # end while

    if (nbad!=0) {
        printlog ("ERROR - GNSSKYSUB: "//nbad//" input files do not exist, \
            or are the wrong type", l_logfile, verbose+)
        goto crash
    }
    if (ninp==0) {
        printlog ("ERROR - GNSSKYSUB: No input images meet wildcard criteria",
            l_logfile,verbose+)
        goto crash
    }
    # At the end of all this inimg[ninp] now contains the input images incl.
    # the directory path

    # Output name list: can be empty if prefix is defined, @list or
    # comma-separated list
    if (stridx(" ",l_outimages)>0) l_outimages = ""

    if ((substr(l_outimages,1,1) == "@") \
        && !access(substr(l_outimages,2,strlen(l_outimages))) ) {

        printlog ("ERROR - GNSSKYSUB: File "//\
            substr(l_outimages,2,strlen(l_outimages))//" does not exist",
            l_logfile,verbose+)
        goto crash
    }
    files (l_outimages, sort-, > newfilelist)
    scanfile = newfilelist
    nout = 0

    # If empty string, prefix must be defined
    if (l_outimages == "") {
        if (l_outpref == "" || stridx(" ",l_outpref)>0 ) {
            printlog ("ERROR - GNSSKYSUB: Neither output name nor output \
                prefix is defined", l_logfile, verbose+)
            goto crash
        }
        nout = ninp
        i = 1

        while (i<=nout) {
            fparse (inimg[i],verbose=no)
            outimg[i] = l_outpref//fparse.root//".fits"

            gimverify (outimg[i])
            if (gimverify.status!=1) {
                printlog("ERROR - GNSSKYSUB: Output image "//outimg[i]//\
                    " already exists", l_logfile, verbose+)
                nbad+=1
            }
            i = i+1
        } # end of while(i) loop

    } else {
        while (fscan(scanfile, img) != EOF) {
            gimverify (img)
            if (gimverify.status!=1) {
                printlog ("ERROR - GNSSKYSUB: Output image "//img//\
                    " already exists", l_logfile, verbose+)
                nbad+=1
            }
            nout = nout+1
            if (nout > maxfiles) {
                printlog ("ERROR - GNSSKYSUB: Maximum number of output \
                    images exceeded", l_logfile, verbose+)
                goto crash
            }
            outimg[nout] = gimverify.outname//".fits"
        } # end of while (img) loop

    } # end else

    scanfile = ""
    delete (newfilelist, verify-, >& "dev$null")

    if (nbad>0)
        goto crash

    # At the end of all this outimg[nout] now contains the output images
    # Input and output number must be the same, if output names are defined
    if (ninp!=nout && l_outimages!="") {
        printlog ("ERROR - GNSSKYSUB: Number of input and output images \
            are not the same", l_logfile, verbose+)
        goto crash
    }

    # Run a quick check to see if the input images have been gprepared
    i = 1
    while (i<=ninp) {
        imgets (inimg[i]//"[0]", "NSCIEXT", >& "dev$null")
        nsci = int(imgets.value)
        if (nsci==0) {
            printlog ("ERROR - GNSSKYSUB: Keyword NSCIEXT not found in image",
                l_logfile,l_verbose)
            printlog ("ERROR - GNSSKYSUB: File "//inimg[i]//" is not \
                gprepared", l_logfile, l_verbose)
            printlog ("ERROR - GNSSKYSUB: Please run gprepare on all input \
                images", l_logfile, l_verbose)
            goto crash
        }
        i = i+1
    }

    # Start making the output images. This is where all the action really is!
    i = 1
    while (i<=ninp) {
        printlog ("Operating on "//inimg[i]//" to create "//outimg[i],
            l_logfile, l_verbose)

        #get header information
        imgets (inimg[i]//"[0]","NODPIX")
        shuffle = int(imgets.value)

        imgets (inimg[i]//"["//l_sci_ext//",1]","CCDSUM")
        print (imgets.value) | scan(dummy,binning)
        shuffle = shuffle/binning

        imgets (inimg[i]//"[0]","ANODCNT")
        anodcnt = int(imgets.value)

        imgets (inimg[i]//"[0]","BNODCNT")
        bnodcnt = int(imgets.value)

        # Set the way gemarith creates the output extension
        gouttype = "ref"
        keypar (inimg[i]//"["//l_sci_ext//",1]", "i_pixtype", silent+)
        if (keypar.found) {
            l_pixtype = int(keypar.value)
            # If it's a shout or ushort change to int to stop wrap around of
            # negative numbers

            if (l_pixtype == 3 || l_pixtype == 11) {
                gouttype = "int"
            }
        } else {
            printlog ("ERROR - GNSSKYSUB: Cannot determine pixel type", \
                l_logfile, verbose+)
            goto crash
        }

        #check whether the nod counts match. if not we may want renormalize the
        #input later
        if (anodcnt != bnodcnt){
            printlog ("WARNING - GNSSKYSUB: Nod counts for A and B positions \
                do not match.", l_logfile, verbose+)
            if (l_fl_fixnc)
                printlog ("WARNING - GNSSKYSUB: Normalizing images to fix nod \
                    count mismatch.", l_logfile, verbose+)
            else
                printlog ("WARNING - GNSSKYSUB: Ignoring the nod count \
                    mismatch.", l_logfile,verbose+)
        }

        inadjusted = no

        # A and B nod counts mismatch and user wants to adjust them
        if (l_fl_fixnc && anodcnt != bnodcnt) {

            temp1 = mktemp("tmp1")//".fits"
            temp2 = mktemp("tmp2")//".fits"

            gouttype = "ref"
            keypar (inimg[i]//"["//l_sci_ext//",1]", "i_pixtype", silent+)
            if (keypar.found) {
                l_pixtype = int(keypar.value)
                # If it's a shout or ushort change to int to stop wrap around
                # of negative numbers

                if (l_pixtype == 7 || l_pixtype == 3) {
                    gouttype = "real"
                }
            } else {
                printlog ("ERROR - GNSSKYSUB: Cannot determine pixel type", \
                    l_logfile, verbose+)
                goto crash
            }

            # Divide inout images by the nodcount b
            gemarith (inimg[i], op="/", operand2=bnodcnt, \
                result=temp1, sci_ext=l_sci_ext, var_ext=l_var_ext, \
                dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, fl_vardq=l_fl_vardq, \
                dims="default", \
                intype="double", outtype=gouttype, refim="default", \
                rangecheck=yes, verbose=no, logfile=l_logfile)

            # Divide inout images by the nodcount a
            gemarith (inimg[i], op="/", operand2=anodcnt, \
                result=temp2, sci_ext=l_sci_ext, var_ext=l_var_ext, \
                dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, fl_vardq=l_fl_vardq, \
                dims="default", \
                intype="double", outtype=gouttype, refim="default", \
                rangecheck=yes, verbose=no, logfile=l_logfile)

            # Set flag for later use and output file for subtraction
            inadjusted = yes
            temp4 = mktemp("tmp4")
        } else {
            # Set temp names to use later on
            temp1 = inimg[i]
            temp2 = inimg[i]
            # Set output name for subtraction
            temp4 = outimg[i]
        }

        # imshifted temp1 - output file
        temp3 = mktemp("tmp3")//".fits"

        # Set up extensions to copy
        extn[1] = l_sci_ext
        extn[2] = l_var_ext
        extn[3] = l_dq_ext

        # Copy over PHU
        imcopy (temp1//"[0]", temp3, verbose-)

        # Check for an MDF and attach it
        gemextn (temp1, process="expand", check="exists", \
            extver="", index="1-", extname=l_mdf_ext, \
            outfile="dev$null", \
            logfile=l_logfile, verbose=l_verbose)

        if (gemextn.count == 1) {
            tcopy (temp1//"["//l_mdf_ext//"]", \
                temp3//"["//l_mdf_ext//"]", verbose=no)
        }

        # Always do SCI
        num_extn = 1

        # Check for VAR and DQ planes and attach them
        gemextn (temp1, process="expand", check="exists,mef", \
            extver="1-", index="", extname=l_var_ext//","//l_dq_ext, \
            outfile="dev$null", omit="index", \
            logfile=l_logfile, verbose=l_verbose)

        # Increase the number of indicies to loop over in extn to include VAR
        # and DQ - MS
        if ((gemextn.count > 0) && (gemextn.count / 2) == nsci) {
            num_extn = 3
        }

        # Shift the required extensions by shuffle
        for (j = 1; j <= num_extn; j += 1) {
            for (k=1; k <= nsci; k += 1) {
                inextn = "["//extn[j]//","//k//"]"
                appendextn = "["//extn[j]//","//k//",append]"

                imshift (input=temp1//inextn, \
                    output=temp3//appendextn, \
                    xshift=0, yshift=shuffle, interp="nearest", \
                    shifts_file="", boundary_type="nearest", constant=0.)
            }
        }

        # Subtract the sky
        gemarith (temp2, op="-", operand2=temp3, \
            result=temp4, sci_ext=l_sci_ext, var_ext=l_var_ext, \
            dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, fl_vardq=l_fl_vardq, \
            dims="default", \
            intype="double", outtype=gouttype, refim="default", \
            rangecheck=yes, verbose=no, logfile=l_logfile)

        # Image a/b nodcnt do not match - created scaled versions for
        # subtraction - rescale sky subtracted image here - MS
        if (inadjusted) {
            # Scale the image back to the greater of anodcount / bnodcount

            gouttype = "ref"
            keypar (inimg[i]//"["//l_sci_ext//",1]", "i_pixtype", silent+)
            if (keypar.found) {
                l_pixtype = int(keypar.value)
                # If it's a shout or ushort change to int to stop wrap around
                # of negative numbers

                if (l_pixtype == 7 || l_pixtype == 3) {
                    gouttype = "real"
                }
            } else {
                printlog ("ERROR - GNSSKYSUB: Cannot determine pixel type", \
                    l_logfile, verbose+)
                goto crash
            }

            gemarith (temp4, op="*", operand2=(max(anodcnt,bnodcnt)), \
                result=outimg[i], sci_ext=l_sci_ext, var_ext=l_var_ext, \
                dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, fl_vardq=l_fl_vardq, \
                dims="default", \
                intype="double", outtype=gouttype, refim="default", \
                rangecheck=yes, verbose=no, logfile=l_logfile)
            imdelete (temp1//", "//temp2//", "//temp4, verify-, >& "dev$null")
        }

        # Delete shifted file
        imdelete (temp3, verify-, >& "dev$null")

        gemdate ()
        gemhedit (outimg[i]//"[0]", "GNSSKYSU", gemdate.outdate,
            "UT Time stamp for GNSSKYSUB", delete-)
        gemhedit (outimg[i]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)

        printlog ("", l_logfile, verbose+)
        i = i+1
    }

    # exit status good
    goto clean

crash:
    # Exit with error subroutine
    status = 1
    goto clean

clean:
    # clean up
    delete (filelist, verify-, >& "dev$null")
    delete (newfilelist, verify-, >& "dev$null")
    scanfile = ""

    #close log file
    if (status == 0)
        printlog ("GNSSKYSUB exit status: good.", l_logfile, l_verbose)
    else
        printlog ("GNSSKYSUB exit status: error.", l_logfile, l_verbose)

    printlog ("---------------------------------------------------------\
        -------------------", l_logfile, l_verbose)

end
