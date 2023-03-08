# Copyright(c) 2012 Association of Universities for Research in Astronomy, Inc.

# F2DISPLAY - Display cut spectroscopic FLAMINGOS-2 data
# Full details of this task are given in the help file

procedure f2display(inimage, frame)

char    inimage     {prompt = "Input cut FLAMINGOS-2 image"}
int     frame       {1, prompt = "Frame to be written to"}
char    extname     {"SCI", prompt = "Extension type to display"}
int     status      {0, prompt = "Exit status (0=good)"}

begin

    char    l_inimage, l_extname, badhdr, extn, section
    int     junk, l_frame, nsciext, ext, nx, ny, x1, x2, y1, y2
    real    xcenter, ycenter

    char    l_key_cut_section = ""

    status = 1

    # Cache tasks used throughout the script
    cache ("keypar", "gemdate")

    # Read input parameters
    junk = fscan (inimage, l_inimage)
    l_frame = frame
    junk = fscan (extname, l_extname)

    # Shared definitions, define parameters from nsheaders
    junk = fscan (nsheaders.key_cut_section, l_key_cut_section)
    if ("" == l_key_cut_section) badhdr = badhdr + " key_cut_section"

    # Check if the input image exists
    gemextn (l_inimage, check="exists,mef,image", process="append", \
        index="", extname=l_extname, extversion="1", ikparams="", \
        omit="", replace="", outfile="dev$null", logfile="", glogpars="",
        verbose-)
    if (gemextn.fail_count != 0 || gemextn.count == 0) {
        print ("ERROR - F2DISPLAY: Cannot display " // l_inimage // ": no " \
            // l_extname // "data.")
        goto clean
    }

    # Only cut data will be displayed
    keypar (l_inimage // "[0]", "F2CUT", silent+)
    if (!keypar.found) {
        print ("ERROR - F2DISPLAY: " // l_inimage // " has not been cut.")
        goto clean
    }

    # Determine the number of science extension in the input image
    keypar (l_inimage // "[0]", "NSCIEXT", silent+)
    if (keypar.found) {
        nsciext = int(keypar.value)
    } else {
        print ("ERROR - F2DISPLAY: Cannot determine the number of science \
            extensions.")
        goto clean
    }
    
    # Set the IRAF stdimage variable
    set stdimage = "imt2048"

    # Loop over science extensions and display each one
    for (ext = 1; ext <= nsciext; ext += 1) {

        # Set the extension to display
        extn = "[" // l_extname // "," // ext // "]"

        # Get the section to display
        keypar (l_inimage // extn, l_key_cut_section, silent+)
        if (keypar.found) {
            section = keypar.value
        } else {
            print ("ERROR - F2DISPLAY: Cannot determine section for \
                extension " // extn)
            goto clean
        }

        print (section)
        print (section) | scanf ("[%d:%d,%d:%d]", x1, x2, y1, y2)

        # Determine the display center for this extension
        xcenter = (((x2 - x1) / 2.) + x1) / 2048
        ycenter = (((y2 - y1) / 2.) + y1) / 2048

        if (ext == 1) {
            # Clear the frame before displaying the first extension
            display (l_inimage // extn, frame=l_frame, bpmask="BPM", \
                bpdisplay="none", bpcolors="red", overlay="", \
                ocolors="green", erase+, border_erase-, select_frame+, \
                repeat-, fill-, zscale+, contrast=0.25, zrange+, zmask="", \
                nsample=1000, xcenter=xcenter, ycenter=ycenter, xsize=1.0, \
                ysize=1.0, xmag=1.0, ymag=1.0, order=0, z1=INDEF, z2=INDEF, \
                ztrans="linear", lutfile="")
        } else {
            # Display the extensions in the current frame
            display (l_inimage // extn, frame=l_frame, bpmask="BPM", \
                bpdisplay="none", bpcolors="red", overlay="", \
                ocolors="green", erase-, border_erase-, select_frame+, \
                repeat-, fill-, zscale+, contrast=0.25, zrange+, zmask="", \
                nsample=1000, xcenter=xcenter, ycenter=ycenter, xsize=1.0, \
                ysize=1.0, xmag=1.0, ymag=1.0, order=0, z1=INDEF, z2=INDEF, \
                ztrans="linear", lutfile="")
        }
    }

    # If we have managed to arrive here, everything has worked correctly
    status = 0

clean:
    if (status == 1) {
        print ("")
        print ("F2DISPLAY exit status: failed.")
    }

end
