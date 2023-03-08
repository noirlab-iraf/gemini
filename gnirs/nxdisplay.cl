# Copyright(c) 2004-2011 Association of Universities for Research in Astronomy, Inc.

procedure nxdisplay (image, frame)
 
# CW - 2004 Jan 30
# Quick and dirty procedure to display all extensions
# of a GNIRS processed image at once

string    image    {prompt="Raw image"}
int       frame    {1, prompt = "Frame to be written to"}

begin

    # local variables
    string  l_image, grating, extn, ext_ver
    bool    erase_frame
    int     junk, l_frame, num_sci_ext, ii, status, x_width, total_width
    int     y_height
    real    current_xcentre, x_offset, xstart, xgap, xboarder

    # Read input parameters
    junk = fscan (image, l_image)
    l_frame = frame

    # Define default status
    status = 0

    # Define the extension type to plot
    ext_ver = "SCI"

    # Read the number of science extensions
    keypar (l_image//"[0]", "NSCIEXT", silent+)
    if (!keypar.found) {
        print ("NXDISPLAY ERROR - "//l_image//" has not been processed."//\
            " Cannot display image. Exiting.")
        goto crash
    } else {
        num_sci_ext = int(keypar.value)
    }
        
    # Set the IRAF stdimage variable
    set stdimage = "imt2048"

    total_width = 0

    # Loop over science extensions to obtain full width off extensions
    for (ii = 1; ii <= num_sci_ext; ii += 1) {

        # Set the current extension to display
        extn = "[SCI,"//ii//"]"
 
        # Read the NAXIS2 keyword
        keypar (l_image//extn, "i_naxis1", silent+)
        x_width = int(keypar.value)

        # Accumilate the total width
        total_width += x_width
    }

    # Perform calculations to centre the image

    # Find initial boarder
    xboarder = ((2048.0 - total_width) / 2048.0) / 2.0

    # Calculate gap between extensions
    xgap = real(xboarder / ( (num_sci_ext - 1)))

    # Calulate the starting point (the new boader)
    xstart = xboarder - (xgap * ((num_sci_ext - 1)/2.0))

    # Initiate the new starting point
    current_xcentre = xstart

    # Loop over SCI extensions and plot 
    for (ii = 1; ii <= num_sci_ext; ii += 1) {

        # Set the current extension to display
        extn = "["//ext_ver//","//ii//"]"

        # Read the width of the extension again
        keypar (l_image//extn, "i_naxis1", silent+)
        x_width = int(keypar.value)

        # Calculate the centering offset for this extension - half
        # the width of the current extension
        x_offset = real (x_width) / (2.0 * 2048.0)

        # Redefine teh centre for this extension
        current_xcentre += x_offset

        # Old xcenter values 0.15, 0.28, 0.41, 0.54, 0.67, 0.80
        
        # Set the display erase parameter according to extension number
        if (ii == 1) {
            erase_frame = yes
        } else {
            erase_frame = no
        }

        # display the current science extension
        display (l_image//extn, l_frame, bpdisplay="none", overlay="", \
            erase=erase_frame, border_erase=no, select_frame=yes, \
            repeat=no, fill=no,  zscale=no, contrast=0.25, zrange=yes, \
            nsample=500, xcenter=current_xcentre, \
            ycenter=0.5, xsize=1., ysize=1., xmag=1.0, ymag=2., order=0, \
            z1=INDEF, z2=INDEF, ztrans="linear", lutfile="")

        # Define shift the xcentre to the left edge of the next extension
        current_xcentre +=  x_offset + xgap
    }

    goto clean

crash:
    status = 1

clean:

    if (status != 1) {
        status = 0
    }

end
