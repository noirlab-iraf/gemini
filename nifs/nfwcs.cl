# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.
#
# NFWCS -- Set WCS information used by TRANSCUBE.
#
# This routine sets a 3D WCS for the IFU slices.  The WCS is a linear system
# in arcsec from some reference point.  It currently uses the MDF (a FITS
# binary table) to define the reference points as well as slit width.
# The MDF table is indexed by the extension version (i.e. EXTVER = row).
# The MDF may be given explicitly or by a header keyword.  The pixel
# scale and any offsets are defined explicitly or by header keywords.
# The convention is that X refers to the the coordinate along the slit and
# Y refers to the coordinate across the slit (i.e. the slice direction).

procedure nfwcs (images)

string	images			{prompt="List of NIFS science slit extensions"}
string	mdf = "!mdf_file"	{prompt="MDF table"}
string	xmdfname = "DEC"	{prompt="MDF column for X position (arcsec)"}
string	ymdfname = "RA"		{prompt="MDF column for Y position (arcsec)"}
string	dymdfname = "slitsize_x"{prompt="MDF column for slit width (arcsec)"}
string	pixscale = "!pixscale"	{prompt="Pixel scale along slit (arcsec/pix)"}
string	xoffset = "0."		{prompt="X offset (arcsec)"}
string	yoffset = "0."		{prompt="Y offset (arcsec)"}
bool	keepwcs = no		{prompt="Keep existing 2D WCS?"}

struct	*scanfile

begin
	file	imlist, im, im1, m
	int	i, extver
	real	x, y, dx, dy, xoff, yoff
	struct	colnames

	# Expand input list and set other input parameters.
	imlist = mktemp ("tmp")
	sections (images, option="fullname", > imlist)

	if (fscan (mdf, m) == 0)
	    error (1, "MDF file not specified")
	if (stridx ("!", pixscale) == 0)
	    dx = real (pixscale)
	if (stridx ("!", xoffset) == 0)
	    xoff = real (xoffset)
	if (stridx ("!", yoffset) == 0)
	    yoff = real (yoffset)
	printf ("%s %s %s\n", xmdfname, ymdfname, dymdfname) | scan (colnames)

	# Loop through list of image extensions.
	scanfile = imlist
	while (fscan (scanfile, im) != EOF) {
	    # Add inherit to get to global header.
	    im1 = im
	    i = stridx ("]", im1)
	    if (i > 0)
	        im1 = substr (im1,1,i-1) // ",inherit" // substr (im1,i,1000)

	    # Get MDF information.
	    # Add fits extension required by TABLES package.
	    i = stridx ("!", mdf)
	    if (i > 0) {
	        m = ""; hselect (im1, substr(mdf,i+1,1000), yes) | scan (m)
		if (m == "")
		    error (1, "MDF file not specified")
	    }
	    i = strstr (".fits", m)
	    if (i == 0) {
		i = strldx ("[", m)
		if (i == 0)
		    m = m // ".fits"
		else
		    m = substr (m,1,i-1) // ".fits" // substr (m,i,1000)
	    }

	    hselect (im1, "extver", yes) | scan (extver)
	    tdump (m, cdfile="", pfile="", datafile="STDOUT",
		columns=colnames, rows=extver, pwidth=-1) |
		scan (x, y, dy)
	    x *= 3600.
	    y *= 3600.

	    # Get pixel scale along the slit.
	    i = stridx ("!", pixscale)
	    if (i > 0) {
	        dx = INDEF; hselect (im1, substr(pixscale,i+1,1000), yes) |
		    scan (dx)
		if (isindef(dx))
		    error (1, "Pixel scale value not found")
	    }

	    # Get offsets.
	    i = stridx ("!", xoffset)
	    if (i > 0) {
	        xoff = INDEF; hselect (im1, substr(xoffset,i+1,1000), yes) |
		    scan (xoff)
		if (isindef(xoff))
		    error (1, "X offset value not found")
	    }
	    x += xoff

	    i = stridx ("!", yoffset)
	    if (i > 0) {
	        yoff = INDEF; hselect (im1, substr(yoffset,i+1,1000), yes) |
		    scan (yoff)
		if (isindef(yoff))
		    error (1, "Y offset value not found")
	    }
	    y += yoff

	    # Clear any previous WCS keywords.
	    if (!keepwcs)
		hedit (im, "crpix*,crval*,cd1*,cd2*,cd3*", del+,
		    add-, addonly-, verify-, show-, update+)
	    hedit (im, "ltv*,ltm*", del+,
		add-, addonly-, verify-, show-, update+)
	    hedit (im, "wat*,waxmap*", del+,
		add-, addonly-, verify-, show-, update+)

	    # Set dimensionality and axis mapping.
        gemhedit (im, "wcsdim", 3, "", delete-)
        gemhedit (im, "waxmap01", "1 0 2 0 0 0", "", delete-)

	    # Set dispersion axis.
        gemhedit (im, "ctype1", "LINEAR", "", delete-)
        gemhedit (im, "WAT1_001", "wtype=linear axtype=wave", "", delete-)
	    if (!keepwcs)
            gemhedit (im, "cd1_1", 1, "", delete-)

	    # Set axis along the slit.
        gemhedit (im, "ctype2", "LINEAR", "", delete-)
	    gemhedit (im, "WAT2_001", "wtype=linear axtype=eta", "", delete-)
	    gemhedit (im, "crval2", x, "", delete-)
        if (!keepwcs)
            gemhedit (im, "cd2_2", dx, "", delete-)
        else
            gemhedit (im, "cd2_2", "("//dx//"*cd2_2)", "", delete-)

	    # Set axis across the slit.
        gemhedit (im, "ctype3", "LINEAR", "", delete-)
        gemhedit (im, "WAT3_001", "wtype=linear axtype=xi", "", delete-)
        gemhedit (im, "crval3", y, "", delete-)
	    gemhedit (im, "cd3_3", dy, "", delete-)
	}
	scanfile = ""; delete (imlist, verify-)
end
