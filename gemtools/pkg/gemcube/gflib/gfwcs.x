# GFWCS -- WCS and STACK geometry function driver for TRANSCUBE.
#
# This driver works with images that are 2D or 3D and that have a
# 2D or 3D WCS.  The case when the images are 2D but the WCS is 3D combines
# the 2D slices into a 3D data cube.  This means it can be used to build data
# cubes from multiple spatial images sampled along some other coordinate
# (typically time or wavelength) or from multiple long slit images with an
# implicit second spatial axis.  The long slit case would typically be data
# from an image slicer IFU.
#
# The distinction between the WCS and STACK driver entry points are
# whether 2D images with 2D WCS are combined into a 2D output or
# "stacked" into a 3D output using the input list index as the
# coordinate for the third dimension.
# 
# This driver also supports FITCOORDS-style distortion functions.

include	<mach.h>
include	<math.h>
include	<error.h>
include	<imhdr.h>
include	<mwset.h>
include	<ctotok.h>
include	<math/gsurfit.h>

# GF structure.
define	GF_LENSHAPE	100			# Length of shape string
define	GF_LEN		66			# Length of structure
define	GF_PP		Memi[$1]		# PSET pointer
define	GF_MW		Memi[$1+1]		# MWCS pointer
define	GF_MWDIM	Memi[$1+2]		# WCS dimensionality
define	GF_CTLW		Memi[$1+3]		# Logical -> world
define	GF_AXES		Memi[$1+3+$2]		# In to out axes mapping (3)
define	GF_WTAXES	Memi[$1+6+$2]		# Wt to out axes mapping (3)
define	GF_UN		Memi[$1+9+$2]		# Units pointers (2)
define	GF_SF		Memi[$1+11+$2]		# Arrays of surface ptrs (2)
define	GF_NSF		Memi[$1+13+$2]		# Number of surface ptrs (2)
define	GF_SHAPE	Memc[P2C($1+16)]	# Shape string

# Axes types.
define	GF_AX1	"|ra|xi|"		# First spatial axis
define	GF_AX2	"|dec|eta|"		# Second spatial axis
define	GF_AX3	"|wave|lambda|"		# Non-spatial axis

# GFWCS_OPEN -- WCS open procedure.

procedure gfwcs_open (gf, images, index, pset)	# PUBLIC

pointer	gf				#I Geometry function pointer
pointer	images				#I List of images
int	index				#I Index of image in list to open
char	pset[ARB]			#I PSET name

begin
	call gfw_open (gf, images, index, pset, INDEFI, 1)
end


# GFWCS_OUT -- WCS output suggestion procedure.

procedure gfwcs_out (images, refmw, pset, mw, imlen)	# PUBLIC

pointer	images				#I List of images
pointer	refmw				#I Reference WCS
char	pset[ARB]			#I PSET name
pointer	mw				#O Suggested WCS
int	imlen[3]			#O Suggested image size

begin
	call gfw_out (images, refmw, pset, mw, imlen, NO)
end


# GFSTACK_OPEN -- Stack open procedure.
#
# This is like the WCS driver functions except that a list of 2D images are
# stacked using their index in the input list.

procedure gfstackwcs_open (gf, images, index, pset)	# PUBLIC

pointer	gf				#I Geometry function pointer
pointer	images				#I List of images
int	index				#I Index of image in list to open
char	pset[ARB]			#I PSET name

begin
	call gfw_open (gf, images, index, pset, 3, index)
end


# GFSTACK_OUT -- STACK output suggestion procedure.
#
# This is like the WCS driver functions except that a list of 2D images are
# stacked using their index in the input list.

procedure gfstack_out (images, refmw, pset, mw, imlen)	# PUBLIC

pointer	images				#I List of images
pointer	refmw				#I Reference WCS
char	pset[ARB]			#I PSET name
pointer	mw				#O Suggested WCS
int	imlen[3]			#O Suggested image size

begin
	call gfw_out (images, refmw, pset, mw, imlen, YES)
end


# GFWCS_CLOSE -- WCS close procedure.

procedure gfwcs_close (gf)	# PUBLIC

pointer	gf				#U Object pointer

int	i, j

begin
	# Close the distortion functions.
	do i = 1, 2 {
	    if (GF_UN(gf,i) != NULL)
	        call un_close (GF_UN(gf,i))
	    do j = 1, GF_NSF(gf,i)
	        call dgsfree (Memi[GF_SF(gf,i)+j-1])
	    call mfree (GF_SF(gf,i), TY_POINTER)
	}

	# Close the WCS.
	if (GF_MW(gf) != NULL)
	    call mw_close (GF_MW(gf))

	# Close the PSET.
	if (GF_PP(gf) != NULL)
	    call clcpset (GF_PP(gf))

	# Free the structure.
	call mfree (gf, TY_STRUCT)
end


# GFWCS_PIXEL -- WCS pixel procedure.
#
# This expects 3D coordinates but it will evaluate 2D WCS.
# It maps the input axes to the output orienation.
# It supports FITCOORDS distortion functions.

procedure gfwcs_pixel (gf, pixel, world)	# PUBLIC

pointer	gf				#I Object pointer
double	pixel[3]			#I Pixel coordinate
double	world[3]			#O World coodinate

int	i, j, n
double	p[3], w[3], dgseval()

begin
	# Compute transformation.

	# Apply distortion functions if defined.
	if (GF_NSF(gf,1) != 0 || GF_NSF(gf,2) != 0) {
	    p[1] = pixel[1]; p[2] = pixel[2]; p[3] = pixel[3]
	    do i = 1, 2 {
		n = GF_NSF(gf,i)
		if (n == 1)
		    p[i] = dgseval (Memi[GF_SF(gf,i)], pixel[1], pixel[2])
		else if (n > 1) {
		    p[i] = 0.
		    do j = 1, n
			p[i] = p[i] +
			    dgseval (Memi[GF_SF(gf,i)+j-1], pixel[1], pixel[2])
		    p[i] = p[i] / n
		}
	    }
	    call mw_ctrand (GF_CTLW(gf), p, w, GF_MWDIM(gf))
	} else
	    call mw_ctrand (GF_CTLW(gf), pixel, w, GF_MWDIM(gf))

	# Map the input axes to the output axes.
	world[1] = w[GF_AXES(gf,1)]
	world[2] = w[GF_AXES(gf,2)]
	if (GF_MWDIM(gf) == 2)
	    world[3] = pixel[3]
	else
	    world[3] = w[GF_AXES(gf,3)]
end


# GFWCS_GEOM -- WCS geometry procedure.

procedure gfwcs_geom (gf, pixel, shape, axmap)	# PUBLIC

pointer	gf				#I Object pointer
double	pixel[3]			#I Pixel coordinate
char	shape[1024]			#O Shape string
int	axmap[3]			#O Axis map

begin
	call strcpy (GF_SHAPE(gf), shape, 1024)
	call amovi (GF_WTAXES(gf,1), axmap, 3)
end


#==============================================================================#


# GFW_OPEN -- WCS and STACK open procedure.
#
# The is the core procedure for setting up the structure for calling the
# other geometry functions.  This procedure accepts 
# This sets up the driver for calls to the other functions.

procedure gfw_open (gf, images, iin, pset, omwdim, iout)	# PRIVATE

pointer	gf				#I Geometry function pointer
pointer	images				#I List of images
int	iin				#I Index of input image in list to open
char	pset[ARB]			#I PSET name
int	omwdim				#I Output dimensionality
int	iout				#I Output index for stacking

int	i, j, mwdim, omwdim1, axis[3]
double	r[3], w[3], cd[9], cd1[9]
pointer	sp, image, key, str, pp, im, mw, ctlw, inmw, tmp

bool	strne()
int	clgeti(), clgpseti(), mw_stati(), imtrgetim()
pointer	clopset(), immap(), mw_openim(), mw_sctran(), mw_open()
errchk	clopset, immap, mw_openim, mw_open, calloc, gfwcs_sf, gfwcs_gwtermd

begin
	call smark (sp)
	call salloc (image, SZ_FNAME, TY_CHAR)
	call salloc (key, SZ_FNAME, TY_CHAR)
	call salloc (str, SZ_FNAME, TY_CHAR)

	iferr {
	    pp = NULL; im = NULL; mw = NULL; gf = NULL

	    # Open the pset.
	    pp = clopset (pset)

	    # Set the WCS from the image.
	    i = imtrgetim (images, iin, Memc[image], SZ_FNAME)
	    tmp = immap (Memc[image], READ_ONLY, 0); im = tmp
	    tmp = mw_openim (im); mw = tmp
	    call mw_seti (mw, MW_USEAXMAP, NO)
	    mwdim = mw_stati (mw, MW_NPHYSDIM)
	    ctlw = mw_sctran (mw, "logical", "world", 0)

	    # For now we only allow input images to be promoted to
	    # a higher dimension.  Some day we might allow reduction.
	    if (IS_INDEFI(omwdim))
	        omwdim1 = mwdim
	    else
	        omwdim1 = omwdim
	    if (omwdim1 < mwdim)
	        call error (2,
		    "Output WCS dimensionality may not be less than the input")

	    # Set the GF structure.
	    call calloc (tmp, GF_LEN, TY_STRUCT); gf = tmp
	    GF_PP(gf) = pp
	    GF_MWDIM(gf) = mwdim
	    GF_MW(gf) = mw
	    GF_CTLW(gf) = ctlw
	    do i = 1, 3
	        GF_AXES(gf,i) = i

	    # Get the distortion functions if defined.
	    call gfwcs_sf (gf, im)

	    # Promote the WCS dimensionality if needed.
	    if (omwdim1 == 3 && mwdim == 2) {

		# Create new WCS.
		inmw = mw
		mw = mw_open (NULL, 3)
		call mw_newsystem (mw, "world", 3)

		# Add axis with defaults.
		call mw_gwtermd (inmw, r, w, cd1, mwdim)
		call aclrd (cd, 9)
		do i = 1, mwdim
		    do j = 1, mwdim
		        cd[(i-1)*3+j] = cd1[(i-1)*mwdim+j]
		r[3] = 1 - (iout - 1)
		w[3] = 1.
		cd[9] = 1.
		j = 0

		# Set the projection.
		Memc[key] = EOS
		do i = 1, mwdim {
		    ifnoerr (call mw_gwattrs (inmw, i, "wtype",
			Memc[str], SZ_FNAME)) {
			if (strne (Memc[str], "linear")) {
			    call strcpy (Memc[str], Memc[key], SZ_FNAME)
			    j = j + 1
			    axis[j] = i
			} else
			    call mw_swtype (mw, i, 1, Memc[str], "")
		    }
		}
		if (j > 0 && Memc[key] != EOS)
		    call mw_swtype (mw, axis, j, Memc[key], "")

		# Set the axes types.
		do i = 1, mwdim {
		    ifnoerr (call mw_gwattrs (inmw, i, "axtype",
			Memc[str], SZ_FNAME))
			call mw_swattrs (mw, i, "axtype", Memc[str])
		}

		# Set the world transformation.
		call mw_swtermd (mw, r, w, cd, 3)

		# Reset elements of the GF structure.
		mwdim = mw_stati (mw, MW_NPHYSDIM)
		ctlw = mw_sctran (mw, "logical", "world", 0)
		GF_MWDIM(gf) = mwdim
		GF_MW(gf) = mw
		GF_CTLW(gf) = ctlw

		call mw_close (inmw)
	    }

	    # Define the mapping between input and output axes.
	    i = clgpseti (pp, "axis1")
	    if (!IS_INDEFI(i))
	        GF_AXES(gf,1) = i
	    i = clgpseti (pp, "axis2")
	    if (!IS_INDEFI(i))
	        GF_AXES(gf,2) = i
	    i = clgpseti (pp, "axis3")
	    if (!IS_INDEFI(i))
	        GF_AXES(gf,3) = i
	    do i = 1, mwdim {
	        do j = 1, mwdim
		    if (i != j && GF_AXES(gf,i) == GF_AXES(gf,j))
		        call error (1, "Bad axis mapping")
	    }

	    # Determine the non-spatial axis, i.e. uncoupled.
	    GF_WTAXES(gf,3) = clgeti ("nonspatial")
	    if (GF_WTAXES(gf,3) == 1) {
		GF_WTAXES(gf,1) = 2
		GF_WTAXES(gf,2) = 3
	    } else if (GF_WTAXES(gf,2) == 1) {
		GF_WTAXES(gf,1) = 1
		GF_WTAXES(gf,2) = 3
	    } else {
		GF_WTAXES(gf,1) = 1
		GF_WTAXES(gf,2) = 2
	    }
		    
	    # Compute the pixel size and orientation in the output world system.
	    call gfwcs_gwtermd (gf, r, w, cd)
	    do i = 1, mwdim {
		w[i] = 0.
		do j = 1, mwdim
		    w[i] = w[i] + cd[(i-1)*mwdim+j]**2
		w[i] = sqrt (w[i])
	    }
	    if (mwdim < 3)
		w[3] = 0.
	    r[1] = RADTODEG(atan2(cd[(GF_WTAXES(gf,2)-1)*mwdim+GF_WTAXES(gf,1)],
	        cd[(GF_WTAXES(gf,2)-1)*mwdim+GF_WTAXES(gf,2)]))

	    # Set the weight information.
	    call sprintf (GF_SHAPE(gf), GF_LENSHAPE, "rectangle %g %g %g %g")
	        call pargd (w[GF_WTAXES(gf,1)])
	        call pargd (w[GF_WTAXES(gf,2)])
	        call pargd (w[GF_WTAXES(gf,3)])
		call pargd (r[1])

	    call imunmap (im)
	} then {
	    if (gf != NULL)
	        call gfwcs_close (gf)
	    if (im != NULL)
	        call imunmap (im)
	    call erract (EA_ERROR)
	}

	call sfree (sp)
end


# GFW_OUT -- WCS and STACK output suggestion procedure.
#
# This function is used to suggest a WCS and return the output pixel size
# that represents a list of input images.  When a reference WCS is provided
# the returned WCS is a copy but with the reference point adjusted for
# an image origin that covers all the input data.

procedure gfw_out (images, refmw, pset, mw, imlen, stack)	# PRIVATE

pointer	images				#I List of images
pointer	refmw				#I Reference WCS
pointer	mw				#O Suggested WCS
char	pset[ARB]			#I PSET name
int	imlen[3]			#O Suggested image size
int	stack				#I Stack 2D images?

int	i, j, k, l, m, n, nim, minlen[3], maxlen[3]
int	mwdim, omwdim, axis[3], axes[3]
double	r[3], w[3], cd[9], x[3]
double	r1[3], w1[3], cd1[9], cmin[3], cmax[3], scale, minscale
pointer	sp, image, key, str, pp, im, inmw, gf, tmp, ctwl

bool	streq(), strne()
int	imtrgetim(), imtlen(), strdic(), nscan(), nowhite(), ctod(), mw_stati()
double	clgpsetd()
pointer	clopset(), immap(), mw_newcopy(), mw_open(), mw_sctran()
errchk	clopset, immap, calloc, gfwcs_gwtermd

begin
	call smark (sp)
	call salloc (image, SZ_FNAME, TY_CHAR)
	call salloc (key, SZ_FNAME, TY_CHAR)
	call salloc (str, SZ_FNAME, TY_CHAR)

	# Open parameter set.
	pp = clopset (pset)

	# Set the suggested output WCS.
	if (refmw != NULL) {
	    mw = mw_newcopy (refmw)
	    ctwl = mw_sctran (mw, "world", "logical", 0)
	} else {
	    # Find the WCS with the minimum scale.
	    # In this version the minimum scale means the norm of the CD matrix.
	    do nim = 1, imtlen (images) {
		if (stack == YES)
		    call gfw_open (gf, images, nim, pset, 3, nim)
		else
		    call gfw_open (gf, images, nim, pset, INDEFI, 1)
		call gfwcs_gwtermd (gf, r1, w1, cd1)
		if (nim == 1) {
		    mwdim = GF_MWDIM(gf)
		    minscale = 0.
		    do i = 1, mwdim*mwdim
			minscale = minscale + cd1[i] * cd1[i]
		    call amovi (GF_AXES(gf,1), axes, 3)
		    inmw = mw_newcopy (GF_MW(gf))
		    call amovd (r1, r, mwdim)
		    call amovd (w1, w, mwdim)
		    call amovd (cd1, cd, mwdim*mwdim)
		} else {
		    if (GF_MWDIM(gf) != mwdim)
			call error (1,
			    "Input images have different WCS dimensions")
		    scale = 0.
		    do i = 1, mwdim*mwdim
			scale = scale + cd1[i] * cd1[i]
		    if (scale < minscale) {
		        minscale = scale
			call mw_close (inmw)
			call amovi (GF_AXES(gf,1), axes, 3)
			inmw = mw_newcopy (GF_MW(gf))
			call amovd (r1, r, mwdim)
			call amovd (w1, w, mwdim)
			call amovd (cd1, cd, mwdim*mwdim)
		    }
		}

		call gfwcs_close (gf)
	    }

	    # Now create the output WCS.
	    # This is complicated to allow any combination of input
	    # logical and world axes.

	    mw = mw_open (NULL, mwdim)
	    call mw_newsystem (mw, "world", mwdim)

	    # Set any projection.
	    k = 0
	    do i = 1, mwdim {
		ifnoerr (call mw_gwattrs (inmw, axes[i], "wtype",
		    Memc[str], SZ_FNAME)) {
		    if (strne (Memc[str], "linear")) {
			call strcpy (Memc[str], Memc[key], SZ_FNAME)
			k = k + 1
			axis[k] = i
		    } else
			call mw_swtype (mw, i, 1, Memc[str], "")
		}
	    }
	    if (k > 0)
		call mw_swtype (mw, axis, k, Memc[key], "")

	    # Set the axes types.
	    k = 0
	    do i = 1, mwdim {
		ifnoerr (call mw_gwattrs (inmw, axes[i], "axtype",
		    Memc[str], SZ_FNAME)) {
		    call mw_swattrs (mw, i, "axtype", Memc[str])
		    if (strdic(Memc[str],Memc[str],SZ_FNAME,GF_AX1)>0)
		        k = k + 1
		    else if (strdic(Memc[str],Memc[str],SZ_FNAME,GF_AX2)>0)
		        k = k + 1
		}
	    }

	    do i = 1, mwdim {
		do j = 1, mwdim {
		    call sprintf (Memc[key], SZ_FNAME, "cd%d_%d")
			call pargi (i)
			call pargi (j)
		    #cd1[1] = clgpsetd (pp, Memc[key])
		    #if (!IS_INDEFD(cd1[1]))
		    #	cd[j+3*(i-1)] = cd1[1]
		    call clgpseta  (pp, Memc[key], Memc[str], SZ_LINE)
		    k = nowhite (Memc[str], Memc[str], SZ_LINE)
		    k = 1
		    if (streq (Memc[str], "INDEF"))
			cd[j+3*(i-1)] = cd[j+3*(i-1)]
		    else if (streq (Memc[str], "+INDEF"))
			cd[j+3*(i-1)] = abs(cd[j+3*(i-1)])
		    else if (streq (Memc[str], "-INDEF"))
			cd[j+3*(i-1)] = -abs(cd[j+3*(i-1)])
		    else if (ctod (Memc[str], k, cd1[1]) > 0)
			cd[j+3*(i-1)] = cd1[1]
		}
	    }

	    # If there are square pixels.
	    call clgpseta  (pp, "square", Memc[str], SZ_LINE)
	    if (Memc[str] != EOS) {
	        call sscan (Memc[str])
		call gargi (axis[1])
		call gargi (axis[2])
		switch (nscan()) {
		case 1:
		    call error (1, "Bad SQUARE parameter")
		case 2:
		    x[1] = sqrt(cd[(axis[1]-1)*mwdim+axis[1]]**2 +
			cd[(axis[1]-1)*mwdim+axis[2]]**2)
		    x[2] = sqrt(cd[(axis[2]-1)*mwdim+axis[1]]**2 +
			cd[(axis[2]-1)*mwdim+axis[2]]**2)
		    if (x[1] > x[2]) {
			scale = x[2] / x[1]
			do i = 1, mwdim
			    cd[(axis[1]-1)*mwdim+i] = scale *
				cd[(axis[1]-1)*mwdim+i]
			r[1] = r[1] / scale
		    } else {
			scale = x[1] / x[2]
			do i = 1, mwdim
			    cd[(axis[2]-1)*mwdim+i] = scale *
				cd[(axis[2]-1)*mwdim+i]
			r[2] = r[2] / scale
		    }
		}
	    }
	    
	    # Update the world terms and set the transformation to
	    # pixels for determining the size of the output image.

	    call mw_swtermd (mw, r, w, cd, mwdim)
	    ctwl = mw_sctran (mw, "world", "logical", 0)

	    call mw_close (inmw)
	}

	# Set output dimensionality.
	omwdim = mw_stati (mw, MW_NPHYSDIM)

	# Find the size of the output that includes the desired input data.
	# Also check the input dimensionality matches the ref dimensionality.
        call amovki (MAX_INT, minlen, 3)
	call amovki (-MAX_INT, maxlen, 3)
	do nim = 1, imtlen (images) {

	    if (stack == YES)
		call gfw_open (gf, images, nim, pset, 3, nim)
	    else
		# This is where we could allow the reference dimensionality
		# to be greater than the input.
		call gfw_open (gf, images, nim, pset, INDEFI, 1)
	    mwdim = GF_MWDIM(gf)
	    inmw = GF_MW(gf)

	    if (mwdim != omwdim)
	        call error (3,
		    "Reference WCS dimensionality different than input")

	    # Accumulate limits by sampling the edges of the data.
	    i = imtrgetim (images, nim, Memc[image], SZ_FNAME)
	    tmp = immap (Memc[image], READ_ONLY, 0); im = tmp
	    switch (mwdim) {
	    case 2:
		cmin[1] = clgpsetd (pp, "cmin1")
		if (IS_INDEFD(cmin[1]))
		    cmin[1] = -MAX_DOUBLE
		cmax[1] = clgpsetd (pp, "cmax1")
		if (IS_INDEFD(cmax[1]))
		    cmax[1] = MAX_DOUBLE
		cmin[2] = clgpsetd (pp, "cmin2")
		if (IS_INDEFD(cmin[2]))
		    cmin[2] = -MAX_DOUBLE
		cmax[2] = clgpsetd (pp, "cmax2")
		if (IS_INDEFD(cmax[2]))
		    cmax[2] = MAX_DOUBLE

		do l = 1, mwdim {
		    do k = 1, 2 {
			if (k == 1)
			    x[l] = 0.5
			else
			    x[l] = IM_LEN(im,l) + 0.5
			m = mod (l, mwdim) + 1
			do i = 1, IM_LEN(im,m),
			    max (1,min(10,IM_LEN(im,m)-1)) {
			    x[m] = i
			    call gfwcs_pixel (gf, x, w)
			    w[1] = max (w[1], cmin[1])
			    w[1] = min (w[1], cmax[1])
			    w[2] = max (w[2], cmin[2])
			    w[2] = min (w[2], cmax[2])
			    call mw_ctrand (ctwl, w, r, mwdim)
			    do n = 1, mwdim {
				minlen[n] = min(minlen[n],
				    int(r[n]+0.501))
				maxlen[n] = max(maxlen[n],
				    int(r[n]+0.499))
			    }
			}
		    }
		}
	    case 3:
		cmin[1] = clgpsetd (pp, "cmin1")
		if (IS_INDEFD(cmin[1]))
		    cmin[1] = -MAX_DOUBLE
		cmax[1] = clgpsetd (pp, "cmax1")
		if (IS_INDEFD(cmax[1]))
		    cmax[1] = MAX_DOUBLE
		cmin[2] = clgpsetd (pp, "cmin2")
		if (IS_INDEFD(cmin[2]))
		    cmin[2] = -MAX_DOUBLE
		cmax[2] = clgpsetd (pp, "cmax2")
		if (IS_INDEFD(cmax[2]))
		    cmax[2] = MAX_DOUBLE
		cmin[3] = clgpsetd (pp, "cmin3")
		if (IS_INDEFD(cmin[3]))
		    cmin[3] = -MAX_DOUBLE
		cmax[3] = clgpsetd (pp, "cmax3")
		if (IS_INDEFD(cmax[3]))
		    cmax[3] = MAX_DOUBLE

		do l = 1, mwdim {
		    do k = 1, 2 {
			if (k == 1)
			    x[l] = 0.5
			else
			    x[l] = IM_LEN(im,l) + 0.5
			do j = 1, 2 {
			    m = mod (l, mwdim) + 1
			    if (j == 1)
				x[m] = 0.5
			    else
				x[m] = IM_LEN(im,m) + 0.5
			    m = mod (l+1, mwdim) + 1
			    do i = 1, IM_LEN(im,m),
				max (1,min(10,IM_LEN(im,m)-1)) {
				x[m] = i
				call gfwcs_pixel (gf, x, w)
				w[1] = max (w[1], cmin[1])
				w[1] = min (w[1], cmax[1])
				w[2] = max (w[2], cmin[2])
				w[2] = min (w[2], cmax[2])
				w[3] = max (w[3], cmin[3])
				w[3] = min (w[3], cmax[3])
				call mw_ctrand (ctwl, w, r, mwdim)
				do n = 1, mwdim {
				    if (r[n] < 0) {
					minlen[n] = min(minlen[n],
					    int(r[n]+0.499))
					maxlen[n] = max(maxlen[n],
					    int(r[n]+0.501))
				    } else {
					minlen[n] = min(minlen[n],
					    int(r[n]+0.501))
					maxlen[n] = max(maxlen[n],
					    int(r[n]+0.499))
				    }
				}
			    }
			}
		    }
		}
	    }
	    call imunmap (im)

	    call gfwcs_close (gf)
	}

	# Update the image size and reference pixel.
	call amovki (1, imlen, 3)
	call mw_gwtermd (mw, r, w, cd, omwdim)
	do i = 1, mwdim {
	    imlen[i] = maxlen[i] - minlen[i] + 1
	    r[i] = r[i] - minlen[i] + 1
	}
	call mw_swtermd (mw, r, w, cd, omwdim)
	call mw_ctfree (ctwl)

	call clcpset (pp)
	call sfree (sp)
end


# GFWCS_SF -- Get FITCOORDS surface fits.

procedure gfwcs_sf (gf, im)			# PRIVATE

pointer	gf			#I GF object
pointer	im			#I Image pointer

int	i, j, rec, ncoeffs, nsf[2]
pointer	dt, un1, sf1, coeffs, un[2], sf[2]
pointer	sp, name, db, fname, key

bool	un_compare()
int	dtlocate(), dtgeti()
pointer	dtmap1(), un_open()
errchk	dtmap1

begin
	call smark (sp)
	call salloc (name, SZ_FNAME, TY_CHAR)
	call salloc (db, SZ_FNAME, TY_CHAR)
	call salloc (fname, SZ_FNAME, TY_CHAR)
	call salloc (key, SZ_FNAME, TY_CHAR)

	un[1] = NULL; un[2] = NULL
	nsf[1] = 0; nsf[2] = 0

	iferr (call imgstr (im, "FCDB", Memc[db], SZ_FNAME))
	    Memc[db] = EOS
	do i = 1, ARB {
	    call sprintf (Memc[key], SZ_FNAME, "FCFIT%d")
	        call pargi (i)
	    iferr (call imgstr (im, Memc[key], Memc[name], SZ_FNAME))
	        break
	    call sprintf (Memc[key], SZ_FNAME, "FCFILE%d")
	        call pargi (i)
	    iferr (call imgstr (im, Memc[key], Memc[fname], SZ_FNAME)) {
	        call sprintf (Memc[fname], SZ_FNAME, "fc%s")
		    call pargstr (Memc[name])
	    }

	    dt = dtmap1 (Memc[db], Memc[fname], READ_ONLY)
	    rec = dtlocate (dt, Memc[name])

	    j = dtgeti (dt, rec, "axis")
	    ifnoerr (call dtgstr (dt, rec, "units", Memc[name], SZ_FNAME))
		un1 = un_open (Memc[name])
	    else
	    	un1 = NULL
	    ncoeffs = dtgeti (dt, rec, "surface")
	    call malloc (coeffs, ncoeffs, TY_DOUBLE)
	    call dtgad (dt, rec, "surface", Memd[coeffs], ncoeffs, ncoeffs)
	    call dgsrestore (sf1, Memd[coeffs])
	    call mfree (coeffs, TY_DOUBLE)
	    call dtunmap (dt)

	    if (un1 != NULL) {
	        if (un[j] == NULL)
		    un[j] = un1
		else if (un_compare (un1, un[j]))
		    call un_close (un1)
		else {
		    call un_close (un1)
		    call un_close (un[j])
		    call error (1, "FITCOORDS units disagree")
		}
	    }

	    if (sf1 != NULL) {
		if (nsf[j] == 0)
		    call malloc (sf[j], 5, TY_POINTER)
		else if (mod (nsf[j],5) == 0)
		    call realloc (sf[j], nsf[j]+5, TY_POINTER)
		Memi[sf[j]+nsf[j]] = sf1
		nsf[j] = nsf[j] + 1
	    }
	}

	do j = 1, 2 {
	    GF_UN(gf,j) = un[j]
	    GF_SF(gf,j) = sf[j]
	    GF_NSF(gf,j) = nsf[j]
	}

	call sfree (sp)
end


# GFWCS_GWTERM -- Get world term in the desired output axes orientation.
#
# This will use the MWCS if there are no distortions otherwise it will
# numerically compute them based on evaluating the distorted WCS.

procedure gfwcs_gwtermd (gf, r, w, cd)		# PRIVATE

pointer	gf			# GF structure
double	r[ARB]			# Reference pixel coordinate
double	w[ARB]			# Reference world coordinate
double	cd[ARB]			# Scale matrix

int	i, j, mwdim, axes[3]
double	dval, r1[3], r2[3], w1[3], cd1[9], ltm[9], iltm[9]
double	dgsgetd()
pointer	pp, mw

bool	clgpsetb()

begin
	pp = GF_PP(gf)
	mw = GF_MW(gf)
	mwdim = GF_MWDIM(gf)
	call amovi (GF_AXES(gf,1), axes, 3)

	if (GF_NSF(gf,1) == 0 && GF_NSF(gf,2) == 0) {
	    # When there are no distortion functions we use
	    # the WCS but rearrange the axes to the standard
	    # order.

	    # Remove any physical transformation.
	    call mw_gltermd (mw, ltm, r1, mwdim)
	    call mw_invertd (ltm, iltm, mwdim)
	    call mw_gwtermd (mw, r, w1, cd, mwdim)
	    call mw_vmuld (ltm, r, w, mwdim) 
	    call aaddd (w, r1, r, mwdim)
	    call mw_mmuld (cd, iltm, cd1, mwdim)

	    # Map the axes.
	    do i = 1, mwdim {
		r[i] = r1[axes[i]]
		w[i] = w1[axes[i]]
		do j = 1, mwdim
		    cd[(i-1)*mwdim+j] = cd1[(axes[i]-1)*mwdim+j]
	    }

	    # Set the logical axes to correspond most closely to the world axes.
	    if (mwdim == 3) {
		do j = 1, 3 {
		    if (cd[2*mwdim+j] != 0.)
			break
		}
		do i = 1, 3 {
		    dval = cd[(i-1)*mwdim+3]
		    cd[(i-1)*mwdim+3] = cd[(i-1)*mwdim+j]
		    cd[(i-1)*mwdim+j] = dval
		}
	    }
	    if (abs(cd[2]) > abs(cd[1])) {
		do i = 1, 2 {
		    dval = cd[(i-1)*mwdim+1]
		    cd[(i-1)*mwdim+1] = cd[(i-1)*mwdim+2]
		    cd[(i-1)*mwdim+2] = dval
		}
	    }
	} else {
	    # Set reference point at first pixel and use
	    # numerical derivatives.  We can't just use the
	    # input tangent point because the distortion
	    # evaluation is only defined within the image.

	    call amovkd (1D0, r, mwdim)
	    call amovkd (2D0, r1, mwdim)
	    r[1] = -MAX_DOUBLE; r[2] = -MAX_DOUBLE
	    r1[1] = MAX_DOUBLE; r1[2] = MAX_DOUBLE
	    do i = 1, 2 {
		do j = 1, GF_NSF(gf,i) {
		    r[1] = max (r[1], dgsgetd (Memi[GF_SF(gf,i)+j-1], GSXMIN))
		    r[2] = max (r[2], dgsgetd (Memi[GF_SF(gf,i)+j-1], GSYMIN))
		    r1[1] = min (r1[1], dgsgetd (Memi[GF_SF(gf,i)+j-1], GSXMAX))
		    r1[2] = min (r1[2], dgsgetd (Memi[GF_SF(gf,i)+j-1], GSYMAX))
		}
	    }
	    
	    call gfwcs_pixel (gf, r, w)
	    call amovd (r, r2, mwdim)
	    do i = 1, mwdim {
		r2[axes[i]] = r1[axes[i]]
		call gfwcs_pixel (gf, r2, w1)
		r2[axes[i]] = r[axes[i]]
		do j = 1, mwdim
		    cd[(j-1)*mwdim+i] = (w1[j] - w[j]) /
		        (r1[axes[i]] - r[axes[i]])
	    }
	}

	if (!clgpsetb (pp, "rotate")) {
	    do i = 1, mwdim {
		do j = 1, mwdim {
		    if (i != j)
			cd[(i-1)*mwdim+j] = 0.
		}
	    }
	}
end
