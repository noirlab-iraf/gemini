# Copyright(c) 2006-2009 Association of Universities for Research in Astronomy, Inc.

include	<math.h>
include	<math/iminterp.h>
include	<mwset.h>
include	"transcube.h"

define	DEBUG	false


# TC_WEIGHTS -- Set the input pixel to output overlap pixels weights.
#
# This involves parsing the shape description and pre-computing the
# weights.  The weights may interpolators as a function of offset
# between the input pixel center and the output pixel center.
# The weights are based on the type of weight chosen.

procedure tc_weights (mw, shape, axmap, wttype, drizscale, wts, no, ns)

pointer	mw			#I Output MWCS
char	shape[ARB]		#I Drop shape description
int	axmap[3]		#I Axis map
int	wttype			#I Weighting type
real	drizscale[3]		#I Drop scale factors
pointer	wts[3]			#O Weight look up array
int	no[3]			#O No. of output pixels covered by input pixel
int	ns[3]			#O Number of output subpixel centers covered

double	dw[3], pa
pointer	sp, str

int	strdic(), nscan()

define	err_	99

begin
	# If using nearest weighting or the drizzle weighting scales are all
	# zero we use nearest weights regardless of the shape of the pixel
	# footprint.

	if (wttype == WT_NEAREST || (wttype == WT_DRIZ &&
	    drizscale[1] <= 0. && drizscale[2] <= 0. && drizscale[3] <= 0.)) {
	    call tc_2dnearest (wts[1], no[1], ns[1])
	    call tc_1dnearest (wts[3], no[3], ns[3])
	    return
	}

	call smark (sp)
	call salloc (str, SZ_FNAME, TY_CHAR)

	# Parse shape string.
	call sscan (shape)
	call gargwrd (Memc[str], SZ_FNAME)
	switch (strdic (Memc[str], Memc[str], SZ_FNAME, SHAPES)) {
	case 1:		# Rectangle
	    call gargd (dw[1])
	    call gargd (dw[2])
	    call gargd (dw[3])
	    call gargd (pa)
	    if (nscan() < 5)
	        goto err_
	    if (drizscale[1] <= 0. && drizscale[2] <= 0.) {
		dw[1] = max (0., drizscale[1]) * dw[1]
		dw[2] = max (0., drizscale[2]) * dw[2]
	    } else {
		dw[1] = max (0.01, drizscale[1]) * dw[1]
		dw[2] = max (0.01, drizscale[2]) * dw[2]
	    }
	    dw[3] = max (0., drizscale[3]) * dw[3]
	    call tc_rectangle (mw, dw, pa, axmap, wttype, wts, no, ns)
	    if (DEBUG) {
		for (dw[2]=-0.5; dw[2]<=0.501; dw[2]=dw[2]+0.5)
		    for (dw[1]=-0.5; dw[1]<=0.501; dw[1]=dw[1]+0.5)
			call tc_2ddp (dw[1], wts[1], no[1], ns[1])
	    }
	case 2:
	    call gargd (dw[1])
	    call gargd (dw[3])
	    if (nscan() < 3)
	        goto err_
	case 3:
	    call gargd (dw[1])
	    call gargd (dw[3])
	    call gargd (pa)
	    if (nscan() < 4)
	        goto err_
	default:
err_
	    call error (1, "Error in shape string")
	}

	call sfree (sp)
end


# TC_FREE -- Free memory.

procedure tc_free (wts, no, ns)

pointer	wts[3]			#U Weight look up array
int	no[3]			#I No. of output pixels covered by input pixel
int	ns[3]			#I Number of output subpixel centers sampled

int	i

begin
	if (ns[1]*ns[2] > 1) {
	    do i = 1, no[1]*no[2]
	        call msifree (Memi[wts[1]+i-1])
	}
	call mfree (wts[1], TY_POINTER)

	if (ns[3] > 1) {
	    do i = 1, no[3]
	        call asifree (Memi[wts[3]+i-1])
	}
	call mfree (wts[3], TY_POINTER)
end


# TC_RECTANGLE -- Compute weight look up array for rectangular input pixels.
#
# This routine assumes that *input* z axis is independent of the x/y axes.
# It also handles the output 2D WCS.

procedure tc_rectangle (mw, dw, pa, axmap, wttype, wts, no, ns)

pointer	mw			#I Output MWCS
double	dw[3]			#I Input pixel size in world coordinates
double	pa			#I Position angle
int	axmap[3]		#I Axis map
int	wttype			#I Weighting type
pointer	wts[3]			#O Weight look up array
int	no[3]			#O No. of output pixels covered by input pixel
int	ns[3]			#O Number of output subpixel centers covered

double	r[3], w[3], cd[3,3], cd2[2,2]

int	mw_stati()

begin
	# Get the CD matrix.
	switch (mw_stati (mw, MW_NPHYSDIM)) {
	case 2:
	    call mw_gwtermd (mw, r, w, cd2, 2)
	    cd[1,1] = cd2[1,1]; cd[1,2] = cd2[1,2]; cd[1,3] = 0.
	    cd[2,1] = cd2[2,1]; cd[2,2] = cd2[2,2]; cd[2,3] = 0.
	    cd[3,1] = 0.;       cd[3,2] = 0.;       cd[3,3] = 1.
	case 3:
	    call mw_gwtermd (mw, r, w, cd, 3)
	default:
	    call error (1, "WCS dimensionality not supported")
	}

	# Map the axes.
	#w[1] = dw[axmap[1]]; w[2] = dw[axmap[2]] w[3] = dw[axmap[3]]
	w[1] = dw[1]; w[2] = dw[2]; w[3] = dw[3]
	cd2[1,1] = cd[axmap[1],axmap[1]]; cd2[1,2] = cd[axmap[1],axmap[2]]
	cd2[2,1] = cd[axmap[2],axmap[1]]; cd2[2,2] = cd[axmap[2],axmap[2]]
	cd2[3,3] = cd[axmap[3],axmap[3]]

	# Compute the weights.
	if (w[1] <= 0. && w[2] <= 0.)
	    call tc_2dnearest (wts[1], no[1], ns[1])
	else
	    call tc_2dwt (wttype, w[1], pa, cd2, wts[1], no[1], ns[1])
	if (w[3] <= 0.)
	    call tc_1dnearest (wts[3], no[3], ns[3])
	else
	    call tc_1dwt (wttype, w[3]/cd[axmap[3],axmap[3]],
	        wts[3], no[3], ns[3])
end


# TC_INTERP -- Convert offset subpixel array to interpolation function.

procedure tc_interp (wts, ns1, ns2, interp)

real	wts[ARB]		#I Weight look up array
int	ns1, ns2		#I Number of output subpixel centers sampled
pointer	interp			#O Pointer to iterpolation functions

begin
	# For a single overlap pixel the weight is the interpolator.
	if (ns1 * ns2 == 1) {
	    interp = wts[1]
	    return
	}

	# Fit interpolator.
	if (ns2 == 1) {
	    call asiinit (interp, II_LINEAR)
	    call asifit (interp, wts, ns1)
	} else {
	    call msiinit (interp, II_BILINEAR)
	    call msifit (interp, wts, ns1, ns2, ns1)
	}
end
