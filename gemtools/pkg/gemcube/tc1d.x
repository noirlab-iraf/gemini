# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.

include	"transcube.h"

define	DEBUG	true


# TC_1DNEAREST -- Set weights for nearest.

procedure tc_1dnearest (wts, no, ns)

pointer	wts			#O Weight look up array
int	no			#0 No. of output pixels covered by input pixel
int	ns			#0 Number of output subpixel centers sampled

begin
	no = 1; ns = 1
	call malloc (wts, 1, TY_POINTER)
	Memr[wts] = 1.
end


# TC_1DWT -- Compute weights for rectangular input pixels.
#
# The weights are an array of interpolator pointers for each overlap pixel.
# The interpolator is indexed by subpixel offset centers.  A special case
# is if no or ns are 1 and then the returned pointer is for a single
# real weight value.

procedure tc_1dwt (wttype, dw, interp, no, ns)

int	wttype			#I Weighting type
double	dw			#I Pixel size (units of output pixels)
pointer	interp			#O Interpolation or weight pointers
int	no			#O No. of output pixels covered by input pixel
int	ns			#O No. of output subpixel centers sampled

int	i, no1
pointer	wts, ptr

begin
	# Loop through the overlap pixels, compute weights at subpixel
	# centers, and fit an interpolation function.  Note that for
	# drizzle weights we can do each offset pixel independently to
	# save on working memory.  For interpolation weights all the
	# weights are computed at once and then the sample interpolators
	# are fit.

	switch (wttype) {
	case WT_DRIZ:
	    # Allocate memory and determine number of overlap and offset pixels.
	    call tc_1ddriz (dw, INDEFI, wts, no, ns)

	    call malloc (interp, no, TY_POINTER)
	    if (no == 1) {
		Memr[interp] = 1.
		return
	    }

	    # For each offset pixel compute the weights and fit interpolator.
	    no1 = (no - 1) / 2
	    ptr = interp
	    do i = -no1, no1 {
		call tc_1ddriz (dw, i, wts, no, ns)
		if (ns == 1)
		    Memr[ptr] = Memr[wts]
		else
		    call tc_interp (Memr[wts], ns, 1, Memi[ptr])
		ptr = ptr + 1
	    }
#	case WT_LIN:
#	    # Allocate memory, determine no and ns, and compute weights.
#	    call tc_1dlin (dw, wts, no, ns)
#
#	    call malloc (interp, no, TY_POINTER)
#	    if (no == 1) {
#		Memr[interp] = 1.
#		return
#	    }
#
#	    no1 = (no - 1) / 2
#	    wt = wts; ptr = interp
#	    do i = -no1, no1 {
#		if (ns == 1)
#		    Memr[ptr] = Memr[wt]
#		else
#		    call tc_interp (Memr[wt], ns, 1, Memi[ptr])
#		wt = wt + ns; ptr = ptr + 1
#	    }
	}

	call mfree (wts, TY_REAL)
end


# TC_1DDRIZ -- Compute overlap weights for a range of subpixel centers
# at a given output pixel.  In this routine the rectangular pixels are
# aligned so that the overlaps are simply linearly related to the offsets.

procedure tc_1ddriz (dw, i, wts, no, ns)

double	dw			#I Input pixel size (units of output pixels)
int	i			#I Output overlapped pixel
pointer	wts			#U Weights
int	no			#U Number of overlap pixels
int	ns			#U Number of subpixel samples

int	j, nx2
real	stepx, dx
double	dw2
pointer	wt

begin
	# Initialize if the overlap pixel is INDEF.
	if (IS_INDEFI(i)) {
	    no = 2 * int (dw/2 + 1.99) - 1
	    if (no == 1)
		ns = 1
	    else
		ns = NSAMPLEZ

	    call malloc (wts, ns, TY_REAL)
	    return
	}

	# If there is only one overlap the weight is always 1.
	if (no == 1) {
	    Memr[wts] = 1.
	    return
	}

	# Compute drizzle weights for pixel.
	dw2 = dw / 2
	nx2 = (ns - 1) / 2 + 1
	stepx = 1. / max (1, ns - 1)
	
	wt = wts
	do j = 1, ns {
	    dx = (j - nx2) * stepx
	    if (i < 0)
		Memr[wt] = min (1., max (0., dw2-dx+i+0.5))
	    else if (i > 0)
		Memr[wt] = min (1., max (0., dw2+dx-i+0.5))
	    else
		Memr[wt] = min (0.5-dx, dw2) + min (0.5+dx, dw2)
	    wt = wt + 1
	}
end


# TC_1DDP -- Dump weight look up array.

procedure tc_1ddp (dw, wts, no, ns)

double	dw			#I Offset from center of pixel
pointer	wts			#I Weight look up array
int	no			#I No. of output pixels covered by input pixel
int	ns			#I Number of output subpixel centers sampled

int	i
int	no1
real	wt, sum
real	dx, asieval()
pointer	asi

begin
	call eprintf ("\n(%.3f)\n")
	    call pargd (dw)

	dx = (dw + 0.5) * (ns - 1) + 1

	no1 = (no - 1) / 2
	asi = wts
	sum = 0.
	do i = -no1, no1 {
	    if (ns > 1)
		wt = asieval(Memi[asi], dx)
	    else
		wt = Memr[asi]
	    call eprintf (" %5.3f")
		call pargr (wt)
	    sum = sum + wt
	    asi = asi + 1
	}
	call eprintf ("\n%7.4f\n")
	    call pargr (sum)
end
