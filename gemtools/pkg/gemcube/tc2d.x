# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.

include	<math.h>
include	"transcube.h"

define	DEBUG	true


# TC_2DNEAREST -- Set weights for nearest.

procedure tc_2dnearest (wts, no, ns)

pointer	wts[2]			#O Weight look up array
int	no[2]			#0 No. of output pixels covered by input pixel
int	ns[2]			#0 Number of output subpixel centers sampled

begin
	no[1] = 1; no[2] = 1
	ns[1] = 1; ns[2] = 1
	call malloc (wts[1], 1, TY_POINTER)
	Memr[wts[1]] = 1.
end


# TC_2DWT -- Compute weights for rectangular input pixels.
#
# The weights are an array of interpolator pointers for each overlap pixel.
# The interpolator is indexed by subpixel offset centers.  A special case
# is if no or ns are 1 and then the returned pointer is for a single
# real weight value.

procedure tc_2dwt (wttype, dw, pa, cd, interp, no, ns)

int	wttype			#I Weighting type
double	dw[2]			#I Pixel size (world coordinates)
double	pa			#I Position angle
double	cd[2,2]			#I Transformation matrix (p -> w)
pointer	interp			#O Interpolation or weight pointers
int	no[2]			#O No. of output pixels covered by input pixel
int	ns[2]			#O No. of output subpixel centers sampled

int	i, j, n, no1, no2
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
	    call tc_2ddriz (dw, cd, pa, INDEFI, INDEFI, wts, no, ns)

	    n = no[1] * no[2]
	    call malloc (interp, n, TY_POINTER)
	    if (n == 1) {
		Memr[interp] = 1.
		return
	    }
	    n = ns[1] * ns[2]

	    # For each offset pixel compute the weights and fit interpolator.
	    no1 = (no[1] - 1) / 2
	    no2 = (no[2] - 1) / 2
	    ptr = interp
	    do j = -no2, no2 {
		do i = -no1, no1 {
		    call tc_2ddriz (dw, cd, pa, i, j, wts, no, ns)
		    if (n == 1)
			Memr[ptr] = Memr[wts]
		    else
			call tc_interp (Memr[wts], ns[1], ns[2], Memi[ptr])
		    ptr = ptr + 1
		}
	    }
#	case WT_LIN:
#	    # Allocate memory, determine no and ns, and compute weights.
#	    call tc_2dlin (dw, wts, no, ns)
#
#	    n = no[1] * no[2]
#	    call malloc (interp, n, TY_POINTER)
#	    if (n == 1) {
#		Memr[interp] = 1.
#		return
#	    }
#	    n = ns[1] * ns[2]
#
#	    no1 = (no[1] - 1) / 2
#	    no2 = (no[2] - 1) / 2
#	    wt = wts; ptr = interp
#	    do j = -no2, no2 {
#		do i = -no1, no1 {
#		    if (n == 1)
#			Memr[ptr] = Memr[wt]
#		    else
#			call tc_interp (Memr[wt], ns[1], ns[2], Memi[ptr])
#		    wt = wt + ns; ptr = ptr + 1
#		}
#	    }
	}

	call mfree (wts, TY_REAL)
end


# TC_2DDRIZ -- Compute overlap weights for a range of subpixel centers
# at a given output pixel.

procedure tc_2ddriz (dw, cd, pa, i, j, wts, no, ns)

double	dw[2]			#I Input pixel size (world coordinates)
double	cd[2,2]			#I Transformation matrix (p -> w)
double	pa			#I Position angle
int	i, j			#I Output overlapped pixel
pointer	wts			#O Weights
int	no[2]			#O No. of overlap pixels
int	ns[2]			#O No. of sample centers

int	k, l, m, n, nx, ny, nx2, ny2, sum
real	stepx, stepy, subpix[2], subpix2[2], val, dx, dy, dx1, dy1
double	c, s, ratio, cdpw[2,2], cdwp[2,2], dw2[2], p[2], w[2]
pointer	wt

begin
	# Initialize if the overlap pixel is INDEF.
	if (IS_INDEFI(i)) {

	    c = cos (DEGTORAD(pa))
	    s = sin (DEGTORAD(pa))

	    # Determine number of overlap pixels.
	    call mw_invertd (cd, cdwp, 2)

	    # Footprint of input pixel on output pixels.
	    dw2[1] =  dw[1]/2 * c + dw[2]/2 * s
	    dw2[2] = -dw[1]/2 * s + dw[2]/2 * c
	    call mw_vmuld (cdwp, dw2, w, 2)

	    dw2[1] = -dw[1]/2 * c + dw[2]/2 * s
	    dw2[2] =  dw[1]/2 * s + dw[2]/2 * c
	    call mw_vmuld (cdwp, dw2, p, 2)

	    p[1] = max (abs(w[1]), abs (p[1]))
	    p[2] = max (abs(w[2]), abs (p[2]))

	    no[1] = 2 * int (p[1] + 1.99) - 1
	    no[2] = 2 * int (p[2] + 1.99) - 1
	    if (no[1] == 1 && no[2] == 1) {
		ns[1] = 1; ns[2] = 1
	    } else {
		ns[1] = NSAMPLEX; ns[2] = NSAMPLEY
	    }

	    call malloc (wts, ns[1]*ns[2], TY_REAL)

	    # Now set the transformations between input and output frames.
	    # If the input pixel PA in world coordinates matches the
	    # the PA of the output world coordinate system then will
	    # result in aligned input and output pixels.

	    cdpw[1,1] = cd[1,1] * c - cd[1,2] * s
	    cdpw[1,2] = cd[1,1] * s + cd[1,2] * c
	    cdpw[2,1] = cd[2,1] * c - cd[2,2] * s
	    cdpw[2,2] = cd[2,1] * s + cd[2,2] * c
	    call mw_invertd (cdpw, cdwp, 2)

	    # Determine which is bigger, the input or output pixel.
	    ratio = abs (cdpw[1,1]*cdpw[2,2] + cdpw[1,2]*cdpw[2,1]) /
	       (dw[1]*dw[2])

	    return
	}

	dw2[1] = dw[1] / 2; dw2[2] = dw[2] / 2
	nx = ns[1]; ny = ns[2]
	nx2 = (nx - 1) / 2 + 1; ny2 = (ny - 1) / 2 + 1
	stepx = 1. / max (1, nx - 1); stepy = 1. / max (1, ny - 1)
	
	if (ratio <= 1.) {
	    subpix[1] = 1. / NSUBPIX; subpix[2] = 1. / NSUBPIX
	    subpix2[1] = subpix[1] / 2; subpix2[2] = subpix[2] / 2
	    val = subpix[1] * subpix[2]
	    
	    wt = wts
	    do l = 1, ny {
		dy = j - (l - ny2) * stepy
		dy1 = dy - 0.5 - subpix2[2]
		do k = 1, nx {
		    dx = i - (k - nx2) * stepx
		    dx1 = dx - 0.5 - subpix2[1]

		    # Check if the output pixel is completely in the input pixel.
		    sum = 0.
		    do n = 1, NSUBPIX, NSUBPIX-1 {
			p[2] = dy1 + n * subpix[2]
			do m = 1, NSUBPIX, NSUBPIX-1 {
			    p[1] = dx1 + m * subpix[1]
			    call mw_vmuld (cdpw, p, w, 2)
			    if (abs(w[1]) <= dw2[1] && abs(w[2]) <= dw2[2])
				sum = sum + 1
			}
		    }
		    if (sum == 4) {
			Memr[wt] = 1
			wt = wt + 1
			next
		    }

		    # Subsample the pixel.
		    do n = 2, NSUBPIX-1 {
			p[2] = dy1 + n * subpix[2]
			do m = 2, NSUBPIX-1 {
			    p[1] = dx1 + m * subpix[1]
			    call mw_vmuld (cdpw, p, w, 2)
			    if (abs(w[1]) <= dw2[1] && abs(w[2]) <= dw2[2])
				sum = sum + 1
			}
		    }
		    Memr[wt] = sum * val
		    wt = wt + 1
		}
	    }


	} else {
	    subpix[1] = 2 * dw2[1] / NSUBPIX; subpix[2] = 2 * dw2[2] / NSUBPIX
	    subpix2[1] = subpix[1] / 2; subpix2[2] = subpix[2] / 2
	    val = (subpix[1] * subpix[2]) / (4 * dw2[1] * dw2[2])
	    
	    wt = wts
	    do l = 1, ny {
		dy = j + (l - ny2) * stepy
		dy1 = -dw2[2] - subpix2[2]
		do k = 1, nx {
		    dx = i + (k - nx2) * stepx
		    dx1 =  -dw2[1] - subpix2[1]

		    # Check if the output pixel is completely in the input pixel.
		    sum = 0.
		    do n = 1, NSUBPIX, NSUBPIX-1 {
			w[2] = dy1 + n * subpix[2]
			do m = 1, NSUBPIX, NSUBPIX-1 {
			    w[1] = dx1 + m * subpix[1]
			    call mw_vmuld (cdwp, w, p, 2)
			    if (abs(p[1]-dx) <= 0.5 && abs(p[2]-dy) <= 0.5)
				sum = sum + 1
			}
		    }
		    if (sum == 4) {
			Memr[wt] = 1
			wt = wt + 1
			next
		    }

		    # Subsample the pixel.
		    do n = 2, NSUBPIX-1 {
			w[2] = dy1 + n * subpix[2]
			do m = 2, NSUBPIX-1 {
			    w[1] = dx1 + m * subpix[1]
			    call mw_vmuld (cdwp, w, p, 2)
			    if (abs(p[1]-dx) <= 0.5 && abs(p[2]-dy) <= 0.5)
				sum = sum + 1
			}
		    }
		    Memr[wt] = sum * val
		    wt = wt + 1
		}
	    }
	}
end


# TC_2DDP -- Dump weight look up array.

procedure tc_2ddp (dw, wts, no, ns)

double	dw[2]			#I Offset from center of pixel
pointer	wts			#I Weight look up array
int	no[2]			#I No. of output pixels covered by input pixel
int	ns[2]			#I Number of output subpixel centers sampled

int	i, j
int	nsfit, no1, no2
real	wt, sum
real	dx, dy, msieval()
pointer	msi

begin
	call eprintf ("\n(%.3f, %.3f)\n")
	    call pargd (dw[1])
	    call pargd (dw[2])

	dx = (dw[1] + 0.5) * (ns[1] - 1) + 1
	dy = (dw[2] + 0.5) * (ns[2] - 1) + 1
	nsfit = ns[1] * ns[2]

	no1 = (no[1] - 1) / 2; no2 = (no[2] - 1) / 2
	msi = wts
	sum = 0.
	do j = -no2, no2 {
	    do i = -no1, no1 {
		if (nsfit > 1)
		    wt = msieval(Memi[msi], dx, dy)
		else
		    wt = Memr[msi]
		#call eprintf (" %5.3f")
		    #call pargr (wt)
		call eprintf (" %2d")
		    call pargr (wt*10)
		sum = sum + wt
		msi = msi + 1
	    }
	    call eprintf ("\n")
	}
	call eprintf ("%7.4f\n")
	    call pargr (sum)
end


## TC_2DLINEAR -- Compute overlap weights for a range of subpixel centers
## at a given output pixel.  In this routine the rectangular pixels are
## aligned so that the overlaps are simply linearly related to the offsets.
#
#procedure tc_2dlinear (dw, i, j, wt, nx, ny)
#
#double	dw[2]				#I Input pixel size (output pixels)
#int	i, j				#I Output overlapped pixel
#real	wt[nx, ny]			#O Weights
#int	nx, ny				#I Number of subpixel samples
#
#int	k, l, nx2, ny2
#real	stepx, stepy, dx, dy, wt2
#
#begin
#	nx2 = (nx - 1) / 2 + 1; ny2 = (ny - 1) / 2 + 1
#	stepx = 1. / max (1, nx - 1); stepy = 1. / max (1, ny - 1)
#	
#	do l = 1, ny {
#	    dy = (l - ny2) * stepy
#	    if (j < 0)
#		wt2 = min (1., max (0., dw[2]-dy+j+0.5))
#	    else if (j > 0)
#		wt2 = min (1., max (0., dw[2]+dy-j+0.5))
#	    else
#		wt2 = min (0.5-dy, dw[2]) + min (0.5+dy, dw[2])
#
#	    do k = 1, nx {
#		dx = (k - nx2) * stepx
#		if (i < 0)
#		    wt[k,l] = min (1., max (0., dw[1]-dx+i+0.5)) * wt2
#		else if (i > 0)
#		    wt[k,l] = min (1., max (0., dw[1]+dx-i+0.5)) * wt2
#		else
#		    wt[k,l] = (min (0.5-dx, dw[1]) + min (0.5+dx, dw[1])) * wt2
#	    }
#	}
#end
