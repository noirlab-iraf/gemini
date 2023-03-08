# Copyright(c) 2006-2009 Association of Universities for Research in Astronomy, Inc.

include	<imhdr.h>
include	<mwset.h>
include	<mach.h>
include	<plset.h>
include	"transcube.h"

define	DEBUG		false

define	TOL1		1E-3		# Tolerance on equality of pixel size
define	TOL2		1		# Tolerance on equality ofposition angle


# TRANSCUBE -- Cube transform procedure.
#
# This is the core task that operates on a single input. 
#
# There are two complexities here.  One is handling I/O since the coordinate
# transformations may be arbitrary.  One optimization is to try and buffer
# as much of the output in memory as possible.  When the output is too
# large then multiple passes are needed to build up the output in subcubes.
#
# The other complexity is adding the input pixel values on to the output
# pixel grid.  The transformation geometry is handled by a geometry function
# interface.  This task only does pixel I/O from the input image and does not
# itself deal with the input WCS.  It deals with the output WCS and
# transformations between world coordinates and pixel coordinate in the
# input image.  The other efficiency is that the pixel weights are precomputed
# assuming the input pixel sizes do not change significantly.
#
# There are various extra features including input and output masks and
# weights, scaling, and log output.
#
# Currently this task is limited to celestial or coupled axes on the
# first two output axes and an independent third coordinate (such as
# wavelength) on the third axes.

procedure transcube (input, nin, output, masks, weights, bpm, scale, wt,
	wcsreference, wttype, drizscale, blank, geofunc, memalloc,
	logfiles, nlogs) 

int	input				#I List of input images
int	nin				#I Index of input image
char	output[ARB]			#I Output image
char	masks[ARB]			#I Output mask
char	weights[ARB]			#I Output weights
char	bpm[ARB]			#I Input bad pixel
char	scale[ARB]			#I Input scale
char	wt[ARB]				#I Input weight
char	wcsreference[ARB]		#I WCS reference
char	wttype[ARB]			#I Weighting type
real	drizscale[3]			#I Drizzle scale factors
real	blank				#I Blank value
char	geofunc[ARB]			#I Geometry function
real	memalloc			#I Memory to alloc in MB
int	logfiles[ARB]			#I Output logfile descriptors
int	nlogs				#I Number of logfile descriptors

int	i, j, k, l, mwdim, ioff[3], mode, ns12, v2, v3
int	wtt, max_npix, npix, blk
int	v[3], imlen[3], buflen[3], vo1[3], vo2[3], nvo[3], nblk[3]
int	axmap[3], ioutpix[3], ioutmin[3], ioutmax[3]
int	no[3], ns[3], nohalf[3], nshalf[3]
real	s, w, pixval, pixwt, wtz, routpix[3]
double	inpix[3], outpix[3], outworld[3], wtpars[4]
pointer	in, out, pl, bp, outwt, inwt, refmw, outmw, gf, ctwl
pointer	inbuf, bpbuf, outbuf, outwtbuf, inwtbuf, outline, outwtline
pointer	sp, image, bpmfname, wtstr, str, shape, shape1
pointer	wtfname, plbuf, wts[3], asi, msi

bool	streq(), strne(), pl_empty(), pl_linenotempty()
int	imtrgetim(), imtlen()
int	mw_stati(), ctor(), strlen(), strdic(), imaccess()
real	asieval(), msieval(), imgetr()
pointer	immap(), imgl3r(), imgl3s(), imps3r(), impl3r(), impl3s()
pointer	pl_create()
pointer	mw_openim(), mw_sctran()

begin
	call smark (sp)
	call salloc (image, SZ_FNAME, TY_CHAR)
	call salloc (bpmfname, SZ_FNAME, TY_CHAR)
	call salloc (wtstr, SZ_FNAME, TY_CHAR)
	call salloc (str, SZ_LINE, TY_CHAR)
	call salloc (shape, 1024, TY_CHAR)
	call salloc (shape1, 1024, TY_CHAR)

	# First check for valid enumerated parameters.
	wtt = strdic (wttype, Memc[str], SZ_LINE, WTTYPES)
	if (wtt == 0)
	    call error (WT_ERR, "Unknown weighting type")

	# Set memory to alloc for output image and weights.
	max_npix = max (65000., memalloc / 2 / 4 * 1000000)

	# Open input image.
	i = imtrgetim (input, nin, Memc[image], SZ_LINE)
	in = immap (Memc[image], READ_ONLY, 0)

	# Open input mask.
	if (bpm[1] != EOS) {
	    if (bpm[1] == '!')
	        call imgstr (in, bpm[2], Memc[bpmfname], SZ_FNAME)
	    else if (streq (bpm,"BPM"))
	        call imgstr (in, bpm, Memc[bpmfname], SZ_FNAME)
	    else
	        call strcpy (bpm, Memc[bpmfname], SZ_FNAME)
	    bp = immap (Memc[bpmfname], READ_ONLY, 0)
	    do i = 1, IM_NDIM(in) {
		if (IM_LEN(in,i) != IM_LEN(bp,i))
		    call error (1, "Bad pixel mask does not match input")
	    }
	} else
	    bp = NULL

	# Set scale.
	if (scale[1] != EOS) {
	    if (scale[1] == '!')
	        s = imgetr (in, scale[2])
	    else {
	        call sscan (scale)
		call gargr (s)
	    }
	} else
	    s = 1.

	# Set input weights.
	if (wt[1] == EOS) {
	    inwt = NULL
	    w = 1.
	} else {
	    if (wt[1] == '!')
	        call imgstr (in, wt[2], Memc[wtstr], SZ_FNAME)
	    else
	        call strcpy (wt, Memc[wtstr], SZ_FNAME)

	    ifnoerr (inwt = immap (Memc[wtstr], READ_ONLY, 0)) {
		do i = 1, IM_NDIM(in) {
		    if (IM_LEN(in,i) != IM_LEN(inwt,i))
			call error (1,
			   "Input weight image does not match input image")
		}
	    } else {
	       inwt = NULL
	       i = 1
	       if (ctor (Memc[wtstr], i, w) != strlen(Memc[wtstr]))
	           call error (1, "Syntax error in input weights")
	    }
	}

	# Open output images and set the MWCS transformation on the first input.
	if (nin == 1) {
	    if (nlogs > 0) {
		do i = 1, nlogs {
		    call fprintf (logfiles[i],
		        "  Geometry function driver: %s\n")
		        call pargstr (geofunc)
		    call fprintf (logfiles[i], "  Weight type: %s\n")
		        call pargstr (wttype)
		    if (wtt == WT_DRIZ) {
		        call fprintf (logfiles[i], "    Drizzle scale(s): %.5g")
			    call pargr (drizscale[1])
			if (drizscale[2] != drizscale[1] ||
			    drizscale[3] != drizscale[1]) {
			    call fprintf (logfiles[i], " %.5g %.5g\n")
			        call pargr (drizscale[2])
				call pargr (drizscale[3])
			} else
			    call fprintf (logfiles[i], "\n")
		    }
		    call fprintf (logfiles[i], "  Blank value: %.5g\n")
		        call pargr (blank)
		    call fprintf (logfiles[i], "  Memory allocation: %.5g Mb\n")
		        call pargr (memalloc)
		    if (masks[1] != EOS) {
			call fprintf (logfiles[i], "  Output mask: %s\n")
			    call pargstr (masks)
		    }
		    if (weights[1] != EOS) {
			call fprintf (logfiles[i], "  Output weights: %s\n")
			    call pargstr (weights)
		    }
		}
	    }

	    mode = READ_WRITE
	    iferr (out = immap (output, mode, 0)) {
		# Determine output WCS and size.
		if (wcsreference[1] != EOS) {
		    if (nlogs > 0) {
		        do i = 1, nlogs {
			    call fprintf (logfiles[i],
			        "  Output WCS reference: %s\n")
			    call pargstr (wcsreference)
			}
		    }
		    outwt = immap (wcsreference, READ_ONLY, 0)
		    refmw = mw_openim (outwt)
		    call imunmap (outwt)
		    call gf_out (geofunc, input, refmw, outmw, imlen)
		    call mw_close (refmw)
		} else {
		    refmw = NULL
		    call gf_out (geofunc, input, refmw, outmw, imlen)
		}

		# Set new output image.
		mode = NEW_COPY
		out = immap (output, mode, in)
		IM_PIXTYPE(out) = TY_REAL
		IM_NDIM(out) =  mw_stati (outmw, MW_NPHYSDIM)
		call amovi (imlen, IM_LEN(out,1), 3)
		call mw_saveim (outmw, out)

		npix = imlen[1] * imlen[2] * imlen[3]
		if (npix > max_npix) {
		    do k = 1, imlen[3] {
			do j = 1, imlen[2]
			    call aclrr (Memr[impl3r(out,j,k)], imlen[1])
		    }
		    call imunmap (out)
		    out = immap (output, READ_WRITE, 0)
		}
	    } else {
		call amovi (IM_LEN(out,1), imlen, 3)
		outmw = mw_openim (out)
	    }
	    npix = imlen[1] * imlen[2] * imlen[3]


	    # We need to use a temporary weight image if the output
	    # is larger than we can buffer in memory.

	    call malloc (wtfname, SZ_FNAME, TY_CHAR)
	    call strcpy (weights, Memc[wtfname], SZ_FNAME)
	    if (Memc[wtfname] == EOS && npix > max_npix)
		call mktemp ("tmp", Memc[wtfname], SZ_FNAME)

	    if (Memc[wtfname] != EOS) {
		if (mode == NEW_COPY) {
		    outwt = immap (Memc[wtfname], mode, in)
		    IM_PIXTYPE(outwt) = TY_REAL
		    IM_NDIM(outwt) = IM_NDIM(out)
		    call amovi (imlen, IM_LEN(outwt,1), 3)
		    call mw_saveim (outmw, outwt)
		    if (npix > max_npix) {
			if (DEBUG)
			    call eprintf ("Creating empty weight map\n")
		        do k = 1, imlen[3] {
			    do j = 1, imlen[2]
			        call aclrr (Memr[impl3r(outwt,j,k)], imlen[1])
			}
			call imunmap (outwt)
			outwt = immap (Memc[wtfname], READ_WRITE, 0)
		    }
		} else {
		    outwt = immap (weights, mode, 0)
		    if (DEBUG)
			call eprintf ("Unweighting output\n")
		    do k = 1, imlen[3] {
			do j = 1, imlen[2] {
			    call amulr (Memr[imgl3r(out,j,k)],
			        Memr[imgl3r(outwt,j,k)],
				Memr[impl3r(out,j,k)], imlen[1])
			}
		    }
		}

		if (DEBUG) {
		    call eprintf ("Using weight image: %s\n")
		        call pargstr (Memc[wtfname])
		}
	    }

	    # Set the coordinate transformations.
	    mwdim = mw_stati (outmw, MW_NPHYSDIM)
	    ctwl = mw_sctran (outmw, "world", "logical", 0)

	    # Set the geometry function and initialize the weights.
	    call gf_open (geofunc, gf, input, nin)
	    call amovkd (1D0, inpix, 3)
	    call gf_geom (gf, inpix, Memc[shape], axmap)
	    wt[1] = NULL
	    wts[1] = NULL

	    if (nlogs > 0) {
		do i = 1, nlogs {
		    call fprintf (logfiles[i], "  Output image: %s[%d")
			call pargstr (output)
			call pargi (IM_LEN(out,1))
		    do j = 2, IM_NDIM(out) {
		        call fprintf (logfiles[i], ",%d")
			    call pargi (IM_LEN(out,j))
		    }
		    call fprintf (logfiles[i], "]\n")
		}
	    }
	        
	} else {
	    # Initialize from the already open output image.
	    call amovi (IM_LEN(out,1), imlen, 3)
	    npix = imlen[1] * imlen[2] * imlen[3]

	    # Set the geometry function and check if previous weights are OK.
	    call gf_open (geofunc, gf, input, nin)
	    call amovkd (1D0, inpix, 3)
	    call gf_geom (gf, inpix, Memc[shape], axmap)
	    call sscan (Memc[shape])
	    call gargwrd (Memc[str], SZ_LINE)
	    if (streq (Memc[str], Memc[shape1])) {
	        call gargd (outworld[1])
	        call gargd (outworld[2])
	        call gargd (outworld[3])
		call gargd (outpix[1])
		outworld[1] = abs (outworld[1] - wtpars[1])
		outworld[2] = abs (outworld[2] - wtpars[2])
		outworld[3] = abs (outworld[3] - wtpars[3])
		outpix[1] = abs (outpix[1] - wtpars[4])
		if (outpix[1] > 350.)
		    outpix[1] = 360. - outpix[1]

		if (wtpars[1] != 0.)
		    outworld[1] = outworld[1] / wtpars[1]
		if (wtpars[2] != 0.)
		    outworld[2] = outworld[2] / wtpars[2]
		if (wtpars[3] != 0.)
		    outworld[3] = outworld[3] / wtpars[3]

		if (outworld[1] > TOL1 || outworld[2] > TOL1 ||
		    outworld[3] > TOL1 || outpix[1] > TOL2)
		    call tc_free (wts, no, ns)
	    } else
		call tc_free (wts, no, ns)
	}

	# Log output.
	if (nlogs > 0) {
	    do i = 1, nlogs {
	        call fprintf (logfiles[i], "  %s -> %s\n")
		    call pargstr (Memc[image])
		    call pargstr (output)
		if (bpm[1] != EOS) {
		    call fprintf (logfiles[i], "    Input mask: %s\n")
			call pargstr (Memc[bpmfname])
		}
		if (scale[1] != EOS) {
		    call fprintf (logfiles[i], "    Scale: %.4g\n")
			call pargr (s)
		}
		if (wt[1] != EOS) {
		    call fprintf (logfiles[i], "    Input weights: %s\n")
			call pargstr (Memc[wtstr])
		}
	    }
	}

	# Compute the weights.
	if (wts[1] == NULL) {
	    if (DEBUG) {
	        call eprintf ("Compute new overlap weights: %s\n")
		    call pargstr (Memc[shape])
	    }
	    call tc_weights (outmw, Memc[shape], axmap, wtt, drizscale,
	        wts, no, ns)
	    
	    do j = 1, 3 {
		nohalf[j] = (no[j] - 1) / 2
		nshalf[j] = (ns[j] - 1) / 2
	    }
	    ns12 = ns[1] * ns[2]

	    # Save parameters to decide whether weights need to be recomputed.
	    call sscan (Memc[shape])
	    call gargwrd (Memc[shape1], SZ_LINE)
	    call gargd (wtpars[1])
	    call gargd (wtpars[2])
	    call gargd (wtpars[3])
	    call gargd (wtpars[4])
	} else {
	    if (DEBUG) {
	        call eprintf ("Reuse overlap weights: %s\n")
		    call pargstr (Memc[shape])
	    }
	}

	# Create pixel list flags for the input data.  The pixel list is used
	# to keep track of the input pixels that have been used and we make
	# use of efficient routines that test for empty list or list lines.
	# This is only needed when the output has to be built up in blocks
	# rather than all in memory at once.
	#
	# Note that there is a first pass through the coordinate mapping
	# to eliminate data which falls outside the output or is in the
	# input bad pixel mask.  This is inefficient in having to evaluate
	# coordinates but it avoids I/O from the input when the output has
	# smaller coverage.

	call amovki (MAX_INT, ioutmin, 3)
	call amovki (-MAX_INT, ioutmax, 3)
	call salloc (plbuf, IM_LEN(in,1), TY_SHORT)
	pl = pl_create (3, IM_LEN(in,1), 1)
	v[1] = 1
	do v3 = 1, IM_LEN(in,3) {
	    inpix[3] = v3
	    do v2 = 1, IM_LEN(in,2) {
		inpix[2] = v2
		if (bp != NULL)
		    bpbuf = imgl3s (bp, v2, v3)
		do i = 1, IM_LEN(in,1) {
		    Mems[plbuf+i-1] = 0
		    if (bp != NULL) {
			if (Mems[bpbuf+i-1] != 0)
			    next
		    }

		    # Check for out of bounds pixels.
		    inpix[1] = i; outpix[3] = 1
		    call gf_pixel (gf, inpix, outworld)
		    call mw_ctrand (ctwl, outworld, outpix, mwdim)

		    do j = 1, 3 {
		        ioutpix[j] = nint (outpix[j])
			if (ioutpix[j] < 1 || ioutpix[j] > imlen[j])
			    break
		    }
		    if (j <= 3)
		        next

		    # Set limits.
		    call amini (ioutpix, ioutmin, ioutmin, 3)
		    call amaxi (ioutpix, ioutmax, ioutmax, 3)

		    # Flag input pixel to be used.
		    Mems[plbuf+i-1] = 1
		}
		v[2] = v2; v[3] = v3
		call plplps (pl, v, Mems[plbuf], 0, IM_LEN(in,1), PIX_SRC)
	    }
	}

	if (ioutmin[1] > ioutmax[1]) {
	    call amovki (1, ioutmin, 3)
	    call amovi (imlen, ioutmax, 3)
	}

	if (DEBUG) {
	    call eprintf ("[%d:%d,%d:%d,%d:%d]\n")
	        call pargi (ioutmin[1])
	        call pargi (ioutmax[1])
	        call pargi (ioutmin[2])
	        call pargi (ioutmax[2])
	        call pargi (ioutmin[3])
	        call pargi (ioutmax[3])
	}

	# We don't need the input bad pixel mask because we have already
	# flagged the bad input pixels in the pixel list.

	if (bp != NULL)
	    call imunmap (bp)

	# I/O is a big challenge since the relative pixel orders in the input
	# and output can be quite different.  The simplest, though not optimal,
	# algorithm is to divide the data into cubes.  We give preference
	# to longer lines by subdividing the higher dimensions first.

	i = 1
	if (npix <= max_npix) {
	    call amovki (1, ioutmin, 3)
	    call amovi (imlen, ioutmax, 3)
	}
	nvo[1] = ioutmax[1] - ioutmin[1] + 1
	nvo[2] = ioutmax[2] - ioutmin[2] + 1
	nvo[3] = ioutmax[3] - ioutmin[3] + 1
	call amovki (1, nblk, 3)
	while (nvo[1]*nvo[2]*nvo[3] > max_npix) {
	    i = mod (i+1, 3) + 1
	    if (nvo[i] >= 2 * no[i]) {
		nblk[i] = nblk[i] + 1
		nvo[i] = (ioutmax[i] - ioutmin[i] + nblk[i]) / nblk[i]
	    }
	}

	if (DEBUG) {
	    do i = 1, 3 {
		call eprintf (
		    "%d: len=%d ioutmin=%d ioutmax=%d nblk=%d nvo=%d\n")
		    call pargi (i)
		    call pargi (imlen[i])
		    call pargi (ioutmin[i])
		    call pargi (ioutmax[i])
		    call pargi (nblk[i])
		    call pargi (nvo[i])
	    }
	}

	# Loop through input and write to output.
	# The pixel list is used to keep track of the input pixels
	# that have been used and we make use of routines that test for empty
	# list or list lines.

	do blk = 0, nblk[1]*nblk[2]*nblk[3]-1 {

	    # Work in memory if small enough.
	    if (nin == 1 || npix > max_npix) {
		# Set block of output.  Include overlap pixels if not at edge.
		i = mod (blk, nblk[1])
		j = mod (blk, nblk[1]*nblk[2]) / nblk[1]
		k = blk / (nblk[1] * nblk[2])
		vo1[1] = max (1, i * nvo[1] + ioutmin[1] - nohalf[1])
		vo2[1] = min (imlen[1], (i+1) * nvo[1] + ioutmin[1] - 1 +
		    nohalf[1])
		vo1[2] = max (1, j * nvo[2] + ioutmin[2] - nohalf[2])
		vo2[2] = min (imlen[2], (j+1) * nvo[2] + ioutmin[2] - 1 +
		    nohalf[2])
		vo1[3] = max (1, k * nvo[3] + ioutmin[3] - nohalf[3])
		vo2[3] = min (imlen[3], (k+1) * nvo[3] + ioutmin[3] - 1 +
		    nohalf[3])

		buflen[1] = vo2[1] - vo1[1] + 1
		buflen[2] = vo2[2] - vo1[2] + 1
		buflen[3] = vo2[3] - vo1[3] + 1

		if (DEBUG) {
		    call eprintf ("OUT: %d=[%d:%d] %d=[%d:%d] %d=[%d:%d]\n")
		        call pargi (i+1)
			call pargi (vo1[1])
			call pargi (vo2[1])
		        call pargi (j+1)
			call pargi (vo1[2])
			call pargi (vo2[2])
		        call pargi (k+1)
			call pargi (vo1[3])
			call pargi (vo2[3])
		}

		outbuf = imps3r (out, vo1[1], vo2[1], vo1[2], vo2[2],
		    vo1[3], vo2[3])
		
		if (mode == READ_WRITE || npix > max_npix) {
		    outline = outbuf
		    do v3 = vo1[3], vo2[3] {
			do v2 = vo1[2], vo2[2] {
			    call amovr (Memr[imgl3r(out,v2,v3)+vo1[1]-1],
				Memr[outline], buflen[1])
			    outline = outline + buflen[1]
			}
		    }
		} else
		    call aclrr (Memr[outbuf], buflen[1]*buflen[2]*buflen[3])

		# Read and write weight array.
		if (Memc[wtfname] != EOS) {
		    outwtbuf = imps3r (outwt, vo1[1], vo2[1], vo1[2],
		        vo2[2], vo1[3], vo2[3])

		    if (mode == READ_WRITE || npix > max_npix) {
			if (DEBUG)
			    call eprintf ("Read weights from disk\n")
			outwtline = outwtbuf
			do v3 = vo1[3], vo2[3] {
			    do v2 = vo1[2], vo2[2] {
				call amovr (Memr[imgl3r(outwt,v2,v3)+vo1[1]-1],
				    Memr[outwtline], buflen[1])
				outwtline = outwtline + buflen[1]
			    }
			}
		    } else {
			if (DEBUG)
			    call eprintf ("Clear output weights in disk image\n")
			call aclrr (Memr[outwtbuf],
			    buflen[1]*buflen[2]*buflen[3])
		    }
		} else {
		    if (DEBUG)
		        call eprintf ("Clear in memory weights\n")
		    call calloc (outwtbuf, buflen[1]*buflen[2]*buflen[3],
		        TY_REAL)
		}
	    }

	    # Loop through the input pixels and set the output pixels which
	    # are currently in memory.

	    if (pl_empty (pl))
	        next

	    do v3 = 1, IM_LEN(in,3) {
	        v[3] = v3
		inpix[3] = v3
		do v2 = 1, IM_LEN(in,2) {
		    v[2] = v2
		    if (!pl_linenotempty (pl, v))
			next
		    call plglps (pl, v, Mems[plbuf], 0, IM_LEN(in,1), 0)

		    inbuf = imgl3r (in, v[2], v[3])
		    if (inwt != NULL)
			inwtbuf = imgl3r (inwt, v[2], v[3])

		    inpix[2] = v2
		    do i = 1, IM_LEN(in,1) {
			if (Mems[plbuf+i-1] == 0)
			    next

			pixval = s * Memr[inbuf+i-1]
			if (inwt != NULL)
			    w = Memr[inwtbuf+i-1]

			# Determine nearest output pixel for input pixel.
			inpix[1] = i; outpix[3] = 1
			call gf_pixel (gf, inpix, outworld)
			call mw_ctrand (ctwl, outworld, outpix, mwdim)

			do j = 1, 3 {
			    # Round to minimize some precision problems.
			    outpix[j] = nint (100*outpix[j]) / 100.

			    # It is more convenient to be zero indexed.
			    outpix[j] = outpix[j] - vo1[j]

			    # Fractional coordinate relative to pixel center.
			    ioutpix[j] = nint (outpix[j])
			    outpix[j] = outpix[j] - ioutpix[j]

			    # Index in lookup weights.
			    routpix[j] = (outpix[j] + 0.5) * (ns[j] - 1) + 1

			    # Offset to origin at first overlap pixel.
			    ioutpix[j] = ioutpix[j] - nohalf[j]

			    # Check for inbounds data.
			    if (vo1[j] == 1) {
				if (ioutpix[j]+nohalf[j] < 0 ||
				    ioutpix[j] >= nvo[j])
				    break
			    } else if (ioutpix[j] < 0 || ioutpix[j] >= nvo[j])
				break
			}
			if (j < 4)
			    next

			# Loop through the overlap pixels.
			# The axis mapping determines the use of the 2D and
			# 1D weights.
			do l = 0, no[3]-1 {
			    ioff[axmap[3]] = l + ioutpix[axmap[3]]
			    if (ioff[axmap[3]] < 0 ||
			        ioff[axmap[3]] >= buflen[axmap[3]])
				next
			    asi = wts[3] + l
			    if (ns[3] > 1)
				wtz = asieval (Memi[asi], routpix[axmap[3]])
			    else
				wtz = Memr[asi]
			    if (wtz == 0.)
				next

			    do k = 0, no[2]-1 {
				ioff[axmap[2]] = k + ioutpix[axmap[2]]
				if (ioff[axmap[2]] < 0 ||
				    ioff[axmap[2]] >= buflen[axmap[2]])
				    next
				msi = wts[1] + k * no[1] - 1
				do j = 0, no[1]-1 {
				    msi = msi + 1
				    ioff[axmap[1]] = j + ioutpix[axmap[1]]
				    if (ioff[axmap[1]] < 0 ||
				        ioff[axmap[1]] >= buflen[axmap[1]])
					next
				    
				    if (ns12 > 1)
					pixwt = w * msieval (Memi[msi],
					    routpix[axmap[1]],
					    routpix[axmap[2]]) * wtz
				    else
					pixwt = w * Memr[msi] * wtz
				    if (pixwt == 0.)
					next

				    outline = outbuf + (ioff[3]*buflen[2]+
					ioff[2])*buflen[1] + ioff[1]
				    Memr[outline] = Memr[outline] +
				        pixwt * pixval
				    outwtline = outwtbuf + (ioff[3]*buflen[2]+
				        ioff[2])*buflen[1] + ioff[1]
				    Memr[outwtline] = Memr[outwtline] + pixwt
				}
			    }
			}

			# Flag that we are done with the pixel.
			Mems[plbuf+i-1] = 0
		    }

		    # Update pixel flags.
		    call plplps (pl, v, Mems[plbuf], 0, IM_LEN(in,1), PIX_SRC)
		}
	    }
	}

	if (DEBUG) {
	    # Sanity check that all input pixels have been accounted for.
	    do v3 = 1, IM_LEN(in,3) {
		v[3] = v3
		do v2 = 1, IM_LEN(in,2) {
		    v[2] = v2
		    if (!pl_linenotempty (pl, v))
			next
		    call plglps (pl, v, Mems[plbuf], 0,
			IM_LEN(in,1), 0)
		    do i = 1, IM_LEN(in,1) {
			if (Mems[plbuf+i-1] == 0)
			    next
			call eprintf ("(%d,%d,%d)\n")
			    call pargi (i)
			    call pargi (v2)
			    call pargi (v3)
			break
		    }
		    if (i <= IM_LEN(in,1))
			break
		}
		if (i <= IM_LEN(in,1))
		    break
	    }
	}

	# Close the pixel list and input weights.
	call pl_close (pl)
	if (inwt != NULL)
	    call imunmap (inwt)

	# Finish up the output if needed.
	if (nin == imtlen(input)) {

	    # Normalize by the weights.  Note we only do this after the last
	    # input so if there is a crash the output will not be normalized.

	    if (npix <= max_npix) {
		inbuf = outbuf
		outline = outbuf
		outwtline = outwtbuf
	    } else {
	        call imflush (out)
		call imflush (outwt)
	    }
	    do v3 = 1, imlen[3] {
		do v2 = 1, imlen[2] {
		    if (npix > max_npix) {
			outline = impl3r (out, v2, v3)
			inbuf = imgl3r (out, v2, v3)
			outwtbuf = imgl3r (outwt, v2, v3)
			outwtline = outwtbuf
		    }
		    do i = 1, imlen[1] {
			if (outwtbuf == NULL)
			    pixwt = 0.
			else
			    pixwt = Memr[outwtline]
			if (pixwt > 0.)
			    Memr[outline] = Memr[inbuf] / pixwt
			else
			    Memr[outline] = blank
			inbuf = inbuf + 1
			outline = outline + 1
			outwtline = outwtline + 1
		    }
		}
	    }

	    # Create an output mask if desired.
	    if (masks[1] != EOS) {
		if (imaccess (masks, 0) == YES)
		    call imdelete (masks)
		bp = immap (masks, NEW_COPY, out)

		if (npix <= max_npix)
		    outwtline = outwtbuf
		do v3 = 1, imlen[3] {
		    do v2 = 1, imlen[2] {
			bpbuf = impl3s (bp, v2, v3)
			if (npix > max_npix) {
			    outwtbuf = imgl3s (outwt, v2, v3)
			    outwtline = outwtbuf
			}
			do i = 1, imlen[1] {
			    if (outwtbuf == NULL)
			        pixwt = 0.
			    else
				pixwt = Memr[outwtline]
			    if (pixwt > 0.)
				Mems[bpbuf] = 0
			    else
				Mems[bpbuf] = 1
			    bpbuf = bpbuf + 1
			    outwtline = outwtline + 1
			}
		    }
		}
		call imunmap (bp)

		call imastr (out, "BPM", masks)
	    } else if (nin == imtlen(input)) {
		iferr (call imdelf (out, "BPM"))
		    ;
	    }

	    # Close or free the output weight data.  If no output weight
	    # image was desired but we had to buffer the weights on disk
	    # then delete the temporary weight image.

	    if (Memc[wtfname] != EOS) {
		call imunmap (outwt)
		if (nin == imtlen(input) && strne (weights, Memc[wtfname])) {
		    iferr (call imdelete (Memc[wtfname]))
		        ;
		    call mfree (wtfname, TY_CHAR)
		}
	    } else
		call mfree (outwtbuf, TY_REAL)

	    # Close the output and MWCS.
	    call mw_ctfree (ctwl)
	    call mw_close (outmw)
	    call imunmap (out)

	    # Free weights.
	    call tc_free (wts, no, ns)
	}

	# Close things dealing with the input.
	call imunmap (in)
	call gf_close (gf)

	call sfree (sp)
end


# TRANSCUBE_LISTS -- Cube transform procedure with lists.
#
# This procedure expands the lists and call TRANSCUBE for each input.

procedure transcube_list (input, output, masks, weights, bpm, scale, wt,
	wcsreference, wttype, drizscale, blank, geofunc, memalloc,
	logfiles, nlogs)

int	input				#I Input image list
int	output				#I Output image list
int	masks				#I Output mask list
int	weights				#I Output weights list
int	bpm				#I Input bad pixel list
int	scale				#I Input scale list
int	wt				#I Input weight list
char	wcsreference[ARB]		#I WCS reference
char	wttype[ARB]			#I Weighting type
real	drizscale[3]			#I Drizzle scale factors
real	blank				#I Blank value
char	geofunc[ARB]			#I Geometry function
real	memalloc			#I Memory to alloc in Mb
int	logfiles[ARB]			#I Output logfile descriptors
int	nlogs				#I Number of logfile descriptors

int	i, nin, nout
pointer	sp, in, out, outmask, outwt, inmask, scl, inwt, str

int	imtlen(), imtrgetim()

begin
	call smark (sp)
	call salloc (in, SZ_FNAME, TY_CHAR)
	call salloc (out, SZ_FNAME, TY_CHAR)
	call salloc (outmask, SZ_FNAME, TY_CHAR)
	call salloc (outwt, SZ_FNAME, TY_CHAR)
	call salloc (inmask, SZ_FNAME, TY_CHAR)
	call salloc (scl, SZ_FNAME, TY_CHAR)
	call salloc (inwt, SZ_FNAME, TY_CHAR)
	call salloc (str, SZ_LINE, TY_CHAR)

	# First check of input parameters.
	nin = imtlen (input)
	nout = imtlen (output)
	if (nout != 1 && nout != nin)
	    call error (1, "Invalid output list")
	i = imtlen (bpm)
	if (i > 1 && i != nin)
	    call error (1, "Invalid input bad pixel list")
	i = imtlen (scale)
	if (i > 1 && i != nin)
	    call error (1, "Invalid input scale list")
	i = imtlen (wt)
	if (i > 1 && i != nin)
	    call error (1, "Invalid input weight list")
	i = imtlen (masks)
	if (i > 0 && i != nout)
	    call error (1, "Invalid output mask list")
	i = imtlen (weights)
	if (i > 0 && i != nout)
	    call error (1, "Invalid output weight list")

	# Initialize.
	Memc[outmask] = EOS
	Memc[outwt] = EOS
	Memc[inmask] = EOS
	Memc[scl] = EOS
	Memc[inwt] = EOS

	# Loop through the input list.
	do i = 1, nin {
	    if (imtrgetim (bpm, i, Memc[str], SZ_LINE) != EOF)
	        call strcpy (Memc[str], Memc[inmask], SZ_FNAME)
	    if (imtrgetim (output, i, Memc[str], SZ_LINE) != EOF)
	        call strcpy (Memc[str], Memc[out], SZ_FNAME)
	    if (imtrgetim (weights, i, Memc[str], SZ_LINE) != EOF)
	        call strcpy (Memc[str], Memc[outwt], SZ_FNAME)
	    if (imtrgetim (masks, i, Memc[str], SZ_LINE) != EOF)
	        call strcpy (Memc[str], Memc[outmask], SZ_FNAME)
	    if (imtrgetim (scale, i, Memc[str], SZ_LINE) != EOF)
	        call strcpy (Memc[str], Memc[scl], SZ_FNAME)
	    if (imtrgetim (wt, i, Memc[str], SZ_LINE) != EOF)
	        call strcpy (Memc[str], Memc[inwt], SZ_FNAME)

	    call transcube (input, i, Memc[out],
		Memc[outmask], Memc[outwt], Memc[inmask], Memc[scl],
		Memc[inwt], wcsreference, wttype, drizscale,
		blank, geofunc, memalloc, logfiles, nlogs)
	}

	call sfree (sp)
end
