# Copyright(c) 2005-2006 Association of Universities for Research in Astronomy, Inc.

include <imhdr.h>

# A very basic routine that modifies the input by scanning along rows
# replacing values according to whether or not they are above or below 
# a certain threshold:
# - Points below the threshold are replaced by zero.
# - Contiguous sections above the threshold are replaced by the square 
#   of the distance from the nearest section edge.
#
# For example, with a threshold value of 3
#
# Data     0  1  -1  5  6  4  4 -1 -2  2  9  8  7 
# Above?   n  n   n  y  y  y  y  n  n  n  y  y  y
# Distance -  -   -  1  2  2  1  -  -  -  1  2  1
# Result   0  0   0  1  4  4  1  0  0  0  1  4  1
#
# This is rather arbitrary, but gives converts a flat into something
# vaguely like a pinhole observation, as required by the IFU.
#
# The details of calculating the threshold and determining the
# dispersion axis (with rotation if necessary) are left to the calling
# CL script (I'm pushed for time and am more used to CL than SPP at
# the moment; this is in SPP only because it process data pixel by
# pixel and would be too slow in CL).

procedure t_peakhelper ()

char    inimage[SZ_FNAME]       # I Image to process
char    outimage[SZ_FNAME]      # O Modified image
real    threshold               # I Cutoff threshold

pointer inimg, outimg, inline, outline
int     nx, ny, j
real    clgetr()
pointer immap(), imgl2r(), impl2r()

begin
	
        call clgstr ("inimage", inimage, SZ_FNAME)
        call clgstr ("outimage", outimage, SZ_FNAME)
        threshold = clgetr ("threshold")

        inimg = immap (inimage, READ_ONLY, NULL)
        outimg = immap (outimage, NEW_COPY, inimg)
        nx = IM_LEN (inimg, 1)
        ny = IM_LEN (inimg, 2)

        for (j = 1; j <= ny; j = j + 1) {
            inline = imgl2r (inimg, j)
            outline = impl2r (outimg, j)
            call triangles (Memr[inline], Memr[outline], nx, threshold)
        }

        call imunmap (inimg)
        call imunmap (outimg)

end


# Implement the algorithm described at the head of this file,
# replacing values above a threshold with "triangles".

procedure triangles (inline, outline, nx, threshold)

real    inline[ARB]     # I existing pixels
real    outline[ARB]    # O triangle pixels
int     nx              # I number of pixels
real    threshold       # I cutoff threshold

int     i
int     left, right
real    val

begin

        right = 0
        left = 0

        for (i = 1; i <= nx; i = i + 1) {

            # If i is to the right of right, then we are looking
            # for a new region
            if (i > right) {
                if (inline[i] < threshold) {
                    # Set low values to zero
                    outline[i] = 0.0
                } else {
                    # Initialise a new region
                    left = i
                    right = i
                    while (inline[right] >= threshold && right < nx) {
                        right = right + 1
                    }
                    # line will be modified below
                }
            }

            # If i is to the left of (or at) right, then we are
            # building a triangle
            if (i <= right) {
                val = min (i - left, right - i)
                val = val + 1.0
                val = val * val
                outline[i] = val
            }
        }

end
