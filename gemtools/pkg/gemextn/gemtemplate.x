# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	<error.h>
include	<syserr.h>
include	"gemerrors.h"

# This is a limited version of the imt routines which supports one level
# of @-file indirection.  This version does not provide getting an image
# by index or the length of the list.  If an application wants the length
# it should loop through all the names using gemtgetim and then call
# gemtrew.  This should be done immediately after opening or a rewind
# because the position in the list would be lost.


# GEMTOPEN -- Open a Gemini template.
# Raises error if file not opened

pointer procedure gemtopen (template)

char	template[ARB]		#I Gemini template
pointer	twoimt			#O Gemini list

pointer imt
pointer	imtopen()

errchk	imtopen()

begin
	imt = imtopen (template)
	call calloc (twoimt, 2, TY_POINTER)
	Memi[twoimt] = imt
	Memi[twoimt+1] = NULL
	return (twoimt)
end


# GEMTGETIM -- Get the next image in the template.

int procedure gemtgetim (imt, image, maxch)

pointer	imt			#I Gemini list
char	image[maxch]		#O Image name
int	maxch			#O Maximum length of image name
int	stat			#R Return status

pointer	tmp

int	imtopen(), imtgetim(), stridx()

begin
	repeat {
	    # If a secondary list is open read from it.  At the end of
	    # the secondary list close it and proceed to the next image
	    # in the primary list.

	    if (Memi[imt+1] != NULL) {
		stat = imtgetim (Memi[imt+1], image, maxch)
		if (stat != EOF)
		    return (stat)
		call imtclose (Memi[imt+1])
		Memi[imt+1] = NULL
	    }

	    # Read the next element in the primary list.
	    # Check for an @-file and open a secondary list.

	    stat = imtgetim (Memi[imt], image, maxch)
	    if (stat != EOF) {
		if (stridx('@', image) > 0) {
		    iferr (tmp = imtopen (image)) {
			call gemtclose (imt)
			call erract (EA_ERROR)
		    } else {
			Memi[imt+1] = tmp
		    }
		    next
		}
	    }

	    return (stat)
	}
end


# GEMTCLOSE -- Close the Gemini template.

procedure gemtclose (imt)

pointer	imt			#I Gemini list

begin
	if (Memi[imt+1] != NULL)
	    call imtclose (Memi[imt+1])
	call imtclose (Memi[imt])
	call mfree (imt, TY_POINTER)
end


# GEMTREW -- Rewind the template.

procedure gemtrew (imt)

pointer	imt			#I Gemini list

begin
	if (Memi[imt+1] != NULL) {
	    call imtclose (Memi[imt+1])
	    Memi[imt+1] = NULL
	}
	call imtrew (Memi[imt])
end

