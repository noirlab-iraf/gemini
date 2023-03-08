# Copyright(c) 2004-2006 Association of Universities for Researchh in Astronomy, Inc.

# Original author: Kathleen Labrie   17-Aug-2004

include "glog.h"
include "gemerrors.h"

.help
.nf
This file contains procedures to format the parameter string to be 
fed to gloginit().  Each routine append a new line to the parameter
list.

   glogpstring - String value
        glogpc - Character value
        glogps - Short value
        glogpi - Integer value
        glogpl - Long value
        glogpr - Real value
        glogpd - Double value
        glogpb - Boolean value
        glogpx - Complex value

GLOGPSTRING -- Append info for a string value.  All the other GLOGPt 
               procedures eventually call this routine once the value has
	       been converted to a string.

	call glogpstring (paramstr, param, valuestr)
	
	paramstr	: Formatted parameter string	[input/output, (string)]
	param		: Name of the parameter		[input, (string)]
	valuestr	: Value, as a string		[input, (string)]

GLOGPt -- Append info for a non-string value.  Replace the 't' by 
          first character of the data type.  For example, glogpi() for
	  integers, glogpb() for boolean, etc.  The routines are simple
	  wrappers around glogpstring().

	call glogpT (paramstr, param, value)
	
	paramstr	: Formatted parameter string	[input/output, (string)]
	param		: Name of the parameter		[input, (string)]
	value		: Value (non-string)		[input, (non-string)]
.fi
.endhelp


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------


# GLOGPSTRING -- Append info for a string value.  All the other GLF 
#                procedures eventually call this routine once the value has
#	         been converted to a string.

procedure glogpstring (paramstr, param, valuestr)

char	paramstr[G_SZ_PARAMSTR]	#IO The formatted parameter string
char	param[ARB]		#I  The name of the parameter
char	valuestr[ARB]		#I  The parameter value, as a string

# Other variables
int	len

# Gemini functions
bool	g_whitespace()

# IRAF functions
int	strlen()

begin
	# First, make sure the string won't get too long
	len = strlen (paramstr) + strlen (valuestr) + 24
	if ( len >= G_SZ_PARAMSTR )
	   call error (G_INTERNAL_ERROR, "Parameter string is too long.")
	
	# Append the new parameter/value pair
	if ( ! g_whitespace (paramstr) ) {
	    call sprintf (paramstr, G_SZ_PARAMSTR, "%s\n")
	        call pargstr (paramstr)
	}
	call sprintf (paramstr, G_SZ_PARAMSTR, "%s%-20s = %s")
	    call pargstr (paramstr)
	    call pargstr (param)
	    call pargstr (valuestr)

	return
end

#--------------------------------------------------------------------------

# GLOGPI -- Append info for an integer value.

procedure glogpi (paramstr, param, value)

char	paramstr[G_SZ_PARAMSTR]	#IO The formatted parameter string
char	param[ARB]		#I  The name of the parameter
int	value			#I  The integer value

# Other variables
char	valuestr[SZ_LINE]

errchk	glogpstring()

begin
	call sprintf (valuestr, SZ_LINE, "%d")
	    call pargi (value)
	call glogpstring (paramstr, param, valuestr)
end

#--------------------------------------------------------------------------

# GLOGPB -- Append info for a boolean value.

procedure glogpb (paramstr, param, value)

char	paramstr[G_SZ_PARAMSTR]	#IO The formatted parameter string
char	param[ARB]		#I  The name of the parameter
bool	value			#I  The integer value

# Other variables
char	valuestr[SZ_LINE]

errchk	glogpstring()

begin
	call sprintf (valuestr, SZ_LINE, "%b")
	    call pargb (value)
	call glogpstring (paramstr, param, valuestr)
end

#--------------------------------------------------------------------------

# GLOGPR -- Append info for a real value.

procedure glogpr (paramstr, param, value)

char	paramstr[G_SZ_PARAMSTR]	#IO The formatted parameter string
char	param[ARB]		#I  The name of the parameter
real	value			#I  The integer value

# Other variables
char	valuestr[SZ_LINE]

errchk	glogpstring()

begin
	call sprintf (valuestr, SZ_LINE, "%g")
	    call pargr (value)
	call glogpstring (paramstr, param, valuestr)
end

#--------------------------------------------------------------------------

# GLOGPC -- Append info for a character value.

procedure glogpc (paramstr, param, value)

char	paramstr[G_SZ_PARAMSTR]	#IO The formatted parameter string
char	param[ARB]		#I  The name of the parameter
char	value			#I  The integer value

# Other variables
char	valuestr[SZ_LINE]

errchk	glogpstring()

begin
	call sprintf (valuestr, SZ_LINE, "%c")
	    call pargc (value)
	call glogpstring (paramstr, param, valuestr)
end

#--------------------------------------------------------------------------

# GLOGPS -- Append info for a short integer value.

procedure glogps (paramstr, param, value)

char	paramstr[G_SZ_PARAMSTR]	#IO The formatted parameter string
char	param[ARB]		#I  The name of the parameter
short	value			#I  The integer value

# Other variables
char	valuestr[SZ_LINE]

errchk	glogpstring()

begin
	call sprintf (valuestr, SZ_LINE, "%d")
	    call pargs (value)
	call glogpstring (paramstr, param, valuestr)
end

#--------------------------------------------------------------------------

# GLOGPL -- Append info for a long integer value.

procedure glogpl (paramstr, param, value)

char	paramstr[G_SZ_PARAMSTR]	#IO The formatted parameter string
char	param[ARB]		#I  The name of the parameter
long	value			#I  The integer value

# Other variables
char	valuestr[SZ_LINE]

errchk	glogpstring()

begin
	call sprintf (valuestr, SZ_LINE, "%d")
	    call pargl (value)
	call glogpstring (paramstr, param, valuestr)
end

#--------------------------------------------------------------------------

# GLOGPD -- Append info for a double precision value.

procedure glogpd (paramstr, param, value)

char	paramstr[G_SZ_PARAMSTR]	#IO The formatted parameter string
char	param[ARB]		#I  The name of the parameter
double	value			#I  The integer value

# Other variables
char	valuestr[SZ_LINE]

errchk	glogpstring()

begin
	call sprintf (valuestr, SZ_LINE, "%g")
	    call pargd (value)
	call glogpstring (paramstr, param, valuestr)
end

#--------------------------------------------------------------------------

# GLOGPX -- Append info for a complex value.

procedure glogpx (paramstr, param, value)

char	paramstr[G_SZ_PARAMSTR]	#IO The formatted parameter string
char	param[ARB]		#I  The name of the parameter
complex	value			#I  The integer value

# Other variables
char	valuestr[SZ_LINE]

errchk	glogpstring()

begin
	call sprintf (valuestr, SZ_LINE, "%z")
	    call pargx (value)
	call glogpstring (paramstr, param, valuestr)
end
