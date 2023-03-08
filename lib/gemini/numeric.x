# Copyright(c) 2002-2005 Association of Universities for Research in Astronomy, Inc.
#
# Some subroutines for numerical operations, supporting the gemini tasks
#
# Version     Sep  3, 2002  JT
# ----

include <mach.h>


int procedure gnu_equalr(a, b)

# Determine whether two real numbers are equal within some rounding err

real a, b, thresh, diff

begin

  thresh = EPSILONR
  diff = a-b
	
  if (diff < thresh && diff > -thresh) return YES
  else return NO
	
end
