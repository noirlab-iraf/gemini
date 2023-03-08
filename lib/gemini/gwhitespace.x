# Copyright(c) 2004-2005 Association of Universities for Research in Astronomy, Inc.

include <ctype.h>

# G_WHITESPACE -- Is string empty?

bool procedure g_whitespace(str)

char	str[ARB]		# I The string to test

int	i, len
int	strlen()

begin
	len = strlen (str)
	do i = 1, len {
	    if (! IS_WHITE(str[i]))
		return false

	}
	return true
end
