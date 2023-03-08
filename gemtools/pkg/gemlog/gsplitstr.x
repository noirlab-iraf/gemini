# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

.help
.nf
G_SPLITSTR -- Given a token, split a string.
	nstr  = g_splitstr( strbuf, token, strptr )
	
	nstr		: Number of strings  [return value, (int)]
	strbuf		: String buffer  [input/output, (char[])]
	token		: Token  [input (char)]
	strptr		: Array of indices pointing into strbuf [output, (int[])]
.fi
.endhelp

# Kathleen Labrie   29-Apr-2004

int procedure g_splitstr (strbuf, token, strptr)

char	strbuf[ARB]	#IO String buffer
char	token		#I  Character token
int	strptr[ARB]	#O  Array of indices pointing into strbuf

int	nstr		#O  Number of strings found

# Other variables
int	len, i, ntoken

# IRAF functions
int strlen()

begin
	len = strlen ( strbuf )

	ntoken = 0
	strptr[1] = 1
	for (i = 1; i <= len; i = i+1) {
	    if (strbuf[i] == token) {
		strbuf[i] = EOS
		ntoken = ntoken + 1
		strptr[ntoken+1] = i+1
	    }
	}

	nstr = ntoken+1
	return (nstr)
end
