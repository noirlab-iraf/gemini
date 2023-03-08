# Copyright(c) 1986 Association of Universities for Research in Astronomy Inc.

include	<syserr.h>
include	"idb.h"

# IMGCOM -- Get the comment field for a keyword.

procedure imgcom (im, key, comment)

pointer	im			#I image descriptor
char	key[ARB]		#I parameter to be set
char	comment[ARB]		#O comment string

bool	string_valued
int	ch, i, n, j
pointer	rp, ip, sp, val
int	idb_findrecord()
errchk	syserrs

begin
	call smark (sp)
	call salloc (val, SZ_LINE, TY_CHAR)

	# Find the record.
	if (idb_findrecord (im, key, rp) == 0)
	    call syserrs (SYS_IDBKEYNF, key)

	# Determine the actual datatype of the parameter.  String valued
	# parameters will have an apostrophe in the first nonblank column
	# of the value field.

	string_valued = false
	for (ip=IDB_STARTVALUE;  ip <= IDB_ENDVALUE;  ip=ip+1) {
	    # Skip leading whitespace.
	    for (; Memc[rp+ip-1] == ' '; ip=ip+1) 
	        ;

	    if (Memc[rp+ip-1] == '\'') {
		# Get string value.
		do i = ip, IDB_RECLEN {
		    ch = Memc[rp+i]
		    if (ch == '\n')
			break
		    if (ch == '\'')
			break
		}
		i = i + 2
		break
	    } else {
		# Numeric value.
		do i = ip, IDB_RECLEN {
		    ch = Memc[rp+i-1]
		    if (ch == '\n' || ch == ' ' || ch == '/')
			break
		}
		break
	    }
	}

	for (; Memc[rp+ip-1] == ' '; ip=ip+1) 
	        ;
	n = 0
	do j = IDB_RECLEN, i ,-1 {
	   ch = Memc[rp+j]
	   n = n + 1
	   if (ch != ' ' && ch != '\n') 
	      break
	}
	for (j=0; Memc[rp+i+1+j] == ' '; j=j+1)
	        ;
	call strcpy (Memc[rp+i+1+j], comment, IDB_RECLEN-i-n-j+1)
call eprintf("COMMDD: <%s>, ip: %d, i: %d, IDB_RECLEN-i-n: %d\n")
   call pargstr(comment)
   call pargi(ip)
   call pargi(i)
   call pargi(IDB_RECLEN-i-n)

	call sfree (sp)
end
