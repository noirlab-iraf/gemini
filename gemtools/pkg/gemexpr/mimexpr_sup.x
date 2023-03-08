# Copyright(c) 2004-2009 Association of Universities for Research in Astronomy, Inc.

include	<ctotok.h>
include	<imhdr.h>
include	<ctype.h>
include	<mach.h>
include	<imset.h>
include	<fset.h>
include	<lexnum.h>
include	<evvexpr.h>
include "../../../lib/mefio/mefio.h"
include	"gettok.h"
include "gemexpr.h"


# IE_EXPANDTEXT -- Scan an expression, performing macro substitution on the
# contents and returning a fully expanded string.

pointer procedure ie_expandtext (st, expr)

pointer	st			#I symbol table (macros)
char	expr[ARB]		#I input expression

pointer	buf, gt
int	buflen, nchars
pointer	locpr()
int	gt_expand()
pointer	gt_opentext()
pointer ie_gsym()
extern	ie_gsym()

begin
	buflen = SZ_COMMAND
	call malloc (buf, buflen, TY_CHAR)

	gt = gt_opentext (expr, locpr(ie_gsym), st, 0, GT_NOFILE)
	nchars = gt_expand (gt, buf, buflen)
	call gt_close (gt)

	return (buf)
end


# IE_GETOPS -- Parse the expression and generate a list of input operands.
# The output operand list is returned as a sequence of EOS delimited strings.

int procedure ie_getops (st, expr, oplist, maxch)

pointer	st			#I symbol table
char	expr[ARB]		#I input expression
char	oplist[ARB]		#O operand list
int	maxch			#I max chars out

int	noperands, ch, i
int	ops[MAX_OPERANDS]
pointer	gt, sp, tokbuf, op
bool ldebug
pointer ie_gsym()
extern	ie_gsym()
pointer	gt_opentext()
pointer	locpr(), gt_rawtok(), gt_nexttok()
errchk	gt_opentext, gt_rawtok

begin
    ldebug = false
    
	call smark (sp)
	call salloc (tokbuf, SZ_LINE, TY_CHAR)

    
	call aclri (ops, MAX_OPERANDS)
	gt = gt_opentext (expr, locpr(ie_gsym), st, 0, GT_NOFILE+GT_NOCOMMAND)

	# This assumes that operand names are the letters "a" to "z".
	while (gt_rawtok (gt, Memc[tokbuf], SZ_LINE) != EOF) {
	    ch = Memc[tokbuf]
        if (ldebug) {
            call printf("tokbuf=%s\n")
            call pargstr(Memc[tokbuf])
            call flush(STDOUT)
        }    
	    if (IS_LOWER(ch) && ((Memc[tokbuf+1] == EOS) || Memc[tokbuf+1] == '.')) {
		    if (gt_nexttok (gt) != '(') {
		        ops[ch-'a'+1] = 1
            }
        }
	}

	call gt_close (gt)

	op = 1
	noperands = 0
	do i = 1, MAX_OPERANDS
	    if (ops[i] != 0 && op < maxch) {
		oplist[op] = 'a' + i - 1
		op = op + 1
		oplist[op] = EOS
		op = op + 1
		noperands = noperands + 1
	    }

	oplist[op] = EOS
	op = op + 1

	call sfree (sp)
	return (noperands)
end


# IE_GSYM -- Get symbol routine for the gettok package.

pointer procedure ie_gsym (st, symname, nargs)

pointer	st			#I symbol table
char	symname[ARB]		#I symbol to be looked up
int	nargs			#O number of macro arguments

pointer	sym
pointer	strefsbuf(), stfind()

begin
	sym = stfind (st, symname)
	if (sym == NULL)
	    return (NULL)

	nargs = SYM_NARGS(sym)
	return (strefsbuf (st, SYM_TEXT(sym)))
end
