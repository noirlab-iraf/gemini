# Copyright(c) 2004-2013 Association of Universities for Research in Astronomy, Inc.

include <evvexpr.h>
include <lexnum.h>
include <ctype.h>
include <mach.h>

include "gemexpr.h"
include "../gemlog/glog.h"

define GEMEXPR_SUBTASK "gemexprpars"

define EXPRTYPE_SCI 0
define EXPRTYPE_VAR 1
define EXPRTYPE_DQ  2

# This file contains:
#
#       t_gemarith(...) - gemarith routine in SPP
#       t_gemexpr(...)  - gemexpr routing in SPP, called by gemarith
#

define CARD_LEN 72  #room for full card plus EOS and an extra byte
                                       #  for word     alignment
define OPSTRLEN 256
define GE_OPNAM_LEN 32
define FULLOPNAM_LEN 32
define GE_MAPSCIOP_LEN 32

procedure t_gemarith()
# Local variables for task parameters
pointer op1,op2
pointer operator
pointer exprstr
pointer sp
pointer outtype
bool ldebug
pointer mepa, mepb
bool aisim, bisim
char aopstr[OPSTRLEN], bopstr[OPSTRLEN]
pointer memap()
char tmpstr[SZ_LINE], errmsg[SZ_LINE]
char outfn[SZ_FNAME]
int tmpi
double tmpd
pointer tmpp

int		l_status

int errget()
int strcmp(), strncmp()
int strlen()
bool streq()
int ctod()
int mehct() # mef find Highest Common Type among extensions
int typea, typeb, intype


begin
    
    call log_init("gemarith","gemtools")

    l_status = 1;    

    mepa = 0
    mepb = 0
    ldebug = false
    call log_info("begin (this task relies on gemexpr)")
    
    call smark(sp)
    
    call salloc (op1,SZ_LINE, TY_CHAR)
    call salloc (op2,SZ_LINE, TY_CHAR)
    call salloc (operator, SZ_LINE, TY_CHAR)
    call salloc (exprstr, SZ_LINE, TY_CHAR)
    call salloc (outtype, SZ_LINE, TY_CHAR)
    
    Memc[exprstr]=0
    call clgstr("operand1", Memc[op1], SZ_LINE)
    call clgstr("op", Memc[operator], SZ_LINE)
    call clgstr("operand2", Memc[op2], SZ_LINE)
    
    # transmit output to gemexpr
    call clgstr("result", outfn, SZ_FNAME)
    call clpstr("gemexprpars.output", outfn)
    
# gemexpr bug workaround
# check type of operand
    iferr (mepa = memap(Memc[op1]))
    {
        tmpp = op1
        if ((strncmp ("operand", Memc[op1], 7) == 0) || (ctod(Memc,tmpp, tmpd) > 0))  {
            aisim = false
            call strcpy (Memc[op1],aopstr,OPSTRLEN)
        } else {
            call sprintf(tmpstr, SZ_LINE, "operand1 (%s) does not exist as a MEF\n  (and is not recognized as a constant\n      or FITS header reference)")
            call pargstr(Memc[op1])
            call task_error(GEERR_INVALOPERANDS, tmpstr)
            return
        }
    }
    else
    {
        aisim = true
        call flush(STDOUT)
        call clpstr("gemexprpars.a", Memc[op1])
        call flush(STDOUT)
        call strcpy ("a",aopstr,OPSTRLEN)
        typea = mehct(mepa)
        
        call meunmap(mepa)
    }
    
    iferr ( mepb = memap(Memc[op2]))
    {
        tmpp = op2
        if ((strncmp ("operand", Memc[op2], 7) == 0) || (ctod(Memc,tmpp, tmpd) > 0))  {
            bisim = false
            call strcpy (Memc[op2],bopstr,OPSTRLEN)
        } else {
            call sprintf(tmpstr, SZ_LINE, "operand2 (%s) does not exist as a MEF\n  (and is not recognized as a constant\n      or FITS header reference)")
            call pargstr(Memc[op2])
            call task_error(GEERR_INVALOPERANDS, tmpstr)
            return
        }
     }
    else
    {
        bisim = true
        call clpstr("gemexprpars.b", Memc[op2])
        call strcpy ("b",bopstr,OPSTRLEN)
        
        typeb = mehct(mepb)
        
        call meunmap(mepb)
    }
    
    #final input type
    if (typea >= typeb) {
        intype = typea
    } else {
        intype = typeb
    }

    if (ldebug) {
        call printf ("%d %d --> %d\n")
        call pargi (typea)
        call pargi (typeb)
        call pargi (intype)
        call flush (STDOUT)
    }

    #check operands for operand references, 
    # brute force permuted conditionals 
    if (strncmp(aopstr,"operand1", 8)==0) {
        call task_error(GEERR_INVALOPERANDS, 
           "self reference in parameter operand1")
        return
    }
    if (strncmp(aopstr,"operand2", 8)==0) {
        aopstr[1] = 'b'
        call amovc(aopstr[9], aopstr[2], strlen(aopstr)-7)
        if (ldebug) {
            call printf("operand1 now xformed: %s\n")
            call pargstr(aopstr)
            call flush(STDOUT)
        }
    }
    if (strncmp(bopstr,"operand2", 8)==0) {
        call task_error(GEERR_INVALOPERANDS, 
           "self reference in parameter operand2")
        return
    }
    if (strncmp(bopstr,"operand1", 8)==0) {
        bopstr[1] = 'a'
        call amovc(bopstr[9], bopstr[2], strlen(bopstr)-7)

        if (ldebug) {
            call printf("operand2 now xformed: %s\n")
            call pargstr(bopstr)
            call flush(STDOUT)
        }
    }
    
    # now check/map refim
    # expected to be default, operand1, or operand2
    call clgstr("refim", tmpstr, SZ_LINE)
    if (false) {
        call printf("refim=%s\n")
        call pargstr(tmpstr)
        call flush(STDOUT)
    }
    if (streq(tmpstr,"default") == false) {
        # it's not default, make sure it is an operand
        if (streq("operand1",tmpstr) == true) {
            call clpstr("gemexprpars.refim","a")       
        } else {
            if (streq("operand2",tmpstr) == true) {
                call clpstr("gemexprpars.refim","b")       
            } else {
                call task_error(GEERR_BADREFIM, "refim must be \"default\", \"operand1\", or \"operand2\"")
                return
            }
        }       
    } else {
        call clpstr("gemexprpars.refim","default")
    }
        
        
    call sprintf(Memc[exprstr],SZ_LINE,"%s %s %s")
    call pargstr(aopstr)
    call pargstr(Memc[operator])
    call pargstr(bopstr)
    call clpstr("gemexprpars.sci_expr", Memc[exprstr])
    
    if (ldebug)
        {
        call printf("sci_expr=%s\n")
        call pargstr(Memc[exprstr])
        call flush(STDOUT)
    }
    
    if (strcmp("+", Memc[operator]) == 0 || strcmp("-", Memc[operator])==0)
        {
        if (aisim && bisim)
            {
            call clpstr("gemexprpars.var_expr", "a[VAR]+b[VAR]")
        }
        else if (aisim)
        {
            call clpstr("gemexprpars.var_expr", "a[VAR]")
        }
        else if (bisim)
        {
            call clpstr("gemexprpars.var_expr", "b[VAR]")
        }
        else
        {
            call task_error(GEERR_INVALOPERANDS, "Neither operand is an image")
            return
        }
    }
    else if (strcmp("*", Memc[operator]) == 0)
    {
        if (aisim && bisim) {
            call sprintf(Memc[exprstr], SZ_LINE, "(a[VAR]*b[SCI]**2)+(b[VAR]*a[SCI]**2)")
            call clpstr("gemexprpars.var_expr", Memc[exprstr])
        }
        else if (aisim) # b is constant
        {
            call sprintf(Memc[exprstr], SZ_LINE, "(a[VAR]*(%s)**2)")
            call pargstr(bopstr)
            call clpstr("gemexprpars.var_expr", Memc[exprstr])
        }
        else if (bisim) # a is constant
        {
            call sprintf(Memc[exprstr], SZ_LINE, "(b[VAR]*(%s)**2)")
            call pargstr(aopstr)
            call clpstr("gemexprpars.var_expr", Memc[exprstr])
        }
        else
        {
            call task_error(GEERR_NEEDIMFORVAR, "Neither operand is an image" )
            return
        }
    } else if (strcmp("/", Memc[operator])==0) {
        if (aisim && bisim) {
            call sprintf(Memc[exprstr], SZ_LINE, "( a[VAR] + a[SCI]**2 * b[VAR] / b[SCI]**2 ) / b[SCI]**2")
            call clpstr("gemexprpars.var_expr", Memc[exprstr])
        }
        else if (aisim) # b is constant
        {
            call sprintf(Memc[exprstr], SZ_LINE,  " a[VAR] / (%s)**2")
            call pargstr(bopstr)
            call clpstr("gemexprpars.var_expr", Memc[exprstr])
        }
        else if (bisim) # a is constant
        {
            call sprintf(Memc[exprstr], SZ_LINE,  "( (%s)**2 * b[VAR] / b[SCI]**2 ) / b[SCI]**2")
            call pargstr(aopstr)
            call clpstr("gemexprpars.var_expr", Memc[exprstr])
        }
        else {
            call task_error(GEERR_NEEDIMFORVAR, "Neither operand is an image" )
            return
        }        
    } else {
        call sprintf(tmpstr, SZ_LINE, "Bad Operator \"%s\"")
        call pargstr(Memc[operator])
        call task_error(GEERR_BADOPERATOR, tmpstr)
        return
    }
    
    # JUGGLE OUTTYPE
    call clgstr("outtype", outtype, SZ_LINE)
    
    if (strcmp(outtype, "default") == 0) {
        if (intype >= TY_DOUBLE) {
            call clpstr("outtype", "double")
        } else {
            call clpstr("outtype", "real")
        }
    }
    
    if (ldebug)
        {
        call printf("var_expr=%s\n")
        call pargstr(Memc[exprstr])
        call flush(STDOUT)
    }
    iferr( call w_gemexpr(GEMEXPR_SUBTASK)) {
        tmpi = errget(errmsg, SZ_LINE)
        call sprintf(tmpstr, SZ_LINE, "(%d) %s\nGEMEXPR subsystem reports error")
        call pargi(tmpi)
        call pargstr(errmsg)
        call task_error(GEERR_W_GEMEXPRFAILED, tmpstr)
        call log_close()
        call sfree(sp)
        call clputi("status",1)
        if (tmpi != GEERR_OUTPUTEXISTS) {
            iferr (call imdelete(outfn));
        }
        return
    }
    # fix the outtype juggle above so it's right in uparm
    call clpstr("outtype", outtype)
    call mtstamp(outfn, "GEMARITH")
    

# Close the logfile
	#     The memory for gl is freed in glogclose (whether glogclose is
	#     successful or not).
	# Free memory
    call log_info("gemarith: done")
    call log_close()
    call sfree(sp)
    
    # if I got here it was success.
    call clputi("status",0)
end
	
procedure t_gemexpr()
# local vars
int errnum
char errmsg[SZ_LINE]
char outfn[SZ_FNAME]
int errget()
begin
	
    call log_init("gemexpr", "gemtools")
    call log_info("begin")
    iferr( call w_gemexpr(""))
    {
        errnum = errget(errmsg, SZ_LINE)
        call task_error(errnum, errmsg)
        call log_close()
        call clputi("status", 1)
        
        if (errnum != GEERR_OUTPUTEXISTS) {
            call clgstr("output", outfn, SZ_FNAME)
            iferr (call imdelete(outfn)) ;
        }
            
        return
    }

    call clgstr("output", outfn, SZ_FNAME)
    call mtstamp(outfn, "GEMEXPR")
   
    call log_info("done")
    call log_close()
	call clputi("status",0)
end
	

procedure w_gemexpr(namespace)
char namespace[ARB]
#local variables
pointer sp
pointer sci_ext, var_ext, dq_ext
pointer sci_expr, var_expr, dq_expr
pointer ge_output
pointer oplist, opnam, opval
pointer fname,st
pointer mdfsource       # filename for MEF with MDF in it
bool    ldebug
int     i
char    ge_opnam[GE_OPNAM_LEN]
char	fullopnam[FULLOPNAM_LEN]
char    ge_mapsciop[GE_MAPSCIOP_LEN]
char	refimnam[SZ_FNAME]
char	tmpstr[SZ_LINE] # reusable tmpstr
int 	tmpi
char	tc	#tmp char
bool    fl_vardq        #propagate VAR and DQ?
bool    verbose
int     noperands
pointer ie_getexprdb()
int 	strcmp()
int		errget()
int 	ie_getops()
bool 	strne(), streq()
bool 	clgetb()
int		clgeti()
int 	access(),imaccess()
int		strlen()
bool 	allwhite()
#---
include "../../../lib/mefio/mefiocommon.h"
include "nscommon.h"
#---
errchk mimexpr, createDQexpr, t_mimexpr,expr_process
begin
    ldebug = false

    call smark(sp)
    
    call salloc (sci_ext, CARD_LEN, TY_CHAR)
    call salloc (var_ext, CARD_LEN, TY_CHAR)
    call salloc (dq_ext, CARD_LEN, TY_CHAR)
    call salloc (sci_expr, SZ_COMMAND, TY_CHAR)
    call salloc (var_expr, SZ_COMMAND, TY_CHAR)
    call salloc (dq_expr,  SZ_COMMAND, TY_CHAR)
    call salloc (ge_output, SZ_LINE, TY_CHAR)
    call salloc (oplist, SZ_LINE, TY_CHAR)
    call salloc (opval, SZ_LINE, TY_CHAR)
    call salloc (fname, SZ_PATHNAME, TY_CHAR)
    call salloc (mdfsource, SZ_PATHNAME, TY_CHAR)
    
    verbose  = clgetb("verbose")

    # put namesapce in common
    call strcpy(namespace, g_ns, SZ_LINE)
    
    call clgstr("sci_ext", Memc[sci_ext],CARD_LEN)
    if (allwhite(Memc[sci_ext])) {
        call sprintf(tmpstr, SZ_LINE, "sci_ext cannot be an empty string")
        call error(GEERR_BADPARM, tmpstr)
    }
    call strupr(Memc[sci_ext])
    call strcpy( Memc[sci_ext], opt_sci_ext, SZ_LINE) 
    # note for above: CARD_LEN doesn't work in place of SZ_LINE
    #   for opt_sci_ext error common where opt_sci_ext resides
    #   NOTE: this common var helps fmrelate refer to sci frames
    #   for VAR and DQ frame types
    call clgstr("var_ext", Memc[var_ext],CARD_LEN)
    if (allwhite(Memc[var_ext])) {
        call sprintf(tmpstr, SZ_LINE, "var_ext cannot be an empty string")
        call error(GEERR_BADPARM, tmpstr)
    }
    call strupr(Memc[var_ext])
    call clgstr("dq_ext", Memc[dq_ext],CARD_LEN)
    if (allwhite(Memc[dq_ext])) {
        call sprintf(tmpstr, SZ_LINE, "dq_ext cannot be an empty string")
        call error(GEERR_BADPARM, tmpstr)
    }
    call strupr(Memc[dq_ext])
    
	# get the gemexpr output parm
    call mkfullopname(namespace, "output", fullopnam, FULLOPNAM_LEN)
    call clgstr(fullopnam, Memc[ge_output], SZ_LINE)
    
    if (ldebug) {
        call printf("----\noutput = %s\n")
        call pargstr(Memc[ge_output])
        call flush(STDOUT)
    }

    ###NOTE!!!! This file check wants to be the first error... why?
    # because if any other error fires the caller will probably try to
    # delete the output thinking it's half-formed/malformed!
    # we don't like to delete the output BECAUSE it already exists
    if (imaccess(Memc[ge_output],0) == YES || access(Memc[ge_output],0,0) == YES)
    {
        call sprintf(tmpstr, SZ_LINE, "File, \"%s\" already exists")
        call pargstr(Memc[ge_output])
        call sfree(sp)
        call error(GEERR_OUTPUTEXISTS, tmpstr)
    }
    
    call mkfullopname(namespace, "sci_expr", fullopnam, FULLOPNAM_LEN)
    call clgstr(fullopnam, Memc[sci_expr],SZ_COMMAND)
    
    fl_vardq = clgetb("fl_vardq")
    if (fl_vardq) {
		call mkfullopname(namespace, "var_expr", fullopnam, FULLOPNAM_LEN)
        call clgstr(fullopnam, Memc[var_expr],SZ_COMMAND)
		call mkfullopname(namespace, "dq_expr", fullopnam, FULLOPNAM_LEN)
        call clgstr(fullopnam, Memc[dq_expr],SZ_COMMAND)
        call strupr(Memc[sci_ext])
        call strupr(Memc[var_ext])
        call strupr(Memc[dq_ext])
        if (false) {
            call printf("test - %s %s %s\n")
            call pargstr(Memc[sci_ext])
            call pargstr(Memc[var_ext])
            call pargstr(Memc[dq_ext])
            call flush(STDOUT)
        }
        
        if (streq(Memc[sci_ext], Memc[var_ext])){
            call sprintf(tmpstr, SZ_LINE, "PARMS: sci_ext(%s) cannot be var_ext(%s)")
            call pargstr(Memc[sci_ext])
            call pargstr(Memc[var_ext])
            call error(GEERR_PARMCOLLIDE, tmpstr)
        }
        if (streq(Memc[sci_ext],Memc[dq_ext])){
            call sprintf(tmpstr, SZ_LINE, "PARMS: sci_ext(%s) cannot be dq_ext(%s)")
            call pargstr(Memc[sci_ext])
            call pargstr(Memc[dq_ext])
            call error(GEERR_PARMCOLLIDE, tmpstr)
        }        
        if (streq(Memc[var_ext],Memc[dq_ext])){
            call sprintf(tmpstr, SZ_LINE, "PARMS: var_ext(%s) cannot be dq_ext(%s)")
            call pargstr(Memc[var_ext])
            call pargstr(Memc[dq_ext])
            call error(GEERR_PARMCOLLIDE, tmpstr)
        }
    }
    	
	# Prepare for going through the list of operands...
	
	# load expresion database.   
    st = NULL
    call mkfullopname(namespace, "exprdb", fullopnam, FULLOPNAM_LEN)
    call clgstr (fullopnam, Memc[fname], SZ_PATHNAME)
    call clpstr ("mimexprpars.exprdb", Memc[fname])
    if (strne (Memc[fname], "none"))
        st = ie_getexprdb (Memc[fname])
    
	# Parse the expression and generate a list of input operands.
    noperands = ie_getops (st, Memc[sci_expr], Memc[oplist], SZ_LINE)
    
    if (ldebug) {
        call printf("noperands=%d\n")
        call pargi(noperands)
        call flush(STDOUT)
    }

    #get refim
    call mkfullopname(namespace, "refim", fullopnam, FULLOPNAM_LEN)
    call clgstr(fullopnam, refimnam, SZ_FNAME)
    call clpstr("mimexprpars.refim", refimnam)
    if (false) {
        call printf("%s=%s\n")
        call pargstr(fullopnam)
        call pargstr(refimnam)
        call flush(STDOUT)
    }
    
    # Process the list of input operands and initialize each operand.
	# Note, these are the gemexpr operands so that mimexpr can get them
	# without asking for each extname
	
    opnam = oplist
    do i = 1, noperands 
	{
		if (Memc[opnam] == EOS)
            {
            call error (1, "malformed operand list")
        }
        
		# create fully qualified parameter name
		call mkfullopname(namespace, Memc[opnam], ge_opnam, GE_OPNAM_LEN)
        
        # check if operand is greater than "m", which is not allowed as
        #  this range is used for science frame reference when the work
        #  is handed to mimexpr
        if (Memc[opnam] > 'm')
        {
            call sprintf(tmpstr, SZ_LINE, "exprproc: \"%s\" is not a legal operator name\n")
            call pargstr(Memc[opnam])
            call log_err(GEERR_OPERANDRANGE, tmpstr )
            call error(GEERR_OPERANDRANGE, "OPERATOR OUT OF RANGE, 'm' is greatest legal operator name")
        }

        #pre-refactor ge_mapsciop[1] = Memc[opnam] + 13
        #pre-refactor ge_mapsciop[2] = EOS
	call sprintf(ge_mapsciop, GE_MAPSCIOP_LEN, "mimexprpars.%c")
	tc = Memc[opnam] + 13
	call pargc(tc) 
# note: 13 is added to provide access to science frames in operands n-z (user doesn't see
#       this... they use a[sci] rather than "n", but gemexpr does this
#       mapping before passing work on to mimexpr)
        if (false)
            {
            call printf("a-%s|\n")
            call pargstr(ge_opnam)
            call flush(STDOUT)
        }
		
	# get the parameter... fully qualified parameter name built above
        call clgstr (ge_opnam, Memc[opval], SZ_LINE)
		
        if (false)
        {
            call printf("opval-%s\n")
            call pargstr(Memc[opval])
            call flush(STDOUT)
            call printf("b-%s|\n")
            call pargstr(ge_opnam)
            call flush(STDOUT)
        }

		call sprintf(fullopnam, FULLOPNAM_LEN, "mimexprpars.%c")
		call pargc(Memc[opnam])
		
        if (ldebug) {
            call printf("fullopnam for put: %s\nge_mapsciop for put:%s\n")
            call pargstr(fullopnam)
            call pargstr(ge_mapsciop)
            call flush (STDOUT)
        }

        # original version: call clpstr (Memc[opnam], Memc[opval])
        call clpstr (fullopnam, Memc[opval])
        call clpstr (ge_mapsciop, Memc[opval])
#call printf("here3\n")
#call flush (STDOUT)
        if (ldebug)
        {
            call printf("(%s, %s)\n")
            call pargstr(Memc[opnam])
            call pargstr(Memc[opval])
            
            call printf("ge_opname,opval = %s,%s\n");
            call pargstr(ge_opnam[1])
            call pargstr(Memc[opval])
            call flush(STDOUT)
        }
        
        if (ldebug){
            call printf("refimnam = %s,%d\n")
            call pargstr(refimnam)
            call pargi(i)
            call flush(STDOUT)
        }

        if (streq("default", refimnam) && (i == 1)) {
            call strcpy(Memc[opval], Memc[mdfsource], SZ_PATHNAME )
        } else {  
#            if (Memc[opnam] == refimnam[1]) {
                if (streq(Memc[opnam], refimnam)) {
                # refim is the MDF source
                    call strcpy(Memc[opval], Memc[mdfsource], SZ_PATHNAME )
            }
        }

           # Get next operand name.
        while (Memc[opnam] != EOS) {
            opnam = opnam + 1
        }
        opnam = opnam + 1
    }
	
    if (ldebug) {
        call printf("Process Science Frames (EXTNAME=%s)\n")
        call pargstr(Memc[sci_ext])
        call flush(STDOUT)
    }
    
    iferr{
        # just in case... because there is an warning that we want to make a error
        call sciexpr_process(Memc[sci_expr])
        # SCI frames
        call clpstr("mimexprpars.expr", Memc[sci_expr])
        call clpstr("mimexprpars.extname", Memc[sci_ext])
        call clpstr("mimexprpars.output", Memc[ge_output])
        

        call clputi("status",0) 
        call t_mimexpr()
        tmpi = clgeti("status")
        if (tmpi == GEWARN_NOFRAMES) {
            # this is an error then, must have SCI frames
            call sprintf(tmpstr, SZ_LINE, "No Science Frames (EXTNAME = %s)")
            call pargstr(Memc[sci_ext])
            call error(GEERR_NOSCIENCE, tmpstr)
        }

        
        if (fl_vardq) {
            # VAR frames
            call expr_process(Memc[var_expr], EXPRTYPE_VAR)
            call clpstr("mimexprpars.expr", Memc[var_expr])
            call clpstr("mimexprpars.extname", Memc[var_ext])
            
            if (ldebug) {
                call printf("Process Variance Frames (EXTNAME=%s)\n")
                call pargstr(Memc[var_ext])
                call flush(STDOUT)
            }
            call t_mimexpr()
                
            # DQ frames
            if (ldebug) {
                call printf("After Variance\n")
                call flush(STDOUT)
            }
            if (strcmp("default", Memc[dq_expr]) == 0 ) {
                call createDQexpr(Memc[sci_expr], Memc[dq_expr], noperands)
            }
            if (ldebug) {
                call printf("dq_expr = %s, dq_ext=%s\n")
                call pargstr(Memc[dq_expr])
                call pargstr(Memc[dq_ext])
                call flush(STDOUT)
            }
            call expr_process(Memc[dq_expr], EXPRTYPE_DQ)
            call clpstr("mimexprpars.expr", Memc[dq_expr])
            call clpstr("mimexprpars.extname", Memc[dq_ext])

            if (ldebug) {
                call printf("Process Data Quality Frames (EXTNAME=%s)\n")
                call pargstr(Memc[dq_ext])
                call flush(STDOUT)
            }
            call t_mimexpr()
        }
    } then {
        tmpi = errget(tmpstr, SZ_LINE)
        call log_err(tmpi, tmpstr)
        call error(GEERR_GENERR, "mimexpr subsystem failed")
    }
    
    iferr {
        call gdpropagate(Memc[mdfsource],Memc[ge_output],"AUTO")
    } then {
        call error(GEERR_GENERR, "mimexpr subsystem failed because gdpropagate failed")
    }

        #set the NEXTEND to however many extensions are really in output
        call setNEXTEND(Memc[ge_output])
        
       call sfree(sp)
end

procedure createDQexpr(sci_expr, dq_expr, noperands)
char sci_expr[SZ_COMMAND]
char dq_expr[SZ_COMMAND]
int noperands

#------
pointer ie, io
pointer sp
char    opnam[30]
int i
int ie_optypes()
bool ldebug
bool wrotefirstop
begin
    ldebug = false
    
    # make sure dq_expr starts empty if were in here
    dq_expr[1]=EOS
    
    wrotefirstop = false
    call smark(sp)
    call salloc (ie, LEN_IMEXPR, TY_STRUCT)
    noperands = ie_optypes(sci_expr, ie)
    
    opnam[2]=EOS
    if (ldebug) {
        call printf ("dq_expr=%s\n")
        call pargstr(dq_expr)
        call flush(STDOUT)
    }

    do i = 1, noperands
    {
        io = IE_IMOP(ie,i)
        if (IO_TYPE(io) == IMAGE) {
            opnam[1] = 'a' - 1 + i
            if (wrotefirstop) {
                call strcat(" | ", dq_expr, SZ_COMMAND)
            }
            call strcat("int(", dq_expr, SZ_COMMAND)
            call strcat(opnam, dq_expr, SZ_COMMAND)
            call strcat("[DQ])", dq_expr, SZ_COMMAND)
            wrotefirstop = true
        }
    }
    if (ldebug) {    
        call printf ("dq_expr=%s\n")
        call pargstr(dq_expr)
        call flush(STDOUT)
    }
    call sfree(sp)
    return

end

# clean up expression if it has "[sci]" in it
# emit errors for other image values
procedure sciexpr_process(input)
char input[ARB]

# local functions below
bool done
int stridx(), strncmp()
int ind
int len, strlen()
bool ldebug

begin
    ldebug = false

    if (ldebug) {
        call printf ("sciexpr_process:input expr = %s\n")
        call pargstr(input)
        call flush(STDOUT)
    }

    done = false
    ind = 1
    len = strlen(input)
    while (done == false) {
        ind = stridx("[", input[1])

        #skip "[0]" syntax
        while (strncmp(input[ind],"[0]",3)==0) {
            ind = stridx("[",input[ind+3])
        }
        
        if (ind == 0) {
            done = true
        } else {
            if (ldebug) {
                call printf("B:%c ind=%d input=%s\n")
                call pargc(input[ind])
                call pargi(ind)
                call pargstr(input[1])
                call flush(STDOUT)
            }

            if (   (strncmp(input[ind],"[sci]",5) == 0) 
                || (strncmp(input[ind],"[SCI]",5) == 0)) {
                # then this is the right place
                # mask out the qualifier with no null
                input[ind] = ' '
                input[ind+1] = ' '
                input[ind+2] = ' '
                input[ind+3] = ' '
                input[ind+4] = ' '
                call amovc(input[ind+5], input[ind], len-ind)
                # map the letter
            } else if (  (strncmp(input[ind],"[VAR]",5) == 0)
                    || (strncmp(input[ind],"[var]",5) == 0)) {
                call error(GEERR_INVALIDEXPR, "[VAR] not allowed in sci_expr")
            } else if (  (strncmp(input[ind],"[DQ]",4) == 0)
                    || (strncmp(input[ind],"[dq]",4) == 0)) {
                call error(GEERR_INVALIDEXPR, "[DQ] not allowed in sci_expr")
            } else {
                call error(GEERR_INVALIDEXPR, "invalid use of \"[\" or \"]\"")
            }
        }    
    }
    if(ldebug) {
        call printf("\nsciexpr_process output expr=%s\n")
        call pargstr(input)
        call flush(STDOUT)
    }
end


procedure expr_process(input,type)
char input[ARB]
int type
# local functions below
bool done
int stridx(), strncmp(), stridxs()
int ind, tind,cind
int len, strlen()
bool ldebug
int test

begin
       ldebug = false

    if (ldebug) {
        call printf ("input expr = %s\n")
        call pargstr(input)
        call flush(STDOUT)
    }

    # this first loops transforms letters that are not qualified 
    # with [...], a second loop takes care of those that are
    done = false
    ind = 1
    tind = ind
    cind = ind
    len = strlen(input)
    test =0
    while (done == false) {

        test = test+1
        ind = stridxs("abcdefghijklm",input[cind])
        if (ind == 0) {
            done = true
        } else {
            cind = cind+ind-1
            tind = cind+1
            
            if (ldebug) {
                call printf("A(%d):%c cind=%d %d %d input=%s\n")
                call pargi(len)
                call pargc(input[cind])
                call pargi(cind)
                call pargi(ind)
                call pargi(tind)
                call pargstr(input)
                call flush(STDOUT)
            }

            if (input[tind] == '.') {
                # then transform letter
                input[cind] = input[cind]+13            
            } else {
                # If the next character is not "]", "[", an alphabetic
                # character and "(", then it's an input operand
                # The '(' is to allow macro functions
                # The remaining two allow exponentials
                if (   (input[tind] != '[') 
                    && !IS_ALPHA(input[tind])
                    && (input[tind] != ']')
                    && (input[tind] != '(')
                    && !(IS_DIGIT(input[cind-1]) && input[tind] == '+' && IS_DIGIT(input[tind+1]))
                    && !(IS_DIGIT(input[cind-1]) && input[tind] == '-' && IS_DIGIT(input[tind+1]))
                ) {
                    input[cind] = input[cind]+13
                }          
            }

            cind = cind+1
        }
        if (cind>=len) done=true
    }
    if (ldebug) {
        call printf("Pass One Complete...\n")
        call flush(STDOUT)
    }
    
    done = false
    ind = 1
    while (done == false) {
        ind = stridx("[", input[1])

        #skip "[0]" syntax
        while (strncmp(input[ind],"[0]",3)==0) {
            ind = stridx("[",input[ind+3])
        }
        
        if (ind == 0) {
            done = true
        } else {
            if (ldebug) {
                call printf("B:%c ind=%d input=%s\n")
                call pargc(input[ind])
                call pargi(ind)
                call pargstr(input[1])
                call flush(STDOUT)
            }

            if (   (strncmp(input[ind],"[sci]",5) == 0) 
                || (strncmp(input[ind],"[SCI]",5) == 0)) {
                # then this is the right place
                # mask out the qualifier with no null
                input[ind] = ' '
                input[ind+1] = ' '
                input[ind+2] = ' '
                input[ind+3] = ' '
                input[ind+4] = ' '
                call amovc(input[ind+5], input[ind], len-ind)
                # map the letter
                if (false) { #if (input[ind-1]+13 > 'z') {
                    call log_err(GEERR_EXPRPROC,"ILLEGAL OPERAND OUT OF RANGE!")
                    call error(GEERR_EXPRPROC,"ILLEGAL OPERAND OUT OF RANGE!")
                }
                input[ind-1] = input[ind-1]+13
            } else if (  (strncmp(input[ind],"[VAR]",5) == 0)
                    || (strncmp(input[ind],"[var]",5) == 0)) {
                if (type != EXPRTYPE_VAR) {
                    call error(GEERR_INVALIDEXPR, "[VAR] qualifier allowed only in var_expr")
                }
                call amovc(input[ind+5],input[ind], len-ind)
            } else if (  (strncmp(input[ind],"[DQ]",4) == 0)
                    || (strncmp(input[ind],"[dq]",4) == 0)) {
                if (type != EXPRTYPE_DQ) {
                    call error(GEERR_INVALIDEXPR, "[DQ] qualifier allowed only in dq_expr")
                }
                call amovc(input[ind+4],input[ind], len-ind)
            } else {
                call error(GEERR_INVALIDEXPR, "invalid use of '[' or ']'")
            }
        }    
    }
    if(ldebug) {
        call printf("\nexpr_process output expr=%s\n")
        call pargstr(input)
        call flush(STDOUT)
    }
end

int procedure ie_optypes(expr, ie)
char expr[ARB]  # input expression (mimexpr ready, e.g. after expr_process(..)
pointer ie      # pointer to fill with operand information, should be
                       #  LEN_IMEXPR sized (see gemexpr.h)

# local variable below
int noperands
pointer sp,opnam,oplist,opval

pointer io,ip,o
pointer st
pointer fname,imname, section, cluster

bool ldebug
int dtype, nchars, i
double dval
int ctod()
int lexnum(),stridxs()
int int()
bool streq()
# pointer ie_getexprdb()
int ie_getops()
char fullopnam[SZ_LINE]
include "nscommon.h"

define  numeric_ 91
define  image_ 92
define  extiter_ 93

begin
    
    ldebug = false
    if (ldebug) {
        call printf("ie_optypes():expr=%s\n")
        call pargstr(expr[1])
        call flush(STDOUT)
    }
    
    call smark  (sp)
    
    call aclri (Memi[ie], LEN_IMEXPR)
    
    call salloc (fname, SZ_PATHNAME, TY_CHAR)
    call salloc (imname, SZ_PATHNAME, TY_CHAR)
    call salloc (section, SZ_FNAME, TY_CHAR)
    call salloc (cluster, SZ_FNAME, TY_CHAR)
    call salloc (oplist, SZ_LINE, TY_CHAR)
    call salloc (opval, SZ_LINE, TY_CHAR)
    
#    call printf("two\n")
#    call flush(STDOUT)
    
       # Load the expression database, if any.
    st = NULL
#       call clgstr ("exprdb", Memc[fname], SZ_PATHNAME)
#       if (strne (Memc[fname], "none"))
#           st = ie_getexprdb (Memc[fname])
    IE_ST(ie) = st

#       call printf("three\n")
#       call flush(STDOUT)

# Parse the expression and generate a list of input operands.
    noperands = ie_getops (st, expr, Memc[oplist], SZ_LINE)
    IE_NOPERANDS(ie) = noperands

    # Process the list of input operands and initialize each operand.
    # This means fetch the value of the operand from the CL, determine
    # the operand type, and initialize the image operand descriptor.
    # The operand list is returned as a sequence of EOS delimited strings.

    opnam = oplist
    do i = 1, noperands {
        io = IE_IMOP(ie,i)
        if (Memc[opnam] == EOS)
            call error (1, "ie_fill(): invalid operand list")

	    call mkfullopname(g_ns, Memc[opnam], fullopnam, FULLOPNAM_LEN)        
        call clgstr (fullopnam, Memc[opval], SZ_LINE)
        
        if (ldebug)
            {
            call printf("t_g(): opname, opval = %s,%s\n")
            call pargstr(Memc[opnam])
            call pargstr(Memc[opval])
            call flush(STDOUT)
        }
        IO_OPNAME(io) = Memc[opnam]
        ip = opval
        
# Initialize the input operand; these values are overwritten below.
        o = IO_OP(io)
        call aclri (Memi[o], LEN_OPERAND)
        
        if (Memc[ip] == '.' && (Memc[ip+1] == EOS || Memc[ip+1] == '['))
            {
# A "." is shorthand for the last output image.
            call strcpy (Memc[ip+1], Memc[section], SZ_FNAME)
            call clgstr ("lastout", Memc[opval], SZ_LINE)
            call strcat (Memc[section], Memc[opval], SZ_LINE)
            goto image_
        }
        else if (IS_LOWER(Memc[ip]) && Memc[ip+1] == '.')
        {
               # "a.foo" refers to parameter foo of image A.  Mark this as
               # a parameter operand for now, and patch it up later.
            
            IO_TYPE(io) = PARAMETER
            IO_DATA(io) = ip
            
            if(ldebug) {
                call printf("PARAMETER\n")
                call flush(STDOUT)
            }
        }
        else if (ctod (Memc, ip, dval) > 0)
        {
            if (Memc[ip] != EOS)
                goto image_
            
# A numeric constant.
            numeric_        IO_TYPE(io) = NUMERIC
            if(ldebug) {
                call printf("NUMERIC\n")
                call flush(STDOUT)
            }
            
            ip = opval
            switch (lexnum (Memc, ip, nchars))
            {
                case LEX_REAL:
                dtype = TY_REAL
                if (stridxs("dD",Memc[opval]) > 0 || nchars > NDIGITS_RP+3)
                    dtype = TY_DOUBLE
                O_TYPE(o) = dtype
                if (dtype == TY_REAL)
                    O_VALR(o) = dval
                else
                O_VALD(o) = dval
                
                default:
                O_TYPE(o) = TY_INT
                O_LEN(o)  = 0
                O_VALI(o) = int(dval)
            }
        }
        else
        {
# Anything else is assumed to be an image name.
            image_
            if(ldebug)
                {
                call printf("IMAGE\n")
                call flush(STDOUT)
            }
            
            ip = opval
            call imgimage (Memc[ip], Memc[fname], SZ_PATHNAME)
            if (streq (Memc[fname], Memc[imname])) {
                call error (2, "input and output images cannot be the same")
            }
            call imgimage(Memc[ip], Memc[cluster], SZ_FNAME)
            IO_MEF(io) = NULL
            
            IO_TYPE(io) = IMAGE
            
            # If one dimensional image read in data and be done with it.
            # Get next operand name.
            while (Memc[opnam] != EOS)
            {
                opnam = opnam + 1
            }
            opnam = opnam + 1
        }
    }
    call sfree(sp)
    
    return noperands
end

procedure mkfullopname(namespace, opnam, out_fullname, ofnlen)
char namespace[ARB]
char opnam[ARB]
char out_fullname[ARB]
int ofnlen #out_fullname buffer length

# --- functions
bool streq()

begin
	if ( streq (namespace, ""))
	{
		call strcpy (opnam, out_fullname, ofnlen)
	}	
	else
	{
		call sprintf(out_fullname, ofnlen, "%s.%s")
		call pargstr(namespace)
		call pargstr(opnam)
	}
		
	return
end
    
bool procedure  allwhite(str)
char str[ARB]
int len, i    
int strlen()
bool ldebug

begin
        
        ldebug=false

        len = strlen(str)
        
     
        for (i = 1; i <= len ; i= i+1) {
             if (ldebug) {   
                call printf("str[%d]=%c\n")
                call pargi(i)
                call pargc(str[i])
            }
            if (! IS_WHITE(str[i]))
                break;   
       }
 
        if (ldebug) {
            call printf("<- allwhite(%s), len=%d, i=%d\n")
            call pargstr(str)
            call pargi(len)
            call pargi(i)
        }
    
        if (i>len) return true
        else return false
        
end
    

