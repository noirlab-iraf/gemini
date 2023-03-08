# Copyright(c) 2005 Association of Universities for Research in Astronomy, Inc.

# Test routines for the GEMINI library

task    demostrtools         = t_demostrtools

# DEMOSTRTOOLS -- Demo/Test of the STRTOOLS functions
#
# The parameter file for DEMOSTRTOOLS is:
#
# demo,s,a,"all","|all|numrecog|",,"Which demo?"
# status,i,h,0,,,"Exit status (0=good)"
#
# To compile and run DEMOSTRTOOLS:
#   cl> gemini
#   ge> cd gemini/lib/gemini
#   ge> mkpkg -p gemini zzdebug
#   ge> cd
#   ge> task demostrtools="gemini$lib/gemini/zzdebug.e"
#   ge> demostrtools

# Author: Kathleen Labrie 16-May-2005

procedure t_demostrtools ()

# Local variables for task parameters
char	l_demo[SZ_LINE]
int     l_status

# Gemini functions
int     d_numrecog()

# IRAF functions
int     strcmp()

begin
        l_status = 0
        
        # Get task parameter values
        call clgstr ("demo", l_demo, SZ_LINE)
        
        if ( strcmp (l_demo, "numrecog") == 0 )
            l_status = l_status + d_numrecog()
        else {      # must be 'all'
            l_status = l_status + d_numrecog()
        }
        
        call clputi ("status", l_status)
        return
end


# D_NUMRECOG -- Function to test the number recognition tools
#       status = d_numrecog ()
#
#       status      : Exit status       [Return value, (int)]

int procedure d_numrecog ()

char    str[SZ_LINE,20]
int     i, l_status

# Gemini functions
bool    g_isDecimal(), g_isInteger(), g_isReal(), g_isExpNumber()
bool    g_isOctal(), g_isHex()

begin
        l_status = 0
        
        call strcpy("a true string", str[1,1], SZ_LINE)
        call strcpy("59", str[1,2], SZ_LINE)
        call strcpy("-432", str[1,3], SZ_LINE)
        call strcpy("3.14159", str[1,4], SZ_LINE)
        call strcpy("-98.2", str[1,5], SZ_LINE)
        call strcpy("6.626e-34", str[1,6], SZ_LINE)
        call strcpy("-0.5e-10", str[1,7], SZ_LINE)
        call strcpy("1e3", str[1,8], SZ_LINE)
        call strcpy(".5", str[1,9], SZ_LINE)
        call strcpy("2005-05-13", str[1,10], SZ_LINE)
        call strcpy("834.43.32", str[1,11], SZ_LINE)
        call strcpy("INDEF", str[1,12], SZ_LINE)
        call strcpy("45b", str[1,13], SZ_LINE)
        call strcpy("-323B", str[1,14], SZ_LINE)
        call strcpy("0-b", str[1,15], SZ_LINE)
        call strcpy("-398B", str[1,16], SZ_LINE)
        call strcpy("0ffx", str[1,17], SZ_LINE)
        call strcpy("F8a0eX", str[1,18], SZ_LINE)
        call strcpy("-578x", str[1,19], SZ_LINE)
        call strcpy("4fgx", str[1,20], SZ_LINE)
        
        for (i=1; i<=20; i=i+1) {
            call printf ("'%s'  ")
                call pargstr (str[1,i])
            if (g_isDecimal(str[1,i]))
                call printf ("Decimal, ")
            else
                call printf ("NOT Decimal, ")

            if (g_isInteger(str[1,i]))
                call printf ("Integer, ")
            else
                call printf ("NOT Integer, ")

            if (g_isReal(str[1,i]))
                call printf ("Real, ")
            else
                call printf ("NOT Real, ")

            if (g_isExpNumber(str[1,i]))
                call printf ("Exponential, ")
            else
                call printf ("NOT Exponential, ")

            if (g_isHex(str[1,i]))
                call printf ("Hex, ")
            else
                call printf ("NOT Hex, ")
            
            if (g_isOctal(str[1,i]))
                call printf ("Octal\n")
            else
                call printf ("NOT Octal\n")
        }
        
        # might want to print that stuff to a file and compare the
        # output with the expected one.  Any difference would show
        # a problem with the code.
        
        return (l_status)
end
