# Copyright(c) 2005-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen labrie 17-May-2005

# T_GTISNUMBER -- Check if a string is a number, and if so, which type.

procedure t_gisnumber ()

#char   instring[SZ_LINE]           #I String to check
#char   ttest[SZ_LINE]              #I Type to test for
#bool   fl_istype                   #O Result of the test
#bool   verbose                     #I Verbose
#int    status                      #O Exit status

char    l_instring[SZ_LINE], l_ttest[SZ_LINE]
bool    l_fl_istype, l_verbose
int     l_status

# GEMINI functions
bool    g_isDecimal(), g_isInteger(), g_isReal(), g_isExpNumber()
bool    g_isHex(), g_isOctal()

# IRAF functions
int     strcmp()
bool    clgetb()

begin
        l_status = 0
        l_fl_istype = FALSE
        
        # Get task parameter values
        call clgstr ("instring", l_instring, SZ_LINE)
        call clgstr ("ttest", l_ttest, SZ_LINE)
        l_verbose = clgetb ("verbose")
        
        if ( strcmp (l_ttest, "decimal") == 0 )
            l_fl_istype = g_isDecimal (l_instring)
        else if ( strcmp (l_ttest, "integer") == 0 )
            l_fl_istype = g_isInteger (l_instring)
        else if ( strcmp (l_ttest, "real") == 0 )
            l_fl_istype = g_isReal (l_instring)
        else if ( strcmp (l_ttest, "exponential") == 0 )
            l_fl_istype = g_isExpNumber (l_instring)
        else if ( strcmp (l_ttest, "hex") == 0 )
            l_fl_istype = g_isHex (l_instring)
        else if ( strcmp (l_ttest, "octal") == 0 )
            l_fl_istype = g_isOctal (l_instring)
        else {
            l_status = 99
            call printf ("GEMISNUMBER ERROR 99 Internal Error\n")
        }
      
       if (l_verbose) {
            call printf ("%s\n")
                call pargb (l_fl_istype)
        }
      
        call clputb ("fl_istype", l_fl_istype)
        call clputi ("status", l_status)
        return
end
