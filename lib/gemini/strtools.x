# Copyright(c) 2005-2009 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie   16-May-2005

include <ctype.h>

.help
.nf
This file contains functions to verify whether a string variable actually
contains a string, a integer, a real, etc.

      g_isDecimal - Is the value numeric?
      g_isInteger - Is the value an integer?
         g_isReal - Is the value a floating point?
    g_isExpNumber - Is the value in exponential notation?
          g_isHex - Is the value in hexadecimal notation?
        g_isOctal - Is the value in octal notation?

G_ISDECIMAL -- Returns true if the value in a string variable is decimal.
        bool = g_isDecimal( str )
        
        bool        : Test result           [return value, (bool)]
        str         : String to test        [input, (string)]

G_ISINTEGER -- Returns true if the value in a string variable is an integer.
        bool = g_isInteger( str )
        
        bool        : Test result           [return value, (bool)]
        str         : String to test        [input, (string)]

G_ISREAL -- Returns true if the value in a string variable is a floating point
          number.
        bool = g_isReal( str )
        
        bool        : Test result           [return value, (bool)]
        str         : String to test        [input, (string)]

G_ISEXPNUMBER -- Returns true if the value in a string variable is a floating
               point number written using the exponential format.
        bool = g_isExpNumber( str )
        
        bool        : Test result           [return value, (bool)]
        str         : String to test        [input, (string)]

G_ISHEX -- Returns true if the value in a string variable is a hexadecimal number
        bool = g_isHex( str )
        
        bool        : Test result           [return value, (bool)]
        str         : String to test        [input, (string)]

G_ISOCTAL -- Returns true if the value in a string variable is an octal number
        bool = g_isOctal( str )
        
        bool        : Test result           [return value, (bool)]
        str         : String to test        [input, (string)]
.fi
.endhelp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

#G_ISDECIMAL -- Returns true if the value in a string variable is numeric.
#        bool = g_isDecimal( str )
#        
#        bool        : Test result           [return value, (bool)]
#        str         : String to test        [input, (string)]

bool procedure g_isDecimal (str)

char    str[ARB]        #I String to test

bool    isDecimal

int     Ecount, Pcount
bool    isINDEF

bool    g_dectest()

begin
        isDecimal = g_dectest (str, Ecount, Pcount, isINDEF)
            
        return (isDecimal)
end


#G_ISINTEGER -- Returns true if the value in a string variable is an integer.
#        bool = g_isInteger( str )
#        
#        bool        : Test result           [return value, (bool)]
#        str         : String to test        [input, (string)]

bool procedure g_isInteger (str)

char    str[ARB]        #I String to test

bool    isInteger

int     Ecount, Pcount
bool    isdec, isINDEF

bool    g_dectest

begin
        isdec = g_dectest (str, Ecount, Pcount, isINDEF)

        if (! isdec)
            isInteger=FALSE
        else if (isINDEF)
            isInteger=TRUE
        else if ((Ecount != 0) || (Pcount != 0))
            isInteger=FALSE
        else
            isInteger=TRUE

        return (isInteger)
end


#G_ISREAL -- Returns true if the value in a string variable is floating point
#            number.
#        bool = g_isReal( str )
#        
#        bool        : Test result           [return value, (bool)]
#        str         : String to test        [input, (string)]

bool procedure g_isReal (str)

char    str[ARB]        #I String to test

bool    isReal

int     Ecount,Pcount
bool    isdec, isINDEF

bool    g_dectest()

begin
        isdec = g_dectest (str, Ecount, Pcount, isINDEF)
        
        if (! isdec)                             # not a number
            isReal=FALSE
        else if ((Pcount == 0) && (Ecount == 0) && (!isINDEF)) # this is an integer
            isReal=FALSE
        else
            isReal=TRUE

        return (isReal)
end


#G_ISEXPNUMBER -- Returns true if the value in a string variable is a floating
#                 point number written using the exponential format.
#        bool = g_isExpNumber( str )
#        
#        bool        : Test result           [return value, (bool)]
#        str         : String to test        [input, (string)]

bool procedure g_isExpNumber (str)

char    str[ARB]        #I String to test

bool    isExpNumber

int     Ecount, Pcount
bool    isdec, isINDEF

bool    g_dectest()

begin
        isdec = g_dectest( str, Ecount, Pcount, isINDEF )
        
        if (! isdec)                            # not a number
            isExpNumber=FALSE
        else if ((Ecount == 0) && (!isINDEF))   # not exponential and not INDEF
            isExpNumber=FALSE
        else
            isExpNumber=TRUE

        return (isExpNumber)
end


#G_DECTEST -- Returns true if the value in a string variable is not a string.
#             Also returns the info about format.
#        bool = g_dectest( str, Ecount, Pcount, isINDEF )
#        
#        test        : Test result           [return value, (bool)]
#        str         : String to test        [input, (string)]
#        Ecount      : Number of 'e' or 'E'  [output, (int)]
#        Pcount      : Number of '.'         [output, (int)]
#        isINDEF     : Flag for INDEF values    [output, (bool)]

bool procedure g_dectest (str, Ecount, Pcount, isINDEF)

char    str[ARB]        #I String to test
int     Ecount          #O Number of 'e' or 'E' (1 or 0 in numeric values)
int     Pcount          #O Number of '.'  (1 or 0 in numeric values)
bool    isINDEF         #O Flag for INDEF values

bool    isdec

char    c
int     len, i
int     Epos

# Gemini functions
bool    g_whitespace()

# IRAF functions
int     strlen(), strncmp()

begin
        isdec = TRUE
        Ecount = 0
        Epos = 0
        Pcount = 0
        isINDEF = FALSE

        if ( g_whitespace (str) ) {     # An empty string
        
            isdec = FALSE
            
        } else {
        
            len = strlen (str)
            for (i=1; i<=len; i=i+1) {

                c = str[i]
                if ( !IS_DIGIT(c) ) {
                    switch (c) {
                    case '-':
                        if (i == len)             #a character must follow
                            isdec=FALSE
                        else if (i != 1) {        #not for negative value
                            if (Ecount == 0)      #not exponential notation
                                isdec=FALSE
                            else if (Epos != i-1) #'-' not after e or E
                                isdec=FALSE
                        }
                    case '.':
                        if (Pcount != 0)      #only one '.' allowed
                            isdec=FALSE
                        else if ((Ecount!=0) && (i > Epos))  #no '.' allowed in the exp part
                            isdec=FALSE
                        else
                            Pcount = Pcount+1
                    case 'e','E':
                        if (Ecount != 0)      #only one 'e' or 'E' allowed
                            isdec=FALSE
                        else if (i == len)  #a character must follow
                            isdec=FALSE
                        else {
                            Ecount = Ecount+1
                            Epos = i
                        }
                    case 'I':   # check for an INDEF value
                        if (i != 1)
                            isdec=FALSE
                        else if ( strncmp( str, "INDEF", len ) != 0 )
                            isdec=FALSE
                        else {
                            isINDEF = TRUE
                            return (isdec)   #is INDEF => isdec=TRUE
                        }
                    default:
                        isdec=FALSE
                    }               
                }

                if (isdec == FALSE)
                    return(isdec)
            }
        }

        return (isdec)
end


#G_ISHEX -- Returns true if the value in a string variable is a hexadecimal number
#        bool = g_isHex( str )
#        
#        bool        : Test result           [return value, (bool)]
#        str         : String to test        [input, (string)]

bool procedure g_isHex (str)

char    str[ARB]        #I String to test

bool    isHex

char    c
int     len, i

# Gemini functions
bool    g_whitespace()

# IRAF functions
int     strlen(), strncmp()

begin
        isHex = TRUE

        if ( g_whitespace (str) ) {     # An empty string
        
            isHex = FALSE
        
        } else {
        
            len = strlen (str)

            # Hexadecimal numbers _must_ have their last character be 'x' or 'X'
            if ((str[len] != 'x') && (str[len] != 'X')) {
                if ( strncmp (str, "INDEF", len) == 0)      #is INDEF
                    isHex=TRUE
                else
                    isHex = FALSE
                return (isHex)
            }

            # isHex assumed to be true, unless one of the conditions below
            # show otherwise.

            for (i=1; i<=len-1; i=i+1) {     #we've already checked the last char.

                c = str[i]
                if ( !IS_DIGIT(c) && !( c>='a' && c<='f') && !(c>='A' && c<='F')) {
                    if (c == '-') {
                        if (i == len-1)         #a valid character must follow
                            isHex=FALSE
                        else if (i != 1)        # not for negative value
                            isHex=FALSE
                    } else
                        isHex=FALSE
                }

                if (isHex == FALSE)
                    return(isHex)          
            }
        }

        return (isHex)
end


#G_ISOCTAL -- Returns true if the value in a string variable is an octal number
#        bool = g_isOctal( str )
#        
#        bool        : Test result           [return value, (bool)]
#        str         : String to test        [input, (string)]

bool procedure g_isOctal (str)

char    str[ARB]        #I String to test

bool    isOctal

char    c
int     len, i

# Gemini functions
bool    g_whitespace()

# IRAF functions
int     strlen(), strncmp()

begin
        isOctal = TRUE
        
        if ( g_whitespace(str) ) {      # An empty string
        
            isOctal = FALSE
        
        } else {
        
            len = strlen (str)

            # Octal numbers _must_ have their last character be 'b' or 'B'
            if ((str[len] != 'b') && (str[len] != 'B')) {
                if ( strncmp (str, "INDEF", len) == 0)      #is INDEF
                    isOctal=TRUE
                else
                    isOctal = FALSE
                return (isOctal)
            }

            # isOctal assumed to be true, unless one of the conditions below
            # show otherwise.

            for (i=1; i<=len-1; i=i+1) {     #we've already checked the last char.

                c = str[i]
                if ( ! (c >= '0' && c <= '7') ) {
                    if (c == '-') {
                        if (i == len-1)         #a valid character must follow
                            isOctal=FALSE
                        else if (i != 1)        # not for negative value
                            isOctal=FALSE
                    } else
                        isOctal=FALSE
                }

                if (isOctal == FALSE)
                    return(isOctal)          
            }
        }

        return (isOctal)
end
