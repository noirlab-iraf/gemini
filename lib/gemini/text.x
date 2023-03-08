# Copyright(c) 2002-2005 Association of Universities for Research in Astronomy, Inc.
#
# Some subroutines for manipulating text strings, supporting the gemini tasks
#
# Version     Sep  5, 2002  JT
# ----

include <ctype.h>


int procedure gte_blankpar(str)

# Determine whether input string is blank except for quotation
# (returns YES or NO)

char str[ARB]                           # I: input string

int n

begin

  for (n=1; str[n]!=EOS; n=n+1)
    if (!IS_WHITE(str[n]) && str[n]!='\'' && str[n]!='\"') return NO

  return YES
  
end


int procedure gte_isblank(str)

# Determine whether input string contains only white space
# (returns YES or NO)

char str[ARB]                           # I: input string

int n

begin

  for (n=1; str[n]!=EOS; n=n+1) if (!IS_WHITE(str[n])) return NO
  return YES
  
end


procedure gte_unpad(str, outstr, maxlen)

# Remove any leading or trailing whitespace from a string

# Str and outstr can be the same.

char str[ARB]                           # I: string to strip padding from
char outstr[ARB]                        # O: string without padding
int maxlen                              # I: maximum output length

int n

int gstrcpy()

begin
  
  # Skip leading whitespace:
  for (n=1; IS_WHITE(str[n]); n=n+1);
  
  # Copy remaining input to output:
  n = gstrcpy(str[n], outstr, maxlen)

  # Remove trailing whitespace:
  for (; n>0 && IS_WHITE(outstr[n]); n=n-1);
  outstr[n+1] = EOS
  
end


int procedure gte_indict(instr, dict)

# Determine whether input string exists in a dictionary string
# (returns YES or NO)

# This differs from the IRAF routine strdic() in that the match must be
# exact, including capitalization and white space. A blank string never
# matches the dictionary.

char instr[ARB]                         # I: input string
char dict[ARB]                          # I: dictionary str to check against

int nd, nw

int gte_isblank()

begin

  # Indicate no match if input string is blank:
  if (gte_isblank(instr)==YES) return NO
  
  # Loop through dictionary entries:
  for (nd=1; dict[nd]!=EOS; nd=nd+1) {
  
    # Match word against current dict entry:
    for (nw=1; instr[nw]==dict[nd] && dict[nd]!='|' && dict[nd]!=EOS; nw=nw+1)
      nd = nd+1

    # If word chars match dict & entry has ended, indicate match:
    if (instr[nw]==EOS && (dict[nd]=='|' || dict[nd]==EOS)) return YES
  
    # Otherwise, forward to end of entry:
    for (; dict[nd]!='|' && dict[nd]!=EOS; nd=nd+1);
  
    # If dictionary finished, indicate no match:
    if (dict[nd]==EOS) return NO
    
  }

  # No match if last entry blank:
  return NO

end
