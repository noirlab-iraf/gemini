# Copyright(c) 2000-2006 Association of Universities for Research in Astronomy, Inc.

procedure printlog (text,logfile,verbose)

# Script to send text to logfile and terminal
# 
# Version  April 24, 2000  BM
#          Sept 20,  2002  BM v1.4 release

string text      {"",prompt="Text string to log"}
string logfile   {"",prompt="Name of log file"}
bool   verbose   {yes,prompt="Output to screen"}

begin

string l_text,l_logfile
bool l_verbose

l_text=text
l_logfile=logfile
l_verbose=verbose

if (l_logfile == "STDOUT") {
   l_logfile = ""
   l_verbose=yes
}

if ((substr(l_text,1,5) == "ERROR") || (substr(l_text,1,7) == "WARNING")) {
    l_verbose=yes
}
if (l_logfile != "" && l_logfile != " ") {
    print(l_text, >> l_logfile)
}
if (l_verbose) {
    print(l_text)
}

end
