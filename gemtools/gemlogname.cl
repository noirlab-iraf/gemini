# Copyright(c) 2002-2009 Association of Universities for Research in Astronomy, Inc.

procedure gemlogname()

# Determine log file name, by parsing param string / package setting
#
# Parameters:
# logpar   - input: log name parameter from calling task
#            output: same as input if non-blank, otherwise package setting
# package  - input: name of current Gemini sub-package
#            output: same as input
# logname  - output: valid log-file name to be used (blank if undetermined)
# status   - output: result of execution; see below
#
# Status:
# 0 - logfile name determined successfully from task or package param
# 1 - logfile name unspecified by task or package; using package default
# 2 - bad logfile name; using package setting or default
# 3 - logfile name bad or unspecified and no valid package name
#
# NB. The check for invalid chars in the filename may not be exhaustive
#
# Version  Sept 20, 2002 JT v1.4 release

string logpar  {"", prompt="Input log name"}
string package {"", prompt="Name of current (sub-) package"}
string logname {"", prompt="Actual log name to use"}
int status {0, prompt="Exit status"}

begin

string l_logpar, l_package, l_logname, pkgset, tstr
int l_status
bool bad_logpar, bad_package, bad_pkgset

# Read the parameters:
l_logpar = logpar
l_package = package

# Default state:
bad_logpar = no
bad_package = no
bad_pkgset = no

tstr = ""

# Strip any whitespace from package name:	
print(l_package) | scanf("%s %s", l_package, tstr)

# Determine package setting:
#   (can't think how to get value <pkg>.logfile with <pkg> as a variable)
if (tstr != "") {
  pkgset = ""
  bad_package = yes
}
else if (l_package=="gmos")
  pkgset = gmos.logfile
else if (l_package=="niri")
  pkgset = niri.logfile
else if (l_package=="agcam")
  pkgset = agcam.logfile
else if (l_package=="flamingos")
  pkgset = flamingos.logfile
else if (l_package=="oscir")
  pkgset = oscir.logfile
else if (l_package=="quirc")
  pkgset = quirc.logfile
else {
  pkgset = ""
  bad_package = yes
}

# Strip any whitespace from package setting:
print(pkgset) | scanf("%s %s", pkgset, tstr)

# Check validity:
if (tstr != "" || stridx("~!$&*?<>'\"", pkgset) != 0 ||
    substr(pkgset, 1, 1) == "-")
  bad_pkgset = yes

# Strip any whitespace from log name:
print(l_logpar) | scanf("%s %s", l_logpar, tstr)

# Check validity:
if (tstr != "" || stridx("~!$&*?<>'\"", l_logpar) != 0 ||
    substr(l_logpar, 1, 1) == "-")
  bad_logpar = yes

# Determine name to use:
if (l_logpar == "") {
  l_logpar = pkgset
  if (bad_package) {
    l_logname = ""
    l_status = 3
  }
  else if (pkgset == "") {
    l_logname = l_package//".log"
    l_status = 1
  }
  else if (bad_pkgset) {
    l_logname = l_package//".log"
    l_status = 2
  }
  else {
    l_logname = pkgset
    l_status = 0
  }
}
else if (bad_logpar) {
  if (bad_package) {
    l_logname = ""
    l_status = 3
  }
  else if (pkgset == "" || bad_pkgset) {
    l_logname = l_package//".log"
    l_status = 2
  }
  else {
    l_logname = pkgset
    l_status = 2
  }
}
else {
  l_logname = l_logpar
  l_status = 0
}

# Set output parameters:
logpar = l_logpar
package = l_package  # so IRAF remembers it
logname = l_logname
status = l_status

end
