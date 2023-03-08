#{ Copyright(c) 2000-2017 Association of Universities for Research in Astronomy, Inc.
#
# Package script for the gemini package
#
# Version: Sept 14, 2002 BR,IJ  Release v1.4
#          Jan   9, 2004 KL     Release v1.5
#          Apr  19, 2004 KL     Release v1.6
#          Oct  25, 2004 KL     Release v1.7
#          May   6, 2005 KL     Release v1.8
#          Mar  24, 2006 KL     Beta release v1.9beta
#          Jul  28, 2006 KL     Release v1.9
#          Jul  28, 2009 JH     Release v1.10
#          Jan  13, 2011 EH,KL  Beta release v1.11beta
#          Aug  19, 2011 EH     Beta release v1.11beta2
#          Dec  30, 2011 EH     Release v1.11
#          Mar  28, 2012 EH     Release v1.11.1
#          Dec  13, 2012 EH     Beta release v1.12beta
#          May  14, 2013 EH     Beta release v1.12beta2
#          Oct  11, 2013 EH     Release v1.12
#          Jul  21, 2014 KL     Commissioning release v1.13 GMOS Ham
#          Oct  10, 2014 MS     Patch for Commissioning release v1.13 GMOS Ham
#          Jan  30, 2015 KL     Release v1.13
#          Dec   7, 2015 KL	Release v1.13.1
#          Apr  21, 2017 KL     Commissioning release v1.14comm
#          Jul  20, 2017 KL     Release v1.14
#
# The latest version of the Gemini IRAF package is only compatible with
# IRAF v2.16 from AstroConda (or Ureka)
if (defpar ("release")) {
    if (release >= "2.16") {
        # There appear to be several issues related to the new image template
        # code in IRAF v2.16, so do not use the new image templates 
        if (defvar("use_new_imt")) {
            set use_new_imt = "no"
        }
        if (!access("iraf$ur/") && !access("iraf$conda_build.sh")) {
            printf ("WARNING: The Gemini IRAF package is not compatible \n")
            printf ("         with IRAF v2.16, unless installed using AstroConda \n")
            printf ("         or Ureka \n")
            printf ("Tested with IRAF 2.16 from AstroConda\n")
            sleep 10
        }
    } else {
        printf ("WARNING: The Gemini IRAF pacakge is only compatible with \n")
        printf ("         IRAF v2.16 and above.\n")
        printf ("Tested with IRAF 2.16 from AstroConda\n")
        sleep 10
    }
    #} else if ((release < "2.14.1") || (release > "2.15.1a")) {
    #    printf ("WARNING: The Gemini IRAF package is only compatible with\n")
    #    printf ("         versions of IRAF between v2.14.1 and v2.15.1a\n")
    #    printf ("         and with v2.16 in Ureka\n")
    #    printf ("Tested with IRAF 2.16 from Ureka/AstroConda\n")
    #    sleep 10
    #}
} else {
        printf ("WARNING: The Gemini IRAF package is only compatible with\n")
        printf ("         IRAF v2.16 and above, as distributed through AstroConda\n")
        printf ("         or Ureka\n")
        printf ("Tested with IRAF 2.16 from AstroConda\n")
        sleep 10
}
;

# Load necessary packages - only those that are used by all packages
images
plot
dataio
lists
tv
utilities
noao
proto
imred
ccdred
bias
ccdtest
utilities
astutil
digiphot
apphot
ptools
if(defpac("stsdas") == no) stsdas (motd-)
;
tools
fitting
artdata

# Load the FITSUTIL package.
if (deftask("fitsutil") == no) {
  print("ERROR: The FITSUTIL package is required but not defined.")
  bye()
}
;
fitsutil

nproto

reset imtype = "fits"

flpr

cl < "gemini$lib/zzsetenv.def"

set quirc       = "gemini$quirc/"
set nifs        = "gemini$nifs/"
set niri        = "gemini$niri/"
set gnirs       = "gemini$gnirs/"
set gmos        = "gemini$gmos/"
set midir       = "gemini$midir/"
set mostools    = "gemini$gmos/mostools/"
set oscir       = "gemini$oscir/"
set gemtools    = "gemini$gemtools/"
set gcal        = "gemini$gcal/"
set flamingos   = "gemini$flamingos/"
set f2          = "gemini$f2/"
set gsaoi       = "gemini$gsaoi/"


package gemini, bin = gembin$

task quirc.pkg    = quirc$quirc.cl
task nifs.pkg     = nifs$nifs.cl
task niri.pkg     = niri$niri.cl
task gnirs.pkg    = gnirs$gnirs.cl
task gmos.pkg     = gmos$gmos.cl
task oscir.pkg    = oscir$oscir.cl
task gemtools.pkg = gemtools$gemtools.cl
task flamingos.pkg  = flamingos$flamingos.cl
task f2.pkg       = f2$f2.cl
task midir.pkg    = midir$midir.cl
task gsaoi.pkg    = gsaoi$gsaoi.cl


task $sed = $foreign
hidetask sed

if (motd) {
print(" ")
print("     +------------------- Gemini IRAF Package -------------------+")
print("     |                 Version 1.14, July 20, 2017               |")
print("     |                                                           |")
print("     |               Requires IRAF v2.16 or greater              |")
print("     |              Tested with AstroConda IRAF v2.16            |") 
print("     |             Gemini Observatory, Hilo, Hawaii              |")
print("     |    Please use the help desk for submission of questions   |")
print("     |  http://www.gemini.edu/sciops/helpdesk/helpdeskIndex.html |")
print("     |     You can also check the Data Reduction User Forum      |")
print("     |                 http://drforum.gemini.edu                 |")
print("     +-----------------------------------------------------------+")
print(" ")
print("     Warning setting imtype=fits")
if (defvar("use_new_imt")) {
    print("     Warning setting use_new_imt=no")
}
print(" ")

}
;
clbye()
