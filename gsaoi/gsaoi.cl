# Copyright(c) 2012-2017 Association of Universities for Research in Astronomy, Inc.
#
# Package script for the gsaoi package
#
# Version: Dec  13, 2012 EH     Beta release v1.12beta
#          May  14, 2013 EH     Beta release v1.12beta2
#          Oct  11, 2013 EH     Release v1.12
#          Jan  30, 2015 KL     Release v1.13
#	   Dec   7, 2015 KL	Release v1.13.1
#          Jul  20, 2017 KL     Release v1.14
#
# load necessary packages - gemini loads most of the packages
gemtools

package gsaoi

# Generic preparations
task gacalfind=gsaoi$gacalfind.cl
task gacaltrim=gsaoi$gacaltrim.cl
task gadark=gsaoi$gadark.cl
task gadimschk=gsaoi$gadimschk.cl
task gadisplay=gsaoi$gadisplay.cl
task gafastsky=gsaoi$gafastsky.cl
task gaflat=gsaoi$gaflat.cl
task gaimchk=gsaoi$gaimchk.cl
task gamosaic=gsaoi$gamosaic.cl
task gaprepare=gsaoi$gaprepare.cl
task gareduce=gsaoi$gareduce.cl
task gasky=gsaoi$gasky.cl
task gastat=gsaoi$gastat.cl
hidetask gacaltrim
hidetask gadimschk
hidetask gaimchk
hidetask gastat

# Information
task gsaoiinfo=gsaoi$gsaoiinfo.cl

# Examples
task gsaoiexamples=gsaoi$gsaoiexamples.cl

clbye()
