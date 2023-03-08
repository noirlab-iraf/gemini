# Copyright(c) 2002-2015 Association of Universities for Research in Astronomy, Inc.
# Copyright(c) 2009      James E.H. Turner

include <gemini.h>
include <error.h>
include <imhdr.h>
include <tbset.h>
include <smw.h>
include <time.h>

procedure t_gfcube()

# Reconstruct x/y/lambda datacube from extracted IFU spectra
# Version     Sep  9, 2002  JT  (preliminary release in v1.4)
#             Mar 26, 2004  JT  (update for GMOS-S IFU mapping)
#               - new version requires INSTRUME keyword to be GMOS-N/S
#             Mar 06, 2009  JT  (compensate for atmospheric dispersion)
#             Apr 17, 2009  JT  (add data quality output)
#             May 07, 2009  JT  (add data quality input & simple log file)
#             Jun 19, 2009  JT  (make arrays bigger for optical distortion)
#             Sep 17, 2013  JT  (bitmask controls use of input data quality)
#             Mar 28, 2014  JT  (propagate variance & output flux/arcsec^2)
#             Apr 25, 2014  JT  (conserve variance overall, not in 1 pix)
# ----

char  inimage[SZ_FNAME], outimage[SZ_FNAME], outprefix[SZ_FNAME]#, field[2]
char  logfile[SZ_FNAME]
real  ssample
int   bitmask
bool  fl_atmdisp, fl_flux, fl_var, fl_dq
#bool  verbose
#int   status

# Local variables
char sciname[SZ_FNAME], mdfname[SZ_FNAME], outname[SZ_FNAME],
     iphuname[SZ_FNAME], ophuname[SZ_FNAME], imsuffix[SZ_FNAME],
     inst_name[SZ_KWVAL], varname[SZ_FNAME], dqname[SZ_FNAME],
     outvarname[SZ_FNAME], outdqname[SZ_FNAME], msg[SZ_LINE], 
     timestr[SZ_TIME], unit[SZ_KWVAL]
pointer sciext, mdfext, varext, dqext, outext, inphu, outphu
pointer outvarext, outdqext, outarr, outvararr, outdqarr, glog
int fnum[1500], nc, xdim, ydim, ldim, lstart, noutpix, imexists, stat
real xc[1500], yc[1500], xmin, xmax, ymin, ymax, xsample, ysample
real lcrpix, lcrval, lsample, flip
real xadisp[6500], yadisp[6500], lc[6500]
real elevation, PA, parA, height, latitude, temperature, pressure, humidity
long clkval

real  clgetr()
int   clgeti()
bool  clgetb()
pointer immap(), tbtopn(), open()
int imgeti(), imaccess(), strcmp(), gte_blankpar(), gin_ismef(),
    gin_multype(), gin_valid()
long clktime()

define  degtorad    0.017453292519943295d0
define  radtodeg    57.295779513082323d0
define  hrstorad    0.26179938779914941d0
define  fibpsqas    28.8675                # for GMOS-IFUs at f/16

begin
	
  # Read task parameters:
  call clgstr("inimage", inimage, SZ_FNAME)
  call clgstr("outimage", outimage, SZ_FNAME)
  call clgstr("outprefix", outprefix, SZ_FNAME)
  ssample = clgetr("ssample")
  bitmask = clgeti("bitmask")
  fl_atmdisp = clgetb("fl_atmdisp")
  fl_flux = clgetb("fl_flux")
  fl_var = clgetb("fl_var")
  fl_dq = clgetb("fl_dq")
  call clgstr("logfile", logfile, SZ_FNAME)

  # If the log file name is blank, get it from the GMOS parameters:
  if (gte_blankpar(logfile)==YES)
    call clgstr("gmos.logfile", logfile, SZ_FNAME)

  # Open the (currently basic) log file if specified, otherwise
  # disable logging:
  if (gte_blankpar(logfile)==NO) {
    iferr(glog = open(logfile, APPEND, TEXT_FILE))
      call error(1, "gfcube: cannot open log file for writing")
  }
  else
    glog = NULL

  # Get readable time string for log:
  clkval = clktime(0)
  call cnvtime(clkval, timestr, SZ_TIME)

  # Print opening log entry:
  call printlog(glog, "----------------------------------------------------------------------------\n")
  call sprintf(msg, SZ_LINE, "GFCUBE -- %s\n\n")
    call pargstr(timestr)
  call printlog(glog, msg)

  # Check that the spatial output sampling is non-zero:
  if (ssample==0.0)
    call error(1, "gfcube: ssample must be > 0")
  
  # Set sampling in x/y individually:
  xsample = ssample
  ysample = ssample
  
  # Make sure a unique input image exists:
  iferr(imexists = imaccess(inimage, READ_ONLY))
    call error(1, "gfcube: inimage is not unique - specify extension")
  if (imexists==NO) call error(1, "gfcube: inimage does not exist")

  # Make sure input name contains file extension:
  call gin_afext(inimage, GEM_DMULTEXT, inimage, YES)
  
  # Complain if input name contains image ext or sec:
  call gin_gsuf(inimage, imsuffix)
  if (strcmp(imsuffix, "")!=0)
	call error(1, "gfcube: inimage contains image extension or section")

  # Make sure input image is a MEF file:
  if (gin_ismef(inimage)==NO)
    call error(1, "gfcube: inimage is not a valid multi-extension FITS file")

  # Determine output name and ensure it is a complete FITS file:
  if (gte_blankpar(outimage)==YES) {
    if (gte_blankpar(outprefix)==YES)
      call error(1, "gfcube: both outimage and outprefix are blank")
    call gin_make(outprefix, inimage, GEM_DMULTEXT, outname)
  }
  else {
    call gin_gsuf(outimage, imsuffix)
    if (strcmp(imsuffix, "")!=0)
      call error(1, "gfcube: outimage contains image extension or section")
    call gin_afext(outimage, GEM_DMULTEXT, outname, NO)
    if(gin_multype(outname)==NO)
      call error(1, "gfcube: outimage type is not valid (try FITS)")
  }

  # Check for problematic output name:
  if (imaccess(outname, READ_ONLY)==YES)
    call error(1, "gfcube: outimage already exists")
  if (gin_valid(outname)==NO)
    call error(1, "gfcube: outimage name is not valid")

  # Log input and output filenames & parameters that aren't otherwise obvious:
  call sprintf(msg, SZ_LINE, "Input file:  %s\n")
    call pargstr(inimage)
  call printlog(glog, msg)
  call sprintf(msg, SZ_LINE, "Output file: %s\n\n")
    call pargstr(outname)
  call printlog(glog, msg)
  call sprintf(msg, SZ_LINE, "bitmask     = %d\n")
    call pargi(bitmask)
  call printlog(glog, msg)
  call sprintf(msg, SZ_LINE, "fl_atmdisp  = %b\n")
    call pargb(fl_atmdisp)
  call printlog(glog, msg)
  call sprintf(msg, SZ_LINE, "fl_flux     = %b\n\n")
    call pargb(fl_flux)
  call printlog(glog, msg)

  # Generate extension names from file names:
  call gin_psuf(inimage, "[SCI,1]", sciname)
  call gin_psuf(inimage, "[VAR,1]", varname)
  call gin_psuf(inimage, "[DQ,1]", dqname)
  call gin_psuf(inimage, "[MDF]", mdfname)
  call gin_psuf(inimage, "[0]", iphuname)
  call gin_psuf(outname, "[0,append]", ophuname)
  call gin_psuf(outname, "[SCI,1,append]", outname)
  call gin_psuf(outname, "[VAR,1,append]", outvarname)
  call gin_psuf(outname, "[DQ,1,append]", outdqname)

  # Open sci extension and, if present, VAR/DQ extensions (whether or not
  # propagating them to the output):
  sciext = immap(sciname, READ_ONLY, NULL)
  iferr (varext = immap(varname, READ_ONLY, NULL)) varext = NULL
  iferr (dqext = immap(dqname, READ_ONLY, NULL)) dqext = NULL

  # Check that it's 2D and matches the VAR/DQ:
  if (imgeti(sciext, "i_naxis") != 2)
    call error(1, "gfcube: input extension [SCI,1] is not 2 dimensional")

  if (varext != NULL)
    if (((IM_LEN(varext, 1) != IM_LEN(sciext, 1)) ||
	(IM_LEN(varext, 2) != IM_LEN(sciext, 2))))
      call error(1, "gfcube: input extension [VAR,1] does not match [SCI,1]")
  if (dqext != NULL)
    if (((IM_LEN(dqext, 1) != IM_LEN(sciext, 1)) ||
        (IM_LEN(dqext, 2) != IM_LEN(sciext, 2))))
      call error(1, "gfcube: input extension [DQ,1] does not match [SCI,1]")

  # Open PHU and get instrument name:
  inphu = immap(iphuname, READ_ONLY, NULL)
  iferr(call imgstr(inphu, "INSTRUME", inst_name, SZ_KWVAL))
    inst_name[1] = EOS
  else call gte_unpad(inst_name, inst_name, SZ_KWVAL)
  if (strcmp(inst_name, "GMOS") == 0)
    call strcpy ("GMOS-N", inst_name, SZ_KWVAL)
  if (strcmp(inst_name, "GMOS-N") != 0 && strcmp(inst_name, "GMOS-S") != 0)
    call error(1, "gfcube: header keyword INSTRUME must be GMOS-N or GMOS-S")

  # Open MDF:
  mdfext = tbtopn(mdfname, READ_ONLY, NULL)

  # Read MDF co-ordinates & boundaries:
  call getmap(mdfext, inst_name, fnum, xc, yc, xmin, xmax, ymin, ymax, nc)
	
  # Close MDF ext:
  call tbtclo(mdfext)

  # Check that nrows == number of good fibres in MDF:
  if (IM_LEN(sciext, 2) != nc)
    call error(1, "gfcube: num image rows != num good fibres in MDF")

  # If correcting atmospheric dispersion, get environmental and pointing
  # data from the headers & observatory database:
  if (fl_atmdisp)
    call getpointenv(inphu, elevation, PA, parA, flip, height, latitude,
                     temperature, pressure, humidity, glog)

  # Copy primary header to output:
  outphu = immap(ophuname, NEW_COPY, inphu)
  call imunmap(inphu)
  call imunmap(outphu)

  # Calculate limits of output cube:
  call calcxydim(xmin, xmax, ymin, ymax, xsample, ysample, xdim, ydim)

  # Determine where to truncate spectra:
  call getllim(sciext, lstart, ldim)

  # Get spectral WCS from the input header:
  call getspecwcs(sciext, lstart, ldim, lcrpix, lcrval, lsample, lc)
  
  # Get flux units from header & update them if scaling to flux/arcsec^2:
  call getunit(sciext, unit, fl_flux)

  # Create temporary array(s) for output cube(s) (can't open IMIO output
  # buffers for SCI, VAR and DQ in different extensions at the same time):
  noutpix = xdim*ydim*ldim
  call malloc(outarr, noutpix, TY_REAL)
  if (fl_var && varext != NULL)  # Only use VAR if requested & present in input
    call malloc(outvararr, noutpix, TY_REAL)
  else
    outvararr = NULL
  if (fl_dq)
    call malloc(outdqarr, noutpix, TY_INT)
  else
    outdqarr = NULL

  # If correcting atmospheric dispersion, get the spectral WCS values
  # from the header and derive the atmospheric dispersion vector for
  # each wavelength (relative to the starting wavelength):
  if (fl_atmdisp)
    call calcdispoff(lc, ldim, elevation, PA, parA, flip, height, latitude,
      temperature, pressure, humidity, xadisp, yadisp, glog)

  # Resample input spectra into datacube:
  call printlog(glog, "Interpolating...\n")
  call recon(sciext, varext, dqext, outarr, outvararr, outdqarr, nc,
    xsample, ysample, xmin, ymin, xdim, ydim, ldim, lstart, fnum, xc, yc,
    xadisp, yadisp, fl_atmdisp, bitmask, glog)

  # Scale output to flux per arcsec^2 if requested:
  if (fl_flux) call fscale(outarr, outvararr, noutpix)

  # Write the data extension(s) to file, first setting their dimensions &
  # WCS in the header (NB. opening more than one extension at a time for
  # writing seems to cause corruption, maybe because appending doesn't work
  # reliably until the previous extension is finalized):
  call writecube(outarr, sciext, outname, xdim, ydim, ldim, xmin, xmax,
    ymin, ymax, xsample, ysample, lsample, lcrpix, lcrval, unit, TY_REAL)
  if (outvararr != NULL)
    call writecube(outvararr, NULL, outvarname, xdim, ydim, ldim, xmin, xmax,
      ymin, ymax, xsample, ysample, lsample, lcrpix, lcrval, "", TY_REAL)
  if (fl_dq)
    call writecube(outdqarr, NULL, outdqname, xdim, ydim, ldim, xmin, xmax,
      ymin, ymax, xsample, ysample, lsample, lcrpix, lcrval, "", TY_USHORT)

  # Free output arrays:
  call mfree(outarr, TY_REAL)
  if (outvararr != NULL) call mfree(outvararr, TY_REAL) 
  if (outdqarr != NULL) call mfree(outdqarr, TY_INT)

  # Close input file:
  call imunmap(sciext)
  if (varext != NULL) call imunmap(varext)
  if (dqext != NULL) call imunmap(dqext)

  # Print closing log entry:
  call printlog(glog, "----------------------------------------------------------------------------\n")

  # Close the log file:
  if (glog != NULL) call close(glog)
	
end


procedure getmap(mdfext, inst_name, fnum, xc, yc, xmin, xmax, ymin, ymax, nc)

# Get co-ordinates & field num for good fibres and return co-ordinate
# boundaries for the big field

# Input instrument name should be GMOS-N or GMOS-S

pointer mdfext
char inst_name[SZ_KWVAL]
int fnum[ARB], nc
real xc[ARB], yc[ARB], xmin, xmax, ymin, ymax

pointer xcol, ycol, becol
int n, nrows, beam, inst
bool initlim

define  GMOS_N        0
define  GMOS_S        1

int tbpsta(), strcmp()

begin

  # Get num rows and col nums:
  nrows = tbpsta(mdfext, TBL_NROWS)
  call tbcfnd(mdfext, "xinst", xcol, 1)
  call tbcfnd(mdfext, "yinst", ycol, 1)
  call tbcfnd(mdfext, "beam", becol, 1)

  # Check instrument name:
  if (strcmp(inst_name, "GMOS-S")==0) inst = GMOS_S
  else inst = GMOS_N
	
  # Loop through table rows & get co-ords for good fibres:
  nc=0
  initlim = true
  for (n=1; n <= nrows; n=n+1) {
    call tbegti(mdfext, becol, n, beam)
    if (beam!=-1) {
      nc = nc+1
      call tbegtr(mdfext, xcol, n, xc[nc])
      call tbegtr(mdfext, ycol, n, yc[nc])
      # Determine field from xc:
      if (inst == GMOS_N && xc[nc] < 30.0 ||
          inst == GMOS_S && xc[nc] > 30.0) {
        fnum[nc] = 1
        # Derive field limits:
        if (initlim) {
          xmin = xc[nc]; xmax = xmin
          ymin = yc[nc]; ymax = ymin
          initlim = false
        } else {
          if (xc[nc] < xmin) xmin = xc[nc]
          else if (xc[nc] > xmax) xmax = xc[nc]
          if (yc[nc] < ymin) ymin = yc[nc]
          else if (yc[nc] > ymax) ymax = yc[nc]
        }
      } # if (field==1)
      else fnum[nc] = 2
    } # end if (beam==1)
  } # end for n

end


procedure setdim(outext, xmin, xmax, ymin, ymax, xsample, ysample, lsample,
  xdim, ydim, ldim, lcrpix, lcrval, unit, dtype)

# Set output image dimensions & WCS

pointer outext
real xmin, xmax, ymin, ymax, xsample, ysample, lsample, lcrpix, lcrval
int xdim, ydim, ldim, dtype
char unit[SZ_KWVAL]

int n
char kwbuf[SZ_KWVAL]

int imaccf(), strcmp()

begin

  # Set output data type:
  call imputi(outext, "i_pixtype", dtype)

  # Set dimensionality:
  call imputi(outext, "i_naxis", 3)
  call imastr(outext, "CTYPE1", "LINEAR")
  call imastr(outext, "CTYPE2", "LINEAR")
  call imastr(outext, "CTYPE3", "LAMBDA")

  call imputi(outext, "i_naxis1", xdim)
  call imputi(outext, "i_naxis2", ydim)

  # Set spatial increment:
  call imaddr(outext, "CD1_1", xsample)
  call imaddr(outext, "CD2_2", ysample)
  call imaddr(outext, "CDELT1", xsample)
  call imaddr(outext, "CDELT2", ysample)
  
  # Set starting co-ord to centre of first pix:
  call imaddr(outext, "CRPIX1", 1.0)
  call imaddr(outext, "CRVAL1", xmin)
  call imaddr(outext, "CRPIX2", 1.0)
  call imaddr(outext, "CRVAL2", ymin)

  # Set wavelength dim:
  call imaddi(outext, "i_naxis3", ldim)
  call imaddi(outext, "DISPAXIS", 3)
	
  # Set wavelength ref & increment:
  if (lcrval != 0.0) {
    call imaddr(outext, "CRPIX3", lcrpix)
    call imaddr(outext, "CRVAL3", lcrval)
    call imaddr(outext, "CD3_3", lsample)
    call imaddr(outext, "CDELT3", lsample)
  }

  # Make sure LTM values are 1:
  call imaddr(outext, "LTM1_1", 1.0)
  call imaddr(outext, "LTM2_2", 1.0)
  call imaddr(outext, "LTM3_3", 1.0)
	
  # Delete any extra WCS keywords:
  iferr(call imdelf(outext, "WCSDIM"));
  iferr(call imdelf(outext, "LTM1_2"));
  iferr(call imdelf(outext, "LTM2_1"));
  iferr(call imdelf(outext, "CD1_2"));
  iferr(call imdelf(outext, "CD2_1"));
  iferr(call imdelf(outext, "WAT0_001"));
  iferr(call imdelf(outext, "WAT1_001"));
  iferr(call imdelf(outext, "WAT2_001"));
  iferr(call imdelf(outext, "REFPIX1"));
  iferr(call imdelf(outext, "REFPIX2"));

  # Update output units if known from the input:
  if (strcmp(unit, "") != 0)
    call imastr(outext, "BUNIT", unit)

  # Delete apextract keywords:
  # (need a more efficient way of doing this in final code)
  for (n=1; n < 1000; n=n+1) {
    call sprintf(kwbuf, SZ_KWVAL, "APNUM%d")
     call pargi(n)
    if (imaccf(outext, kwbuf) == YES) call imdelf(outext, kwbuf)
	else break
  }
  for (n=1; true; n=n+1) {
    call sprintf(kwbuf, SZ_KWVAL, "APID%d")
     call pargi(n)
    if (imaccf(outext, kwbuf) == YES) call imdelf(outext, kwbuf)
	else break
  }
  for (n=1000; true; n=n+1) {
    call sprintf(kwbuf, SZ_KWVAL, "AP%d")
     call pargi(n)
    if (imaccf(outext, kwbuf) == YES) call imdelf(outext, kwbuf)
	else break
  }
	
end


procedure calcxydim(xmin, xmax, ymin, ymax, xsample, ysample, xdim, ydim)

# Calculate output image dimensions

pointer sciext
real xmin, xmax, ymin, ymax, xsample, ysample
int xdim, ydim

real cenx, ceny

begin

  # Calc / set spatial dimensions from sampling:
  xdim = nint((xmax-xmin) / xsample)
  ydim = nint((ymax-ymin) / ysample)
	
  # Calc spatial centre co-ords:
  cenx = 0.5 * (xmin + xmax)
  ceny = 0.5 * (ymin + ymax)
	
  # Round co-ordinate boundaries to match output grid:
  # (and reverse x co-ordinate, to get correct orientation)
#  xmin = cenx - 0.5 * (xdim-1) * xsample
#  xmax = cenx + cenx - xmin
  xmax = cenx - 0.5 * (xdim-1) * xsample  # really xmin
  xmin = cenx + cenx - xmax               # really xmax
  ymin = ceny - 0.5 * (ydim-1) * ysample
  ymax = ceny + ceny - ymin

  # Reverse the sample increment in x & y:
  # (was x only until rev. 1.7, but GMOS images have North down at PA=0)
  xsample = -xsample
  ysample = -ysample
  
end


procedure getpointenv(phu, elevation, PA, parA, flip, height, latitude,
                      temperature, pressure, humidity, glog)

# Gather the observatory, pointing and environmental data needed to
# calculate atmospheric dispersion.

pointer phu, glog
real elevation, PA, parA, height, latitude, temperature, pressure, humidity
real flip

# The units are as follows, similar to SLALIB but with angles kept in
# degrees at this stage, mainly for ergonomic reasons:
#
#   elevation  degrees
#          PA  degrees
#        parA  degrees
#      height  metres above sea level
#    latitude  degrees
# temperature  Kelvin
#    pressure  millibars
#    humidity  fraction (0-1)

char observatory[SZ_LINE]
pointer obsdb
double sla_HA, sla_dec, sla_latitude
real dispA
char msg[SZ_LINE]
int port

real imgetr()
double imgetd()
int imgeti(), imaccf(), strcmp()
pointer obsopen()
real obsgetr()
double slPA()

begin

  # Get the observatory from the header if the keyword exists, otherwise
  # retrieve it from the gmos package parameters:
  if (imaccf(phu, "OBSERVAT") == YES) {
    call imgstr(phu, "OBSERVAT", observatory, SZ_LINE)
  }
  else {
    call clgstr("observatory", observatory, SZ_LINE)
  }
  call strlwr(observatory)

  # Read the latitude & height from the observatory database with xtools:
  obsdb = obsopen(observatory)
  latitude = obsgetr(obsdb, "latitude")
  height = obsgetr(obsdb, "altitude")
  call obsclose(obsdb)

  # Read header pointing keyword values, generating an error if missing:
  if (imaccf(phu, "ELEVATIO") == YES) {
    elevation = imgetr(phu, "ELEVATIO")
  }
  else {
    call error(1, "gfcube: fl_atmdisp+ and no ELEVATIO header keyword")
  }
  if (imaccf(phu, "DEC") == YES) {
    sla_dec = degtorad*imgetd(phu, "DEC")
  }
  else {
    call error(1, "gfcube: fl_atmdisp+ and no DEC header keyword")
  }
  if (imaccf(phu, "PA") == YES) {
    PA = imgetr(phu, "PA")
  }
  else {
    call error(1, "gfcube: fl_atmdisp+ and no PA header keyword")
  }
  if (imaccf(phu, "HA") == YES) {
    # Convert hour angle and dec to parallactic angle (-180>180deg):
    sla_HA = hrstorad*imgetd(phu, "HA")
    sla_latitude = degtorad*double(latitude)
    parA = radtodeg*real(slPA(sla_HA, sla_dec, sla_latitude))
  }
  else {
    call error(1, "gfcube: fl_atmdisp+ and no HA header keyword")
  }

  # Log the angle of dispersion WRT the field:
  dispA = parA-PA
  if (dispA < -180.) dispA = dispA + 360.
  if (dispA > 180.) dispA = dispA - 360.
  call sprintf(msg, SZ_LINE, "Parallactic angle minus PA = %.1f deg\n\n")
    call pargr(dispA)
  call printlog(glog, msg)

  # Determine science fold mirror flip from instrument port:
  if (imaccf(phu, "INPORT") == YES)
    port = imgeti(phu, "INPORT")
  else {
    port = 0
    call printlog(glog,
      "WARNING: no port number (INPORT) in header; assuming side port\n")
  }
  if (port == 1)
    flip = 1.0
  else
    flip = -1.0

  # Read environmental header keyword values, generating a warning and
  # default value if missing:
  #
  # From a quick inspection of GEA plots, the following values seem
  # fairly close to average: MK 3deg C, RH 30%, P 615. They could be
  # improved, but the calculation is not highly sensitive to these
  # environmental parameters anyway.
  #
  if (imaccf(phu, "TAMBIENT") == YES) {
    temperature = 273.15 + imgetr(phu, "TAMBIENT")
    call sprintf(msg, SZ_LINE, "T = %.0f K\n")
      call pargr(temperature)
    call printlog(glog, msg)
  }
  else {
    if (strcmp(observatory, "gemini-north") == 0) temperature = 276.
    else temperature = 283.
    call sprintf(msg, SZ_LINE,
      "WARNING: no T (TAMBIENT) in header: using %.0f K\n")
      call pargr(temperature)
    call printlog(glog, msg)
  }
  if (imaccf(phu, "PRESSUR2") == YES) {
    pressure = 0.01 * imgetr(phu, "PRESSUR2")
    call sprintf(msg, SZ_LINE, "P = %.0f mb\n")
      call pargr(pressure)
    call printlog(glog, msg)
  }
  else {
    if (strcmp(observatory, "gemini-north") == 0) pressure = 615.
    else pressure = 730.
    call sprintf(msg, SZ_LINE,
      "WARNING: no P (PRESSUR2) in header: using %.0f mb\n")
      call pargr(pressure)
    call printlog(glog, msg)
  }
  if (imaccf(phu, "HUMIDITY") == YES) {
    humidity = 0.01 * imgetr(phu, "HUMIDITY")
    call sprintf(msg, SZ_LINE, "RH = %.0f %%%%\n")
      call pargr(100.*humidity)
    call printlog(glog, msg)
  }
  else {
    humidity = 0.3
    call printlog(glog, "WARNING: no RH (HUMIDITY) in header: using 30 %%\n")

    # Note that in practice, the humidity seems to make little or no
    # difference to the atmospheric dispersion results (more of an
    # achromatic effect?).
  }

  call printlog(glog, "")

end


procedure getspecwcs(img, lstart, ldim, lcrpix, lcrval, lsample, lc)

# Get spectral co-ordinate value for each wavelength plane of output data
# and simple spectral WCS from header
#
# Come back later and see if I can figure out how to apply a shift to the
# SMW WCS itself and write that out directly in setdim

pointer img
int lstart, ldim
real lcrpix, lcrval, lsample, lc[ARB]

pointer smw_wcs, smw_ct
real x1, x2
int l
bool gotkw

pointer smw_openim(), smw_sctran()
real smw_c1tranr(), imgetr()

begin

  # Open NOAO SMW WCS structure based on image header:
  smw_wcs = smw_openim(img)

  # Set up transformation:
  smw_ct = smw_sctran(smw_wcs, "logical", "world", 1)
  
  # Populate array of WCS values corresponding to output pixel values:
  # (values are in Angstroms; does SMW guarantee that (has no docs)?)
  for (l=0; l < ldim; l=l+1)
    lc[l+1] = smw_c1tranr(smw_ct, real(lstart)+l)

  # Close SMW structure:  
  call smw_close(smw_wcs)

  # Get wavelength increment:
  ifnoerr(lsample = imgetr(img, "CD1_1"));
  else ifnoerr(lsample = imgetr(img, "CDELT1"));
  else lsample = 1.0

  # Get wavelength ref pix:
  gotkw = false
  ifnoerr(lcrval = imgetr(img, "CRVAL1"))
    ifnoerr(lcrpix = imgetr(img, "CRPIX1")) {
      lcrpix = lcrpix-real(lstart)+1.0
      gotkw = true
    }
  if (!gotkw) {
    lcrpix = 1.0
    lcrval = 0.0
  }

end


procedure getunit(sciext, unit, flux)

pointer sciext
char unit[SZ_KWVAL]
bool flux

begin

  # Get the units from the SCI header, or a blank string if undefined:
  iferr(call imgstr(sciext, "BUNIT", unit, SZ_KWVAL))
    unit[1] = EOS
  else if (flux)
    call strcat("/arcsec2", unit, SZ_KWVAL)

end


procedure calcdispoff(lc, ldim, elevation, PA, parA, flip, height, latitude,
  temperature, pressure, humidity, xadisp, yadisp, glog)

# Derive offsets with wavelength due to atmospheric dispersion, in
# arcseconds relative to the mean, using pointing, environmental data etc.
# (via SLALIB). The output vectors are the displacement of the target with
# respect to the IFU field, assuming x increases to the left(?) and y to
# the top. [CHECK]

real lc[ARB], xadisp[ARB], yadisp[ARB]
real elevation, PA, parA, height, latitude, temperature, pressure, humidity
real flip
int ldim
pointer glog

double ZD, dZD, sla_height, sla_temperature, sla_pressure, sla_humidity
double sla_wavelength, sla_latitude, dZD_ref
double as_sintheta, as_costheta
real   xavg, yavg
int l
char msg[SZ_LINE]

define  Tlapserate  0.0065d0
define  epsilon     1d-8

begin

  # Convert angles to radians and all variables to double for SLALIB:
  ZD = degtorad*(90.0d0-double(elevation))
  sla_height = double(height)
  sla_temperature = double(temperature)
  sla_pressure = double(pressure)
  sla_humidity = double(humidity)
  sla_latitude = degtorad*double(latitude)

  # Precalculate the conversion from refraction in elevation in radians to
  # displacement along each IFU axis in arcseconds:
  as_sintheta = 3600.0d0*radtodeg*sin(degtorad*double(parA-PA))
  as_costheta = 3600.0d0*radtodeg*cos(degtorad*double(parA-PA))

  # Initialize sum of offsets for calculating the mean WRT first wavelength:
  xavg=0.0; yavg=0.0

  # Loop over output wavelength planes, calculating an offset in arcseconds
  # relative to the first wavelength, in the co-ordinate system of the MDF:
  for (l=1; l <= ldim; l=l+1) {

    # Convert Angstroms to Microns:
    sla_wavelength = 0.0001d0*double(lc[l])

    # Call SLALIB to calculate refraction at the relevant wavelength:
    # (could also use slRFCO(), slATMD())
    call slRFRO(ZD, sla_height, sla_temperature, sla_pressure,
      sla_humidity, sla_wavelength, sla_latitude, Tlapserate, epsilon, dZD)

    # Get the starting refraction on first iteration:
    if (l==1) dZD_ref = dZD

    # Convert the N-S dispersion WRT the first plane into an offset vector
    # in arcseconds:
    xadisp[l] = real(as_sintheta*(dZd-dZD_ref))
    yadisp[l] = real(as_costheta*(dZD-dZD_ref))*flip

    # Accumulate the sum of offsets to derive an average below:
    xavg = xavg + xadisp[l]
    yavg = yavg + yadisp[l]

  } # end (loop over wavelength planes)

  # Convert sum of offsets to mean:
  xavg = xavg / real(ldim)
  yavg = yavg / real(ldim)

  # Subtract the mean from all the offsets to get corrections WRT the
  # average image position (for a flat spectrum, this would match the image
  # formed by collapsing the cube in wavelength):
  for (l=1; l <= ldim; l=l+1) {
    xadisp[l] = xadisp[l] - xavg
    yadisp[l] = yadisp[l] - yavg
  }

  # Log the range of atmospheric dispersion:
  call printlog (glog,
    "\nCalculated atmospheric dispersion relative to the mean:\n")
  call sprintf(msg, SZ_LINE, "    %.1fA : x=%6.2f, y=%6.2f arcsec\n")
    call pargr(lc[1])
    call pargr(xadisp[1])
    call pargr(yadisp[1])
  call printlog(glog, msg)
  call sprintf(msg, SZ_LINE, "    %.1fA : x=%6.2f, y=%6.2f arcsec\n\n")
    call pargr(lc[ldim])
    call pargr(xadisp[ldim])
    call pargr(yadisp[ldim])
  call printlog(glog, msg)

end


procedure recon(sciext, varext, dqext, outarr, outvararr, outdqarr, nc,
  xsample, ysample, x1, y1, xdim, ydim, ldim, lstart, fnum, xc, yc,
  xadisp, yadisp, atmdisp, bitmask, glog)

# Resample the input data onto the datacube

pointer sciext, varext, dqext, outarr, outvararr, outdqarr
int nc, xdim, ydim, ldim, lstart, fnum[ARB], bitmask
real xsample, ysample, x1, y1, xc[ARB], yc[ARB], xadisp[ARB], yadisp[ARB]
bool atmdisp
pointer glog

# Local variables
long nspix, ninpix
pointer iwk, wk, xout, yout, vout, varout, dqout, xarr, yarr, varr, varvarr
pointer indat, invar, indq, inbuf, invarbuf, indqbuf
int x, y, l, n, ncf, istat, stat, ltmp
real xadj, yadj, sfact
char msg[SZ_LINE]

pointer imgs2r(), imgs2i()

begin

  nspix = xdim * ydim
  ninpix = ldim * nc

  # Factor by which we're upsampling the fibre pattern (for var. calc.):
  sfact = 1.0 / (xsample * ysample * fibpsqas)

  # Get memory for workspace arrays:
  call malloc(iwk, 31*nc+nspix, TY_INT)
  call malloc(wk, 8*nc, TY_REAL)
  call malloc(xout, nspix, TY_REAL)
  call malloc(yout, nspix, TY_REAL)
  call malloc(vout, nspix, TY_REAL)
  call malloc(varout, nspix, TY_REAL)
  call malloc(dqout, nspix, TY_INT)
  call malloc(xarr, nc, TY_REAL)
  call malloc(yarr, nc, TY_REAL)
  call malloc(varr, nc, TY_REAL)
  call malloc(varvarr, nc, TY_REAL)
  call malloc(inbuf, ninpix, TY_REAL)
  if (outvararr != NULL) call malloc(invarbuf, ninpix, TY_REAL)
  if (dqext != NULL) call malloc(indqbuf, ninpix, TY_INT)

  # Generate output grid co-ords:
  for (y=0; y < ydim; y=y+1)
    for (x=0; x < xdim; x=x+1) {
      Memr[xout+y*xdim+x] = x1 + x*xsample
      Memr[yout+y*xdim+x] = y1 - y*ysample  # Reversed ysample May 2009
    }

  # Get input array(s) and copy to safe buffer(s):
  ltmp = lstart+ldim-1
  indat = imgs2r(sciext, lstart, ltmp, 1, nc)
  call amovr(Memr[indat], Memr[inbuf], ldim*nc)
  if (outvararr != NULL) {
    invar = imgs2r(varext, lstart, ltmp, 1, nc)
    call amovr(Memr[invar], Memr[invarbuf], ldim*nc)
  }
  if (dqext != NULL) {
    indq = imgs2i(dqext, lstart, ltmp, 1, nc)
    call amovi(Memi[indq], Memi[indqbuf], ldim*nc)
  }
	
  # Work through datacube, one plane at a time:
  for (l=0; l < ldim; l=l+1) {
    
    # Clear output plane buffer(s):
    call aclrr(Memr[vout], nspix)
    call aclrr(Memr[varout], nspix)
    call aclri(Memi[dqout], nspix)

    # Use supplied adjustment for atmospheric dispersion or zero if
    # correction not requested:
    if (atmdisp) {
      xadj = xadisp[l+1]
      yadj = yadisp[l+1]
    } else {
      xadj = 0.0
      yadj = 0.0
    }

    # Extract co-ords & values for field 1 into arrays for interp:
    ncf = 0
    for(n=1; n <= nc; n=n+1) {
      if (fnum[n]==1) {
        # Copy pixel values to array (where DQ is good or unknown):
        ltmp = (n-1)*ldim+l
        Memr[varr+ncf] = Memr[inbuf+ltmp]
        Memr[xarr+ncf] = xc[n] - xadj
        Memr[yarr+ncf] = yc[n] - yadj
        if (outvararr != NULL) Memr[varvarr+ncf] = Memr[invarbuf+ltmp]
        # Keep this point if there's no DQ or the DQ is good (do the same
        # for VAR so it matches the SCI values):
        if (dqext == NULL) ncf = ncf+1
        else if (and(Memi[indqbuf+ltmp], bitmask) == 0) ncf = ncf+1
      }
    } # end for (n <= nc)

    # call PDA_IDBVIP(md,ncp,ndp,xd,yd,zd,nip,xi,yi,zi,iwk,wk,istat,status
    istat = 0; stat = 0
    call JDBVIP(1, 4, ncf, Memr[xarr], Memr[yarr], Memr[varr], nspix,
                Memr[xout], Memr[yout], Memr[vout], Memi[dqout],
                Memi[iwk], Memr[wk], istat, stat)

    if (istat != 0) {
      # call eprintf("\nWARNING: gfcube: interpolation failed at plane %d\n")
      call sprintf(msg, SZ_LINE,
        "WARNING: gfcube: interpolation failed at plane %d\n")
       call pargi(l+1)
      call printlog(glog, msg)

      # Set DQ bad if interpolation failed, since still uninitialized:
      call amovki(1, Memi[dqout], nspix)
    }

    # Copy output plane buffer to output array:
    call amovr(Memr[vout], Memr[outarr+l*nspix], nspix)
    if (outdqarr != NULL)
      call amovi(Memi[dqout], Memi[outdqarr+l*nspix], nspix)

    # Also resample the variance if applicable & copy it to the output:
    #
    # Since we don't propagate covariance, this variance is upweighted here by
    # the same factor as the fibre pattern is upsampled onto the cube, so as to
    # conserve the right errors overall when combining pixel values over
    # spatial regions during analysis. The error on any individual pixel value
    # will therefore be overestimated by 1.86x with the default ssample=0.1.
    if (outvararr != NULL) {
      if (istat == 0)  # Only interpolate if it succeeded for SCI
        # Setting first parameter to 3 re-uses previous interpolation weights:
        call JDBVIP(3, 4, ncf, Memr[xarr], Memr[yarr], Memr[varvarr], nspix,
                    Memr[xout], Memr[yout], Memr[varout], Memi[dqout],
                    Memi[iwk], Memr[wk], istat, stat)
        call amulkr(Memr[varout], sfact, Memr[outvararr+l*nspix], nspix)
    }

  } # end (for l < ldim)

  # Free buffers:	
  call mfree(iwk, TY_INT)
  call mfree(wk, TY_REAL)
  call mfree(xout, TY_REAL)
  call mfree(yout, TY_REAL)
  call mfree(vout, TY_REAL)
  call mfree(varout, TY_REAL)
  call mfree(xarr, TY_REAL)
  call mfree(yarr, TY_REAL)
  call mfree(varr, TY_REAL)
  call mfree(varvarr, TY_REAL)
  call mfree(inbuf, TY_REAL)
  if (outvararr != NULL) call mfree(invarbuf, TY_REAL)
  if (dqext != NULL) call mfree(indqbuf, TY_INT)

end


procedure fscale(sciarr, vararr, npix)

pointer sciarr, vararr
int npix

begin

  # Ideally one would derive the scaling from the MDF co-ordinate density
  # rather than the hard-wired scale factor fibpsqas but the MDF parsing in
  # getmap is already GMOS-specific so just keep this simple until such
  # time as more flexibility is needed.
  call amulkr(Memr[sciarr], fibpsqas, Memr[sciarr], npix)
  if (vararr != NULL)
    call amulkr(Memr[vararr], fibpsqas*fibpsqas, Memr[vararr], npix)

end


procedure getllim(sciext, lstart, ldim)

# Calculate useful wavelength limits for output cube
# This method turns out to be a bit hit-and-miss, pending DQ propagation

pointer sciext
int lstart, ldim

# Local variables
pointer indat
int nrows, l, n, ngood, minspec

pointer imgs2r()
int gnu_equalr()

begin

  # Get input dimension in lambda:
  ldim = IM_LEN(sciext, 1)
  nrows = IM_LEN(sciext, 2)

  # Set threshold for inclusion:
  minspec = nint(0.99 * nrows)
	
  # Get input image data:
  indat = imgs2r(sciext, 1, ldim, 1, nrows)

  # Find first wavelength with >99% of spectra:
  for (l=0; l < ldim; l=l+1) {

    # Count non-blank spectra in both fields:
    ngood = 0
    for(n=0; n < nrows; n=n+1)
      if (gnu_equalr(Memr[indat+n*ldim+l], 0.0)==NO)
        ngood = ngood+1
		
    if (ngood > minspec) break
				
  } # end for (l < ldim)

  if (l==ldim) lstart = 1
  else lstart = l+1


  # Find last wavelength element with >99% of spectra:
  for (l=ldim-1; l > -1; l=l-1) {

    # Count non-blank spectra in both fields:
    ngood = 0
    for(n=0; n < nrows; n=n+1)
      if (gnu_equalr(Memr[indat+n*ldim+l], 0.0)==NO)
        ngood = ngood+1
		
    if (ngood > minspec) break
				
  } # end for (l < ldim)

  if (l==-1) ldim = ldim-lstart+1
  else ldim = l-lstart+2

end


procedure writecube(arr, inheadext, outname, xdim, ydim, ldim, xmin, xmax,
  ymin, ymax, xsample, ysample, lsample, lcrpix, lcrval, unit, dtype)

pointer arr, inheadext
char outname[SZ_FNAME], unit[SZ_KWVAL]
int xdim, ydim, ldim, dtype, buftype
real xmin, xmax, ymin, ymax, xsample, ysample, lsample, lcrpix, lcrval

pointer outext, outdat

pointer immap(), imps3i(), imps3r(), imps3d()

begin

  # Expect a real, double or integer input array, depending on the
  # output image type:
  if (dtype == TY_INT || dtype == TY_SHORT || dtype == TY_USHORT)
    buftype = TY_INT
  else if (dtype == TY_REAL)
    buftype = TY_REAL
  else if (dtype == TY_DOUBLE)
    buftype = TY_DOUBLE
  else
    call error(1, "writecube(): unsupported dtype")

  # Open output extension:
  if (inheadext != NULL) 
    outext = immap(outname, NEW_COPY, inheadext)
  else
    outext = immap(outname, NEW_IMAGE, NULL)

  # Set dimensions in the header:
  call setdim(outext, xmin, xmax, ymin, ymax, xsample, ysample, lsample,
    xdim, ydim, ldim, lcrpix, lcrval, unit, dtype)

  # Get output buffer and copy datacube to output:
  if (buftype == TY_INT) {
    outdat = imps3i(outext, 1, xdim, 1, ydim, 1, ldim)
    call amovi(Memi[arr], Memi[outdat], xdim*ydim*ldim)
  }
  else if (buftype == TY_REAL) {
    outdat = imps3r(outext, 1, xdim, 1, ydim, 1, ldim)
    call amovr(Memr[arr], Memr[outdat], xdim*ydim*ldim)
  }
  else if (buftype == TY_DOUBLE) {
    outdat = imps3d(outext, 1, xdim, 1, ydim, 1, ldim)
    call amovd(Memd[arr], Memd[outdat], xdim*ydim*ldim)
  }

  # Close the extension:
  call imunmap(outext)

end


procedure printlog(glog, str)

pointer glog
char str[ARB]

begin

  # Write message to stdout:
  call printf(str)
  call flush(STDOUT)

  # Write to log if open:
  if (glog != NULL) {
    call fprintf(glog, str)
    call flush(glog)
  }

end

