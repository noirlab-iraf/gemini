# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
#
# Version    Sep 2014  JT  Initial version for Hamamatsu QE correction
#            Nov 2014  JT  Deal with binned data

include <imhdr.h>
include	<pkg/gtools.h>  # for identify
include <gemini.h>
include <ctype.h>
include <math/gsurfit.h>
include "../lib/irafdb/apertures.h"
include "../lib/irafdb/identify.h"

procedure t_gfunexwl()

# Determine wavelength for each pixel in a 2D mosaicked GMOS IFU spectrum,
# based on already-determined fibre wavelength solutions & traces in the
# IRAF database. The output is later used, for example, to apply a QE
# correction as a function of wavelength for each chip before doing things
# like continuum fitting & normalization, to help avoid discontinuities.

# Parameters:
char waveref[SZ_FNAME], outimage[SZ_FNAME]
char database[SZ_FNAME]
int  ncols, xord, yord

# Local variables:
char wavesciname[SZ_FNAME], wavephuname[SZ_FNAME] 
char tracedbname[SZ_FNAME], wavedbname[SZ_FNAME]
char outwavname[SZ_FNAME], dbsuff[SZ_FNAME], msg[SZ_LINE], ccdsum[SZ_KWVAL]
char sciext[SZ_FNAME], wavext[SZ_FNAME], detsec[SZ_KWVAL], maskname[SZ_KWVAL]
pointer wavephu, wavesci, wavedbim, outphu, outwav, outdat
pointer ap, aps, cf, id, sf, iarr, jarr
int col, line, line1, line2, nap, naps, firstline[1000], lastline[1000], axis
int ext, nsci, nslit, lap, tap, block, nlim, nlimits, c1, c2, bin[2]
int  i, j, il, ih, jl, jh, ni, nj, nio, njo, njr, npix, st, idx
real coord, step, y, limit[1000], ulim, llim, wl
bool debug, gendbim

int clgeti(), imgeti(), sscan(), parsesec(), parseints(), imaccess()
int gstrmatch(), strsearch()
pointer immap(), gt_init(), imgl2r(), impl2r(), imps2r()
real apcveval()
double id_fitpt()
bool new_block()

begin

  debug = false

  # Get task parameters:
  call clgstr("waveref", waveref, SZ_FNAME)
  call clgstr("outimage", outimage, SZ_FNAME)
  call clgstr("database", database, SZ_FNAME)
  ncols = clgeti("ncols")
  xord = clgeti("xorder")
  yord = clgeti("yorder")

  if (xord >= ncols)
    call error(1, "gfunexwl: xorder must be < ncols\n")

  # Strip any .fits extensions in order to construct the database names:
  if(gstrmatch (waveref, "{.fit?$", c1, c2) > 0)
    waveref[c1] = EOS

  # Check we can open the arc and get the number of sci extensions / slits.
  call gin_psuf(waveref, "[0]", wavephuname)
  iferr(wavephu = immap(wavephuname, READ_ONLY, NULL)) {
    call error(1, "gfunexwl: cannot open image specified by waveref\n")
  }
  iferr(nsci = imgeti(wavephu, "NSCIEXT")) {
    call imunmap(wavephu)
    call error(1, "gfunexwl: failed to read NSCIEXT keyword\n")
  }
  iferr(call imgstr(wavephu, "MASKNAME", maskname, SZ_KWVAL)) {
    call imunmap(wavephu)
    call error(1, "gfunexwl: failed to read MASKNAME keyword\n")
  }
  call imunmap(wavephu)

  # Create empty output image to hold output SCI extensions:
  if (imaccess(outimage, READ_ONLY)==YES)
    call error(1, "gfunexwl: outimage already exists")
  outphu = immap(outimage, NEW_IMAGE, NULL)
  call imunmap(outphu)

  # Determine the starting slit number from the mask name:
  if(strsearch(maskname, "_slitb_mdf") > 0)
    nslit = 2
  else
    nslit = 1

  # Buffer for aperture structs, starting at 100 as in apall (gets realloced):
  call malloc (aps, 100, TY_POINTER)

  # Loop over the slits:
  for (ext=1; ext <= nsci; ext=ext+1) {

    # Construct the database filenames corresponding to this slit:
    # (currently the GMOS scripts seem to number the apextract & identify
    # database files differently for IFU-B data, one ap..._2 and the other
    # id..._001 for the same slit!)
    call sprintf(dbsuff, SZ_FNAME, "_%d")
      call pargi(nslit)
    call gin_psuf(waveref, dbsuff, tracedbname)
    call sprintf(dbsuff, SZ_FNAME, "_%03d")
      call pargi(ext)
    call gin_psuf(waveref, dbsuff, wavedbname)

    # Get dimensions of original image section from which the apertures were
    # extracted, from the arc's DETSEC keyword, and recreate the per-slit
    # simple FITS file matching the database, to run code from onedspec:
    call sprintf(sciext, SZ_FNAME, "[SCI,%d]")
      call pargi(ext)
    call gin_psuf(waveref, sciext, wavesciname)
    wavesci = immap(wavesciname, READ_ONLY, NULL)  # default err names SCI ext
    ni = IM_LEN(wavesci, 1)
    nj = IM_LEN(wavesci, 2)
    call imgstr(wavesci, "DETSEC", detsec, SZ_KWVAL)
    if (parsesec(detsec, il, ih, jl, jh) != 0)
      call error(1, "gfunexwl: invalid DETSEC value\n")
    call imgstr(wavesci, "CCDSUM", ccdsum, SZ_KWVAL)
    if (parseints(ccdsum, bin, 2, st) != 0 || st != 2)
      call error(1, "gfunexwl: invalid CCDSUM value\n")
    if (imaccess(wavedbname, READ_ONLY) == YES) gendbim = false
    else {
      gendbim = true
      wavedbim = immap(wavedbname, NEW_COPY, wavesci)
      do line = 1, nj
        call amovr(Memr[imgl2r(wavesci, line)],
                   Memr[impl2r(wavedbim, line)], ni)
      call imunmap(wavedbim)
    }
    call imunmap(wavesci)

    # Check that DETSEC was parsed as expected:
    if (debug) {
      call printf("size %d %d, detsec %d %d %d %d\n")
        call pargi(ni)
        call pargi(nj)
        call pargi(il)
        call pargi(ih)
        call pargi(jl)
        call pargi(jh)
      call printf("xbin %d, ybin %d\n")
        call pargi(bin[1])
        call pargi(bin[2])
    }

    # Make sure the number of DETSEC columns is consistent with the actual
    # extracted image:
    nio = (ih-il+1) / bin[1]   # divides with no remainder for binned data
    njo = (jh-jl+1) / bin[2]   # ditto
    if (ni != nio) {
      if (gendbim) call imdelete(wavedbname)
      call error(1, "gfunexwl: DETSEC doesn't match number of columns\n")
    }

    # Create output SCI extension for this slit:
    call sprintf(wavext, SZ_FNAME, "[WAV,%d,append]")
      call pargi(ext)
    call gin_psuf(outimage, wavext, outwavname)
    outwav = immap(outwavname, NEW_IMAGE, NULL)
    call setdim2d(outwav, nio, njo)
    call imastr(outwav, "DETSEC", detsec)  # propagate original image sec.
    call imastr(outwav, "CCDSUM", ccdsum)

    # Get the output buffer:
    outdat = imps2r(outwav, 1, nio, 1, njo)

    # Ensure that any unwritten pixels will default to zero:
    call aclrr(Memr[outdat], nio*njo)

    # Read aperture database for this slit (returns populated aps & naps):
    naps = 0
    # call ap_alloc(ap);
    call ap_dbread(tracedbname, aps, naps)
    # call ap_free(ap);

    # Get dispaxis from first aperture:
    axis = AP_AXIS(Memi[aps])

    # Initialize identify data structure that we can read the DB into:
    call id_init(id)
    call ic_open (ID_IC(id))
    call strcpy(database, ID_DATABASE(id), SZ_FNAME)
    call strcpy(wavedbname, ID_IMAGE(id), SZ_FNAME)

    # Identify calls these to read image headers etc. corresponding to the
    # database entries. These seem partly to populate the data structure with
    # info. for the reference aperture, which we've hard wired to "first line".
    call id_map(id)       # obscure error if file missing but checked above
    call id_gdata(id)
    
    # Make sure the solution matches the image dimensions (shouldn't fail):
    if (ID_NPTS(id) != ni || ID_MAXLINE(id,1) != nj)
      call error(1, "gfunexwl: wavelength solution doesn't match image dims\n")

    # Determine control columns:
    step = (ni-1) / real(ncols-1)
    
    # Check (some of) what we got back from the database:
    if (debug) {
      call eprintf("ID_IMAGE %s\n")
      call pargstr(ID_IMAGE(id))
      call eprintf("ID_SECTION %s\n")
      call pargstr(ID_SECTION(id))
      call eprintf("ID_DATABASE %s\n")
      call pargstr(ID_DATABASE(id))
      call eprintf("ID_COORDLIST %s\n")
      call pargstr(ID_COORDLIST(id))
      call eprintf("ID_SAVEID %s\n")
      call pargstr(ID_SAVEID(id))
      call eprintf("ID_LINE %d %d\n")
      call pargi(ID_LINE(id,1))
      call pargi(ID_LINE(id,2))
      call eprintf("ID_MAXLINE %d %d\n")
      call pargi(ID_MAXLINE(id,1))
      call pargi(ID_MAXLINE(id,2))
      call eprintf("ID_AP %d %d\n")
      call pargi(ID_AP(id,1))   # This defaults to the middle line
      call pargi(ID_AP(id,2))   # This seems always just to be set to 1
      # APS is a pointer
      call eprintf("ID_NSUM %d %d\n")
      call pargi(ID_NSUM(id,1))
      call pargi(ID_NSUM(id,2))
      call eprintf("ID_MAXFEATURES %d\n")
      call pargi(ID_MAXFEATURES(id))
      call eprintf("ID_FTYPE %d\n")
      call pargi(ID_FTYPE(id))
      call eprintf("ID_NPTS %d\n")
      call pargi(ID_NPTS(id))
      call eprintf("\n")
    }

    # Loop over the apertures and identify image regions corresponding to
    # distinct blocks. This makes the whole thing a bit GMOS-IFU-specific, but
    # it's a lot easier than handling one block at a time and then trying to
    # reconstruct which image rows match which fits most closely after the fact.
    lap = Memi[aps]
    block = 0   # NB. this is defined by aperture ID across slits, not from 1
    firstline[1] = 1
    nlimits = 1; limit[1] = 0.5
    for (nap=1; nap <= naps; nap=nap+1) {

      ap = Memi[aps+nap-1]

      if (new_block(ap, block)) {

        # Increment block gap counter and note bounding input rows for later:
        lastline[nlimits] = nap-1
        nlimits = nlimits + 1
        firstline[nlimits] = nap

        # Average the nominal aperture centres (at the ref. column) above/below
        # the gap to get the boundary line between blocks before extraction:
        limit[nlimits] = 0.5 * (AP_CEN(ap, axis) + AP_CEN(lap, axis))

      }

      lap = ap

    }
    lastline[nlimits] = nap-1    # remember final loop increment gives naps+1
    nlimits = nlimits + 1
    limit[nlimits] = njo + 0.5   # last limit is top of image

    # Check that block gaps were identified correctly:
    if (debug)
      for (line=1; line <= nlimits; line=line+1) {
        call eprintf("%d %.3f \n")
          call pargi(line)
         call pargr(limit[line])
      }

    # # Check block identifications visually:
    # nlim=2
    # for (line=1; line <= njo; line=line+1) {
    #   if (real(line) > limit[nlim]) nlim=nlim+1
    #   call amovkr(real(nlim-1), Memr[outdat+(line-1)*nio], nio)
    # }

    # Ensure realloc for co-ord arrays below gets null pointer 1st time:
    iarr = NULL; jarr = NULL

    # Iterate over the blocks, fit a surface lambda(i,j) to control points
    # and evaluate the fit onto each output pixel nearest the block:
    for (nlim=1; nlim < nlimits; nlim=nlim+1) {

      # Range of row co-ordinates over which this fit is defined:
      llim = limit[nlim]
      ulim = limit[nlim+1]

      # Initialize gsurfit data structure:
      call gsinit (sf, GS_CHEBYSHEV, xord, yord, NO, 1., real(nio), llim, ulim)

      if (debug) {
        call printf("%d - %d:  %.3f - %.3f\n")
          call pargi(firstline[nlim])
          call pargi(lastline[nlim])
          call pargr(limit[nlim])
          call pargr(limit[nlim+1])
      }

      # Iterate over input arc rows & evaluate aperture fits at control
      # columns to get the samples to interpolate:
      line1 = firstline[nlim]
      line2 = lastline[nlim]
      nap = 1   # count aps separately in DB (allows inimage to be a subset)
      for (line=line1; line <= line2; line=line+1) {

        # Read this aperture's wavelength solution from the database:
        ID_AP(id,1)=Memi[ID_APS(id)+line-1]
        call id_dbread (id, wavedbname, ID_AP(id,1), NO, NO)

        # Find next apextract structure with same AP_ID reported by identify:
        # - apextract appears to sort these by increasing AP_ID, so at each
        #   iteration we scan forward to find the next one.
        for (; nap <= naps; nap=nap+1) {
          ap = Memi[aps+nap-1]
          if (AP_ID(ap) == ID_AP(id,1)) break;
            else if (nap == naps)
              call error(1, "gfunexwl: apextract database missing aperture from identify")
        }

        # call printf("%d\n")
        #   call pargi(AP_ID(ap))

        # Loop over control columns:
        for (coord=1.0; coord < real(ni+1); coord=coord+step) {

          # This looks to be how identify rounds to an integer ref. column:
          col = int(coord)

          # Evaluate fibre trace & wavelength at this column:
          y = AP_CEN(ap, axis) + apcveval(AP_CV(ap), real(col))
          wl = real(id_fitpt(id, double(col)))

          # # Test evaluated fit prior to interpolating:
          # # Memr[outdat+(int(y)-1)*nio+col-1] = wl
          # call printf("block %d ap %d col %d line %.3f wl %.3f\n")
          #  call pargi(nlim)
          #  call pargi(AP_ID(ap))
          #  call pargi(col)
          #  call pargr(y)
          #  call pargr(wl)

          # Accumulate points to be fitted:
          call gsaccum(sf, real(col), y, wl, 1.0, WTS_UNIFORM)

        } # End loop over control columns.

        # At the next line, start database AP_ID search from next aperture
        # since we've already done the current one:
        nap = nap+1

      }  # End iteration over arc rows / apertures for this block

      # Fit interpolant surface to the control points for this block:
      call gssolve(sf, st)

      # Clean up if the fit somehow failed (eg. ncols too small):
      if (st != OK) {
        call gsfree(sf)
        call imunmap(outwav)
        if (gendbim) call imdelete(wavedbname)
        call error(1, "gfunexwl: surface fit failed")
      }

      # Derive region dims & construct output co-ordinate arrays:
      njr = int(ulim) - int(llim)   # round limits inwards (+1-1 cancel)
      npix = nio * njr
      call realloc(iarr, npix, TY_REAL)
      call realloc(jarr, npix, TY_REAL)
      idx=0
      for (j=int(llim)+1; j <= int(ulim); j=j+1)
        for (i=1; i <= nio; i=i+1) {
          Memr[iarr+idx] = real(i)
          Memr[jarr+idx] = real(j)
          idx = idx+1
        }

      # Evaluate the fit onto the applicable output region:
      call gsvector(sf, Memr[iarr], Memr[jarr],
        Memr[outdat+int(llim)*nio], npix)

      # Free the gsurfit struct after interpolating this co-ord region:
      call gsfree(sf)

    }  # end iteration over blocks & wavelength interpolation onto output

    # Free co-ordinate arrays & IRAF database structs:
    call mfree(iarr, TY_REAL)
    call mfree(jarr, TY_REAL)
    call id_free(id)
    call mfree(aps, TY_POINTER)

    # Close output SCI extension:
    call imunmap(outwav)

    # Remove the temporary simple FITS file matching the onedspec database:
    if (gendbim) call imdelete(wavedbname)

    # Increment slit number for identify database matching next SCI extension:
    nslit = nslit+1

  }  # end loop over slits

end


# Function to decide how apertures are grouped into blocks:
#
# The first time this is called, lastblock should be supplied with a value of
# 0 so it gets initialized correctly.
#
bool procedure new_block(ap, lastblock)

pointer ap
int block, lastblock

begin

  block = (AP_ID(ap)-1) / 50 + 1

  if (lastblock == 0) lastblock = block

  if (block > lastblock) {
    lastblock = block
    return true
  }
  return false

end


# Parse a 2D IRAF image section in the specific format "[x1:x2,y1:y2]" used
# for detector section keywords:
int procedure parsesec(str, il, ih, jl, jh)

char str[ARB]
int il, ih, jl, jh

int i, len, nch
int strlen(), ctoi()

begin

  # This is the brute force parsing method, as it's reliable and I don't have
  # time to start fixing up missing SPP functionality 30 yrs down the line...

  len = strlen(str)
  i = 1

  # Skip any leading whitespace:
  while (i <= len && IS_WHITE(str[i])) i=i+1

  # Skip any opening square bracket:
  if (str[i] == '[') i=i+1

  # First numeric limit:
  if (ctoi(str, i, il) < 1) return 1

  # Require the next character to be a colon:
  if (str[i] == ':') i=i+1
  else return 1

  # Second numeric limit:
  if (ctoi(str, i, ih) < 1) return 1

  # Require the next character to be a comma:
  if (str[i] == ',') i=i+1
  else return 1

  # Third numeric limit:
  if (ctoi(str, i, jl) < 1) return 1

  # Require the next character to be a colon:
  if (str[i] == ':') i=i+1
  else return 1

  # Fourth numeric limit:
  if (ctoi(str, i, jh) < 1) return 1

  # Skip any closing square bracket:
  if (str[i] == ']') i=i+1

  # Skip any trailing whitespace:
  while (i <= len && IS_WHITE(str[i])) i=i+1

  # Success==0 if we've now reached the end of the string:
  if (str[i] == EOS) return 0
  return 1

end


# Parse space-separated integers from a string (CCDSUM):
int procedure parseints(str, iv, nmax, n)

char str[ARB]
int iv[ARB]    # must be hold as many vals as found, up to nmax
int nmax, n

int i, len
int strlen(), ctoi()

begin

  len = strlen(str)
  i = 1
  n = 1

  # Skip any leading whitespace:
  while (i <= len && IS_WHITE(str[i])) i=i+1

  # Loop over space-separated integers:
  while (i <= len && n <= nmax) {

    # Get next integer (anything else is an error):
    if (ctoi(str, i, iv[n]) < 1) return 1
    n = n+1

    # Skip any trailing whitespace:
    while (i <= len && IS_WHITE(str[i])) i=i+1

  }

  # Decrement count started at 1 instead of 0:
  n = n-1

  # Success==0 if we've now reached the end of the string
  # (it's an error if nmax was reached with chars remaining):
  if (str[i] == EOS) return 0
  return 1

end


procedure setdim2d(outext, ni, nj)

pointer outext
int ni, nj

begin

  # Set output data type:
  call imputi(outext, "i_pixtype", TY_REAL)

  # Set dimensionality:
  call imputi(outext, "i_naxis", 2)
  call imputi(outext, "i_naxis1", ni)
  call imputi(outext, "i_naxis2", nj)

end

