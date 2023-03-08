# Copyright(c) 2002-2005 Association of Universities for Research in Astronomy, Inc.
#
# Subroutines for checking and manipulating image name strings, consistent
# with the behaviour of the gemini tasks.
#
# Version     Sep  5, 2002  JT
# ----

include <gemini.h>
include <ctype.h>


procedure gin_gfext(name, fext)

# Extract the file type extension from an image name string

# The input may optionally include an image extension or section suffix.
# If there is no file extension or it is not recognized by IRAF as an image
# format, an empty string is returned.

char name[SZ_FNAME]                     # I: full image name
char fext[SZ_FNAME]                     # O: file extension

char cluster[SZ_FNAME], dummyname[SZ_FNAME]
int n, clustlen, fextlen

int strlen(), imaccess()

begin

  # Get the cluster (file) name and length:
  call imgcluster(name, cluster, SZ_FNAME)
  clustlen = strlen(cluster)

  # Find last dot in cluster name:
  for (n=clustlen; n > 0; n=n-1)
    if (cluster[n]=='.') break
  
  # If there is a dot, check whether it begins a known image ext:
  if (n > 0) {
	
    # Copy extension name, from last dot to the end:
    call strcpy(cluster[n+1], fext, SZ_FNAME)
	fextlen = strlen(fext)
	
    # Remove any trailing whitespace after extension:
    for (n=fextlen; n>0 && IS_WHITE(fext[n]); n=n-1);
    if (n < fextlen) fext[n+1] = EOS
	
    # Check the extension name against those recognized by IRAF and
	# reset it to an empty string if not found:
    call strcpy("name.", dummyname, SZ_FNAME)  # (dummy root)
    call strcat(fext, dummyname, SZ_FNAME)
    if (imaccess(dummyname, NEW_IMAGE)==NO) fext[1]=EOS

  } 
  
  # Otherwise (no dot), set the file extension name to an empty string:
  else fext[1]=EOS
  
end


procedure gin_gsuf(name, suffix)

# Extract the image extension/section suffix from an image name string

# The suffix may contain a cluster index or kernel section (specifying an
# image extension within the host file) and/or an image section. If none
# of these are present, the string returned is blank. See imio$imparse.x
# for more info. The suffix is not checked here for validity.

char name[SZ_FNAME]                        # I: full image name
char suffix[SZ_FNAME]                      # O: suffix (if present)

char cluster[SZ_FNAME]
int clustlen

int strlen()

begin

  # Get the cluster (file) name and length:
  call imgcluster(name, cluster, SZ_FNAME)
  clustlen = strlen(cluster)
  
  # Copy any image extension/section after the cluster name:
  call strcpy(name[clustlen+1], suffix, SZ_FNAME)

end


procedure gin_groot(name, root)

# Extract the root filename from an image name string

# Any file extension and image extension or section suffix are removed,
# but the full path is preserved.

# Name and root can be the same.

char name[SZ_FNAME]                     # I: full image name
char root[SZ_FNAME]                     # O: root filename (cluster-extn)

char fext[SZ_FNAME]
int n, clustlen

int strlen(), strcmp()

begin
  
  # Get the full cluster name:
  call imgcluster(name, root, SZ_FNAME)
  clustlen = strlen(root)

  # Check whether the cluster name includes a recognized file ext:
  call gin_gfext(root, fext)
  
  # If there is an extension, delete everything after the last dot:
  if (strcmp(fext, "")!=0) {
    for (n=clustlen; n > 0; n=n-1) if (root[n]=='.') break
	root[n]=EOS
  }
  
end


procedure gin_afext(name, defext, outname, find)

# Make sure a filename string includes a file type extension.

# If no image extension name recognized by IRAF is already present, one is
# inserted between the root filename and any image extension/section suffix.
# If find==YES, the extension is set by searching for the image on disk.
# If find==NO, or if the image is not found or ambiguous, the extension
# defaults to "defext", which should be set to GEM_DEFEXT or GEM_DMULTEXT by
# the calling procedure. In the cases of (unexpected) non-existence and
# ambiguity, the main program should generate a separate error, since the
# default is really just a placeholder which could be confusing for the user.

# Name and outname can be the same.

char name[SZ_FNAME]                     # I: Filename, with/without ext
char defext[SZ_FNAME]                   # I: Default ext to use if none
char outname[SZ_FNAME]                  # O: Filename with ext
int find                                # I: Look for extension on disk?

char cluster[SZ_FNAME], fext[SZ_FNAME], suffix[SZ_FNAME], tempname[SZ_FNAME]
int n, ikistat, clustlen

int strlen(), iki_access(), strcmp(), gte_isblank()
bool itob()

begin

  # We have to use the low-level iki_access() function to get the image
  # name extension corresponding to a file on disk. If this changes, the
  # procedure will break when find==YES.

  # Get the cluster name (filename, including file ext):
  call imgcluster(name, cluster, SZ_FNAME)
  clustlen = strlen(cluster)

  # Shift any trailing whitespace from cluster name to suffix:
  for (n=clustlen; n>0 && IS_WHITE(cluster[n]); n=n-1);
  if (n < clustlen) cluster[n+1] = EOS
  clustlen = n

  # Copy any image extension/section after the cluster name:
  call strcpy(name[clustlen+1], suffix, SZ_FNAME)
 
  # Get any existing file extension:
  call gin_gfext(cluster, fext)

  # If there is no extension, add one:
  if (strcmp(fext, "")==0) {

    # If find==YES, try to determine extension by looking on disk:
	if (itob(find)) {

      # Determine image (not ext) existence and file extension:
	  ikistat = iki_access(cluster, tempname, fext, READ_ONLY)

	  # Append actual extension if found, otherwise (image doesn't
	  # exist or name ambiguous) append default extension:
	  if (ikistat > 0) {
		call strcat(".", cluster, SZ_FNAME)
		call strcat(fext, cluster, SZ_FNAME)
      }
	  else if (gte_isblank(defext)==NO) {
		call strcat(".", cluster, SZ_FNAME)
        call strcat(defext, cluster, SZ_FNAME)
      }
	  
	} # end if (find==YES)
	
	# Otherwise (find==NO), set extension to the default
	else if (gte_isblank(defext)==NO) {
	  call strcat(".", cluster, SZ_FNAME)
	  call strcat(defext, cluster, SZ_FNAME)
    }

  } # end if (no existing extension)

  # Reconstruct the full image name:
  call strcpy(cluster, outname, SZ_FNAME)
  call strcat(suffix, outname, SZ_FNAME)
  
end


procedure gin_psuf(name, suffix, outname)

# Add an image extension or section suffix to an image name, replacing any
# existing suffix.

# Name and outname can be the same.

char name[SZ_FNAME]                     # I: image name
char suffix[SZ_FNAME]                   # I: suffix to append
char outname[SZ_FNAME]                  # O: new image name with suffix

begin

  # Get the cluster name and copy to output:
  call imgcluster(name, outname, SZ_FNAME)
  
  # Append the new suffix:
  call strcat(suffix, outname, SZ_FNAME)
  
end


int procedure gin_ismef(name)

# Determine whether a specified file exists in multi-extension FITS format
# (returns YES or NO)

# The file is considered to be MEF if all these conditions are met:
#   1. the name does not specify a particular image extension
#   2. the file contains a readable extension [0]
#   3. the header contains EXTEND=="T"

# NB. it is not possible to test conformance with the NOST 100-2.0 standard
# rigorously without accessing the header directly. For the time being we
# shall avoid this, since it makes things complicated and it is not clear
# that it is necessary. Note that Space Telescope format has EXTEND=='F'.
# It IS possible, though, to make a bad MEF file with EXTEND=='F' in IRAF.

char name[SZ_FNAME]                        # I: filename to test

char imextname[SZ_FNAME], clustname[SZ_FNAME]

pointer phu
int ismef

int strlen()
pointer immap()
bool imgetb()

begin

  # If name specifies an image extension, return NO:
  call imgimage(name, imextname, SZ_FNAME)
  call imgcluster(name, clustname, SZ_FNAME)
  if (strlen(imextname) != strlen(clustname)) return NO

  # Otherwise, append [0] to cluster name:
  call strcpy(clustname, imextname, SZ_FNAME)
  call strcat("[0]", imextname, SZ_FNAME)
  
  # Set ismef==YES until incompatibility found:
  ismef = YES
  
  # Attempt to open the file (will open the PHU if multi-image cluster):
  iferr(phu = immap(imextname, READ_ONLY, NULL)) return NO

  # Check that EXTEND keyword exists and is true:
  iferr(if(!imgetb(phu, "EXTEND")) ismef = NO) ismef = NO
  
  # Close input file:  
  call imunmap(phu)
  
  return ismef
  
end


procedure gin_make(prefix, name, fext, outname)

# Construct an output filename, based on an input image name

# The input "name" is modified with "prefix" and has the file extension,
# "fext" appended, replacing any existing file extension. Any existing
# image extension or section suffix is deleted.

# Name and outname can be the same.

char prefix[SZ_FNAME]                   # I: prefix to add
char name[SZ_FNAME]                     # I: base root or full filename
char fext[SZ_FNAME]                     # I: file extension to use
char outname[SZ_FNAME]                  # O: derivative filename

char root[SZ_FNAME]

begin

  # Get the root of the base name:
  call gin_groot(name, root)

  # Copy prefix to outname and append the root name:
  call strcpy(prefix, outname, SZ_FNAME)
  call strcat(root, outname, SZ_FNAME)
  
  # Append the file extension:
  call strcat(".", outname, SZ_FNAME)
  call strcat(fext, outname, SZ_FNAME)
  
end


int procedure gin_multype(name)

# Check whether the file extension of an image name corresponds to a
# multi-image format that we can write to (initially just MEF).

char name[SZ_FNAME]                     # I: filename

char fext[SZ_FNAME]

int gte_indict()

begin

  # Get file extension:
  call gin_gfext(name, fext)
  
  # Check against dictionary of valid types:
  return gte_indict(fext, GEM_MULTYPES)
  
end


int procedure gin_valid(name)

# Decide whether input string is a valid name for a new image

# To avoid trouble, wildcards and some other special characters are not
# allowed in the filename, even though IRAF can handle them. Likewise, the
# root name must not be blank (don't want hidden files). The main program
# must check separately whether the file already exists.

char name[SZ_FNAME]                     # I: filename

char clustname[SZ_FNAME], root[SZ_FNAME]

int strcmp(), stridxs()

begin

  # Get filename with and without file extension:
  call imgcluster(name, clustname, SZ_FNAME)
  call gin_groot(name, root)
  
  # Strip surrounding spaces from root name:
  call gte_unpad(root, root, SZ_FNAME)
  
  # Return invalid if root name blank or begins with '-':
  if(strcmp(root, "")==0 || root[1]=='.' || root[1]=='-') return NO
  
  # Search for some other problematic characters in root+ext:
  if(stridxs("@#%&*?<>|,=\"\'", clustname) > 0) return NO
  
  return YES
  
end
