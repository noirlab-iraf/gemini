# Copyright(c) 2000-2012 Association of Universities for Research in Astronomy, Inc.

procedure gemhedit (image, field, value, comment)

# Edit FITS header keyword
# Should be safe and never return status = 1 unless image does not exist
#
# Version  May 31, 2001   IJ, v1.2 release
#          Mar 03, 2002   IJ, v1.3 release
#          Sep 20, 2002   IJ, v1.4 release
#          Feb 03, 2003   IJ  fix segmentation violation bug under Linux7.3
#          Aug 11, 2003   KL, Added 'addonly-' parameter to hedit calls
#          Jan 09, 2009   JH, Propagate comment update for existing keyword
#          Feb 03, 2010   EH, Remove original comment before writing new one
#                             (longer comments remain), add delete option

char    image   {prompt = "Input image"}
char    field   {prompt = "Keyword name"}
char    value   {prompt = "Value"}
char    comment {prompt = "Comment to add for keyword"}
bool    delete  {no, prompt = "Delete keyword(s)?"}
char    upfile  {"", prompt = "Update Command file"}
int     status  {0,prompt = "Exit status"}

begin

    char    l_image, l_field, l_upfile, l_after, l_before, test_str
    bool    l_keyexist, l_delete, l_add, l_addonly, l_update, l_show
    int     tlen, test_fn, test_len, comma_pos
    real    test_val
    string  l_value
    struct l_comment, sdate

    cache ("keypar", "nhedit")

    l_image = image ; l_field = field ; l_value = value ;
    l_comment = comment ; l_delete = delete
    l_upfile = upfile

    # Initialize parameters
    status = 0

    # Default nhedit parameters - that don't get changed
    l_after = ""
    l_before = ""
    l_update = yes
    l_add = yes # gets changed if delete+
    l_addonly = no
    l_show = no
    l_keyexist = no

    # Check Image exists
    if (no == imaccess(l_image)) {
        print("ERROR - GEMHEDIT: Image "//l_image//" not found")
        goto crash
    }

    # Command file has been supplied
    if (l_upfile != "") {
        # No checks of the file and any other parameters are made
        if(access(l_upfile)) {

            # These parameters need to reset for a command file
            l_field = ""
            l_value = ""
            l_comment = ""

        } else {
            print ("ERROR - GEMHEDIT: upfile "//l_upfile//" not found")
            goto crash
        }
    }

    # Make sure the field name is in upper case. Don't capitalise
    # keywords beginning with i_ EH

    if (substr(l_field, 1, 2) != "i_") {
        print (l_field) | \
            translit ("STDIN", "a-z", "A-Z", delete-, collapse=no) | \
            scan (l_field)
    }

    # Cut keyword to eight if it is longer - only if not deleteing (to allow
    # comma separated keyword lists to be deleted - MS
    if (!l_delete && (strlen(l_field) > 8)) {
        l_field = substr(l_field,1,8)

        # Check for any commas - shouldn't really get this unless it's a list
        # which is likely only to be used for deleting! - MS
        comma_pos = stridx(",",l_field)
        if (comma_pos > 0) {
            l_field = substr (l_field,1,comma_pos-1)
        }
    }

    # Check if keyword already exists
    #
    # A keyword can exist *and* have an empty value. The original call to
    # hselect didn't work properly if a keyword existed but had an empty
    # value. Previously, in this case, it was assumed that the keyword
    # did not exist, which meant that the nhedit call with add+ was used.
    # This, unfortunately, wrote the same keyword to the header for a
    # second time. Using keypar instead allows you to test whether the
    # keyword actually exists, rather than whether it has a value. EH
    #hselect(l_image,l_field,yes) | scan(l_keycheck)

    # If deleteing change some parameters - ignore if it exists or not - this
    # allows comma separated lists of keywords to be deleted. Add must be
    # switched off. - MS
    if (l_delete) {

        l_value = "."
        l_comment = "."
        l_add = no

    } else {

        keypar(l_image, l_field, silent+)
        if (keypar.found) {
            l_keyexist = yes
        }
    }

    if ((l_keyexist == yes) && (l_comment == "")) {
        l_comment = "." # leave original comment in place

    } else if (l_keyexist && l_comment != ".") {

        # Keyword exists and comment is not an empty string or a dot
        # Update comment to contain extra spaces if needed to overwrite the
        # original comment

        test_str = str(l_value)
        test_fn = fscan (test_str,test_val)

        # The general space for a value in a card is 20 characters - at least
        if (test_fn == 0) {
            # It's a string (include quotes in length)
            test_len = 18
        } else {
            # Number
            test_len = 20
        }

        # Test if value length is greater than the default value length
        tlen = strlen(test_str)
        if (tlen > test_len) {
            # Card length is 80
            tlen = 80 - tlen
        } else {
            tlen = 80 - test_len
        }
        # Remove the remaing length of standard characters:
        #     keyword(8), "=", " / "
        tlen = tlen - 12

        if (strlen(l_comment) < tlen) {
            printf ("%-"//tlen//"s\n", l_comment) | scan (l_comment)
        }
    }

SKIP_CHECKS:

    # Because nhedit was removed from fitsutil and we do want
    # to get backward compatibility for a little while, the first
    # block checks if fitsutil.nhedit exists and if so use it,
    # if not, the following block checks if imutil.nhedit exists and if
    # so use that one.  The final block just catches error cases where
    # no known nhedit is found.
    #
    # It is a lot of code duplication, but hopefully in a couple release
    # we will be able to remove the backward compatible block to fitsutil.
    # When that day comes it will be easy to remove it all in one shot.
    # -- KL,  13 Jan 2010

    if (deftask("fitsutil.nhedit")) {

        if (l_upfile == "") {
            l_upfile = "NULL"
        }

        # This is to keep fits standard if a value of INDEF is used
        # Can't seem to get this to work any other way - MS
        if ((isindef(l_value) || l_value == "INDEF") && !l_keyexist) {
            # Have to add a string vale to
            fitsutil.nhedit (l_image, l_field, "CHANGEME", comment=l_comment, \
                comfile=l_upfile, after=l_after, before=l_before, \
                update=l_update, add=l_add, addonly=l_addonly, \
                delete=l_delete, verify-, show=l_show)
            l_comment = "."
        }

        fitsutil.nhedit (l_image, l_field, l_value, comment=l_comment, \
            comfile=l_upfile, after=l_after, before=l_before, \
            update=l_update, add=l_add, addonly=l_addonly, delete=l_delete, \
            verify-, show=l_show)

    } else if (deftask("imutil.nhedit")) {

        # This is to keep fits standard if a value of INDEF is used
        # Can't seem to get this to work any other way - MS
        if ((isindef(l_value) || l_value == "INDEF") && !l_keyexist) {
            # Have to add a string vale to
            imutil.nhedit (l_image, l_field, "CHANGEME", comment=l_comment, \
                comfile=l_upfile, after="", before="", update=l_update, \
                add=l_add, addonly=l_addonly, delete=l_delete, \
                verify-, show=l_show)
            l_comment = "."
        }

        imutil.nhedit (l_image, l_field, l_value, comment=l_comment, \
            comfile=l_upfile, after="", before="", update=l_update, \
            add=l_add, addonly=l_addonly, delete=l_delete, \
            verify-, show=l_show)

    } else {
        print ("ERROR - GEMHEDIT: No NHEDIT found")
        goto crash
    }

    goto clean

crash:

    status = 1

clean:

    # Reset upfile to ""
    upfile = ""

end
