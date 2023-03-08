# Copyright(c) 2011-2017 Association of Universities for Research in Astronomy, Inc.

# The Gemini IRAF package requires a number of other external IRAF packages.
# The Gemini IRAF package will not load if any of these packages are missing.
# This script checks if external package versions are installed and up to date.
#
# To use this task to check a Gemini IRAF installation, define and execute the
# task by typing the following at the IRAF or PyRAF prompt:
#
#     ecl> task $install_check=gemini$install_check.cl
#     ecl> install_check

procedure install_check

begin

    struct l_check
    bool   isiraf215plus, is_in_ureka
    bool   isiraf216, is_in_astroconda
    string pkg_path, ureka_file, test_IMT
    string astroconda_file
    real   min_version, test_version

    delete ("extern_check.log", verify-, >& "dev$null")
    show iraf | tee ("extern_check.log", out_type="text", append+)
    show IRAFARCH | tee ("extern_check.log", out_type="text", append+)
    show version | tee ("extern_check.log", out_type="text", append+)
    show gemini | tee ("extern_check.log", out_type="text", append+)
    if (defvar("extern")) {
        show extern | tee ("extern_check.log", out_type="text", append+)
    }

    # Check to see if using Ureka - valid use of IRAF 2.16
    is_in_ureka = no
    show iraf | scan (ureka_file)
    ureka_file = ureka_file//"../bin/ur-setup-real"
    if (access(ureka_file)) {
        is_in_ureka = yes
    }
    
    # Check to see if using AstroConda - valid use of IRAF 2.16
    is_in_astroconda = no
    show iraf | scan (astroconda_file)
    astroconda_file = astroconda_file//"../bin/conda"
    if (access(astroconda_file)) {
        is_in_astroconda = yes
    }

    # Check IRAF. The Gemini IRAF package requires IRAF v2.14.1.
    show version | scan (line)
    i = stridx ("V", line)
    j = stridx ("-", line)
    s1 = substr (line, i+1, j-1)
    s2 = substr (line, i+1, strlen(line))
    l_check = ""
    show IRAFARCH | scan (line)
    isiraf215plus = no
    isiraf216 = no

    if (s2 == "2.14.1") {
        print ("IRAF installation 2.14 is no longer supported. FAIL") | \
            tee ("extern_check.log", out_type="text", append+)
#        # Make sure the critical patch is installed
#        dir ("iraf$bin."//line//"/x_tv.e", long+, ncols=0, \
#            maxch=18, sort+, all-) | match ("Jan") | scan (l_check)
#        dir ("iraf$bin."//line//"/x_tv.e", long+, ncols=0, \
#            maxch=18, sort+, all-) | match ("Aug") | scan (l_check)
#        dir ("iraf$bin."//line//"/x_tv.e", long+, ncols=0, \
#            maxch=18, sort+, all-) | match ("Sep") | scan (l_check)
#        if (l_check == "") {
#            print ("Installation of IRAF 2.14.1 patch has \
#                not been installed. FAIL") | \
#                tee ("extern_check.log", out_type="text", append+)
#            dir ("iraf$bin."//line//"/x_tv.e", long+, ncols=0, \
#                maxch=18, sort+, all-) | \
#                tee ("extern_check.log", out_type="text", append+)
#        } else {
#            print ("IRAF installation 2.14.1 with patch is up to date. \
#                PASS") | tee ("extern_check.log", out_type="text", append+)
#        }
    } else if (substr(s2,1,4) == "2.15") {
        print ("IRAF installation 2.15 is no longer supported. FAIL") | \
            tee ("extern_check.log", out_type="text", append+)
#        # iraf 2.15 onwards
#        isiraf215plus = yes
#        if (s2 == "2.15.1a") {
#            print ("IRAF installation 2.15.1a is up to date. PASS") | \
#                tee ("extern_check.log", out_type="text", append+)
#        } else {
#            print ("IRAF v2.15.1a is required. FAIL") | \
#                tee ("extern_check.log", out_type="text", append+)
#        }
    } else if (substr(s2,1,4) == "2.16" && is_in_ureka) {
        isiraf216 = yes

        print ("IRAF installation v2.16 in Ureka valid. PASS") | \
            tee ("extern_check.log", out_type="text", append+)
        print ("  Though we recommend that start using AstroConda") | \
            tee ("extern_check.log", out_type="text", append+)

        show ("use_new_imt") | scan (test_IMT)

        print ("use_new_imt is set to "//test_IMT)
    } else if (substr(s2,1,4) == "2.16" && is_in_astroconda) {
        isiraf216 = yes

        print ("IRAF installation v2.16 in AstroConda is up-to-date. PASS") | \
            tee ("extern_check.log", out_type="text", append+)

        show ("use_new_imt") | scan (test_IMT)

        print ("use_new_imt is set to "//test_IMT)

    } else {
        # IRAF version is neither 2.14.1, 2.15.1a or 2.16 in Ureka
        print ("Incompatible version of IRAF installed. Please install \
            AstroConda (containing v2.16). FAIL")
    }

    # Check to see if the external packages that are required by the Gemini
    # IRAF package are installed

    # Check stsdas. v3.15 is required, since it contains an improject bug fix
    # required for the Gemini IRAF package
    min_version = 3.15
    if (access("stsdas$stsdas.par")) {
        test_version = real(substr(stsdas.version,2,strlen(stsdas.version)))
        if (test_version >= min_version)
            print ("stsdas v"//test_version//" installed. PASS") | \
                tee ("extern_check.log", out_type="text", append+)
        else {
            print ("STSDAS version: "//stsdas.version) | \
                tee ("extern_check.log", out_type="text", append+)
            print ("stsdas less than v"//min_version//". Need to update. \
                FAIL") | tee ("extern_check.log", out_type="text", append+)
        }
    } else
        print ("stsdas package not correctly installed! FAIL") | \
            tee ("extern_check.log", out_type="text", append+)

    # Check tables
    if (access("tables$tables.par")) {
        test_version = real(substr(tables.version,2,strlen(tables.version)))
        if (test_version >= min_version)
            print ("tables v"//test_version//" installed. PASS") | \
                tee ("extern_check.log", out_type="text", append+)
        else {
            print ("TABLES version: "//tables.version) | \
                tee ("extern_check.log", out_type="text", append+)
            print ("tables less than v"//min_version//". Need to update. \
                FAIL") | tee ("extern_check.log", out_type="text", append+)
        }
    } else
        print ("tables package not correctly installed! FAIL") | \
            tee ("extern_check.log", out_type="text", append+)

    # Check fitsutil
    if (access("fitsutil$fitsutil.par")) {
        l_check = ""
        # Set the correct path according to iraf version
        if (isiraf216) {
            pkg_path = "src"
            dir ("fitsutil$src", long+, ncols=0, maxch=18, sort+, all-) | \
                scan (l_check)
            if (l_check == "no") {
                # The src directory was not found. The most likely reason is
                # that the 2.14 version of fitsutil is installed.
                print ("2.14 version of fitsutil package installed, but IRAF \
                    v2.16 detected. Update fitsutil. FAIL.") | \
                    tee ("extern_check.log", out_type="text", append+)
            }
        } else {
            pkg_path = "pkg"
            dir ("fitsutil$pkg", long+, ncols=0, maxch=18, sort+, all-) | \
                scan (l_check)
            if (l_check == "no") {
                # The pkg directory was not found. The most likely reason is
                # that the 2.15 version of fitsutil is installed.
                print ("2.15 version of fitsutil package installed, but IRAF \
                    v2.14 detected. Update fitsutil. FAIL.") | \
                    tee ("extern_check.log", out_type="text", append+)
            }
        }

        l_check = ""
        # The getcmd.x file is dated Aug 26, 2005 for both the 2.14 and 2.15
        # versions of the fitsutil package
        dir ("fitsutil$"//pkg_path//"/getcmd.x", long+, ncols=0, maxch=18, \
            sort+, all-)| match 'Aug 26  2005' | scan (l_check)
        if (l_check != "")
            print ("fitsutil package up to date. PASS") | \
                tee ("extern_check.log", out_type="text", append+)
        else {
            dir ("fitsutil$"//pkg_path//"/fxinsert.x", long+, ncols=0, \
                maxch=18, sort+, all-) | match 'Mar 30  2005' | scan (l_check)
            if (l_check != "")
                print ("fitsutil package v.Mar 30 2005. Consider upgrading. \
                    PASS.") | tee ("extern_check.log", out_type="text", \
                    append+)
            else {
                dir ("fitsutil$"//pkg_path//"/nhedit.x", long+, ncols=0, \
                    maxch=18, sort+, all-) | match 'Mar 11  2005' | \
                    scan (l_check)
                if (l_check != "")
                    print ("fitsutil package v.Mar 11 2005. Consider \
                        upgrading. PASS.") | tee ("extern_check.log", \
                        out_type="text", append+)
                else {
                    dir ("fitsutil$"//pkg_path//"/getcmd.x", long+, ncols=0, \
                        maxch=18, sort+, all-)
                    print ("fitsutil package not up to date. FAIL.") | \
                        tee ("extern_check.log", out_type="text", append+)
                }
            }
        }
    } else
        print ("fitsutil package not correctly installed! FAIL") | \
            tee ("extern_check.log", out_type="text", append+)

    # Check gemini
    if (access("gemini$gemini.par")) {
        if (gemini.verno == "v1.14")
            print ("gemini v1.14 installed. PASS") | \
                tee ("extern_check.log", out_type="text", append+)
        else {
            print ("gemini version: "//gemini.verno) | \
                tee ("extern_check.log", out_type="text", append+)
            print ("gemini package less than v1.14. Need to update. \
                FAIL") | tee ("extern_check.log", out_type="text", append+)
        }
    } else
        print ("gemini package not correctly installed! FAIL") | \
            tee ("extern_check.log", out_type="text", append+)

    # Check machine name
    !uname -a | tee -a extern_check.log

end
