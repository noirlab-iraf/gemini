
               Installation of the Gemini IRAF package v1.14
               ---------------------------------------------

[0] All users are requested to use AstroConda 
    (https://astroconda.readthedocs.io/) to install IRAF and the Gemini IRAF 
    package. Full installation instructions are available on the AstroConda
    website.  Quick installation instructions can be found at
    http://www.gemini.edu/node/12665

    The Gemini IRAF package v1.14 should run fine under Ureka.  We still 
    encourage users to switch to AstroConda at their earliest convenience.

    The Gemini IRAF package v1.14 should run under both IRAF and PyRAF
    environments, however it has been tested only on PyRAF.

    WARNING!! The Gemini IRAF package v1.14 does NOT work reliably with the
    32-bit or 64-bit version of IRAF v2.16 as distributed by NOAO.

    The Gemini IRAF package was tested under IRAF v2.16 from AstroConda
    (as distributed in June 2017)


AstroConda Installation
-----------------------
[1] The Gemini IRAF package v1.14 is distributed with AstroConda.  Following
    the "AstroConda with IRAF installation" instructions will install 
    Gemini IRAF v1.14.
    
    Remember to run "mkiraf" to configure IRAF before using it for the 
    first time.

    No further action is required.
    
    The package can also be installed manually from the tarball.  The 
    installation instructions for that follow below.


Manual Installation from tarball
--------------------------------
[2] The Gemini IRAF v1.14 tarball is distributed via the Gemini web pages:

        http://www.gemini.edu/sciops/data-and-results/processing-software

    The file gemini_v114.tar.gz is to be downloaded
    
    The gemini_readme.txt file contains these instructions.

[3] Create a directory to contain the Gemini IRAF package files. This 
    directory should be outside the IRAF directory tree.  We will call
    that directory "my_gemini_v114".
    
        % mkdir my_gemini_v114

[4] Unpack the tar file in your installation directory

	% cd my_gemini_v114
	% tar xvzf <path_to>/gemini_v114.tar.gz

[5] No need to compile anything

    32-bit Linux and Mac OS X binaries, compatible with AstroConda are
    already included with the package.  You do not need to compile anything.

[6] Configure IRAF to pick up the new package.

    In your iraf home directory (commonly ~/iraf), create a "loginuser.cl"
    file that will contain the package configuration.  The iraf home directory 
    is where you ran "mkiraf", and where the login.cl is located.
	
	% vi loginuser.cl  (or use your favorite editor)
	
    Add the following lines:
        reset gemini=<full_path_to>/my_gemini_v114/	
        task gemini.pkg=gemini$gemini.cl
        keep
    
    The trailing "/" on the first line is important, don't forget it.

[7] If updating an older installation of the Gemini IRAF package, it is
    recommended that users initialize their uparm directories by typing
    "rm uparm/*" in their iraf home directory. NOTE: if you wish to 
    make a note of any stored parameters, please do so before running this command.

[8] Please use the Gemini HelpDesk for submitting questions
    http://www.gemini.edu/sciops/helpdesk/helpdeskIndex.html
