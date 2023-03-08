# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
#
# JT: Prototype version without some of the usual sanity checks (it can
# nevertheless save quite a bit of time & verbosity in processing scripts).

procedure addbpm(inimages, bpm)

string  inimages   {prompt="MEF file(s) to add BPM to"}
string  bpm        {prompt="BPM file"}
string  logfile    {"gemtools.log", prompt = "Log file"}
bool    verbose    {yes, prompt = "Verbose output?"}
struct  *scanfile  {"", prompt="Internal use only"}

begin

    char    l_inimages = ""
    char    l_bpm = ""
    char    l_logfile = ""
    bool    l_verbose

    char    paramstr, tmplist, tmplist2, filename, tmpbpm, dqim
    int     n, nsci, nsci2, nimages

#    paramstr =  "inimages       = "//inimages.p_value//"\n"
#    paramstr += "bpm           = "//bpm.p_value//"\n"

    l_inimages = inimages
    l_bpm = bpm
    l_logfile = logfile
    l_verbose = verbose

    cache ("gloginit", "gemextn", "gemdate")

    tmplist  = mktemp("tmplist")
    tmplist2 = mktemp("tmplist2")
    tmpbpm   = mktemp("tmpbpm")//".fits"

    gemextn (l_inimages, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension", replace="",
        outfile=tmplist, logfile=l_logfile, glogpars="", verbose=l_verbose)
    gemextn ("@"//tmplist, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmplist2, logfile=l_logfile, glogpars="", verbose=l_verbose)
    nimages = gemextn.count
    delete (tmplist, verify-, >& "dev$null")

    scanfile = tmplist2
    while (fscan (scanfile, filename) != EOF) {

        hselect (filename//"[0]", "NSCIEXT", yes) | scan (nsci)

        # Need to check that the number of BPM extensions matches NSCIEXT here
        # This should eventually be updated later to match the CCD identifier
        # in each extension rather than assuming nsciext matches.

        # This should really check whether the science data have been trimmed
        # (if applicable) and either pad the input BPM or at least throw an
        # error if not, to ensure all the dimensions match, but for now it
        # just appends what it's given.

        for (n=1; n <= nsci; n+=1) {

            dqim=filename//"[DQ,"//n//"]"

            # Combine the BPM with any existing DQ in the input file:
            if (imaccess(dqim)) {
                imexpr("a | b", tmpbpm, l_bpm//"[DQ,"//n//"]", dqim, \
                  outtype="ushort", verbose-)
                imcopy(tmpbpm, filename//"[DQ,"//n//",overwrite]")
                imdelete(tmpbpm, verify-, >& "dev$null")
            } else {
                imcopy(l_bpm//"[DQ,"//n//"]", filename//"[DQ,"//n//",append]")
            }
        } # end n <= nsci

    } # end loop over input files

    delete(tmplist2, verify-, >& "dev$null")

end
