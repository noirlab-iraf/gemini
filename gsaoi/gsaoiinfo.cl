# Copyright(c) 2012-2017 Association of Universities for Research in Astronomy, Inc.

procedure gsaoiinfo

char version {"20Jul2017", prompt="Package version date"}

begin

    help("gsaoiinfo", file_templat-, all-, parameter="all", section="all",
        option="help", page+, nlpp=59, lmargin=1, rmargin=72, curpack="AskCL",
        device="terminal", helpdb="helpdb")

end
