# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure flamingosinfo

char version {"20Jul2017",prompt="Package version date"}

begin

    help("flamingosinfo", file_templat-, all-, parameter="all", section="all",
        option="help", page+, nlpp=59, lmargin=1, rmargin=72, curpack="AskCL",
        device="terminal", helpdb="helpdb")

end
