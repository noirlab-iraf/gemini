# Copyright(c) 2000-2017 Association of Universities for Research in Astronomy, Inc.

procedure niriinfo

char version {"20Jul2017",prompt="Package version date"}

begin

    help("niriinfo", file_templat-, all-, parameter="all", section="all",
        option="help", page+, nlpp=59, lmargin=1, rmargin=72, curpack="AskCL",
        device="terminal", helpdb="helpdb")

end
