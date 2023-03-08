# Copyright(c) 2010-2017 Association of Universities for Research in Astronomy, Inc.

# F2INFO - Generic information on FLAMINGOS-2 data and reduction

procedure f2info

char version {"20Jul2017", prompt = "Package version date"}

begin

help("f2info", file_templat-, all-, parameter="all", section="all", 
    option="help", page+, nlpp=59, lmargin=1, rmargin=72, curpack="AskCL", 
    device="terminal", helpdb="helpdb")

end
