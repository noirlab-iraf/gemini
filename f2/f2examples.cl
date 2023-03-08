# Copyright(c) 2010-2012 Association of Universities for Research in Astronomy, Inc.

procedure f2examples (example)

char example {enum="imaging|longslit|MOS", prompt="Example to print"}

begin

    char l_example
    l_example = example

    if (l_example == "imaging") 
        page("f2$doc/f2_imaging_example.cl")

    if (l_example == "longslit")
        page("f2$doc/f2_longslit_example.cl")

    if (l_example == "MOS")
        page("f2$doc/f2_mos_example.cl")

end
