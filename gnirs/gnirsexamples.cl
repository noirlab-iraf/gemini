# Copyright(c) 2004-2012 Association of Universities for Research in Astronomy, Inc.

procedure gnirsexamples (example)

char    example    {enum = "longslit|XD|IFU", prompt="Example to print"}

begin

    char l_example
    l_example = example

    if (l_example == "longslit") 
        page("gnirs$doc/gnirs_longslit_example.cl")

    if (l_example == "XD")
        page("gnirs$doc/gnirs_xd_example.cl")

    if (l_example == "IFU")
        page("gnirs$doc/gnirs_ifu_example.cl")

end
