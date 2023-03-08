# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

procedure nifsexamples (example)

char    example    {enum = "calibration|science|telluric", prompt="Example to print."}

begin

    char    l_example

    l_example = example

    if (l_example == "calibration") 
        page("nifs$doc/nifs_ifu_cal_example.cl")
    else if (l_example == "science")
        page("nifs$doc/nifs_ifu_science_example.cl")
    else if (l_example == "telluric")
        page("nifs$doc/nifs_ifu_telluric_example.cl")

end
