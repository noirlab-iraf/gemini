# Copyright(c) 2003-2011 Association of Universities for Research in Astronomy, Inc.

procedure midirexamples (example)

char    example {enum="michelle_imaging|trecs_imaging|michelle_longslit|midir_longslit|polarimetry", prompt="Example to print"}

begin

    char    l_example
    
    l_example=example

    if (l_example == "michelle_imaging")
        page ("midir$doc/michelle_imaging_example.cl")
    else if (l_example == "trecs_imaging")
        page ("midir$doc/trecs_imaging_example.cl")
    else if (l_example == "michelle_longslit")
        page ("midir$doc/michelle_longslit_example.cl")
    else if (l_example == "midir_longslit")
        page ("midir$doc/midir_longslit_example.cl")
    else if (l_example == "polarimetry")
        page ("midir$doc/michelle_polarimetry_example.cl")
    else
        printf ("MIDIREXAMPLES - ERROR: Invalid example ID\n")

end
