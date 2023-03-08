# Copyright(c) 2002-2012 Association of Universities for Research in Astronomy, Inc.

procedure niriexamples(example)

char example {enum="imaging|longslit",prompt="Example to print"}

begin

    char l_example
    l_example = example

    if (l_example=="imaging") 
        page("niri$doc/niri_imaging_example.cl")

    if (l_example=="longslit")
        page("niri$doc/niri_longslit_example.cl")

end
