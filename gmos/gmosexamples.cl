# Copyright(c) 2002-2012 Association of Universities for Research in Astronomy, Inc.

procedure gmosexamples(example)

char example {enum="imaging|longslit|standard_longslit|MOS|NS_longslit|NS_MOS|IFU_1slit|IFU_2slit|standard_IFU", prompt="Example to print"}

begin

    char l_example
    l_example = example

    if (l_example == "imaging")
        page("gmos$doc/gmos_imaging_example.cl")

    if (l_example == "longslit")
        page("gmos$doc/gmos_longslit_example.cl")

    if (l_example == "standard_longslit")
        page("gmos$doc/gmos_longslit_standard_example.cl")

    if (l_example == "MOS")
        page("gmos$doc/gmos_mos_example.cl")

    if (l_example == "NS_longslit")
        page("gmos$doc/gmos_ns_longslit_example.cl")

    if (l_example == "NS_MOS")
        page("gmos$doc/gmos_ns_mos_example.cl")

    if (l_example == "IFU_1slit")
        page("gmos$doc/gmos_ifu_1slit_example.cl")

    if (l_example == "IFU_2slit")
        page("gmos$doc/gmos_ifu_2slit_example.cl")

    if (l_example == "standard_IFU")
        page("gmos$doc/gmos_ifu_standard_example.cl")

end
