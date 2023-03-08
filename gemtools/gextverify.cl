# Copyright(c) 2002-2006 Association of Universities for Research in Astronomy, Inc.

procedure gextverify(ext)

# Version Feb 28, 2002  MT v1.3 release
#         Sept 20,2002  MT v1.4 release

string  ext              {prompt="Extension name"}
string  outext           {"", prompt="Output extension name"}
int     status       	 {0, prompt="Status"}

# status:
# 0 - extension name is not an empty string
# 1 - extension is an empty string

begin

string l_ext
int    n

l_ext = ext

status = 0
# Get rid of blank spaces
print (l_ext) | translit ("STDIN", " ", "", delete+, collapse-) | scan (l_ext)

if(l_ext=="" || l_ext=="[]") {
     status = 1
     print("ERROR GEXTVERIFY: Empty extension.")
} else outext = l_ext

end
