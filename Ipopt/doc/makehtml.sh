#!/bin/bash

# http://static.usenix.org/events/samples/conversion.html

# get rid of underlines in formulas
cat > l2hinit <<EOF
\$DVIPSOPT = '-E';
return 1
EOF

export HOME=`pwd`
export L2HINIT_NAME=l2hinit

cat $HOME/$L2HINIT_NAME

latex2html -nofootnode -antialias -local_icons -noaddress -noinfo -html_version 4.0 -split +2 documentation.tex

rm l2hinit
