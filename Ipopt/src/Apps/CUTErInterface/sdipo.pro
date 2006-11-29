#!/bin/csh -f
# sdipopt: script to decode a sif file and then run IPOPT on the output
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
#
#{version}
#

#{args}

#
#  N. Gould, D. Orban & Ph. Toint, November 7th, 2000
#

#{cmds}

#
#  define a short acronym for the package to which you wish to make an interface
#

setenv caller sdipo
setenv PAC ipo

#
#  Check the arguments
#

set PRECISION = "double"
@ decode_problem = 0
@ last=$#argv
@ i=1

while ($i <= $last)
  set opt=$argv[$i]
  if("$opt" == '-s')then
    set PRECISION = "single"
  else if("$opt" == '-decode' ) then
    @ decode_problem = 1
  else if("$opt" == '-h' )then
    $MYCUTER/bin/helpmsg
    exit 0
  endif
  @ i++
end

if( $decode_problem == 0 ) then
  set arguments = ('-decode' "$argv[1-$last]")
else
  set arguments = ("$argv[1-$last]")
endif

# call main script

${PAC} ${arguments}
