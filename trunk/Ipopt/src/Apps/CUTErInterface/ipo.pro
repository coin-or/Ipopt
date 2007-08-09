#!/bin/csh -f
#  ( Last modified on 23 Dec 2000 at 17:29:56 )
#
# ipopt: apply IPOPT to a problem and delete the executable after use.
#
#{version}
#

#{args}

#
#  N. Gould, D. Orban & Ph. Toint, November 7th, 2000
#

#{cmds}

#
# Environment check
#

envcheck
if( $status != 0 ) exit $status

#
#  define a short acronym for the package to which you wish to make an interface
#

setenv caller ipo
setenv PAC ipo

#
#  define the name of the subdirectory of $CUTER/common/src/pkg
#  in which the package lies
#

setenv PACKAGE ipo

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
    echo 'IPOPT is only available in double precision'
    echo 'rerun without -s'
    exit 1
  else if("$opt" == '-decode' ) then
    @ decode_problem = 1
  else if("$opt" == '-h' )then
    $MYCUTER/bin/helpmsg
    exit 0
  endif
  @ i++
end

#
#  define the system libraries needed by the package
#  using the format -lrary to include library.a
#

setenv SYSLIBS "`cat ${MYCUTER}/double/lib/ipo.liblist`"

#  define the name(s) of the object file(s) (files of the form *.o)
#  and/or shared-object libraries (files of the form lib*.so)
#  and/or libraries (files of the form lib*.a) for the IPOPT package.
#  Object files must lie in 
#    $MYCUTER/(precision)/bin  
#  while libraries and/or shared-object libraries must be in 
#    $MYCUTER/(precision)/lib 
#  where (precision) is either single (when the -s flag is present) or 
#  double (when there is no -s flag)
#

setenv PACKOBJ ""

#
#  define the name of the package specification file (if any)
#  (this possibly precision-dependent file must either lie in
#  the current directory or in $MYCUTER/common/src/pkg/$PACKAGE/ )
#

setenv SPECS ""

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted prior to 
#  decoding the problem file


#  decode the problem file

#{sifdecode}
if( $decode_problem == 1 ) then
  if( $?MYSIFDEC ) then
     $MYSIFDEC/bin/sifdecode $argv
  else
     echo " ${caller} : environment variable MYSIFDEC not set"
     echo "      Either SifDec is not installed or you"
     echo "      should properly set MYSIFDEC"
     exit 7
  endif
endif

if ( $status != 0 ) exit $status

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted prior to 
#  running the package


#  run the package, removing -decode option if present

if( $decode_problem == 1 ) then
  @ n = ${#argv} - 1
  @ i = 1
  set arguments = ''
  while( $i <= $n )
    if( "$argv[$i]" != '-decode' ) then
      set arguments = ( "$arguments" "$argv[$i]" )
    endif
    @ i++
  end
else
  set arguments = "$argv"
endif

$MYCUTER/bin/runpackage $arguments -n

if ( $status != 0 ) exit $status

#  this is where package-dependent commands (using the commands defined
#  in the current system.(your system) file) may be inserted to clean
#  up after running the package


