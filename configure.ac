# Copyright (C) 2004, 2005 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AC_PREREQ(2.59)

AC_INIT([Ipopt],[3.0.0],[coin-ipopt@www-124.ibm.com])

AC_COPYRIGHT([
Copyright 2004, 2005 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package IPOPT which is distributed
under the Common Public License.])

AC_REVISION(3.0.0)

AC_CONFIG_SRCDIR(Common/IpDebug.hpp)

# Get the system type
AC_CANONICAL_BUILD

# Check if user wants to produce debugging code
AC_ARG_ENABLE([debug],
  AC_HELP_STRING([--enable-debug], [compile with debug options and runtime tests]),
  [case "${enableval}" in
    yes) debugcompile=true
         AC_DEFINE([IP_DEBUG],[1],
                   [If defined, debug sanity checks are performed during runtime])
         ;;
    no)  debugcompile=false ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
  esac],[debugcompile=false])

########################
# C++ Compiler Options #
########################

AC_LANG_PUSH(C++)
SAVE_CXXFLAGS="$CXXFLAGS"
case $build in
  *-cygwin) comps="g++" ;;
  *)        comps="xlC aCC CC g++ c++ icpc" ;;
esac
AC_PROG_CXX([$comps])
CXXFLAGS="$SAVE_CXXFLAGS"

warn_cxxflags=

if test -z "$CXXFLAGS"; then
  add_cxxflags=
  opt_cxxflags=
  dbg_cxxflags=

  if test "$GXX" = "yes"; then
    case "$CXX" in
      icpc | */icpc) ;;
      *)
# ToDo decide about unroll-loops
        opt_cxxflags="-O3"
#        case $build in
#          x86_64-* | i686-* | i586-* )
#            # Check for enabling special instructions
#            myflags=`grep flags /proc/cpuinfo | tail -1 | cut -d: -f2-`
#            have_sse=no
#            for i in $myflags; do
#              case $i in
#                sse3) CXXFLAGS="-mfpmath=sse -msse3"
#                     AC_MSG_CHECKING([whether flag \"$CXXFLAGS\" works])
#                     AC_TRY_LINK([#include <cmath>],[double s = sin(1.0);],
#	                [AC_MSG_RESULT([yes])
#                         opt_cxxflags="$opt_cxxflags $CXXFLAGS"],
#                        [AC_MSG_RESULT([no])])
#                     have_sse=yes;
#                     CXXFLAGS=
#                     ;;
#                sse2)
#                     CXXFLAGS="-mfpmath=sse -msse2"
#                     AC_MSG_CHECKING([whether flag \"$CXXFLAGS\" works])
#                     AC_TRY_LINK([#include <cmath>],[double s = sin(1.0);],
#	                [AC_MSG_RESULT([yes])
#                         if test $have_sse = no; then
#                           opt_cxxflags="$opt_cxxflags $CXXFLAGS"
#			 else
#			   opt_cxxflags="$opt_cxxflags -msse2"
#                         fi],
#                        [AC_MSG_RESULT([no])])
#                     have_sse=yes;
#                     CXXFLAGS=
#                     ;;
#                sse)
#                     CXXFLAGS="-mfpmath=sse -msse"
#                     AC_MSG_CHECKING([whether flag \"$CXXFLAGS\" works])
#                     AC_TRY_LINK([#include <cmath>],[double s = sin(1.0);],
#	                [AC_MSG_RESULT([yes])
#                         if test $have_sse = no; then
#                           opt_cxxflags="$opt_cxxflags $CXXFLAGS"
#			 else
#			   opt_cxxflags="$opt_cxxflags -msse"
#                         fi],
#                        [AC_MSG_RESULT([no])])
#                     have_sse=yes;
#                     CXXFLAGS=
#                     ;;
#                mmx) CXXFLAGS="-mmmx"
#                     AC_MSG_CHECKING([whether flag \"$CXXFLAGS\" works])
#                     AC_TRY_LINK([#include <cmath>],[double s = sin(1.0);],
#	                [AC_MSG_RESULT([yes])
#                         opt_cxxflags="$opt_cxxflags $CXXFLAGS"],
#                        [AC_MSG_RESULT([no])])
#                     CXXFLAGS= ;;
#                3dnow) CXXFLAGS="-m3dnow"
#                     AC_MSG_CHECKING([whether flag \"$CXXFLAGS\" works])
#                     AC_TRY_LINK([#include <cmath>],[double s = sin(1.0);],
#	                [AC_MSG_RESULT([yes])
#                         opt_cxxflags="$opt_cxxflags $CXXFLAGS"],
#                        [AC_MSG_RESULT([no])])
#                     CXXFLAGS= ;;
#              esac
#            done
#          ;;
#        esac
        add_cxxflags="-pipe"
        dbg_cxxflags="-g"
        warn_cxxflags="-pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Werror -Wconversion"
        case $build in
          *-cygwin*) AC_MSG_CHECKING([whether flag \" -mno-cygwin\" works])
                    CXXFLAGS="-mno-cygwin"
                    AC_TRY_LINK(,[int i=0; i++;],
                                    [AC_MSG_RESULT([yes]); add_cxxflags="-mno-cygwin $add_cxxflags"; use_mingw=yes],
                                    [AC_MSG_RESULT([no])])
                    CXXFLAGS= ;;
        esac ;;
    esac
  fi
  if test -z "$opt_cxxflags"; then
    case $build in
     *-cygwin*)
      case "$CXX" in
      cl | */cl)
             opt_cxxflags="-Ot1"
             add_cxxflags="-nologo -EHsc -GR"
             dbg_cxxflags="-Yd";;
      esac;;
     *-linux-*)
      case "$CXX" in
      icpc | */icpc)
#             AC_MSG_CHECKING([whether flag \"-fast\" works])
#               CXXFLAGS="-fast"
#               AC_TRY_RUN(,[int i=0; i++;],
#                               [AC_MSG_RESULT([yes]); opt_cxxflags="-fast"],
#                               [AC_MSG_RESULT([no]); opt_cxxflags="-O3 -ip -static"])
#               CXXFLAGS=
             opt_cxxflags="-O3 -ip"
             add_cxxflags="-mp -pc64"
             dbg_cxxflags="-g"
             # Check if -i_dynamic is necessary (for new glibc library)
             AC_MSG_CHECKING([whether flag \"-i_dynamic\" is needed])
             CXXFLAGS=
             AC_TRY_LINK(,[int i=0; i++;],
                         [AC_MSG_RESULT([apparently not])],
                         [AC_MSG_RESULT([seems so]);
                               add_cxxflags="-i_dynamic $add_cxxflags"])
             ;;
      esac;;
    *-ibm-*)
      case "$CXX" in
      xlC* | */xlC* | mpxlC* | */mpxlC*)
             opt_cxxflags="-O3 -qarch=auto -qcache=auto -qhot -qtune=auto -qmaxmem=-1"
             add_cxxflags="-bmaxdata:0x80000000 -qrtti=dyna"
             dbg_cxxflags="-g";;
      esac;;
    *-hp-*)
      case "$CXX" in
      aCC | */aCC )
             opt_cxxflags="-O"
             add_cxxflags="-AA"
             dbg_cxxflags="-g";;
      esac;;
    *-sun-*)
             opt_cxxflags="-O5 -xtarget=native"
             dbg_cxxflags="-g";;
    esac
  fi

  if test "$ac_cv_prog_cxx_g" = yes && test -z "$dbg_cxxflags" ; then
    dbg_cxxflags="-g"
  fi

  if test "$debugcompile" = "true"; then
    CXXFLAGS="$dbg_cxxflags $add_cxxflags"
  else
    if test -z "$opt_cxxflags"; then
# Try if -O option works if nothing else is set
      CXXFLAGS="-O"
      AC_TRY_LINK([],[int i=0; i++;],[opt_cxxflags="-O"],[])
    fi
    CXXFLAGS="$opt_cxxflags $add_cxxflags"
  fi
fi

# Try if CXXFLAGS works
AC_MSG_CHECKING([whether CXXFLAGS=\"$CXXFLAGS\" works])
AC_TRY_LINK([],[int i=0; i++;],[],[CXXFLAGS=])
if test -z "$CXXFLAGS"; then
  AC_MSG_RESULT([no])
  AC_MSG_WARN([This value for CXXFLAGS does not work.  I will now just try '-O', but you might want to set CXXFLAGS manually.])
  CXXFLAGS='-O'
  AC_TRY_LINK([],[int i=0; i++;],[],[CXXFLAGS=])
  if test -z "$CXXFLAGS"; then
    AC_MSG_RESULT([no])
    AC_MSG_WARN([This value for CXXFLAGS does not work.  I will continue with empty CXXFLAGS, but you might want to set CXXFLAGS manually.])
  fi
else
  AC_MSG_RESULT([yes])
fi

# C++ link libraries #
AC_ARG_VAR(CXXLIBS,[Libraries necessary for linking C++ code with Fortran compiler])
if test -z "$CXXLIBS"; then
  if test "$GXX" = "yes"; then
    case "$CXX" in
      icpc | */icpc)
        CXXLIBS=""
        ;;
      *)
        CXXLIBS="-lstdc++ -lm -lgcc"
        ;;
    esac
  else
    case $build in
     *-linux-*)
      case "$CXX" in
      icpc | */icpc)
        CXXLIBS=""
             ;;
      esac;;
    *-ibm-*)
      CXXLIBS="-lC -lc"
      ;;
    *-hp-*)
      CXXLIBS="-L/opt/aCC/lib -l++ -lstd_v2 -lCsup_v2 -lm -lcl -lc"
      ;;
    *-sun-*)
      CXXLIBS="-lCstd -lCrun"
    esac
  fi
fi
if test -z "$CXXLIBS"; then
  AC_MSG_WARN([Could not automatically determine CXXLIBS (C++ link libraries; necessary if main program is in Fortran of C).])
else
  AC_MSG_NOTICE([Assuming that CXXLIBS is \"$CXXLIBS\".])
fi
AC_LANG_POP(C++)

###############################################
# Set variables for the archiver              #
# (for AIX in 64bit mode also on SGI, SUN...) #
###############################################
if test -z "$AR"; then
  if test ! "$GXX" = "yes"; then
    case $build in
      *-sun-*)
        AR="CC -xar"
        ;;
      *-sgi-*)
        AR="CC -ar"
        ;;
    esac
  fi
  if test -z "$AR"; then
    AR="ar"
  fi
fi
AC_ARG_VAR(AR,[Command for handling archives (default is `ar')])
AC_MSG_NOTICE([Using AR=\"$AR\".])

if test -z "$ARFLAGS"; then
  if test ! "$GXX" = "yes"; then
    case $build in
      *-sun-*)
        ARFLAGS="-o"
        ;;
      *-sgi-*)
        ARFLAGS="-o"
        ;;
    esac
  fi
  if test -z "$ARFLAGS"; then
    ARFLAGS="crv"
  fi
fi
AC_ARG_VAR(ARFLAGS,[Flags for adding file for archive command (default is `crv')])
AC_MSG_NOTICE([Using ARFLAGS=\"$ARFLAGS\".])

if test -z "$AR_X"; then
  AR_X="ar x"
fi
AC_ARG_VAR(AR_X,[Command for extracing an archive (default is `ar x')])
AC_MSG_NOTICE([Using AR_X=\"$AR_X\".])

######################
# C Compiler Options #
######################

AC_LANG_PUSH(C)
SAVE_CFLAGS="$CFLAGS"
case $build in
  *-cygwin) comps="gcc" ;;
  *)        comps="xlc cc gcc icc" ;;
esac
AC_PROG_CC([$comps])
CFLAGS="$SAVE_CFLAGS"

warn_cflags=

if test -z "$CFLAGS"; then
  add_cflags=
  opt_cflags=
  dbg_cflags=

  if test "$GCC" = "yes"; then
    case "$CC" in
      icc | */icc) ;;
      *)
        opt_cflags="-O3"
#        case $build in
#          x86_64-* | i686-* | i586-* )
#            # Check for enabling special instructions
#            myflags=`grep flags /proc/cpuinfo | tail -1 | cut -d: -f2-`
#            have_sse=no
#            SAVE_LIBS="$LIBS"
#            LIBS="$LIBS -lm"
#            for i in $myflags; do
#              case $i in
#                sse3) CFLAGS="-mfpmath=sse -msse3"
#                     AC_MSG_CHECKING([whether flag \"$CFLAGS\" works])
#                     AC_TRY_LINK([#include <math.h>],[double s = sin(1.0);],
#	                [AC_MSG_RESULT([yes])
#                         opt_cflags="$opt_cflags $CFLAGS"],
#                        [AC_MSG_RESULT([no])])
#                     have_sse=yes;
#                     CFLAGS=
#                     ;;
#                sse2)
#                     CFLAGS="-mfpmath=sse -msse2"
#                     AC_MSG_CHECKING([whether flag \"$CFLAGS\" works])
#                     AC_TRY_LINK([#include <math.h>],[double s = sin(1.0);],
#	                [AC_MSG_RESULT([yes])
#                         if test $have_sse = no; then
#                           opt_cflags="$opt_cflags $CFLAGS"
#			 else
#			   opt_cflags="$opt_cflags -msse2"
#                         fi],
#                        [AC_MSG_RESULT([no])])
#                     have_sse=yes;
#                     CFLAGS=
#                     ;;
#                sse)
#                     CFLAGS="-mfpmath=sse -msse"
#                     AC_MSG_CHECKING([whether flag \"$CFLAGS\" works])
#                     AC_TRY_LINK([#include <math.h>],[double s = sin(1.0);],
#	                [AC_MSG_RESULT([yes])
#                         if test $have_sse = no; then
#                           opt_cflags="$opt_cflags $CFLAGS"
#			 else
#			   opt_cflags="$opt_cflags -msse"
#                         fi],
#                        [AC_MSG_RESULT([no])])
#                     have_sse=yes;
#                     CFLAGS=
#                     ;;
#                mmx) CFLAGS="-mmmx"
#                     AC_MSG_CHECKING([whether flag \"$CFLAGS\" works])
#                     AC_TRY_LINK([#include <math.h>],[double s = sin(1.0);],
#	                [AC_MSG_RESULT([yes])
#                         opt_cflags="$opt_cflags $CFLAGS"],
#                        [AC_MSG_RESULT([no])])
#                     CFLAGS= ;;
#                3dnow) CFLAGS="-m3dnow"
#                     AC_MSG_CHECKING([whether flag \"$CFLAGS\" works])
#                     AC_TRY_LINK([#include <math.h>],[double s = sin(1.0);],
#	                [AC_MSG_RESULT([yes])
#                         opt_cflags="$opt_cflags $CFLAGS"],
#                        [AC_MSG_RESULT([no])])
#                     CFLAGS= ;;
#              esac
#            done
#            LIBS="$SAVE_LIBS"
#          ;;
#        esac
        add_cflags="-pipe"
        dbg_cflags="-g"
	warn_cflags="-pedantic-errors -Wimplicit -Wparentheses -Wsequence-point -Wreturn-type -Wcast-qual -Wall -Werror "
        case $build in
          *-cygwin*) AC_MSG_CHECKING([whether flag \"-mno-cygwin\" works])
                    CFLAGS="-mno-cygwin"
                    AC_TRY_LINK(,[int i=0; i++;],
                                    [AC_MSG_RESULT([yes]); add_cflags="-mno-cygwin $add_cflags"; use_mingw=yes],
                                    [AC_MSG_RESULT([no])])
                    CFLAGS= ;;
        esac ;;
    esac
  fi
  if test -z "$opt_cflags"; then
    case $build in
     *-cygwin*)
      case "$CC" in
      cl | */cl)
             opt_cflags="-Ot1"
             add_cflags="-nologo"
             dbg_cflags="-Yd";;
      esac;;
     *-linux-*)
      case "$CC" in
      icc | */icc)
#             AC_MSG_CHECKING([whether flag \"-fast\" works])
#               CFLAGS="-fast"
#               AC_TRY_RUN(,[int i=0; i++;],
#                               [AC_MSG_RESULT([yes]); opt_cflags="-fast"],
#                               [AC_MSG_RESULT([no]); opt_cflags="-O3 -ipo -static"])
#               CFLAGS=
             opt_cflags="-O3 -ip"
             add_cflags="-mp -pc64"
             dbg_cflags="-g"
             # Check if -i_dynamic is necessary (for new glibc library)
             AC_MSG_CHECKING([whether flag \"-i_dynamic\" is needed])
             CFLAGS=
             AC_TRY_LINK(,[int i=0; i++;],
                         [AC_MSG_RESULT([apparently not])],
                         [AC_MSG_RESULT([seems so]);
                               add_cflags="-i_dynamic $add_cflags"])
             ;;
      pgcc | */pgcc)
             opt_cflags="-fast"
             add_cflags="-Kieee -Mstabs -pc 64"
             dbg_cflags="-g" ;;
      esac;;
    *-ibm-*)
      case "$CC" in
      xlc* | */xlc* | mpxlc* | */mpxlc*)
             opt_cflags="-O3 -qarch=auto -qcache=auto -qhot -qtune=auto -qmaxmem=-1"
             add_cflags="-bmaxdata:0x80000000"
             dbg_cflags="-g";;
      esac;;
    *-hp-*)
             opt_cflags="-O"
             add_cflags="-Ae"
             dbg_cflags="-g";;
    *-sun-*)
             opt_cflags="-xO5 -xtarget=native"
             dbg_cflags="-g";;
    *-sgi-*)
             opt_cflags="-O -OPT:Olimit=0"
             dbg_cflags="-g";;
    esac
  fi

  if test "$ac_cv_prog_cc_g" = yes && test -z "$dbg_cflags" ; then
    dbg_cflags="-g"
  fi

  if test "$debugcompile" = "true"; then
    CFLAGS="$dbg_cflags $add_cflags"
  else
    if test -z "$opt_cflags"; then
# Try if -O option works if nothing else is set
      CFLAGS="-O"
      AC_TRY_LINK([],[int i=0; i++;],[opt_cflags="-O"],[])
    fi
    CFLAGS="$opt_cflags $add_cflags"
  fi
fi

# Try if CFLAGS works
AC_MSG_CHECKING([whether CFLAGS=\"$CFLAGS\" works])
AC_TRY_LINK([],[int i=0; i++;],[],[CFLAGS=])
if test -z "$CFLAGS"; then
  AC_MSG_RESULT([no])
  AC_MSG_WARN([This value for CFLAGS does not work.  I will now just try '-O', but you might want to set CFLAGS manually.])
  CFLAGS='-O'
  AC_TRY_LINK([],[int i=0; i++;],[],[CFLAGS=])
  if test -z "$CFLAGS"; then
    AC_MSG_RESULT([no])
    AC_MSG_WARN([This value for CFLAGS does not work.  I will continue with empty CFLAGS, but you might want to set CFLAGS manually.])
  fi
else
  AC_MSG_RESULT([yes])
fi
AC_LANG_POP(C)

###################################
#   Fortran 77 compiler options   #
###################################

AC_LANG_PUSH(Fortran 77)
SAVE_FFLAGS="$FFLAGS"
case $build in
  *-cygwin) comps="gfortran g77" ;;
  *)        comps="xlf fort77 gfortran f77 g77 ifort ifc" ;;
esac
AC_PROG_F77($comps)
FFLAGS="$SAVE_FFLAGS"

if test -z "$FFLAGS"; then
  add_fflags=
  opt_fflags=
  dbg_fflags=

  if test "$G77" = "yes"; then
    opt_fflags="-O3"
#        case $build in
#          x86_64-* | i686-* | i586-* )
#            # Check for enabling special instructions
#            myflags=`grep flags /proc/cpuinfo | tail -1 | cut -d: -f2-`
#            have_sse=no
#            for i in $myflags; do
#              case $i in
#                sse3) FFLAGS="-mfpmath=sse -msse3"
#                     AC_MSG_CHECKING([whether flag \"$FFLAGS\" works])
#                     AC_TRY_LINK(,[      double precision x
#      x = sin(1.d0)],
#	                [AC_MSG_RESULT([yes])
#                         opt_fflags="$opt_fflags $FFLAGS"],
#                        [AC_MSG_RESULT([no])])
#                     have_sse=yes;
#                     FFLAGS=
#                     ;;
#                sse2)
#                     FFLAGS="-mfpmath=sse -msse2"
#                     AC_MSG_CHECKING([whether flag \"$FFLAGS\" works])
#                     AC_TRY_LINK(,[      double precision x
#      x = sin(1.d0)],
#	                [AC_MSG_RESULT([yes])
#                         if test $have_sse = no; then
#                           opt_fflags="$opt_fflags $FFLAGS"
#			 else
#			   opt_fflags="$opt_fflags -msse2"
#                         fi],
#                        [AC_MSG_RESULT([no])])
#                     have_sse=yes;
#                     FFLAGS=
#                     ;;
#                sse)
#                     FFLAGS="-mfpmath=sse -msse"
#                     AC_MSG_CHECKING([whether flag \"$FFLAGS\" works])
#                     AC_TRY_LINK(,[      double precision x
#      x = sin(1.d0)],
#	                [AC_MSG_RESULT([yes])
#                         if test $have_sse = no; then
#                           opt_fflags="$opt_fflags $FFLAGS"
#			 else
#			   opt_fflags="$opt_fflags -msse"
#                         fi],
#                        [AC_MSG_RESULT([no])])
#                     have_sse=yes;
#                     FFLAGS=
#                     ;;
#                mmx) FFLAGS="-mmmx"
#                     AC_MSG_CHECKING([whether flag \"$FFLAGS\" works])
#                     AC_TRY_LINK(,[      double precision x
#      x = sin(1.d0)],
#	                [AC_MSG_RESULT([yes])
#                         opt_fflags="$opt_fflags $FFLAGS"],
#                        [AC_MSG_RESULT([no])])
#                     FFLAGS= ;;
#                3dnow) FFLAGS="-m3dnow"
#                     AC_MSG_CHECKING([whether flag \"$FFLAGS\" works])
#                     AC_TRY_LINK(,[      double precision x
#      x = sin(1.d0)],
#	                [AC_MSG_RESULT([yes])
#                         opt_fflags="$opt_fflags $FFLAGS"],
#                        [AC_MSG_RESULT([no])])
#                     FFLAGS= ;;
#              esac
#            done
#          ;;
#        esac

#    case $build in
#      i686-* )
#        FFLAGS="-march=i686 -mfpmath=sse -msse"
#        AC_MSG_CHECKING([whether flags \"$FFLAGS\" work])
#        AC_TRY_LINK(,[
#      double precision x
#      x = sin(1.d0)],
#                       [AC_MSG_RESULT([yes]); opt_fflags="$opt_fflags $FFLAGS"],
#                       [AC_MSG_RESULT([no])])
#        FFLAGS=
#      ;;
#    esac
    add_fflags="-pipe"
    dbg_fflags="-g"
    case $build in
      *-cygwin*) AC_MSG_CHECKING([whether flag \"-mno-cygwin\" works])
                FFLAGS="-mno-cygwin"
                AC_TRY_LINK(,[      write(*,*) 'Hello world'],
                                [AC_MSG_RESULT([yes]); add_fflags="-mno-cygwin $add_fflags"],
                                [AC_MSG_RESULT([no])])
                FFLAGS= ;;
    esac
  else
    case $build in
    *-cygwin*)
      case $F77 in
      ifort | */ifort)
             opt_fflags="-O3"
             add_fflags="-nologo -Qprec -fpp"
             dbg_fflags="-debug" ;;
      esac;;
    *-linux-*)
      case $F77 in
      ifc | */ifc | ifort | */ifort)
#             AC_MSG_CHECKING([whether flag \"-fast\" works])
#               FFLAGS="-fast"
#               AC_TRY_RUN(,[      write(*,*) 'Hello world'],
#                               [AC_MSG_RESULT([yes]); opt_fflags="-fast"],
#                               [AC_MSG_RESULT([no]); opt_fflags="-O3 -ipo -static"])
#               FFLAGS= 
             opt_fflags="-O3 -ip"
             add_fflags="-mp -pc64 -cm -w90 -w95"
             dbg_fflags="-g -CA -CB -CS -CU"
             # Check if -i_dynamic is necessary (for new glibc library)
             AC_MSG_CHECKING([whether flag \"-i_dynamic\" is needed])
             FFLAGS=
             AC_TRY_LINK(,[      write(*,*) 'Hello world'],
                         [AC_MSG_RESULT([apparently not])],
                         [AC_MSG_RESULT([seems so]);
                               add_fflags="-i_dynamic $add_fflags"])
             ;;
      pgf77 | */pgf77)
             opt_fflags="-fast"
             add_fflags="-Kieee -Mstabs -pc 64"
             dbg_fflags="-g" ;;
      esac;;
    *-ibm-*)
      case $F77 in
      xlf* | */xlf* | mpxlf* | */mpxlf* )
             opt_fflags="-O3 -qarch=auto -qcache=auto -qhot -qtune=auto -qmaxmem=-1"
             add_fflags="-bmaxdata:0x80000000"
             dbg_fflags="-g -C";;
      esac;;
    *-hp-*)
             opt_fflags="+O3"
             add_fflags="+U77"
             dbg_fflags="-C -g";;
    *-sun-*)
             opt_fflags="-O5 -xtarget=native"
             dbg_fflags="-g";;
    *-sgi-*)
             opt_fflags="-O5 -OPT:Olimit=0"
             dbg_fflags="-g";;
    esac
  fi

  if test "$ac_cv_prog_f77_g" = yes && test -z "$dbg_fflags" ; then
    dbg_fflags="-g"
  fi

  if test "$debugcompile" = true; then
    FFLAGS="$dbg_fflags $add_fflags"
  else
    if test -z "$opt_fflags"; then
# Try if -O option works if nothing else is set
      AC_MSG_CHECKING([whether FFLAGS=\"-O\" works])
      FFLAGS=-O
      AC_TRY_LINK(,[      integer i],[AC_MSG_RESULT([yes]); opt_fflags="-O"],[AC_MSG_RESULT([no])])
    fi
    FFLAGS="$opt_fflags $add_fflags"
  fi
fi

# Try if FFLAGS works
AC_MSG_CHECKING([whether FFLAGS=\"$FFLAGS\" works])
AC_TRY_LINK([],[      integer i],[],[FFLAGS=])
if test -z "$FFLAGS"; then
  AC_MSG_RESULT([no])
  AC_MSG_WARN([This value for FFLAGS does not work.  I will now just try '-O', but you might want to set FFLAGS manually.])
  FFLAGS='-O'
  AC_TRY_LINK([],[      integer i],[],[FFLAGS=])
  if test -z "$FFLAGS"; then
    AC_MSG_RESULT([no])
    AC_MSG_WARN([This value for FFLAGS does not work.  I will continue with empty FFLAGS, but you might want to set FFLAGS manually.])
  fi
else
  AC_MSG_RESULT([yes])
fi

AC_F77_WRAPPERS

case $build in
# The following is a fix to define FLIBS for ifort on Windows
   *-cygwin*)
     case $F77 in
       ifort | */ifort)
           FLIBS="libifcorert.lib libcmt.lib /link /NODEFAULTLIB:libc.lib $LIBS";;
     esac;;
   *-hp-*)
       FLIBS="$FLIBS -lm";;
   *-ibm-*)
       FLIBS=`echo $FLIBS | sed 's/-lc)/-lc/g'` ;;
esac

AC_LANG_POP(Fortran 77)

AC_PROG_RANLIB

AC_PREFIX_DEFAULT([`pwd`])

AC_PROG_LN_S 

# Stuff for automake
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

############################################################################
############################################################################
#                   Stuff that we need for C++ programs                    #
############################################################################
############################################################################

AC_LANG_PUSH(C++)

################
# Header Files #
################

AC_DEFUN([MY_CHECK_CXX_HEADER],
         [AC_CHECK_HEADERS([c$1])
          if test ! "$ac_cv_header_c$1" = "yes"; then
            AC_CHECK_HEADERS([$1.h])
          fi])

MY_CHECK_CXX_HEADER(math)
MY_CHECK_CXX_HEADER(ctype)
MY_CHECK_CXX_HEADER(assert)
MY_CHECK_CXX_HEADER(stdio)
MY_CHECK_CXX_HEADER(stdarg)
MY_CHECK_CXX_HEADER(float)

###################
# Function finite #
###################

# First check for a few headers

case $build in
  *-cygwin*) ;;
  *) MY_CHECK_CXX_HEADER([ieeefp]) ;;
esac

# now try to guess the correct name of the finite function

MY_C_FINITE=
AC_CHECK_DECL([finite],[MY_C_FINITE=finite],,[
#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# endif
#endif
#ifdef HAVE_CFLOAT
# include <cfloat>
#else
# ifdef HAVE_FLOAT_H
#  include <float.h>
# endif
#endif
#ifdef HAVE_CIEEEFP
# include <cieeefp>
#else
# ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
# endif
#endif])
if test -z "$MY_C_FINITE"; then
  AC_CHECK_DECL([_finite],[MY_C_FINITE=_finite],,[
#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# endif
#endif
#ifdef HAVE_CFLOAT
# include <cfloat>
#else
# ifdef HAVE_FLOAT_H
#  include <float.h>
# endif
#endif
#ifdef HAVE_CIEEEFP
# include <cieeefp>
#else
# ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
# endif
#endif])
  if test -z "$MY_C_FINITE"; then
    AC_CHECK_DECL([isfinite],[MY_C_FINITE=isfinite],,[
#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# endif
#endif
#ifdef HAVE_CFLOAT
# include <cfloat>
#else
# ifdef HAVE_FLOAT_H
#  include <float.h>
# endif
#endif
#ifdef HAVE_CIEEEFP
# include <cieeefp>
#else
# ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
# endif
#endif])
  fi
fi
if test -z "$MY_C_FINITE"; then
  AC_MSG_WARN(Cannot find C-function for checking NaN and Inf.  Check will be disabled)
else
  AC_DEFINE_UNQUOTED([MY_C_FINITE],[$MY_C_FINITE],
                     [Define to be the name of C-function for NaNInf check])
fi

###########
# va_copy #
###########

AC_CHECK_DECL([va_copy],[AC_DEFINE([HAVE_VA_COPY],[1],
                                   [Define to 1 if va_copy is avaliable])],,[
#ifdef HAVE_CSTDARG
# include <cstdarg>
#else
# ifdef HAVE_STDARG_H
#  include <stdarg.h>
# endif
#endif])

##########################################################################

# Additional libraries and packages

IPADDLIBS=
AC_SUBST(IPADDLIBS)

##########
#  BLAS  #
##########

ipblaslibdir="\$(top_builddir)/Extern/blas"

AC_ARG_WITH([blas],
            AC_HELP_STRING([--with-blas],
                           [specify BLAS library (or BUILD for compilation)]),
            [use_blas=$withval], [use_blas=])
# Check if user supplied option makes sense
if test ! -z "$use_blas"; then
  if test "$use_blas" = "BUILD"; then
    AC_CHECK_FILE([$srcdir/Extern/blas/dgemm.f],[],
                  [AC_MSG_ERROR([option \"BUILD\" specified for blas, but sources are missing])])
  else
    AC_MSG_CHECKING([whether user supplied BLASLIB=\"$use_blas\" works])
    LIBS="$use_blas $LIBS"
    IPADDLIBS="$use_blas $IPADDLIBS"
    # Probably the blas library requires the Fortran libraries, so add
    # them for this test.  They are added at the end anyway
    SAVE_LIBS="$LIBS"
    LIBS="$LIBS $FLIBS"
    AC_F77_FUNC(daxpy,cfunc)
    AC_TRY_LINK([extern "C" {void $cfunc();}],[$cfunc()],
                [AC_MSG_RESULT([yes])],
                [AC_MSG_RESULT([no])
                 AC_MSG_ERROR([user supplied BLAS library \"$use_blas\" does not work])])
    LIBS="$SAVE_LIBS"
  fi
else
# Try to autodetect the library for blas based on build system
  case $build in
    *-sgi-*) AC_F77_FUNC(daxpy,cfunc)
             AC_CHECK_LIB([complib.sgimath],[$cfunc],
                          [use_blas=-lcomplib.sgimath;
                           LIBS="-lcomplib.sgimath $LIBS";
                           IPADDLIBS="-lcomplib.sgimath $IPADDLIBS"])
             ;;
    *-sun-*) SAVE_LIBS="$LIBS"
             AC_MSG_CHECKING([whether -xlic_lib=sunperf has BLAS])
             LIBS="-xlic_lib=sunperf $LIBS"
             AC_F77_FUNC(daxpy,cfunc)
             AC_TRY_LINK([extern "C" {void $cfunc();}],[$cfunc()],
                         [AC_MSG_RESULT([yes])
			  IPADDLIBS="-xlic_lib=sunperf $IPADDLIBS"
			  use_blas=-xlic_lib=sunperf],
                         [AC_MSG_RESULT([no])
			  LIBS="$SAVE_LIBS"])
             ;;
  esac
  if test -z "$use_blas"; then
    # Probably the blas library requires the Fortran libraries, so add
    # them for this test.  They are added at the end anyway
    AC_F77_FUNC(daxpy,cdaxpy)
    SAVE_LIBS="$LIBS"
    LIBS="$LIBS $FLIBS"
    AC_CHECK_LIB(blas,$cdaxpy,[use_blas=-lblas;
                               LIBS="-lblas $SAVE_LIBS";
                               IPADDLIBS="-lblas $IPADDLIBS"],
                              [LIBS="$SAVE_LIBS"])
  fi
  if test -z "$use_blas"; then
    AC_CHECK_FILE([$srcdir/Extern/blas/dgemm.f],[use_blas=BUILD],
                  [AC_MSG_ERROR([cannot find BLAS, and sources are missing in $srcdir/Extern/blas/])])
  fi
fi

# ##########
# # LAPACK #
# ##########
# 
# iplapacklibdir="\$(top_builddir)/Extern/lapack"
# lapackdir=Extern/lapack/LAPACK/SRC
# lapacksrcdir=$srcdir/$lapackdir
# lapackobjdir=Extern/lapack
# 
# AC_ARG_WITH([lapack],
#             AC_HELP_STRING([--with-lapack],
#                            [specify LAPACK library (or BUILD for compilation)]),
#             [use_lapack=$withval], [use_lapack=])
# # Check if user supplied option makes sense
# if test ! -z "$use_lapack"; then
#   if test "$use_lapack" = "BUILD"; then
#     AC_CHECK_FILE([$lapacksrcdir/dlarf.f],[],
#                   [AC_MSG_ERROR([option \"BUILD\" specified for lapack, but sources are missing])])
#   else
#     AC_MSG_CHECKING([whether user supplied LAPACKLIB=\"$use_lapack\" works])
#     LIBS="$use_lapack $LIBS"
#     IPADDLIBS="$use_lapack $IPADDLIBS"
#     AC_TRY_LINK([],[      call dlasr],[AC_MSG_RESULT([yes])],
#                 [AC_MSG_RESULT([no]); AC_MSG_ERROR([user supplied LAPACK library \"$use_lapack\" does not work])])
#   fi
# else
# # First try to see if LAPACK is already available with BLAS library
#   AC_MSG_CHECKING([whether LAPACK is already available with BLAS library])
#   AC_TRY_LINK([],[      call dgetrf()],
#             [AC_MSG_RESULT([yes]); use_lapack=ok],
#             [AC_MSG_RESULT([no])])
#   if test -z "$use_lapack"; then
#     # Try to autodetect the library for lapack based on build system
#     case $build in
#       *-sgi-*)   AC_CHECK_LIB([complib.sgimath],[dgetrf],
#                             [use_lapack=-lcomplib.sgimath;
#                              LIBS="-lcomplib.sgimath $LIBS";
#                              IPADDLIBS="-lcomplib.sgimath $IPADDLIBS"]);;
#       *-sun-*) SAVE_LDFLAGS="$LDFLAGS"
#                AC_MSG_CHECKING([whether -xlic_lib=sunperf has LAPACK])
#                LDFLAGS="-xlic_lib=sunperf $LDFLAGS"
#                AC_TRY_LINK([],[      call dgetrf],
#                            [AC_MSG_RESULT([yes]); use_lapack=-xlic_lib=sunperf],
#                            [AC_MSG_RESULT([no]); LDFLAGS="$SAVE_LDFLAGS"]) ;;
#     esac
#   fi
#   if test -z "$use_lapack"; then
#     AC_CHECK_LIB(lapack,dgetrf,[use_lapack=-llapack;
#                                 LIBS="-llapack $LIBS";
#                                 IPADDLIBS="-llapack $IPADDLIBS"],[])
#   fi
#   if test -z "$use_lapack"; then
#     AC_CHECK_FILE([$lapacksrcdir/dlarf.f],[use_lapack=BUILD],
#                   [AC_MSG_ERROR([cannot find LAPACK, and sources are missing in $lapacksrcdir])])
#   fi
# fi

#######################
# HARWELL SUBROUTINES #
#######################

hslsrcdir="$srcdir/Extern/HSL"
hsllibdir="\$(top_builddir)/Extern/HSL"

AC_DEFUN([MY_CHECK_HSL],
         [AC_MSG_CHECKING([for $1 in $use_hsl])
          use_$2=
          AC_F77_FUNC($1,cfunc)
          AC_TRY_LINK([extern "C" {void $cfunc();}],[$cfunc()],
                      [AC_DEFINE(HAVE_$2,1,[Define to 1 if $2 is available])
		       AC_MSG_RESULT([yes]); use_$2=yes],
                      [AC_MSG_RESULT([no]); use_$2=no])])

AC_DEFUN([MY_CHECK_HSLFILE],
         [AC_CHECK_FILE([$hslsrcdir/$1.f],
                        [AC_DEFINE(HAVE_$2,1,[Define to 1 if $2 is available]) use_$2=BUILD],
                        [use_$2=no])])

AC_MSG_NOTICE([Checking for Harwell subroutine library])

AC_ARG_WITH([hsl],
            AC_HELP_STRING([--with-hsl],
                           [specify Harwell library (or BUILD for compilation)]),
            [use_hsl=$withval], [use_hsl=BUILD])

if test "$use_hsl" = BUILD; then
# Check which files are available
  MY_CHECK_HSLFILE(ma27ad,MA27)
  if test "$use_MA27" = BUILD; then
    AC_MSG_CHECKING([whether MA27 is threadsafe])
    if $EGREP MA27UD $hslsrcdir/ma27ad.f >&5 ; then
      AC_MSG_RESULT(yes)
    else
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([Please obtain the threadsafe (newer) version of MA27])
    fi
  fi
  MY_CHECK_HSLFILE(mc19ad,MC19)
else
  SAVE_LIBS="$LIBS"
  LIBS="$use_hsl $LIBS $FLIBS"
  IPADDLIBS="$use_hsl $IPADDLIBS"
  MY_CHECK_HSL([ma27ad],[MA27])
  if test "$use_MA27" = "yes"; then
    AC_MSG_CHECKING([whether MA27 is threadsafe])
    AC_F77_FUNC(ma27ud, cfunc)
    AC_TRY_LINK([extern "C" {void $cfunc();}],[$cfunc],
                [AC_MSG_RESULT(yes)],
                [AC_MSG_RESULT(no)
		 AC_MSG_ERROR([Please obtain the threadsafe (newer) version of MA27])])
  fi
  MY_CHECK_HSL([mc19ad],[MC19])
  LIBS="$use_hsl $SAVE_LIBS"
fi

###########
# PARDISO #
###########

AC_ARG_WITH([pardiso],
            AC_HELP_STRING([--with-pardiso],
                           [specify Pardiso library (unsupported)]),
            [use_pardiso=$withval], [use_pardiso=no])

if test ! "$use_pardiso" = "no"; then
  LIBS="$use_pardiso $LIBS"
  IPADDLIBS="$use_pardiso $IPADDLIBS"

# ToDo put this back in
#  AC_MSG_CHECKING([whether user-supplied Pardiso library \"$use_pardiso\" works])
#  AC_TRY_LINK([],[      call PARDISO()],
#              [AC_MSG_RESULT(yes)],
#              [AC_MSG_ERROR([Pardiso library $use_pardiso does not seem to work])])

  AC_DEFINE(HAVE_PARDISO,1,
            [Define to 1 if Pardiso is available (unsupported)])
fi

#########
# TAUCS #
#########

AC_ARG_WITH([taucs],
            AC_HELP_STRING([--with-taucs],
                           [specify TAUCS library and dependencies (unsupported)]),
            [use_taucs=$withval], [use_taucs=no])

if test ! "$use_taucs" = "no"; then
  LIBS="$use_taucs $LIBS"
  IPADDLIBS="$use_taucs $IPADDLIBS"

  AC_DEFINE(HAVE_TAUCS,1,
            [Define to 1 if TAUCS is available (unsupported)])

  AC_ARG_VAR(TAUCS_INCFLAGS,[Flags for including taucs header files])
  CXXFLAGS="$CXXFLAGS $TAUCS_INCFLAGS"
fi

# Check if at least one linear solver was found

AC_ARG_ENABLE([distconf],,
  [case "${enableval}" in
    yes) distconf=true ;;
    no)  distconf=false ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-distconf) ;;
  esac],[distconf=false])

if test "$use_MA27" = "no" -a "$use_pardiso" = "no" -a "$use_taucs" = no; then
  if test "$distconf" = false; then
    AC_MSG_ERROR("No linear solver (such as MA27) available")
  fi
fi

# For the remaining tests we need to include FLIBS into LIBS, so that
# the C compiler can link programs
LIBS="$LIBS $FLIBS"

#######################
# AMPL SOLVER LIBRARY #
#######################

asldir=Extern/ASL/solvers
aslsrcdir=$srcdir/$asldir

AC_ARG_WITH([asldir],
            AC_HELP_STRING([--with-asldir],
                           [specify absolute path to AMPL solver directory (or BUILD for compilation, or "no" for disabling AMPL)]),
            [use_asldir=$withval], [use_asldir=])

if test "$use_asldir" = BUILD; then
  AC_CHECK_FILE([$aslsrcdir/asl.h],[],
                [AC_MSG_ERROR([option \"BUILD\" specified for asldir, but sources are missing])])
elif test -z "$use_asldir"; then
 # try to find sources - if not given don't compile
  AC_CHECK_FILE([$aslsrcdir/asl.h],[use_asldir=BUILD],[use_asldir=no])
elif test ! "$use_asldir" = "no"; then
  AC_CHECK_FILE([$use_asldir/amplsolver.a],[],[AC_MSG_ERROR([ASL directory \"$use_asldir\" specified, but library missing])])
  AC_CHECK_FILE([$use_asldir/asl.h],[],[AC_MSG_ERROR([ASL directory \"$use_asldir\" specified, but header files are missing])])
  use_asldir=`cd $use_asldir; pwd`
  case $build in
    *-cygwin*) use_asldir=`cygpath -w $use_asldir` ;;
  esac
fi

if test "$use_asldir" = BUILD; then
  if test -z "$ASLMAKEFLAGS" ; then
    asladdcflags=
    case $build in
      *-linux*) asladdcflags="-DNON_STDIO"
                case $CC in
                  pgcc | */pgcc)
                     asladdcflags="$asladdcflags -DNO_fpu_control";;
                esac
                  ;;
      *-cygwin*) asladdcflags="-DNON_STDIO"
                if test "$use_mingw" = yes; then
                  asladdcflags="$asladdcflags -DNo_Control87 -DNO_RUSAGE"
                fi
                ;;
    esac
#    ASLMAKEFLAGS="CFLAGS=\"\$(CFLAGS) $asladdcflags\""
# Do the following for now to make sure we are not interupting compilation
# of ASL because of warning message:
    bla=
    for i in $CFLAGS; do
      if test ! "$i" = "-Werror"; then
        bla="$bla $i"
      fi
    done
    ASLMAKEFLAGS="CFLAGS=\"$bla $asladdcflags\""

    AC_MSG_CHECKING([whether $am_make is the GNU make])
    rm -f conftest.mak
    ($am_make -v -f conftest.mak > conftest.grp) 2>&5
    if $EGREP GNU conftest.grp >&5 ; then
      AC_MSG_RESULT([yes])
      ASLMAKEFLAGS="-j 1 $ASLMAKEFLAGS"
    else
      AC_MSG_RESULT([no])
    fi
    rm -f conftest.mak conftest.grp

#    if test $USE_MAINTAINER_MODE = 'yes'; then
#      AC_MSG_WARN([Assuming that maintainer mode is only used for GNU make])
#      ASLMAKEFLAGS="-j 1 $ASLMAKEFLAGS"
#    else
#      AC_LANG_PUSH(Fortran 77)
#      AC_MSG_CHECKING([whether $am_make understands the `-j 1' flag])
#      AC_LANG_CONFTEST([[bladfaoekja:]])
#      if ($am_make -f conftest.c -j 1 bladfaoekja >&5) 2>&5 ; then
#        AC_MSG_RESULT([yes])
#        ASLMAKEFLAGS="-j 1 $ASLMAKEFLAGS"
#      else
#        AC_MSG_RESULT([no])
#        AC_MSG_CHECKING([whether $am_make understands the `-J 1' flag])
#        if ($am_make -f conftest.c -J 1 bladfaoekja >&5) 2>&5 ; then
#          AC_MSG_RESULT([yes])
#          ASLMAKEFLAGS="-J 1 $ASLMAKEFLAGS"
#        else
#          AC_MSG_RESULT([no])
#        fi
#      fi
#      AC_LANG_POP(Fortran 77)
#    fi
  fi
  ASLMAKEFILECMDS="sed 's/ar \$\$(ARFLAGS)/$AR $ARFLAGS/' makefile.u | sed 's/ar ruv/$AR $ARFLAGS/' "
  if test "$EXEEXT" = ".exe"; then
    ASLMAKEFILECMDS="$ASLMAKEFILECMDS | sed 's/a.out/a.exe/'"
  fi
  case $build in
    *-cygwin*) ASLMAKEFILECMDS="$ASLMAKEFILECMDS | sed 's%s/System_details/\`uname -sr\`/%s|System_details|\`uname -sr\`|%'" ;;
  esac
  ASLMAKEFILECMDS="$ASLMAKEFILECMDS > makefile"
fi
AC_ARG_VAR(ASLMAKEFLAGS,[Compilation flags necessary for compiling the AMPL solver library])
AC_SUBST(ASLMAKEFILECMDS)

if test ! "$use_asldir" = no; then
  AC_CHECK_LIB(dl,[dlopen],[LIBS="$LIBS -ldl"],[])
fi

AC_LANG_POP(C++)

######################################
# Equivalent int Fortran and C types #
######################################

AC_DEFINE([FORTRAN_INTEGER_TYPE],[int],[Define to the C type corresponding to Fortran INTEGER])
#AC_CHECK_SIZEOF([long])
#AC_CHECK_SIZEOF([int])
#AC_CHECK_SIZEOF([double])
#AC_MSG_CHECKING([for C type corresponding to Fortran INTEGER])
#if test $ac_cv_sizeof_long = $ac_cv_sizeof_double; then
#  AC_DEFINE([FORTRAN_INTEGER_TYPE],[int],[Define to the C type corresponding to Fortran INTEGER])
#  AC_MSG_RESULT([int])
#else
#  AC_DEFINE([FORTRAN_INTEGER_TYPE],[long])
#  AC_MSG_RESULT([long])
#fi


IPADDLIBS="$IPADDLIBS"

##########
# ASTYLE #
##########
ASTYLE=astyle
AC_SUBST(ASTYLE)
ASTYLEFLAGS="--mode=c --indent=spaces=2 --indent-cases --indent-namespaces --min-conditional-indent=1 --brackets=linux --brackets=break-closing-headers --max-instatement-indent=2"
AC_SUBST(ASTYLEFLAGS)

# The following one is for the example makefiles
case $build in
  *-cygwin*) ipoptlib='`cygpath -w \$(IPOPTLIBDIR)`\\libipopt.a';;
  *)         ipoptlib="-L\$(IPOPTLIBDIR) -lipopt";;
esac
AC_SUBST(ipoptlib)

# For ASL
AC_SUBST(AMPL_CPPFLAGS)
AMPL_CPPFLAGS="-I\$(top_builddir)/Extern/ASL -I\$(top_srcdir)/Extern/ASL/solvers"

# Add the compiler warning flags for GNU compilers
# (we can't use them during tests)
case $build in
# Cygwin seems to be giving too many error messages
  *-cygwin*) warn_cflags=""
             warn_cxxflags=""
             ;;
esac
CFLAGS="$CFLAGS $warn_cflags"
CXXFLAGS="$CXXFLAGS $warn_cxxflags"

#########################
# Makefile conditionals #
#########################

# The following variable collects the names of libraries that should be
# included into libipopt.a
AC_SUBST(IPALLLIBS)
IPALLLIBS="Algorithm/libipoptalg.a Algorithm/LinearSolvers/liblinsolvers.a Interfaces/libinterfaces.a Common/libcommon.a LinAlg/liblinalg.a LinAlg/TMatrices/libtmatrices.a"

AM_CONDITIONAL([BUILDBLAS],[test "$use_blas" = "BUILD"])
if test "$use_blas" = BUILD; then
  LIBS="$ipblaslibdir/libipblas.a $LIBS"
  IPALLLIBS="Extern/blas/libipblas.a $IPALLLIBS"
fi

# AM_CONDITIONAL([BUILDLAPACK],[test "$use_lapack" = "BUILD"])
# if test "$use_lapack" = BUILD; then
#   LIBS="$iplapacklibdir/libiplapack.a $LIBS"
#   IPALLLIBS="Extern/lapack/libiplapack.a $IPALLLIBS"
#   for file in dgesv.f dgetf2.f dgetrf.f dgetrs.f dlabad.f dlacon.f dlacpy.f dlae2.f dlaebz.f dlaev2.f dlagtf.f dlagts.f dlamch.f dlansp.f dlanst.f dlapy2.f dlaqsp.f dlarf.f dlarfg.f dlarnv.f dlartg.f dlaruv.f dlascl.f dlaset.f dlasr.f dlasrt.f dlassq.f dlaswp.f dlatps.f dopgtr.f dopmtr.f dorg2l.f dorg2r.f dppcon.f dppequ.f dpprfs.f dppsv.f dppsvx.f dpptrf.f dpptrs.f drscl.f dspev.f dspevx.f dsptrd.f dstebz.f dstein.f dsteqr.f dsterf.f ieeeck.f ilaenv.f; do
#     AC_CONFIG_LINKS($lapackobjdir/$file:$lapackdir/$file)
#   done
# fi

AM_CONDITIONAL([BUILDHSL],[test "$use_hsl" = "BUILD"])
if test "$use_hsl" = BUILD ; then
#  LDFLAGS="-L$hsllibdir $LDFLAGS"
  LIBS="$hsllibdir/libiphsl.a $LIBS"
  IPALLLIBS="Extern/HSL/libiphsl.a $IPALLLIBS"
fi
AM_CONDITIONAL([BUILDMA27],[test "$use_MA27" = BUILD])
AM_CONDITIONAL([BUILDMC19],[test "$use_MC19" = BUILD])

AM_CONDITIONAL([BUILDBLASDIR],[test "$use_blas" = BUILD])
# AM_CONDITIONAL([BUILDLAPACKDIR],[test "$use_lapack" = BUILD])
AM_CONDITIONAL([BUILDHSLDIR],[test "$use_hsl" = BUILD])

AM_CONDITIONAL([HAVE_MA27],[test ! "$use_MA27" = no])
AM_CONDITIONAL([HAVE_MC19],[test ! "$use_MC19" = no])
AM_CONDITIONAL([HAVE_PARDISO],[test ! "$use_pardiso" = no])
AM_CONDITIONAL([HAVE_TAUCS],[test ! "$use_taucs" = no])

AM_CONDITIONAL([BUILDASLDIR],[test "$use_asldir" = BUILD])
AM_CONDITIONAL([BUILDAMPLINTERFACE],[test ! "$use_asldir" = "no"])
AC_SUBST(ASLLIB)
AC_SUBST(AMPL_CPPFLAGS)
if test "$use_asldir" = BUILD; then
  ASLLIB="\$(top_builddir)/Extern/ASL/amplsolver.a"
  AMPL_CPPFLAGS="-I\$(top_builddir)/Extern/ASL -I\$(top_srcdir)/Extern/ASL/solvers"
  AC_CONFIG_LINKS(Extern/ASL/ASLsolvers:$asldir)
elif test ! "$use_asldir" = no; then
  case $build in
    *-cygwin*) AMPL_CPPFLAGS="-I\"$use_asldir\""
               ASLLIB="\"$use_asldir\amplsolver.a\"" ;;
    *)         AMPL_CPPFLAGS="-I$use_asldir"
               ASLLIB="$use_asldir/amplsolver.a";;
  esac
fi


AC_CONFIG_FILES([Makefile Common/Makefile
                 LinAlg/Makefile LinAlg/TMatrices/Makefile 
		 Interfaces/Makefile Algorithm/Makefile
		 Algorithm/LinearSolvers/Makefile
		 Apps/Makefile Apps/AmplSolver/Makefile
		 Apps/CUTErInterface/Makefile
		 Examples/hs071_f/Makefile
		 Examples/Cpp_example/Makefile
		 Examples/hs071_cpp/Makefile
		 Examples/hs071_c/Makefile
		 Extern/Makefile Extern/blas/Makefile
		 Extern/HSL/Makefile
		 Extern/ASL/Makefile])

AC_CONFIG_HEADER([inc/config.h])

AC_OUTPUT

AC_MSG_NOTICE([Configuration successful])
