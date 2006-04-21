# Copyright (C) 2006 International Business Machines..
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id: Makefile.am 669 2006-03-14 03:53:48Z andreasw $
#
# Author: Andreas Wachter    IBM      2006-04-14

# This file defines the common autoconf macros for COIN
#

# Check requirements
AC_PREREQ(2.59)

###########################################################################
#                          COIN_DEBUG_COMPILE                             #
###########################################################################

# enable the configure flag --enable-debug and set the variable
# coin_debug_compile to true or false
# This is used by COIN_PROG_CXX, COIN_PROG_CC and COIN_PROG_F77
# to determine the compilation flags.
# It define the macro COIN_DEBUG if it is chosen, and the makefile
# conditional COIN_DEBUG_COMPILE is defined

AC_DEFUN([AC_COIN_DEBUG_COMPILE],
[AC_BEFORE([$0],[AC_COIN_PROG_CXX])dnl
AC_BEFORE([$0],[AC_COIN_PROG_CC])dnl
AC_BEFORE([$0],[AC_COIN_PROG_F77])dnl

AC_MSG_CHECKING([whether we want to compile in debug mode])

AC_ARG_ENABLE([debug],
[AC_HELP_STRING([--enable-debug],
                [compile with debug options and runtime tests])],
                [case "${enableval}" in
                   yes) coin_debug_compile=true
                     AC_DEFINE([COIN_DEBUG],[1],
                               [If defined, debug sanity checks are performed during runtime])
                     ;;
                   no)  coin_debug_compile=false
                     ;;
                    *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug)
                     ;;
                 esac],
                [coin_debug_compile=false])

if test $coin_debug_compile = true; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

AM_CONDITIONAL([COIN_DEBUG_COMPILE],[test "$coin_debug_compile" = true])
]) # AC_COIN_DEBUG_COMPILE


###########################################################################
#                             COIN_PROG_CXX                               #
###########################################################################

# Find the compile command by running AC_PROG_CXX (with compiler names
# for different operating systems) and put it into CXX (unless it was
# given my the user), and find an appropriate value for CXXFLAGS

AC_DEFUN([AC_COIN_PROG_CXX],
[AC_LANG_PUSH(C++)

coin_use_cxx=yes

save_cxxflags="$CXXFLAGS"
case $build in
  *-cygwin*) comps="g++ cl" ;;
  *-mingw*)  comps="g++ cl" ;;
  *-darwin*) comps="g++ c++ CC" ;;
          *) comps="xlC aCC CC g++ c++ pgCC icpc" ;;
esac
AC_PROG_CXX([$comps])
CXXFLAGS="$save_cxxflags"

if test x"$CXXFLAGS" != x; then
  coin_user_set_cxxflags=yes
fi

AC_CACHE_CHECK([for C++ compiler options],[coin_cv_cxxflags],
[if test "$coin_user_set_cxxflags" != yes; then

  coin_add_cxxflags=
  coin_opt_cxxflags=
  coin_dbg_cxxflags=

  if test "$GXX" = "yes"; then
    case "$CXX" in
      icpc | */icpc)
        ;;
      *)
# ToDo decide about unroll-loops
        coin_opt_cxxflags="-O3"
        coin_add_cxxflags="-pipe"
        coin_dbg_cxxflags="-g"
 
        case $build in
          *-cygwin*)
            CXXFLAGS="-mno-cygwin"
            AC_TRY_LINK([],[int i=0; i++;],
                        [coin_add_cxxflags="-mno-cygwin $coin_add_cxxflags"])
            CXXFLAGS=
            ;;
        esac
        ;;
    esac
  fi
  if test -z "$coin_opt_cxxflags"; then
    case $build in
      *-cygwin*)
        case "$CXX" in
          cl | */cl)
            coin_opt_cxxflags='/Ot1'
            coin_add_cxxflags='/nologo /EHsc /GR /MT'
            coin_dbg_cxxflags='/-Yd'
            ;;
        esac
        ;;
      *-linux-*)
        case "$CXX" in
          icpc | */icpc)
            coin_opt_cxxflags="-O3 -ip"
            coin_add_cxxflags=""
            coin_dbg_cxxflags="-g"
            # Check if -i_dynamic is necessary (for new glibc library)
            CXXFLAGS=
            AC_TRY_LINK([],[int i=0; i++;],[],
                        [coin_add_cxxflags="-i_dynamic $coin_add_cxxflags"])
            ;;
          pgCC | */pgCC)
            coin_opt_cxxflags="-fast"
            coin_add_cxxflags="-Kieee -pc 64"
            coin_dbg_cxxflags="-g"
            ;;
        esac
        ;;
      *-ibm-*)
        case "$CXX" in
          xlC* | */xlC* | mpxlC* | */mpxlC*)
            coin_opt_cxxflags="-O3 -qarch=auto -qcache=auto -qhot -qtune=auto -qmaxmem=-1"
            coin_add_cxxflags="-bmaxdata:0x80000000 -qrtti=dyna"
            coin_dbg_cxxflags="-g"
            ;;
        esac
        ;;
      *-hp-*)
        case "$CXX" in
          aCC | */aCC )
            coin_opt_cxxflags="-O"
            coin_add_cxxflags="-AA"
            coin_dbg_cxxflags="-g"
            ;;
        esac
        ;;
      *-sun-*)
          coin_opt_cxxflags="-O4 -xtarget=native"
          coin_dbg_cxxflags="-g"
        ;;
    esac
  fi

  if test "$ac_cv_prog_cxx_g" = yes && test -z "$coin_dbg_cxxflags" ; then
    coin_dbg_cxxflags="-g"
  fi

  if test "$coin_debug_compile" = "true"; then
    CXXFLAGS="$coin_dbg_cxxflags $coin_add_cxxflags"
  else
    if test -z "$coin_opt_cxxflags"; then
      # Try if -O option works if nothing else is set
      CXXFLAGS="-O"
      AC_TRY_LINK([],[int i=0; i++;],[coin_opt_cxxflags="-O"])
    fi
    CXXFLAGS="$coin_opt_cxxflags $coin_add_cxxflags"
  fi
fi

# Try if CXXFLAGS works
AC_TRY_LINK([],[int i=0; i++;],[],[CXXFLAGS=])
if test -z "$CXXFLAGS"; then
  AC_MSG_WARN([The flags CXXFLAGS="$CXXFLAGS" do not work.  I will now just try '-O', but you might want to set CXXFLAGS manually.])
  CXXFLAGS='-O'
  AC_TRY_LINK([],[int i=0; i++;],[],[CXXFLAGS=])
  if test -z "$CXXFLAGS"; then
    AC_MSG_WARN([This value for CXXFLAGS does not work.  I will continue with empty CXXFLAGS, but you might want to set CXXFLAGS manually.])
  fi
fi
coin_cv_cxxflags="$CXXFLAGS"
]) # AC_CACHE_CHECK([for C++ compiler options CXXFLAGS]
CXXFLAGS="$coin_cv_cxxflags"

AC_LANG_POP(C++)
]) # AC_COIN_PROG_CXX


###########################################################################
#                             COIN_CXXLIBS                                #
###########################################################################

# Determine the C++ runtime libraries required for linking a C++ library
# with a Fortran or C compiler.  The result is available in CXXLIBS.

AC_DEFUN([AC_COIN_CXXLIBS],
[AC_REQUIRE([AC_PROG_CXX])dnl
AC_LANG_PUSH(C++)
AC_ARG_VAR(CXXLIBS,[Libraries necessary for linking C++ code with Fortran compiler])
if test -z "$CXXLIBS"; then
  if test "$GXX" = "yes"; then
    case "$CXX" in
      icpc | */icpc)
        CXXLIBS=""
        ;;
      *)
        CXXLIBS="-lstdc++ -lm" # -lgcc"
        ;;
    esac
  else
    case $build in
     *-linux-*)
      case "$CXX" in
      icpc | */icpc)
        CXXLIBS=""
             ;;
      pgCC | */pgCC)
        CXXLIBS="-lstd -lC -lc"
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
]) # AC_COIN_CXXLIBS

###########################################################################
#                           COIN_CHECK_HEADER                             #
###########################################################################

# This macro checks for a header file, but it does so without the
# standard header.  This avoids warning messages like:
#
# configure: WARNING: dlfcn.h: present but cannot be compiled
# configure: WARNING: dlfcn.h:     check for missing prerequisite headers?
# configure: WARNING: dlfcn.h: see the Autoconf documentation
# configure: WARNING: dlfcn.h:     section "Present But Cannot Be Compiled"
# configure: WARNING: dlfcn.h: proceeding with the preprocessor's result
# configure: WARNING: dlfcn.h: in the future, the compiler will take precedence

AC_DEFUN([AC_COIN_CHECK_HEADER],
[if test x"$4" = x; then
  hdr="#include <$1>"
else
  hdr="$4"
fi
AC_CHECK_HEADERS([$1],[$2],[$3],[$hdr])
]) # AC_COIN_CHECK_HEADER

###########################################################################
#                       COIN_CHECK_CXX_CHEADER                             #
###########################################################################

# This macro checks for C header files that are used from C++.  For a give
# stub (e.g., math), it first checks if the C++ library (cmath) is available.
# If it is, it defines HAVE_CMATH (or whatever the stub is).  If it is not
# available, it checks for the old C head (math.h) and defines HAVE_MATH_H
# if that one exists.

AC_DEFUN([AC_COIN_CHECK_CXX_CHEADER],
[AC_LANG_PUSH(C++)
AC_COIN_CHECK_HEADER([c$1],[$2],[$3],[$4])
if test "$ac_cv_header_c$1" != "yes"; then
  AC_COIN_CHECK_HEADER([$1.h],[$2],[$3],[$4])
fi
AC_LANG_POP(C++)
]) # AC_COIN_CHECK_CXX_CHEADER

###########################################################################
#                             COIN_PROG_CC                                #
###########################################################################

# Find the compile command by running AC_PROG_CC (with compiler names
# for different operating systems) and put it into CC (unless it was
# given my the user), and find an appropriate value for CFLAGS

AC_DEFUN([AC_COIN_PROG_CC],
[AC_LANG_PUSH(C)

coin_use_cc=yes

save_cflags="$CFLAGS"
case $build in
  *-cygwin*) comps="gcc cl" ;;
  *-mingw*)  comps="gcc cl" ;;
  *-linux-*) comps="xlc gcc cc pgcc icc" ;;
  *)         comps="xlc cc gcc pgcc icc" ;;
esac
AC_PROG_CC([$comps])
CFLAGS="$save_cflags"

if test x"$CFLAGS" != x; then
  coin_user_set_cflags=yes
fi

AC_CACHE_CHECK([for C compiler options],[coin_cv_cflags],
[if test "$coin_user_set_cflags" != yes; then

  coin_add_cflags=
  coin_opt_cflags=
  coin_dbg_cflags=

  if test "$GCC" = "yes"; then
    case "$CC" in
      icc | */icc)
        ;;
      *)
        coin_opt_cflags="-O3"
        coin_add_cflags="-pipe"
        coin_dbg_cflags="-g"

        case $build in
          *-cygwin*)
            CFLAGS="-mno-cygwin"
            AC_TRY_LINK([],[int i=0; i++;],
                        [coin_add_cflags="-mno-cygwin $coin_add_cflags"])
            CFLAGS=
          ;;
        esac
        ;;
    esac
  fi
  if test -z "$coin_opt_cflags"; then
    case $build in
      *-cygwin*)
        case "$CC" in
          cl | */cl)
            coin_opt_cflags='/Ot1'
            coin_add_cflags='/nologo'
            coin_dbg_cflags='/Yd'
            ;;
        esac
        ;;
      *-linux-*)
        case "$CC" in
          icc | */icc)
            coin_opt_cflags="-O3 -ip"
            coin_add_cflags=""
            coin_dbg_cflags="-g"
            # Check if -i_dynamic is necessary (for new glibc library)
            CFLAGS=
            AC_TRY_LINK([],[int i=0; i++;],[],
                        [coin_add_cflags="-i_dynamic $coin_add_cflags"])
            ;;
          pgcc | */pgcc)
            coin_opt_cflags="-fast"
            coin_add_cflags="-Kieee -pc 64"
            coin_dbg_cflags="-g"
            ;;
        esac
        ;;
      *-ibm-*)
        case "$CC" in
          xlc* | */xlc* | mpxlc* | */mpxlc*)
            coin_opt_cflags="-O3 -qarch=auto -qcache=auto -qhot -qtune=auto -qmaxmem=-1"
            coin_add_cflags="-bmaxdata:0x80000000"
            coin_dbg_cflags="-g"
          ;;
        esac
        ;;
      *-hp-*)
        coin_opt_cflags="-O"
        coin_add_cflags="-Ae"
        coin_dbg_cflags="-g"
        ;;
      *-sun-*)
        coin_opt_cflags="-xO4 -xtarget=native"
        coin_dbg_cflags="-g"
        ;;
      *-sgi-*)
        coin_opt_cflags="-O -OPT:Olimit=0"
        coin_dbg_cflags="-g"
        ;;
    esac
  fi

  if test "$ac_cv_prog_cc_g" = yes && test -z "$coin_dbg_cflags" ; then
    coin_dbg_cflags="-g"
  fi

  if test "$coin_debug_compile" = "true"; then
    CFLAGS="$coin_dbg_cflags $coin_add_cflags"
  else
    if test -z "$coin_opt_cflags"; then
      # Try if -O option works if nothing else is set
      CFLAGS="-O"
      AC_TRY_LINK([],[int i=0; i++;],[coin_opt_cflags="-O"],[])
    fi
    CFLAGS="$coin_opt_cflags $coin_add_cflags"
  fi
fi

# Try if CFLAGS works
AC_TRY_LINK([],[int i=0; i++;],[],[CFLAGS=])
if test -z "$CFLAGS"; then
  AC_MSG_WARN([The value CFLAGS="$CFLAGS" do not work.  I will now just try '-O', but you might want to set CFLAGS manually.])
  CFLAGS='-O'
  AC_TRY_LINK([],[int i=0; i++;],[],[CFLAGS=])
  if test -z "$CFLAGS"; then
    AC_MSG_WARN([This value for CFLAGS does not work.  I will continue with empty CFLAGS, but you might want to set CFLAGS manually.])
  fi
fi
coin_cv_cflags="$CFLAGS"
]) # AC_CACHE_CHECK([for C compiler options CXXFLAGS]
CFLAGS="$coin_cv_cflags"

AC_LANG_POP(C)
]) # AC_COIN_PROG_CC

###########################################################################
#                             COIN_PROG_F77                               #
###########################################################################

# Find the compile command by running AC_PROG_F77 (with compiler names
# for different operating systems) and put it into F77 (unless it was
# given my the user), and find an appropriate value for FFLAGS

AC_DEFUN([AC_COIN_PROG_F77],
[AC_LANG_PUSH([Fortran 77])

coin_use_f77=yes

save_fflags="$FFLAGS"
case $build in
  *-cygwin*) comps="gfortran g77 ifort" ;;
  *-mingw*)  comps="gfortran g77 ifort" ;;
  *)         comps="xlf fort77 gfortran f77 g77 pgf90 pgf77 ifort ifc" ;;
esac
AC_PROG_F77($comps)
FFLAGS="$save_fflags"

if test x"$FFLAGS" != x; then
  coin_user_set_fflags=yes
fi

AC_CACHE_CHECK([for Fortran compiler options],[coin_cv_fflags],
[if test "$coin_user_set_fflags" != yes; then

  coin_add_fflags=
  coin_opt_fflags=
  coin_dbg_fflags=

  if test "$G77" = "yes"; then
    coin_opt_fflags="-O3"
    coin_add_fflags="-pipe"
    coin_dbg_fflags="-g"
    case $build in
      *-cygwin*)
        FFLAGS="-mno-cygwin"
        AC_TRY_LINK([],[      write(*,*) 'Hello world'],
                    [coin_add_fflags="-mno-cygwin $coin_add_fflags"])
        FFLAGS=
      ;;
    esac
  else
    case $build in
      *-cygwin*)
        case $F77 in
          ifort | */ifort)
            coin_opt_fflags='/O3'
            coin_add_fflags='/nologo'
            coin_dbg_fflags='/debug'
          ;;
        esac
        ;;
      *-linux-*)
        case $F77 in
          ifc | */ifc | ifort | */ifort)
            coin_opt_fflags="-O3 -ip"
            coin_add_fflags="-cm -w90 -w95"
            coin_dbg_fflags="-g -CA -CB -CS"
            # Check if -i_dynamic is necessary (for new glibc library)
            FFLAGS=
            AC_TRY_LINK([],[      write(*,*) 'Hello world'],[],
                        [coin_add_fflags="-i_dynamic $coin_add_fflags"])
            ;;
          pgf77 | */pgf77 | pgf90 | */pgf90)
            coin_opt_fflags="-fast"
            coin_add_fflags="-Kieee -pc 64"
            coin_dbg_fflags="-g"
          ;;
        esac
        ;;
      *-ibm-*)
        case $F77 in
          xlf* | */xlf* | mpxlf* | */mpxlf* )
            coin_opt_fflags="-O3 -qarch=auto -qcache=auto -qhot -qtune=auto -qmaxmem=-1"
            coin_add_fflags="-bmaxdata:0x80000000"
            coin_dbg_fflags="-g -C"
            ;;
        esac
        ;;
      *-hp-*)
        coin_opt_fflags="+O3"
        coin_add_fflags="+U77"
        coin_dbg_fflags="-C -g"
        ;;
      *-sun-*)
        coin_opt_fflags="-O4 -xtarget=native"
        coin_dbg_fflags="-g"
        ;;
      *-sgi-*)
        coin_opt_fflags="-O5 -OPT:Olimit=0"
        coin_dbg_fflags="-g"
        ;;
    esac
  fi

  if test "$ac_cv_prog_f77_g" = yes && test -z "$coin_dbg_fflags" ; then
    coin_dbg_fflags="-g"
  fi

  if test "$coin_debug_compile" = true; then
    FFLAGS="$coin_dbg_fflags $coin_add_fflags"
  else
    if test -z "$coin_opt_fflags"; then
      # Try if -O option works if nothing else is set
      FFLAGS=-O
      AC_TRY_LINK([],[      integer i],
                  [coin_opt_fflags="-O"])
    fi
    FFLAGS="$coin_opt_fflags $coin_add_fflags"
  fi
fi

# Try if FFLAGS works
AC_TRY_LINK([],[      integer i],[],[FFLAGS=])
if test -z "$FFLAGS"; then
  AC_MSG_WARN([The flags FFLAGS="$FFLAGS" do not work.  I will now just try '-O', but you might want to set FFLAGS manually.])
  FFLAGS='-O'
  AC_TRY_LINK([],[      integer i],[],[FFLAGS=])
  if test -z "$FFLAGS"; then
    AC_MSG_WARN([This value for FFLAGS does not work.  I will continue with empty FFLAGS, but you might want to set FFLAGS manually.])
  fi
fi
coin_cv_fflags="$FFLAGS"
]) # AC_CACHE_CHECK([for Fortran compiler options FFLAGS]
FFLAGS="$coin_cv_fflags"

AC_LANG_POP([Fortran 77])
]) # AC_COIN_PROG_F77

###########################################################################
#                           COIN_F77_WRAPPERS                             #
###########################################################################

# Calls autoconfs AC_F77_WRAPPERS and does additional corrections to FLIBS

AC_DEFUN([AC_COIN_F77_WRAPPERS],
[AC_BEFORE([AC_COIN_PROG_F77],[$0])dnl
AC_BEFORE([AC_PROG_F77],[$0])dnl

AC_LANG_PUSH([Fortran 77])

AC_F77_WRAPPERS

# This is to correct a missing exclusion in autoconf 2.59
if test x"$FLIBS" != x; then
  my_flibs=
  for flag in $FLIBS; do
    case flag in
      -lcrt*.o) ;;
             *) my_flibs="$my_flibs $flag" ;;
    esac
  done
  FLIBS="$my_flibs"
fi

case $build in
# The following is a fix to define FLIBS for ifort on Windows
   *-cygwin*)
     case $F77 in
       ifort | */ifort)
           FLIBS="/link libifcorert.lib $LIBS /NODEFAULTLIB:libc.lib";;
     esac;;
   *-hp-*)
       FLIBS="$FLIBS -lm";;
   *-ibm-*)
       FLIBS=`echo $FLIBS | sed 's/-lc)/-lc/g'` ;;
   *-linux-*)
     case "$F77" in
       pgf77 | */pgf77 | pgf90 | */pgf90)
# ask linker to go through the archives multiple times
# (the Fortran compiler seems to do that automatically...
         FLIBS="-Wl,--start-group $FLIBS -Wl,--end-group" ;;
     esac
esac

]) # AC_COIN_F77_WRAPPERS

###########################################################################
#                          COIN_ADD_WARN_FLAGS                            #
###########################################################################

# This macro can be used at the end of the configure.ac script to add
# warning options to the compiler flags.  This has to be done so late,
# because otherwise the tests don't work.  It is assumed that
# AC_COIN_PROG_... has been used before to set up the compilers.
# (so far, we don't set warning flags for Fortran compilers...)

AC_DEFUN([AC_COIN_ADD_WARN_FLAGS],
[AC_BEFORE([AC_COIN_PROG_CXX],[$0])dnl
AC_BEFORE([AC_COIN_PROG_CC],[$0])dnl
AC_BEFORE([AC_COIN_PROG_F77],[$0])dnl

# ToDo find more compiler options
if test x"$coin_use_cxx" = xyes && test x"$coin_user_set_cxxflags" != xyes; then

  AC_MSG_CHECKING([for warning flags for C++ compiler])
  coin_warn_cxxflags=

  if test "$GXX" = "yes"; then
    case "$CXX" in
      icpc | */icpc)
        ;;
      *)
        coin_warn_cxxflags="-pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wconversion"
        ;;
    esac
  fi

  AC_MSG_RESULT([$coin_warn_cxxflags])
  CXXFLAGS="$CXXFLAGS $coin_warn_cxxflags"
fi

if test x"$coin_use_cc" = xyes && test x"$coin_user_set_cflags" != xyes; then

  AC_MSG_CHECKING([for warning flags for C compiler])
  coin_warn_cflags=

  if test "$GCC" = "yes"; then
    case "$CC" in
      icc | */icc)
        ;;
      *)
	coin_warn_cflags="-pedantic-errors -Wimplicit -Wparentheses -Wsequence-point -Wreturn-type -Wcast-qual -Wall"
        ;;
    esac
  fi

  AC_MSG_RESULT([$coin_warn_cflags])
  CFLAGS="$CFLAGS $coin_warn_cflags"
fi

if test x"$coin_use_f77" = xyes && test x"$coin_user_set_fflags" != xyes; then
  AC_MSG_CHECKING([for warning flags for Fortran compiler])
  coin_warn_fflags=
  AC_MSG_RESULT([$coin_warn_fflags])
fi

]) # AC_COIN_ADD_WARN_FLAGS

###########################################################################
#                         COIN_INIT_AUTO_TOOLS                            #
###########################################################################

# Initialize the auto tools automake and libtool, with all
# modifications we want for COIN packages.
#
# This also defines the AC_SUBST variables:
# pkg_source_dir     absolute path to source code for this package
# pkg_bin_dir        absolute path to the directory where binaries are
#                    going to be installed (prefix/bin)
# pkg_lib_dir        absolute path to the directory where libraries are
#                    going to be installed (prefix/lib)
# pkg_include_dir    absolute path to the directory where the header files
#                    are installed (prefix/include)

AC_DEFUN([AC_COIN_INIT_AUTO_TOOLS],
[AC_BEFORE([AC_COIN_PROG_CXX],[$0])dnl
AC_BEFORE([AC_COIN_PROG_CC],[$0])dnl
AC_BEFORE([AC_COIN_PROG_F77],[$0])dnl

# Stuff for automake
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

# Stuff for libtool
AC_COIN_PROG_LIBTOOL

# If maintainer mode is chosen, we make sure that the correct versions
# of the tools are used, and that we know where libtoo.m4 is (to
# recreate acinclude.m4)

AC_SUBST(LIBTOOLM4)
LIBTOOLM4=

if test "$enable_maintainer_mode" = yes; then

  # Check if we have autoconf
  AC_CHECK_PROG([have_autoconf],[autoconf],[yes],[no])
  if test $have_autoconf = no; then
    AC_MSG_ERROR([You specified you want to use maintainer mode, but I cannot find autoconf in your path.])
  fi

  # Check whether autoconf is the correct version
  correct_version='2.59'
  grep_version=`echo  $correct_version | sed -e 's/\\./\\\\\\./g'`
  AC_MSG_CHECKING([whether we are using the correct version ($correct_version) of autoconf])
  autoconf --version > confauto.out 2>&1
  if $EGREP $grep_version confauto.out >/dev/null 2>&1; then
    AC_MSG_RESULT([yes])
  else
    rm -f confauto.out
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([You don't have the correct version of autoconf as the first one in your path.])
  fi
  rm -f confauto.out

  # Check if we have automake
  AC_CHECK_PROG([have_automake],[automake],[yes],[no])
  if test $have_automake = no; then
    AC_MSG_ERROR([You specified you want to use maintainer mode, but I cannot find automake in your path.])
  fi
  
  # Check whether automake is the correct version
  correct_version='1.9.6'
  grep_version=`echo  $correct_version | sed -e 's/\\./\\\\\\./g'`
  AC_MSG_CHECKING([whether we are using the correct version ($correct_version) of automake])
  automake --version > confauto.out 2>&1
  if $EGREP $grep_version confauto.out >/dev/null 2>&1; then
    AC_MSG_RESULT([yes])
  else
    rm -f confauto.out
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([You don't have the correct version of automake as the first one in your path.])
  fi
  rm -f confauto.out

  # Check if we can find the libtool file
  if test "${LIBTOOLPREFIX:+set}" != set; then
    for p in $HOME ; do
      AC_CHECK_FILE([$p/share/aclocal/libtool.m4],
                    [LIBTOOLM4="$p/share/aclocal/libtool.m4"
                     LIBTOOLPREFIX="$p"],)
      if test x"$LIBTOOLM4" != x; then
        break;
      fi
    done
    if test x"$LIBTOOLM4" = x; then
      AC_MSG_ERROR([You specified you want to use maintainer mode, but I cannot find the file libtool.m4 on your system.  Please set the prefix of the location of the correct file with the LIBTOOLPREFIX variable, so that it is in $LIBTOOLPREFIX/share/aclocal.  We assume here that it is the plain version obtained from the GNU tarball.])
    fi
  else
    AC_CHECK_FILE([$LIBTOOLPREFIX/share/aclocal/libtool.m4],
                  [LIBTOOLM4="$LIBTOOLPREFIX/share/aclocal/libtool.m4"],
                  [AC_MSG_ERROR([You specified LIBTOOLPREFIX, but I cannot find the file libtool.m4 in $LIBTOOLPREFIX/share/aclocal.])])
  fi

  # Check if this is the correct version of libtool (with escaped dots)
  correct_version='1.5.22'
  grep_version=`echo  $correct_version | sed -e 's/\\./\\\\\\./g'`
  AC_CHECK_FILE([$LIBTOOLPREFIX/share/libtool/ltmain.sh],
	        [have_ltmain=yes],
                [have_ltmain=no])
  AC_MSG_CHECKING([whether we are using the correct version ($correct_version) of libtool.])
  if test $have_ltmain = yes; then
    if $EGREP $grep_version $LIBTOOLPREFIX/share/libtool/ltmain.sh >/dev/null 2>&1; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([You don't have the correct version of libtool.  Please set LIBTOOLPREFIX to the correct installation prefix, so that the correct version of ltmain.sh is in $LIBTOOLPREFIX/share/libtool.])
    fi
  else
    AC_MSG_ERROR([I cannot find the file ltmain.sh in $LIBTOOLPREFIX/share/libtool])
  fi  
fi

# helpful variable for the base directory of this package
pkg_source_dir=`cd $srcdir; pwd`

# Stuff for example Makefiles
full_prefix=`echo $exec_prefix | pwd`
AC_SUBST(pkg_lib_dir)
pkg_lib_dir=$full_prefix/lib
AC_SUBST(pkg_include_dir)
pkg_include_dir=$full_prefix/include
AC_SUBST(pkg_bin_dir)
pkg_bin_dir=$full_prefix/bin

]) # AC_COIN_INIT_AUTO_TOOLS

###########################################################################
#                           COIN_PROG_LIBTOOL                             #
###########################################################################

# Setup the libtool stuff together with any modifications to make it
# work on additional platforms

AC_DEFUN([AC_COIN_PROG_LIBTOOL],
[AC_REQUIRE([AC_COIN_DLFCN_H])

# We check for this header here in a non-standard way to avoid warning
# messages
AC_PROG_LIBTOOL

# Fix bugs in libtool script for Windows native compilation:
# - cygpath is not correctly quoted in fix_srcfile_path
# - paths generated for .lib files is not run through cygpath -w


# - lib includes subdirectory information; we want to replace
#
# old_archive_cmds="lib /OUT:\$oldlib\$oldobjs\$old_deplibs"
#
# by
#
# old_archive_cmds="echo \$oldlib | grep .libs >/dev/null; if test \$? = 0; then cd .libs; lib /OUT:\`echo \$oldlib\$oldobjs\$old_deplibs | sed -e s@\.libs/@@g\`; cd .. ; else lib /OUT:\$oldlib\$oldobjs\$old_deplibs ; fi"
#
#          -e 's%old_archive_cmds="lib /OUT:\\\$oldlib\\\$oldobjs\\\$old_deplibs"%old_archive_cmds="echo \\\$oldlib \| grep .libs >/dev/null; if test \\\$? = 0; then cd .libs; lib /OUT:\\\`echo \\\$oldlib\\\$oldobjs\\\$old_deplibs \| sed -e s@\\.libs/@@g\\\`; cd .. ; else lib /OUT:\\\$oldlib\\\$oldobjs\\\$old_deplibs; fi"%' \

# The following was a hack for chaniing @BACKSLASH to \
#          -e 'sYcompile_command=`\$echo "X\$compile_command" | \$Xsed -e '"'"'s%@OUTPUT@%'"'"'"\$output"'"'"'%g'"'"'`Ycompile_command=`\$echo "X\$compile_command" | \$Xsed -e '"'"'s%@OUTPUT@%'"'"'"\$output"'"'"'%g'"'"' | \$Xsed -e '"'"'s%@BACKSLASH@%\\\\\\\\\\\\\\\\%g'"'"'`Y' \

case $build in
  *-cygwin*)
  case "$CXX" in
    cl | */cl) 

      sed -e 's|fix_srcfile_path=\"`cygpath -w \"\$srcfile\"`\"|fix_srcfile_path=\"\\\`cygpath -w \\\"\\$srcfile\\\"\\\`\"|' \
          -e 's|fix_srcfile_path=\"\"|fix_srcfile_path=\"\\\`cygpath -w \\\"\\$srcfile\\\"\\\`\"|' \
          -e 's%compile_deplibs=\"\$dir/\$old_library \$compile_deplibs\"%compile_deplibs="`cygpath -w \$dir/\$old_library | sed -e '"'"'sY\\\\\\\\Y/Yg'"'"'` \$compile_deplibs\"'% \
          -e 's%compile_deplibs=\"\$dir/\$linklib \$compile_deplibs\"%compile_deplibs="`cygpath -w \$dir/\$linklib | sed -e '"'"'sY\\\\\\\\Y/Yg'"'"'` \$compile_deplibs\"'% \
      libtool > conftest.bla

      mv conftest.bla libtool
      chmod 755 libtool
      ;;
  esac
esac
]) # AC_COIN_PROG_LIBTOOL

# This is a trick to force the check for the dlfcn header to be done before
# the checks for libtool
AC_DEFUN([AC_COIN_DLFCN_H],
[AC_LANG_PUSH(C)
AC_COIN_CHECK_HEADER([dlfcn.h])
AC_LANG_POP(C)
]) # AC_COIN_DLFCN_H

###########################################################################
#                            COIN_RPATH_FLAGS                             #
###########################################################################

# This macro, in case shared objects are used, defines a variable
# RPATH_FLAGS that can be used by the linker to hardwire the library
# search path for the given directories.  This is useful for example
# Makefiles

AC_DEFUN([AC_COIN_RPATH_FLAGS],
[RPATH_FLAGS=

if test "$GXX" = "yes"; then
  RPATH_FLAGS=
  for dir in $1; do
    RPATH_FLAGS="$RPATH_FLAGS -Wl,--rpath -Wl,$dir"
  done
else
  case $build in
    *-linux-*)
      case "$CXX" in
      icpc | */icpc)
        RPATH_FLAGS=
        for dir in $1; do
          RPATH_FLAGS="$RPATH_FLAGS -Wl,--rpath -Wl,$dir"
        done
      esac ;;
    *-ibm-*)
      case "$CXX" in
      xlC* | */xlC* | mpxlC* | */mpxlC*)
        RPATH_FLAGS=nothing ;;
      esac ;;
    *-hp-*)
        RPATH_FLAGS=nothing ;;
    *-mingw32)
        RPATH_FLAGS=nothing ;;
    *-sun-*)
        RPATH_FLAGS=
        for dir in $1; do
          RPATH_FLAGS="$RPATH_FLAGS -R$dir"
        done
   esac
fi

if test "$RPATH_FLAGS" = ""; then
  AC_MSG_WARN([Could not automatically determine how to tell the linker about automatic inclusion of the path for shared libraries.  The test examples might not work if you link against shared objects.  You will need to set the LD_LIBRARY_PATH or LIBDIR variable manually.])
fi
if test "$RPATH_FLAGS" = "nothing"; then
  RPATH_FLAGS=
fi

AC_SUBST(RPATH_FLAGS)
]) # AC_COIN_RPATH_FLAGS
