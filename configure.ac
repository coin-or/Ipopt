# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

## SGI: ./configure CC=cc F77=f77 CXX=CC CXXFLAGS="-g -LANG:std -DOLD_C_HEADERS" AR="CC -ar" ARFLAGS="-o "
## AIX: ./configure CC=xlc F77=xlf CXX=xlC --with-lapack=$HOME/lib/liblapack.a CXXFLAGS="-g -qrtti=dyna" CXXLIBS="-lC -lc"

AC_PREREQ(2.59)

AC_INIT([ipopt],[3.0beta],[coin-ipopt@www-124.ibm.com])

AC_COPYRIGHT([
Copyright 2004 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package IPOPT which is distributed
under the Common Public License.])

AC_REVISION(beta)

AC_CONFIG_SRCDIR(Interfaces/IpVector.hpp)

# Get the system type
AC_CANONICAL_BUILD

if test -z "$AR"; then
  AR="ar"
fi
AC_ARG_VAR(AR,[Command for handling archives (default is `ar')])

# Check if user wants to produce debugging code
AC_ARG_ENABLE([debug],
  AC_HELP_STRING([--enable-debug], [compile with debug options and runtime tests]),
  [case "${enableval}" in
    yes) debugcompile=true ;;
    no)  debugcompile=false ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
  esac],[debugcompile=false])

######################
# C Compiler Options #
######################

if test "$debugcompile" = "true"; then
  AC_DEFINE([IP_DEBUG],[1],[If defined, debug sanity checks are performed during runtime])
  if test -z "$CXXFLAGS"; then
    CXXFLAGS="-g -pedantic-errors -Wimplicit -Wparentheses -Wsequence-point -Wreturn-type -Wcast-qual -Werror -pipe"
  fi

  if test -z "$CFLAGS"; then
    CFLAGS="-g -pedantic-errors -Wimplicit -Wparentheses -Wsequence-point -Wreturn-type -Wcast-qual -Werror -pipe"
  fi

  if test -z "$FFLAGS"; then
    FFLAGS="-g"
  fi
else
  if test -z "$CXXFLAGS"; then
    CXXFLAGS="-O3 -pipe"
  fi

  if test -z "$CFLAGS"; then
    CFLAGS="-O3 -pipe"
  fi

  if test -z "$FFLAGS"; then
    FFLAGS="-O3"
  fi
fi

SAVE_CFLAGS="$CFLAGS"
#AC_PROG_CC([icc cl pgcc xlc cc gcc cl])
AC_PROG_CC
CFLAGS="$SAVE_CFLAGS"

SAVE_CXXFLAGS="$CXXFLAGS"
#AC_PROG_CXX([icpc cl KCC CC cxx cc++ xlC aCC c++ g++ egcs gcc])
AC_PROG_CXX
CXXFLAGS="$SAVE_CXXFLAGS"

AC_PROG_RANLIB

AC_PREFIX_DEFAULT([`pwd`])

AC_PROG_LN_S 

# Stuff for automake
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

###################################
#   Fortran 77 compiler options   #
###################################

AC_LANG_PUSH(Fortran 77)
SAVE_FFLAGS="$FFLAGS"
#AC_PROG_F77(ifort ifc pgf77 xlf fort77 f77 g77 df)
AC_PROG_F77
FFLAGS="$SAVE_FFLAGS"

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
esac

##########################################################################

# Additional libraries and packages

IPADDLIBS=
AC_SUBST(IPADDLIBS)

##########
#  BLAS  #
##########

AC_ARG_WITH([blas],
            AC_HELP_STRING([--with-blas],
                           [specify BLAS library (or BUILD for compilation)]),
            [use_blas=$withval], [use_blas=])
# Check if user supplied option makes sense
if test ! -z "$use_blas"; then
  if test "$use_blas" = "BUILD"; then
    AC_CHECK_FILE([$srcdir/OTHERS/blas/dgemm.f],[],
                  [AC_MSG_ERROR([option \"BUILD\" specified for blas, but sources are missing])])
  else
    AC_MSG_CHECKING([whether user supplied BLASLIB=\"$use_blas\" works])
    LIBS="$use_blas $LIBS"
    IPADDLIBS="$use_blas $IPADDLIBS"
    AC_TRY_LINK([],[      call dgemm],[AC_MSG_RESULT([yes])],
                [AC_MSG_RESULT([no]); AC_MSG_ERROR([user supplied BLAS library \"$use_blas\" does not work])])
  fi
else
# Try to autodetect the library for blas based on build system
  case $build in
    *-sgi-*)   AC_CHECK_LIB([complib.sgimath],[dgemm],
                          [use_blas=-lcomplib.sgimath;
                           LIBS="-lcomplib.sgimath $LIBS";
                           IPADDLIBS="-lcomplib.sgimath $IPADDLIBS"]);;
    *-sun-*) SAVE_LDFLAGS="$LDFLAGS"
             AC_MSG_CHECKING([whether -xlic_lib=sunperf has BLAS])
             LDFLAGS="-xlic_lib=sunperf $LDFLAGS"
             AC_TRY_LINK([],[      call dgemm],
                         [AC_MSG_RESULT([yes]); use_blas=-xlic_lib=sunperf],
                         [AC_MSG_RESULT([no]); LDFLAGS="$SAVE_LDFLAGS"]) ;;
  esac
  if test -z "$use_blas"; then
    AC_CHECK_LIB(blas,dgemm,[use_blas=-lblas;
                             LIBS="-lblas $LIBS";
                             IPADDLIBS="-lblas $IPADDLIBS"],[])
  fi
  if test -z "$use_blas"; then
    AC_CHECK_FILE([$srcdir/OTHERS/blas/dgemm.f],[use_blas=BUILD],
                  [AC_MSG_ERROR([cannot find BLAS, and sources are missing in $srcdir/OTHERS/blas/])])
  fi
fi

##########
# LAPACK #
##########

AC_ARG_WITH([lapack],
            AC_HELP_STRING([--with-lapack],
                           [specify LAPACK library (or BUILD for compilation)]),
            [use_lapack=$withval], [use_lapack=])
# Check if user supplied option makes sense
if test ! -z "$use_lapack"; then
  if test "$use_lapack" = "BUILD"; then
    AC_CHECK_FILE([$lapacksrcdir/dlarf.f],[],
                  [AC_MSG_ERROR([option \"BUILD\" specified for lapack, but sources are missing])])
  else
    AC_MSG_CHECKING([whether user supplied LAPACKLIB=\"$use_lapack\" works])
    LIBS="$use_lapack $LIBS"
    IPADDLIBS="$use_lapack $IPADDLIBS"
    AC_TRY_LINK([],[      call dlasr],[AC_MSG_RESULT([yes])],
                [AC_MSG_RESULT([no]); AC_MSG_ERROR([user supplied LAPACK library \"$use_lapack\" does not work])])
  fi
else
# First try to see if LAPACK is already available with BLAS library
  AC_MSG_CHECKING([whether LAPACK is already available with BLAS library])
  AC_TRY_LINK([],[      call dgetrf()],
            [AC_MSG_RESULT([yes]); use_lapack=ok],
            [AC_MSG_RESULT([no])])
  if test -z "$use_lapack"; then
    # Try to autodetect the library for lapack based on build system
    case $build in
      *-sgi-*)   AC_CHECK_LIB([complib.sgimath],[dgetrf],
                            [use_lapack=-lcomplib.sgimath;
                             LIBS="-lcomplib.sgimath $LIBS";
                             IPADDLIBS="-lcomplib.sgimath $IPADDLIBS"]);;
      *-sun-*) SAVE_LDFLAGS="$LDFLAGS"
               AC_MSG_CHECKING([whether -xlic_lib=sunperf has LAPACK])
               LDFLAGS="-xlic_lib=sunperf $LDFLAGS"
               AC_TRY_LINK([],[      call dgetrf],
                           [AC_MSG_RESULT([yes]); use_lapack=-xlic_lib=sunperf],
                           [AC_MSG_RESULT([no]); LDFLAGS="$SAVE_LDFLAGS"]) ;;
    esac
  fi
  if test -z "$use_lapack"; then
    AC_CHECK_LIB(lapack,dgetrf,[use_lapack=-llapack;
                                LIBS="-llapack $LIBS";
                                IPADDLIBS="-llapack $IPADDLIBS"],[])
  fi
  if test -z "$use_lapack"; then
    AC_CHECK_FILE([$lapacksrcdir/dlarf.f],[use_lapack=BUILD],
                  [AC_MSG_ERROR([cannot find LAPACK, and sources are missing in $lapacksrcdir])])
  fi
fi

#######
# HSL #
#######

LIBS="-L${HOME}/lib -lhsl $LIBS"
IPADDLIBS="-L${HOME}/lib -lhsl $IPADDLIBS"



AC_LANG_POP(Fortran 77)

# For the remaining tests we need to include FLIBS into LIBS, so that
# the C compiler can link programs
LIBS="$LIBS $FLIBS"


######################################
# Equivalent int Fortran and C types #
######################################

AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([double])
AC_MSG_CHECKING([for C type corresponding to Fortran INTEGER])
if test $ac_cv_sizeof_long = $ac_cv_sizeof_double; then
  AC_DEFINE([INTEGER_TYPE],[int],[Define to the C type corresponding to Fortran INTEGER])
  AC_MSG_RESULT([int])
else
  AC_DEFINE([INTEGER_TYPE],[long])
  AC_MSG_RESULT([long])
fi

######################
# C++ link libraries #
######################
if test -z "$CXXLIBS"; then
  CXXLIBS="-L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2/../../.. -lstdc++ -lm -lgcc_s -lgcc"
fi
AC_ARG_VAR(CXXLIBS,[Libraries necessary for linking C++ code with Fortran compiler])

IPADDLIBS="$IPADDLIBS $CXXLIBS"

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


AC_CONFIG_FILES([Makefile Common/Makefile Algorithm/Makefile
                 LinAlgImpl/Makefile LinAlgImpl/Serial/Makefile 
                 PDSystemImpl/Makefile PDSystemImpl/PDFullSpace/Makefile
		 Apps/Makefile Apps/AmplSolver/Makefile
		 Apps/StdInterface/Makefile
		 Apps/CompositeInterface/Makefile
		 Apps/CUTErInterface/Makefile
		 Examples/Fortran_example/Makefile
		 Examples/Cpp_example/Makefile])

AC_CONFIG_HEADER([inc/config.h])

AC_OUTPUT

AC_MSG_NOTICE([Configuration successful])
