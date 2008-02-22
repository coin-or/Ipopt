# Copyright (C) 2004, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Ipopt],[3.3dev],[http://projects.coin-or.org/Ipopt/newticket])

AC_COPYRIGHT(
[Copyright 2004, 2007 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package IPOPT which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/Common/IpDebug.hpp)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Ipopt)

# We need egrep and ln -s somewherw
AC_PROG_EGREP
AC_PROG_LN_S

# Get the name of the C compiler and appropriate compiler options
AC_COIN_PROG_CC

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX
# We add one more definition in if we are using g++ compilers and are
# in debug mode
if test "$coin_debug_compile" = true; then
  if test "$GXX" = yes; then
    CXXFLAGS="$CXXFLAGS -DCHECK_SMARTPTR"
  fi
fi

# Get the C++ runtime libraries in case we want to link with a C or
# Fortran compiler
AC_COIN_CXXLIBS

# Get the name of the Fortran compiler and appropriate compiler options
AC_COIN_PROG_F77

# Find out how to call Fortran from C and determine Fortran runtime libraries
AC_COIN_F77_WRAPPERS

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

# Some extra stuff for blue gene on Laci's request
AC_SUBST(BLUEGENEEXTRA)
BLUEGENEEXTRA="echo"
AC_ARG_ENABLE([bg],,
  [case "${enableval}" in
    yes)  BLUEGENEEXTRA="chmod 644 *.o; cp ../Extern/ASL/amplsolver.a ../lib/libamplsolver.a"
          is_bg=yes
          ;;
  esac],[is_bg=no])

# A Fix for the AIX compiler
AM_CONDITIONAL([AIX_FIX_COMPILER_BUG], test "x$AIX_FIX_COMPILER_BUG" != x)

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(Ipopt)

#############################################################################
#                             Stuff for examples                            #
#############################################################################

# find out how long an int pointer is to know if we need INTEGER*4 or
# INTEGER*8 in Fortran to capture pointers.
#
# This test doesn't work on BlueGene
if test "$is_bg" = no; then
  AC_LANG_PUSH(C)
  AC_CHECK_SIZEOF(int *)
  AC_LANG_POP(C)
else
  ac_cv_sizeof_int_p=4
fi
AC_SUBST(BITS_PER_POINTER)
AC_SUBST(BIT32FCOMMENT)
AC_SUBST(BIT64FCOMMENT)
case "$ac_cv_sizeof_int_p" in
  4 | 4?) BITS_PER_POINTER=32
     BIT32FCOMMENT=''
     BIT64FCOMMENT='C' ;;
  8 | 8?) BITS_PER_POINTER=64
     BIT32FCOMMENT='C'
     BIT64FCOMMENT='' ;;
  *) AC_MSG_ERROR([Unknown length of int *]);;
esac

########################################################
## Set variables for the archiver for C++ object files #
## (for AIX in 64bit mode also on SGI, SUN...)         #
########################################################

if test -z "$CXXAR"; then
  if test ! "$GXX" = "yes"; then
    case $build in
      *-sun-*)
        CXXAR="CC -xar -o"
        ;;
      *-sgi-*)
        CXXAR="CC -ar -o"
        ;;
    esac
  fi
  if test -z "$CXXAR"; then
    CXXAR="$AR cru"
  fi
fi
AC_ARG_VAR(CXXAR,[Command for making archives from C++ objects (default is AR crv) - only needed for Example])
AC_MSG_NOTICE([Using CXXAR=\"$CXXAR\".])

############################################################################
############################################################################
#                   Stuff that we need for C++ programs                    #
############################################################################
############################################################################

AC_LANG_PUSH(C++)

################
# Header Files #
################


AC_COIN_CHECK_CXX_CHEADER(math)
AC_COIN_CHECK_CXX_CHEADER(ctype)
AC_COIN_CHECK_CXX_CHEADER(assert)
AC_COIN_CHECK_CXX_CHEADER(stdio)
AC_COIN_CHECK_CXX_CHEADER(stdlib)
AC_COIN_CHECK_CXX_CHEADER(stdarg)
AC_COIN_CHECK_CXX_CHEADER(string)
AC_COIN_CHECK_CXX_CHEADER(float)
AC_COIN_CHECK_CXX_CHEADER(time)

###################
# Function finite #
###################

# First check for a few headers

case $build in
  *-cygwin*) ;;
  *) AC_COIN_CHECK_CXX_CHEADER([ieeefp]) ;;
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

###########################
# Random number generator #
###########################

AC_CHECK_DECL([drand48],[AC_DEFINE([HAVE_DRAND48],[1],
                                   [Define to 1 if function drand48 is available])],,[
#ifdef HAVE_CSTDLIB
# include <cstdlib>
#else
# ifdef HAVE_STDLIB_H
#  include <stdlib.h>
# endif
#endif
])

AC_CHECK_DECL([rand],[AC_DEFINE([HAVE_RAND],[1],
                                [Define to 1 if function rand is available])],,[
#ifdef HAVE_CSTDLIB
# include <cstdlib>
#else
# ifdef HAVE_STDLIB_H
#  include <stdlib.h>
# endif
#endif
])

AC_CHECK_DECL([srand],[AC_DEFINE([HAVE_STD__RAND],[1],
                                 [Define to 1 if function std::rand is available])],,[
#ifdef HAVE_CSTDLIB
# include <cstdlib>
#else
# ifdef HAVE_STDLIB_H
#  include <stdlib.h>
# endif
#endif
using namespace std;])


##########################################################################

# Additional libraries and packages

# We need to restore the value of LIBS before those test, otherwise libtool
# will hick-up
LIBS_without_ADDLIBS="$LIBS"

#######################
#        BLAS         #
#######################

AC_COIN_HAS_BLAS
if test $coin_has_blas = no; then
  AC_MSG_ERROR([cannot find BLAS])
fi

#######################
#       LAPACK        #
#######################

AC_COIN_HAS_LAPACK
if test $coin_has_lapack = no; then
  AC_MSG_WARN([Compiling code without LAPACK.  Certain options (e.g., quasi-Newton) will not work.])
fi

#######################
# HARWELL SUBROUTINES #
#######################

hslobjdir=../ThirdParty/HSL
hslsrcdir=$srcdir/$hslobjdir

AC_DEFUN([MY_CHECK_HSL],
         [AC_MSG_CHECKING([for $1 in $use_hsl])
          use_$2=
          AC_F77_FUNC($1,cfunc)
          AC_TRY_LINK([extern "C" {void $cfunc();}],[$cfunc()],
                      [AC_DEFINE(HAVE_$2,1,[Define to 1 if $2 is available])
		       AC_MSG_RESULT([yes]); use_$2=yes],
                      [AC_MSG_RESULT([no]); use_$2=no])])

AC_DEFUN([MY_CHECK_HSLFILE],
         [AC_COIN_CHECK_FILE([$hslsrcdir/$1.f],
                        [AC_DEFINE(HAVE_$2,1,[Define to 1 if $2 is available]) use_$2=BUILD],
                        [use_$2=no])])

AC_MSG_NOTICE([Checking for Harwell subroutine library])

AC_ARG_WITH([hsl],
            AC_HELP_STRING([--with-hsl],
                           [specify Harwell library (or BUILD for compilation)]),
            [use_hsl=$withval], [use_hsl=BUILD])

# First check if there is a working HSL ThirdParty project
if test "$use_hsl" = BUILD; then
  AC_COIN_CHECK_FILE([$hslobjdir/Makefile],[],[use_hsl=])
fi

if test "$use_hsl" = BUILD; then
  # First check if there is a working HSL ThirdParty project
  AC_COIN_CHECK_FILE([$hslobjdir/Makefile],[],[use_hsl=unavailable])

  if test "$use_hsl" = BUILD; then
    # Check which files are available
    MY_CHECK_HSLFILE(ma27ad,MA27)
    ## This test is already done in ThirdParty configure
    #  if test "$use_MA27" = BUILD; then
    #    AC_MSG_CHECKING([whether MA27 is threadsafe])
    #    if $EGREP MA27UD $hslsrcdir/ma27ad.f >/dev/null 2>&1; then
    #      AC_MSG_RESULT(yes)
    #    else
    #      AC_MSG_RESULT(no)
    #      AC_MSG_ERROR([Please obtain the threadsafe (newer) version of MA27])
    #    fi
    #  fi
    MY_CHECK_HSLFILE(mc19ad,MC19)
    MY_CHECK_HSLFILE(ma57ad,MA57)
    MY_CHECK_HSLFILE(ma28ad,MA28)
    coin_need_flibs=yes
  fi
else
  if test "$use_hsl" != no; then
    SAVE_LIBS="$LIBS"
    LIBS="$use_hsl $LIBS $FLIBS"
    ADDLIBS="$use_hsl $ADDLIBS"
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
    MY_CHECK_HSL([ma57ad],[MA57])
    MY_CHECK_HSL([ma28ad],[MA28])
    LIBS="$use_hsl $SAVE_LIBS"
  fi
fi

###########
# PARDISO #
###########

AC_LANG_PUSH([Fortran 77])
AC_ARG_WITH([pardiso],
            AC_HELP_STRING([--with-pardiso],
                           [specify Pardiso library]),
            [use_pardiso=$withval], [use_pardiso=no])

if test "$use_pardiso" != "no"; then
  LIBS="$use_pardiso $LIBS"
  ADDLIBS="$use_pardiso $ADDLIBS"

# ToDo put this back in
  AC_MSG_CHECKING([whether user-supplied Pardiso library \"$use_pardiso\" works])
  AC_TRY_LINK([],[      call PARDISO()],
              [AC_MSG_RESULT(yes)],
              [AC_MSG_ERROR([Pardiso library $use_pardiso does not seem to work])])

# Test if this is the parallel version of Pardiso
  AC_DEFINE(HAVE_PARDISO,1,[Define to 1 if Pardiso is available])
  
  AC_MSG_CHECKING([whether Pardiso library is the parallel version])
  AC_TRY_LINK([],[      call pardiso_exist_parallel()],
              [AC_MSG_RESULT(yes)
               AC_DEFINE(HAVE_PARDISO_PARALLEL,1,
                [Define to 1 if you are using the parallel version of Pardiso])],
              [AC_MSG_RESULT(no)])
fi
AC_LANG_POP([Fortran 77])

#########
# TAUCS #
#########

AC_ARG_WITH([taucs],
            AC_HELP_STRING([--with-taucs],
                           [specify TAUCS library and dependencies]),
            [use_taucs=$withval], [use_taucs=no])

if test "$use_taucs" != "no"; then
  LIBS="$use_taucs $LIBS"
  ADDLIBS="$use_taucs $ADDLIBS"

  AC_DEFINE(HAVE_TAUCS,1,
            [Define to 1 if TAUCS is available])

  AC_ARG_VAR(TAUCS_INCFLAGS,[Flags for including taucs header files])
  CXXFLAGS="$CXXFLAGS $TAUCS_INCFLAGS"
fi

#########
# MUMPS #
#########

AC_COIN_HAS_MUMPS

########
# WSMP #
########

AC_ARG_WITH([wsmp],
            AC_HELP_STRING([--with-wsmp],
                           [specify WSMP library]),
            [use_wsmp=$withval], [use_wsmp=no])

if test "$use_wsmp" != "no"; then
  LIBS="$use_wsmp $LIBS"
  ADDLIBS="$use_wsmp $ADDLIBS"

  AC_LANG_PUSH([Fortran 77])
  AC_MSG_CHECKING([whether user-supplied WSMP library \"$use_wsmp\" works])
  AC_TRY_LINK([],[      call WSSMP()],
              [AC_MSG_RESULT(yes)],
              [AC_MSG_ERROR([WSMP library $use_wsmp does not seem to work])])
  AC_LANG_POP([Fortran 77])

  AC_DEFINE(HAVE_WSMP,1,
            [Define to 1 if WSMP is available])
fi

# Check if at least one linear solver was found

AC_ARG_ENABLE([distconf],,
  [case "${enableval}" in
    yes) distconf=true ;;
    no)  distconf=false ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-distconf) ;;
  esac],[distconf=false])

if test "$use_MA27" = "no" -a "$use_MA57" = "no" -a "$use_pardiso" = "no" -a "$use_taucs" = no -a "$use_wsmp" = no -a "$coin_has_mumps" = no; then
  if test "$distconf" = false; then
    AC_MSG_ERROR("No linear solver (such as MUMPS) available.  You need to provide at least one.")
  fi
fi

# For the remaining tests we need to include FLIBS into LIBS, so that
# the C compiler can link programs
#LIBS="$LIBS $FLIBS"
# ToDo should this go into BLAS and LAPACK test????

#######################
# AMPL SOLVER LIBRARY #
#######################

LIBS="$LIBS_without_ADDLIBS"

AC_COIN_HAS_ASL

AC_LANG_POP(C++)

######################################
# Equivalent int Fortran and C types #
######################################

# FIXME: The following test should be active, but this requires chance in
#        code to copy Index* to ipfint* arrays...
if test "$cross_compiling" = no && test "$is_bg" != yes; then
  AC_LANG_PUSH(C)
  AC_DEFINE([FORTRAN_INTEGER_TYPE],[int],[Define to the C type corresponding to Fortran INTEGER])
#  AC_CHECK_SIZEOF([long])
#  AC_CHECK_SIZEOF([int])
#  AC_CHECK_SIZEOF([double])
#  AC_MSG_CHECKING([for C type corresponding to Fortran INTEGER])
#  if test $ac_cv_sizeof_long = $ac_cv_sizeof_double; then
#    AC_DEFINE([FORTRAN_INTEGER_TYPE],[int],[Define to the C type corresponding to Fortran INTEGER])
#    AC_MSG_RESULT([int])
#  else
#    AC_DEFINE([FORTRAN_INTEGER_TYPE],[long])
#    AC_MSG_RESULT([long])
#  fi
  AC_LANG_POP(C)
else
  AC_MSG_WARN([We are cross compiling, assuming Fortran `INTEGER' type corresponds to C `int' type])
  AC_DEFINE([FORTRAN_INTEGER_TYPE],[int])
fi

############# JAVA

if test "${JAVA_HOME+set}" = set; then
  AC_MSG_NOTICE([JAVA_HOME = $JAVA_HOME])
fi
AC_SUBST(JAVA_HOME)
case $build in
  *-cygwin* | *-mingw*)
    JIPOPTDLL=jipopt.dll
    ;;
  *)
    JIPOPTDLL=libjipopt.so
    ;;
esac
AC_SUBST(JIPOPTDLL)

# Compiler flag for linking DLLs/shared objects (not yet portable!)
SHARED_FLAGS=-shared
AC_SUBST(SHARED_FLAGS)

AC_COIN_LINK_FROM_FILELIST(contrib/JavaInterface,contrib/JavaInterface/javafiles.lst)

##########
# ASTYLE #
##########
ASTYLE=astyle
AC_SUBST(ASTYLE)
ASTYLEFLAGS="--mode=c --indent=spaces=2 --indent-cases --indent-namespaces --min-conditional-indent=1 --brackets=linux --brackets=break-closing --convert-tabs"
AC_SUBST(ASTYLEFLAGS)

# The following one is for the example makefiles
case $CXX in
  cl*)
    ipoptlib='`$(CYGPATH_W) $(IPOPTLIBDIR)`\\libipopt.lib';;
  *)
    ipoptlib="-L\$(IPOPTLIBDIR) -lipopt";;
esac
AC_SUBST(ipoptlib)

#########################
# Makefile conditionals #
#########################

# The following variable collects the names of libraries that should
# be included into libipopt.a (relative to subdir Interfaces, where it
# is made)

AC_SUBST(IPALLLIBS)
IPALLLIBS="../contrib/CGPenalty/libcgpenalty.la ../contrib/LinearSolverLoader/libLinearSolverLoader.la ../Algorithm/libipoptalg.la ../Algorithm/LinearSolvers/liblinsolvers.la ../Common/libcommon.la ../LinAlg/liblinalg.la ../LinAlg/TMatrices/libtmatrices.la"

if test "$use_blas" = BUILD; then
  IPALLLIBS="../../../ThirdParty/Blas/libcoinblas.la $IPALLLIBS"
fi

if test "$use_lapack" = BUILD; then
  IPALLLIBS="../../../ThirdParty/Lapack/libcoinlapack.la $IPALLLIBS"
fi

if test "$use_hsl" = BUILD ; then
  IPALLLIBS="../../../ThirdParty/HSL/libcoinhsl.la $IPALLLIBS"
fi

if test "$use_mumps" = BUILD ; then
  IPALLLIBS="../../../ThirdParty/Mumps/libcoinmumps.la $IPALLLIBS"
fi

AM_CONDITIONAL([HAVE_MA27],[test ! "$use_MA27" = no])
AM_CONDITIONAL([HAVE_MA57],[test ! "$use_MA57" = no])
AM_CONDITIONAL([HAVE_MA28],[test ! "$use_MA28" = no])
AM_CONDITIONAL([HAVE_MC19],[test ! "$use_MC19" = no])
AM_CONDITIONAL([HAVE_PARDISO],[test ! "$use_pardiso" = no])
AM_CONDITIONAL([HAVE_TAUCS],[test ! "$use_taucs" = no])
AM_CONDITIONAL([HAVE_WSMP],[test ! "$use_wsmp" = no])

AM_CONDITIONAL([BUILD_HSLLOADER],[test "$use_MA27" = no -o "$use_MA28" = no -o "$use_MA57" = no -o "$use_MC19" = no])
AM_CONDITIONAL([BUILD_PARDISOLOADER],[test "$use_pardiso" = no])

########################################################################
##                 Linear solver loader stuff                         ##
########################################################################

AC_COIN_CHECK_HEADER([windows.h])
AC_COIN_CHECK_HEADER([dlfcn.h])
AC_CHECK_LIB(dl,[dlopen],[ADDLIBS="$ADDLIBS -ldl"; LIBS="$LIBS -ldl"],[])

case $build in
  *-cygwin* |  *-mingw*)
     AC_DEFINE([SHAREDLIBEXT], ["dll"], [Set to extension for shared libraries in quotes.])
     ;;
  *-darwin*)
     AC_DEFINE([SHAREDLIBEXT], ["dylib"], [Set to extension for shared libraries in quotes.])
     ;;
  *)
     AC_DEFINE([SHAREDLIBEXT], ["so"], [Set to extension for shared libraries in quotes.])
esac


########################################################################
##             Create links for the test source files                 ##
########################################################################

AC_CONFIG_LINKS([test/hs071_main.cpp:examples/hs071_cpp/hs071_main.cpp
                 test/hs071_nlp.cpp:examples/hs071_cpp/hs071_nlp.cpp
                 test/hs071_nlp.hpp:examples/hs071_cpp/hs071_nlp.hpp
                 test/hs071_c.c:examples/hs071_c/hs071_c.c])

########################################################################
##          Create links for VPATH config of certain files            ##
########################################################################

AC_COIN_VPATH_LINK(contrib/MatlabInterface/ipopt.m)

AC_CONFIG_FILES([Makefile
		 src/Common/Makefile
                 src/LinAlg/Makefile
		 src/LinAlg/TMatrices/Makefile 
		 src/Interfaces/Makefile
		 src/Interfaces/ipopt_addlibs_f.txt
		 src/Interfaces/ipopt_addlibs_cpp.txt
		 src/Algorithm/Makefile
		 src/Algorithm/LinearSolvers/Makefile
		 src/contrib/CGPenalty/Makefile
		 src/contrib/LinearSolverLoader/Makefile
		 src/Apps/Makefile
		 src/Apps/AmplSolver/Makefile
		 src/Apps/CUTErInterface/Makefile
		 examples/hs071_f/Makefile
		 examples/hs071_f/hs071_f.f
		 examples/Cpp_example/Makefile
		 examples/hs071_cpp/Makefile
		 examples/hs071_c/Makefile
		 examples/ScalableProblems/Makefile
		 test/Makefile
		 test/run_unitTests
		 contrib/MatlabInterface/src/Makefile
		 contrib/JavaInterface/Makefile])

AC_CONFIG_HEADER([inc/config_ipopt.h])

AC_COIN_FINALIZE
