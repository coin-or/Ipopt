# Copyright (C) 2004, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Ipopt],[3.9trunk],[http://projects.coin-or.org/Ipopt/newticket])

AC_COPYRIGHT(
[Copyright 2004, 2010 International Business Machines and others.
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
AC_COIN_PROJECTDIR_INIT(Ipopt)

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

# Also set up the FC compiler variable, required to compile .f90
# files.  We use the same compile and flags as F77 and FFLAGS, unless
# the user has set FC and FCFLAGS explicitly.
if test "$FC" = ""; then
  FC="$F77"
fi
if test "$FCFLAGS" = ""; then
  FCFLAGS="$FFLAGS"
fi
AC_PROG_FC

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

AC_COIN_CHECK_LIBM(IpoptLib)

AC_COIN_CHECK_PACKAGE_BLAS(IpoptLib)
if test $coin_has_blas != yes; then
  AC_MSG_ERROR([Required package BLAS not found.])
fi

AC_COIN_CHECK_PACKAGE_LAPACK(IpoptLib)
if test $coin_has_lapack != yes; then
  AC_MSG_WARN([Compiling code without LAPACK.  Certain options (e.g., quasi-Newton) will not work.])
fi

AC_COIN_CHECK_PACKAGE(ASL,   [coinasl])
AC_COIN_CHECK_PACKAGE(HSL,   [coinhsl],   [IpoptLib])
AC_COIN_CHECK_PACKAGE(Mumps, [coinmumps], [IpoptLib])

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

#######
# HSL #
#######

# check HSL_CFLAGS to see which HSL functions are available in the HSL lib

have_ma27=no
have_ma28=no
have_ma57=no
have_mc19=no
have_ma77=no
if test $coin_has_hsl = yes ; then
  case $HSL_CFLAGS in 
    *HAVE_MA27*)
       have_ma27=yes
       ;;
  esac
  case $HSL_CFLAGS in 
    *HAVE_MA28*)
       have_ma28=yes
       ;;
  esac
  case $HSL_CFLAGS in 
    *HAVE_MA57*)
       have_ma57=yes
       ;;
  esac
  case $HSL_CFLAGS in 
    *HAVE_MC19*)
       have_mc19=yes
       ;;
  esac
  case $HSL_CFLAGS in 
    *HAVE_MA77*)
       have_ma77=yes
       ;;
  esac
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
  SAVE_LIBS="$LIBS"
  PARDISO_LIB="$use_pardiso"
  LIBS="$PARDISO_LIB $LIBS $LAPACK_LIBS $BLAS_LIBS"

# ToDo put this back in
  AC_MSG_CHECKING([whether user-supplied Pardiso library \"$use_pardiso\" works])
  AC_TRY_LINK([],[      call PARDISO()],
              [AC_MSG_RESULT(yes)],
              [AC_MSG_ERROR([Pardiso library $use_pardiso does not seem to work])])

  AC_DEFINE(HAVE_PARDISO,1,[Define to 1 if Pardiso is available])
  
# Test if this is the parallel version of Pardiso
  AC_MSG_CHECKING([whether Pardiso library is the parallel version])
  AC_TRY_LINK([],[      call pardiso_exist_parallel()],
              [AC_MSG_RESULT(yes)
               AC_DEFINE(HAVE_PARDISO_PARALLEL,1,[Define to 1 if you are using the parallel version of Pardiso])],
              [AC_MSG_RESULT(no)])
  
  AC_MSG_CHECKING([whether Pardiso library is at least version 4.0])
  AC_TRY_LINK([],[      call pardiso_ipopt_newinterface()],
              [AC_MSG_RESULT(yes)
               AC_DEFINE(HAVE_PARDISO_NEWINTERFACE,1,[Define to 1 if you are using at least a 4.0 version of Pardiso])],
              [AC_MSG_RESULT(no)])
  LIBS="$SAVE_LIBS"
  IPOPTLIB_LIBS="$PARDISO_LIB $IPOPTLIB_LIBS"
  IPOPTLIB_PCLIBS="$PARDISO_LIB $IPOPTLIB_PCLIBS"
fi
AC_LANG_POP([Fortran 77])

########
# WSMP #
########

AC_ARG_WITH([wsmp],
            AC_HELP_STRING([--with-wsmp],
                           [specify WSMP library]),
            [use_wsmp=$withval], [use_wsmp=no])

if test "$use_wsmp" != "no"; then
  SAVE_LIBS="$LIBS"
  LIBS="$LIBS $use_wsmp $BLAS_LIBS"

  AC_LANG_PUSH([Fortran 77])
  AC_MSG_CHECKING([whether user-supplied WSMP library \"$use_wsmp\" works])
  AC_TRY_LINK([],[      call WSSMP()],
              [AC_MSG_RESULT(yes)],
              [AC_MSG_ERROR([WSMP library $use_wsmp does not seem to work])])
  AC_LANG_POP([Fortran 77])

  AC_DEFINE(HAVE_WSMP,1,
            [Define to 1 if WSMP is available])
  LIBS="$SAVE_LIBS"
  IPOPTLIB_LIBS="$use_wsmp $IPOPTLIB_LIBS"
  IPOPTLIB_PCLIBS="$use_wsmp $IPOPTLIB_PCLIBS"
fi

# Check if at least one linear solver was found

if test "$have_ma27" != "yes" -a "$have_ma57" != "yes" -a "$have_ma77" != "yes" -a "$use_pardiso" = "no" -a "$use_taucs" = no -a "$use_wsmp" = no -a "$coin_has_mumps" != yes; then
  AC_MSG_WARN("You are compiling Ipopt without any linear solver (such as MUMPS).  You will only be able to run Ipopt if you compile a shared library with a linear solver.")
fi

###################################################
# Check if user wants inexact algorithm available #
###################################################

AC_ARG_ENABLE([inexact-solver],
  [AC_HELP_STRING([--enable-inexact-solver],
     [enable inexact linear solver version EXPERIMENTAL! (default: no)])],
  [case "$enableval" in
     no | yes) ;;
     *)
       AC_MSG_ERROR([invalid argument for --enable-inexact-solver: $enableval]);;
   esac
   use_inexact=$enableval],
  [use_inexact=no])

if test $use_inexact = yes; then
  if test "$use_pardiso" = "no"; then
    AC_MSG_ERROR([The inexact solver option is currently only available with Pardiso])
  fi
  # Check if the global function pointer variable is defined in the Pardiso library
  SAVE_LIBS="$LIBS"
  LIBS="$LIBS $PARDISO_LIB $LAPACK_LIBS $BLAS_LIBS $FLIBS"
  AC_LANG_PUSH(C)
  AC_MSG_CHECKING([if Pardiso library has global variable IpoptFunction])
  AC_LINK_IFELSE(
[extern int (*IpoptFunction)(int n, double* x,  double* r, int k, double b);

int main() {
  IpoptFunction = 0;
  return 0;
}],
    [AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)
     AC_MSG_ERROR([Pardiso library not recent enough for the inexact solver option])])
  AC_LANG_POP(C)
  LIBS="$SAVE_LIBS"
fi


# For the remaining tests we need to include FLIBS into LIBS, so that
# the C compiler can link programs
#LIBS="$LIBS $FLIBS"
# ToDo should this go into BLAS and LAPACK test????

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

#doing this check here avoids its repetition for each file in the following loop
AC_COIN_CHECK_VPATH
for file in `cat $srcdir/contrib/JavaInterface/javafiles.lst`
do
  AC_COIN_VPATH_LINK(contrib/JavaInterface/$file)
done

############ Matlab

# Try to detect if there is matlab installed on the system and set the
# MATLAB_HOME variable

AC_MSG_CHECKING([for matlab installation])

AC_ARG_WITH([matlab-home],
            AC_HELP_STRING([--with-matlab-home],
                           [specify location of matlab installation for the matlab Ipopt interface.  The "mex" compiler should be in MATLAB_HOME/bin]),
            [use_matlab_home=$withval], [use_matlab_home=])


if test x"$use_matlab_home" = x; then
  # Try to detect matlab home automatically
  mexpath=`which mexext 2>/dev/null`
  if test -z "$mexpath" ||
     expr "$mexpath" : '.*no mexext.*' >/dev/null 2>&1 ; then
    AC_MSG_RESULT([no])
  else
    MATLAB_HOME=`echo $mexpath | sed -e 's=/bin/mexext=='`
    AC_MSG_RESULT([$MATLAB_HOME])
  fi
else
  MATLAB_HOME="$use_matlab_home"
  AC_MSG_RESULT([$MATLAB_HOME])
  if test ! -x $MATLAB_HOME/bin/mex; then
    AC_MSG_ERROR([Matlab compiler $MATLAB_HOME/bin/mex not found or not executable])
  fi
fi

MEXSUFFIX=
if test x"$MATLAB_HOME" != x; then
  AC_MSG_CHECKING([for Matlab executable extension (MEXSUFFIX)])

##   cat >conftest.c <<_ACEOF
## #include "mex.h"
## void mexFunction(int nlhs, mxArray *plhs@<:@@:>@,
##   int nrhs, const mxArray *prhs@<:@@:>@) {
##   int i=0;
## }
## _ACEOF
##   if $MATLAB_HOME/bin/mex conftest.c >5 2>&5; then
##    MEXSUFFIX=`ls conftest.mex* 5>/dev/null | sed -e 's=conftest.=='`
##    AC_MSG_RESULT([$MEXSUFFIX])
##  else
##    AC_MSG_RESULT([error running MEX compiler])
##  fi

  if ! test -x $MATLAB_HOME/bin/mexext; then
    AC_MSG_WARN([Cannot find $MATLAB_HOME/bin/mexext executable to determine MEXSUFFIX])
  else
    MEXSUFFIX=`$MATLAB_HOME/bin/mexext`
    AC_MSG_RESULT([$MEXSUFFIX])
  fi
fi

AC_SUBST(MATLAB_HOME)
AC_SUBST(MEXSUFFIX)

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
IPALLLIBS="../contrib/CGPenalty/libcgpenalty.la ../Algorithm/libipoptalg.la ../Algorithm/LinearSolvers/liblinsolvers.la ../Common/libcommon.la ../LinAlg/liblinalg.la ../LinAlg/TMatrices/libtmatrices.la"

AM_CONDITIONAL([HAVE_MA27],[test "$have_ma27" = yes])
AM_CONDITIONAL([HAVE_MA57],[test "$have_ma57" = yes])
AM_CONDITIONAL([HAVE_MA28],[test "$have_ma28" = yes])
AM_CONDITIONAL([HAVE_MA77],[test "$have_ma77" = yes])
AM_CONDITIONAL([HAVE_MC19],[test "$have_mc19" = yes])
AM_CONDITIONAL([HAVE_PARDISO],[test ! "$use_pardiso" = no])
AM_CONDITIONAL([HAVE_WSMP],[test ! "$use_wsmp" = no])

AC_ARG_ENABLE([linear-solver-loader],
  [AC_HELP_STRING([--enable-linear-solver-loader],
     [compile linear solver loader (default: yes)])],
  [case "$enableval" in
     no | yes) ;;
     *)
       AC_MSG_ERROR([invalid argument for --enable-linear-solver-loader: $enableval]);;
   esac
   use_linearsolverloader=$enableval],
  [use_linearsolverloader=yes])

AC_MSG_CHECKING([whether the linear solver loader should be compiled])
AC_MSG_RESULT([$use_linearsolverloader])

if test $use_linearsolverloader = yes; then
  AC_DEFINE([HAVE_LINEARSOLVERLOADER],[1],
            [Define to 1 if the linear solver loader should be compiled to allow dynamic loading of shared libaries with linear solvers])
  IPALLLIBS="../contrib/LinearSolverLoader/libLinearSolverLoader.la $IPALLLIBS"
fi

AM_CONDITIONAL([BUILD_LINEARSOLVERLOADER],[test $use_linearsolverloader = yes])
AM_CONDITIONAL([BUILD_HSLLOADER],[test "$have_ma27" != yes -o "$have_ma28" != yes -o "$have_ma57" != yes -o "$have_ma77" != yes -o "$have_mc19" != yes])
AM_CONDITIONAL([BUILD_PARDISOLOADER],[test "$use_pardiso" = no])

AM_CONDITIONAL([BUILD_INEXACT], [test $use_inexact = yes])
if test $use_inexact = yes; then
  IPALLLIBS="../Algorithm/Inexact/libinexact.la $IPALLLIBS"
  AC_DEFINE([BUILD_INEXACT],[1],
            [Define to 1 if the inexact linear solver option is included])
fi

########################################################################
##                 Linear solver loader stuff                         ##
########################################################################

AC_LANG_PUSH(C)
AC_COIN_CHECK_HEADER([windows.h])
AC_COIN_CHECK_HEADER([dlfcn.h])
AC_CHECK_LIB(dl,[dlopen],[IPOPTLIB_LIBS="$IPOPTLIB_LIBS -ldl"; IPOPTLIB_PCLIBS="$IPOPTLIB_PCLIBS -ldl"],[])

case $build in
  *-cygwin* |  *-mingw*)
     AC_DEFINE([SHAREDLIBEXT], ["dll"], [Set to extension for shared libraries in quotes.])
     # Take -pedantic-errors out of CFLAGS
     AC_MSG_NOTICE([Removing flag -pedantic-errors from CFLAGS])
     CFLAGS=`echo $CFLAGS | sed -e s/-pedantic-errors//`
     ;;
  *-darwin*)
     AC_DEFINE([SHAREDLIBEXT], ["dylib"], [Set to extension for shared libraries in quotes.])
     ;;
  *)
     AC_DEFINE([SHAREDLIBEXT], ["so"], [Set to extension for shared libraries in quotes.])
esac
AC_CHECK_FUNCS([snprintf _snprintf],[break])
AC_CHECK_FUNCS([vsnprintf _vsnprintf],[break])
AC_LANG_POP(C)

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

AC_COIN_VPATH_LINK(contrib/MatlabInterface/examples/examplehs038.m)
AC_COIN_VPATH_LINK(contrib/MatlabInterface/examples/examplehs051.m)
AC_COIN_VPATH_LINK(contrib/MatlabInterface/examples/examplehs071.m)
AC_COIN_VPATH_LINK(contrib/MatlabInterface/examples/examplelasso.m)
AC_COIN_VPATH_LINK(contrib/MatlabInterface/examples/lasso.m)

AC_COIN_VPATH_LINK(examples/Cpp_example/cpp_example.cpp)
AC_COIN_VPATH_LINK(examples/Cpp_example/MyNLP.cpp)
AC_COIN_VPATH_LINK(examples/Cpp_example/MyNLP.hpp)
AC_COIN_VPATH_LINK(examples/hs071_cpp/hs071_main.cpp)
AC_COIN_VPATH_LINK(examples/hs071_cpp/hs071_nlp.cpp)
AC_COIN_VPATH_LINK(examples/hs071_cpp/hs071_nlp.hpp)
AC_COIN_VPATH_LINK(examples/hs071_c/hs071_c.c)

AC_COIN_VPATH_LINK(tutorial/AmplExperiments/hs71.mod)
AC_COIN_VPATH_LINK(tutorial/AmplExperiments/infeasible.mod)
AC_COIN_VPATH_LINK(tutorial/AmplExperiments/MoreAmplModels.txt)
AC_COIN_VPATH_LINK(tutorial/AmplExperiments/car1.run)
AC_COIN_VPATH_LINK(tutorial/AmplExperiments/car1.gp)

AC_COIN_VPATH_LINK(tutorial/Modeling/bad1.mod)
AC_COIN_VPATH_LINK(tutorial/Modeling/bad1-fix1.mod)
AC_COIN_VPATH_LINK(tutorial/Modeling/bad1-fix2.mod)

AC_COIN_VPATH_LINK(tutorial/CodingExercise/exercise_example.mod)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/C/1-skeleton/TutorialC.c)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/C/2-mistake/TutorialC.c)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/C/3-solution/TutorialC.c)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Cpp/1-skeleton/TutorialCpp_main.cpp)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Cpp/1-skeleton/TutorialCpp_nlp.hpp)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Cpp/1-skeleton/TutorialCpp_nlp.cpp)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Cpp/2-mistake/TutorialCpp_main.cpp)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Cpp/2-mistake/TutorialCpp_nlp.hpp)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Cpp/2-mistake/TutorialCpp_nlp.cpp)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Cpp/3-solution/TutorialCpp_main.cpp)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Cpp/3-solution/TutorialCpp_nlp.hpp)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Cpp/3-solution/TutorialCpp_nlp.cpp)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Matlab/1-skeleton/TutorialMatlab.m)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Matlab/2-mistake/TutorialMatlab.m)
AC_COIN_VPATH_LINK(tutorial/CodingExercise/Matlab/3-solution/TutorialMatlab.m)

AC_COIN_VPATH_LINK(contrib/AsNMPC/examples/parametric_cpp/parametric_driver.cpp)
AC_COIN_VPATH_LINK(contrib/AsNMPC/examples/parametric_cpp/parametricTNLP.cpp)
AC_COIN_VPATH_LINK(contrib/AsNMPC/examples/parametric_cpp/parametricTNLP.hpp)
AC_COIN_VPATH_LINK(contrib/AsNMPC/examples/redhess_cpp/MyNmpcTNLP.cpp)
AC_COIN_VPATH_LINK(contrib/AsNMPC/examples/redhess_cpp/MyNmpcTNLP.hpp)
AC_COIN_VPATH_LINK(contrib/AsNMPC/examples/redhess_cpp/redhess_cpp.cpp)

########################################################################
##          Create Makefiles and other stuff                          ##
########################################################################

AC_CONFIG_FILES([Makefile
		 src/Common/Makefile
		 src/LinAlg/Makefile
		 src/LinAlg/TMatrices/Makefile 
		 src/Interfaces/Makefile
		 src/Interfaces/ipopt_addlibs_f.txt
		 src/Algorithm/Makefile
		 src/Algorithm/LinearSolvers/Makefile
		 src/Algorithm/Inexact/Makefile
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
		 tutorial/CodingExercise/C/1-skeleton/Makefile
		 tutorial/CodingExercise/C/2-mistake/Makefile
		 tutorial/CodingExercise/C/3-solution/Makefile
		 tutorial/CodingExercise/Cpp/1-skeleton/Makefile
		 tutorial/CodingExercise/Cpp/2-mistake/Makefile
		 tutorial/CodingExercise/Cpp/3-solution/Makefile
		 tutorial/CodingExercise/Matlab/1-skeleton/startup.m
		 tutorial/CodingExercise/Matlab/2-mistake/startup.m
		 tutorial/CodingExercise/Matlab/3-solution/startup.m
		 tutorial/CodingExercise/Fortran/1-skeleton/TutorialFortran.f
		 tutorial/CodingExercise/Fortran/2-mistake/TutorialFortran.f
		 tutorial/CodingExercise/Fortran/3-solution/TutorialFortran.f
		 tutorial/CodingExercise/Fortran/1-skeleton/Makefile
		 tutorial/CodingExercise/Fortran/2-mistake/Makefile
		 tutorial/CodingExercise/Fortran/3-solution/Makefile
		 test/Makefile
		 test/run_unitTests
		 contrib/MatlabInterface/src/Makefile
		 contrib/MatlabInterface/examples/startup.m
		 contrib/MatlabInterface-Old/src/Makefile
		 contrib/MatlabInterface-Old/examples/startup.m
		 contrib/JavaInterface/Makefile
		 contrib/AsNMPC/Makefile
		 contrib/AsNMPC/AmplSolver/Makefile
		 contrib/AsNMPC/examples/parametric_cpp/Makefile
		 contrib/AsNMPC/examples/redhess_cpp/Makefile
		 contrib/AsNMPC/src/Makefile
		 ipopt.pc ipopt-uninstalled.pc])

AC_CONFIG_HEADER([inc/config_ipopt.h])

AC_COIN_FINALIZE
