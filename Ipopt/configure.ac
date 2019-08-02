# Copyright (C) 2004, 2011 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Ipopt],[trunk],[https://github.com/coin-or/Ipopt/issues/new],,[https://github.com/coin-or/Ipopt])

AC_COPYRIGHT(
[Copyright 2004, 2011 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package IPOPT which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/Common/IpDebug.hpp)

# Do some project-level initialization work (version numbers, change prefix
# default, ...)
AC_COIN_INITIALIZE(Ipopt)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C and C++ compilers and appropriate compiler options.
AC_COIN_PROG_CC
AC_COIN_PROG_CXX
AM_CONDITIONAL([COIN_HAS_F77],false)

# This is a C++ package, set the language accordingly.
#AC_LANG_PUSH(C++)

# Initialize libtool
AC_COIN_PROG_LIBTOOL



# We need egrep and ln -s somewherw
#AC_PROG_EGREP
#AC_PROG_LN_S

# Get the C++ runtime libraries in case we want to link with a C or
# Fortran compiler
#AC_COIN_CXXLIBS

# Get the name of the Fortran compiler and appropriate compiler options
#AC_COIN_PROG_F77
#AM_CONDITIONAL(COIN_HAS_F77, test "$F77" != "unavailable")

# Find out how to call Fortran from C and determine Fortran runtime libraries
#if test "$F77" != "unavailable" ; then
#  AC_COIN_F77_WRAPPERS
#fi

# Initialize automake and libtool
#AC_COIN_INIT_AUTO_TOOLS

# A Fix for the AIX compiler
#AM_CONDITIONAL([AIX_FIX_COMPILER_BUG], test "x$AIX_FIX_COMPILER_BUG" != x)

# Doxygen
AC_COIN_DOXYGEN


#############################################################################
#                              COIN-OR components                           #
#############################################################################

AC_COIN_CHK_LIBM(IpoptLib)

#AC_COIN_CHECK_PACKAGE_BLAS(IpoptLib)
#if test $coin_has_blas != yes; then
#  AC_MSG_ERROR([Required package BLAS not found.])
#fi

AC_COIN_CHK_LAPACK(IpoptLib)
#if test $coin_has_lapack != yes; then
#  AC_MSG_WARN([Compiling code without LAPACK.  Certain options (e.g., quasi-Newton) will not work.])
#fi

AC_COIN_CHK_PKG(ASL,[IpoptAmplInterfaceLib],[coinasl],[build])
#AC_COIN_CHECK_PACKAGE(HSL,   [coinhsl >= 1.3],   [IpoptLib])
#AC_COIN_CHECK_PACKAGE(Mumps, [coinmumps], [IpoptLib])

AM_CONDITIONAL([COIN_HAS_HSL],[false])
AM_CONDITIONAL([COIN_HAS_MUMPS],[false])
AM_CONDITIONAL([HAVE_PARDISO],[false])
AM_CONDITIONAL([HAVE_MA28],[false])
AM_CONDITIONAL([HAVE_WSMP],[false])


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

############################################################################
############################################################################
#                   Stuff that we need for C++ programs                    #
############################################################################
############################################################################

AC_LANG_PUSH(C++)

################
# Header Files #
################

#AC_COIN_CHECK_CXX_CHEADER(math)
#AC_COIN_CHECK_CXX_CHEADER(ctype)
#AC_COIN_CHECK_CXX_CHEADER(assert)
#AC_COIN_CHECK_CXX_CHEADER(stdio)
#AC_COIN_CHECK_CXX_CHEADER(stdlib)
#AC_COIN_CHECK_CXX_CHEADER(stdarg)
#AC_COIN_CHECK_CXX_CHEADER(string)
#AC_COIN_CHECK_CXX_CHEADER(float)
#AC_COIN_CHECK_CXX_CHEADER(time)
#AC_COIN_CHECK_CXX_CHEADER(stddef)

###################
# Function finite #
###################

#AC_COIN_CHECK_ISFINITE

###########
# va_copy #
###########

AC_CHECK_DECL([va_copy],[AC_DEFINE([HAVE_VA_COPY],[1],
                                   [Define to 1 if va_copy is available])],,[#include <cstdarg>])

###########################
# Random number generator #
###########################

AC_CHECK_DECL([drand48],[AC_DEFINE([HAVE_DRAND48],[1],
                                   [Define to 1 if function drand48 is available])],,[#include <cstdlib>])

AC_CHECK_DECL([rand],[AC_DEFINE([HAVE_RAND],[1],
                                [Define to 1 if function rand is available])],,[#include <cstdlib>])

AC_CHECK_DECL([srand],[AC_DEFINE([HAVE_STD__RAND],[1],
                                 [Define to 1 if function std::rand is available])],,[
#include <cstdlib>
using namespace std;])

##########################################################################

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
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [extern int (*IpoptFunction)(int n, double* x,  double* r, int k, double b);],
      [IpoptFunction = 0;])],
    [AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)
     AC_MSG_ERROR([Pardiso library not recent enough for the inexact solver option])])
  AC_LANG_POP(C)
  LIBS="$SAVE_LIBS"
fi

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
  AC_MSG_WARN([We are cross compiling, assuming Fortran 'INTEGER' type corresponds to C 'int' type])
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
  *-darwin* )
    JIPOPTDLL=libjipopt.dylib
    ;;  
  *)
    JIPOPTDLL=libjipopt.so
    ;;
esac
AC_SUBST(JIPOPTDLL)

# Compiler flag for linking DLLs/shared objects (not yet portable!)
#if test $coin_cxx_is_cl = true; then
#  SHARED_FLAGS=-LD
#else 
#  SHARED_FLAGS=-shared
#fi
#AC_SUBST(SHARED_FLAGS)

#doing this check here avoids its repetition for each file in the following loop
AC_COIN_CHECK_VPATH
for file in `cat $srcdir/contrib/JavaInterface/javafiles.lst`
do
  AC_COIN_VPATH_LINK(contrib/JavaInterface/$file)
done


##########
# ASTYLE #
##########

ASTYLE=astyle
AC_SUBST(ASTYLE)
ASTYLEFLAGS="--mode=c -A1 --indent=spaces=3 --indent-switches --min-conditional-indent=1 --convert-tabs --align-pointer=type --pad-oper --add-brackets"
AC_SUBST(ASTYLEFLAGS)


#########################
# Makefile conditionals #
#########################

# The following variable collects the names of libraries that should
# be included into libipopt.a (relative to subdir Interfaces, where it
# is made)

AC_SUBST(IPALLLIBS)
IPALLLIBS="../contrib/CGPenalty/libcgpenalty.la ../Algorithm/libipoptalg.la ../Algorithm/LinearSolvers/liblinsolvers.la ../Common/libcommon.la ../LinAlg/liblinalg.la ../LinAlg/TMatrices/libtmatrices.la"

#AM_CONDITIONAL([HAVE_PARDISO],[test ! "$use_pardiso" = no])
#AM_CONDITIONAL([HAVE_WSMP],[test ! "$use_wsmp" = no])

AC_ARG_ENABLE([linear-solver-loader],
  [AC_HELP_STRING([--enable-linear-solver-loader],
     [compile linear solver loader (default: yes, if F77 available, otherwise no)])],
  [case "$enableval" in
     no | yes) ;;
     *)
       AC_MSG_ERROR([invalid argument for --enable-linear-solver-loader: $enableval]);;
   esac
   use_linearsolverloader=$enableval],
  [if test "$F77" != unavailable ; then use_linearsolverloader=yes; else use_linearsolverloader=no; fi])

AC_MSG_CHECKING([whether the linear solver loader should be compiled])
AC_MSG_RESULT([$use_linearsolverloader])

if test $use_linearsolverloader = yes; then
  AC_DEFINE([HAVE_LINEARSOLVERLOADER],[1],
            [Define to 1 if the linear solver loader should be compiled to allow dynamic loading of shared libraries with linear solvers])
  IPALLLIBS="../contrib/LinearSolverLoader/libLinearSolverLoader.la $IPALLLIBS"
fi

AM_CONDITIONAL([BUILD_LINEARSOLVERLOADER],[test $use_linearsolverloader = yes])
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
AC_CHECK_HEADER([windows.h])
AC_CHECK_HEADER([dlfcn.h])
AC_CHECK_LIB(dl,[dlopen],[
  IPOPTLIB_LIBS="$IPOPTLIB_LIBS -ldl"
  IPOPTLIB_PCLIBS="$IPOPTLIB_PCLIBS -ldl"
  IPOPTLIB_LIBS_INSTALLED="$IPOPTLIB_LIBS_INSTALLED -ldl"
 ],[])

case $host in
  *-cygwin* |  *-mingw*)
     AC_DEFINE([SHAREDLIBEXT], ["dll"], [Set to extension for shared libraries in quotes.])
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

AC_COIN_VPATH_LINK(contrib/RInterface/DESCRIPTION)
AC_COIN_VPATH_LINK(contrib/RInterface/NAMESPACE)
AC_COIN_VPATH_LINK(contrib/RInterface/inst/CITATION)
AC_COIN_VPATH_LINK(contrib/RInterface/inst/doc/ipoptr.Rnw)
AC_COIN_VPATH_LINK(contrib/RInterface/inst/doc/ipoptr.pdf)
AC_COIN_VPATH_LINK(contrib/RInterface/inst/doc/reflist.bib)
AC_COIN_VPATH_LINK(contrib/RInterface/man/ipoptr-package.Rd)
AC_COIN_VPATH_LINK(contrib/RInterface/man/ipoptr.Rd)
AC_COIN_VPATH_LINK(contrib/RInterface/man/is.ipoptr.Rd)
AC_COIN_VPATH_LINK(contrib/RInterface/man/make.sparse.Rd)
AC_COIN_VPATH_LINK(contrib/RInterface/man/plot.sparseness.Rd)
AC_COIN_VPATH_LINK(contrib/RInterface/man/print.ipoptr.Rd)
AC_COIN_VPATH_LINK(contrib/RInterface/man/print.sparseness.Rd)
AC_COIN_VPATH_LINK(contrib/RInterface/R/get.option.types.R)
AC_COIN_VPATH_LINK(contrib/RInterface/R/ipoptr.R)
AC_COIN_VPATH_LINK(contrib/RInterface/R/is.ipoptr.R)
AC_COIN_VPATH_LINK(contrib/RInterface/R/make.sparse.R)
AC_COIN_VPATH_LINK(contrib/RInterface/R/plot.sparseness.R)
AC_COIN_VPATH_LINK(contrib/RInterface/R/print.ipoptr.R)
AC_COIN_VPATH_LINK(contrib/RInterface/R/print.sparseness.R)
AC_COIN_VPATH_LINK(contrib/RInterface/tests/approx_banana.R)
AC_COIN_VPATH_LINK(contrib/RInterface/tests/banana.R)
AC_COIN_VPATH_LINK(contrib/RInterface/tests/hs071_nlp.R)
AC_COIN_VPATH_LINK(contrib/RInterface/tests/lasso.R)
AC_COIN_VPATH_LINK(contrib/RInterface/tests/mynlp.R)
AC_COIN_VPATH_LINK(contrib/RInterface/tests/parameters.R)
AC_COIN_VPATH_LINK(contrib/RInterface/tests/sparseness.R)
AC_COIN_VPATH_LINK(contrib/RInterface/src/ipoptr.cpp)
AC_COIN_VPATH_LINK(contrib/RInterface/src/IpoptRJournal.cpp)
AC_COIN_VPATH_LINK(contrib/RInterface/src/IpoptRNLP.cpp)

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

AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_cpp/parametric_driver.cpp)
AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_cpp/parametricTNLP.cpp)
AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_cpp/parametricTNLP.hpp)
AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_dsdp_cpp/parametric_dsdp_driver.cpp)
AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_dsdp_cpp/parametricTNLP.cpp)
AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_dsdp_cpp/parametricTNLP.hpp)
AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/redhess_cpp/MySensTNLP.cpp)
AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/redhess_cpp/MySensTNLP.hpp)
AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/redhess_cpp/redhess_cpp.cpp)

########################################################################
##          Create Makefiles and other stuff                          ##
########################################################################

AC_CONFIG_FILES([Makefile
		 src/Common/Makefile
		 src/LinAlg/Makefile
		 src/LinAlg/TMatrices/Makefile 
		 src/Interfaces/Makefile
		 src/Algorithm/Makefile
		 src/Algorithm/LinearSolvers/Makefile
		 src/Algorithm/Inexact/Makefile
		 src/contrib/CGPenalty/Makefile
		 src/contrib/LinearSolverLoader/Makefile
		 src/Apps/Makefile
		 src/Apps/AmplSolver/Makefile
		 test/Makefile
		 test/run_unitTests
		 ipopt.pc ipopt-uninstalled.pc
		 doc/Doxyfile])
#		 src/Apps/CUTErInterface/Makefile
#		 examples/hs071_f/Makefile
#		 examples/hs071_f/hs071_f.f
#		 examples/Cpp_example/Makefile
#		 examples/hs071_cpp/Makefile
#		 examples/hs071_c/Makefile
#		 examples/ScalableProblems/Makefile
#		 tutorial/CodingExercise/C/1-skeleton/Makefile
#		 tutorial/CodingExercise/C/2-mistake/Makefile
#		 tutorial/CodingExercise/C/3-solution/Makefile
#		 tutorial/CodingExercise/Cpp/1-skeleton/Makefile
#		 tutorial/CodingExercise/Cpp/2-mistake/Makefile
#		 tutorial/CodingExercise/Cpp/3-solution/Makefile
#		 tutorial/CodingExercise/Matlab/1-skeleton/startup.m
#		 tutorial/CodingExercise/Matlab/2-mistake/startup.m
#		 tutorial/CodingExercise/Matlab/3-solution/startup.m
#		 tutorial/CodingExercise/Fortran/1-skeleton/TutorialFortran.f
#		 tutorial/CodingExercise/Fortran/2-mistake/TutorialFortran.f
#		 tutorial/CodingExercise/Fortran/3-solution/TutorialFortran.f
#		 tutorial/CodingExercise/Fortran/1-skeleton/Makefile
#		 tutorial/CodingExercise/Fortran/2-mistake/Makefile
#		 tutorial/CodingExercise/Fortran/3-solution/Makefile
#		 contrib/MatlabInterface/src/Makefile
#		 contrib/MatlabInterface/examples/startup.m
#		 contrib/JavaInterface/Makefile
#		 contrib/sIPOPT/Makefile
#		 contrib/sIPOPT/AmplSolver/Makefile
#		 contrib/sIPOPT/examples/parametric_cpp/Makefile
#		 contrib/sIPOPT/examples/parametric_dsdp_cpp/Makefile	
#		 contrib/sIPOPT/examples/redhess_cpp/Makefile
#		 contrib/sIPOPT/src/Makefile

# under Windows, the Makevars file for the R Interface need to be named Makevars.win
case $build in 
  *-cygwin* | *-mingw*)
    AC_CONFIG_FILES([contrib/RInterface/src/Makevars.win:contrib/RInterface/src/Makevars.in])
    ;;
  *)
    AC_CONFIG_FILES([contrib/RInterface/src/Makevars])
    ;;
esac

if test $coin_has_asl = yes ; then
  AC_CONFIG_FILES([ipoptamplinterface.pc:src/Apps/AmplSolver/ipoptamplinterface.pc.in
                   ipoptamplinterface-uninstalled.pc:src/Apps/AmplSolver/ipoptamplinterface-uninstalled.pc.in
                  ])
fi

AC_CONFIG_HEADER([src/Common/config.h src/Common/config_ipopt.h examples/ScalableProblems/config.h])

AC_COIN_FINALIZE
