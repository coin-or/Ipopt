# Copyright (C) 2004, 2011 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Ipopt],[3.13.2],[https://github.com/coin-or/Ipopt/issues/new],,[https://github.com/coin-or/Ipopt])

AC_COPYRIGHT(
[Copyright 2004, 2011 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package IPOPT which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/Common/IpDebug.hpp)

# Do some project-level initialization work (version numbers, ...)
AC_COIN_INITIALIZE(Ipopt)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C, C++, and Fortran compilers and appropriate compiler options.
AC_COIN_PROG_CC
AC_COIN_PROG_CXX
AC_COIN_PROG_F77

# If there is a Fortran compiler, then setup everything to use it, including F77_FUNC
if test -n "$F77" ; then
  AC_COIN_F77_SETUP
fi

# This is a C++ package, set the language accordingly.
#AC_LANG_PUSH(C++)

# Initialize libtool
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
AC_COIN_RPATH_FLAGS([$expanded_libdir])

# Get the C++ runtime libraries in case we want to link a static Ipopt library
# with a C or Fortran compiler
AC_COIN_CXXLIBS

# Doxygen
AC_COIN_DOXYGEN

# COIN_IPOPT_VERBOSITY and COIN_IPOPT_DEBUGLEVEL
AC_COIN_DEBUGLEVEL(IPOPT)

#############################################################################
#                                 Dependencies                              #
#############################################################################

AC_COIN_CHK_LIBM(IpoptLib)

#AC_COIN_CHK_BLAS(IpoptLib)
#if test $coin_has_blas != yes; then
#  AC_MSG_ERROR([Required package BLAS not found.])
#fi

AC_COIN_CHK_LAPACK(IpoptLib)
if test $coin_has_lapack != yes; then
  AC_MSG_ERROR([Required package LAPACK not found.])
#  AC_MSG_WARN([Compiling code without LAPACK.  Certain options (e.g., quasi-Newton) will not work.])
fi

AC_COIN_CHK_PKG(ASL,[IpoptAmplInterfaceLib],[coinasl],[build])

#########
# MUMPS #
#########

AC_COIN_CHK_PKG(Mumps,[IpoptLib],[coinmumps],[build])

# Check whether MPI_Initialized is available
#   we assume that MPI_Finalized is present if MPI_Initialized is present
AC_CHECK_FUNCS([MPI_Initialized])

#######
# HSL #
#######

AC_COIN_CHK_PKG(HSL,[IpoptLib HSLLib],[coinhsl],[build])

have_ma28=no
if test "$coin_has_hsl" = yes ; then
  AC_COIN_FINALIZE_FLAGS([HSLLib])
  AC_COIN_NAMEMANGLING([HSL],[ma27ad],[$HSLLIB_LFLAGS])

  # extra check for MA28, since that decides whether we have to build IpMa28Partition.F
  if test -n "$F77" ; then
    AC_COIN_TRY_LINK([ma28ad],[$HSLLIB_LFLAGS],,[
      if test "$ma28ad_namemangling" != "$ac_cv_f77_mangling" ; then
        AC_MSG_WARN([Name mangling of MA28 different than Fortran. This will not link. Disabling MA28])
      else
        have_ma28=yes
      fi
    ])
  fi
fi
AM_CONDITIONAL([HAVE_MA28],[test "$have_ma28" = yes])

###########
# PARDISO #
###########

AC_ARG_WITH([pardiso],
            AC_HELP_STRING([--with-pardiso],[specify Pardiso library (>= 4.0) from pardiso-project.org]),
            [have_pardiso_project=yes; pardiso_lflags=$withval],
            [have_pardiso_project=no])

have_pardiso_mkl=no
if test "$have_pardiso_project" = yes ; then
  # check whether flags from --with-pardiso work and figure out name mangling
  # if so, define PARDISO_FUNC and keep lflags
  AC_COIN_TRY_LINK([pardiso_ipopt_newinterface],[$pardiso_lflags $lapack_lflags],[$lapack_pcfiles],
    [AC_COIN_DEFINENAMEMANGLING([PARDISO],[$pardiso_ipopt_newinterface_namemangling])
     IPOPTLIB_LFLAGS="$pardiso_lflags $IPOPTLIB_LFLAGS"
     AC_DEFINE(HAVE_PARDISO,1,[Define to 1 if Pardiso is available])
    ],
    [AC_MSG_ERROR([Symbol pardiso_ipopt_newinterface not found with Pardiso flags $pardiso_lflags and Lapack. Require Pardiso >= 4.0.])])

  AC_COIN_TRY_LINK([pardiso_exist_parallel],[$pardiso_lflags $lapack_lflags],[$lapack_pcfiles],
    [AC_DEFINE(HAVE_PARDISO_PARALLEL,1,[Define to 1 if you are using the parallel version of Pardiso])])

else
  # check whether Pardiso is available via Lapack, which should then be MKL
  # figure out name mangling and define PARDISO_FUNC
  AC_COIN_TRY_LINK([pardiso],[$lapack_lflags],[$lapack_pcfiles],
    [AC_COIN_DEFINENAMEMANGLING([PARDISO],[$pardiso_namemangling])
     have_pardiso_mkl=yes
     AC_DEFINE(HAVE_PARDISO,1,[Define to 1 if Pardiso is available])
     # assume MKL Pardiso is parallel (it never has pardiso_exist_parallel)
     # TODO does sequential MKL also have a parallel Pardiso? do we need some check here?
     AC_DEFINE(HAVE_PARDISO_PARALLEL,1,[Define to 1 if you are using the parallel version of Pardiso])
     AC_DEFINE(HAVE_PARDISO_MKL,1,[Define to 1 if you are using Pardiso from MKL])
    ])
fi

AM_CONDITIONAL([HAVE_PARDISO],[test "$have_pardiso_mkl$have_pardiso_project" != nono])

########
# WSMP #
########

AC_ARG_WITH([wsmp],
            AC_HELP_STRING([--with-wsmp],[specify WSMP library]),
            [have_wsmp=yes; wsmp_lflags=$withval],
            [have_wsmp=no])

if test "$have_wsmp" = "yes"; then
  AC_COIN_TRY_LINK([wssmp],[$wsmp_lflags],[],
    [AC_COIN_DEFINENAMEMANGLING([WSMP],[$wssmp_namemangling])
     IPOPTLIB_LFLAGS="$wsmp_lflags $IPOPTLIB_LFLAGS"
     AC_DEFINE(HAVE_WSMP,1,[Define to 1 if WSMP is available])
    ],
    [AC_MSG_ERROR([Symbol wssmp not found with WSMP flags $wsmp_lflags.])])
fi

AM_CONDITIONAL([HAVE_WSMP],[test $have_wsmp = yes])

#############################################################################
#                             Stuff for examples                            #
#############################################################################

# find out how long an int pointer is to know if we need INTEGER*4 or
# INTEGER*8 in Fortran to capture pointers.
AC_LANG_PUSH(C)
AC_CHECK_SIZEOF(int *)
AC_LANG_POP(C)
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

#####################
# Function isfinite #
#####################

AC_COIN_CHECK_ISFINITE

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

AC_COIN_CHECK_NAMESPACE_DECL([std::srand],[1],
  [AC_DEFINE([HAVE_STD__RAND],[1],[Define to 1 if function std::rand is available])],[],
  [#include <cstdlib>])

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
  if test $have_pardiso_project = no; then
    AC_MSG_ERROR([The inexact solver option is currently only available with Pardiso from pardiso-project.org])
  fi
  AC_DEFINE([BUILD_INEXACT],[1],[Define to 1 if the inexact linear solver option is included])
fi
AM_CONDITIONAL([BUILD_INEXACT], [test $use_inexact = yes])

AC_LANG_POP(C++)

######################################
# Equivalent int Fortran and C types #
######################################

# FIXME: The following test should be active, but this requires change in
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

AC_ARG_ENABLE([java],
  [AC_HELP_STRING([--disable-java],[disable building of Java interface])],
  [enable_java="$enableval"],
  [case "$JAVA_HOME" in
     *\ * ) enable_java=no ;;   # do not enable java-check by default, if there are spaces in JAVA_HOME - that causes trouble
     * )    enable_java="$enable_shared" ;;
   esac
  ])

if test "$enable_java" != no ; then
  # look for javac: required to compile Java code and build C-header
  # this is a modified version of AX_PROG_JAVAC
  m4_define([m4_ax_prog_javac_list],["gcj -C" guavac jikes javac])dnl
  AS_IF([test "x$JAVAPREFIX" = x],
        [test "x$JAVAC" = x && AC_CHECK_PROGS([JAVAC], [m4_ax_prog_javac_list])],
        [test "x$JAVAC" = x && AC_CHECK_PROGS([JAVAC], [m4_ax_prog_javac_list], [], [$JAVAPREFIX/bin])])
  m4_undefine([m4_ax_prog_javac_list])dnl

  if test -z "$JAVAC" ; then
    AC_MSG_NOTICE([No JAVA compiler. Disabling build of Java interface.])
    enable_java=no
  else
    AX_PROG_JAVAC_WORKS
  fi
fi

if test "$enable_java" != no ; then
  # look for jni header: required to compile C++ part of Java interface
  AX_JNI_INCLUDE_DIR
  if test -z "$JNI_INCLUDE_DIRS" ; then
    AC_MSG_NOTICE([No JNI header directory. Disabling build of Java interface.])
    enable_java=no
  else
    for JNI_INCLUDE_DIR in $JNI_INCLUDE_DIRS ; do
      CPPFLAGS="$CPPFLAGS -I$JNI_INCLUDE_DIR"
    done
  fi
fi

if test "$enable_java" != no ; then
  # the following macros can make configure stop with an error
  # we could work around that, but having javac and no jar, java, or javadoc would be odd anyway
  
  # look for jar: required to pack Java interface
  AX_PROG_JAR
  # look for more java to run tests and examples and do documentation
  AX_PROG_JAVA
  AX_PROG_JAVADOC
fi

AM_CONDITIONAL([BUILD_JAVA], [test "$enable_java" != no])

#########################
# Makefile conditionals #
#########################

# The following variable collects the names of libraries that should
# be included into libipopt.a (relative to subdir Interfaces, where it
# is made)

AC_SUBST(IPALLLIBS)
IPALLLIBS="../contrib/CGPenalty/libcgpenalty.la ../Algorithm/libipoptalg.la ../Algorithm/LinearSolvers/liblinsolvers.la ../Common/libcommon.la ../LinAlg/liblinalg.la ../LinAlg/TMatrices/libtmatrices.la"

AM_CONDITIONAL([BUILD_INEXACT], [test $use_inexact = yes])
if test $use_inexact = yes; then
  IPALLLIBS="../Algorithm/Inexact/libinexact.la $IPALLLIBS"
fi

########################################################################
##                      Linear solver loader                          ##
########################################################################

AC_ARG_ENABLE([linear-solver-loader],
  [AC_HELP_STRING([--disable-linear-solver-loader],[disable build of linear solver loader])],
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
            [Define to 1 if the linear solver loader should be compiled to allow dynamic loading of shared libraries with linear solvers])
  IPALLLIBS="../contrib/LinearSolverLoader/libLinearSolverLoader.la $IPALLLIBS"
fi

AM_CONDITIONAL([BUILD_LINEARSOLVERLOADER],[test $use_linearsolverloader = yes])

AC_LANG_PUSH(C)
AC_CHECK_HEADER([windows.h])
AC_CHECK_HEADER([dlfcn.h])
AC_CHECK_LIB(dl,[dlopen],[IPOPTLIB_LFLAGS="$IPOPTLIB_LFLAGS -ldl"],[])
AC_CHECK_FUNCS([snprintf _snprintf],[break])
AC_CHECK_FUNCS([vsnprintf _vsnprintf],[break])
AC_LANG_POP(C)

########################################################################
##                             sIpopt                                 ##
########################################################################

AC_ARG_ENABLE([sipopt],
  [AC_HELP_STRING([--disable-sipopt],[disable build of sIpopt])],
  [case "$enableval" in
     no | yes) ;;
     *)
       AC_MSG_ERROR([invalid argument for --enable-sipopt: $enableval]);;
   esac
   use_sipopt=$enableval],
  [use_sipopt=yes])
AM_CONDITIONAL(BUILD_SIPOPT, [test "$use_sipopt" = yes])

########################################################################
##             Create links for the test source files                 ##
########################################################################

AC_CONFIG_LINKS([test/hs071_main.cpp:examples/hs071_cpp/hs071_main.cpp
                 test/hs071_nlp.cpp:examples/hs071_cpp/hs071_nlp.cpp
                 test/hs071_nlp.hpp:examples/hs071_cpp/hs071_nlp.hpp
                 test/hs071_c.c:examples/hs071_c/hs071_c.c])

if test "$use_sipopt" = yes ; then
AC_CONFIG_LINKS([test/parametric_driver.cpp:contrib/sIPOPT/examples/parametric_cpp/parametric_driver.cpp
  test/parametricTNLP.cpp:contrib/sIPOPT/examples/parametric_cpp/parametricTNLP.cpp
  test/MySensTNLP.cpp:contrib/sIPOPT/examples/redhess_cpp/MySensTNLP.cpp
  test/redhess_cpp.cpp:contrib/sIPOPT/examples/redhess_cpp/redhess_cpp.cpp])
fi

########################################################################
##          Create links for VPATH config of certain files            ##
########################################################################

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
if test "$enable_java" != no ; then
  AC_COIN_VPATH_LINK(examples/hs071_java/HS071.java)
fi

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

if test "$use_sipopt" = yes ; then
  AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_cpp/parametric_driver.cpp)
  AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_cpp/parametricTNLP.cpp)
  AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_cpp/parametricTNLP.hpp)
  AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_dsdp_cpp/parametric_dsdp_driver.cpp)
  AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_dsdp_cpp/parametricTNLP.cpp)
  AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/parametric_dsdp_cpp/parametricTNLP.hpp)
  AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/redhess_cpp/MySensTNLP.cpp)
  AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/redhess_cpp/MySensTNLP.hpp)
  AC_COIN_VPATH_LINK(contrib/sIPOPT/examples/redhess_cpp/redhess_cpp.cpp)
fi

########################################################################
##          Create Makefiles and other stuff                          ##
########################################################################

AC_COIN_FINALIZE_FLAGS([IpoptLib IpoptAmplInterfaceLib])

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
		 ipopt.pc
		 doc/Doxyfile
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
])

if test -n "$F77" ; then
  AC_CONFIG_FILES([examples/hs071_f/hs071_f.f examples/hs071_f/Makefile
    tutorial/CodingExercise/Fortran/1-skeleton/TutorialFortran.f
    tutorial/CodingExercise/Fortran/2-mistake/TutorialFortran.f
    tutorial/CodingExercise/Fortran/3-solution/TutorialFortran.f
    tutorial/CodingExercise/Fortran/1-skeleton/Makefile
    tutorial/CodingExercise/Fortran/2-mistake/Makefile
    tutorial/CodingExercise/Fortran/3-solution/Makefile])
fi

if test "$enable_java" != no ; then
  AC_CONFIG_FILES([examples/hs071_java/Makefile examples/ScalableProblems_java/Makefile])
fi

if test "$use_sipopt" = yes ; then
  AC_COIN_FINALIZE_FLAGS([sIpoptLib])
  AC_CONFIG_FILES([
		 contrib/sIPOPT/Makefile
		 contrib/sIPOPT/src/Makefile
		 contrib/sIPOPT/AmplSolver/Makefile
		 contrib/sIPOPT/examples/parametric_cpp/Makefile
		 contrib/sIPOPT/examples/parametric_dsdp_cpp/Makefile
		 contrib/sIPOPT/examples/redhess_cpp/Makefile
  ])
fi


# under Windows, the Makevars file for the R Interface need to be named Makevars.win
case $build in 
  *-cygwin* | *-mingw* | *-msys* )
    AC_CONFIG_FILES([contrib/RInterface/src/Makevars.win:contrib/RInterface/src/Makevars.in])
    ;;
  *)
    AC_CONFIG_FILES([contrib/RInterface/src/Makevars])
    ;;
esac

if test $coin_has_asl = yes ; then
  AC_CONFIG_FILES([ipoptamplinterface.pc:src/Apps/AmplSolver/ipoptamplinterface.pc.in])
fi

AC_CONFIG_HEADER([src/Common/config.h src/Common/config_ipopt.h])

AC_COIN_FINALIZE
