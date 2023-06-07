# Copyright (C) 2004, 2011 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([Ipopt],[3.14.13],[https://github.com/coin-or/Ipopt/issues/new],[],[https://github.com/coin-or/Ipopt])

AC_COPYRIGHT([
Copyright 2004, 2011 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package IPOPT which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/Common/IpDebug.hpp)

# Do some project-level initialization work (version numbers, ...)
AC_COIN_INITIALIZE

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

# IPOPT_VERBOSITY and IPOPT_DEBUGLEVEL
AC_COIN_DEBUGLEVEL

########################################################################
#      Floating-point Precision   and    Integer size                  #
########################################################################

AC_ARG_WITH([precision],
            [AS_HELP_STRING([--with-precision],[floating-point precision to use: single or double (default)])],
            [precision=$withval],
            [precision=double])
case "$precision" in
  single ) AC_DEFINE([IPOPT_SINGLE],[1],[Define to 1 if using single precision floating point]) ;;
  double ) ;;
  *) AC_MSG_ERROR([unsupported value $precision for option --with-precision])
esac
AM_CONDITIONAL([IPOPT_SINGLE],[test "$precision" = single])

# check which integer size to build
AC_ARG_WITH([intsize],
            [AS_HELP_STRING([--with-intsize],[integer type to use: specify 32 for int or 64 for int64_t])],
            [intsize=$withval],
            [intsize=32])
case "$intsize" in
  32 ) ;;
  64 ) AC_DEFINE(IPOPT_INT64, [1], [Define to 1 if Ipopt index type is int64_t]) ;;
  *) AC_MSG_ERROR([unsupported value $intsize for option --with-intsize]) ;;
esac
AM_CONDITIONAL([IPOPT_INT64],[test "$intsize" = 64])

# Specify 8-byte integers if intsize=64 (TODO This should become a macro in BuildTools that figures out the right flag to use)
if test $intsize = 64 ; then
  case "$F77" in
    *gfortran*) FFLAGS="$FFLAGS -fdefault-integer-8" ;;
    *ifort*)
      case $build in
        *-cygwin* | *-mingw* | *-msys* ) FFLAGS="$FFLAGS -integer-size:64" ;;
        *) FFLAGS="$FFLAGS -integer-size 64" ;;
      esac
      ;;
    "" ) ;; # no Fortran compiler, then don't do anything
    *) AC_MSG_ERROR([Do not know how to select 8-byte integers for Fortran compiler $FC]) ;;
  esac
fi

#############################################################################
#                                 Dependencies                              #
#############################################################################

AC_COIN_CHK_LIBM(IpoptLib)

AC_COIN_CHK_LAPACK(IpoptLib, int$intsize)
if test $coin_has_lapack != yes; then
  AC_MSG_ERROR([Required package LAPACK not found.])
#  AC_MSG_WARN([Compiling code without LAPACK.  Certain options (e.g., quasi-Newton) will not work.])
fi

# check for ASL only if building for double precision
if test $precision = double ; then
  AC_COIN_CHK_PKG(ASL,[IpoptAmplInterfaceLib SIpoptAmplInterfaceLib],[coinasl])
else
  AM_CONDITIONAL([COIN_HAS_ASL],[false])
fi

#########
# MUMPS #
#########

AC_COIN_CHK_PKG(Mumps,[IpoptLib],[coinmumps])

# Check whether MPI_Initialized is available
#   we assume that MPI_Finalized is present if MPI_Initialized is present
AC_CHECK_FUNCS([MPI_Initialized])

AC_ARG_ENABLE([mpiinit],
  [AS_HELP_STRING([--disable-mpiinit],[disable that (un)loading the Ipopt library initalizes (finalizes) MPI if the MPI version of MUMPS is linked])],
  [case "$enableval" in
     no | yes) ;;
     *)
       AC_MSG_ERROR([invalid argument for --enable-mpiinit: $enableval]);;
   esac
   mpiinit=$enableval],
  [mpiinit=yes])
if test $mpiinit = yes ; then
  AC_DEFINE(IPOPT_MPIINIT,1,[Define to 1 if MPI should be initialized (finalized) when Ipopt library is (un)loaded.])
fi

#######
# HSL #
#######

if test $intsize = 32 ; then
  AC_COIN_CHK_PKG(HSL,[IpoptLib HSLLib],[coinhsl])
else
  coin_has_hsl=no
  AM_CONDITIONAL([COIN_HAS_HSL],[false])
fi

have_ma28=no
if test "$coin_has_hsl" = yes ; then
  AC_COIN_FINALIZE_FLAGS([HSLLib])
  case $precision in
    single ) ma27func=ma27a  MA27FUNC=MA27A  ;;
    double ) ma27func=ma27ad MA27FUNC=MA27AD ;;
    * ) AC_MSG_ERROR([HSL cannot be used with --with-precision=$precision]) ;;
  esac
  AC_COIN_TRY_LINK(ma27ad,[$HSLLIB_LFLAGS],,
    [AC_COIN_DEFINENAMEMANGLING(IPOPT_HSL,$ma27ad_namemangling)],
    [AC_MSG_ERROR([Provided package HSL is not working or does not contain MA27. See config.log for details on failed checks.])])

  # extra check for MA28(s), since that decides whether we have to build IpMa28(s)Partition.F
  if test -n "$F77" ; then
    case $precision in
      single )
        AC_COIN_TRY_LINK([ma28a],[$HSLLIB_LFLAGS],,[
          if test "$ma28a_namemangling" != "$ac_cv_f77_mangling" ; then
            AC_MSG_WARN([Name mangling of MA28 different than Fortran. This will not link. Disabling MA28])
          else
            have_ma28=yes
          fi
        ]) ;;
      double )
        AC_COIN_TRY_LINK([ma28ad],[$HSLLIB_LFLAGS],,[
          if test "$ma28ad_namemangling" != "$ac_cv_f77_mangling" ; then
            AC_MSG_WARN([Name mangling of MA28 different than Fortran. This will not link. Disabling MA28])
          else
            have_ma28=yes
          fi
        ]) ;;
      * ) AC_MSG_ERROR([HSL cannot be used with --with-precision=$precision]) ;;
    esac
  fi
fi
AM_CONDITIONAL([HAVE_MA28],[test "$have_ma28" = yes])

###########
# PARDISO #
###########

if test $intsize = 32 ; then
AC_ARG_WITH([pardiso],
  [AS_HELP_STRING([--with-pardiso],[specify Pardiso library (>= 4.0) from pardiso-project.org])],
  [case "$withval" in
     no)  have_pardiso_project=no ;;  # so with_pardiso=no
     yes) AC_MSG_ERROR(["expected filename as argument for --with-pardiso"]) ;;
     *)   AC_DEFINE_UNQUOTED([PARDISO_LIB],["$withval"],[Name of Pardiso library from pardiso-project.org])
          AC_MSG_NOTICE([Pardiso from pardiso-project.org specified as $withval])
          have_pardiso_project=yes
          ;;
   esac],
   [have_pardiso_project=no])
fi

# check whether Pardiso is available via Lapack, which should then be MKL
AC_ARG_ENABLE([pardisomkl],
  [AS_HELP_STRING([--disable-pardisomkl],[disable check for MKL version of Pardiso in Lapack])],
  [case "$enableval" in
     no | yes) ;;
     *)
       AC_MSG_ERROR([invalid argument for --enable-pardisomkl: $enableval]);;
   esac
   checkpardisomkl=$enableval],
  [checkpardisomkl=yes])
have_pardiso_mkl=no
if test $checkpardisomkl = yes ; then
  AC_COIN_TRY_LINK([pardiso],[$lapack_lflags],[$lapack_pcfiles],
    [have_pardiso_mkl=yes
     AC_DEFINE(IPOPT_HAS_PARDISO_MKL,1,[Define to 1 if Pardiso from MKL is available])
    ])
fi
AM_CONDITIONAL([HAVE_PARDISO_MKL],[test "$have_pardiso_mkl" = yes])

#########
# SPRAL #
#########

# check for Spral only if building for double precision and 32-bit integers
# spral_ssids seems to support more than 2^31 nonzeros in a matrix, but not dimensions > 2^31
if test $precision$intsize = double32 ; then
  AC_LANG_PUSH(C++)
  AC_COIN_CHK_LIBHDR(Spral,[IpoptLib],[],[],[],
    [spral_ssids_solve(0, 0, (double*)0, 0, (void*)0, (void*)0, (const struct spral_ssids_options*)0, (struct spral_ssids_inform*)0);],
    [#include "spral_ssids.h"],
    [default_skip])
  AC_LANG_POP(C++)
else
  AM_CONDITIONAL([COIN_HAS_SPRAL],[false])
fi

########
# WSMP #
########

AC_ARG_WITH([wsmp],
            AS_HELP_STRING([--with-wsmp],[specify WSMP library]),
            [have_wsmp=yes; wsmp_lflags=$withval],
            [have_wsmp=no])

if test $have_wsmp$precision$intsize = yesdouble32 ; then
  AC_COIN_TRY_LINK([wssmp],[$wsmp_lflags],[],
    [AC_COIN_DEFINENAMEMANGLING([IPOPT_WSMP],[$wssmp_namemangling])
     IPOPTLIB_LFLAGS="$wsmp_lflags $IPOPTLIB_LFLAGS"
     AC_DEFINE(IPOPT_HAS_WSMP,1,[Define to 1 if WSMP is available])
    ],
    [AC_MSG_ERROR([Symbol wssmp not found with WSMP flags $wsmp_lflags. See config.log for details on failed checks.])])
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
  *) AC_MSG_ERROR([Unknown length of int*]);;
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

#########################
# va_copy and vsnprintf #
#########################

AC_CHECK_DECL([va_copy],[AC_DEFINE([IPOPT_HAS_VA_COPY],[1],
                                   [Define to 1 if va_copy is available])],,[#include <cstdarg>])
AC_CHECK_FUNCS([vsnprintf _vsnprintf],[break])

###########################
# Random number generator #
###########################

AC_CHECK_DECL([drand48],[AC_DEFINE([IPOPT_HAS_DRAND48],[1],
                                   [Define to 1 if function drand48 is available])],,[#include <cstdlib>])

AC_CHECK_DECL([rand],[AC_DEFINE([IPOPT_HAS_RAND],[1],
                                [Define to 1 if function rand is available])],,[#include <cstdlib>])

AC_COIN_CHECK_NAMESPACE_DECL([std::srand],[1],
  [AC_DEFINE([IPOPT_HAS_STD__RAND],[1],[Define to 1 if function std::rand is available])],[],
  [#include <cstdlib>])

##################
# feenableexcept #
##################

AC_CHECK_DECL([feenableexcept],
  [AC_DEFINE([IPOPT_HAS_FEENABLEEXCEPT],[1],[Define to 1 if function feenableexcept is available])],,
  [#include <cfenv>])

###################################################
# Check if user wants inexact algorithm available #
###################################################

AC_ARG_ENABLE([inexact-solver],
  [AS_HELP_STRING([--enable-inexact-solver],[enable inexact linear solver version EXPERIMENTAL! (default: no)])],
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

########
# Java #
########

AC_ARG_ENABLE([java],
  [AS_HELP_STRING([--disable-java],[disable building of Java interface])],
  [enable_java="$enableval"],
  [case "$JAVA_HOME" in
     *\ * ) enable_java=no ;;   # do not enable java-check by default, if there are spaces in JAVA_HOME - that causes trouble
     * )    enable_java="$enable_shared" ;;
   esac])

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
  AC_MSG_CHECKING([if javac supports -h])
  echo "public abstract class conftest { private native boolean test(); }" > conftest.java
  AS_ECHO(["$as_me:${as_lineno-$LINENO}: $JAVAC conftest.java -h conftest.header"]) >&AS_MESSAGE_LOG_FD
  "$JAVAC" conftest.java -h conftest.header >&AS_MESSAGE_LOG_FD
  if test -e conftest.header/conftest.h ; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
    cat conftest.java >&AS_MESSAGE_LOG_FD
    enable_java=no
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
  # look for jar: required to pack Java interface
  dnl copied from AX_PROG_JAR
  AS_IF([test "x$JAVAPREFIX" = x],
        [test "x$JAR" = x && AC_CHECK_PROGS([JAR], [jar])],
        [test "x$JAR" = x && AC_CHECK_PROGS([JAR], [jar], [], [$JAVAPREFIX/bin])])
  if test -z "$JAR" ; then
    AC_MSG_WARN([JAVA compiler found, but no acceptable jar program found in \$PATH. Disabling build of Java interface.])
    enable_java=no
  fi
fi

if test "$enable_java" != no ; then
  # look for java to run tests, but skip check whether java works (so just the test may fail)
  dnl copied from AX_PROG_JAVA
  AS_IF([test "x$JAVAPREFIX" = x],
        [test x$JAVA = x && AC_CHECK_PROGS([JAVA], [kaffe java])],
        [test x$JAVA = x && AC_CHECK_PROGS([JAVA], [kaffe java], [], [$JAVAPREFIX/bin])])
  if test -z "$JAVA" ; then
    AC_MSG_WARN([JAVA compiler found, but no acceptable JAVA executable found in \$PATH. Disabling tests of Java interface. Java example will not run.])
  fi

  # look for javadoc to build documentation
  dnl copied from AX_PROG_JAVADOC
  AS_IF([test "x$JAVAPREFIX" = x],
        [test "x$JAVADOC" = x && AC_CHECK_PROGS([JAVADOC], [javadoc])],
        [test "x$JAVADOC" = x && AC_CHECK_PROGS([JAVADOC], [javadoc], [], [$JAVAPREFIX/bin])])
  if test -z "$JAVADOC" ; then
    AC_MSG_WARN([JAVA compiler found, but no javadoc found in \$PATH. Building the Java documentation (make javadoc) will fail.])
  fi
fi

AM_CONDITIONAL([BUILD_JAVA], [test "$enable_java" != no])

########################################################################
##                      Linear solver loader                          ##
########################################################################

AC_ARG_ENABLE([linear-solver-loader],
  [AS_HELP_STRING([--disable-linear-solver-loader],[disable use of linear solver loader])],
  [case "$enableval" in
     no | yes) ;;
     *)
       AC_MSG_ERROR([invalid argument for --enable-linear-solver-loader: $enableval]);;
   esac
   use_linearsolverloader=$enableval],
  [use_linearsolverloader=yes])

AC_MSG_CHECKING([whether the linear solver loader should be used])
AC_MSG_RESULT([$use_linearsolverloader])

if test $use_linearsolverloader = yes; then
  AC_DEFINE([IPOPT_HAS_LINEARSOLVERLOADER],[1],
            [Define to 1 if the linear solver loader should be used to load libraries with linear solvers at runtime])

  AC_LANG_PUSH(C)
  AC_CHECK_LIB(dl,[dlopen],[IPOPTLIB_LFLAGS="$IPOPTLIB_LFLAGS -ldl"],[])
  AC_LANG_POP(C)
fi

AC_CHECK_HEADER([windows.h],AC_DEFINE(HAVE_WINDOWS_H,[1],[Define to 1 if windows.h is available.]))

########################################################################
##                             sIpopt                                 ##
########################################################################

AC_ARG_ENABLE([sipopt],
  [AS_HELP_STRING([--disable-sipopt],[disable build of sIpopt])],
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
  case $precision in
    single ) AC_COIN_VPATH_LINK(examples/hs071_java/HS071s.java) ;;
    double ) AC_COIN_VPATH_LINK(examples/hs071_java/HS071.java) ;;
  esac
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

# the SIpoptLib is here to get SIPOPTLIB_EXPORT
AC_COIN_FINALIZE_FLAGS([IpoptLib IpoptAmplInterfaceLib SIpoptLib SIpoptAmplInterfaceLib])

AC_CONFIG_FILES([Makefile
  src/ipopt.pc
  src/Makefile
  src/Apps/AmplSolver/Makefile
  src/Apps/AmplSolver/ipoptamplinterface.pc
  test/Makefile
  test/run_unitTests
  doc/Doxyfile
  examples/Cpp_example/Makefile
  examples/recursive_nlp/Makefile
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
  AC_CONFIG_FILES([examples/hs071_f/Makefile])
  case $precision in
    single )
      AC_CONFIG_FILES([examples/hs071_f/hs071_fs.f])
      ;;
    double )
      AC_CONFIG_FILES([examples/hs071_f/hs071_f.f
        tutorial/CodingExercise/Fortran/1-skeleton/TutorialFortran.f
        tutorial/CodingExercise/Fortran/2-mistake/TutorialFortran.f
        tutorial/CodingExercise/Fortran/3-solution/TutorialFortran.f
        tutorial/CodingExercise/Fortran/1-skeleton/Makefile
        tutorial/CodingExercise/Fortran/2-mistake/Makefile
        tutorial/CodingExercise/Fortran/3-solution/Makefile])
      ;;
  esac
fi

if test "$enable_java" != no ; then
  AC_CONFIG_FILES([examples/hs071_java/Makefile examples/ScalableProblems_java/Makefile])
fi

if test "$use_sipopt" = yes ; then
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

AC_CONFIG_HEADERS([src/Common/config.h src/Common/config_ipopt.h])

AC_COIN_FINALIZE
