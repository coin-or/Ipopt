# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AC_PREREQ(2.59)

AC_INIT([Ipopt],[3.1.0dev],[http://projects.coin-or.org/Ipopt/newticket])

AC_COPYRIGHT([
Copyright 2004, 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package IPOPT which is distributed
under the Common Public License.])

AC_REVISION(3.1.0dev)

AC_CONFIG_SRCDIR(Common/IpDebug.hpp)

# Get the system type
AC_CANONICAL_BUILD

AC_PREFIX_DEFAULT([`pwd`])

AC_PROG_LN_S 
AC_PROG_EGREP

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE

if test "$coin_debug_compile" = true; then
  AC_DEFINE([IP_DEBUG],[1],
            [If defined, debug sanity checks are performed in Ipopt during runtime])
fi

# Some extra stuff for blue gene on Laci's request
AC_SUBST(BLUEGENEEXTRA)
BLUEGENEEXTRA="echo"
AC_ARG_ENABLE([bg],,
  [case "${enableval}" in
    yes)  BLUEGENEEXTRA="chmod 644 *.o; cp ../Extern/ASL/amplsolver.a ../lib/libamplsolver.a"
          is_bg=yes
          ;;
  esac],[is_bg=no])

# Get the name of the C++ compiler and appropriate options
AC_COIN_PROG_CXX
AC_COIN_CXXLIBS

# We add one more definition in if we are using g++ compilers and are
# in debug mode
if test "$coin_debug_compile" = true; then
  if test "$GXX" = yes; then
    CXXFLAGS="$CXXFLAGS -DCHECK_SMARTPTR"
  fi
fi

# C compiler
AC_COIN_PROG_CC

# Fortran compiler
AC_COIN_PROG_F77

# Find out how to call Fortran from C and determine Fortran runtime libraries
AC_COIN_F77_WRAPPERS

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
case $ac_cv_sizeof_int_p in
  4) BITS_PER_POINTER=32
     BIT32FCOMMENT=''
     BIT64FCOMMENT='C' ;;
  8) BITS_PER_POINTER=64
     BIT32FCOMMENT='C'
     BIT64FCOMMENT='' ;;
  *) AC_MSG_ERROR([Unknown length of int *]);;
esac

# On Cygwin, building DLLs doesn't work
case $build in
  *-cygwin*)
    if test x"$enable_shared" = xyes; then
      AC_MSG_WARN([On Cygwin, shared objects are not supported. I'm disabling your choice.])
    fi
    enable_shared=no
  ;;
esac

# By default, we only want the shared objects to be compiled
AC_DISABLE_STATIC

# The codes is not yet "DLL clean" (ToDo?)
AC_LIBTOOL_WIN32_DLL
# Stuff for libtool
AC_COIN_INIT_AUTO_TOOLS

# Find out which options are needed to tell the linker about the
# shared objects path (for example Makefiles)
AC_COIN_RPATH_FLAGS(['$(IPOPTLIBDIR)'])

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

##########################################################################

# Additional libraries and packages

# IPADDLIBS will contain all libraries for additional packages
IPADDLIBS=
AC_SUBST(IPADDLIBS)

# We need to restore the value of LIBS before those test, otherwise libtool
# will hick-up
LIBS_without_IPADDLIBS="$LIBS"

##########
#  BLAS  #
##########

ipblaslibdir="\$(top_builddir)/Extern/blas"

AC_ARG_WITH([blas],
            AC_HELP_STRING([--with-blas],
                           [specify BLAS library (or BUILD for compilation)]),
            [use_blas=$withval], [use_blas=])

# On cygwin, unless otherwise specified, recompile blas because it
# otherwise links with the cygwin blas which doesn't run under DOS>
if test -z "$use_blas"; then
  case $build in
    *-cygwin*) use_blas="BUILD" ;;
  esac
fi

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

##########
# LAPACK #
##########

AC_LANG_PUSH([Fortran 77])

iplapacklibdir="\$(top_builddir)/Extern/lapack"
lapackdir=Extern/lapack/LAPACK/SRC
lapacksrcdir=$srcdir/$lapackdir
lapackobjdir=Extern/lapack

AC_ARG_WITH([lapack],
            AC_HELP_STRING([--with-lapack],
                           [specify LAPACK library (or BUILD for compilation, NONE for unavailable)]),
            [use_lapack=$withval], [use_lapack=])

# On cygwin, unless otherwise specified, recompile lapack because it
# otherwise links with the cygwin lapack which doesn't run under DOS>
if test -z "$use_lapack"; then
  case $build in
    *-cygwin*) use_lapack="BUILD" ;;
  esac
fi

# Check if user supplied option makes sense
if test ! -z "$use_lapack"; then
  if test "$use_lapack" = "BUILD"; then
    AC_CHECK_FILE([$lapacksrcdir/dlarf.f],[],
                  [AC_MSG_ERROR([option \"BUILD\" specified for lapack, but sources are missing])])
  else
    if test ! "$use_lapack" = "NONE"; then
      AC_MSG_CHECKING([whether user supplied LAPACKLIB=\"$use_lapack\" works])
      LIBS="$use_lapack $LIBS"
      IPADDLIBS="$use_lapack $IPADDLIBS"
      AC_TRY_LINK([],[      call dsyev],[AC_MSG_RESULT([yes])],
                  [AC_MSG_RESULT([no]); AC_MSG_ERROR([user supplied LAPACK library \"$use_lapack\" does not work])])
    fi
  fi
else
# First try to see if LAPACK is already available with BLAS library
  AC_MSG_CHECKING([whether LAPACK is already available with BLAS library])
  AC_TRY_LINK([],[      call dsyev()],
            [AC_MSG_RESULT([yes]); use_lapack=ok],
            [AC_MSG_RESULT([no])])
  if test -z "$use_lapack"; then
    # Try to autodetect the library for lapack based on build system
    case $build in
      *-sgi-*) AC_CHECK_LIB([complib.sgimath],[dsyev],
                            [use_lapack=-lcomplib.sgimath;
                             LIBS="-lcomplib.sgimath $LIBS";
                             IPADDLIBS="-lcomplib.sgimath $IPADDLIBS"]);;
      *-sun-*) SAVE_LIBS="$LIBS"
               AC_MSG_CHECKING([whether -xlic_lib=sunperf has LAPACK])
               LIBS="-xlic_lib=sunperf $LIBS"
               AC_TRY_LINK([],[      call dsyev],
                           [AC_MSG_RESULT([yes])
		  	    IPADDLIBS="-xlic_lib=sunperf $IPADDLIBS"
			    use_blas=-xlic_lib=sunperf],
                           [AC_MSG_RESULT([no])
			    LIBS="$SAVE_LIBS"])
             ;;
    esac
  fi
  if test -z "$use_lapack"; then
    AC_CHECK_LIB(lapack,dpotrs,[use_lapack=-llapack;
                                LIBS="-llapack $LIBS";
                                IPADDLIBS="-llapack $IPADDLIBS"],[])
  fi
  if test -z "$use_lapack"; then
    AC_CHECK_FILE([$lapacksrcdir/dlarf.f],[use_lapack=BUILD],
                  [AC_MSG_WARN([cannot find LAPACK, and sources are missing in $lapacksrcdir]); use_lapack=NONE])
  fi
fi

if test "$use_lapack" = "NONE"; then
  AC_MSG_WARN([Compiling code without LAPACK.  Certain options (e.g., quasi-Newton) will not work.])
else
  AC_DEFINE(HAVE_LAPACK, 1, [Define to 1 if LAPACK is avaliable])
fi

AC_LANG_POP([Fortran 77])

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
    if $EGREP MA27UD $hslsrcdir/ma27ad.f >/dev/null 2>&1; then
      AC_MSG_RESULT(yes)
    else
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([Please obtain the threadsafe (newer) version of MA27])
    fi
  fi
  MY_CHECK_HSLFILE(mc19ad,MC19)
  MY_CHECK_HSLFILE(ma57ad,MA57)
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
  MY_CHECK_HSL([ma57ad],[MA57])
  LIBS="$use_hsl $SAVE_LIBS"
fi

###########
# PARDISO #
###########

AC_LANG_PUSH([Fortran 77])
AC_ARG_WITH([pardiso],
            AC_HELP_STRING([--with-pardiso],
                           [specify Pardiso library (unsupported)]),
            [use_pardiso=$withval], [use_pardiso=no])

if test ! "$use_pardiso" = "no"; then
  LIBS="$use_pardiso $LIBS"
  IPADDLIBS="$use_pardiso $IPADDLIBS"

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

if test ! "$use_taucs" = "no"; then
  LIBS="$use_taucs $LIBS"
  IPADDLIBS="$use_taucs $IPADDLIBS"

  AC_DEFINE(HAVE_TAUCS,1,
            [Define to 1 if TAUCS is available])

  AC_ARG_VAR(TAUCS_INCFLAGS,[Flags for including taucs header files])
  CXXFLAGS="$CXXFLAGS $TAUCS_INCFLAGS"
fi

#########
# MUMPS #
#########

AC_ARG_WITH([mumps],
            AC_HELP_STRING([--with-mumps],
                           [specify MUMPS library and dependencies]),
            [use_mumps="$withval"], [use_mumps=no])

if test ! "$use_mumps" = "no"; then
  LIBS="$use_mumps $LIBS"
  IPADDLIBS="$use_mumps $IPADDLIBS"

  AC_DEFINE(HAVE_MUMPS,1,
            [Define to 1 if MUMPS is available])

  AC_ARG_VAR(MUMPS_INCFLAGS,[Flags for including mumps and (psydo) mpi header files.])
  if test -z "$MUMPS_INCFLAGS"; then
    AC_MSG_ERROR([You need to specify MUMPS_INCFLAGS if you want to use MUMPS.])
  fi
  CXXFLAGS="$CXXFLAGS $MUMPS_INCFLAGS"
fi

########
# WSMP #
########

AC_ARG_WITH([wsmp],
            AC_HELP_STRING([--with-wsmp],
                           [specify WSMP library]),
            [use_wsmp=$withval], [use_wsmp=no])

if test ! "$use_wsmp" = "no"; then
  LIBS="$use_wsmp $LIBS"
  IPADDLIBS="$use_wsmp $IPADDLIBS"

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

if test "$use_MA27" = "no" -a "$use_MA57" = "no" -a "$use_pardiso" = "no" -a "$use_taucs" = no -a "$use_wsmp" = no -a "$use_mumps" = no; then
  if test "$distconf" = false; then
    AC_MSG_ERROR("No linear solver (such as MA27) available.  You need to provide at least one.")
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

case $CXX in
  cl | */cl)   ampllib=amplsolv.lib ;;
  *)           ampllib=amplsolver.a ;;
esac

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
  AC_CHECK_FILE([$use_asldir/$ampllib],[],[AC_MSG_ERROR([ASL directory \"$use_asldir\" specified, but library missing])])
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
                if test "$is_bg" = "yes"; then
                  asladdcflags="$asladdcflags -DNo_Control87 -DBLUE_GENE"
                fi
                ;;
      *-cygwin*) # ToDo: If on Cygwin but not using -mno-cygwin, we should not define -DNo_Control87 -DNO_RUSAGE
                asladdcflags="-DNON_STDIO -DNo_Control87 -DNO_RUSAGE"
                ;;
      *-mingw32)
                asladdcflags="-DNON_STDIO -DNo_Control87 -DNO_RUSAGE"
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

#   Somehow, -qtune=auto makes ASL segfault on AIX compiler, so take them out
    case $CC in
      xlc* | */xlc* | mpxlc* | */mpxlc*)
        ASLMAKEFLAGS=`echo $ASLMAKEFLAGS | sed 's/-qtune=auto//g' | sed 's/-qarch=auto//g' | sed 's/-qcache=auto//g' | sed 's/-qhot//g' | sed 's/-O3/-O/g'`
      ;;
      pgcc | */pgcc)
        ASLMAKEFLAGS=`echo $ASLMAKEFLAGS | sed 's/-fast/-O1/g'`
      ;;
    esac

    AC_MSG_CHECKING([whether $am_make is the GNU make])
    rm -rf conftest.mak
    ($am_make -v -f conftest.mak > conftest.grp) 2>&5
    if $EGREP GNU conftest.grp >&5 ; then
      AC_MSG_RESULT([yes])
      ASLMAKEFLAGS="-j 1 $ASLMAKEFLAGS"
    else
      AC_MSG_RESULT([no])
    fi
    rm -f conftest.mak conftest.grp
  fi
  ASLMAKEFILECMDS="sed 's/ar \$\$(ARFLAGS)/$AR \$\$(ARFLAGS)/' makefile.u | sed 's/ar ruv/$AR ruv/' "
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

# For AmplApplication Makefile
ASLLIB=
AMPL_CPPFLAGS="-I\$(top_builddir)/Extern/ASL -I\$(top_srcdir)/Extern/ASL/solvers"

if test "$use_asldir" = BUILD; then
  ASLLIB="\$(top_builddir)/Extern/ASL//$ampllib"
  case $build in
    *-mingw32)
      # This is a hack for mingw where AC_CONFIG_LINKS doesn't work for directories
      rm -rf Extern/ASL/ASLsolvers
      $LN_S $srcdir/$asldir Extern/ASL/ASLsolvers
      ;;
    *)
      AC_CONFIG_LINKS(Extern/ASL/ASLsolvers:$asldir)
      ;;
  esac
elif test ! "$use_asldir" = no; then
  case "$CXX" in
    cl | */cl)
      use_asldir="`cygpath -w $use_asldir | sed -e sX\\\\\\\\X/Xg`"
      AMPL_CPPFLAGS="-I$use_asldir"
      ASLLIB="$use_asldir/$ampllib"
    ;;
  esac
fi
if test x"$ASLLIB" = x; then
  ASLLIB="$use_asldir/$ampllib"
fi

if test ! "$use_asldir" = no; then
  AC_CHECK_LIB(dl,[dlopen],[ASLLIB="$ASLLIB -ldl"],[])
fi

AC_SUBST(ASLLIB)
AC_SUBST(AMPL_CPPFLAGS)

AC_LANG_POP(C++)

######################################
# Equivalent int Fortran and C types #
######################################

if test "$cross_compiling" = no && test "$is_bg" != yes; then
  AC_LANG_PUSH(C)
  AC_DEFINE([FORTRAN_INTEGER_TYPE],[int],[Define to the C type corresponding to Fortran INTEGER])
  AC_CHECK_SIZEOF([long])
  AC_CHECK_SIZEOF([int])
  AC_CHECK_SIZEOF([double])
  AC_MSG_CHECKING([for C type corresponding to Fortran INTEGER])
  if test $ac_cv_sizeof_long = $ac_cv_sizeof_double; then
    AC_DEFINE([FORTRAN_INTEGER_TYPE],[int],[Define to the C type corresponding to Fortran INTEGER])
    AC_MSG_RESULT([int])
  else
    AC_DEFINE([FORTRAN_INTEGER_TYPE],[long])
    AC_MSG_RESULT([long])
  fi
  AC_LANG_POP(C)
else
  AC_MSG_WARN([We are cross compiling, assuming Fortran `INTEGER' type corresponds to C `int' type])
  AC_DEFINE([FORTRAN_INTEGER_TYPE],[int])
fi
IPADDLIBS="$IPADDLIBS"
LIBS="$LIBS_without_IPADDLIBS"

##########
# ASTYLE #
##########
ASTYLE=astyle
AC_SUBST(ASTYLE)
ASTYLEFLAGS="--mode=c --indent=spaces=2 --indent-cases --indent-namespaces --min-conditional-indent=1 --brackets=linux --brackets=break-closing-headers --max-instatement-indent=2"
AC_SUBST(ASTYLEFLAGS)

# The following one is for the example makefiles
case $build in
#  *-cygwin*) ipoptlib='`cygpath -w \$(IPOPTLIBDIR)`\\libipopt.a';;
  *)         ipoptlib="-L\$(IPOPTLIBDIR) -lipopt";;
esac
AC_SUBST(ipoptlib)

AC_COIN_ADD_WARN_FLAGS

#########################
# Makefile conditionals #
#########################

# The following variable collects the names of libraries that should
# be included into libipopt.a (relative to subdir Interfaces, where it
# is made)

AC_SUBST(IPALLLIBS)
IPALLLIBS="../Algorithm/libipoptalg.la ../Algorithm/LinearSolvers/liblinsolvers.la ../Common/libcommon.la ../LinAlg/liblinalg.la ../LinAlg/TMatrices/libtmatrices.la"

AM_CONDITIONAL([BUILDBLAS],[test "$use_blas" = "BUILD"])
if test "$use_blas" = BUILD; then
#  LIBS="$ipblaslibdir/libipblas.a $LIBS"
  IPALLLIBS="../Extern/blas/libipblas.la $IPALLLIBS"
fi

AM_CONDITIONAL([BUILDLAPACK],[test "$use_lapack" = "BUILD"])
if test "$use_lapack" = BUILD; then
#  LDFLAGS="-L$iplapacklibdir $LDFLAGS"
#  LIBS="$iplapacklibdir/libiplapack.a $LIBS"
  IPALLLIBS="../Extern/lapack/libiplapack.la $IPALLLIBS"
  for file in dlae2.f dlaev2.f dlamch.f dlanst.f dlansy.f dlapy2.f dlarf.f dlarfb.f dlarfg.f dlarft.f dlartg.f dlascl.f dlaset.f dlasr.f dlasrt.f dlassq.f dlatrd.f dorg2l.f dorg2r.f dorgql.f dorgqr.f dorgtr.f dpotf2.f dpotrf.f dpotrs.f dsteqr.f dsterf.f dsyev.f dsytd2.f dsytrd.f ieeeck.f ilaenv.f; do
    AC_CONFIG_LINKS($lapackobjdir/$file:$lapackdir/$file)
  done
fi

AM_CONDITIONAL([BUILDHSL],[test "$use_hsl" = "BUILD"])
if test "$use_hsl" = BUILD ; then
#  LDFLAGS="-L$hsllibdir $LDFLAGS"
#  LIBS="$hsllibdir/libiphsl.la $LIBS"
  IPALLLIBS="../Extern/HSL/libiphsl.la $IPALLLIBS"
fi
AM_CONDITIONAL([BUILDMA27],[test "$use_MA27" = BUILD])
AM_CONDITIONAL([BUILDMA57],[test "$use_MA57" = BUILD])
AM_CONDITIONAL([BUILDMC19],[test "$use_MC19" = BUILD])

AM_CONDITIONAL([BUILDBLASDIR],[test "$use_blas" = BUILD])
AM_CONDITIONAL([BUILDLAPACKDIR],[test "$use_lapack" = BUILD])
AM_CONDITIONAL([BUILDHSLDIR],[test "$use_hsl" = BUILD])

AM_CONDITIONAL([HAVE_MA27],[test ! "$use_MA27" = no])
AM_CONDITIONAL([HAVE_MA57],[test ! "$use_MA57" = no])
AM_CONDITIONAL([HAVE_MC19],[test ! "$use_MC19" = no])
AM_CONDITIONAL([HAVE_PARDISO],[test ! "$use_pardiso" = no])
AM_CONDITIONAL([HAVE_TAUCS],[test ! "$use_taucs" = no])
AM_CONDITIONAL([HAVE_MUMPS],[test ! "$use_mumps" = no])
AM_CONDITIONAL([HAVE_WSMP],[test ! "$use_wsmp" = no])

AM_CONDITIONAL([BUILDASLDIR],[test "$use_asldir" = BUILD])
AM_CONDITIONAL([BUILDAMPLINTERFACE],[test ! "$use_asldir" = "no"])

AC_CONFIG_FILES([Makefile
		 Common/Makefile
                 LinAlg/Makefile
		 LinAlg/TMatrices/Makefile 
		 Interfaces/Makefile
		 Interfaces/ipopt_addlibs_f.txt
		 Interfaces/ipopt_addlibs_cpp.txt
		 Algorithm/Makefile
		 Algorithm/LinearSolvers/Makefile
		 Apps/Makefile
		 Apps/AmplSolver/Makefile
		 Apps/CUTErInterface/Makefile
		 Examples/hs071_f/Makefile
		 Examples/hs071_f/hs071_f.f
		 Examples/Cpp_example/Makefile
		 Examples/hs071_cpp/Makefile
		 Examples/hs071_c/Makefile
		 Examples/ScalableProblems/Makefile
		 Extern/Makefile
		 Extern/blas/Makefile
		 Extern/lapack/Makefile
		 Extern/HSL/Makefile
		 Extern/ASL/Makefile
		 Test/run_unitTests])

AC_CONFIG_HEADER([inc/config_ipopt.h])

AC_OUTPUT

AC_MSG_NOTICE([Ipopt configuration successful])
