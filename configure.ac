# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

AC_PREREQ(2.59)
#
# NOTE: The configure file in this distribution has been created with
# autoconf 2.59 with the modification that the line 650 in fortran.m4
# is replace by
#         -lang* | -lcrt[[012]].o | -lcrtbegin.o | -lc | -lgcc | -libmil | -LANG:=
*)
#

AC_INIT([ipopt],[2.2.1e],[coin-ipopt@www-124.ibm.com])

AC_COPYRIGHT([
Copyright 2003 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package IPOPT which is distributed
under the Common Public License.])

AC_REVISION($Revision: 676 $)

AC_CONFIG_SRCDIR(IPOPT/ipopt/ipopt.F)

# Get the system type
AC_CANONICAL_BUILD

if test -z "$AR"; then
  AR="ar"
fi
AC_ARG_VAR(AR,[Command for handling archives (default is `ar')])

# Check if user wants to produce debugging code
AC_ARG_ENABLE([debug],
  AC_HELP_STRING([--enable-debug], [compile with debug options]),
  [case "${enableval}" in
    yes) debugcompile=true ;;
    no)  debugcompile=false ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
  esac],[debugcompile=false])

######################
# C Compiler Options #
######################

SAVE_CFLAGS="$CFLAGS"
AC_PROG_CC([icc cl pgcc xlc cc gcc cl])
CFLAGS="$SAVE_CFLAGS"

# A few fixes for different systems
AC_AIX
AC_GNU_SOURCE
AC_ISC_POSIX
AC_MINIX

if test -z "$CFLAGS"; then
  add_cflags=
  opt_cflags=
  dbg_cflags=

  if test "$GCC" = "yes"; then
    case "$CC" in
      icc | */icc) ;;
      *)
        opt_cflags="-O3 -funroll-loops"
        add_cflags=
# The -mieee-fp flag doesn't seem to appear in the GCC documentation anymore,
# and it doesn't work
#        case $build in
#          *-linux*) add_cflags="$add_cflags -mieee-fp" ;;
#        esac
        dbg_cflags="-g"
        case $build in
          *-cygwin*) AC_MSG_CHECKING([whether -mno-cygwin works])
                    CFLAGS="-mno-cygwin"
                    AC_TRY_COMPILE(,[int i=0;],
                                    [AC_MSG_RESULT([yes]); add_cflags="-mno-cygwin $add_cflags"; use_mingw=yes],
                                    [AC_MSG_RESULT([no])])
                    CFLAGS= ;;
        esac
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
             opt_cflags="-O3"
             add_cflags="-mp -pc64"
             dbg_cflags="-g"
             # Check if -i_dynamic is necessary (for new glibc library)
             AC_MSG_CHECKING([whether flag `-i_dynamic' is needed])
             CFLAGS=
             AC_TRY_LINK(,[int i=0;],
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
             opt_cflags="-O -qmaxmem=-1"
             add_cflags="-bmaxdata:0x80000000"
             dbg_cflags="-g";;
      esac;;
    *-hp-*)
             opt_cflags="-O"
             add_cflags="-Ae"
             dbg_cflags="-g";;
    *-sun-*)
             opt_cflags="-O"
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
      AC_TRY_COMPILE([],[int i;],[opt_cflags="-O"],[])
    fi
    CFLAGS="$opt_cflags $add_cflags"
  fi
fi

# Try if CFLAGS works
AC_MSG_CHECKING([whether CFLAGS=\"$CFLAGS\" works])
AC_TRY_COMPILE([],[int i;],[],[CFLAGS=])
if test -z "$CFLAGS"; then
  AC_MSG_RESULT([no])
  AC_MSG_WARN([This value for CFLAGS does not work.  I will now just try '-O', but you might want to set CFLAGS manually.])
  CFLAGS='-O'
  AC_TRY_COMPILE([],[int i;],[],[CFLAGS=])
  if test -z "$CFLAGS"; then
    AC_MSG_RESULT([no])
    AC_MSG_WARN([This value for CFLAGS does not work.  I will continue with empty CFLAGS, but you might want to set CFLAGS manually.])
  fi
else
  AC_MSG_RESULT([yes])
fi

AC_PROG_RANLIB

AC_PREFIX_DEFAULT([`pwd`])

AC_PROG_LN_S 

# Stuff for automake
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

# This used to be further down, but it should be fine here (and not
# cause problems with missing Fortran runtime libraries
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

###################################
#   Fortran 77 compiler options   #
###################################

AC_LANG_PUSH(Fortran 77)
SAVE_FFLAGS="$FFLAGS"
AC_PROG_F77(ifort ifc pgf77 xlf fort77 f77 g77 df)
FFLAGS="$SAVE_FFLAGS"

if test -z "$FFLAGS"; then
  add_fflags=
  opt_fflags=
  dbg_fflags=

  if test "$G77" = "yes"; then
    opt_fflags="-O3 -funroll-loops"
    add_fflags=
# The -mieee-fp flag doesn't seem to appear in the GCC documentation anymore,
# and it doesn't work
#    case $build in
#      *-linux*) add_fflags="$add_fflags -mieee-fp" ;;
#    esac
    dbg_fflags="-g"
    case $build in
      *-cygwin*) AC_MSG_CHECKING([whether -mno-cygwin works])
                FFLAGS="-mno-cygwin"
                AC_TRY_COMPILE(,[      write(*,*) 'Hello world'],
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
             opt_fflags="-O3"
             add_fflags="-mp -pc64 -cm -w90 -w95"
             dbg_fflags="-g -CA -CB -CS -CU"
             # Check if -i_dynamic is necessary (for new glibc library)
             AC_MSG_CHECKING([whether flag `-i_dynamic' is needed])
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
             opt_fflags="-O -qmaxmem=-1"
             add_fflags="-bmaxdata:0x80000000"
             dbg_fflags="-g -C";;
      esac;;
    *-hp-*)
             opt_fflags="+O2"
             add_fflags="+U77"
             dbg_fflags="-C -g";;
    *-sun-*)
             opt_fflags="-O"
             dbg_fflags="-g";;
    *-sgi-*)
             opt_fflags="-O -OPT:Olimit=0"
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
      AC_TRY_COMPILE(,[      integer i],[AC_MSG_RESULT([yes]); opt_fflags="-O"],[AC_MSG_RESULT([no])])
    fi
    FFLAGS="$opt_fflags $add_fflags"
  fi
fi

# Try if FFLAGS works
AC_MSG_CHECKING([whether FFLAGS=\"$FFLAGS\" works])
AC_TRY_COMPILE([],[      integer i],[],[FFLAGS=])
if test -z "$FFLAGS"; then
  AC_MSG_RESULT([no])
  AC_MSG_WARN([This value for FFLAGS does not work.  I will now jst try '-O', but you might want to set FFLAGS manually.])
  FFLAGS='-O'
  AC_TRY_COMPILE([],[      integer i],[],[FFLAGS=])
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
esac

#Find out if compiler understands %VAL
AC_MSG_CHECKING([whether Fortran compiler understands \"%VAL\"])
AC_TRY_COMPILE([],[
      integer i
      call sub(%VAL(i))],
      [AC_MSG_RESULT([yes]) AC_DEFINE(USE_MALLOC,1,
                                      [Define to 1 if Fortran compiler understand %VAL and malloc should be used]) use_malloc=yes],
      [AC_MSG_RESULT([no]) use_malloc=no])
if test "$use_malloc" = yes; then
  AC_LANG_PUSH(C)
  AC_CHECK_SIZEOF(int *)
  AC_LANG_POP(C)
  case $ac_cv_sizeof_int_p in
    4) AC_DEFINE(SIZEOF_INT_P_IS_4,1,
                 [Define to 1 if pointers consists of 4 bytes]) ;;
    8) AC_DEFINE(SIZEOF_INT_P_IS_8,1,
                 [Define to 1 if pointers consists of 8 bytes]) ;;
    *) AC_MSG_ERROR([Unknown length of int *]);;
  esac
fi

# Fortran Proprocessor
AC_ARG_VAR(FPP,[Preprocessor for Fortran code (is used as `$FPP input.F >output.f'), e.g. `cpp -E -P'])
case $CPP in
"icc -E") FPP="icc -EP" ;;
*)   FPP="$CPP -P" ;;
esac

##########################################################################

# Additional libraries and packages

IPADDLIBS=
AC_SUBST(IPADDLIBS)

##########
#  BLAS  #
##########

ipblaslibdir="\$(top_builddir)/OTHERS/blas"

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
    *-sun-*) SAVE_LIBS="$LIBS"
             AC_MSG_CHECKING([whether -xlic_lib=sunperf has BLAS])
             LIBS="-xlic_lib=sunperf $LIBS"
             AC_TRY_LINK([],[      call dgemm],
                         [AC_MSG_RESULT([yes])
			  IPADDLIBS="-xlic_lib=sunperf $IPADDLIBS"
			  use_blas=-xlic_lib=sunperf],
                         [AC_MSG_RESULT([no])
			  LIBS="$SAVE_LIBS"])
             ;;
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
# D1MACH #
##########

# Is this already in blas library?
AC_MSG_CHECKING([whether d1mach is already available])
AC_TRY_LINK([],[
      double precision macheps, d1mach
      macheps = d1mach(4)],
          [AC_MSG_RESULT([yes]); use_d1mach=ok],
          [AC_MSG_RESULT([no]); use_d1mach=])
if test -z "$use_d1mach"; then
  AC_CHECK_FILE([$srcdir/OTHERS/blas/d1mach.f],[use_d1mach=BUILD],
                [AC_MSG_ERROR([cannot find source for d1mach subroutine])])
fi

##########
# LAPACK #
##########

iplapacklibdir="\$(top_builddir)/OTHERS/lapack"
lapackdir=OTHERS/lapack/LAPACK/SRC
lapacksrcdir=$srcdir/$lapackdir
lapackobjdir=OTHERS/lapack

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
      *-sgi-*) AC_CHECK_LIB([complib.sgimath],[dgetrf],
                            [use_lapack=-lcomplib.sgimath;
                             LIBS="-lcomplib.sgimath $LIBS";
                             IPADDLIBS="-lcomplib.sgimath $IPADDLIBS"]);;
      *-sun-*) SAVE_LIBS="$LIBS"
               AC_MSG_CHECKING([whether -xlic_lib=sunperf has LAPACK])
               LIBS="-xlic_lib=sunperf $LIBS"
               AC_TRY_LINK([],[      call dgetrf],
                           [AC_MSG_RESULT([yes])
		  	    IPADDLIBS="-xlic_lib=sunperf $IPADDLIBS"
			    use_blas=-xlic_lib=sunperf],
                           [AC_MSG_RESULT([no])
			    LIBS="$SAVE_LIBS"])
             ;;
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

#######################
# HARWELL SUBROUTINES #
#######################

hslsrcdir="$srcdir/OTHERS/HSL"
hsllibdir="\$(top_builddir)/OTHERS/HSL"

AC_DEFUN([MY_CHECK_HSL],
         [AC_MSG_CHECKING([for $1 in $use_hsl])
          use_$2=
          AC_TRY_LINK([],[      call $1()],
                         [AC_DEFINE(HAVE_$2,1,[Define to 1 if $2 is available]) AC_MSG_RESULT([yes]); use_$2=yes],
                         [AC_MSG_RESULT([no])])])
AC_DEFUN([MY_CHECK_HSLFILE],
         [AC_CHECK_FILE([$hslsrcdir/$1.f],
                        [AC_DEFINE(HAVE_$2,1,[Define to 1 if $2 is available]) use_$2=BUILD],
                        [use_$2=])])

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
    if grep MA27UD $hslsrcdir/ma27ad.f >&5 ; then
      AC_MSG_RESULT(yes)
    else
      AC_MSG_RESULT(no)
      AC_DEFINE(OLD_MA27,[],
                [Define to 1 if MA27 is given in the old, non-threadsafe version])
    fi
  fi
  MY_CHECK_HSLFILE(ma28ad,MA28)
  MY_CHECK_HSLFILE(ma47ad,MA47)
  MY_CHECK_HSLFILE(ma48ad,MA48)
  MY_CHECK_HSLFILE(ma57ad,MA57)
  MY_CHECK_HSLFILE(mc19ad,MC19)
  MY_CHECK_HSLFILE(mc29ad,MC29)
  MY_CHECK_HSLFILE(mc30ad,MC30)
  MY_CHECK_HSLFILE(mc35ad,MC35)
  MY_CHECK_HSLFILE(mc39ad,MC39)
else
  LIBS="$use_hsl $LIBS"
  IPADDLIBS="$use_hsl $IPADDLIBS"
  MY_CHECK_HSL([ma27ad],[MA27])
  if test "$use_MA27" = "yes"; then
    AC_MSG_CHECKING([whether MA27 is threadsafe])
    AC_TRY_LINK([],[      call MA27UD()],
                   [AC_MSG_RESULT(yes)],
                   [AC_MSG_RESULT(no)
                    AC_DEFINE(OLD_MA27,[],
                        [Define to 1 if MA27 is given in the old, non-threadsafe version])])
  fi
  MY_CHECK_HSL([ma28ad],[MA28])
  MY_CHECK_HSL([ma47ad],[MA47])
  MY_CHECK_HSL([ma48ad],[MA48])
  MY_CHECK_HSL([ma57ad],[MA57])
  MY_CHECK_HSL([mc19ad],[MC19])
  MY_CHECK_HSL([mc29ad],[MC29])
  MY_CHECK_HSL([mc30ad],[MC30])
  MY_CHECK_HSL([mc35ad],[MC35])
  MY_CHECK_HSL([mc39ad],[MC39])
fi

# Try to make sure the correct FD05 routine is there (if needed)
if test "$use_MA48" = BUILD || test "$use_MA57" = BUILD; then
  AC_CHECK_FILE([$hslsrcdir/fd05ad.f],[use_FD05="fd05ad.f"],[use_FD05=])
  if test -z "$use_FD05"; then
    case $build in
      rs6000-* | powerpc*)
                use_FD05="fd05adj.f";;
      vax-*)    use_FD05="fd05adv.f";;
      hp*)      use_FD05="fd05adh.f";;
      *-osf*)   use_FD05="fd05ada.f";;
      *ymp*)    use_FD05="fd05adc.f";;
      *)        use_FD05="fd05ads.f";;
    esac
    AC_MSG_WARN([The file `fd05ad.f' is not provided.  My guess for this system is that you should use `$use_FD05', but you better verify this!!!])
    AC_CHECK_FILE([$hslsrcdir/$use_FD05],
                  [],[AC_MSG_ERROR([File $use_FD05 not found])])
    AC_CONFIG_LINKS([OTHERS/HSL/fd05ad.f:OTHERS/HSL/$use_FD05])
  fi
fi

AC_LANG_POP(Fortran 77)

# For the remaining tests we need to include FLIBS into LIBS, so that
# the C compiler can link programs
LIBS="$LIBS $FLIBS"

#######################
# AMPL SOLVER LIBRARY #
#######################

AC_LANG_PUSH(C)

asldir=OTHERS/ASL/solvers
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
    ASLMAKEFLAGS="CFLAGS=\"\$(CFLAGS) $asladdcflags\""
    if test $USE_MAINTAINER_MODE = 'yes'; then
      AC_MSG_WARN([Assuming that maintainer mode is only used for GNU make])
      ASLMAKEFLAGS="-j 1 $ASLMAKEFLAGS"
    else
      AC_MSG_CHECKING([whether $am_make understands the `-j 1' flag])
      AC_LANG_CONFTEST([[bla:]])
      if ($am_make -f conftest.c -j 1 bla >&5) 2>&5 ; then
        AC_MSG_RESULT([yes])
        ASLMAKEFLAGS="-j 1 $ASLMAKEFLAGS"
      else
        AC_MSG_RESULT([no])
        AC_MSG_CHECKING([whether $am_make understands the `-J 1' flag])
        if ($am_make -f conftest.c -J 1 bla >&5) 2>&5 ; then
          AC_MSG_RESULT([yes])
          ASLMAKEFLAGS="-J 1 $ASLMAKEFLAGS"
        else
          AC_MSG_RESULT([no])
        fi
      fi
    fi
  fi
  ASLMAKEFILECMDS="sed 's/ar ruv/$AR ruv/' makefile.u"
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

AC_ARG_ENABLE([mpcc],
  AC_HELP_STRING([--enable-mpcc], [compile AMPL executable with support for complementarity constraints]),
  [case "$enableval" in
    yes) mpcccompile=true;;
    no)  mpcccompile=false;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpcc) ;;
  esac],[mpcccompile=false])
if test "$mpcccompile" = true; then
  AC_DEFINE(INCLUDE_CC,[],[Define to 1 if you want to compile the AMPL executable with support for MPECs])
fi

AC_LANG_POP(C)

########
# TRON #
########

tronsrcdir=OTHERS/TRON/tron
tronbuilddir=OTHERS/TRON
iptronlibdir="\$(top_builddir)/OTHERS/TRON"

PATCH="patch"
AC_SUBST(PATCH)

AC_ARG_ENABLE([tron],
  AC_HELP_STRING([--enable-tron], [include TRON for restoration phase]),
  [case "$enableval" in
    yes) use_tron=BUILD;;
    no)  use_tron=no;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-tron) ;;
  esac],[use_tron=check])
if test ! "$use_tron" = no; then
  AC_CHECK_FILE([$srcdir/$tronsrcdir/src/tron/dtron.f],
                [use_tron=BUILD],[have_tron=false])
  if test "$use_tron" = BUILD && test "$have_tron" = false; then
    AC_MSG_ERROR([TRON source files are missing])
  fi
fi
if test "$use_tron" = BUILD; then
  AC_DEFINE(HAVE_TRON,[],[Define to 1 if TRON is available])
fi


# Stuff that we need for C programs

###################
# Function finite #
###################

# First check for a few headers

AC_CHECK_HEADERS([math.h])
AC_CHECK_HEADERS([float.h])
case $build in
  *-cygwin*) ;;
  *) AC_CHECK_HEADERS([ieeefp.h]) ;;
esac

# now try to guess the correct name of the finite function

MY_C_FINITE=
AC_CHECK_DECL([finite],[MY_C_FINITE=finite],,[
#ifdef HAVE_MATH_H
# include <math.h>
#endif
#ifdef HAVE_FLOAT_H
# include <float.h>
#endif
#ifdef HAVE_IEEEFP_H
# include <ieeefp.h>
#endif])
if test -z "$MY_C_FINITE"; then
  AC_CHECK_DECL([_finite],[MY_C_FINITE=_finite],,[
#ifdef HAVE_MATH_H
# include <math.h>
#endif
#ifdef HAVE_FLOAT_H
# include <float.h>
#endif
#ifdef HAVE_IEEEFP_H
# include <ieeefp.h>
#endif])
  if test -z "$MY_C_FINITE"; then
    AC_CHECK_DECL([isfinite],[MY_C_FINITE=isfinite],,[
#ifdef HAVE_MATH_H
# include <math.h>
#endif
#ifdef HAVE_FLOAT_H
# include <float.h>
#endif
#ifdef HAVE_IEEEFP_H
# include <ieeefp.h>
#endif])
  fi
fi
if test -z "$MY_C_FINITE"; then
  AC_MSG_WARN(Cannot find C-function for checking NaN and Inf.  Check will be disabled)
else
  AC_DEFINE_UNQUOTED([MY_C_FINITE],[$MY_C_FINITE],
                     [Define to be the name of C-function for NaNInf check])
fi

################################
# Measuring CPU time (Fortran) #
################################

AC_LANG_PUSH(Fortran 77)
have_timer=
AC_MSG_CHECKING([whether ETIME is available])
AC_TRY_LINK([],[
      real t, tarray(2)
      real ETIME
      t = ETIME(tarray)],
               [have_timer=yes; AC_MSG_RESULT([yes]); AC_DEFINE([HAVE_ETIME],[],
               [Define to 1 if ETIME can be use to measure CPU time])],
               [AC_MSG_RESULT([no])])
if test -z "$have_timer"; then
  AC_MSG_CHECKING([whether CPU_TIME is available])
  AC_TRY_LINK([],[
      real t
      call CPU_TIME(time)],
                 [have_timer=yes; AC_MSG_RESULT([yes]); AC_DEFINE([HAVE_CPU_TIME],[],
                 [Define to 1 if CPU_TIME can be use to measure CPU time])],
                 [AC_MSG_RESULT([no])])
fi
if test -z "$have_timer"; then
  AC_MSG_WARN([Cannot figure out how to measure CPU time - this will be disabled])
fi

AC_LANG_POP(Fortran 77)

# The following variable collects the names of libraries that should be
# included into libipopt.a
AC_SUBST(IPALLLIBS)
IPALLLIBS=

# Makefile conditionals

AM_CONDITIONAL([USE_MALLOC],[test "$use_malloc" = "yes"])
AM_CONDITIONAL([BUILDBLAS],[test "$use_blas" = "BUILD"])
AM_CONDITIONAL([BUILDD1MACH],[test "$use_d1mach" = "BUILD"])
if test "$use_blas" = BUILD || test "$use_d1mach" = "BUILD"; then
#  LDFLAGS="-L$ipblaslibdir $LDFLAGS"
  LIBS="$ipblaslibdir/libipblas.a $LIBS"
  IPALLLIBS="OTHERS/blas/libipblas.a $IPALLLIBS"
fi

AM_CONDITIONAL([BUILDLAPACK],[test "$use_lapack" = "BUILD"])
if test "$use_lapack" = BUILD; then
#  LDFLAGS="-L$iplapacklibdir $LDFLAGS"
  LIBS="$iplapacklibdir/libiplapack.a $LIBS"
  IPALLLIBS="OTHERS/lapack/libiplapack.a $IPALLLIBS"
  for file in dgesv.f dgetf2.f dgetrf.f dgetrs.f dlabad.f dlacon.f dlacpy.f dlae2.f dlaebz.f dlaev2.f dlagtf.f dlagts.f dlamch.f dlansp.f dlanst.f dlapy2.f dlaqsp.f dlarf.f dlarfg.f dlarnv.f dlartg.f dlaruv.f dlascl.f dlaset.f dlasr.f dlasrt.f dlassq.f dlaswp.f dlatps.f dopgtr.f dopmtr.f dorg2l.f dorg2r.f dppcon.f dppequ.f dpprfs.f dppsv.f dppsvx.f dpptrf.f dpptrs.f drscl.f dspev.f dspevx.f dsptrd.f dstebz.f dstein.f dsteqr.f dsterf.f ieeeck.f ilaenv.f; do
    AC_CONFIG_LINKS($lapackobjdir/$file:$lapackdir/$file)
  done
fi

AM_CONDITIONAL([BUILDHSL],[test "$use_hsl" = "BUILD"])
if test "$use_hsl" = BUILD ; then
#  LDFLAGS="-L$hsllibdir $LDFLAGS"
  LIBS="$hsllibdir/libiphsl.a $LIBS"
  IPALLLIBS="OTHERS/HSL/libiphsl.a $IPALLLIBS"
fi

AM_CONDITIONAL([BUILDMA27],[test "$use_MA27" = BUILD])
AM_CONDITIONAL([BUILDMA28],[test "$use_MA28" = BUILD])
AM_CONDITIONAL([BUILDMA47],[test "$use_MA47" = BUILD])
AM_CONDITIONAL([BUILDMA48],[test "$use_MA48" = BUILD])
AM_CONDITIONAL([BUILDMA57],[test "$use_MA57" = BUILD])
AM_CONDITIONAL([BUILDMC19],[test "$use_MC19" = BUILD])
AM_CONDITIONAL([BUILDMC29],[test "$use_MC29" = BUILD])
AM_CONDITIONAL([BUILDMC30],[test "$use_MC30" = BUILD])
AM_CONDITIONAL([BUILDMC35],[test "$use_MC35" = BUILD])
AM_CONDITIONAL([BUILDMC39],[test "$use_MC39" = BUILD])

AM_CONDITIONAL([BUILDFD05],[test "$use_MA48" = BUILD ||
                            test "$use_MA57" = BUILD])
AM_CONDITIONAL([BUILDMC13],[test "$use_MA28" = BUILD ||
                            test "$use_MA48" = BUILD])
AM_CONDITIONAL([BUILDMC20],[test "$use_MA28" = BUILD ||
                            test "$use_MC35" = BUILD])
AM_CONDITIONAL([BUILDMC21],[test "$use_MA28" = BUILD ||
                            test "$use_MA48" = BUILD])

AM_CONDITIONAL([BUILDBLASDIR],[test "$use_blas" = BUILD ||
                               test "$use_d1mach" = BUILD])
AM_CONDITIONAL([BUILDLAPACKDIR],[test "$use_lapack" = BUILD])
AM_CONDITIONAL([BUILDHSLDIR],[test "$use_hsl" = BUILD])

AM_CONDITIONAL([BUILDASLDIR],[test "$use_asldir" = BUILD])
AM_CONDITIONAL([BUILDAMPLINTERFACE],[test ! "$use_asldir" = "no"])
AC_SUBST(ASLLIB)
AC_SUBST(AMPL_CPPFLAGS)
if test "$use_asldir" = BUILD; then
  ASLLIB="\$(top_builddir)/OTHERS/ASL/amplsolver.a"
  AMPL_CPPFLAGS="-I\$(top_builddir)/OTHERS/ASL -I\$(top_srcdir)/OTHERS/ASL/solvers"
  AC_CONFIG_LINKS(OTHERS/ASL/ASLsolvers:OTHERS/ASL/solvers)
elif test ! "$use_asldir" = no; then
  case $build in
    *-cygwin*) AMPL_CPPFLAGS="-I\"$use_asldir\""
               ASLLIB="\"$use_asldir\amplsolver.a\"" ;;
    *)         AMPL_CPPFLAGS="-I$use_asldir"
               ASLLIB="$use_asldir/amplsolver.a";;
  esac
fi

AM_CONDITIONAL([BUILDTRONDIR],[test "$use_tron" = BUILD])
if test "$use_tron" = BUILD; then
#  for file in dtron.f dspcg.f dcauchy.f dtrpcg.f dgpstep.f dprsrch.f dtrqsol.f dmid.f dgpnrm2.f; do
  for file in dtrqsol.f dmid.f dgpnrm2.f dbreakpt.f; do
    AC_CONFIG_LINKS($tronbuilddir/$file:$tronsrcdir/src/tron/$file)
  done
  for file2 in dicfs.f dssyax.f dstrsol.f dicf.f ihsort.f insort.f dsel2.f; do
    AC_CONFIG_LINKS($tronbuilddir/$file2:$tronsrcdir/src/icf/$file2)
  done
#  LDFLAGS="-L$iptronlibdir $LDFLAGS"
  LIBS="$iptronlibdir/libiptron.a \$(top_builddir)/IPOPT/ipopt/libipcore.a $LIBS"
  IPALLLIBS="OTHERS/TRON/libiptron.a $IPALLLIBS"
fi

AM_CONDITIONAL([USE_WSMP],[test "$use_wsmp" = "yes"])

AM_CONDITIONAL([DUMMYCUTER],[test 1 = 0])

LIBS="\$(top_builddir)/IPOPT/ipopt/libipcore.a $LIBS"
#LDFLAGS="-L\$(top_builddir)/IPOPT/ipopt $LDFLAGS"
IPALLLIBS="IPOPT/ipopt/libipcore.a $IPALLLIBS"

# # Now make sure that the Fortran runtime libraries are in LIBS
# LIBS="$LIBS $FLIBS"

# The following is necessary for HP
case $build in
   *-hp-*)
       IPADDLIBS="$IPADDLIBS -lm";;
esac


IPINCLUDEDIR="\$(top_builddir)/IPOPT/include"
AC_SUBST(IPINCLUDEDIR)

# The following one is for the example makefiles
case $build in
  *-cygwin*) ipoptlib='`cygpath -w \$(IPOPTLIBDIR)`\\libipopt.a';;
  *)         ipoptlib="-L\$(IPOPTLIBDIR) -lipopt";;
esac
AC_SUBST(ipoptlib)

# The following is a fix for the AIX Fortran compiler (preprocessing)
AC_SUBST(SETPPF77COMPILE)
case $F77 in
  xlf* | */xlf* | mpxlf* | */mpxlf* ) 
                 SETPPF77COMPILE="PPF77COMPILE = \$(F77) -WF,\"-I$IPINCLUDEDIR\" -WF,\"\$(DEFS)\" \$(DEFAULT_INCLUDES) \$(INCLUDES) \$(AM_CPPFLAGS) \$(CPPFLAGS) \$(AM_FFLAGS) \$(FFLAGS)" ;;
esac

AC_CONFIG_FILES([Makefile OTHERS/Makefile OTHERS/blas/Makefile
                 OTHERS/lapack/Makefile OTHERS/HSL/Makefile
                 OTHERS/ASL/Makefile OTHERS/TRON/Makefile
                 IPOPT/Makefile
                 IPOPT/ipopt/Makefile IPOPT/AMPL_interface/Makefile
                 IPOPT/CUTEr_interface/Makefile
                 IPOPT/include/Makefile EXAMPLES/Fortran_example/Makefile
                 EXAMPLES/C_example/Makefile])

AC_CONFIG_HEADER([IPOPT/include/config.h IPOPT/include/Ipopt.h IPOPT/include/Ipopt_Interface.h])

AC_OUTPUT

AC_MSG_NOTICE([Configuration successful])
