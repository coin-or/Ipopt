# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

## SGI: ./configure CC=cc F77=f77 CXX=CC CXXFLAGS="-g -LANG:std -DOLD_C_HEADERS" AR="CC -ar" ARFLAGS="-o " FFLAGS=-g CFLAGS=-g CXXLIBS='-lCsup -lC -lCio'
## AIX: ./configure CC=xlc F77=xlf CXX=xlC --with-lapack=$HOME/lib/liblapack.a CFLAGS=-g FFLAGS=-g CXXFLAGS="-g -qrtti=dyna" CXXLIBS="-lC -lc"
## HP1: configure --with-asl=$HOME/COIN/OBJhp1/OTHERS/ASL/amplsolver.a --with-blas=/opt/fortran/lib/libblas.a --with-hsl=$HOME/COIN/OBJhp1/OTHERS/HSL/libiphsl.a --with-lapack=$HOME/COIN/OBJhp1/OTHERS/lapack/libiplapack.a F77=/opt/fortran/bin/fort77 CC=/opt/ansic/bin/cc CXX=/opt/aCC/bin/aCC FFLAGS='-g' CFLAGS='-g' CXXFLAGS='-AA -DMY_C_FINITE=isfinite -g' CXXLIBS='-L/opt/aCC/lib -l++ -lstd_v2 -lCsup_v2 -lm -lcl -lc'

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
if test -z "$ARFLAGS"; then
  ARFLAGS="cru"
fi
AC_ARG_VAR(ARFLAGS,[Flags for archive command (default is `crv')])

# Check if user wants to produce debugging code
AC_ARG_ENABLE([debug],
  AC_HELP_STRING([--enable-debug], [compile with debug options and runtime tests]),
  [case "${enableval}" in
    yes) debugcompile=true ;;
    no)  debugcompile=false ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
  esac],[debugcompile=false])

####################
# Compiler Options #
####################

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
  MY_CHECK_HSLFILE(mc19ad,MC19)
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
  MY_CHECK_HSL([mc19ad],[MC19])
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

###########
# PARDISO #
###########

AC_ARG_WITH([pardiso],
            AC_HELP_STRING([--with-pardiso],
                           [specify Pardiso library]),
            [use_pardiso=$withval], [use_pardiso=no])

if test ! "$use_pardiso" = "no"; then
  LIBS="$use_pardiso $LIBS"
  IPADDLIBS="$use_pardiso $IPADDLIBS"
  AC_MSG_CHECKING([whether user supplied Pardiso library \"$use_pardiso\" works])
  AC_TRY_LINK([],[      call PARDISO()],
              [AC_MSG_RESULT(yes)],
              [AC_MSG_ERROR([Pardiso subroutine not found in $use_pardiso])])

  AC_DEFINE(HAVE_PARDISO,1,
            [Define to 1 if Pardiso is available])
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

AC_LANG_POP(C)

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
AMPL_CPPFLAGS="-I\$(top_builddir)/OTHERS/ASL -I\$(top_srcdir)/OTHERS/ASL/solvers"

#########################
# Makefile conditionals #
#########################

# The following variable collects the names of libraries that should be
# included into libipopt.a
AC_SUBST(IPALLLIBS)
IPALLLIBS="Algorithm/libipoptalg.a Apps/StdInterface/libstdinterface.a Common/libcommon.a LinAlgImpl/Serial/liblaserial.a PDSystemImpl/PDFullSpace/libpdfullspace.a"

AM_CONDITIONAL([BUILDBLAS],[test "$use_blas" = "BUILD"])
if test "$use_blas" = BUILD; then
  LIBS="$ipblaslibdir/libipblas.a $LIBS"
  IPALLLIBS="OTHERS/blas/libipblas.a $IPALLLIBS"
fi

AM_CONDITIONAL([BUILDLAPACK],[test "$use_lapack" = "BUILD"])
if test "$use_lapack" = BUILD; then
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
AM_CONDITIONAL([BUILDMC19],[test "$use_MC19" = BUILD])

AM_CONDITIONAL([BUILDBLASDIR],[test "$use_blas" = BUILD])
AM_CONDITIONAL([BUILDLAPACKDIR],[test "$use_lapack" = BUILD])
AM_CONDITIONAL([BUILDHSLDIR],[test "$use_hsl" = BUILD])

AM_CONDITIONAL([HAVE_PARDISO],[test ! "$use_pardiso" = no])

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


AC_CONFIG_FILES([Makefile Common/Makefile Algorithm/Makefile
                 LinAlgImpl/Makefile LinAlgImpl/Serial/Makefile 
                 PDSystemImpl/Makefile PDSystemImpl/PDFullSpace/Makefile
		 Apps/Makefile Apps/AmplSolver/Makefile
		 Apps/StdInterface/Makefile
		 Apps/CompositeInterface/Makefile
		 Apps/CUTErInterface/Makefile
		 Examples/Fortran_example/Makefile
		 Examples/Cpp_example/Makefile
		 OTHERS/Makefile OTHERS/blas/Makefile
		 OTHERS/lapack/Makefile OTHERS/HSL/Makefile
		 OTHERS/ASL/Makefile])

AC_CONFIG_HEADER([inc/config.h])

AC_OUTPUT

AC_MSG_NOTICE([Configuration successful])
