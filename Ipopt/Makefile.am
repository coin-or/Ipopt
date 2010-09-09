# Copyright (C) 2004, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = src/Common src/LinAlg src/Algorithm src/contrib/CGPenalty

if BUILD_LINEARSOLVERLOADER
  SUBDIRS += src/contrib/LinearSolverLoader
endif

SUBDIRS += src/Interfaces src/Apps

if ALWAYS_FALSE
  SUBDIRS += test
endif

EXTRA_DIST = \
	examples/Cpp_example/cpp_example.cpp \
	examples/Cpp_example/Makefile.in \
	examples/Cpp_example/MyNLP.cpp \
	examples/Cpp_example/MyNLP.hpp \
	examples/hs071_f/Makefile.in \
	examples/hs071_f/hs071_f.f.in \
	examples/hs071_c/Makefile.in \
	examples/hs071_c/hs071_c.c \
	examples/hs071_cpp/Makefile.in \
	examples/hs071_cpp/hs071_main.cpp \
	examples/hs071_cpp/hs071_nlp.cpp \
	examples/hs071_cpp/hs071_nlp.hpp \
	examples/ScalableProblems/LuksanVlcek1.cpp \
	examples/ScalableProblems/LuksanVlcek1.hpp \
	examples/ScalableProblems/LuksanVlcek2.cpp \
	examples/ScalableProblems/LuksanVlcek2.hpp \
	examples/ScalableProblems/LuksanVlcek3.cpp \
	examples/ScalableProblems/LuksanVlcek3.hpp \
	examples/ScalableProblems/LuksanVlcek4.cpp \
	examples/ScalableProblems/LuksanVlcek4.hpp \
	examples/ScalableProblems/LuksanVlcek5.cpp \
	examples/ScalableProblems/LuksanVlcek5.hpp \
	examples/ScalableProblems/LuksanVlcek6.cpp \
	examples/ScalableProblems/LuksanVlcek6.hpp \
	examples/ScalableProblems/LuksanVlcek7.cpp \
	examples/ScalableProblems/LuksanVlcek7.hpp \
	examples/ScalableProblems/Makefile.in \
	examples/ScalableProblems/MittelmannBndryCntrlDiri.cpp \
	examples/ScalableProblems/MittelmannBndryCntrlDiri.hpp \
	examples/ScalableProblems/MittelmannBndryCntrlNeum.cpp \
	examples/ScalableProblems/MittelmannBndryCntrlNeum.hpp \
	examples/ScalableProblems/MittelmannDistCntrlDiri.cpp \
	examples/ScalableProblems/MittelmannDistCntrlDiri.hpp \
	examples/ScalableProblems/MittelmannDistCntrlNeumA.cpp \
	examples/ScalableProblems/MittelmannDistCntrlNeumA.hpp \
	examples/ScalableProblems/MittelmannDistCntrlNeumB.cpp \
	examples/ScalableProblems/MittelmannDistCntrlNeumB.hpp \
	examples/ScalableProblems/MittelmannParaCntrl.hpp \
	examples/ScalableProblems/README \
	examples/ScalableProblems/RegisteredTNLP.cpp \
	examples/ScalableProblems/RegisteredTNLP.hpp \
	examples/ScalableProblems/solve_problem.cpp \
	README \
	INSTALL \
	LICENSE \
	AUTHORS \
	doc/documentation.bbl \
	doc/documentation.tex \
	doc/documentation.pdf \
	doc/options.tex \
	test/run_unitTests.in \
	test/mytoy.nl \
	Windows/VisualStudio_dotNET/README \
	Windows/VisualStudio_dotNET/CppExample/CppExample.vcproj \
	Windows/VisualStudio_dotNET/hs071_cpp/hs071_cpp.vcproj \
	Windows/VisualStudio_dotNET/hs071_c/hs071_c.vcproj \
	Windows/VisualStudio_dotNET/hs071_f/hs071_f.vfproj \
	Windows/VisualStudio_dotNET/hs071_f/hs071_f.f \
	Windows/VisualStudio_dotNET/FortranLinAlg/FortranLinAlg.vfproj \
	Windows/VisualStudio_dotNET/Ipopt/Ipopt.sln \
	Windows/VisualStudio_dotNET/Ipopt/Ipopt.vcproj \
	Windows/VisualStudio_dotNET/ScalableProblems/ScalableProblems.vcproj

.PHONY: cuter test unitTest astyle

# Stuff to clean up from "make test"

DISTCLEANFILES = \
	examples/Cpp_example/cpp_example@EXEEXT@ \
	examples/Cpp_example/cpp_example.@OBJEXT@ \
	examples/Cpp_example/MyNLP.@OBJEXT@ \
	examples/hs071_cpp/hs071_cpp@EXEEXT@ \
	examples/hs071_cpp/hs071_main.@OBJEXT@ \
	examples/hs071_cpp/hs071_nlp.@OBJEXT@ \
	examples/hs071_cpp/ipopt.out \
	examples/hs071_c/hs071_c@EXEEXT@ \
	examples/hs071_c/hs071_c.@OBJEXT@ \
	examples/hs071_c/ipopt.out \
	examples/hs071_f/hs071_f@EXEEXT@ \
	examples/hs071_f/hs071_f.@OBJEXT@ \
	examples/hs071_f/IPOPT.OUT \
	examples/ScalableProblems/*.@OBJEXT@ \
	examples/ScalableProblems/solve_problem@EXEEXT@ \
	examples/ScalableProblems/solution.txt \
	tutorial/CodingExercise/*/*/*.@OBJEXT@ \
	tutorial/CodingExercise/C/*/TutorialC@EXEEXT@ \
	tutorial/CodingExercise/Cpp/*/TutorialCpp@EXEEXT@ \
	tutorial/CodingExercise/Fortran/*/TutorialFortran@EXEEXT@

tmpdir = tmpdir

libipopt = libipopt.a

cuter: install
	cd src/Apps/CUTErInterface; $(MAKE) cuter

test: unitTest

unitTest: all
	cd test; $(MAKE) test

install-exec-local: install-doc

uninstall-local: uninstall-doc

astyle:
	cd src/Algorithm; $(MAKE) astyle
	cd src/Algorithm/LinearSolvers; $(MAKE) astyle
	cd src/Algorithm/Inexact; $(MAKE) astyle
	cd src/contrib/CGPenalty; $(MAKE) astyle
	cd src/Apps/AmplSolver; $(MAKE) astyle
	cd src/Common; $(MAKE) astyle
	cd src/LinAlg; $(MAKE) astyle
	cd src/LinAlg/TMatrices; $(MAKE) astyle
	cd src/Interfaces; $(MAKE) astyle

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = ipopt.pc

addlibsdir = $(DESTDIR)$(datadir)/coin/doc/Ipopt

if COIN_HAS_PKGCONFIG
install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@ \
	$(PKG_CONFIG) --libs ipopt > $(addlibsdir)/ipopt_addlibs_cpp.txt
else
#temporary hack to get addlibs file also if no pkg-config is available (points to build directories instead of installation directories)
install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
	echo $(MUMPS_LIBS) $(HSL_LIBS) $(ADDLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) > $(addlibsdir)/ipopt_addlibs_cpp.txt
endif

uninstall-hook:
	rm -f $(addlibsdir)/ipopt_addlibs_cpp.txt

include BuildTools/Makemain.inc
