# Copyright (C) 2004, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

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

.PHONY: cuter test unitTest astyle doxydoc

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
	
doxydoc :
	doxygen doxydoc/doxygen.conf
	
clean-doxydoc:
	cd doxydoc && rm -rf html *.log *.tag

install-exec-local: install-doc

uninstall-local: uninstall-doc

clean-local : clean-doxydoc

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

# setup addlibs files
# the addlibs_cpp file is created by running pkg-config on the installed ipopt.pc, if pkg-config available,
#   or by using variables setup in configure
# the addlibs_c file is created by adding @CXXLIBS@ to the content of the addlibs_cpp file
# the addlibs_f file is creating by taking the content of the addlibs_cpp file and removing the FLIBS
#   since pkg-config may have reordered the flags from FLIBS, we take them out one by one
install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
if COIN_HAS_PKGCONFIG
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@:$(DESTDIR)$(pkgconfiglibdir) \
	$(PKG_CONFIG) --libs ipopt > $(addlibsdir)/ipopt_addlibs_cpp.txt
else
if COIN_CXX_IS_CL
	echo "-libpath:`$(CYGPATH_W) @abs_lib_dir@` libipopt.lib @IPOPTLIB_LIBS_INSTALLED@" > $(addlibsdir)/ipopt_addlibs_cpp.txt
else
	echo "-L@abs_lib_dir@ -lipopt @IPOPTLIB_LIBS_INSTALLED@" > $(addlibsdir)/ipopt_addlibs_cpp.txt
endif
endif
	addlibs=`cat $(addlibsdir)/ipopt_addlibs_cpp.txt` ; \
	echo "$$addlibs @CXXLIBS@" > $(addlibsdir)/ipopt_addlibs_c.txt ; \
	for i in @FLIBS@ coin_dummy ; do \
		addlibs=`echo -n " $$addlibs " | sed -e "s! $$i ! !g"` ; \
	done ; \
	echo "$$addlibs @CXXLIBS@" > $(addlibsdir)/ipopt_addlibs_f.txt

uninstall-hook:
	rm -f $(addlibsdir)/ipopt_addlibs_{cpp,c,f}.txt

include BuildTools/Makemain.inc
