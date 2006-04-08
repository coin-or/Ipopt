# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = Extern Common LinAlg Algorithm Interfaces Apps

IPALLLIBS = @IPALLLIBS@

EXTRA_DIST = \
	Examples/Cpp_example/cpp_example.cpp \
	Examples/Cpp_example/Makefile.in \
	Examples/Cpp_example/MyNLP.cpp \
	Examples/Cpp_example/MyNLP.hpp \
	Examples/hs071_f/Makefile.in \
	Examples/hs071_f/hs071_f.f.in \
	Examples/hs071_c/Makefile.in \
	Examples/hs071_c/hs071_c.c \
	Examples/hs071_cpp/Makefile.in \
	Examples/hs071_cpp/hs071_main.cpp \
	Examples/hs071_cpp/hs071_nlp.cpp \
	Examples/hs071_cpp/hs071_nlp.hpp \
	Examples/ScalableProblems/LuksanVlcek1.cpp \
	Examples/ScalableProblems/LuksanVlcek1.hpp \
	Examples/ScalableProblems/LuksanVlcek2.cpp \
	Examples/ScalableProblems/LuksanVlcek2.hpp \
	Examples/ScalableProblems/LuksanVlcek3.cpp \
	Examples/ScalableProblems/LuksanVlcek3.hpp \
	Examples/ScalableProblems/LuksanVlcek4.cpp \
	Examples/ScalableProblems/LuksanVlcek4.hpp \
	Examples/ScalableProblems/LuksanVlcek5.cpp \
	Examples/ScalableProblems/LuksanVlcek5.hpp \
	Examples/ScalableProblems/LuksanVlcek6.cpp \
	Examples/ScalableProblems/LuksanVlcek6.hpp \
	Examples/ScalableProblems/LuksanVlcek7.cpp \
	Examples/ScalableProblems/LuksanVlcek7.hpp \
	Examples/ScalableProblems/Makefile.in \
	Examples/ScalableProblems/MittelmannBndryCntrlDiri.cpp \
	Examples/ScalableProblems/MittelmannBndryCntrlDiri.hpp \
	Examples/ScalableProblems/MittelmannBndryCntrlNeum.cpp \
	Examples/ScalableProblems/MittelmannBndryCntrlNeum.hpp \
	Examples/ScalableProblems/MittelmannDistCntrlDiri.cpp \
	Examples/ScalableProblems/MittelmannDistCntrlDiri.hpp \
	Examples/ScalableProblems/MittelmannDistCntrlNeumA.cpp \
	Examples/ScalableProblems/MittelmannDistCntrlNeumA.hpp \
	Examples/ScalableProblems/MittelmannDistCntrlNeumB.cpp \
	Examples/ScalableProblems/MittelmannDistCntrlNeumB.hpp \
	Examples/ScalableProblems/MittelmannParaCntrl.hpp \
	Examples/ScalableProblems/README \
	Examples/ScalableProblems/RegisteredTNLP.cpp \
	Examples/ScalableProblems/RegisteredTNLP.hpp \
	Examples/ScalableProblems/solve_problem.cpp \
	README \
	INSTALL \
	LICENSE \
	AUTHORS \
	Docs/documentation.bbl \
	Docs/documentation.tex \
	Docs/documentation.pdf \
	Docs/options.tex \
	Test/run_unitTests.in \
	Test/mytoy.nl \
	Windows/VisualStudio_dotNET/README \
	Windows/VisualStudio_dotNET/CppExample/CppExample.vcproj \
	Windows/VisualStudio_dotNET/hs071_cpp/hs071_cpp.vcproj \
	Windows/VisualStudio_dotNET/hs071_c/hs071_c.vcproj \
	Windows/VisualStudio_dotNET/hs071_f/hs071_f.vfproj \
	Windows/VisualStudio_dotNET/hs071_f/hs071_f.f \
	Windows/VisualStudio_dotNET/FortranLinAlg/FortranLinAlg.vfproj \
	Windows/VisualStudio_dotNET/include/config.h \
	Windows/VisualStudio_dotNET/Ipopt/Ipopt.sln \
	Windows/VisualStudio_dotNET/Ipopt/Ipopt.vcproj \
	Windows/VisualStudio_dotNET/ScalableProblems/ScalableProblems.vcproj


.PHONY: doc

# Stuff to clean up from "make test"

DISTCLEANFILES = \
	Examples/hs071_cpp/hs071_cpp \
	Examples/hs071_cpp/hs071_cpp.exe \
	Examples/hs071_cpp/hs071_main.o \
	Examples/hs071_cpp/hs071_main.obj \
	Examples/hs071_cpp/hs071_nlp.o \
	Examples/hs071_cpp/hs071_nlp.obj \
	Examples/hs071_cpp/ipopt.out \
	Examples/hs071_c/hs071_c \
	Examples/hs071_c/hs071_c.exe \
	Examples/hs071_c/hs071_c.o \
	Examples/hs071_c/hs071_c.obj \
	Examples/hs071_c/ipopt.out \
	Examples/hs071_f/hs071_f \
	Examples/hs071_f/hs071_f.exe \
	Examples/hs071_f/hs071_f.o \
	Examples/hs071_f/hs071_f.obj

doc:
	doxygen

tmpdir = tmpdir

libipopt = libipopt.a

ipopt_addlibs_cpp_txt = ipopt_addlibs_cpp.txt
ipopt_addlibs_f_txt = ipopt_addlibs_f.txt

install-exec-local:
	set -x; \
	$(mkinstalldirs) $(DESTDIR)$(libdir); \
	rm -rf $(tmpdir); \
	mkdir $(tmpdir); \
	basdir=`pwd`; \
	cd $(tmpdir); \
	for lib in $(IPALLLIBS); \
	do \
	  $(AR_X) $$basdir/$$lib; \
	done; \
	$(BLUEGENEEXTRA); \
	$(AR) $(ARFLAGS) $(libipopt) *.$(OBJEXT)
	$(INSTALL_DATA) $(tmpdir)/$(libipopt) $(DESTDIR)$(libdir)
	echo "$(IPADDLIBS) $(FLIBS)" >$(tmpdir)/$(ipopt_addlibs_cpp_txt)
	echo "$(IPADDLIBS)" >$(tmpdir)/$(ipopt_addlibs_f_txt)
	$(INSTALL_DATA) $(tmpdir)/$(ipopt_addlibs_cpp_txt) $(DESTDIR)$(libdir)
	$(INSTALL_DATA) $(tmpdir)/$(ipopt_addlibs_f_txt) $(DESTDIR)$(libdir)
	rm -rf $(tmpdir)
	$(RANLIB) $(DESTDIR)$(libdir)/$(libipopt)

uninstall-local:
	echo DESTDIR $(DESTDIR)
	echo libdir $(libdir)
	echo libipopt $(libipopt)
	rm -f $(DESTDIR)$(libdir)/$(libipopt)
	rm -f $(DESTDIR)$(libdir)/$(ipopt_addlibs_cpp_txt)
	rm -f $(DESTDIR)$(libdir)/$(ipopt_addlibs_f_txt)

cuter: install
	cd Apps/CUTErInterface; $(MAKE) cuter

test: unitTest

unitTests: unitTest

unitTest:
	cd Examples/hs071_cpp; make clean; make
	cd Examples/hs071_c; make clean; make
	cd Examples/hs071_f; make clean; make
	cd Test; chmod u+x ./run_unitTests; ./run_unitTests

astyle:
	cd Algorithm; $(MAKE) astyle
	cd Algorithm/LinearSolvers; $(MAKE) astyle
	cd Apps/AmplSolver; $(MAKE) astyle
	cd Common; $(MAKE) astyle
	cd LinAlg; $(MAKE) astyle
	cd LinAlg/TMatrices; $(MAKE) astyle
	cd Interfaces; $(MAKE) astyle

