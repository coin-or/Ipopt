# Copyright (C) 2004, 2005 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = Extern Common LinAlg Algorithm Interfaces Apps

IPALLLIBS = @IPALLLIBS@

EXTRA_DIST = \
	Examples/Cpp_example/cpp_example.cpp \
	Examples/Cpp_example/MyNLP.cpp \
	Examples/Cpp_example/MyNLP.hpp \
	Examples/hs071_f/hs071_f.f \
	Examples/hs071_c/hs071_c.c \
	Examples/hs071_cpp/hs071_main.cpp \
	Examples/hs071_cpp/hs071_nlp.cpp \
	Examples/hs071_cpp/hs071_nlp.hpp \
	README \
	Docs/documentation.tex \
	Windows/VisualStudio_dotNET/README \
	Windows/VisualStudio_dotNET/CppExample/CppExample.vcproj \
	Windows/VisualStudio_dotNET/FortranExample/FortranExample.vfproj \
	Windows/VisualStudio_dotNET/FortranLinAlg/FortranLinAlg.vfproj \
	Windows/VisualStudio_dotNET/include/config.h \
	Windows/VisualStudio_dotNET/Ipopt/Ipopt.sln \
	Windows/VisualStudio_dotNET/Ipopt/Ipopt.vcproj

.PHONY: doc

doc:
	doxygen

tmpdir = tmpdir

libipopt = libipopt.a

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
	$(AR) $(ARFLAGS) $(libipopt) *.$(OBJEXT)
	$(INSTALL_DATA) $(tmpdir)/$(libipopt) $(DESTDIR)$(libdir)
	rm -rf $(tmpdir)
	$(RANLIB) $(DESTDIR)$(libdir)/$(libipopt)

uninstall-local:
	rm -f $(DESTDIR)$(libdir)/$(libipopt)

cuter: install
	cd Apps/CUTErInterface; $(MAKE) cuter

astyle:
	cd Algorithm; $(MAKE) astyle
	cd Algorithm/LinearSolvers; $(MAKE) astyle
	cd Apps/AmplSolver; $(MAKE) astyle
	cd Common; $(MAKE) astyle
	cd LinAlg; $(MAKE) astyle
	cd LinAlg/TMatrices; $(MAKE) astyle
	cd Interfaces; $(MAKE) astyle

