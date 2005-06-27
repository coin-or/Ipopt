# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = OTHERS Common LinAlgImpl Algorithm PDSystemImpl Apps

IPALLLIBS = @IPALLLIBS@

EXTRA_DIST = Interfaces/IpEqMultCalculator.hpp \
	Interfaces/IpMatrix.hpp \
	Interfaces/IpNLP.hpp \
	Interfaces/IpPDSystemSolver.hpp \
	Interfaces/IpSymMatrix.hpp \
	Interfaces/IpTNLP.hpp \
	Interfaces/IpVector.hpp \
	Interfaces/IpStdCInterface.h \
	Examples/Fortran_example/example.f \
	README

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
	  ar x $$basdir/$$lib; \
	done; \
	$(AR) $(ARFLAGS) $(libipopt) *.$(OBJEXT)
	$(INSTALL_DATA) $(tmpdir)/$(libipopt) $(DESTDIR)$(libdir)
	rm -rf $(tmpdir)
	$(RANLIB) $(DESTDIR)$(libdir)/$(libipopt)

uninstall-local:
	rm -f $(DESTDIR)$(libdir)/$(libipopt)

cuter: install
	cd Apps/CUTErInterface; $(MAKE) cuter

# Astyle stuff

ASTYLE_FILES = Interfaces/IpEqMultCalculator.hppbak \
	Interfaces/IpMatrix.hppbak \
	Interfaces/IpNLP.hppbak \
	Interfaces/IpPDSystemSolver.hppbak \
	Interfaces/IpSymMatrix.hppbak \
	Interfaces/IpTNLP.hppbak \
	Interfaces/IpVector.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES = $(ASTYLE_FILES)

SUFFIXES = .cppbak .hppbak

astyle: astyle-local
	cd Algorithm; $(MAKE) astyle
	cd Apps/AmplSolver; $(MAKE) astyle
	cd Apps/CompositeInterface; $(MAKE) astyle
	cd Apps/StdInterface; $(MAKE) astyle
	cd Common; $(MAKE) astyle
	cd LinAlgImpl/Serial; $(MAKE) astyle
	cd PDSystemImpl/PDFullSpace; $(MAKE) astyle

astyle-local: $(ASTYLE_FILES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@
