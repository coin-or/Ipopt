# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = Common LinAlgImpl Algorithm PDSystemImpl Apps

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

IPALLLIBS = Algorithm/libipoptalg.a \
	Apps/StdInterface/libstdinterface.a \
	Common/libcommon.a \
	LinAlgImpl/Serial/liblaserial.a \
	PDSystemImpl/PDFullSpace/libpdfullspace.a

install-exec-local:
	set -x; \
	$(mkinstalldirs) $(DESTDIR)$(libdir); \
	rm -rf $(tmpdir); \
	mkdir $(tmpdir); \
	basdir=`pwd`; \
	cd $(tmpdir); \
	for lib in $(IPALLLIBS); \
	do \
	  $(AR) x $$basdir/$$lib; \
	done; \
	$(AR) cur $(libipopt) *.$(OBJEXT)
	$(INSTALL_DATA) $(tmpdir)/$(libipopt) $(DESTDIR)$(libdir)
	rm -rf $(tmpdir)
	$(RANLIB) $(DESTDIR)$(libdir)/$(libipopt)

uninstall-local:
	rm -f $(DESTDIR)$(libdir)/$(libipopt)


# Astyle stuff

ASTYLE_SOURCES = Interfaces/IpEqMultCalculator.hppbak \
	Interfaces/IpMatrix.hppbak \
	Interfaces/IpNLP.hppbak \
	Interfaces/IpPDSystemSolver.hppbak \
	Interfaces/IpSymMatrix.hppbak \
	Interfaces/IpTNLP.hppbak \
	Interfaces/IpVector.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES = $(ASTYLE_SOURCES)

SUFFIXES = .cppbak .hppbak

astyle: $(ASTYLE_SOURCES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@
