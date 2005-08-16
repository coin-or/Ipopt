# Copyright (C) 2004, 2005 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = Extern Common LinAlg Algorithm Interfaces Apps

IPALLLIBS = @IPALLLIBS@

EXTRA_DIST = Examples/Fortran_example/example.f \
	README \
	Docs/documentation.tex

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
	cd Apps/CompositeInterface; $(MAKE) astyle
	cd Common; $(MAKE) astyle
	cd LinAlg; $(MAKE) astyle
	cd LinAlg/TMatrices; $(MAKE) astyle
	cd Interfaces; $(MAKE) astyle

