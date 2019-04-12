# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 672 2005-03-12 20:44:24Z andreasw $

AUTOMAKE_OPTIONS = foreign

SUBDIRS = OTHERS IPOPT

IPALLLIBS = @IPALLLIBS@

IPADDLIBS = @IPADDLIBS@
FLIBS = @FLIBS@

EXTRA_DIST = EXAMPLES/Fortran_example/Makefile.in \
	EXAMPLES/Fortran_example/example.f \
	EXAMPLES/C_example/Makefile.in \
	EXAMPLES/C_example/example.c \
	README \
	doc/CHANGELOG \
	doc/INSTALL \
	doc/LICENSE \
	doc/PARAMS.DAT.example \
	doc/QUICKINSTALL \
	doc/README.IPOPT \
	doc/README.Windows \
	VisualStudio/Ipopt.sln \
	VisualStudio/include/Ipopt.h \
	VisualStudio/include/Ipopt_Interface.h \
	VisualStudio/include/config.h \
	VisualStudio/include/config_f.h \
	VisualStudio/Fortran_example/Fortran_example.vfproj \
	VisualStudio/Ipopt/Ipopt.vfproj \
	VisualStudio/Ipopt-C/Ipopt-C.vcproj \
	VisualStudio/C_example/C_example.vcproj

tmpdir = tmpdir

libipopt = libipopt.a

ipopt_addlibs_c_txt = ipopt_addlibs_c.txt
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
	  $(AR) x $$basdir/$$lib; \
	done; \
	$(AR) cur $(libipopt) *.$(OBJEXT)
	$(INSTALL_DATA) $(tmpdir)/$(libipopt) $(DESTDIR)$(libdir)
	echo "$(IPADDLIBS) $(FLIBS)" >$(tmpdir)/$(ipopt_addlibs_c_txt)
	echo "$(IPADDLIBS)" >$(tmpdir)/$(ipopt_addlibs_f_txt)
	$(INSTALL_DATA) $(tmpdir)/$(ipopt_addlibs_c_txt) $(DESTDIR)$(libdir)
	$(INSTALL_DATA) $(tmpdir)/$(ipopt_addlibs_f_txt) $(DESTDIR)$(libdir)
	rm -rf $(tmpdir)
	$(RANLIB) $(DESTDIR)$(libdir)/$(libipopt)

uninstall-local:
	rm -f $(DESTDIR)$(libdir)/$(libipopt)
	rm -f $(DESTDIR)$(libdir)/$(ipopt_addlibs_c_txt)
	rm -f $(DESTDIR)$(libdir)/$(ipopt_addlibs_f_txt)

cuter: install
	cd IPOPT/CUTEr_interface; $(MAKE) cuter
