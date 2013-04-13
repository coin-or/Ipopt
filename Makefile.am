# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = doxydoc/doxygen.conf

########################################################################
#                          Subdirectories                              #
########################################################################

# subdirs is set by configure as the list of all subdirectories to recurse
# into
SUBDIRS = $(subdirs)

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd Ipopt; $(MAKE) test

unitTest: test

tests: all
	for dir in $(subdirs); do \
	  if test -r $$dir/test/Makefile; then \
	    (cd $$dir; $(MAKE) test) \
	  fi; \
	done

unitTests: tests

doxydoc:
	doxygen doxydoc/doxygen.conf

clean-doxydoc :
	cd doxydoc && rm -rf html *.log *.tag
	
install-doxydoc : doxydoc
	$(mkdir_p) $(DocInstallDir)/doxydoc
	cp -R doxydoc/html doxydoc/ipopt_doxy.tag $(DocInstallDir)/doxydoc

uninstall-doxydoc:
	if test -e $(DocInstallDir)/doxydoc ; then rm -rf $(DocInstallDir)/doxydoc ; fi

cuter:
	cd Ipopt; $(MAKE) cuter

.PHONY: test unitTest tests unitTests doxydoc cuter

clean-local: clean-doxydoc
uninstall-local: uninstall-doxydoc
install-exec-local:
	@if test -e doxydoc/html ; then \
	  $(MAKE) install-doxydoc ; \
	else \
	  echo "Run make install-doxydoc to generate and install Doxygen documentation." ; \
	fi

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES = coin_subdirs.txt

include BuildTools/Makemain.inc
