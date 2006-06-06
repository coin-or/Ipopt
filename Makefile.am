# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 56 2006-06-01 15:53:57Z andreasw $

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

test:
	cd Ipopt; make test

unitTest: test

tests: all
	for dir in $(subdirs); do \
	  if test -r $$dir/test/Makefile; then \
	    (cd $$dir; make test) \
	  fi; \
	done

unitTests: tests

doxydoc:
	cd $(srcdir); doxygen doxydoc/doxygen.conf

.PHONY: test unitTest tests unitTests doxydoc

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
