# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libdummy.a

libdummy_a_SOURCES = ipoma.f

ipoma.f: CUTErInterface.f
	if test -z "$$MYCUTER" -o -z "$$CUTER"; then \
	  touch ipoma.f; \
	else \
	  sed "/#{sizing}/r $$MYCUTER/double/config/tools.siz" $< | sed -f $$MYCUTER/double/config/cast.sed > ipoma.f; \
	fi

ipo: ipo.pro
	sed -f $$MYCUTER/double/config/script.sed $< > $@
	chmod a+x $@

sdipo: sdipo.pro
	sed -f $$MYCUTER/double/config/script.sed $< > $@
	chmod a+x $@

IPADDLIBS = @IPADDLIBS@

CLEANFILES = ipoma.f ipo sdipo

cuter:	ipoma.o ipo sdipo
	if test -z "$$MYCUTER" -o -z "$$CUTER"; then \
	  echo CUTEr environment variables not set - cannot install ipoma.f; \
	else \
	  echo "-L$(libdir) -lipopt $(IPADDLIBS)" >$$MYCUTER/double/lib/ipo.liblist ;\
	  cp ipoma.o $(MYCUTER)/double/bin/ipoma.o; \
	  cp ipo $(MYCUTER)/bin/ipo ;\
	  cp sdipo $(MYCUTER)/bin/sdipo ;\
	fi
