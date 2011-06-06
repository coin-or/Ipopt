# Copyright (C) 2003, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = \
	CUTErInterface.f \
	ipoma-dummy.f

noinst_LTLIBRARIES = libdummy.la

nodist_libdummy_la_SOURCES = ipoma.f

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Common

ipoma.f: CUTErInterface.f
	if test -z "$$MYCUTER" -o -z "$$CUTER"; then \
	  rm -f  ipoma.f; \
	  cp $(srcdir)/ipoma-dummy.f ipoma.f; \
	else \
	  sed "/#{sizing}/r $$MYCUTER/double/config/tools.siz" $< | sed -f $$MYCUTER/double/config/cast.sed | sed "s/C$(BITS_PER_POINTER)BIT//" > ipoma.f; \
	fi

ipo:
	if test -r $$CUTER/build/prototypes/gen.pro ; then \
	  sed -e 's/gen/ipo/g' $$CUTER/build/prototypes/gen.pro \
	    -e 's/"ipo.o"/""/' \
	    -e 's!SYSLIBS ""!SYSLIBS "\`cat $${MYCUTER}/double/lib/ipo.liblist\`"!'\
            -e "s!set arguments = ''!set arguments = \"--blas -lm --lapack -lm\"!" > ipo.pro; \
	  sed -f $$MYCUTER/double/config/script.sed ipo.pro > $@; \
	  chmod a+x $@; \
	else \
	  touch ipo; \
	fi

sdipo:
	if test -r $$CUTER/build/prototypes/gen.pro ; then \
	  sed -e 's/gen/ipo/g' $$CUTER/build/prototypes/sdgen.pro > sdipo.pro; \
	  sed -f $$MYCUTER/double/config/script.sed sdipo.pro > $@ ; \
	  chmod a+x $@; \
	else \
	  touch sdipo; \
	fi

CLEANFILES = ipoma.f ipo sdipo ipo.pro sdipo.pro

cuter:	ipoma.$(OBJEXT) ipo sdipo
	if test -z "$$MYCUTER" -o -z "$$CUTER"; then \
	  echo CUTEr environment variables not set - cannot install ipoma.f; \
	else \
	  echo "-L$(libdir) -lipopt $(IPOPTLIB_LIBS_INSTALLED) $(CXXLIBS) $(RPATH_FLAGS)" >$$MYCUTER/double/lib/ipo.liblist ;\
	  if test -r $$CUTER/build/prototypes/gen.pro ; then \
	    cp ipoma.o $(MYCUTER)/double/bin/ipoma.$(OBJEXT); \
	    cp ipo $(MYCUTER)/bin/ipo ;\
	    cp sdipo $(MYCUTER)/bin/sdipo ;\
	  else \
	    cp ipoma.o $(MYCUTER)/double/bin/ipoptma.$(OBJEXT); \
	  fi; \
	fi
