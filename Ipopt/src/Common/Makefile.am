# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

includeipoptdir = $(includedir)/coin
includeipopt_HEADERS = \
	IpCachedResults.hpp \
	IpDebug.hpp \
	IpException.hpp \
	IpJournalist.hpp \
	IpObserver.hpp \
	IpOptionsList.hpp \
	IpReferenced.hpp \
	IpRegOptions.hpp \
	IpSmartPtr.hpp \
	IpTaggedObject.hpp \
	IpTimedTask.hpp \
	IpTypes.hpp \
	IpUtils.hpp

noinst_LTLIBRARIES = libcommon.la

libcommon_la_SOURCES = \
	IpoptConfig.h \
	IpCachedResults.hpp \
	IpDebug.cpp IpDebug.hpp \
	IpException.hpp \
	IpJournalist.cpp IpJournalist.hpp \
	IpObserver.cpp IpObserver.hpp \
	IpOptionsList.cpp IpOptionsList.hpp \
	IpRegOptions.cpp IpRegOptions.hpp \
	IpReferenced.hpp \
	IpSmartPtr.hpp \
	IpTaggedObject.cpp IpTaggedObject.hpp \
	IpTimedTask.hpp \
	IpTypes.hpp \
	IpUtils.cpp IpUtils.hpp

libcommon_la_LDFLAGS = $(LT_LDFLAGS)

install-exec-local:
	$(install_sh_DATA) config_ipopt.h $(DESTDIR)$(includeipoptdir)/IpoptConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includeipoptdir)/IpoptConfig.h

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

# Astyle stuff

ASTYLE_FILES = \
	IpCachedResults.hppbak \
	IpDebug.cppbak IpDebug.hppbak \
	IpException.hppbak \
	IpJournalist.cppbak IpJournalist.hppbak \
	IpObserver.cppbak IpObserver.hppbak \
	IpOptionsList.cppbak IpOptionsList.hppbak \
	IpRegOptions.cppbak IpRegOptions.hppbak \
	IpReferenced.hppbak \
	IpSmartPtr.hppbak \
	IpTaggedObject.cppbak IpTaggedObject.hppbak \
	IpTimedTask.hppbak \
	IpTypes.hppbak \
	IpUtils.cppbak IpUtils.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

DISTCLEANFILES = $(ASTYLE_FILES)

SUFFIXES = .cppbak .hppbak

astyle: $(ASTYLE_FILES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@
