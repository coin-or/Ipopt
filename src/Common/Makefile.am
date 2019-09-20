# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

includeipoptdir = $(includedir)/coin-or
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
	IpDebug.cpp \
	IpJournalist.cpp \
	IpObserver.cpp \
	IpOptionsList.cpp \
	IpRegOptions.cpp \
	IpTaggedObject.cpp \
	IpUtils.cpp

#libcommon_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPT_BUILD

install-exec-local:
	$(install_sh_DATA) config_ipopt.h $(DESTDIR)$(includeipoptdir)/IpoptConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includeipoptdir)/IpoptConfig.h
