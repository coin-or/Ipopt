# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

includeipoptdir = $(includedir)/ipopt
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
	IpUtils.hpp

noinst_LIBRARIES = libcommon.a

libcommon_a_SOURCES = \
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

CLEANFILES = $(ASTYLE_FILES)

SUFFIXES = .cppbak .hppbak

install-exec-local:
	echo "#ifndef __IPCONFIG__" >bla
	echo "#define __IPCONFIG__" >>bla
	cat ../inc/config.h >> bla
	echo "#endif" >>bla
	sed -e 's/#include "config.h"//' <$(srcdir)/IpTypes.hpp >> bla
	$(install_sh_DATA) bla $(includeipoptdir)/IpTypes.hpp
	rm bla

uninstall-local:
	rm -f $(includeipoptdir)/IpTypes.hpp

astyle: $(ASTYLE_FILES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@
