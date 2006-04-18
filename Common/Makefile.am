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
	IpTypes.hpp \
	IpUtils.hpp

nodist_includeipopt_HEADERS = ../inc/IpConfig.hpp

noinst_LTLIBRARIES = libcommon.la

libcommon_la_SOURCES = \
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

libcommon_la_LDFLAGS = -no-undefined

$(libcommon_la_OBJECTS): ../inc/IpConfig.hpp

# Here we automatically generate IpConfig.hpp so that we can distribute it
# without having unnecessary stuff in it (e.g., PACKAGE_NAME will conflict
# with another packages that uses autoconf)

../inc/IpConfig.hpp: ../inc/config_ipopt.h
	echo "/** Required defines from config_ipopt.h */" >IpConfig.tmp
	echo "#ifndef __IPCONFIG_HPP__" >>IpConfig.tmp
	echo "#define __IPCONFIG_HPP__" >>IpConfig.tmp
	$(EGREP) 'ASSERT|MATH|STDARG|STDLIB|TIME|IP_DEBUG' $< >>IpConfig.tmp
	echo " " >>IpConfig.tmp
	echo "/** Type of Fortran integer translated into C */" >>IpConfig.tmp
	echo "typedef `$(EGREP) 'FORTRAN_INTEGER_TYPE' $< | sed -e 's/#define FORTRAN_INTEGER_TYPE//'` ipfint;" >>IpConfig.tmp
	echo "#endif" >>IpConfig.tmp
	if test -r ../inc/IpConfig.hpp; then \
	  if diff ./inc/IpConfig.hpp IpConfig.tmp >/dev/null 2>1 ; then \
	    rm -f IpConfig.tmp; \
	  else \
	    rm -f ../inc/IpConfig.hpp; \
	    mv IpConfig.tmp ../inc/IpConfig.hpp; \
	  fi; \
	else \
	  mv IpConfig.tmp ../inc/IpConfig.hpp; \
	fi

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

CLEANFILES = $(ASTYLE_FILES) IpConfig.tmp

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
