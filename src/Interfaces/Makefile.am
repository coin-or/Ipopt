# Copyright (C) 2004, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

includeipoptdir = $(includedir)/coin-or
includeipopt_HEADERS = \
	IpAlgTypes.hpp \
	IpIpoptApplication.hpp \
	IpNLP.hpp \
	IpReturnCodes.h \
	IpReturnCodes.hpp \
	IpReturnCodes_inc.h \
	IpReturnCodes.inc \
	IpSolveStatistics.hpp \
	IpStdCInterface.h \
	IpTNLP.hpp \
	IpTNLPAdapter.hpp \
	IpTNLPReducer.hpp

lib_LTLIBRARIES = libipopt.la

libipopt_la_SOURCES = \
	IpInterfacesRegOp.cpp \
	IpIpoptApplication.cpp \
	IpSolveStatistics.cpp \
	IpStdCInterface.cpp \
	IpStdInterfaceTNLP.cpp \
	IpStdFInterface.c \
	IpTNLPAdapter.cpp \
	IpTNLPReducer.cpp

if BUILD_JAVA
libipopt_la_SOURCES += IpStdJInterface.cpp org_coinor_Ipopt.h
BUILT_SOURCES = org_coinor_Ipopt.h

org.coinor.ipopt.jar : org/coinor/Ipopt.class
	$(JAR) cf $@ $<

# had  -target 1.8 -source 1.8, but that produces a warning if Java 8 RT is not setup
org/coinor/Ipopt.class org_coinor_Ipopt.h : $(srcdir)/Ipopt.java
	$(JAVAC) $< -d $(PWD) -h $(PWD)

all-local : org.coinor.ipopt.jar

CLEANFILES = org.coinor.ipopt.jar org/coinor/Ipopt.class org_coinor_Ipopt.h

install-exec-local :
	$(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL) $(INSTALL_STRIP_FLAG) org.coinor.ipopt.jar "$(DESTDIR)$(libdir)"
uninstall-local :
	$(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=uninstall rm -f "$(DESTDIR)$(libdir)/org.coinor.ipopt.jar"

endif

libipopt_la_LIBADD = $(IPALLLIBS) $(IPOPTLIB_LFLAGS)
# probably we need to make this explicit because we set use a variable $(IPALLLIBS) in LIBADD
libipopt_la_DEPENDENCIES = $(IPALLLIBS)

libipopt_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPT_BUILD \
	-I$(srcdir)/../Common \
	-I$(srcdir)/../LinAlg \
	-I$(srcdir)/../Algorithm/LinearSolvers \
	-I$(srcdir)/../LinAlg/TMatrices \
	-I$(srcdir)/../Algorithm \
	-I$(srcdir)/../Algorithm/Inexact \
	-I$(srcdir)/../contrib/CGPenalty \
	-I$(srcdir)/../contrib/LinearSolverLoader \
	$(IPOPTLIB_CFLAGS)
