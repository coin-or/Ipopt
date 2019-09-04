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

org/coinor/Ipopt.class org_coinor_Ipopt.h : $(srcdir)/Ipopt.java
	$(JAVAC) $< -d $(PWD) -h $(PWD) -target 1.8 -source 1.8

all-local : org.coinor.ipopt.jar

CLEANFILES = org.coinor.ipopt.jar org/coinor/Ipopt.class org_coinor_Ipopt.h

install-exec-local :
	$(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL) $(INSTALL_STRIP_FLAG) org.coinor.ipopt.jar "$(DESTDIR)$(libdir)"
uninstall-local :
	$(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=uninstall rm -f "$(DESTDIR)$(libdir)/org.coinor.ipopt.jar"

endif

libipopt_la_LIBADD = $(IPALLLIBS) $(IPOPTLIB_LFLAGS)

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

# Astyle stuff

ASTYLE_FILES = \
	IpAlgTypes.hppbak \
	IpInterfacesRegOp.cppbak IpInterfacesRegOp.hppbak \
	IpIpoptApplication.cppbak IpIpoptApplication.hppbak \
	IpNLP.hppbak \
	IpReturnCodes_inc.hbak IpReturnCodes.hbak IpReturnCodes.hppbak \
	IpSolveStatistics.cppbak IpSolveStatistics.hppbak \
	IpStdCInterface.cppbak IpStdCInterface.hbak \
	IpStdFInterface.cbak \
	IpStdInterfaceTNLP.cppbak IpStdInterfaceTNLP.hppbak \
	IpTNLP.hppbak \
	IpTNLPAdapter.cppbak IpTNLPAdapter.hppbak \
	IpTNLPReducer.cppbak IpTNLPReducer.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

DISTCLEANFILES = $(ASTYLE_FILES)

SUFFIXES = .cppbak .hppbak .cbak .hbak

astyle: $(ASTYLE_FILES)

.h.hbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ | sed -e 's/  *$$//' > $<
	touch $@

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ | sed -e 's/  *$$//' > $<
	touch $@

.c.cbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ | sed -e 's/  *$$//' > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ | sed -e 's/  *$$//' > $<
	touch $@
