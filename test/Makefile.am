# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author:  Andreas Waechter           IBM    2006-04-13

########################################################################
#                        unitTest for Ipopt                            #
########################################################################

noinst_PROGRAMS = hs071_cpp hs071_c emptynlp getcurr

if COIN_HAS_F77
noinst_PROGRAMS += hs071_f
endif

if BUILD_SIPOPT
noinst_PROGRAMS += parametric_cpp redhess_cpp
endif

nodist_hs071_cpp_SOURCES = hs071_main.cpp hs071_nlp.cpp hs071_nlp.hpp
hs071_cpp_LDADD = ../src/libipopt.la

nodist_hs071_c_SOURCES = hs071_c.c
hs071_c_LDADD = ../src/libipopt.la $(CXXLIBS)

nodist_emptynlp_SOURCES = emptynlp.cpp
emptynlp_LDADD = ../src/libipopt.la

nodist_getcurr_SOURCES = getcurr.cpp
getcurr_LDADD = ../src/libipopt.la

if !IPOPT_SINGLE
  nodist_hs071_f_SOURCES = hs071_f.f
else
  nodist_hs071_f_SOURCES = hs071_fs.f
endif
hs071_f_LDADD = ../src/libipopt.la $(CXXLIBS)

if BUILD_JAVA
if !IPOPT_SINGLE
  HS071J=HS071
else
  HS071J=HS071s
endif

$(HS071J).class : $(top_srcdir)/examples/hs071_java/$(HS071J).java
	$(JAVAC) -classpath ../src/org.coinor.ipopt.jar $< -d $(PWD)

all-local : $(HS071J).class

CLEANFILES = $(HS071J).class

endif

nodist_parametric_cpp_SOURCES = parametricTNLP.cpp parametric_driver.cpp
parametric_cpp_LDADD = ../contrib/sIPOPT/src/libsipopt.la

nodist_redhess_cpp_SOURCES = MySensTNLP.cpp redhess_cpp.cpp
redhess_cpp_LDADD = ../contrib/sIPOPT/src/libsipopt.la

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
  -I$(srcdir)/../src/Common \
  -I$(srcdir)/../src/LinAlg \
  -I$(srcdir)/../src/LinAlg/TMatrices \
  -I$(srcdir)/../src/Algorithm \
  -I$(srcdir)/../src/Interfaces \
  -I$(srcdir)/../contrib/sIPOPT/src \
  -I$(srcdir)/../contrib/sIPOPT/examples/parametric_cpp \
  -I$(srcdir)/../contrib/sIPOPT/examples/redhess_cpp

AM_FFLAGS = -I$(srcdir)/../src/Interfaces

hs071_f.f hs071_fs.f:
	$(LN_S) ../examples/hs071_f/$@ $@

if BUILD_JAVA
test: $(noinst_PROGRAMS) $(HS071J).class
else
test: $(noinst_PROGRAMS)
endif
	chmod u+x ./run_unitTests
	./run_unitTests

unitTest: test

.PHONY: test unitTest
