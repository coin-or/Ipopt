# Copyright (C) 2004, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

SUBDIRS = src/Common src/LinAlg src/Algorithm src/contrib/CGPenalty src/Interfaces src/Apps

if BUILD_SIPOPT
  SUBDIRS += contrib/sIPOPT
endif

doc_DATA = README.md AUTHORS LICENSE

.PHONY: test unitTest doc javadoc astyle cppcheck cppcheck-config

test: unitTest

unitTest: all
	cd test; $(MAKE) test

doc :
	cd doc && doxygen

if BUILD_JAVA
javadoc :
	$(JAVADOC) -d javadoc -windowtitle "JIpopt API documentation" $(srcdir)/src/Interfaces/Ipopt.java
endif

clean-doc:
	cd doc && rm -rf html *.log *.tag
	rm -rf javadoc

clean-local : clean-doc

astyle:
	cd $(srcdir) && astyle --mode=c -A1 --indent=spaces=3 --indent-switches --min-conditional-indent=1 --convert-tabs --align-pointer=type --pad-oper --add-brackets -n -r "*.hpp" "*.h"  "*.cpp" "*.c"
	cd $(srcdir) && astyle --mode=java -A1 --indent=spaces=3 --indent-switches --min-conditional-indent=1 --convert-tabs --align-pointer=type --pad-oper --add-brackets -n -r "*.java"

cppcheck:
	cppcheck -q -j4 --enable=all --suppress=uninitMemberVar --suppress=noExplicitConstructor --suppress=cstyleCast:src/Apps/AmplSolver/AmplTNLP.cpp --inline-suppr \
	  --addon=threadsafety.py --std=c++03 --std=c99 --error-exitcode=2 --relative-paths=$(srcdir) \
	  $(IPOPTLIB_CFLAGS) $(IPOPTAMPLINTERFACELIB_CFLAGS) $(addprefix -I,$(shell find $(srcdir)/src -type d)) -Isrc/Common \
	  -DHAVE_CONFIG_H $(srcdir)/src $(srcdir)/test

cppcheck-config:
	cppcheck -q --check-config --suppress=missingIncludeSystem --std=c++03 --std=c99 --error-exitcode=2 \
	  $(IPOPTLIB_CFLAGS) $(IPOPTAMPLINTERFACELIB_CFLAGS) $(addprefix -I,$(shell find $(srcdir)/src -type d)) -Isrc/Common \
	  -DHAVE_CONFIG_H $(srcdir)/src $(srcdir)/test

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = ipopt.pc
