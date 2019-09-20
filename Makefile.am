# Copyright (C) 2004, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = src/Common src/LinAlg src/Algorithm src/contrib/CGPenalty

if BUILD_LINEARSOLVERLOADER
  SUBDIRS += src/contrib/LinearSolverLoader
endif

SUBDIRS += src/Interfaces src/Apps

if BUILD_SIPOPT
  SUBDIRS += contrib/sIPOPT
endif

if ALWAYS_FALSE
  SUBDIRS += test
endif

doc_DATA = README.md AUTHORS LICENSE

.PHONY: test unitTest astyle doc javadoc

# Stuff to clean up from "make test"

DISTCLEANFILES = \
	examples/Cpp_example/cpp_example@EXEEXT@ \
	examples/Cpp_example/cpp_example.@OBJEXT@ \
	examples/Cpp_example/MyNLP.@OBJEXT@ \
	examples/hs071_cpp/hs071_cpp@EXEEXT@ \
	examples/hs071_cpp/hs071_main.@OBJEXT@ \
	examples/hs071_cpp/hs071_nlp.@OBJEXT@ \
	examples/hs071_cpp/ipopt.out \
	examples/hs071_c/hs071_c@EXEEXT@ \
	examples/hs071_c/hs071_c.@OBJEXT@ \
	examples/hs071_c/ipopt.out \
	examples/hs071_f/hs071_f@EXEEXT@ \
	examples/hs071_f/hs071_f.@OBJEXT@ \
	examples/hs071_f/IPOPT.OUT \
	examples/ScalableProblems/*.@OBJEXT@ \
	examples/ScalableProblems/solve_problem@EXEEXT@ \
	examples/ScalableProblems/solution.txt \
	tutorial/CodingExercise/*/*/*.@OBJEXT@ \
	tutorial/CodingExercise/C/*/TutorialC@EXEEXT@ \
	tutorial/CodingExercise/Cpp/*/TutorialCpp@EXEEXT@ \
	tutorial/CodingExercise/Fortran/*/TutorialFortran@EXEEXT@

tmpdir = tmpdir

libipopt = libipopt.a

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

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = ipopt.pc
