# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

# $Id$

# Directory to the sources for the (example) problem definition files
SRCDIR = .


# Compiler options
CC = gcc
CFLAGS = -g

# Directory with libipopt.a
IPOPTLIBDIR = ${exec_prefix}/lib
exec_prefix = ${prefix}
prefix = /home/claird/devel/NewIpopt

# Libraries necessary to link with IPOPT
#LIBS = -L$(IPOPTLIBDIR) -lipopt -L/home/claird/lib -lhsl -llapack -lblas  -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2/../../.. -lstdc++ -lm -lgcc_s -lgcc
LIBS = -L$(IPOPTLIBDIR) -lipopt -L/home/claird/lib -lhsl -llapack -lblas  -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2/../../.. -lstdc++ -lm -lgcc_s -lgcc

# Necessary include directories
INCL = 	-I ${prefix}/Interfaces \
	-I ${prefix}/Common \
	-I ${prefix}/inc \
	-I ${prefix}/LinAlgImpl/Serial \
	-I ${prefix}/Algorithm

# CHANGEME: This should be the name of your executable
EXE = cpp_example

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS = 	MyNLP.o 	\
	cpp_example.o

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .cpp .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CC) $(CFLAGS) -o $@ $$bla $(LIBS)

clean:
	rm -rf $(EXE) $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<


.cpp.obj:
	$(CC) $(CFLAGS) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`
