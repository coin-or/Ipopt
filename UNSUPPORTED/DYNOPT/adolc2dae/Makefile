# Copyright (C) 2002, Carnegie Mellon University and others.
# All Rights Reserved.
# This code is published under the Common Public License.

# $Id: Makefile 531 2004-03-11 01:31:07Z andreasw $
THIS_DIST_DIR   := DYNOPT/adolc2dae
#
#  The next lines automatically determine the root dir for the objects
#
slash := /
empty :=
space := $(empty) $(empty)
dirs  := $(subst /, ,$(THIS_DIST_DIR))
dots  := $(foreach dir,$(dirs), .. )
OBJROOTDIR := $(subst $(space),$(slash),$(strip $(dots)))

MAKEFILEDEFS = $(OBJROOTDIR)/Makefile.defs
MAKEFILETRGS = $(OBJROOTDIR)/Makefile.trgs

include $(MAKEFILEDEFS)

objs = \
	cprogs.o \
	fprogs.o \
	model_init.o \
	model_init_add.o

f2cpost     = f2c_post
f2cpostmain = $(F2CPOST).o

CLEANOBJS = $(objs) $(f2cpostmain)
CLEANLIBS = 
CLEANEXES = $(f2cpost)

all: $(f2cpost)

$(f2cpost): $(f2cpostmain)
	$(CC) $(CFLAGS) -o $@ $(f2cpostmain) $(LNLIBS)

#
#  Dependencies
#
cprogs.o \
model_init.o \
model_init_add.o: f2c_adolc.h

incdirs     = $(ADOLCSRCDIR)
CXXINCFLAGS = $(patsubst %, -I%, $(incdirs))


include $(MAKEFILETRGS)
