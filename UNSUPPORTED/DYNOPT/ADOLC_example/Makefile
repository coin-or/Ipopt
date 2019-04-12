# Copyright (C) 2002, Carnegie Mellon University and others.
# All Rights Reserved.
# This code is published under the Common Public License.

# $Id: Makefile 531 2004-03-11 01:31:07Z andreasw $
THIS_DIST_DIR   := DYNOPT/ADOLC_example
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

model    = tank
exe      = $(model)_ADOLC
execheck = $(exe)_check
modelinc = $(datadir)/TANK_CONSTS.INC

mainobj  = $(OBJROOTDIR)/$(DIST_DYNOPT_AUXIL)/ip2dae_main.o 
maincheckobj = $(OBJROOTDIR)/$(DIST_DYNOPT_AUXIL)/ip2dae_check_derivs.o

modelf   = $(model).f
modelc   = $(model).c
modelcxx = $(model).cc
modelobj = $(model).o
modelso  = model.so

hereobjs =  $(modelobj) eval.o

ipoptauxilobjfiles = \
	c_out.o \
	iter_out.o \
	eval_h_dummy.o \
	eval_hesslag_v_dummy.o
auxilobjs1 = $(addprefix  $(OBJROOTDIR)/$(DIST_IPOPT_AUXIL)/, \
	$(ipoptauxilobjfiles))

dynoptauxilobjfiles = \
	constr.o \
	daemodel_df.o \
	daemodel_df_struc.o \
	daemodel_f.o \
	daemodel_hv.o \
	get_iele4t.o \
	get_start.o
auxilobjs2 = $(addprefix  $(OBJROOTDIR)/$(DIST_DYNOPT_AUXIL)/, \
	$(dynoptauxilobjfiles))

dynoptadolc2daefiles = \
	cprogs.o \
	fprogs.o
auxilobjs3 = $(addprefix  $(OBJROOTDIR)/$(DIST_DYNOPT_ADOLC2DAE)/, \
	$(dynoptadolc2daefiles))

modelinitobj = $(OBJROOTDIR)/$(DIST_DYNOPT_ADOLC2DAE)/model_init.o
modelinitaddobj = $(OBJROOTDIR)/$(DIST_DYNOPT_ADOLC2DAE)/model_init_add.o

auxilobjs4 = $(modelinitobj) $(modelinitaddobj)

AUXILOBJS = $(MAINOBJ) $(MAINCHECKOBJ) $(auxilobjs1) $(auxilobjs2) \
	 $(auxilobjs3) $(auxilobjs4)

data    = $(notdir $(wildcard $(datadir)/*.DAT))
datadir = $(SRCROOTDIR)/$(DIST_DYNOPT_EXAMPLEDATA)

CLEANOBJS = $(hereobjs) $(mainobj) $(maincheckobj) $(modelc) $(modelcxx)
CLEANEXES = $(exe) $(execheck) $(modelso)
CLEANDATA = $(data) IPOPT.OUT X.sol

objs    = $(hereobjs) $(auxilobjs1) $(auxilobjs2) $(auxilobjs3)
LIBS    = $(IPOPTLIB) $(TRONLIB) $(IPOPTLIB) $(DAE2NLPLIB) $(LINALGLIBS) $(ADOLCLIB)

#ifdef SOLIBDIR
#  depso = $(modelso)
#  objs += $(modelinitaddobj)
#else
  objs += $(modelinitobj)
#endif

all: $(exe) $(execheck) $(data)

$(exe): $(mainobj) $(objs) $(DEPLIBS) $(depso)
	$(CXX) $(LNFLAGS) -o $@ $(mainobj) $(objs) $(LNKLIBS) $(LNLIBS) \
	$(CXXLNLIBS)

$(execheck): $(maincheckobj) $(objs) $(DEPLIBS) $(depso)
	$(CXX) $(LNFLAGS) -o $@ $(maincheckobj) $(objs) $(LNKLIBS) $(LNLIBS) \
	$(CXXLNLIBS)

$(modelso): $(modelinitobj) $(modelobj)
	$(SOAR) $(SOARFLAGS) $@ $(modelinitobj) $(modelobj)

$(data):
	cp $(datadir)/$@ .
#
#  Special Building rules using f2c to obtain C++ model from Fortran
#
.INTERMEDIATE: $(modelc)
$(modelc): $(modelf) $(modelinc)
	$(F2C) $(F2CFLAGS) $(FINCFLAGS) $<

$(modelcxx): $(modelc) $(F2CPOST)
	$(F2CPOST) $< $@

$(modelobj): $(modelcxx) $(f2cadolc_h)
	$(CXX) -c $(CXXFLAGS) $(CXXINCFLAGS) $< -o $@

#
#  Mainly for development:
#  Force check of libraries and object files in other directories
#
$(patsubst %.a,%.so, $(DAE2NLPLIB) $(IPOPTLIB)) \
$(AUXILOBJS) : FORCE

#
#  Dependencies
#
fincdirs   = $(SRCROOTDIR)/$(DIST_DYNOPT_AUXIL) $(datadir)
FINCFLAGS  = $(patsubst %, -I%, $(fincdirs))

eval.o: $(SRCROOTDIR)/$(DIST_DYNOPT_AUXIL)/DYNAUX.INC

f2cadol_hdir = $(SRCROOTDIR)/$(DIST_DYNOPT_ADOLC2DAE)
cxxincdirs   = $(f2cadol_hdir) $(ADOLCSRCDIR)
CXXINCFLAGS  = $(patsubst %, -I%, $(cxxincdirs))

f2cadolc_h = $(f2cadol_hdir)/f2c_adolc.h
$(modelcc): $(f2cadolc_h)

#
#  Need to establish some order for parallelism in case FLAGS changed
#
break := $(firstword $(auxilobj1))
$(break): $(filter-out $(break), $(auxilobj1))

$(mainobj): $(auxilobj2)
$(maincheckobj) : $(mainobj)

$(modelcxx): $(ADOLCLIB)

include $(MAKEFILETRGS)
