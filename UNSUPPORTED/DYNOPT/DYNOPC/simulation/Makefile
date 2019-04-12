# Copyright (C) 2002, Carnegie Mellon University and others.
# All Rights Reserved.
# This code is published under the Common Public License.

# $Id: Makefile 531 2004-03-11 01:31:07Z andreasw $
THIS_DIST_DIR   := DYNOPT/DYNOPC/simulation
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

## Uncomment following line to compile ESO version
#ESOVERSION := 1    # This is not part of the distribution

#
#  Need to redefine search path for source code
#
exe      = DYNOPC_simulation

ifndef PROCESS_DATA_DIR
  dummy := $(error You need to set the environment variable PROCESS_DATA_DIR for DYNOPC)
endif

# modellist := $(shell cd $(PROCESS_DATA_DIR); ls -d

ifndef ESOVERSION
  model    = $(MODEL)
  modelf   = $(PROCESS_DATA_DIR)/$(MODEL)/$(model).f
  modelc   = $(model).c
  modelcxx = $(model).cc
  modelobj = $(model).o
  ifdef SOLIBDIR
    modelso  = $(PROCESS_DATA_DIR)/$(MODEL)/$(model).so
  else
    modelso = $(error DYNOPC with ADOL-C only works if you specify SOLIBDIR.)
  endif
endif

mainobj  = $(OBJROOTDIR)/$(DIST_DYNOPC_UTILS)/main.o

hereobjs =  simu_fmain.o discret_ipopt.o initcond.o const_prof.o

checkhereobjs = DYNOPC_check_derivs.o \

dynopcutilsobjfiles = \
	calc_nx.o \
	cprintf.o \
	daspk_call.o \
	eval_f.o \
	eval_g.o \
	eval_hessobj_v.o \
	get_filename.o \
	get_time.o \
	read_cmp.o \
	write_stp.o

auxilobjs1 = $(addprefix  $(OBJROOTDIR)/$(DIST_DYNOPC_UTILS)/, \
	$(dynopcutilsobjfiles))

dynoptauxilobjfiles = \
	constr.o \
        daemodel_df.o \
        daemodel_df_struc.o \
        daemodel_f.o \
        daemodel_hv.o \
	get_iele4t.o \
	myalloc.o

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

ipoptauxilobjfiles = \
	c_out.o \
	iter_out.o \
	eval_h_dummy.o \
	eval_hesslag_v_dummy.o
auxilobjs5 = $(addprefix  $(OBJROOTDIR)/$(DIST_IPOPT_AUXIL)/, \
	$(ipoptauxilobjfiles))

AUXILOBJS = $(auxilobjs1) $(auxilobjs2) $(auxilobjs3) $(auxilobjs4) \
	$(auxilobjs5) $(mainobj)

CLEANOBJS = $(hereobjs) $(modelc) $(modelcxx) $(modelobj) \
	$(wildcard *.o) $(wildcard *.c) $(wildcard *.cc)
CLEANEXES = $(exe)
CLEANDATA = $(data)

ifndef ESOVERSION   # This is for ADOL-C models (default)
  objs    = $(hereobjs) $(auxilobjs1) $(auxilobjs2) $(auxilobjs3) \
	$(auxilobjs5) $(NLEQ1SOBJ)
  LIBS    = $(DASPKLIB) $(IPOPTLIB) $(TRONLIB) $(IPOPTLIB) $(DAE2NLPLIB) \
	$(ADOLCLIB) $(LINALGLIBS) $(DLLIB)
  checkobjs = $(checkhereobjs) $(auxilobjs1) $(auxilobjs2) $(auxilobjs3) \
	$(auxilobjs5) $(NLEQ1SOBJ)
#  ifdef SOLIBDIR
    objs += $(modelinitaddobj)
    checkobjs += $(modelinitaddobj)
#  else
#    objs += $(modelinitobj) $(modelobj)
#    checkobjs += $(modelinitobj) $(modelobj)
#  endif
else
  dynoptauxilobjfiles += model_hv_dummy.o
  objs    = $(hereobjs) $(auxilobjs1) $(auxilobjs2) \
	$(auxilobjs5) $(NLEQ1SOBJ)
  LIBS    = $(DASPKLIB) $(IPOPTLIB) $(TRONLIB) $(IPOPTLIB) $(DAE2NLPLIB) \
	$(ESO2DAELIB) $(IDLGPROMSLIB) $(LINALGLIBS) $(CORBALIBS)
  checkobjs = $(checkhereobjs) $(auxilobjs1) $(auxilobjs2) \
	$(auxilobjs5) $(NLEQ1SOBJ)
endif

#first:
#	echo $(DISTBUILD) $(filter %.so, $(patsubst %.a, %.so,$(DISTBUILD))) $(AUXILOBJS)

all: $(exe)

$(exe): $(mainobj) $(objs) $(DEPLIBS)
	$(CXX)  $(CXXLNFLAGS) -o $@ $(mainobj) $(objs) $(LNKLIBS) $(CXXLNLIBS)

check: $(mainobj) $(checkobjs) $(DEPLIBS)
	$(CXX)  $(CXXLNFLAGS) -o DYNOPC_check_derivs $(mainobj) $(checkobjs) \
	$(LNKLIBS) $(CXXLNLIBS)

$(data):
	cp $(datadir)/$@ .

model: $(modelso)
#
#  Special Building rules using f2c to obtain C++ model from Fortran
#
.INTERMEDIATE: $(modelc) $(modelcxx)
$(modelc): $(modelf) $(modelinc)
	$(F2C) $(F2CFLAGS) $(FINCFLAGS) $<

$(modelcxx): $(modelc) $(F2CPOST)
	$(F2CPOST) $< $@

$(modelobj): $(modelcxx) $(f2cadolc_h)
	$(CXX) -c $(CXXFLAGS) $(CXXINCFLAGS) $< -o $@

ifeq ($(UNAME),AIX)
$(modelso): $(modelinitobj) $(modelobj)
	$(SOAR) $(SOARFLAGS) $@ $(modelinitobj) $(modelobj) $(ADOLCLIB) #-lc -lC
else
$(modelso): $(modelinitobj) $(modelobj)
	$(SOAR) $(SOARFLAGS) $@ $(modelinitobj) $(modelobj)
endif

ifdef SOLIBDIR
.INTERMEDIATE: $(modelobj) #$(modelcxx)
endif
#
#  Mainly for development:
#  Force check of libraries and object files in other directories
#
ifdef SOLIBDIR
$(patsubst %.a,%.so, $(DAE2NLPLIB) $(IPOPTLIB) $(ESO2DAELIB) $(IDLGPROMSLIB)) \
$(AUXILOBJS) $(F2CPOST): FORCE
else
$(DAE2NLPLIB) $(IPOPTLIB) $(ESO2DAELIB) $(IDLGPROMSLIB) \
$(AUXILOBJS) $(F2CPOST): FORCE
endif

#
#  Dependencies
#
incdirs    = $(SRCROOTDIR)/$(DIST_DYNOPC_UTILS) \
	$(SRCROOTDIR)/$(DIST_DYNOPT_AUXIL)
FINCFLAGS  = $(patsubst %, -I%, $(incdirs))

DYNOPC_check_derivs.o const_prof.o discret_ipopt.o initcond.o simu_fmain.o: \
	$(SRCROOTDIR)/$(DIST_DYNOPC_UTILS)/DYNOPC.INC

DYNOPC_check_derivs.o const_prof.o daspk_call.o discret_ipopt.o initcond.o \
simu_fmain.o: \
	$(SRCROOTDIR)/$(DIST_DYNOPT_AUXIL)/DYNAUX.INC

#
f2cadol_hdir = $(SRCROOTDIR)/$(DIST_DYNOPT_ADOLC2DAE)
cxxincdirs   = $(f2cadol_hdir) $(ADOLCSRCDIR)
CXXINCFLAGS  = $(patsubst %, -I%, $(cxxincdirs))

f2cadolc_h = $(f2cadol_hdir)/f2c_adolc.h
$(modelcc): $(f2cadolc_h)

#
#  Need to establish some order for parallelism in case FLAGS changed
#
break := $(firstword $(auxilobjs))
$(break): $(filter-out $(break), $(auxilobjs))

$(mainobj): $(hereobjs) # Need to point these in order to avoid 

include $(MAKEFILETRGS)
