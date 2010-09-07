# The location of the mesh library
LIBMESH_DIR ?= $(HOME)/Develop/libmesh/libmesh-trunk

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include $(LIBMESH_DIR)/Make.common

# Add flags for compilation with Ipopt
#IPOPTDIR ?= $(HOME)/COIN-svn/OBJgcc-debug/Ipopt/branches/parallel-inexact
IPOPTDIR ?= $(HOME)/COIN-svn/OBJgcc/Ipopt/branches/parallel-inexact
libmesh_INCLUDE += -I$(IPOPTDIR)/include/coin

libmesh_LIBS += -Wl,-rpath,$(IPOPTDIR)/lib/coin -Wl,-rpath,$(IPOPTDIR)/lib/coin/ThirdParty -L$(IPOPTDIR)/lib/coin -L$(IPOPTDIR)/lib/coin/ThirdParty -lipopt `cat $(IPOPTDIR)/share/coin/doc/Ipopt/ipopt_addlibs_cpp.txt`


#libmesh_CXXFLAGS += -DCHECK_SMARTPTR
# libmesh_CXXFLAGS = -DDEBUG  -g -ansi -pedantic -DPETSC_USE_DYNAMIC_LIBRARIES -Wno-long-long  -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
# -felide-constructors -llapack_netlib -lblas_netlib  

###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= Airflow.C IpLibMeshPDE.C IpLibMeshPDENLP.C ProblemGeometry.C

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################

.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	   := ./Airflow


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects)
	@echo "Linking "$@"..."
	$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)


# Useful rules.
clean:
	rm -f $(objects) *~
	rm -f $(target)
clobber:
	@$(MAKE) clean
	@rm -f $(target) out_*

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o

# include the dependency list
include .depend


#
# Dependencies
#
.depend:
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend

###############################################################################
