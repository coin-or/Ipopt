# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 632 2004-07-20 23:52:15Z andreasw $

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libipampl.a

libipampl_a_SOURCES = get_amplmu.f ipopt_getdata.f

bin_PROGRAMS = ipopt

ipopt_SOURCES = ipoptAMPL.c

LDADD = libipampl.a $(ASLLIB) $(AM_LDADD)

FFLAGS += -I$(srcdir)/../ipopt

CPPFLAGS += $(AMPL_CPPFLAGS) -I$(srcdir)/../ipopt

ipopt_DEPENDENCIES = ../ipopt/libipcore.a
