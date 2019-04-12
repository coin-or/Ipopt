# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 632 2004-07-20 23:52:15Z andreasw $

AUTOMAKE_OPTIONS = foreign

SUBDIRS =

if BUILDBLASDIR
  SUBDIRS += blas
endif

if BUILDLAPACKDIR
  SUBDIRS += lapack
endif

if BUILDHSLDIR
  SUBDIRS += HSL
endif

if BUILDASLDIR
  SUBDIRS += ASL
endif

if BUILDTRONDIR
  SUBDIRS += TRON
endif
