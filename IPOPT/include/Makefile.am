# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 632 2004-07-20 23:52:15Z andreasw $

AUTOMAKE_OPTIONS = foreign

include_HEADERS = Ipopt_Interface.h

all-local: config_f.h

config_f.h: config.h
	$(EGREP) '^#' config.h > config_f.h

CLEANFILES = config_f.h
