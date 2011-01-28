// Copyright (C) 2009 International Business Machines
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter, Sanjeeb Dash        IBM    2009-07-10

#ifndef __IPMPI_HPP__
#define __IPMPI_HPP__

#include "IpoptConfig.h"

#ifndef HAVE_MPI
# error This header should only be included if Ipopt was configured with MPI
#endif

// This is a wrapper for MPI header file inclusions.  Those migh be
// different on different systems.

//extern "C" {
#define MPICH_SKIP_MPICXX
#include "mpi.h"
//}

#endif
