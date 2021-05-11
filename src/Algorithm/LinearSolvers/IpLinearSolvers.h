// Copyright (C) 2021 COIN-OR Foundation
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef __IPLINEARSOLVERS_H__
#define __IPLINEARSOLVERS_H__

#include "IpoptConfig.h"
#include "IpTypes.h"

#define IPOPTLINEARSOLVER_MA27    0x001u
#define IPOPTLINEARSOLVER_MA57    0x002u
#define IPOPTLINEARSOLVER_MA77    0x004u
#define IPOPTLINEARSOLVER_MA86    0x008u
#define IPOPTLINEARSOLVER_MA97    0x010u
#define IPOPTLINEARSOLVER_MC19    0x020u
#define IPOPTLINEARSOLVER_ALLHSL (IPOPTLINEARSOLVER_MA27 | IPOPTLINEARSOLVER_MA57 | IPOPTLINEARSOLVER_MA77 | IPOPTLINEARSOLVER_MA86 | IPOPTLINEARSOLVER_MA97 | IPOPTLINEARSOLVER_MC19)

#define IPOPTLINEARSOLVER_PARDISO    0x040u
#define IPOPTLINEARSOLVER_PARDISOMKL 0x080u
#define IPOPTLINEARSOLVER_SPRAL   0x100u
#define IPOPTLINEARSOLVER_WSMP    0x200u
#define IPOPTLINEARSOLVER_MUMPS   0x400u

#ifdef __cplusplus
extern "C"
{
#endif

typedef unsigned int IpoptLinearSolver;

/** get bitflags indicating which linear solvers (and mc19) are available to Ipopt
 *
 * If buildinonly if set to a nonzero value, then only report linear solvers that have been linked into Ipopt.
 * Otherwise, also linear solvers are reported which are loaded from a shared library at runtime, if this feature has been compiled.
 * @since 3.14.0
 */
IPOPTLIB_EXPORT IpoptLinearSolver IPOPT_CALLCONV IpoptGetAvailableLinearSolvers(
   int buildinonly
);

#ifdef __cplusplus
}
#endif

#endif /* __IPLINEARSOLVERS_H__ */
