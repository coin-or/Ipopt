// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2005-10-127

#include "LuksanVlcek5.hpp"

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "configall_system.h"
#endif

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

using namespace Ipopt;

inline static Index Sgn(Number a)
{
  if (a>0.) {
    return 1;
  }
  else {
    return -1;
  }
}

LuksanVlcek5::LuksanVlcek5(Number g_l, Number g_u)
{
  g_l_ = g_l;
  g_u_ = g_u;
}

bool LuksanVlcek5::InitializeProblem(Index N)
{
  N_=N;
  if (N_<=4) {
    printf("N needs to be at least 5.\n");
    return false;
  }
  return true;
}

// returns the size of the problem
bool LuksanVlcek5::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in LuksanVlcek5.hpp has 4 variables, x[0] through x[3]
  n = N_+2;

  m = N_-4;

  nnz_jac_g = 5 * m;

  nnz_h_lag = n + n-1 + n-2;

  // use the C style numbering of matrix indices (starting at 0)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool LuksanVlcek5::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                   Index m, Number* g_l, Number* g_u)
{
  // none of the variables have bounds
  for (Index i=1; i<n-1; i++) {
    x_l[i] = -1e20;
    x_u[i] =  1e20;
  }
  // except for the first and last
  x_l[0] = x_u[0] = 0.;
  x_l[n-1] = x_u[n-1] = 0.;

  // Set the bounds for the constraints
  for (Index i=0; i<m; i++) {
    g_l[i] = g_l_;
    g_u[i] = g_u_;
  }

  return true;
}

// returns the initial point for the problem
bool LuksanVlcek5::get_starting_point(Index n, bool init_x, Number* x,
                                      bool init_z, Number* z_L, Number* z_U,
                                      Index m, bool init_lambda,
                                      Number* lambda)
{
  if (!init_x || init_z || init_lambda) {
    return false;
  }

  // set the starting point
  for (Index i=0; i<n; i++) {
    x[i] = -1.;
  }

  /*
  // DELETEME
  for (Index i=0; i<n; i++) {
    x[i] += 0.1*((Number) i);
  }
  */

  return true;
}

// returns the value of the objective function
bool LuksanVlcek5::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  Number p = 7./3.;
  obj_value = 0.;
  for (Index i=1; i<=N_; i++) {
    Number b = (3.-2.*x[i])*x[i] - x[i-1] - x[i+1] + 1.;
    obj_value += pow(fabs(b),p);
  }

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool LuksanVlcek5::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  Number p = 7./3.;
  grad_f[0] = 0.;
  grad_f[1] = 0.;
  for (Index i=1; i<=N_; i++) {
    Number b = (3.-2.*x[i])*x[i] - x[i-1] - x[i+1] + 1.;
    Number pb = pow(fabs(b),p-1.);

    grad_f[i+1] = -p*Sgn(b)*pb;
    grad_f[i-1] += grad_f[i+1];
    grad_f[i] += p*Sgn(b)*(3.-4.*x[i])*pb;
  }

  return true;
}

// return the value of the constraints: g(x)
bool LuksanVlcek5::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  for (Index i=0; i<N_-4; i++) {
    g[i] = 8.*x[i+3]*(x[i+3]*x[i+3]-x[i+2])
           - 2.*(1-x[i+3])
           + 4.*(x[i+3]-x[i+4]*x[i+4])
           + x[i+2]*x[i+2]
           - x[i+1]
           + x[i+4]
           - x[i+5]*x[i+5];
  }
  return true;
}

// return the structure or values of the jacobian
bool LuksanVlcek5::eval_jac_g(Index n, const Number* x, bool new_x,
                              Index m, Index nele_jac, Index* iRow, Index *jCol,
                              Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian

    Index ijac = 0;
    for (Index i=0; i<N_-4; i++) {
      iRow[ijac] = i;
      jCol[ijac] = i+1;
      ijac++;
      iRow[ijac] = i;
      jCol[ijac] = i+2;
      ijac++;
      iRow[ijac] = i;
      jCol[ijac] = i+3;
      ijac++;
      iRow[ijac] = i;
      jCol[ijac] = i+4;
      ijac++;
      iRow[ijac] = i;
      jCol[ijac] = i+5;
      ijac++;
    }
    DBG_ASSERT(ijac == nele_jac);
  }
  else {
    // return the values of the jacobian of the constraints

    Index ijac = 0;
    for (Index i=0; i<N_-4; i++) {
      values[ijac] = -1.;
      ijac++;
      values[ijac] = -8.*x[i+3] + 2.*x[i+2];
      ijac++;
      values[ijac] = 6. - 8.*x[i+2] + 24.*x[i+3]*x[i+3];
      ijac++;
      values[ijac] = -8.*x[i+4] + 1.;
      ijac++;
      values[ijac] = -2.*x[i+5];
      ijac++;
    }
  }

  return true;
}

//return the structure or values of the hessian
bool LuksanVlcek5::eval_h(Index n, const Number* x, bool new_x,
                          Number obj_factor, Index m, const Number* lambda,
                          bool new_lambda, Index nele_hess, Index* iRow,
                          Index* jCol, Number* values)
{
  if (values == NULL) {
    Index ihes=0;
    // First the diagonal
    for (Index i=0; i<n; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i;
      ihes++;
    }
    // Now the first off-diagonal
    for (Index i=0; i<n-1; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i+1;
      ihes++;
    }
    // And finally the second off-diagonal
    for (Index i=0; i<n-2; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i+2;
      ihes++;
    }
    DBG_ASSERT(ihes == nele_hess);
  }
  else {
    Number p = 7./3.;
    // First the objective
    values[0] = 0.;
    values[1] = 0.;
    values[n] = 0.;
    for (Index i=1; i<=N_; i++) {
      Number b = (3.-2.*x[i])*x[i] - x[i-1] - x[i+1] + 1.;
      Number pb1 = pow(fabs(b),p-1.);
      Number pb2 = pow(fabs(b),p-2.);
      Number a1 = 3.-4.*x[i];
      Number a2 = p*(p-1.)*pb2;
      Number a3 = a1*a2;

      // x[i-1] x[i-1]
      values[i-1] += obj_factor*a2;
      // x[i-1] x[i]
      values[n+i-1] += obj_factor*(-a3);
      // x[i-1] x[i+1]
      values[n+(n-1)+i-1] = obj_factor*a2;
      // x[i] x[i]
      values[i] += obj_factor*(a3*a1 - 4.*p*Sgn(b)*pb1);
      // x[i] x[i+1]
      values[n+i] = obj_factor*(-a3);
      // x[i+1] x[i+1]
      values[i+1] = obj_factor*a2;
    }

    // Now the constraints
    for (Index i=0; i<N_-4; i++) {
      // x[i+2] x[i+2]
      values[i+2] += lambda[i]*2.;
      // x[i+3] x{i+3]
      values[i+3] += lambda[i]*48.*x[i+3];
      // x[i+4] x[i+4]
      values[i+4] += lambda[i]*(-8.);
      // x[i+5] x[i+5]
      values[i+5] += lambda[i]*(-2.);
      // x[i+2] x[i+3]
      values[n+i+2] += -lambda[i]*8.;
    }
  }
  return true;
}

void LuksanVlcek5::finalize_solution(SolverReturn status,
                                     Index n, const Number* x, const Number* z_L, const Number* z_U,
                                     Index m, const Number* g, const Number* lambda,
                                     Number obj_value,
				     const IpoptData* ip_data,
				     IpoptCalculatedQuantities* ip_cq)
{}

