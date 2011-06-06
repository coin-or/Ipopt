// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2005-10-127

#include "LuksanVlcek3.hpp"

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

LuksanVlcek3::LuksanVlcek3(Number g_l, Number g_u)
{
  g_l_ = g_l;
  g_u_ = g_u;
}

bool LuksanVlcek3::InitializeProblem(Index N)
{
  N_=N;
  if (N_<=5 || 4*((N_+2)/4)!=N_+2) {
    printf("N needs to be at least 6 and N+2 needs to be a multiple of 4.\n");
    return false;
  }
  return true;
}

// returns the size of the problem
bool LuksanVlcek3::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in LuksanVlcek3.hpp has 4 variables, x[0] through x[3]
  n = N_+2;

  m = 2;

  nnz_jac_g = 4;

  nnz_h_lag = 5*N_/2 + 3;

  // use the C style numbering of matrix indices (starting at 0)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool LuksanVlcek3::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                   Index m, Number* g_l, Number* g_u)
{
  // none of the variables have bounds
  for (Index i=0; i<n; i++) {
    x_l[i] = -1e20;
    x_u[i] =  1e20;
  }

  // Set the bounds for the constraints
  for (Index i=0; i<m; i++) {
    g_l[i] = g_l_;
    g_u[i] = g_u_;
  }

  return true;
}

// returns the initial point for the problem
bool LuksanVlcek3::get_starting_point(Index n, bool init_x, Number* x,
                                      bool init_z, Number* z_L, Number* z_U,
                                      Index m, bool init_lambda,
                                      Number* lambda)
{
  if (!init_x || init_z || init_lambda) {
    return false;
  }

  // set the starting point
  for (Index i=0; i<n/4; i++) {
    x[4*i] = 3.;
    x[4*i+1] = -1.;
    x[4*i+2] = 0.;
    x[4*i+3] = 1.;
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
bool LuksanVlcek3::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  obj_value = 0.;
  for (Index i=0; i<N_/2; i++) {
    Number a1 = x[2*i]+10.*x[2*i+1];
    Number a2 = x[2*i+2] - x[2*i+3];
    Number a3 = x[2*i+1] - 2.*x[2*i+2];
    Number a4 = x[2*i] - x[2*i+3];
    obj_value += a1*a1 + 5.*a2*a2 + pow(a3,4)+ 10.*pow(a4,4);
  }

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool LuksanVlcek3::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  grad_f[0] = 0.;
  grad_f[1] = 0.;
  for (Index i=0; i<N_/2; i++) {
    Number a1 = x[2*i]+10.*x[2*i+1];
    Number a2 = x[2*i+2] - x[2*i+3];
    Number a3 = x[2*i+1] - 2.*x[2*i+2];
    Number a4 = x[2*i] - x[2*i+3];
    grad_f[2*i] += 2.*a1 + 40.*pow(a4,3);
    grad_f[2*i+1] += 20.*a1 + 4.*pow(a3,3);
    grad_f[2*i+2] = 10.*a2 - 8.*pow(a3,3);
    grad_f[2*i+3] = -10.*a2 - 40.*pow(a4,3);
  }

  return true;
}

// return the value of the constraints: g(x)
bool LuksanVlcek3::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  g[0] = 3.*pow(x[0],3) + 2.*x[1] - 5. + sin(x[0]-x[1])*sin(x[0]+x[1]);
  g[1] = 4.*x[n-3] - x[n-4]*exp(x[n-4]-x[n-3]) - 3;
  ;
  return true;
}

// return the structure or values of the jacobian
bool LuksanVlcek3::eval_jac_g(Index n, const Number* x, bool new_x,
                              Index m, Index nele_jac, Index* iRow, Index *jCol,
                              Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian

    Index ijac = 0;
    iRow[ijac] = 0;
    jCol[ijac] = 0;
    ijac++;
    iRow[ijac] = 0;
    jCol[ijac] = 1;
    ijac++;
    iRow[ijac] = 1;
    jCol[ijac] = n-4;
    ijac++;
    iRow[ijac] = 1;
    jCol[ijac] = n-3;
    ijac++;

    DBG_ASSERT(ijac == nele_jac);
  }
  else {
    // return the values of the jacobian of the constraints

    Index ijac = 0;
    values[ijac] = 9.*x[0]*x[0]
                   + cos(x[0]-x[1])*sin(x[0]+x[1])
                   + sin(x[0]-x[1])*cos(x[0]+x[1]);
    ijac++;
    values[ijac] = 2.
                   - cos(x[0]-x[1])*sin(x[0]+x[1])
                   + sin(x[0]-x[1])*cos(x[0]+x[1]);
    ijac++;
    values[ijac] = -(1.+x[n-4])*exp(x[n-4]-x[n-3]);
    ijac++;
    values[ijac] = 4. + x[n-4]*exp(x[n-4]-x[n-3]);
    ijac++;
  }

  return true;
}

//return the structure or values of the hessian
bool LuksanVlcek3::eval_h(Index n, const Number* x, bool new_x,
                          Number obj_factor, Index m, const Number* lambda,
                          bool new_lambda, Index nele_hess, Index* iRow,
                          Index* jCol, Number* values)
{
  if (values == NULL) {
    Index ihes=0;
    for (Index i=0; i<N_/2; i++) {
      iRow[ihes] = 2*i;
      jCol[ihes] = 2*i;
      ihes++;
      iRow[ihes] = 2*i;
      jCol[ihes] = 2*i+1;
      ihes++;
      iRow[ihes] = 2*i;
      jCol[ihes] = 2*i+3;
      ihes++;
      iRow[ihes] = 2*i+1;
      jCol[ihes] = 2*i+1;
      ihes++;
      iRow[ihes] = 2*i+1;
      jCol[ihes] = 2*i+2;
      ihes++;
    }
    iRow[ihes] = N_;
    jCol[ihes] = N_;
    ihes++;
    iRow[ihes] = N_;
    jCol[ihes] = N_+1;
    ihes++;
    iRow[ihes] = N_+1;
    jCol[ihes] = N_+1;
    ihes++;
    DBG_ASSERT(ihes == nele_hess);
  }
  else {
    Index ihes=0;
    values[0] = 0.;
    values[1] = 0.;
    values[3] = 0.;
    for (Index i=0; i<N_/2-1; i++) {
      Number a3 = x[2*i+1] - 2.*x[2*i+2];
      Number a4 = x[2*i] - x[2*i+3];
      // x[2*i] x[2*i]
      values[ihes] += obj_factor*(2. + 120.*a4*a4);
      ihes++;
      // x[2*i] x[2*i+1]
      values[ihes] += obj_factor*20.;
      ihes++;
      // x[2*i] x[2*i+3]
      values[ihes] = -obj_factor*120.*a4*a4;
      ihes++;
      // x[2*i+1] x[2*i+1]
      values[ihes] += obj_factor*(200. + 12.*a3*a3);
      ihes++;
      // x[2*i+1] x[2*i+2]
      values[ihes] = -obj_factor*24.*a3*a3;
      ihes++;
      // x[2*i+2] x[2*i+2]
      values[ihes] = obj_factor*(10. + 48.*a3*a3);
      // x[2*i+2] x[2*i+3]
      values[ihes+1] = -obj_factor*10.;
      // x[2*i+3] x[2*i+3]
      values[ihes+3] = obj_factor*(10. + 120.*a4*a4);
    }
    {
      Index i = N_/2-1;
      Number a3 = x[2*i+1] - 2.*x[2*i+2];
      Number a4 = x[2*i] - x[2*i+3];
      // x[2*i] x[2*i]
      values[ihes] += obj_factor*(2. + 120.*a4*a4);
      ihes++;
      // x[2*i] x[2*i+1]
      values[ihes] += obj_factor*20.;
      ihes++;
      // x[2*i] x[2*i+3]
      values[ihes] = -obj_factor*120.*a4*a4;
      ihes++;
      // x[2*i+1] x[2*i+1]
      values[ihes] += obj_factor*(200. + 12.*a3*a3);
      ihes++;
      // x[2*i+1] x[2*i+2]
      values[ihes] = -obj_factor*24.*a3*a3;
      ihes++;
      // x[2*i+2] x[2*i+2]
      values[ihes] = obj_factor*(10. + 48.*a3*a3);
      // x[2*i+2] x[2*i+3]
      values[ihes+1] = -obj_factor*10.;
      // x[2*i+3] x[2*i+3]
      values[ihes+2] = obj_factor*(10. + 120.*a4*a4);
    }

    // Now the constraints
    ihes = 0;
    Number d1 = x[0] - x[1];
    Number d2 = x[0] + x[1];
    values[ihes] += lambda[0]*(18.*x[0]
                               -2.*sin(d1)*sin(d2)
                               +2.*cos(d1)*cos(d2));
    ihes+=3;
    values[ihes] += lambda[0]*(-2.*sin(d1)*sin(d2)
                               -2.*cos(d1)*cos(d2));

    d1 = x[n-4]-x[n-3];

    // x[n-4] x[n-4]
    ihes = nele_hess - 8;
    values[ihes] += -lambda[1]*(2.+x[n-4])*exp(d1);

    // x[n-4] x[n-3]
    ihes++;
    values[ihes] += lambda[1]*(1.+x[n-4])*exp(d1);

    // x[n-3] x[n-3]
    ihes += 2;
    values[ihes] += -lambda[1]*x[n-4]*exp(d1);
  }
  return true;
}

void LuksanVlcek3::finalize_solution(SolverReturn status,
                                     Index n, const Number* x, const Number* z_L, const Number* z_U,
                                     Index m, const Number* g, const Number* lambda,
                                     Number obj_value,
				     const IpoptData* ip_data,
				     IpoptCalculatedQuantities* ip_cq)
{}

