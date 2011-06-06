// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2005-10-127

#include "LuksanVlcek4.hpp"

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

LuksanVlcek4::LuksanVlcek4(Number g_l, Number g_u)
{
  g_l_ = g_l;
  g_u_ = g_u;
}

bool LuksanVlcek4::InitializeProblem(Index N)
{
  N_=N;
  if (N_<=1 || 4*((N_+2)/4)!=N_+2) {
    printf("N needs to be at least 2 and N+2 needs to be a multiple of 4.\n");
    return false;
  }
  return true;
}

// returns the size of the problem
bool LuksanVlcek4::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in LuksanVlcek4.hpp has 4 variables, x[0] through x[3]
  n = N_+2;

  m = N_-2;

  nnz_jac_g = 3 * m;

  nnz_h_lag = n + n-1;

  // use the C style numbering of matrix indices (starting at 0)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool LuksanVlcek4::get_bounds_info(Index n, Number* x_l, Number* x_u,
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
bool LuksanVlcek4::get_starting_point(Index n, bool init_x, Number* x,
                                      bool init_z, Number* z_L, Number* z_U,
                                      Index m, bool init_lambda,
                                      Number* lambda)
{
  if (!init_x || init_z || init_lambda) {
    return false;
  }

  // set the starting point
  for (Index i=0; i<n/4; i++) {
    x[4*i] = 1.;
    x[4*i+1] = 2.;
    x[4*i+2] = 2.;
    x[4*i+3] = 2.;
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
bool LuksanVlcek4::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  obj_value = 0.;
  for (Index i=0; i<N_/2; i++) {
    Number e0 = exp(x[2*i]);
    Number e0mx1 = e0 - x[2*i+1];
    Number x1mx2 = x[2*i+1] - x[2*i+2];
    Number x2mx3 = x[2*i+2] - x[2*i+3];
    Number t = tan(x2mx3);
    Number x3m1 = x[2*i+3] - 1.;
    obj_value += pow(e0mx1,4) + 100.*pow(x1mx2,6) + pow(t,4)
                 + pow(x[2*i],8) + x3m1*x3m1;
  }

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool LuksanVlcek4::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  grad_f[0] = 0.;
  grad_f[1] = 0.;
  for (Index i=0; i<N_/2; i++) {
    Number e0 = exp(x[2*i]);
    Number e0mx1 = e0 - x[2*i+1];
    Number x1mx2 = x[2*i+1] - x[2*i+2];
    Number x2mx3 = x[2*i+2] - x[2*i+3];
    Number x3m1 = x[2*i+3] - 1.;
    Number dt = 4.*pow(tan(x2mx3),3)/pow(cos(x2mx3),2);

    grad_f[2*i] += 4.*e0*pow(e0mx1,3) + 8.*pow(x[2*i],7);
    grad_f[2*i+1] += -4.*pow(e0mx1,3) + 600.*pow(x1mx2,5);
    grad_f[2*i+2] = -600.*pow(x1mx2,5) + dt;
    grad_f[2*i+3] = -dt + 2.*x3m1;
  }

  return true;
}

// return the value of the constraints: g(x)
bool LuksanVlcek4::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  for (Index i=0; i<N_-2; i++) {
    g[i] = 8.*x[i+1]*(x[i+1]*x[i+1]-x[i]) - 2.*(1-x[i+1])
           + 4.*(x[i+1]-x[i+2]*x[i+2]);
  }
  return true;
}

// return the structure or values of the jacobian
bool LuksanVlcek4::eval_jac_g(Index n, const Number* x, bool new_x,
                              Index m, Index nele_jac, Index* iRow, Index *jCol,
                              Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian

    Index ijac = 0;
    for (Index i=0; i<N_-2; i++) {
      iRow[ijac] = i;
      jCol[ijac] = i;
      ijac++;
      iRow[ijac] = i;
      jCol[ijac] = i+1;
      ijac++;
      iRow[ijac] = i;
      jCol[ijac] = i+2;
      ijac++;
    }
    DBG_ASSERT(ijac == nele_jac);
  }
  else {
    // return the values of the jacobian of the constraints

    Index ijac = 0;
    for (Index i=0; i<N_-2; i++) {
      values[ijac] = -8.*x[i+1];
      ijac++;
      values[ijac] = 6. - 8.*x[i] + 24.*x[i+1]*x[i+1];
      ijac++;
      values[ijac] = -8.*x[i+2];
      ijac++;
    }
  }

  return true;
}

//return the structure or values of the hessian
bool LuksanVlcek4::eval_h(Index n, const Number* x, bool new_x,
                          Number obj_factor, Index m, const Number* lambda,
                          bool new_lambda, Index nele_hess, Index* iRow,
                          Index* jCol, Number* values)
{
  if (values == NULL) {
    Index ihes=0;
    for (Index i=0; i<n-1; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i;
      ihes++;
      iRow[ihes] = i;
      jCol[ihes] = i+1;
      ihes++;
    }
    iRow[ihes] = n-1;
    jCol[ihes] = n-1;
    ihes++;
    DBG_ASSERT(ihes == nele_hess);
  }
  else {
    // First the objective
    Index ihes=0;
    values[0] = 0.;
    values[1] = 0.;
    values[2] = 0.;
    for (Index i=0; i<N_/2; i++) {
      Number e0 = exp(x[2*i]);
      Number e0mx1 = e0 - x[2*i+1];
      Number x1mx2 = x[2*i+1] - x[2*i+2];
      Number x2mx3 = x[2*i+2] - x[2*i+3];
      Number s = sin(x2mx3);
      Number ss = s*s;
      Number c = cos(x2mx3);
      Number ddt = 4.*(3.*ss*c*c + 5.*ss*ss)/pow(c,6);

      // x[2*i] x[2*i]
      values[ihes] += obj_factor*(4.*e0*pow(e0mx1,3)
                                  + 12*e0*e0*e0mx1*e0mx1
                                  + 56.*pow(x[2*i],6));
      ihes++;
      // x[2*i] x[2*i+1]
      values[ihes] += obj_factor*(-12*e0*e0mx1*e0mx1);
      ihes++;
      // x[2*i+1] x[2*i+1]
      values[ihes] += obj_factor*(3000.*pow(x1mx2,4) + 12.*e0mx1*e0mx1);
      ihes++;
      // x[2*i+1] x[2*i+2]
      values[ihes] = -obj_factor*(3000.*pow(x1mx2,4));
      ihes++;
      // x[2*i+2] x[2*i+2]
      values[ihes] = obj_factor*(3000.*pow(x1mx2,4) + ddt);
      // x[2*i+2] x[2*i+3]
      values[ihes+1] = -obj_factor*ddt;
      // x[2*i+3] x[2*i+3]
      values[ihes+2] = obj_factor*(2. + ddt);
    }

    // Now the constraints
    ihes = 0;
    for (Index i=0; i<N_-2; i++) {
      // x[i] x[i+1]
      values[2*i+1] -= lambda[i]*8.;
      // x[i+1] x[i+1]
      values[2*i+2] += lambda[i]*48.*x[i+1];
      // x[i+2] x[i+2]
      values[2*i+4] -= lambda[i]*8.;
    }
  }
  return true;
}

void LuksanVlcek4::finalize_solution(SolverReturn status,
                                     Index n, const Number* x, const Number* z_L, const Number* z_U,
                                     Index m, const Number* g, const Number* lambda,
                                     Number obj_value,
				     const IpoptData* ip_data,
				     IpoptCalculatedQuantities* ip_cq)
{}

