// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2005-10-127

#include "LuksanVlcek7.hpp"

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

LuksanVlcek7::LuksanVlcek7(Number g_l, Number g_u)
{
  g_l_ = g_l;
  g_u_ = g_u;
}

bool LuksanVlcek7::InitializeProblem(Index N)
{
  N_=N;
  if (N_<3) {
    printf("N has to be at least 3.\n");
    return false;
  }
  return true;
}

// returns the size of the problem
bool LuksanVlcek7::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in LuksanVlcek7.hpp has 4 variables, x[0] through x[3]
  n = N_+2;

  m = 4;

  nnz_jac_g = 3 + 4 + 4 + 3;

  nnz_h_lag = n + 3;

  // use the C style numbering of matrix indices (starting at 0)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool LuksanVlcek7::get_bounds_info(Index n, Number* x_l, Number* x_u,
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
bool LuksanVlcek7::get_starting_point(Index n, bool init_x, Number* x,
                                      bool init_z, Number* z_L, Number* z_U,
                                      Index m, bool init_lambda,
                                      Number* lambda)
{
  if (!init_x || init_z || init_lambda) {
    return false;
  }

  // set the starting point
  for (Index i=0; i<n; i++) {
    x[i] = 1.;
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
bool LuksanVlcek7::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  obj_value = 0.;
  for (Index i=1; i<=N_; i++) {
    obj_value += i*((1.-cos(x[i])) + sin(x[i-1]) - sin(x[i+1]));
  }

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool LuksanVlcek7::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  grad_f[0] = 0.;
  grad_f[1] = 0.;
  for (Index i=1; i<=N_; i++) {
    grad_f[i-1] += i*cos(x[i-1]);
    grad_f[i] += i*sin(x[i]);
    grad_f[i+1] = -i*cos(x[i+1]);
  }

  return true;
}

// return the value of the constraints: g(x)
bool LuksanVlcek7::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  g[0] = 4.*(x[1]-x[2]*x[2]) + x[2] - x[3]*x[3];
  g[1] = 8.*x[2]*(x[2]*x[2]-x[1])
         - 2.*(1.-x[2]) + 4.*(x[2]-x[3]*x[3]) + x[3] - x[4]*x[4];
  g[2] = 8.*x[N_-1]*(x[N_-1]*x[N_-1]-x[N_-2])
         - 2.*(1.-x[N_-1])
         + 4.*(x[N_-1]-x[N_]*x[N_])
         + x[N_-2]*x[N_-2]
         - x[N_-3];
  g[3] = 8.*x[N_]*(x[N_]*x[N_]-x[N_-1])
         - 2.*(1.-x[N_]) + x[N_-1]*x[N_-1] - x[N_-2];
  return true;
}

// return the structure or values of the jacobian
bool LuksanVlcek7::eval_jac_g(Index n, const Number* x, bool new_x,
                              Index m, Index nele_jac, Index* iRow, Index *jCol,
                              Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian

    Index ijac = 0;
    // g[0]
    iRow[ijac] = 0;
    jCol[ijac] = 1;
    ijac++;
    iRow[ijac] = 0;
    jCol[ijac] = 2;
    ijac++;
    iRow[ijac] = 0;
    jCol[ijac] = 3;
    ijac++;

    // g[1]
    iRow[ijac] = 1;
    jCol[ijac] = 1;
    ijac++;
    iRow[ijac] = 1;
    jCol[ijac] = 2;
    ijac++;
    iRow[ijac] = 1;
    jCol[ijac] = 3;
    ijac++;
    iRow[ijac] = 1;
    jCol[ijac] = 4;
    ijac++;

    // g[2]
    iRow[ijac] = 2;
    jCol[ijac] = N_-3;
    ijac++;
    iRow[ijac] = 2;
    jCol[ijac] = N_-2;
    ijac++;
    iRow[ijac] = 2;
    jCol[ijac] = N_-1;
    ijac++;
    iRow[ijac] = 2;
    jCol[ijac] = N_;
    ijac++;

    // g[3]
    iRow[ijac] = 3;
    jCol[ijac] = N_-2;
    ijac++;
    iRow[ijac] = 3;
    jCol[ijac] = N_-1;
    ijac++;
    iRow[ijac] = 3;
    jCol[ijac] = N_;
    ijac++;

    DBG_ASSERT(ijac == nele_jac);
  }
  else {
    // return the values of the jacobian of the constraints

    Index ijac = 0;
    // g[0]
    values[ijac] = 4.;
    ijac++;
    values[ijac] = 1. - 8.*x[2];
    ijac++;
    values[ijac] = -2*x[3];
    ijac++;
    // g[1]
    values[ijac] = -8.*x[2];
    ijac++;
    values[ijac] = 24.*x[2]*x[2] - 8.*x[1] + 6.;
    ijac++;
    values[ijac] = -8*x[3] + 1.;
    ijac++;
    values[ijac] = -2.*x[4];
    ijac++;
    // g[2]
    values[ijac] = -1.;
    ijac++;
    values[ijac] = -8.*x[N_-1] + 2.*x[N_-2];
    ijac++;
    values[ijac] = 24.*x[N_-1]*x[N_-1] - 8.*x[N_-2] + 6.;
    ijac++;
    values[ijac] = -8.*x[N_];
    ijac++;
    // g[3]
    values[ijac] = -1.;
    ijac++;
    values[ijac] = -8.*x[N_] + 2.*x[N_-1];
    ijac++;
    values[ijac] = 24.*x[N_]*x[N_] - 8.*x[N_-1] + 2.;
    ijac++;
  }

  return true;
}

//return the structure or values of the hessian
bool LuksanVlcek7::eval_h(Index n, const Number* x, bool new_x,
                          Number obj_factor, Index m, const Number* lambda,
                          bool new_lambda, Index nele_hess, Index* iRow,
                          Index* jCol, Number* values)
{
  if (values == NULL) {
    // The diagonal
    for (Index i=0; i<n; i++) {
      iRow[i] = i;
      jCol[i] = i;
    }
    // x[1] x[2]
    iRow[n] = 1;
    jCol[n] = 2;
    // x[N_-2] x[N_-1]
    iRow[n+1] = N_-2;
    jCol[n+1] = N_-1;
    // x[N_-1] x[N_]
    iRow[n+2] = N_-1;
    jCol[n+2] = N_;
  }
  else {

    // Objective function
    values[0] = 0.;
    values[1] = 0.;
    for (Index i=1; i<=N_; i++) {
      values[i-1] -= obj_factor*(i*sin(x[i-1]));
      values[i] += obj_factor*i*cos(x[i]);
      values[i+1] = obj_factor*i*sin(x[i+1]);
    }

    // g[0]
    values[2] -= lambda[0]*8.;
    values[3] -= lambda[0]*2.;
    // g[1]
    values[2] += lambda[1]*48.*x[2];
    values[3] -= lambda[1]*8.;
    values[4] -= lambda[1]*2.;
    values[n] = -lambda[1]*8.;
    // g[2]
    values[N_-2] += lambda[2]*2.;
    values[N_-1] += lambda[2]*48.*x[N_-1];
    values[N_] -= lambda[2]*8.;
    values[n+1] = -lambda[2]*8.;
    // g[3]
    values[N_-1] += lambda[3]*2.;
    values[N_] += lambda[3]*48.*x[N_];
    values[n+2] = -lambda[3]*8.;

  }
  return true;
}

void LuksanVlcek7::finalize_solution(SolverReturn status,
                                     Index n, const Number* x, const Number* z_L, const Number* z_U,
                                     Index m, const Number* g, const Number* lambda,
                                     Number obj_value,
				     const IpoptData* ip_data,
				     IpoptCalculatedQuantities* ip_cq)
{}

