// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2005-10-127

#include "LuksanVlcek2.hpp"

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

LuksanVlcek2::LuksanVlcek2(Number g_l, Number g_u)
{
  g_l_ = g_l;
  g_u_ = g_u;
}

bool LuksanVlcek2::InitializeProblem(Index N)
{
  N_=N;
  if (N_<=13 || 2*(N_/2)!=N_) {
    printf("N needs to be at least 14 and even.\n");
    return false;
  }
  return true;
}

// returns the size of the problem
bool LuksanVlcek2::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in LuksanVlcek2.hpp has 4 variables, x[0]
  // through x[3]
  n = N_+2;

  m = N_-7;

  nnz_jac_g = 25 + (m-5)*8;

  nnz_h_lag = n + N_ + 1;

  // use the C style numbering of matrix indices (starting at 0)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool LuksanVlcek2::get_bounds_info(Index n, Number* x_l, Number* x_u,
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
bool LuksanVlcek2::get_starting_point(Index n, bool init_x, Number* x,
                                      bool init_z, Number* z_L, Number* z_U,
                                      Index m, bool init_lambda,
                                      Number* lambda)
{
  if (!init_x || init_z || init_lambda) {
    return false;
  }

  // set the starting point
  for (Index i=0; i<n/2; i++) {
    x[2*i] = -2.;
    x[2*i+1] = 1.;
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
bool LuksanVlcek2::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  obj_value = 0.;
  for (Index i=0; i<N_/2; i++) {
    Number a1 = x[2*i]*x[2*i] - x[2*i+1];
    Number a2 = x[2*i] - 1.;
    Number a3 = x[2*i+2]*x[2*i+2] - x[2*i+3];
    Number a4 = x[2*i+2] - 1.;
    Number a5 = x[2*i+1] + x[2*i+3] - 2.;
    Number a6 = x[2*i+1] - x[2*i+3];
    obj_value += 100.*a1*a1 + a2*a2 + 90.*a3*a3 + a4*a4 + 10.*a5*a5 + .1*a6*a6;
  }

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool LuksanVlcek2::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  grad_f[0] = 0.;
  grad_f[1] = 0.;
  for (Index i=0; i<N_/2; i++) {
    grad_f[2*i] += 400.*x[2*i]*(x[2*i]*x[2*i]-x[2*i+1]) + 2.*(x[2*i]-1.);
    grad_f[2*i+1] += -200.*(x[2*i]*x[2*i]-x[2*i+1])
                     + 20*(x[2*i+1] + x[2*i+3] - 2.)
                     + 0.2*(x[2*i+1] - x[2*i+3]);
    grad_f[2*i+2] = 360.*x[2*i+2]*(x[2*i+2]*x[2*i+2] - x[2*i+3])
                    + 2.*(x[2*i+2] -1.);
    grad_f[2*i+3] = -180.*(x[2*i+2]*x[2*i+2] - x[2*i+3])
                    + 20.*(x[2*i+1] + x[2*i+3] -2.)
                    - 0.2*(x[2*i+1] - x[2*i+3]);
  }

  return true;
}

// return the value of the constraints: g(x)
bool LuksanVlcek2::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  for (Index i=0; i<N_-7; i++) {
    g[i] = (2.+5.*x[i+5]*x[i+5])*x[i+5] + 1.;
    for (Index k=Max(0,i-5); k<=i+1; k++) {
      g[i] += x[k]*(x[k]+1.);
    }
  }

  return true;
}

// return the structure or values of the jacobian
bool LuksanVlcek2::eval_jac_g(Index n, const Number* x, bool new_x,
                              Index m, Index nele_jac, Index* iRow, Index *jCol,
                              Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian

    Index ijac=0;
    for (Index i=0; i<N_-7; i++) {
      for (Index k=Max(0,i-5); k<=i+1; k++) {
        iRow[ijac] = i;
        jCol[ijac] = k;
        ijac++;
      }
      iRow[ijac] = i;
      jCol[ijac] = i+5;
      ijac++;
    }
    DBG_ASSERT(ijac == nele_jac);
  }
  else {
    // return the values of the jacobian of the constraints

    Index ijac=0;
    for (Index i=0; i<N_-7; i++) {
      for (Index k=Max(0,i-5); k<=i+1; k++) {
        values[ijac] = 2.*x[k] + 1.;
        ijac++;
      }
      values[ijac] = 2. + 15.*x[i+5]*x[i+5];
      ijac++;
    }
  }

  return true;
}

//return the structure or values of the hessian
bool LuksanVlcek2::eval_h(Index n, const Number* x, bool new_x,
                          Number obj_factor, Index m, const Number* lambda,
                          bool new_lambda, Index nele_hess, Index* iRow,
                          Index* jCol, Number* values)
{
  if (values == NULL) {
    Index ihes=0;
    // First the diagonal elements
    for (Index i=0; i<n; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i;
      ihes++;
    }
    // And now the off-diagonal elements
    for (Index i=0; i<N_/2; i++) {
      iRow[ihes] = 2*i;
      jCol[ihes] = 2*i+1;
      ihes++;
      iRow[ihes] = 2*i+1;
      jCol[ihes] = 2*i+3;
      ihes++;
    }
    iRow[ihes] = n-2;
    jCol[ihes] = n-1;
    ihes++;
    DBG_ASSERT(ihes == nele_hess);
  }
  else {
    // First we take care of the diagonal elements coming from the
    // objective function
    values[0] = 0.;
    values[1] = 0.;
    for (Index i=0; i<N_/2; i++) {
      values[2*i] += obj_factor*(1200.*x[2*i]*x[2*i] - 400.*x[2*i+1] + 2.);
      values[2*i+1] += obj_factor*220.2;
      values[2*i+2] = obj_factor*(1080.*x[2*i+2]*x[2*i+2] - 360*x[2*i+3] + 2.);
      values[2*i+3] = obj_factor*200.2;
    }
    // Now we take care of the off-diagonal elements coming from the
    // objective function
    Index ihes = n;
    values[ihes] = 0.;
    for (Index i=0; i<N_/2; i++) {
      values[ihes] += obj_factor*(-400.*x[2*i]);
      ihes++;
      values[ihes] = obj_factor*19.8;
      ihes++;
      values[ihes] = obj_factor*(-360.*x[2*i+2]);
    }

    // Ok, now the diagonal elements from the constraints
    for (Index i=0; i<N_-7; i++) {
      for (Index k=Max(0,i-5); k<=i+1; k++) {
        values[k] += lambda[i]*2.;
      }
      values[i+5] += lambda[i]*30.*x[i+5];
    }
  }

  return true;
}

void LuksanVlcek2::finalize_solution(SolverReturn status,
                                     Index n, const Number* x, const Number* z_L, const Number* z_U,
                                     Index m, const Number* g, const Number* lambda,
                                     Number obj_value,
				     const IpoptData* ip_data,
				     IpoptCalculatedQuantities* ip_cq)
{}

