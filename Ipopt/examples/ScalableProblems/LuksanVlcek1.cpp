// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2005-10-127

#include "LuksanVlcek1.hpp"

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

LuksanVlcek1::LuksanVlcek1(Number g_l, Number g_u)
{
  g_l_ = g_l;
  g_u_ = g_u;
}

bool LuksanVlcek1::InitializeProblem(Index N)
{
  N_=N;
  if (N_<=2) {
    printf("N needs to be at least 3.\n");
    return false;
  }
  return true;
}

// returns the size of the problem
bool LuksanVlcek1::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in LuksanVlcek1.hpp has 4 variables, x[0] through x[3]
  n = N_;

  m = N_-2;

  nnz_jac_g = m*3;

  nnz_h_lag = n + n-1;

  // use the C style numbering of matrix indices (starting at 0)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool LuksanVlcek1::get_bounds_info(Index n, Number* x_l, Number* x_u,
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
bool LuksanVlcek1::get_starting_point(Index n, bool init_x, Number* x,
                                      bool init_z, Number* z_L, Number* z_U,
                                      Index m, bool init_lambda,
                                      Number* lambda)
{
  if (!init_x || init_z || init_lambda) {
    return false;
  }

  // set the starting point
  for (Index i=0; i<n/2; i++) {
    x[2*i] = -1.2;
    x[2*i+1] = 1.;
  }
  if (n != 2*(n/2)) {
    x[n-1] = -1.2;
  }

  /*
  // DELETEME
  for (Index i=0; i<n; i++) {
    x[i] += 0.001*((Number) i);
  }
  */

  return true;
}

// returns the value of the objective function
bool LuksanVlcek1::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  obj_value = 0.;
  for (Index i=0; i<N_-1; i++) {
    Number a1 = x[i]*x[i]-x[i+1];
    Number a2 = x[i] - 1.;
    obj_value += 100.*a1*a1 + a2*a2;
  }

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool LuksanVlcek1::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  grad_f[0] = 0.;
  for (Index i=0; i<N_-1; i++) {
    grad_f[i] += 400.*x[i]*(x[i]*x[i]-x[i+1]) + 2.*(x[i]-1.);
    grad_f[i+1] = -200.*(x[i]*x[i]-x[i+1]);
  }

  return true;
}

// return the value of the constraints: g(x)
bool LuksanVlcek1::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  for (Index i=0; i<N_-2; i++) {
    g[i] = 3.*pow(x[i+1],3.) + 2.*x[i+2] - 5.
           + sin(x[i+1]-x[i+2])*sin(x[i+1]+x[i+2]) + 4.*x[i+1]
           - x[i]*exp(x[i]-x[i+1]) - 3.;
  }

  return true;
}

// return the structure or values of the jacobian
bool LuksanVlcek1::eval_jac_g(Index n, const Number* x, bool new_x,
                              Index m, Index nele_jac, Index* iRow, Index *jCol,
                              Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian

    Index ijac=0;
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
  }
  else {
    // return the values of the jacobian of the constraints

    Index ijac=0;
    for (Index i=0; i<N_-2; i++) {
      // x[i]
      values[ijac] = -(1.+x[i])*exp(x[i]-x[i+1]);
      ijac++;
      // x[i+1]
      values[ijac] = 9.*x[i+1]*x[i+1]
                     + cos(x[i+1]-x[i+2])*sin(x[i+1]+x[i+2])
                     + sin(x[i+1]-x[i+2])*cos(x[i+1]+x[i+2])
                     + 4. + x[i]*exp(x[i]-x[i+1]);
      ijac++;
      // x[i+2]
      values[ijac] = 2.
                     - cos(x[i+1]-x[i+2])*sin(x[i+1]+x[i+2])
                     + sin(x[i+1]-x[i+2])*cos(x[i+1]+x[i+2]);
      ijac++;
    }
  }

  return true;
}

//return the structure or values of the hessian
bool LuksanVlcek1::eval_h(Index n, const Number* x, bool new_x,
                          Number obj_factor, Index m, const Number* lambda,
                          bool new_lambda, Index nele_hess, Index* iRow,
                          Index* jCol, Number* values)
{
  if (values == NULL) {
    Index ihes=0;
    for (Index i=0; i<N_; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i;
      ihes++;
      if (i<N_-1) {
        iRow[ihes] = i;
        jCol[ihes] = i+1;
        ihes++;
      }
    }
    DBG_ASSERT(ihes == nele_hess);
  }
  else {
    Index ihes=0;
    for (Index i=0; i<N_; i++) {
      // x[i],x[i]
      if (i<N_-1) {
        values[ihes] = obj_factor*(2.+400.*(3.*x[i]*x[i]-x[i+1]));
        if (i<N_-2) {
          values[ihes] -= lambda[i]*(2.+x[i])*exp(x[i]-x[i+1]);
        }
      }
      else {
        values[ihes] = 0.;
      }
      if (i>0) {
        // x[i+1]x[i+1]
        values[ihes] += obj_factor*200.;
        if (i<N_-1) {
          values[ihes] += lambda[i-1]*(18.*x[i]
                                       - 2.*sin(x[i]-x[i+1])*sin(x[i]+x[i+1])
                                       + 2.*cos(x[i]-x[i+1])*cos(x[i]+x[i+1])
                                       - x[i-1]*exp(x[i-1]-x[i]));
        }
      }
      if (i>1) {
        // x[i+2]x[i+2]
        values[ihes] +=
          lambda[i-2]*(- 2.*sin(x[i-1]-x[i])*sin(x[i-1]+x[i])
                       - 2.*cos(x[i-1]-x[i])*cos(x[i-1]+x[i]));
      }
      ihes++;

      if (i<N_-1) {
        // x[i],x[i+1]
        values[ihes] = obj_factor*(-400.*x[i]);
        if (i<N_-2) {
          values[ihes] += lambda[i]*(1.+x[i])*exp(x[i]-x[i+1]);
        }
        /*
        if (i>0) {
        // x[i+1],x[i+2]
        values[ihes] +=
        lambda[i-1]*(  sin(x[i]-x[i+1])*sin(x[i]+x[i+1])
           + cos(x[i]-x[i+1])*cos(x[i]+x[i+1])
           - cos(x[i]-x[i+1])*cos(x[i]+x[i+1])
           - sin(x[i]-x[i+1])*sin(x[i]+x[i+1])
        );
        }
        */
        ihes++;
      }
    }
    DBG_ASSERT(ihes == nele_hess);
  }

  return true;
}

void LuksanVlcek1::finalize_solution(SolverReturn status,
                                     Index n, const Number* x, const Number* z_L, const Number* z_U,
                                     Index m, const Number* g, const Number* lambda,
                                     Number obj_value,
				     const IpoptData* ip_data,
				     IpoptCalculatedQuantities* ip_cq)
{}

