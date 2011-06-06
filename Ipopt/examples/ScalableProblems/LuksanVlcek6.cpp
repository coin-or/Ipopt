// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2005-10-127

#include "LuksanVlcek6.hpp"

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

LuksanVlcek6::LuksanVlcek6(Number g_l, Number g_u)
{
  g_l_ = g_l;
  g_u_ = g_u;
}

bool LuksanVlcek6::InitializeProblem(Index N)
{
  N_=N;
  if (2*(N_/2)!=N_ || N<=1) {
    printf("N needs to be even and at least 2.\n");
    return false;
  }
  return true;
}

// returns the size of the problem
bool LuksanVlcek6::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in LuksanVlcek6.hpp has 4 variables, x[0] through x[3]
  n = N_+1;

  m = N_/2;

  nnz_jac_g = 3 * m;

  nnz_h_lag = n + n-1 + n-2 + n-3 + n-4 + n-5 + n-6;

  // use the C style numbering of matrix indices (starting at 0)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool LuksanVlcek6::get_bounds_info(Index n, Number* x_l, Number* x_u,
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
bool LuksanVlcek6::get_starting_point(Index n, bool init_x, Number* x,
                                      bool init_z, Number* z_L, Number* z_U,
                                      Index m, bool init_lambda,
                                      Number* lambda)
{
  if (!init_x || init_z || init_lambda) {
    return false;
  }

  // set the starting point
  for (Index i=0; i<n; i++) {
    x[i] = 3.;
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
bool LuksanVlcek6::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  Number p = 7./3.;
  obj_value = 0.;
  for (Index i=0; i<N_; i++) {
    Number b = (2.+5.*x[i]*x[i])*x[i] + 1.;
    for (Index j=Max(0,i-5); j<=Min(N_-1,i+1); j++) {
      b += x[j]*(1.+x[j]);
    }
    obj_value += pow(fabs(b),p);
  }

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool LuksanVlcek6::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  Number p = 7./3.;
  grad_f[0] = 0.;
  for (Index i=0; i<N_; i++) {
    Number b = (2.+5.*x[i]*x[i])*x[i] + 1.;
    for (Index j=Max(0,i-5); j<=Min(N_-1,i+1); j++) {
      b += x[j]*(1.+x[j]);
    }
    Number pb1 = pow(fabs(b),p-1.);
    Number apb1 = p*Sgn(b)*pb1;
    for (Index j=Max(0,i-5); j<i; j++) {
      grad_f[j] += (1.+2.*x[j])*apb1;
    }
    grad_f[i] += (3. + 2.*x[i] + 15.*x[i]*x[i])*apb1;
    if (i<N_-1) {
      grad_f[i+1] = (1.+2.*x[i+1])*apb1;
    }
  }
  grad_f[n-1] = 0.;

  return true;
}

// return the value of the constraints: g(x)
bool LuksanVlcek6::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  for (Index i=0; i<N_/2; i++) {
    Number e = exp(x[2*i] - x[2*i+1] - x[2*i+2]);
    g[i] = 4.*x[2*i+1] - (x[2*i] - x[2*i+2])*e - 3;
  }
  return true;
}

// return the structure or values of the jacobian
bool LuksanVlcek6::eval_jac_g(Index n, const Number* x, bool new_x,
                              Index m, Index nele_jac, Index* iRow, Index *jCol,
                              Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian

    Index ijac = 0;
    for (Index i=0; i<N_/2; i++) {
      iRow[ijac] = i;
      jCol[ijac] = 2*i;
      ijac++;
      iRow[ijac] = i;
      jCol[ijac] = 2*i+1;
      ijac++;
      iRow[ijac] = i;
      jCol[ijac] = 2*i+2;
      ijac++;
    }
    DBG_ASSERT(ijac == nele_jac);
  }
  else {
    // return the values of the jacobian of the constraints

    Index ijac = 0;
    for (Index i=0; i<N_/2; i++) {
      Number e = exp(x[2*i] - x[2*i+1] - x[2*i+2]);
      Number a1 = (1. + x[2*i] - x[2*i+2])*e;
      values[ijac] = -a1;
      ijac++;
      values[ijac] = 4. + (x[2*i] - x[2*i+2])*e;
      ijac++;
      values[ijac] = a1;
      ijac++;
    }
  }

  return true;
}

//return the structure or values of the hessian
bool LuksanVlcek6::eval_h(Index n, const Number* x, bool new_x,
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
    // 1st off-diagonal
    for (Index i=0; i<n-1; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i+1;
      ihes++;
    }
    // 2nd off-diagonal
    for (Index i=0; i<n-2; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i+2;
      ihes++;
    }
    // 3rd off-diagonal
    for (Index i=0; i<n-3; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i+3;
      ihes++;
    }
    // 4th off-diagonal
    for (Index i=0; i<n-4; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i+4;
      ihes++;
    }
    // 5th off-diagonal
    for (Index i=0; i<n-5; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i+5;
      ihes++;
    }
    // 6th off-diagonal
    for (Index i=0; i<n-6; i++) {
      iRow[ihes] = i;
      jCol[ihes] = i+6;
      ihes++;
    }
    DBG_ASSERT(ihes == nele_hess);
  }
  else {
    Number p = 7./3.;
    // First the objective
    values[0] = 0.;
    for (Index i=0; i<N_; i++) {
      Number b = (2.+5.*x[i]*x[i])*x[i] + 1.;
      for (Index j=Max(0,i-5); j<=Min(N_-1,i+1); j++) {
        b += x[j]*(1.+x[j]);
      }
      Number pb1 = pow(fabs(b),p-1.);
      Number pb2 = pow(fabs(b),p-2.);
      Number apb1 = p*Sgn(b)*pb1;
      Number apb2 = p*(p-1.)*pb2;
      Number a1 = 3. + 2.*x[i] + 15.*x[i]*x[i];
      Number a1apb2 = apb2*a1;

      // x[i] x[i]
      values[i] += obj_factor*(a1apb2*a1 + apb1*(2. + 30.*x[i]));
      // x[i] x[j] j<i
      Index ih = n;
      Index ip = n-1;
      for (Index j=i-1; j>=Max(0,i-5); j--) {
        values[ih+j] += obj_factor*a1apb2*(1.+2.*x[j]);
        ih += ip;
        ip--;
      }
      for (Index j=Max(0,i-5); j<i; j++) {
        Number axj = (1.+2.*x[j]);
        // x[j] x[j] j<i
        values[j] += obj_factor*(apb2*axj*axj + 2.*apb1);
        ih = n;
        ip = n-1;
        // x[j] x[k] j<i, k<j
        for (Index k=j-1; k>=Max(0,i-5); k--) {
          values[ih+k] += obj_factor*apb2*(1.+2.*x[k])*axj;
          ih += ip;
          ip--;
        }
      }
      if (i<N_-1) {
        Number axj = (1.+2.*x[i+1]);
        // x[i] x[i+1]
        values[n+i] = obj_factor*a1apb2*axj;
        // x[j=i+1] x[k] k<i
        ih = n+n-1;
        ip = n-2;
        for (Index k=i-1; k>=Max(0,i-5); k--) {
          values[ih+k] = obj_factor*apb2*(1.+2.*x[k])*axj;
          ih += ip;
          ip--;
        }
        // x[j=i+1] x[k=i+1]
        values[i+1] = obj_factor*(apb2*axj*axj+ 2.*apb1);
      }
      else {
        // We could just not set those non-zero elements in the
        // structure, but I'm too lazy to do that now
        values[n+i] = 0.;
        ih = n+n-1;
        ip = n-2;
        for (Index k=i-1; k>=Max(0,i-5); k--) {
          values[ih+k] = 0.;
          ih += ip;
          ip--;
        }
        // x[j=i+1] x[k=i+1]
        values[i+1] = 0.;
      }
    }

    // Now the constraints
    for (Index i=0; i<N_/2; i++) {
      Number e = exp(x[2*i] - x[2*i+1] - x[2*i+2]);
      Number a1 = 1. + x[2*i] - x[2*i+2];
      Number a2 = 1. + a1;
      // x[2*i] x[2*i]
      values[2*i] -= lambda[i]*a2*e;
      // x[2*i] x[2*i+1]
      values[n+2*i] += lambda[i]*a1*e;
      // x[2*i] x[2*i+2]
      values[n+(n-1)+2*i] += lambda[i]*a2*e;
      // x[2*i+1] x[2*i+1]
      values[2*i+1] -= lambda[i]*(x[2*i]-x[2*i+2])*e;
      // x[2*i+1] x[2*i+2]
      values[n+2*i+1] -= lambda[i]*a1*e;
      // x[2*i+2] x[2*i+2]
      values[2*i+2] -= lambda[i]*a2*e;
    }
  }
  return true;
}

void LuksanVlcek6::finalize_solution(SolverReturn status,
                                     Index n, const Number* x, const Number* z_L, const Number* z_U,
                                     Index m, const Number* g, const Number* lambda,
                                     Number obj_value,
				     const IpoptData* ip_data,
				     IpoptCalculatedQuantities* ip_cq)
{}

