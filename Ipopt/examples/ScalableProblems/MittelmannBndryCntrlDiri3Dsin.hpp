// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter             IBM    2005-10-18
//           Olaf Schenk   (Univ. of Basel)      2007-08-01
//              modified MittelmannBndryCntrlDiri.hpp for 3-dim problem

#ifndef __MITTELMANNBNDRYCNTRLDIRI3DSIN_HPP__
#define __MITTELMANNBNDRYCNTRLDIRI3DSIN_HPP__

#include "RegisteredTNLP.hpp"

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

/** Base class for boundary control problems with Dirichlet boundary
 *  conditions, as formulated by Hans Mittelmann as Examples 1-4 in
 *  "Optimization Techniques for Solving Elliptic Control Problems
 *  with Control and State Constraints. Part 2: Boundary Control"
 *
 *  Here, the control variables are identical to the values of y on
 *  the boundary, and therefore we don't need any explicit
 *  optimization variables for u.
 */
class MittelmannBndryCntrlDiriBase3Dsin : public RegisteredTNLP
{
public:
  /** Constructor. */
  MittelmannBndryCntrlDiriBase3Dsin();

  /** Default destructor */
  virtual ~MittelmannBndryCntrlDiriBase3Dsin();

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /** Method for returning scaling parameters */
  virtual bool get_scaling_parameters(Number& obj_scaling,
                                      bool& use_x_scaling, Index n,
                                      Number* x_scaling,
                                      bool& use_g_scaling, Index m,
                                      Number* g_scaling);

  /** @name Solution Methods */
  //@{
  /** This method is called after the optimization, and could write an
   *  output file with the optimal profiles */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_valu,
                                 const IpoptData* ip_data,
                                 IpoptCalculatedQuantities* ip_cq);
  //@}

protected:
  /** Method for setting the internal parameters that define the
   *  problem. It must be called by the child class in its
   *  implementation of InitializeParameters. */
  void SetBaseParameters(Index N, Number alpha, Number lb_y,
                         Number ub_y, Number lb_u, Number ub_u,
                         Number d_const);

  /**@name Functions that defines a particular instance. */
  //@{
  /** Target profile function for y */
  virtual Number y_d_cont(Number x1, Number x2, Number x3) const =0;
  //@}

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  MittelmannBndryCntrlDiriBase3Dsin(const MittelmannBndryCntrlDiriBase3Dsin&);
  MittelmannBndryCntrlDiriBase3Dsin& operator=(const MittelmannBndryCntrlDiriBase3Dsin&);
  //@}

  /**@name Problem specification */
  //@{
  /** Number of mesh points in one dimension (excluding boundary) */
  Index N_;
  /** Step size */
  Number h_;
  /** h_ squaredd */
  Number hh_;
  /** overall lower bound on y */
  Number lb_y_;
  /** overall upper bound on y */
  Number ub_y_;
  /** overall lower bound on u */
  Number lb_u_;
  /** overall upper bound on u */
  Number ub_u_;
  /** Constant value of d appearing in elliptical equation */
  Number d_const_;
  /** Weighting parameter for the control target deviation functional
   *  in the objective */
  Number alpha_;
  /** Array for the target profile for y */
  Number* y_d_;
  //@}

  /**@name Auxilliary methods */
  //@{
  /** Translation of mesh point indices to NLP variable indices for
   *  y(x_ijk) */
  inline Index y_index(Index i, Index j, Index k) const
  {
    return k + (N_+2)*j + (N_+2)*(N_+2)*i;
  }
  /** Translation of interior mesh point indices to the corresponding
   *  PDE constraint number */
  inline Index pde_index(Index i, Index j, Index k) const
  {
    return (k-1) + N_*(j-1) + N_*N_*(i-1);
  }
  /** Compute the grid coordinate for given index in x1 direction */
  inline Number x1_grid(Index i) const
  {
    return h_*(Number)i;
  }
  /** Compute the grid coordinate for given index in x2 direction */
  inline Number x2_grid(Index i) const
  {
    return h_*(Number)i;
  }
  /** Compute the grid coordinate for given index in x3 direction */
  inline Number x3_grid(Index i) const
  {
    return h_*(Number)i;
  }
  //@}
};

/** Class implementating Example 1 */
class MittelmannBndryCntrlDiri3Dsin : public MittelmannBndryCntrlDiriBase3Dsin
{
public:
  MittelmannBndryCntrlDiri3Dsin()
  {}

  virtual ~MittelmannBndryCntrlDiri3Dsin()
  {}

  virtual bool InitializeProblem(Index N)
  {
    if (N<1) {
      printf("N has to be at least 1.");
      return false;
    }
    printf("olaf N %d has to be at least 1.", N);
    Number alpha = 0.1;
    Number lb_y = -1e20;
    Number ub_y = 3.5;
    Number lb_u = 0.;
    Number ub_u = 10.;
    Number d_const = -20.;
    SetBaseParameters(N, alpha, lb_y, ub_y, lb_u, ub_u, d_const);
    return true;
  }
protected:
  /** Target profile function for y */
  virtual Number y_d_cont(Number x1, Number x2, Number x3)  const
  {
    return 3. + 5.*(x1*(x1-1.)*x2*(x2-1.)*x3*(x3-1.));
  }
private:
  /**@name hide implicitly defined contructors copy operators */
  //@{
  MittelmannBndryCntrlDiri3Dsin(const MittelmannBndryCntrlDiri3Dsin&);
  MittelmannBndryCntrlDiri3Dsin& operator=(const MittelmannBndryCntrlDiri3Dsin&);
  //@}

};


#endif
