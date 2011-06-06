// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter           IBM    2005-10-18
//                  based on MyNLP.hpp

#ifndef __MITTELMANNDISTRCNTRLDIRI_HPP__
#define __MITTELMANNDISTRCNTRLDIRI_HPP__

#include "IpTNLP.hpp"
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

/** Base class for distributed control problems with Dirichlet
 *  boundary conditions, as formulated by Hans Mittelmann as Examples
 *  1-3 in "Optimization Techniques for Solving Elliptic Control
 *  Problems with Control and State Constraints. Part 2: Distributed
 *  Control"
 */
class MittelmannDistCntrlDiriBase : public RegisteredTNLP
{
public:
  /** Constructor.  N is the number of mesh points in one dimension
   *  (excluding boundary). */
  MittelmannDistCntrlDiriBase();

  /** Default destructor */
  virtual ~MittelmannDistCntrlDiriBase();

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
                                 Number obj_value,
				 const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq);
  //@}

protected:
  /** Method for setting the internal parameters that define the
   *  problem. It must be called by the child class in its
   *  implementation of InitializeParameters. */
  void SetBaseParameters(Index N, Number alpha, Number lb_y,
                         Number ub_y, Number lb_u, Number ub_u,
                         Number u_init);

  /**@name Functions that defines a particular instance. */
  //@{
  /** Target profile function for y */
  virtual Number y_d_cont(Number x1, Number x2) const =0;
  /** Forcing function for the elliptic equation */
  virtual Number d_cont(Number x1, Number x2, Number y, Number u) const =0;
  /** First partial derivative of forcing function w.r.t. y */
  virtual Number d_cont_dy(Number x1, Number x2, Number y, Number u) const =0;
  /** First partial derivative of forcing function w.r.t. u */
  virtual Number d_cont_du(Number x1, Number x2, Number y, Number u) const =0;
  /** Second partial derivative of forcing function w.r.t. y,y */
  virtual Number d_cont_dydy(Number x1, Number x2, Number y, Number u) const =0;
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
  MittelmannDistCntrlDiriBase(const MittelmannDistCntrlDiriBase&);
  MittelmannDistCntrlDiriBase& operator=(const MittelmannDistCntrlDiriBase&);
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
  /** Initial value for the constrols u */
  Number u_init_;
  /** Weighting parameter for the control target deviation functional
   *  in the objective */
  Number alpha_;
  /** Array for the target profile for y */
  Number* y_d_;
  //@}

  /**@name Auxilliary methods */
  //@{
  /** Translation of mesh point indices to NLP variable indices for
   *  y(x_ij) */
  inline Index y_index(Index i, Index j) const
  {
    return j + (N_+2)*i;
  }
  /** Translation of mesh point indices to NLP variable indices for
   *  u(x_ij) */
  inline Index u_index(Index i, Index j) const
  {
    return (N_+2)*(N_+2) + (j-1) + (N_)*(i-1);
  }
  /** Translation of interior mesh point indices to the corresponding
   *  PDE constraint number */
  inline Index pde_index(Index i, Index j) const
  {
    return (j-1) + N_*(i-1);
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
  //@}
};

/** Class implementating Example 1 */
class MittelmannDistCntrlDiri1 : public MittelmannDistCntrlDiriBase
{
public:
  MittelmannDistCntrlDiri1()
  {}

  virtual ~MittelmannDistCntrlDiri1()
  {}

  virtual bool InitializeProblem(Index N)
  {
    if (N<1) {
      printf("N has to be at least 1.");
      return false;
    }
    Number alpha = 0.001;
    Number lb_y = -1e20;
    Number ub_y = 0.185;
    Number lb_u = 1.5;
    Number ub_u = 4.5;
    Number u_init = (ub_u+lb_u)/2.;

    SetBaseParameters(N, alpha, lb_y, ub_y, lb_u, ub_u, u_init);
    return true;
  }
protected:
  /** Target profile function for y */
  virtual Number y_d_cont(Number x1, Number x2)  const
  {
    return 1. + 2.*(x1*(x1-1.)+x2*(x2-1.));
  }
  /** Forcing function for the elliptic equation */
  virtual Number d_cont(Number x1, Number x2, Number y, Number u)  const
  {
    return pow(y,3) - y - u;
  }
  /** First partial derivative of forcing function w.r.t. y */
  virtual Number d_cont_dy(Number x1, Number x2, Number y, Number u)  const
  {
    return 3.*y*y - 1.;
  }
  /** First partial derivative of forcing function w.r.t. u */
  virtual Number d_cont_du(Number x1, Number x2, Number y, Number u)  const
  {
    return -1.;
  }
  /** Second partial derivative of forcing function w.r.t y,y */
  virtual Number d_cont_dydy(Number x1, Number x2, Number y, Number u)  const
  {
    return 6.*y;
  }
private:
  /**@name hide implicitly defined contructors copy operators */
  //@{
  MittelmannDistCntrlDiri1(const MittelmannDistCntrlDiri1&);
  MittelmannDistCntrlDiri1& operator=(const MittelmannDistCntrlDiri1&);
  //@}
};

/** Class implementating Example 2 */
class MittelmannDistCntrlDiri2 : public MittelmannDistCntrlDiriBase
{
public:
  MittelmannDistCntrlDiri2()
  {}
  virtual ~MittelmannDistCntrlDiri2()
  {}
  virtual bool InitializeProblem(Index N)
  {
    if (N<1) {
      printf("N has to be at least 1.");
      return false;
    }
    Number alpha = 0.;
    Number lb_y = -1e20;
    Number ub_y = 0.185;
    Number lb_u = 1.5;
    Number ub_u = 4.5;
    Number u_init = (ub_u+lb_u)/2.;

    SetBaseParameters(N, alpha, lb_y, ub_y, lb_u, ub_u, u_init);
    return true;
  }
protected:
  /** Target profile function for y */
  virtual Number y_d_cont(Number x1, Number x2)  const
  {
    return 1. + 2.*(x1*(x1-1.)+x2*(x2-1.));
  }
  /** Forcing function for the elliptic equation */
  virtual Number d_cont(Number x1, Number x2, Number y, Number u)  const
  {
    return pow(y,3) - y - u;
  }
  /** First partial derivative of forcing function w.r.t. y */
  virtual Number d_cont_dy(Number x1, Number x2, Number y, Number u)  const
  {
    return 3.*y*y - 1.;
  }
  /** First partial derivative of forcing function w.r.t. u */
  virtual Number d_cont_du(Number x1, Number x2, Number y, Number u)  const
  {
    return -1.;
  }
  /** Second partial derivative of forcing function w.r.t y,y */
  virtual Number d_cont_dydy(Number x1, Number x2, Number y, Number u)  const
  {
    return 6.*y;
  }
private:
  /**@name hide implicitly defined contructors copy operators */
  //@{
  MittelmannDistCntrlDiri2(const MittelmannDistCntrlDiri2&);
  MittelmannDistCntrlDiri2& operator=(const MittelmannDistCntrlDiri2&);
  //@}
};

/** Class implementating Example 3 */
class MittelmannDistCntrlDiri3 : public MittelmannDistCntrlDiriBase
{
public:
  MittelmannDistCntrlDiri3()
      :
      pi_(4.*atan(1.))
  {}
  virtual ~MittelmannDistCntrlDiri3()
  {}
  virtual bool InitializeProblem(Index N)
  {
    if (N<1) {
      printf("N has to be at least 1.");
      return false;
    }
    Number alpha = 0.001;
    Number lb_y = -1e20;
    Number ub_y = 0.11;
    Number lb_u = -5;
    Number ub_u = 5.;
    Number u_init = (ub_u+lb_u)/2.;

    SetBaseParameters(N, alpha, lb_y, ub_y, lb_u, ub_u, u_init);
    return true;
  }
protected:
  /** Target profile function for y */
  virtual Number y_d_cont(Number x1, Number x2)  const
  {
    return sin(2.*pi_*x1)*sin(2.*pi_*x2);
  }
  /** Forcing function for the elliptic equation */
  virtual Number d_cont(Number x1, Number x2, Number y, Number u)  const
  {
    return -exp(y) - u;
  }
  /** First partial derivative of forcing function w.r.t. y */
  virtual Number d_cont_dy(Number x1, Number x2, Number y, Number u)  const
  {
    return -exp(y);
  }
  /** First partial derivative of forcing function w.r.t. u */
  virtual Number d_cont_du(Number x1, Number x2, Number y, Number u)  const
  {
    return -1.;
  }
  /** Second partial derivative of forcing function w.r.t y,y */
  virtual Number d_cont_dydy(Number x1, Number x2, Number y, Number u)  const
  {
    return -exp(y);
  }
private:
  /**@name hide implicitly defined contructors copy operators */
  //@{
  MittelmannDistCntrlDiri3(const MittelmannDistCntrlDiri3&);
  MittelmannDistCntrlDiri3& operator=(const MittelmannDistCntrlDiri3&);
  //@}
  /** Value of pi (made available for convenience) */
  const Number pi_;
};

class MittelmannDistCntrlDiri3a : public MittelmannDistCntrlDiriBase
{
public:
  MittelmannDistCntrlDiri3a()
      :
      pi_(4.*atan(1.))
  {}
  virtual ~MittelmannDistCntrlDiri3a()
  {}
  virtual bool InitializeProblem(Index N)
  {
    if (N<1) {
      printf("N has to be at least 1.");
      return false;
    }
    Number alpha = 0.;
    Number lb_y = -1e20;
    Number ub_y = 0.11;
    Number lb_u = -5;
    Number ub_u = 5.;
    Number u_init = (ub_u+lb_u)/2.;

    SetBaseParameters(N, alpha, lb_y, ub_y, lb_u, ub_u, u_init);
    return true;
  }
protected:
  /** Target profile function for y */
  virtual Number y_d_cont(Number x1, Number x2)  const
  {
    return sin(2.*pi_*x1)*sin(2.*pi_*x2);
  }
  /** Forcing function for the elliptic equation */
  virtual Number d_cont(Number x1, Number x2, Number y, Number u)  const
  {
    return -exp(y) - u;
  }
  /** First partial derivative of forcing function w.r.t. y */
  virtual Number d_cont_dy(Number x1, Number x2, Number y, Number u)  const
  {
    return -exp(y);
  }
  /** First partial derivative of forcing function w.r.t. u */
  virtual Number d_cont_du(Number x1, Number x2, Number y, Number u)  const
  {
    return -1.;
  }
  /** Second partial derivative of forcing function w.r.t y,y */
  virtual Number d_cont_dydy(Number x1, Number x2, Number y, Number u)  const
  {
    return -exp(y);
  }
private:
  /**@name hide implicitly defined contructors copy operators */
  //@{
  MittelmannDistCntrlDiri3a(const MittelmannDistCntrlDiri3a&);
  MittelmannDistCntrlDiri3a& operator=(const MittelmannDistCntrlDiri3a&);
  //@}
  /** Value of pi (made available for convenience) */
  const Number pi_;
};

#endif
