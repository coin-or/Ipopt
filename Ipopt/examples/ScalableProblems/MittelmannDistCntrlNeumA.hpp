// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter           IBM    2005-10-18
//                  based on MyNLP.hpp

#ifndef __MITTELMANNDISTRCNTRLNEUMA_HPP__
#define __MITTELMANNDISTRCNTRLNEUMA_HPP__

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

/** Base class for distributed control problems with homogeneous
 *  Neumann boundary conditions, as formulated by Hans Mittelmann as
 *  Examples 4-6 in "Optimization Techniques for Solving Elliptic
 *  Control Problems with Control and State Constraints. Part 2:
 *  Distributed Control"
 */
class MittelmannDistCntrlNeumABase : public RegisteredTNLP
{
public:
  /** Constructor.  N is the number of mesh points in one dimension
   *  (excluding boundary). */
  MittelmannDistCntrlNeumABase();

  /** Default destructor */
  virtual ~MittelmannDistCntrlNeumABase();

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
  void SetBaseParameters(Index N, Number lb_y,
                         Number ub_y, Number lb_u, Number ub_u,
                         Number b_0j, Number b_1j, Number b_i0, Number b_i1,
                         Number u_init);

  /**@name Functions that defines a particular instance. */
  //@{
  /** Target profile function for y (and initial guess function) */
  virtual Number y_d_cont(Number x1, Number x2) const =0;
  /** Integrant in objective function */
  virtual Number fint_cont(Number x1, Number x2, Number y, Number u) const =0;
  /** First partial derivative of fint_cont w.r.t. y */
  virtual Number fint_cont_dy(Number x1, Number x2, Number y, Number u) const =0;
  /** First partial derivative of fint_cont w.r.t. u */
  virtual Number fint_cont_du(Number x1, Number x2, Number y, Number u) const =0;
  /** Second partial derivative of fint_cont w.r.t. y,y */
  virtual Number fint_cont_dydy(Number x1, Number x2, Number y, Number u) const =0;
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. y,y is always zero. */
  virtual bool fint_cont_dydy_alwayszero() const =0;
  /** Second partial derivative of fint_cont w.r.t. u,u */
  virtual Number fint_cont_dudu(Number x1, Number x2, Number y, Number u) const =0;
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. u,u is always zero. */
  virtual bool fint_cont_dudu_alwayszero() const =0;
  /** Second partial derivative of fint_cont w.r.t. y,u */
  virtual Number fint_cont_dydu(Number x1, Number x2, Number y, Number u) const =0;
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. y,u is always zero. */
  virtual bool fint_cont_dydu_alwayszero() const =0;
  /** Forcing function for the elliptic equation */
  virtual Number d_cont(Number x1, Number x2, Number y, Number u) const =0;
  /** First partial derivative of forcing function w.r.t. y */
  virtual Number d_cont_dy(Number x1, Number x2, Number y, Number u) const =0;
  /** First partial derivative of forcing function w.r.t. u */
  virtual Number d_cont_du(Number x1, Number x2, Number y, Number u) const =0;
  /** Second partial derivative of forcing function w.r.t. y,y */
  virtual Number d_cont_dydy(Number x1, Number x2, Number y, Number u) const =0;
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,y is always zero. */
  virtual bool d_cont_dydy_alwayszero() const =0;
  /** Second partial derivative of forcing function w.r.t. u,u */
  virtual Number d_cont_dudu(Number x1, Number x2, Number y, Number u) const =0;
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,y is always zero. */
  virtual bool d_cont_dudu_alwayszero() const =0;
  /** Second partial derivative of forcing function w.r.t. y,u */
  virtual Number d_cont_dydu(Number x1, Number x2, Number y, Number u) const =0;
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,u is always zero. */
  virtual bool d_cont_dydu_alwayszero() const =0;
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
  MittelmannDistCntrlNeumABase(const MittelmannDistCntrlNeumABase&);
  MittelmannDistCntrlNeumABase& operator=(const MittelmannDistCntrlNeumABase&);
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
  /** Value of beta function (in Neumann boundary condition) for
   *  (0,x2) bounray */
  Number b_0j_;
  /** Value of beta function (in Neumann boundary condition) for
   *  (1,x2) bounray */
  Number b_1j_;
  /** Value of beta function (in Neumann boundary condition) for
   *  (x1,0) bounray */
  Number b_i0_;
  /** Value of beta function (in Neumann boundary condition) for
   *  (x1,1) bounray */
  Number b_i1_;
  /** Initial value for the constrols u */
  Number u_init_;
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

/** Class implementating Example 4 */
class MittelmannDistCntrlNeumA1 : public MittelmannDistCntrlNeumABase
{
public:
  MittelmannDistCntrlNeumA1()
      :
      pi_(4.*atan(1.)),
      alpha_(0.001)
  {}

  virtual ~MittelmannDistCntrlNeumA1()
  {}

  virtual bool InitializeProblem(Index N)
  {
    if (N<1) {
      printf("N has to be at least 1.");
      return false;
    }
    Number lb_y = -1e20;
    Number ub_y = 0.371;
    Number lb_u = -8.;
    Number ub_u = 9.;
    Number b_0j = 1.;
    Number b_1j = 1.;
    Number b_i0 = 1.;
    Number b_i1 = 1.;
    Number u_init = (ub_u+lb_u)/2.;

    SetBaseParameters(N, lb_y, ub_y, lb_u, ub_u, b_0j, b_1j, b_i0, b_i1, u_init);
    return true;
  }
protected:
  /** Target profile function for y */
  virtual Number y_d_cont(Number x1, Number x2)  const
  {
    return sin(2.*pi_*x1)*sin(2.*pi_*x2);
  }
  /** Integrant in objective function */
  virtual Number fint_cont(Number x1, Number x2, Number y, Number u) const
  {
    Number diff_y = y-y_d_cont(x1,x2);
    return 0.5*(diff_y*diff_y + alpha_*u*u);
  }
  /** First partial derivative of fint_cont w.r.t. y */
  virtual Number fint_cont_dy(Number x1, Number x2, Number y, Number u) const
  {
    return  y-y_d_cont(x1,x2);
  }

  /** First partial derivative of fint_cont w.r.t. u */
  virtual Number fint_cont_du(Number x1, Number x2, Number y, Number u) const
  {
    return alpha_*u;
  }
  /** Second partial derivative of fint_cont w.r.t. y,y */
  virtual Number fint_cont_dydy(Number x1, Number x2, Number y, Number u) const
  {
    return 1.;
  }
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. y,y is always zero. */
  virtual bool fint_cont_dydy_alwayszero() const
  {
    return false;
  }
  /** Second partial derivative of fint_cont w.r.t. u,u */
  virtual Number fint_cont_dudu(Number x1, Number x2, Number y, Number u) const
  {
    return alpha_;
  }
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. u,u is always zero. */
  virtual bool fint_cont_dudu_alwayszero() const
  {
    return false;
  }
  /** Second partial derivative of fint_cont w.r.t. y,u */
  virtual Number fint_cont_dydu(Number x1, Number x2, Number y, Number u) const
  {
    return 0.;
  }
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. y,u is always zero. */
  virtual bool fint_cont_dydu_alwayszero() const
  {
    return true;
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
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,y is always zero. */
  virtual bool d_cont_dydy_alwayszero() const
  {
    return false;
  }
  /** Second partial derivative of forcing function w.r.t. u,u */
  virtual Number d_cont_dudu(Number x1, Number x2, Number y, Number u) const
  {
    return 0.;
  }
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,y is always zero. */
  virtual bool d_cont_dudu_alwayszero() const
  {
    return true;
  }
  /** Second partial derivative of forcing function w.r.t. y,u */
  virtual Number d_cont_dydu(Number x1, Number x2, Number y, Number u) const
  {
    return 0.;
  }
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,u is always zero. */
  virtual bool d_cont_dydu_alwayszero() const
  {
    return true;
  }
private:
  /**@name hide implicitly defined contructors copy operators */
  //@{
  MittelmannDistCntrlNeumA1(const MittelmannDistCntrlNeumA1&);
  MittelmannDistCntrlNeumA1& operator=(const MittelmannDistCntrlNeumA1&);
  //@}
  /** Value of pi (made available for convenience) */
  const Number pi_;
  /** Value for parameter alpha in objective functin */
  const Number alpha_;
};

/** Class implementating Example 5 */
class MittelmannDistCntrlNeumA2 : public MittelmannDistCntrlNeumABase
{
public:
  MittelmannDistCntrlNeumA2()
      :
      pi_(4.*atan(1.))
  {}

  virtual ~MittelmannDistCntrlNeumA2()
  {}

  virtual bool InitializeProblem(Index N)
  {
    if (N<1) {
      printf("N has to be at least 1.");
      return false;
    }
    Number lb_y = -1e20;
    Number ub_y = 0.371;
    Number lb_u = -8.;
    Number ub_u = 9.;
    Number b_0j = 1.;
    Number b_1j = 1.;
    Number b_i0 = 1.;
    Number b_i1 = 1.;
    Number u_init = (ub_u+lb_u)/2.;

    SetBaseParameters(N, lb_y, ub_y, lb_u, ub_u, b_0j, b_1j, b_i0, b_i1, u_init);
    return true;
  }
protected:
  /** Target profile function for y */
  virtual Number y_d_cont(Number x1, Number x2)  const
  {
    return sin(2.*pi_*x1)*sin(2.*pi_*x2);
  }
  /** Integrant in objective function */
  virtual Number fint_cont(Number x1, Number x2, Number y, Number u) const
  {
    Number diff_y = y-y_d_cont(x1,x2);
    return 0.5*diff_y*diff_y;
  }
  /** First partial derivative of fint_cont w.r.t. y */
  virtual Number fint_cont_dy(Number x1, Number x2, Number y, Number u) const
  {
    return  y-y_d_cont(x1,x2);
  }

  /** First partial derivative of fint_cont w.r.t. u */
  virtual Number fint_cont_du(Number x1, Number x2, Number y, Number u) const
  {
    return 0.;
  }
  /** Second partial derivative of fint_cont w.r.t. y,y */
  virtual Number fint_cont_dydy(Number x1, Number x2, Number y, Number u) const
  {
    return 1.;
  }
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. y,y is always zero. */
  virtual bool fint_cont_dydy_alwayszero() const
  {
    return false;
  }
  /** Second partial derivative of fint_cont w.r.t. u,u */
  virtual Number fint_cont_dudu(Number x1, Number x2, Number y, Number u) const
  {
    return 0.;
  }
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. u,u is always zero. */
  virtual bool fint_cont_dudu_alwayszero() const
  {
    return true;
  }
  /** Second partial derivative of fint_cont w.r.t. y,u */
  virtual Number fint_cont_dydu(Number x1, Number x2, Number y, Number u) const
  {
    return 0.;
  }
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. y,u is always zero. */
  virtual bool fint_cont_dydu_alwayszero() const
  {
    return true;
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
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,y is always zero. */
  virtual bool d_cont_dydy_alwayszero() const
  {
    return false;
  }
  /** Second partial derivative of forcing function w.r.t. u,u */
  virtual Number d_cont_dudu(Number x1, Number x2, Number y, Number u) const
  {
    return 0.;
  }
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,y is always zero. */
  virtual bool d_cont_dudu_alwayszero() const
  {
    return true;
  }
  /** Second partial derivative of forcing function w.r.t. y,u */
  virtual Number d_cont_dydu(Number x1, Number x2, Number y, Number u) const
  {
    return 0.;
  }
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,u is always zero. */
  virtual bool d_cont_dydu_alwayszero() const
  {
    return true;
  }
private:
  /**@name hide implicitly defined contructors copy operators */
  //@{
  MittelmannDistCntrlNeumA2(const MittelmannDistCntrlNeumA2&);
  MittelmannDistCntrlNeumA2& operator=(const MittelmannDistCntrlNeumA2&);
  //@}
  /** Value of pi (made available for convenience) */
  const Number pi_;
};

/** Class implementating Example 6 */
class MittelmannDistCntrlNeumA3 : public MittelmannDistCntrlNeumABase
{
public:
  MittelmannDistCntrlNeumA3()
      :
      pi_(4.*atan(1.)),
      M_(1.),
      K_(0.8),
      b_(1.)
  {}

  virtual ~MittelmannDistCntrlNeumA3()
  {}

  virtual bool InitializeProblem(Index N)
  {
    if (N<1) {
      printf("N has to be at least 1.");
      return false;
    }
    Number lb_y = 3.;//-1e20;
    Number ub_y = 6.09;
    Number lb_u = 1.4;
    Number ub_u = 1.6;
    Number b_0j = 1.;
    Number b_1j = 0.;
    Number b_i0 = 1.;
    Number b_i1 = 0.;
    Number u_init = (ub_u+lb_u)/2.;

    SetBaseParameters(N, lb_y, ub_y, lb_u, ub_u, b_0j, b_1j, b_i0, b_i1, u_init);
    return true;
  }
protected:
  /** Profile function for initial y */
  virtual Number y_d_cont(Number x1, Number x2)  const
  {
    return 6.;
  }
  /** Integrant in objective function */
  virtual Number fint_cont(Number x1, Number x2, Number y, Number u) const
  {
    return u*(M_*u - K_*y);
  }
  /** First partial derivative of fint_cont w.r.t. y */
  virtual Number fint_cont_dy(Number x1, Number x2, Number y, Number u) const
  {
    return -K_*u;
  }

  /** First partial derivative of fint_cont w.r.t. u */
  virtual Number fint_cont_du(Number x1, Number x2, Number y, Number u) const
  {
    return 2.*M_*u - K_*y;
  }
  /** Second partial derivative of fint_cont w.r.t. y,y */
  virtual Number fint_cont_dydy(Number x1, Number x2, Number y, Number u) const
  {
    return 0.;
  }
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. y,y is always zero. */
  virtual bool fint_cont_dydy_alwayszero() const
  {
    return true;
  }
  /** Second partial derivative of fint_cont w.r.t. u,u */
  virtual Number fint_cont_dudu(Number x1, Number x2, Number y, Number u) const
  {
    return 2.*M_;
  }
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. u,u is always zero. */
  virtual bool fint_cont_dudu_alwayszero() const
  {
    return false;
  }
  /** Second partial derivative of fint_cont w.r.t. y,u */
  virtual Number fint_cont_dydu(Number x1, Number x2, Number y, Number u) const
  {
    return -K_;
  }
  /** returns true if second partial derivative of fint_cont
   *  w.r.t. y,u is always zero. */
  virtual bool fint_cont_dydu_alwayszero() const
  {
    return false;
  }
  /** Forcing function for the elliptic equation */
  virtual Number d_cont(Number x1, Number x2, Number y, Number u)  const
  {
    return y*(u + b_*y - a(x1,x2));
  }
  /** First partial derivative of forcing function w.r.t. y */
  virtual Number d_cont_dy(Number x1, Number x2, Number y, Number u)  const
  {
    return (u + 2.*b_*y -a(x1,x2));
  }
  /** First partial derivative of forcing function w.r.t. u */
  virtual Number d_cont_du(Number x1, Number x2, Number y, Number u)  const
  {
    return y;
  }
  /** Second partial derivative of forcing function w.r.t y,y */
  virtual Number d_cont_dydy(Number x1, Number x2, Number y, Number u)  const
  {
    return 2.*b_;
  }
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,y is always zero. */
  virtual bool d_cont_dydy_alwayszero() const
  {
    return false;
  }
  /** Second partial derivative of forcing function w.r.t. u,u */
  virtual Number d_cont_dudu(Number x1, Number x2, Number y, Number u) const
  {
    return 0.;
  }
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,y is always zero. */
  virtual bool d_cont_dudu_alwayszero() const
  {
    return true;
  }
  /** Second partial derivative of forcing function w.r.t. y,u */
  virtual Number d_cont_dydu(Number x1, Number x2, Number y, Number u) const
  {
    return 1.;
  }
  /** returns true if second partial derivative of d_cont
   *  w.r.t. y,u is always zero. */
  virtual bool d_cont_dydu_alwayszero() const
  {
    return false;
  }
private:
  /**@name hide implicitly defined contructors copy operators */
  //@{
  MittelmannDistCntrlNeumA3(const MittelmannDistCntrlNeumA3&);
  MittelmannDistCntrlNeumA3& operator=(const MittelmannDistCntrlNeumA3&);
  //@}
  /** Value of pi (made available for convenience) */
  const Number pi_;
  /*@name constrants appearing in problem formulation */
  //@{
  const Number M_;
  const Number K_;
  const Number b_;
  //@}
  //* Auxiliary function for state equation */
  inline Number a(Number x1, Number x2) const
  {
    return 7. + 4.*sin(2.*pi_*x1*x2);
  }
};

#endif
