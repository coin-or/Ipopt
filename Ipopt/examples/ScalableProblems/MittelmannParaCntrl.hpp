// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter           IBM    2005-10-18

#ifndef __MITTELMANNPARACNTRL_HPP__
#define __MITTELMANNPARACNTRL_HPP__

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

using namespace Ipopt;

/** Base class for parabolic and elliptic control problems, as
 *  formulated by Hans Mittelmann as problem (P) in "Sufficient
 *  Optimality for Discretized Parabolic and Elliptic Control
 *  Problems"
 */
template<class T>
class MittelmannParaCntrlBase : public RegisteredTNLP
{
public:
  /** Constructor. */
  MittelmannParaCntrlBase();

  /** Default destructor */
  virtual ~MittelmannParaCntrlBase();

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

  virtual bool InitializeProblem(Index N);

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
  MittelmannParaCntrlBase(const MittelmannParaCntrlBase<T>&);
  MittelmannParaCntrlBase& operator=(const MittelmannParaCntrlBase<T>&);
  //@}

  /**@name Problem specification */
  //@{
  /** Upper bound on t */
  Number T_;
  /** Upper bound on x */
  Number l_;
  /** Number of mesh points in t direction */
  Index Nt_;
  /** Number of mesh points in x direction */
  Index Nx_;
  /** Step size in t direction*/
  Number dt_;
  /** Step size in x direction*/
  Number dx_;
  /** overall lower bound on y */
  Number lb_y_;
  /** overall upper bound on y */
  Number ub_y_;
  /** overall lower bound on u */
  Number lb_u_;
  /** overall upper bound on u */
  Number ub_u_;
  /** Weighting parameter for the control target deviation functional
   *  in the objective */
  Number alpha_;
  /** Weighting parameter in PDE */
  Number beta_;
  /** Array for the target profile for y in objective */
  Number* y_T_;
  /** Array for weighting function a_y in objective */
  Number* a_y_;
  /** Array for weighting function a_u in objective */
  Number* a_u_;
  //@}

  /**@name Auxilliary methods */
  //@{
  /** Translation of mesh point indices to NLP variable indices for
   *  y(x_ij) */
  inline Index y_index(Index jx, Index it) const
  {
    return jx + (Nx_+1)*it;
  }
  inline Index u_index(Index it) const
  {
    return (Nt_+1)*(Nx_+1) + it - 1;
  }
  /** Compute the grid coordinate for given index in t direction */
  inline Number t_grid(Index i) const
  {
    return dt_*(Number)i;
  }
  /** Compute the grid coordinate for given index in x direction */
  inline Number x_grid(Index j) const
  {
    return dx_*(Number)j;
  }
  //@}
};

template <class T>
MittelmannParaCntrlBase<T>::MittelmannParaCntrlBase()
    :
    y_T_(NULL),
    a_y_(NULL),
    a_u_(NULL)
{}

template <class T>
MittelmannParaCntrlBase<T>::~MittelmannParaCntrlBase()
{
  delete [] y_T_;
  delete [] a_y_;
  delete [] a_u_;
}

template <class T>
bool MittelmannParaCntrlBase<T>::
InitializeProblem(Index N)
{
  typename T::ProblemSpecs p;

  if (N<1) {
    printf("N has to be at least 1.");
    return false;
  }

  T_ = p.T();
  l_ = p.l();
  Nt_ = N;
  Nx_ = N;
  dt_ = T_/Nt_;
  dx_ = l_/Nx_;
  lb_y_ = p.lb_y();
  ub_y_ = p.ub_y();
  lb_u_ = p.lb_u();
  ub_u_ = p.ub_u();
  alpha_ = p.alpha();
  beta_ = p.beta();

  y_T_ = new Number[Nx_+1];
  for (Index j=0; j<=Nx_; j++) {
    y_T_[j] = p.y_T(x_grid(j));
  }
  a_y_ = new Number[Nt_];
  for (Index i=1; i<=Nt_; i++) {
    a_y_[i-1] = p.a_y(t_grid(i));
  }
  a_u_ = new Number[Nt_];
  for (Index i=1; i<=Nt_; i++) {
    a_u_[i-1] = p.a_u(t_grid(i));
  }

  return true;
}

template <class T>
bool MittelmannParaCntrlBase<T>::
get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  typename T::ProblemSpecs p;

  n = (Nt_+1)*(Nx_+1) + Nt_;

  m = Nt_*(Nx_-1) + Nt_ + Nt_;

  nnz_jac_g = 6*Nt_*(Nx_-1) + 3*Nt_ + 4*Nt_;

  nnz_h_lag = Nx_+1 + Nt_;
  if (!p.phi_dydy_always_zero()) {
    nnz_h_lag += Nt_;
  }

  index_style = C_STYLE;

  return true;
}

template <class T>
bool MittelmannParaCntrlBase<T>::
get_bounds_info(Index n, Number* x_l, Number* x_u,
                Index m, Number* g_l, Number* g_u)
{
  typename T::ProblemSpecs p;

  // Set overall bounds on the variables
  for (Index jx=0; jx<=Nx_; jx++) {
    for (Index it=1; it<=Nt_; it++) {
      Index iy = y_index(jx,it);
      x_l[iy] = lb_y_;
      x_u[iy] = ub_y_;
    }
  }
  for (Index i=1; i<=Nt_; i++) {
    Index iu = u_index(i);
    x_l[iu] = lb_u_;
    x_u[iu] = ub_u_;
  }

  /*
  // Boundary condition for t=0
  for (Index it=0; it<=Nt_; it++) {
    Index iy = y_index(0,it);
    x_u[iy] = x_l[iy] = p.a(t_grid(it));
  }
  */
  // Boundary condition for t=0
  for (Index jx=0; jx<=Nx_; jx++) {
    Index iy = y_index(jx,0);
    x_u[iy] = x_l[iy] = p.a(x_grid(jx));
  }

  // all discretized PDE constraints have right hand side zero
  for (Index i=0; i<Nt_*(Nx_-1) + Nt_; i++) {
    g_l[i] = 0.;
    g_u[i] = 0.;
  }

  // but we put b on the right hand side for the x=L boundary condition
  for (Index i=0; i<Nt_; i++) {
    g_l[Nt_*(Nx_-1) + Nt_ + i]
    = g_u[Nt_*(Nx_-1) + Nt_ + i]
      = p.b(t_grid(i+1));
  }

  return true;
}

template <class T>
bool MittelmannParaCntrlBase<T>::
get_starting_point(Index n, bool init_x, Number* x,
                   bool init_z, Number* z_L, Number* z_U,
                   Index m, bool init_lambda,
                   Number* lambda)
{
  DBG_ASSERT(init_x==true && init_z==false && init_lambda==false);

  // Set starting point for y
  for (Index jx=0; jx<=Nx_; jx++) {
    for (Index it=0; it<=Nt_; it++) {
      x[y_index(jx,it)] = 0.;
    }
  }
  // for u
  for (Index i=1; i<=Nt_; i++) {
    x[u_index(i)] = (ub_u_+lb_u_)/2.;
  }

  /*
  // DELETEME
  for (Index i=0; i<n; i++) {
    x[i] += 0.01*i;
  }
  */

  return true;
}

template <class T>
bool MittelmannParaCntrlBase<T>::
get_scaling_parameters(Number& obj_scaling,
                       bool& use_x_scaling,
                       Index n, Number* x_scaling,
                       bool& use_g_scaling,
                       Index m, Number* g_scaling)
{
  obj_scaling = 1./Min(dx_,dt_);
  use_x_scaling = false;
  use_g_scaling = false;
  return true;
}

template <class T>
bool MittelmannParaCntrlBase<T>::
eval_f(Index n, const Number* x,
       bool new_x, Number& obj_value)
{
  // Deviation of y from target
  Number a = x[y_index(0,Nt_)] - y_T_[0];
  Number sum = 0.5*a*a;
  for (Index jx=1; jx<Nx_; jx++) {
    a = x[y_index(jx,Nt_)] - y_T_[jx];
    sum += a*a;
  }
  a = x[y_index(Nx_,Nt_)] - y_T_[Nx_];
  sum += 0.5*a*a;
  obj_value = 0.5*dx_*sum;

  // smoothing for control
  if (alpha_!=.0) {
    sum = 0.5*x[u_index(Nt_)]*x[u_index(Nt_)];
    for (Index it=1; it < Nt_; it++) {
      sum += x[u_index(it)]*x[u_index(it)];
    }
    obj_value += 0.5*alpha_*dt_*sum;
  }

  // third term
  sum = 0.;
  for (Index it=1; it<Nt_; it++) {
    sum += a_y_[it-1]*x[y_index(Nx_,it)] + a_u_[it-1]*x[u_index(it)];
  }
  sum += 0.5*(a_y_[Nt_-1]*x[y_index(Nx_,Nt_)] + a_u_[Nt_-1]*x[u_index(Nt_)]);
  obj_value += dt_*sum;

  return true;
}

template <class T>
bool MittelmannParaCntrlBase<T>::
eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // First set all y entries to zero
  for (Index jx=0; jx<=Nx_; jx++) {
    for (Index it=0; it<=Nt_; it++) {
      grad_f[y_index(jx,it)] = 0.;
    }
  }

  // y entries from y target
  grad_f[y_index(0,Nt_)] = 0.5*dx_*(x[y_index(0,Nt_)] - y_T_[0]);
  for (Index jx=1; jx<Nx_; jx++) {
    grad_f[y_index(jx,Nt_)] = dx_*(x[y_index(jx,Nt_)] - y_T_[jx]);
  }
  grad_f[y_index(Nx_,Nt_)] = 0.5*dx_*(x[y_index(Nx_,Nt_)] - y_T_[Nx_]);

  // y entries from thrid term
  for (Index it=1; it<Nt_; it++) {
    grad_f[y_index(Nx_,it)] = dt_*a_y_[it-1];
  }
  grad_f[y_index(Nx_,Nt_)] += 0.5*dt_*a_y_[Nt_-1];

  // u entries from smoothing and third term
  for (Index it=1; it<Nt_; it++) {
    grad_f[u_index(it)] = alpha_*dt_*x[u_index(it)] + dt_*a_u_[it-1];
  }
  grad_f[u_index(Nt_)] = 0.5*dt_*(alpha_*x[u_index(Nt_)] + a_u_[Nt_-1]);

  return true;
}

template <class T>
bool MittelmannParaCntrlBase<T>::
eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  typename T::ProblemSpecs p;

  Index ig=0;

  Number f = 1./(2.*dx_*dx_);
  for (Index jx=1; jx<Nx_; jx++) {
    for (Index it=0; it<Nt_; it++) {
      g[ig] = (x[y_index(jx,it)]-x[y_index(jx,it+1)])/dt_
              + f*(x[y_index(jx-1,it)] - 2.*x[y_index(jx,it)]
                   + x[y_index(jx+1,it)] + x[y_index(jx-1,it+1)]
                   - 2.*x[y_index(jx,it+1)] + x[y_index(jx+1,it+1)]);
      ig++;
    }
  }

  for (Index it=1; it<=Nt_; it++) {
    g[ig] = (x[y_index(2,it)] - 4.*x[y_index(1,it)] + 3.*x[y_index(0,it)]);
    ig++;
  }

  f = 1./(2.*dx_);
  for (Index it=1; it<=Nt_; it++) {
    g[ig] =
      f*(x[y_index(Nx_-2,it)] - 4.*x[y_index(Nx_-1,it)]
         + 3.*x[y_index(Nx_,it)]) + beta_*x[y_index(Nx_,it)]
      - x[u_index(it)] + p.phi(x[y_index(Nx_,it)]);
    ig++;
  }

  DBG_ASSERT(ig == m);

  return true;
}

template <class T>
bool MittelmannParaCntrlBase<T>::
eval_jac_g(Index n, const Number* x, bool new_x,
           Index m, Index nele_jac, Index* iRow, Index *jCol,
           Number* values)
{
  typename T::ProblemSpecs p;

  Index ijac = 0;

  if (values == NULL) {
    Index ig = 0;
    for (Index jx=1; jx<Nx_; jx++) {
      for (Index it=0; it<Nt_; it++) {
        iRow[ijac] = ig;
        jCol[ijac] = y_index(jx-1,it);
        ijac++;
        iRow[ijac] = ig;
        jCol[ijac] = y_index(jx,it);
        ijac++;
        iRow[ijac] = ig;
        jCol[ijac] = y_index(jx+1,it);
        ijac++;
        iRow[ijac] = ig;
        jCol[ijac] = y_index(jx-1,it+1);
        ijac++;
        iRow[ijac] = ig;
        jCol[ijac] = y_index(jx,it+1);
        ijac++;
        iRow[ijac] = ig;
        jCol[ijac] = y_index(jx+1,it+1);
        ijac++;

        ig++;
      }
    }

    for (Index it=1; it<=Nt_; it++) {
      iRow[ijac] = ig;
      jCol[ijac] = y_index(0,it);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = y_index(1,it);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = y_index(2,it);
      ijac++;

      ig++;
    }

    for (Index it=1; it<=Nt_; it++) {
      iRow[ijac] = ig;
      jCol[ijac] = y_index(Nx_-2,it);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = y_index(Nx_-1,it);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = y_index(Nx_,it);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = u_index(it);
      ijac++;

      ig++;
    }
    DBG_ASSERT(ig == m);
  }
  else {
    Number f = 1./(2.*dx_*dx_);
    Number f2 = 1./dt_;
    for (Index jx=1; jx<Nx_; jx++) {
      for (Index it=0; it<Nt_; it++) {
        values[ijac++] = f;
        values[ijac++] = f2 - 2.*f;
        values[ijac++] = f;
        values[ijac++] = f;
        values[ijac++] = -f2 - 2.*f;
        values[ijac++] = f;
      }
    }

    for (Index it=1; it<=Nt_; it++) {
      values[ijac++] = 3.;
      values[ijac++] = -4.;
      values[ijac++] = 1.;
    }

    f = 1./(2.*dx_);
    for (Index it=1; it<=Nt_; it++) {
      values[ijac++] = f;
      values[ijac++] = -4.*f;
      values[ijac++] = 3.*f + beta_ + p.phi_dy(x[y_index(Nx_,it)]);
      values[ijac++] = -1.;
    }

  }

  return true;
}

template <class T>
bool MittelmannParaCntrlBase<T>::
eval_h(Index n, const Number* x, bool new_x,
       Number obj_factor, Index m, const Number* lambda,
       bool new_lambda, Index nele_hess, Index* iRow,
       Index* jCol, Number* values)
{
  typename T::ProblemSpecs p;

  Index ihes = 0;

  if (values == NULL) {
    // y values from objective
    for (Index jx=0; jx<= Nx_; jx++) {
      iRow[ihes] = y_index(jx,Nt_);
      jCol[ihes] = y_index(jx,Nt_);
      ihes++;
    }
    // u from objective
    for (Index it=1; it<=Nt_; it++) {
      iRow[ihes] = u_index(it);
      jCol[ihes] = u_index(it);
      ihes++;
    }

    // constraint
    if (!p.phi_dydy_always_zero()) {
      for (Index it=1; it<=Nt_; it++) {
        iRow[ihes] = y_index(Nx_,it);
        jCol[ihes] = y_index(Nx_,it);
        ihes++;
      }
    }
  }
  else {
    // y values from objective
    values[ihes++] = obj_factor*0.5*dx_;
    for (Index jx=1; jx<Nx_; jx++) {
      values[ihes++] = obj_factor*dx_;
    }
    values[ihes++] = obj_factor*0.5*dx_;
    // u from objective
    for (Index it=1; it<Nt_; it++) {
      values[ihes++] = obj_factor*alpha_*dt_;
    }
    values[ihes++] = obj_factor*0.5*alpha_*dt_;

    // constrainT
    if (!p.phi_dydy_always_zero()) {
      Index ig = (Nx_-1)*Nt_ + Nt_;
      for (Index it=1; it<=Nt_; it++) {
        values[ihes++] = lambda[ig++]*p.phi_dydy(x[y_index(Nx_,it)]);
      }
    }
  }

  DBG_ASSERT(ihes==nele_hess);

  return true;
}

template <class T>
void MittelmannParaCntrlBase<T>::
finalize_solution(SolverReturn status,
                  Index n, const Number* x, const Number* z_L,
                  const Number* z_U,
                  Index m, const Number* g, const Number* lambda,
                  Number obj_value,
		  const IpoptData* ip_data,
		  IpoptCalculatedQuantities* ip_cq)
{}

class MittelmannParaCntrl5_1
{
public:
  class ProblemSpecs
  {
  public:
    ProblemSpecs ()
        :
        pi_(4.*atan(1.)),
        exp13_(exp(1./3.)),
        exp23_(exp(2./3.)),
        exp1_(exp(1.)),
        expm1_(exp(-1.)),
        sqrt2_(sqrt(2.))
    {}
    Number T()
    {
      return 1.;
    }
    Number l()
    {
      return pi_/4.;
    }
    Number lb_y()
    {
      return -1e20;
    }
    Number ub_y()
    {
      return 1e20;
    }
    Number lb_u()
    {
      return 0.;
    }
    Number ub_u()
    {
      return 1.;
    }
    Number alpha()
    {
      return sqrt2_/2.*(exp23_-exp13_);
    }
    Number beta()
    {
      return 1.;
    }
    inline Number y_T(Number x)
    {
      return (exp1_ + expm1_)*cos(x);
    }
    inline Number a(Number x)
    {
      return cos(x);
    }
    inline Number a_y(Number t)
    {
      return -exp(-2.*t);
    }
    inline Number a_u(Number t)
    {
      return sqrt2_/2.*exp13_;
    }
    inline Number b(Number t)
    {
      return exp(-4.*t)/4.
             - Min(1.,Max(0.,(exp(t)-exp13_)/(exp23_-exp13_)));
    }
    inline Number phi(Number y)
    {
      return y*pow(fabs(y),3);
    }
    inline Number phi_dy(Number y)
    {
      return 4.*pow(fabs(y),3);
    }
    inline Number phi_dydy(Number y)
    {
      return 12.*y*y;
    }
    inline bool phi_dydy_always_zero()
    {
      return false;
    }
  private:
    const Number pi_;
    const Number exp13_;
    const Number exp23_;
    const Number exp1_;
    const Number expm1_;
    const Number sqrt2_;
  };
};

class MittelmannParaCntrl5_2_1
{
public:
  class ProblemSpecs
  {
  public:
    ProblemSpecs ()
    {}
    Number T()
    {
      return 1.58;
    }
    Number l()
    {
      return 1.;
    }
    Number lb_y()
    {
      return -1e20;
    }
    Number ub_y()
    {
      return 1e20;
    }
    Number lb_u()
    {
      return -1.;
    }
    Number ub_u()
    {
      return 1.;
    }
    Number alpha()
    {
      return 0.001;
    }
    Number beta()
    {
      return 1.;
    }
    inline Number y_T(Number x)
    {
      return .5*(1.-x*x);
    }
    inline Number a(Number x)
    {
      return 0.;
    }
    inline Number a_y(Number t)
    {
      return 0.;
    }
    inline Number a_u(Number t)
    {
      return 0.;
    }
    inline Number b(Number t)
    {
      return 0.;
    }
    inline Number phi(Number y)
    {
      return 0.;
    }
    inline Number phi_dy(Number y)
    {
      return 0.;
    }
    inline Number phi_dydy(Number y)
    {
      DBG_ASSERT(false);
      return 0.;
    }
    inline bool phi_dydy_always_zero()
    {
      return true;
    }
  };
};

class MittelmannParaCntrl5_2_2
{
public:
  class ProblemSpecs
  {
  public:
    ProblemSpecs ()
    {}
    Number T()
    {
      return 1.58;
    }
    Number l()
    {
      return 1.;
    }
    Number lb_y()
    {
      return -1e20;
    }
    Number ub_y()
    {
      return 1e20;
    }
    Number lb_u()
    {
      return -1.;
    }
    Number ub_u()
    {
      return 1.;
    }
    Number alpha()
    {
      return 0.001;
    }
    Number beta()
    {
      return 0.;
    }
    inline Number y_T(Number x)
    {
      return .5*(1.-x*x);
    }
    inline Number a(Number x)
    {
      return 0.;
    }
    inline Number a_y(Number t)
    {
      return 0.;
    }
    inline Number a_u(Number t)
    {
      return 0.;
    }
    inline Number b(Number t)
    {
      return 0.;
    }
    inline Number phi(Number y)
    {
      return y*y;
    }
    inline Number phi_dy(Number y)
    {
      return 2.*y;
    }
    inline Number phi_dydy(Number y)
    {
      return 2.;
    }
    inline bool phi_dydy_always_zero()
    {
      return false;
    }
  };
};

class MittelmannParaCntrl5_2_3
{
public:
  class ProblemSpecs
  {
  public:
    ProblemSpecs ()
    {}
    Number T()
    {
      return 1.58;
    }
    Number l()
    {
      return 1.;
    }
    Number lb_y()
    {
      return 0.;
    }
    Number ub_y()
    {
      return 0.675;
    }
    Number lb_u()
    {
      return -1.;
    }
    Number ub_u()
    {
      return 1.;
    }
    Number alpha()
    {
      return 0.001;
    }
    Number beta()
    {
      return 0.;
    }
    inline Number y_T(Number x)
    {
      return .5*(1.-x*x);
    }
    inline Number a(Number x)
    {
      return 0.;
    }
    inline Number a_y(Number t)
    {
      return 0.;
    }
    inline Number a_u(Number t)
    {
      return 0.;
    }
    inline Number b(Number t)
    {
      return 0.;
    }
    inline Number phi(Number y)
    {
      return y*y;
    }
    inline Number phi_dy(Number y)
    {
      return 2.*y;
    }
    inline Number phi_dydy(Number y)
    {
      return 2.;
    }
    inline bool phi_dydy_always_zero()
    {
      return false;
    }
  };
};

class MittelmannParaCntrl5_try
{
public:
  class ProblemSpecs
  {
  public:
    ProblemSpecs ()
        :
        pi_(4.*atan(1.)),
        exp13_(exp(1./3.)),
        exp23_(exp(2./3.)),
        exp1_(exp(1.)),
        expm1_(exp(-1.)),
        sqrt2_(sqrt(2.))
    {}
    Number T()
    {
      return 1.;
    }
    Number l()
    {
      return pi_/4.;
    }
    Number lb_y()
    {
      return -1e20;
    }
    Number ub_y()
    {
      return 1e20;
    }
    Number lb_u()
    {
      return 0.;
    }
    Number ub_u()
    {
      return 1.;
    }
    Number alpha()
    {
      return sqrt2_/2.*(exp23_-exp13_);
    }
    Number beta()
    {
      return 1.;
    }
    inline Number y_T(Number x)
    {
      return (exp1_ + expm1_)*cos(x);
    }
    inline Number a(Number x)
    {
      return cos(x);
    }
    inline Number a_y(Number t)
    {
      return -exp(-2.*t);
    }
    inline Number a_u(Number t)
    {
      return sqrt2_/2.*exp13_;
    }
    inline Number b(Number t)
    {
      return exp(-4.*t)/4.
             - Min(1.,Max(0.,(exp(t)-exp13_)/(exp23_-exp13_)));
    }
    inline Number phi(Number y)
    {
      return -y*sin(y/10.);
    }
    inline Number phi_dy(Number y)
    {
      return -y*cos(y/10.)/10. - sin(y/10.);
    }
    inline Number phi_dydy(Number y)
    {
      return y*sin(y/10.)/100.;
    }
    inline bool phi_dydy_always_zero()
    {
      return false;
    }
  private:
    const Number pi_;
    const Number exp13_;
    const Number exp23_;
    const Number exp1_;
    const Number expm1_;
    const Number sqrt2_;
  };
};

#endif
