// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPORIGIPOPTNLP_HPP__
#define __IPORIGIPOPTNLP_HPP__

#include "IpIpoptNLP.hpp"
#include "IpException.hpp"
#include"IpIpoptType.hpp"

namespace Ipopt
{

  DeclareIpoptType(OrigIpoptNLP);

  /** This class maps the traditional NLP into
   *  something that is more useful by Ipopt.
   *  This class takes care of storing the
   *  calculated model results, handles cacheing,
   *  and (some day) takes care of addition of slacks.
   */
  class OrigIpoptNLP : public IpoptNLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    OrigIpoptNLP(const SmartPtr<const Journalist>& jnlst,
                 const SmartPtr<NLP>& nlp,
                 const SmartPtr<NLPScalingObject>& nlp_scaling);

    /** Default destructor */
    virtual ~OrigIpoptNLP();
    //@}

    /** Initialize - overloaded from IpoptNLP */
    virtual bool Initialize(const Journalist& jnlst,
                            const OptionsList& options,
                            const std::string& prefix);

    /** Initialize (create) structures for
     *  the iteration data */
    virtual bool InitializeStructures(SmartPtr<Vector>& x,
                                      bool init_x,
                                      SmartPtr<Vector>& y_c,
                                      bool init_y_c,
                                      SmartPtr<Vector>& y_d,
                                      bool init_y_d,
                                      SmartPtr<Vector>& z_L,
                                      bool init_z_L,
                                      SmartPtr<Vector>& z_U,
                                      bool init_z_U,
                                      SmartPtr<Vector>& v_L,
                                      bool init_v_L,
                                      SmartPtr<Vector>& v_U,
                                      bool init_v_U
                                     );

    /** Accessor methods for model data */
    //@{
    /** Objective value */
    virtual Number f(const Vector& x);

    /** Gradient of the objective */
    virtual SmartPtr<const Vector> grad_f(const Vector& x);

    /** Equality constraint residual */
    virtual SmartPtr<const Vector> c(const Vector& x);

    /** Jacobian Matrix for equality constraints
     *  (current iteration) */
    virtual SmartPtr<const Matrix> jac_c(const Vector& x);

    /** Inequality constraint residual (reformulated
     *  as equalities with slacks */
    virtual SmartPtr<const Vector> d(const Vector& x);

    /** Jacobian Matrix for inequality constraints
     *  (current iteration) */
    virtual SmartPtr<const Matrix> jac_d(const Vector& x);

    /** Hessian of the lagrangian
     *  (current iteration) */
    virtual SmartPtr<const SymMatrix> h(const Vector& x,
                                        Number obj_factor,
                                        const Vector& yc,
                                        const Vector& yd
                                       );

    /** Lower bounds on x */
    virtual SmartPtr<const Vector> x_L()
    {
      return ConstPtr(x_L_);
    }

    /** Permutation matrix (x_L_ -> x) */
    virtual SmartPtr<const Matrix> Px_L()
    {
      return ConstPtr(Px_L_);
    }

    /** Upper bounds on x */
    virtual SmartPtr<const Vector> x_U()
    {
      return ConstPtr(x_U_);
    }

    /** Permutation matrix (x_U_ -> x */
    virtual SmartPtr<const Matrix> Px_U()
    {
      return ConstPtr(Px_U_);
    }

    /** Lower bounds on d */
    virtual SmartPtr<const Vector> d_L()
    {
      return ConstPtr(d_L_);
    }

    /** Permutation matrix (d_L_ -> d) */
    virtual SmartPtr<const Matrix> Pd_L()
    {
      return ConstPtr(Pd_L_);
    }

    /** Upper bounds on d */
    virtual SmartPtr<const Vector> d_U()
    {
      return ConstPtr(d_U_);
    }

    /** Permutation matrix (d_U_ -> d */
    virtual SmartPtr<const Matrix> Pd_U()
    {
      return ConstPtr(Pd_U_);
    }
    //@}

    /** Accessor method for vector/matrix spaces pointers */
    virtual void GetSpaces(SmartPtr<const VectorSpace>& x_space,
                           SmartPtr<const VectorSpace>& c_space,
                           SmartPtr<const VectorSpace>& d_space,
                           SmartPtr<const VectorSpace>& x_l_space,
                           SmartPtr<const MatrixSpace>& px_l_space,
                           SmartPtr<const VectorSpace>& x_u_space,
                           SmartPtr<const MatrixSpace>& px_u_space,
                           SmartPtr<const VectorSpace>& d_l_space,
                           SmartPtr<const MatrixSpace>& pd_l_space,
                           SmartPtr<const VectorSpace>& d_u_space,
                           SmartPtr<const MatrixSpace>& pd_u_space,
                           SmartPtr<const MatrixSpace>& Jac_c_space,
                           SmartPtr<const MatrixSpace>& Jac_d_space,
                           SmartPtr<const SymMatrixSpace>& Hess_lagrangian_space);

    /** Method for adapting the variable bounds.  This is called if
     *  slacks are becoming too small */
    virtual void AdjustVariableBounds(const Vector& new_x_L,
                                      const Vector& new_x_U,
                                      const Vector& new_d_L,
                                      const Vector& new_d_U);

    /** @name Counters for the number of function evaluations. */
    //@{
    virtual Index f_evals() const
    {
      return f_evals_;
    }
    virtual Index grad_f_evals() const
    {
      return grad_f_evals_;
    }
    virtual Index c_evals() const
    {
      return c_evals_;
    }
    virtual Index jac_c_evals() const
    {
      return jac_c_evals_;
    }
    virtual Index d_evals() const
    {
      return d_evals_;
    }
    virtual Index jac_d_evals() const
    {
      return jac_d_evals_;
    }
    virtual Index h_evals() const
    {
      return h_evals_;
    }
    //@}

    /** Methods for IpoptType */
    //@{
    /** Called by IpoptType to register the options */
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

  private:
    /** Pointer to the NLP */
    SmartPtr<NLP> nlp_;

    /** journalist */
    SmartPtr<const Journalist> jnlst_;

    /** Necessary Vector/Matrix spaces */
    //@{
    SmartPtr<const VectorSpace> x_space_;
    SmartPtr<const VectorSpace> c_space_;
    SmartPtr<const VectorSpace> d_space_;
    SmartPtr<const VectorSpace> x_l_space_;
    SmartPtr<const MatrixSpace> px_l_space_;
    SmartPtr<const VectorSpace> x_u_space_;
    SmartPtr<const MatrixSpace> px_u_space_;
    SmartPtr<const VectorSpace> d_l_space_;
    SmartPtr<const MatrixSpace> pd_l_space_;
    SmartPtr<const VectorSpace> d_u_space_;
    SmartPtr<const MatrixSpace> pd_u_space_;
    SmartPtr<const MatrixSpace> jac_c_space_;
    SmartPtr<const MatrixSpace> jac_d_space_;
    SmartPtr<const SymMatrixSpace> h_space_;

    SmartPtr<const MatrixSpace> scaled_jac_c_space_;
    SmartPtr<const MatrixSpace> scaled_jac_d_space_;
    SmartPtr<const SymMatrixSpace> scaled_h_space_;
    //@}
    /**@name Storage for Model Quantities */
    //@{
    /** Objective function */
    CachedResults<Number> f_cache_;

    /** Gradient of the objective function */
    CachedResults<SmartPtr<const Vector> > grad_f_cache_;

    /** Equality constraint residuals */
    CachedResults<SmartPtr<const Vector> > c_cache_;

    /** Jacobian Matrix for equality constraints
     *  (current iteration) */
    CachedResults<SmartPtr<const Matrix> > jac_c_cache_;

    /** Inequality constraint residual (reformulated
     *  as equalities with slacks */
    CachedResults<SmartPtr<const Vector> > d_cache_;

    /** Jacobian Matrix for inequality constraints
     *  (current iteration) */
    CachedResults<SmartPtr<const Matrix> > jac_d_cache_;

    /** Hessian of the lagrangian
     *  (current iteration) */
    CachedResults<SmartPtr<const SymMatrix> > h_cache_;

    /** Lower bounds on x */
    SmartPtr<Vector> x_L_;

    /** Permutation matrix (x_L_ -> x) */
    SmartPtr<Matrix> Px_L_;

    /** Upper bounds on x */
    SmartPtr<Vector> x_U_;

    /** Permutation matrix (x_U_ -> x */
    SmartPtr<Matrix> Px_U_;

    /** Lower bounds on d */
    SmartPtr<Vector> d_L_;

    /** Permutation matrix (d_L_ -> d) */
    SmartPtr<Matrix> Pd_L_;

    /** Upper bounds on d */
    SmartPtr<Vector> d_U_;

    /** Permutation matrix (d_U_ -> d */
    SmartPtr<Matrix> Pd_U_;
    //@}

    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    OrigIpoptNLP();

    /** Copy Constructor */
    OrigIpoptNLP(const OrigIpoptNLP&);

    /** Overloaded Equals Operator */
    void operator=(const OrigIpoptNLP&);
    //@}

    /** @name auxilliary functions */
    //@{
    /** relax the bounds by a relative move of relax_bound_factor.
     *  Here, relax_bound_factor should be negative (or zero) for
     *  lower bounds, and positive (or zero) for upper bounds.
     */
    void relax_bounds(Number bound_relax_factor, Vector& bounds);
    //@}

    /** @name Algorithmic parameters */
    //@{
    /** relaxation factor for the bounds */
    Number bound_relax_factor_;
    //@}

    /** Flag indicating if initialization method has been called */
    bool initialized_;

    /** @name Counters for the function evaluations */
    //@{
    Index f_evals_;
    Index grad_f_evals_;
    Index c_evals_;
    Index jac_c_evals_;
    Index d_evals_;
    Index jac_d_evals_;
    Index h_evals_;
    //@}
  };

} // namespace Ipopt

#endif
