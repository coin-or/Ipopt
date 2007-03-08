// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPCOMPOSITENLP_HPP__
#define __IPCOMPOSITENLP_HPP__

#include "IpNLP.hpp"

namespace Ipopt
{

  /* forward declaration */
  class SymLinearSolver;

  /** This class creates a composite NLP from a list of NLP's.
   *  This is a Composite class (Design Patterns) that creates a single NLP from
   *  a list of NLP's with some common variables. This allows users to "link" 
   *  different NLP's into a single NLP easily
   */
  class CompositeNLP : public NLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default constructor */
    CompositeNLP(std::vector<SmartPtr<NLP> > nlps, SmartPtr<VectorSpace> q_space,
                 std::vector<SmartPtr<VectorSpace> > linking_eqn_c_spaces,
                 std::vector<SmartPtr<Matrix> > Jx_linking_eqns,
                 std::vector<SmartPtr<Matrix> > Jq_linking_eqns);

    /** Default destructor */
    virtual ~CompositeNLP();
    //@}

    /**@name Exceptions */
    //@{
    DECLARE_STD_EXCEPTION(INVALID_JACOBIAN_DIMENSION_FOR_LINKING_EQUATIONS);
    //@}

    /** @name CompositeNLP Initialization. */
    //@{
    virtual bool ProcessOptions(const OptionsList& options,
                                const std::string& prefix);

    /** Method for creating the derived vector / matrix types
     *  (Do not delete these, the ). */
    virtual bool GetSpaces(SmartPtr<const VectorSpace>& x_space,
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

    /** Method for obtaining the bounds information */
    virtual bool GetBoundsInformation(const Matrix& Px_L,
                                      Vector& x_L,
                                      const Matrix& Px_U,
                                      Vector& x_U,
                                      const Matrix& Pd_L,
                                      Vector& d_L,
                                      const Matrix& Pd_U,
                                      Vector& d_U);

    /** Method for obtaining the starting point
     *  for all the iterates. */
    virtual bool GetStartingPoint(
      SmartPtr<Vector> x,
      bool need_x,
      SmartPtr<Vector> y_c,
      bool need_y_c,
      SmartPtr<Vector> y_d,
      bool need_y_d,
      SmartPtr<Vector> z_L,
      bool need_z_L,
      SmartPtr<Vector> z_U,
      bool need_z_U
    );
    //@}

    /** @name CompositeNLP evaluation routines. */
    //@{
    virtual bool Eval_f(const Vector& x, Number& f);

    virtual bool Eval_grad_f(const Vector& x, Vector& g_f);

    virtual bool Eval_c(const Vector& x, Vector& c);

    virtual bool Eval_jac_c(const Vector& x, Matrix& jac_c);

    virtual bool Eval_d(const Vector& x, Vector& d);

    virtual bool Eval_jac_d(const Vector& x, Matrix& jac_d);

    virtual bool Eval_h(const Vector& x,
                        Number obj_factor,
                        const Vector& yc,
                        const Vector& yd,
                        SymMatrix& h);
    //@}

    /** @name Information about the Composite Structure */
    //@{
    /** returns an appropriate linear solver for the problem structure */
    SmartPtr<SymLinearSolver> CreateLinearSolver();
    //@}
  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    CompositeNLP();

    /** Copy Constructor */
    CompositeNLP(const CompositeNLP&);

    /** Overloaded Equals Operator */
    void operator=(const CompositeNLP&);
    //@}

    /** Journalist */
    SmartPtr<const Journalist> jnlst_;

    /** std::vector of nlps */
    std::vector<SmartPtr<NLP> > nlps_;

    /** vector space for the linking variables, q */
    SmartPtr<VectorSpace> q_space_;

    /** std::vector of VectorSpaces for the linking equations */
    std::vector<SmartPtr<VectorSpace> > linking_eqn_c_spaces_;

    /**@name Data about Linking Equations. For now,
     *  this class only allows linear linking equations. Therefore,
     *  the jacobian of the linking equations is set only once at
     *  the beginning. */
    //@{
    /** std::vector of Jacobian of the linking eqns (one for each nlp)
     *   with respect to x variables */
    std::vector<SmartPtr<Matrix> > Jx_linking_eqns_;

    /** std::vector of Jacobian of the linking eqns (one for each nlp)
     *   with respect to linking variables (q) */
    std::vector<SmartPtr<Matrix> > Jq_linking_eqns_;
    //@}
  };

} // namespace Ipopt

#endif
