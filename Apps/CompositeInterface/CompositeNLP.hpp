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
    virtual bool GetSpaces(SmartPtr<VectorSpace>& x_space,
                           SmartPtr<VectorSpace>& c_space,
                           SmartPtr<VectorSpace>& d_space,
                           SmartPtr<VectorSpace>& x_l_space,
                           SmartPtr<MatrixSpace>& px_l_space,
                           SmartPtr<VectorSpace>& x_u_space,
                           SmartPtr<MatrixSpace>& px_u_space,
                           SmartPtr<VectorSpace>& d_l_space,
                           SmartPtr<MatrixSpace>& pd_l_space,
                           SmartPtr<VectorSpace>& d_u_space,
                           SmartPtr<MatrixSpace>& pd_u_space,
                           SmartPtr<MatrixSpace>& Jac_c_space,
                           SmartPtr<MatrixSpace>& Jac_d_space,
                           SmartPtr<SymMatrixSpace>& Hess_lagrangian_space);

    /** Method for obtaining the bounds information */
    virtual bool GetBoundsInformation(Matrix& Px_L,
                                      Vector& x_L,
                                      Matrix& Px_U,
                                      Vector& x_U,
                                      Matrix& Pd_L,
                                      Vector& d_L,
                                      Matrix& Pd_U,
                                      Vector& d_U);

    /** Method for obtaining the starting point
     *  for all the iterates. */
    virtual bool GetStartingPoint(
      Vector& x,
      bool need_x,
      Vector& y_c,
      bool need_y_c,
      Vector& y_d,
      bool need_y_d,
      Vector& z_L,
      bool need_z_L,
      Vector& z_U,
      bool need_z_U,
      Vector& v_L,
      bool need_v_L,
      Vector& v_U,
      bool need_v_U
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
