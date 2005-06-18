// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPTNLPADAPTER_HPP__
#define __IPTNLPADAPTER_HPP__

#include "IpNLP.hpp"
#include "IpTNLP.hpp"
#include "IpDenseVector.hpp"
#include "IpExpansionMatrix.hpp"
#include "IpGenTMatrix.hpp"
#include "IpSymTMatrix.hpp"
#include "IpIpoptType.hpp"
namespace Ipopt
{

  DeclareIpoptType(TNLPAdapter);

  /** This class Adapts the TNLP interface so it looks like an NLP interface.
   *  This is an Adapter class (Design Patterns) that converts  a TNLP to an
   *  NLP. This allows users to write to the "more convenient" TNLP interface.
   */
  class TNLPAdapter : public NLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default constructor */
    TNLPAdapter(const SmartPtr<TNLP> tnlp);

    /** Default destructor */
    virtual ~TNLPAdapter();
    //@}

    /**@name Exceptions */
    //@{
    DECLARE_STD_EXCEPTION(INVALID_TNLP);
    //@}

    /** @name TNLPAdapter Initialization. */
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

    /** @name TNLPAdapter evaluation routines. */
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

    /** @name Solution Reporting Methods */
    //@{
    virtual void FinalizeSolution(ApplicationReturnStatus status,
                                  const Vector& x, const Vector& z_L, const Vector& z_U,
                                  const Vector& c, const Vector& d,
                                  const Vector& y_c, const Vector& y_d,
                                  Number obj_value);
    //@}

    /** @name Methods for translating data for IpoptNLP into the TNLP
     *  data.  These methods can be used to obtain the current (or
     *  final) data for the TNLP formulation from the IpoptNLP
     *  structure. */
    //@{
    /** Sort the primal variables, and add the fixed values in x */
    void ResortX(const Vector& x, Number* x_orig);
    void ResortG(const Vector& c, const Vector& d, Number *g_orig);
    void ResortBnds(const Vector& x_L, Number* x_L_orig,
                    const Vector& x_U, Number* x_U_orig);
    //@}

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
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
    /** Copy Constructor */
    TNLPAdapter(const TNLPAdapter&);

    /** Overloaded Equals Operator */
    void operator=(const TNLPAdapter&);
    //@}

    /** Journalist */
    SmartPtr<const Journalist> jnlst_;

    /** Pointer to the TNLP class (class specific to Number* vectors and
     *  harwell triplet matrices) */
    SmartPtr<TNLP> tnlp_;

    /**@name Algorithmic parameters */
    //@{
    /** Value for a lower bound that denotes -infinity */
    Number nlp_lower_bound_inf_;
    /** Value for a upper bound that denotes infinity */
    Number nlp_upper_bound_inf_;
    /** Maximal slack for one-sidedly bounded variables.  If a
     *  variable has only one bound, say a lower bound xL, then an
     *  upper bound xL + max_onesided_bound_slack_.  If this value is
     *  zero, no upper bound is added. */
    Number max_onesided_bound_slack_;
    //@}

    /**@name Problem Size Data */
    //@{
    Index n_full_x_; /** full dimension of x (fixed + non-fixed) */
    Index n_full_g_; /** full dimension of g (c + d) */
    Index nz_jac_c_; /** non-zeros of the jacobian of c */
    Index nz_jac_d_; /** non-zeros of the jacobian of d */
    Index nz_full_jac_g_; /** number of non-zeros in full-size Jacobian of g */
    Index nz_full_h_; /** number of non-zeros in full-size Hessian */
    Index nz_h_;     /** number of non-zeros in the non-fixed-size Hessian */
    //@}

    /**@name Local Copy of the Data */
    //@{
    Number* full_x_; /** copy of the full x vector (fixed & non-fixed) */
    Number* full_lambda_; /** copy of lambda (yc & yd) */
    Number* full_g_; /** copy of g (c & d) */
    Number* jac_g_; /** the values for the full jacobian of g */
    Number* c_rhs_; /** the rhs values of c */
    //@}

    /**@name Tags for deciding when to update internal copies of vectors */
    //@{
    TaggedObject::Tag x_tag_for_iterates_;
    TaggedObject::Tag y_c_tag_for_iterates_;
    TaggedObject::Tag y_d_tag_for_iterates_;
    TaggedObject::Tag x_tag_for_g_;
    TaggedObject::Tag x_tag_for_jac_g_;
    //@}

    /**@name Methods to update the values in the local copies of vectors */
    //@{
    bool update_local_x(const Vector& x);
    bool update_local_lambda(const Vector& y_c, const Vector& y_d);
    //@}

    /**@name Internal routines for evaluating g and jac_g (values stored since
     * they are used in both c and d routins */
    //@{
    bool internal_eval_g(bool new_x);
    bool internal_eval_jac_g(bool new_x);
    //@}

    /**@name Internal Permutation Spaces and matrices
     */
    //@{
    /** Expansion from fixed x (ipopt) to full x */
    SmartPtr<ExpansionMatrix> P_x_full_x_;
    SmartPtr<ExpansionMatrixSpace> P_x_full_x_space_;

    /** Expansion from fixed x_L (ipopt) to full x */
    SmartPtr<ExpansionMatrix> P_x_x_L_;
    SmartPtr<ExpansionMatrixSpace> P_x_x_L_space_;

    /** Expansion from fixed x_U (ipopt) to full x */
    SmartPtr<ExpansionMatrix> P_x_x_U_;
    SmartPtr<ExpansionMatrixSpace> P_x_x_U_space_;

    /** Expansion from c only (ipopt) to full ampl c */
    SmartPtr<ExpansionMatrixSpace> P_c_g_space_;
    SmartPtr<ExpansionMatrix> P_c_g_;

    /** Expansion from d only (ipopt) to full ampl d */
    SmartPtr<ExpansionMatrixSpace> P_d_g_space_;
    SmartPtr<ExpansionMatrix> P_d_g_;

    Index* jac_idx_map_;
    Index* h_idx_map_;
    //@}
  };

} // namespace Ipopt

#endif
