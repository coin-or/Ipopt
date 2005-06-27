// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPNLP_HPP__
#define __IPNLP_HPP__

#include "IpUtils.hpp"
#include "IpVector.hpp"
#include "IpSmartPtr.hpp"
#include "IpMatrix.hpp"
#include "IpSymMatrix.hpp"
#include "IpOptionsList.hpp"

namespace Ipopt
{

  /** Brief Class Description.
   *  Detailed Class Description.
   */
  class NLP : public ReferencedObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default constructor */
    NLP()
    {}

    /** Default destructor */
    virtual ~NLP()
    {}
    //@}

    /** Exceptions */
    //@{
    DECLARE_STD_EXCEPTION(USER_SCALING_NOT_IMPLEMENTED);
    //@}

    /** @name NLP Initialization (overload in
     *  derived classes).*/
    //@{
    /** Overload if you want the chance to process options or parameters that
     *  may be specific to the NLP */
    virtual bool ProcessOptions(const OptionsList& options,
                                const std::string& prefix)
    {
      return true;
    }

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
                           SmartPtr<SymMatrixSpace>& Hess_lagrangian_space)=0;

    /** Method for obtaining the bounds information */
    virtual bool GetBoundsInformation(Matrix& Px_L,
                                      Vector& x_L,
                                      Matrix& Px_U,
                                      Vector& x_U,
                                      Matrix& Pd_L,
                                      Vector& d_L,
                                      Matrix& Pd_U,
                                      Vector& d_U)=0;

    /** Method for obtaining the starting point for all the
     *  iterates. ToDo it might not make sense to ask for initial
     *  values for v_L and v_U? */
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
    )=0;
    //@}

    /** @name NLP evaluation routines (overload
     *  in derived classes. */
    //@{
    virtual bool Eval_f(const Vector& x, Number& f) = 0;

    virtual bool Eval_grad_f(const Vector& x, Vector& g_f) = 0;

    virtual bool Eval_c(const Vector& x, Vector& c) = 0;

    virtual bool Eval_jac_c(const Vector& x, Matrix& jac_c) = 0;

    virtual bool Eval_d(const Vector& x, Vector& d) = 0;

    virtual bool Eval_jac_d(const Vector& x, Matrix& jac_d) = 0;

    virtual bool Eval_h(const Vector& x,
                        Number obj_factor,
                        const Vector& yc,
                        const Vector& yd,
                        SymMatrix& h) = 0;
    //@}

    /** @name NLP solution routines. (Overload in derived classes.) */
    //@{
    virtual void FinalizeSolution(ApplicationReturnStatus status,
                                  const Vector& x, const Vector& z_L, const Vector& z_U,
                                  const Vector& c, const Vector& d,
                                  const Vector& y_c, const Vector& y_d,
                                  Number obj_value)
    {}
    //@}

    /** Routines to get the scaling parameters. These do not need to be overloaded
     *  unless the options are set for User scaling
     */
    //@{
    virtual void GetScalingParameters(Number& obj_scaling, Vector& x_scaling,
                                      Vector& c_scaling, Vector& d_scaling) const
    {
      THROW_EXCEPTION(USER_SCALING_NOT_IMPLEMENTED,
                      "You have set options for user provided scaling, but have"
                      " not implemented GetScalingParameters in the NLP interface");
    }
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
    NLP(const NLP&);

    /** Overloaded Equals Operator */
    void operator=(const NLP&);
    //@}
  };

} // namespace Ipopt

#endif
