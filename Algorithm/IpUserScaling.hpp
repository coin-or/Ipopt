// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpNLPScaling.hpp,v 1.3 2004/11/19 00:55:18 andreasw Exp $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPUSERSCALING_HPP__
#define __IPUSERSCALING_HPP__

#include "IpNLPScaling.hpp"
#include "IpNLP.hpp"
#include "IpScaledMatrix.hpp"
#include "IpSymScaledMatrix.hpp"

namespace Ipopt
{
  /** This class does problem scaling by getting scaling parameters
   *  from the user (through the NLP interface).
   */
  class UserScaling : public NLPScalingObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    UserScaling(const SmartPtr<const NLP>& nlp)
        : nlp_(nlp)
    {}

    /** Default destructor */
    virtual ~UserScaling()
    {}
    //@}

    /** Methods to map scaled and unscaled matrices */
    //@{
    /** Returns an obj-scaled version of the given scalar */
    virtual Number apply_obj_scaling(const Number& f);
    /** Returns an obj-unscaled version of the given scalar */
    virtual Number unapply_obj_scaling(const Number& f);
    /** Returns an x-scaled version of the given vector */
    virtual SmartPtr<Vector> apply_vector_scaling_x_NonConst(const SmartPtr<const Vector>& v);
    /** Returns an x-scaled version of the given vector */
    virtual SmartPtr<const Vector> apply_vector_scaling_x(const SmartPtr<const Vector>& v);
    /** Returns an x-unscaled version of the given vector */
    virtual SmartPtr<Vector> unapply_vector_scaling_x_NonConst(const SmartPtr<const Vector>& v);
    /** Returns an x-unscaled version of the given vector */
    virtual SmartPtr<const Vector> unapply_vector_scaling_x(const SmartPtr<const Vector>& v);
    /** Returns an c-scaled version of the given vector */
    virtual SmartPtr<const Vector> apply_vector_scaling_c(const SmartPtr<const Vector>& v);
    /** Returns an c-unscaled version of the given vector */
    virtual SmartPtr<const Vector> unapply_vector_scaling_c(const SmartPtr<const Vector>& v);
    /** Returns an c-scaled version of the given vector */
    virtual SmartPtr<Vector> apply_vector_scaling_c_NonConst(const SmartPtr<const Vector>& v);
    /** Returns an c-unscaled version of the given vector */
    virtual SmartPtr<Vector> unapply_vector_scaling_c_NonConst(const SmartPtr<const Vector>& v);
    /** Returns an d-scaled version of the given vector */
    virtual SmartPtr<const Vector> apply_vector_scaling_d(const SmartPtr<const Vector>& v);
    /** Returns an d-unscaled version of the given vector */
    virtual SmartPtr<const Vector> unapply_vector_scaling_d(const SmartPtr<const Vector>& v);
    /** Returns an d-scaled version of the given vector */
    virtual SmartPtr<Vector> apply_vector_scaling_d_NonConst(const SmartPtr<const Vector>& v);
    /** Returns an d-unscaled version of the given vector */
    virtual SmartPtr<Vector> unapply_vector_scaling_d_NonConst(const SmartPtr<const Vector>& v);
    /** Returns a scaled version of the jacobian for c.
     *  If the overloaded method does not make a new matrix, make sure to set the matrix
     *  ptr passed in to NULL.
     */
    virtual SmartPtr<const Matrix> apply_jac_c_scaling(SmartPtr<const Matrix> matrix);
    /** Returns a scaled version of the jacobian for d
     *  If the overloaded method does not create a new matrix, make sure to set the matrix
     *  ptr passed in to NULL.
     */
    virtual SmartPtr<const Matrix> apply_jac_d_scaling(SmartPtr<const Matrix> matrix);
    /** Returns a scaled version of the hessian of the lagrangian
     *  If the overloaded method does not create a new matrix, make sure to set the matrix
     *  ptr passed in to NULL.
     */
    virtual SmartPtr<const SymMatrix> apply_hessian_scaling(SmartPtr<const SymMatrix> matrix);
    //@}

    /** This method is called by the IpoptNLP's at a convenient time to
     *  compute and/or read scaling factors 
     */
    virtual void DetermineScaling(const SmartPtr<const VectorSpace> x_space,
                                  const SmartPtr<const VectorSpace> c_space,
                                  const SmartPtr<const VectorSpace> d_space,
                                  const SmartPtr<const MatrixSpace> jac_c_space,
                                  const SmartPtr<const MatrixSpace> jac_d_space,
                                  const SmartPtr<const SymMatrixSpace> h_space,
                                  SmartPtr<const MatrixSpace>& new_jac_c_space,
                                  SmartPtr<const MatrixSpace>& new_jac_d_space,
                                  SmartPtr<const SymMatrixSpace>& new_h_space);

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
    UserScaling(const UserScaling&);

    /** Overloaded Equals Operator */
    void operator=(const UserScaling&);
    //@}

    /** pointer to the NLP to get scaling parameters */
    SmartPtr<const NLP> nlp_;

    /** Scaling parameters - we only need to keep copies of
     *  the objective scaling and the x scaling - the others we can
     *  get from the scaled matrix spaces.
     */
    //@{
    /** objective scaling parameter */
    Number df_;
    /** x scaling */
    SmartPtr<Vector> dx_;
    //@}

    /** Scaled Matrix Spaces */
    //@{
    /** Scaled jacobian of c space */
    SmartPtr<ScaledMatrixSpace> scaled_jac_c_space_;
    /** Scaled jacobian of d space */
    SmartPtr<ScaledMatrixSpace> scaled_jac_d_space_;
    /** Scaled hessian of lagrangian spacea */
    SmartPtr<SymScaledMatrixSpace> scaled_h_space_;
    //@}
  };
} // namespace Ipopt
#endif
