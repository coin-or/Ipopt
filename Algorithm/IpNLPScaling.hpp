// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpNLPScaling.hpp,v 1.3 2004/11/19 00:55:18 andreasw Exp $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPNLPSCALING_HPP__
#define __IPNLPSCALING_HPP__

#include "IpSmartPtr.hpp"
#include "IpVector.hpp"
#include "IpMatrix.hpp"
#include "IpSymMatrix.hpp"

namespace Ipopt
{
  /** This is the abstract base class for problem scaling.
   *  It is repsonsible for determining the scaling factors
   *  and mapping quantities in and out of scaled and unscaled
   *  versions 
   */
  class NLPScalingObject : public ReferencedObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    NLPScalingObject()
    {}

    /** Default destructor */
    virtual ~NLPScalingObject()
    {}
    //@}

    /** Methods to map scaled and unscaled matrices */
    //@{
    /** Returns an obj-scaled version of the given scalar */
    virtual Number apply_obj_scaling(const Number& f)=0;
    /** Returns an obj-unscaled version of the given scalar */
    virtual Number unapply_obj_scaling(const Number& f)=0;
    /** Returns an x-scaled version of the given vector */
    virtual SmartPtr<Vector> apply_vector_scaling_x_NonConst(const SmartPtr<const Vector>& v)=0;
    /** Returns an x-scaled version of the given vector */
    virtual SmartPtr<const Vector> apply_vector_scaling_x(const SmartPtr<const Vector>& v)=0;
    /** Returns an x-unscaled version of the given vector */
    virtual SmartPtr<Vector> unapply_vector_scaling_x_NonConst(const SmartPtr<const Vector>& v)=0;
    /** Returns an x-unscaled version of the given vector */
    virtual SmartPtr<const Vector> unapply_vector_scaling_x(const SmartPtr<const Vector>& v)=0;
    /** Returns an c-scaled version of the given vector */
    virtual SmartPtr<const Vector> apply_vector_scaling_c(const SmartPtr<const Vector>& v)=0;
    /** Returns an c-unscaled version of the given vector */
    virtual SmartPtr<const Vector> unapply_vector_scaling_c(const SmartPtr<const Vector>& v)=0;
    /** Returns an c-scaled version of the given vector */
    virtual SmartPtr<Vector> apply_vector_scaling_c_NonConst(const SmartPtr<const Vector>& v)=0;
    /** Returns an c-unscaled version of the given vector */
    virtual SmartPtr<Vector> unapply_vector_scaling_c_NonConst(const SmartPtr<const Vector>& v)=0;
    /** Returns an d-scaled version of the given vector */
    virtual SmartPtr<const Vector> apply_vector_scaling_d(const SmartPtr<const Vector>& v)=0;
    /** Returns an d-unscaled version of the given vector */
    virtual SmartPtr<const Vector> unapply_vector_scaling_d(const SmartPtr<const Vector>& v)=0;
    /** Returns an d-scaled version of the given vector */
    virtual SmartPtr<Vector> apply_vector_scaling_d_NonConst(const SmartPtr<const Vector>& v)=0;
    /** Returns an d-unscaled version of the given vector */
    virtual SmartPtr<Vector> unapply_vector_scaling_d_NonConst(const SmartPtr<const Vector>& v)=0;
    /** Returns a scaled version of the jacobian for c.
     *  If the overloaded method does not make a new matrix, make sure to set the matrix
     *  ptr passed in to NULL.
     */
    virtual SmartPtr<const Matrix> apply_jac_c_scaling(SmartPtr<const Matrix> matrix)=0;
    /** Returns a scaled version of the jacobian for d 
     *  If the overloaded method does not create a new matrix, make sure to set the matrix
     *  ptr passed in to NULL.
     */
    virtual SmartPtr<const Matrix> apply_jac_d_scaling(SmartPtr<const Matrix> matrix)=0;
    /** Returns a scaled version of the hessian of the lagrangian
     *  If the overloaded method does not create a new matrix, make sure to set the matrix
     *  ptr passed in to NULL.
     */
    virtual SmartPtr<const SymMatrix> apply_hessian_scaling(SmartPtr<const SymMatrix> matrix)=0;
    //@}

    /** Methods for scaling bounds - these wrap those above */
    //@{
    /** Returns an x-scaled vector in the x_L space */
    SmartPtr<Vector> apply_vector_scaling_x_L_NonConst(SmartPtr<Matrix> Px_L, const SmartPtr<const Vector>& l, const SmartPtr<VectorSpace> x_space);
    /** Returns an x-scaled vector in the x_U space */
    SmartPtr<Vector> apply_vector_scaling_x_U_NonConst(SmartPtr<Matrix> Px_U, const SmartPtr<const Vector>& u, const SmartPtr<VectorSpace> x_space);    
    /** Returns an d-scaled vector in the d_L space */
    SmartPtr<Vector> apply_vector_scaling_d_L_NonConst(SmartPtr<Matrix> Pd_L, const SmartPtr<const Vector>& l, const SmartPtr<VectorSpace> d_space);
    /** Returns an d-scaled vector in the d_U space */
    SmartPtr<Vector> apply_vector_scaling_d_U_NonConst(SmartPtr<Matrix> Pd_U, const SmartPtr<const Vector>& u, const SmartPtr<VectorSpace> d_space);    
    //@}

    /** Methods for scaling the gradient of the objective - wraps the
     *  virtual methods above
     */
    //@{
    /** Returns a grad_f scaled version (d_f * D_x^{-1}) of the given vector */
    virtual SmartPtr<Vector> apply_grad_obj_scaling_NonConst(const SmartPtr<const Vector>& v);
    /** Returns a grad_f scaled version (d_f * D_x^{-1}) of the given vector */
    virtual SmartPtr<const Vector> apply_grad_obj_scaling(const SmartPtr<const Vector>& v);
    /** Returns a grad_f unscaled version (d_f * D_x^{-1}) of the given vector */
    virtual SmartPtr<Vector> unapply_grad_obj_scaling_NonConst(const SmartPtr<const Vector>& v);
    /** Returns a grad_f unscaled version (d_f * D_x^{-1}) of the given vector */
    virtual SmartPtr<const Vector> unapply_grad_obj_scaling(const SmartPtr<const Vector>& v);
    //@}

    /** This method is called by the IpoptNLP's at a convenient time to
     *  compute and/or read scaling factors 
     */
    virtual void DetermineScaling(const SmartPtr<const VectorSpace> x_space,
				  const SmartPtr<const VectorSpace> c_space,
				  const SmartPtr<const VectorSpace> d_space,
				  const SmartPtr<const MatrixSpace> jac_c_space,
				  const SmartPtr<const MatrixSpace> jac_d_space,
				  const SmartPtr<const SymMatrixSpace> h_space)=0;

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
    NLPScalingObject(const NLPScalingObject&);

    /** Overloaded Equals Operator */
    void operator=(const NLPScalingObject&);
    //@}
  };


  class NoNLPScalingObject : public NLPScalingObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    NoNLPScalingObject()
    {}

    /** Default destructor */
    virtual ~NoNLPScalingObject()
    {}
    //@}

    /** Methods to map scaled and unscaled matrices - overloaded from NLPScalingObject */
    //@{
    /** Returns an obj-scaled version of the given scalar */
    virtual Number apply_obj_scaling(const Number& f) { return f; }
    /** Returns an obj-unscaled version of the given scalar */
    virtual Number unapply_obj_scaling(const Number& f) { return f; }
    /** Returns an x-scaled version of the given vector */
    virtual SmartPtr<Vector> apply_vector_scaling_x_NonConst(const SmartPtr<const Vector>& v) { return v->MakeNewCopy(); }
    /** Returns an x-scaled version of the given vector */
    virtual SmartPtr<const Vector> apply_vector_scaling_x(const SmartPtr<const Vector>& v) { return v; }
    /** Returns an x-unscaled version of the given vector */
    virtual SmartPtr<Vector> unapply_vector_scaling_x_NonConst(const SmartPtr<const Vector>& v) { return v->MakeNewCopy(); }
    /** Returns an x-unscaled version of the given vector */
    virtual SmartPtr<const Vector> unapply_vector_scaling_x(const SmartPtr<const Vector>& v) { return v; }
    /** Returns an c-scaled version of the given vector */
    virtual SmartPtr<const Vector> apply_vector_scaling_c(const SmartPtr<const Vector>& v) { return v; }
    /** Returns an c-unscaled version of the given vector */
    virtual SmartPtr<const Vector> unapply_vector_scaling_c(const SmartPtr<const Vector>& v) { return v; }
    /** Returns an c-scaled version of the given vector */
    virtual SmartPtr<Vector> apply_vector_scaling_c_NonConst(const SmartPtr<const Vector>& v) { return v->MakeNewCopy(); }
    /** Returns an c-unscaled version of the given vector */
    virtual SmartPtr<Vector> unapply_vector_scaling_c_NonConst(const SmartPtr<const Vector>& v) { return v->MakeNewCopy(); }
    /** Returns an d-scaled version of the given vector */
    virtual SmartPtr<const Vector> apply_vector_scaling_d(const SmartPtr<const Vector>& v) { return v; }
    /** Returns an d-unscaled version of the given vector */
    virtual SmartPtr<const Vector> unapply_vector_scaling_d(const SmartPtr<const Vector>& v)  { return v; }
    /** Returns an d-scaled version of the given vector */
    virtual SmartPtr<Vector> apply_vector_scaling_d_NonConst(const SmartPtr<const Vector>& v) { return v->MakeNewCopy(); }
    /** Returns an d-unscaled version of the given vector */
    virtual SmartPtr<Vector> unapply_vector_scaling_d_NonConst(const SmartPtr<const Vector>& v)  { return v->MakeNewCopy(); }
    /** Returns a scaled version of the jacobian for c.
     *  If the overloaded method does not create a new matrix, make sure to set the matrix
     *  ptr passed in to NULL.
     */
    virtual SmartPtr<const Matrix> apply_jac_c_scaling(SmartPtr<const Matrix> matrix) { SmartPtr<const Matrix> ret = matrix; matrix = NULL; return ret; }
    /** Returns a scaled version of the jacobian for d 
     *  If the overloaded method does not create a new matrix, make sure to set the matrix
     *  ptr passed in to NULL.
     */
    virtual SmartPtr<const Matrix> apply_jac_d_scaling(SmartPtr<const Matrix> matrix)  { SmartPtr<const Matrix> ret = matrix; matrix = NULL; return ret; }
    /** Returns a scaled version of the hessian of the lagrangian
     *  If the overloaded method does not create a new matrix, make sure to set the matrix
     *  ptr passed in to NULL.
     */
    virtual SmartPtr<const SymMatrix> apply_hessian_scaling(SmartPtr<const SymMatrix> matrix) { SmartPtr<const SymMatrix> ret = matrix; matrix = NULL; return ret; }
    //@}

    /** Method to determine the scaling - overloaded from NLPScaling.
     *  Since no scaling is done, we do not need to determine anything
     */
    virtual void DetermineScaling(const SmartPtr<const VectorSpace> x_space,
				  const SmartPtr<const VectorSpace> c_space,
				  const SmartPtr<const VectorSpace> d_space,
				  const SmartPtr<const MatrixSpace> jac_c_space,
				  const SmartPtr<const MatrixSpace> jac_d_space,
				  const SmartPtr<const SymMatrixSpace> h_space)
    {}

    /** Methods for scaling the gradient of the objective - made efficient for
     *  the NoScaling implementation.
     */
    //@{
    /** Returns a grad_f scaled version (d_f * D_x^{-1}) of the given vector */
    virtual SmartPtr<Vector> apply_grad_obj_scaling_NonConst(const SmartPtr<const Vector>& v) { return v->MakeNewCopy(); }
    /** Returns a grad_f scaled version (d_f * D_x^{-1}) of the given vector */
    virtual SmartPtr<const Vector> apply_grad_obj_scaling(const SmartPtr<const Vector>& v) { return v; }
    /** Returns a grad_f unscaled version (d_f * D_x^{-1}) of the given vector */
    virtual SmartPtr<Vector> unapply_grad_obj_scaling_NonConst(const SmartPtr<const Vector>& v) { return v->MakeNewCopy(); }
    /** Returns a grad_f unscaled version (d_f * D_x^{-1}) of the given vector */
    virtual SmartPtr<const Vector> unapply_grad_obj_scaling(const SmartPtr<const Vector>& v) { return v; }
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
    NoNLPScalingObject(const NoNLPScalingObject&);

    /** Overloaded Equals Operator */
    void operator=(const NoNLPScalingObject&);
    //@}
  };

} // namespace Ipopt

#endif
