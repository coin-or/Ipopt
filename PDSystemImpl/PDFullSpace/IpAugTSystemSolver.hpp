// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IP_AUGTSYSTEMSOLVER_HPP__
#define __IP_AUGTSYSTEMSOLVER_HPP__

#include "IpAugSystemSolver.hpp"
#include "IpSymTMatrix.hpp"
#include "IpGenTMatrix.hpp"
#include "IpTaggedObject.hpp"
#include "IpDenseVector.hpp"

namespace Ipopt
{

  /** Solver for the augmented system for triple type matrices.
   *
   *  The current implemetation assumes that all matrices are of the
   *  type SymTMatrix, and all vectors are of the type DenseVector.
   */
  class AugTSystemSolver : public AugSystemSolver
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor using only a linear solver object */
    AugTSystemSolver(SymLinearSolver& LinSolver);

    /** Default destructor */
    virtual ~AugTSystemSolver();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);

    /** Set up the augmented system and solve it for a given right hand
     *  side - implementation for GenTMatrices and SymTMatrices.
     */
    virtual SymLinearSolver::ESolveStatus Solve(
      const SymMatrix* W,
      const Vector* D_x,
      double delta_x,
      const Vector* D_s,
      double delta_s,
      const Matrix* J_c,
      const Vector* D_c,
      double delta_c,
      const Matrix* J_d,
      const Vector* D_d,
      double delta_d,
      const Vector& rhs_x,
      const Vector& rhs_s,
      const Vector& rhs_c,
      const Vector& rhs_d,
      Vector& sol_x,
      Vector& sol_s,
      Vector& sol_c,
      Vector& sol_d,
      bool check_NegEVals,
      Index numberOfNegEVals);

    /** Number of negative eigenvalues detected during last
     * solve.  Returns the number of negative eigenvalues of
     * the most recent factorized matrix.  This must not be called if
     * the linear solver does not compute this quantities (see
     * ProvidesInertia).
     */
    virtual Index NumberOfNegEVals() const;

    /** Query whether inertia is computed by linear solver.
     * Returns true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const;

    /** Request to increase quality of solution for next solve.  Ask
     *  underlying linear solver to increase quality of solution for
     *  the next solve (e.g. increase pivot tolerance).  Returns
     *  false, if this is not possible (e.g. maximal pivot tolerance
     *  already used.)
     */
    virtual bool IncreaseQuality();

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default constructor. */
    AugTSystemSolver();
    /** Copy Constructor */
    AugTSystemSolver(const AugTSystemSolver&);

    /** Overloaded Equals Operator */
    void operator=(const AugTSystemSolver&);
    //@}

    /** Initilialize the internal structures of this object.  Allocate
     *	the matrix that stores the augmented system.  This method must
     *  be called before any call to the Set method, and it can only
     *  called once.
     */
    void InternalInitialize(const SymMatrix* W,
                            const Matrix* J_c,
                            const Matrix* J_d);

    /** Set the components of the augmented Matrix.  Here, if W, A,
     *  D_W or D_C are given as NULL, they * are assume to be the null
     *  matrix.
     */
    void Set(const SymMatrix* W,
             const Vector* D_x,
             double delta_x,
             const Vector* D_s,
             double delta_s,
             const Matrix* J_c,
             const Vector* D_c,
             double delta_c,
             const Matrix* J_d,
             const Vector* D_d,
             double delta_d);

    /** The linear solver object that is to be used to solve the
     *  linear systems.
     */
    SmartPtr<SymLinearSolver> linsolver_;

    /** Spaces for creating augmented vectors and matrices */
    SmartPtr<DenseVectorSpace> aug_vec_space_;
    SmartPtr<SymTMatrixSpace> aug_mat_space_;

    /**@name Tags and values to track in order to decide whether the
       matrix has to be updated compared to the most recent call of
       the Set method.
     */
    //@{
    /** Tag for W matrix.  If W has been given to Set as NULL, then
     *  this tag is set to 0
     */
    TaggedObject::Tag w_tag_;
    /** Tag for D_x vector, representing the diagonal matrix D_x.  If
     *  D_x has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_x_tag_;
    /** Most recent value of delta_x from Set method */
    double delta_x_;
    /** Tag for D_s vector, representing the diagonal matrix D_s.  If
     *  D_s has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_s_tag_;
    /** Most recent value of delta_s from Set method */
    double delta_s_;
    /** Tag for J_c matrix.  If J_c has been given to Set as NULL, then
     *  this tag is set to 0
     */
    TaggedObject::Tag j_c_tag_;
    /** Tag for D_c vector, representing the diagonal matrix D_c.  If
     *  D_c has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_c_tag_;
    /** Most recent value of delta_c from Set method */
    double delta_c_;
    /** Tag for J_d matrix.  If J_d has been given to Set as NULL, then
     *  this tag is set to 0
     */
    TaggedObject::Tag j_d_tag_;
    /** Tag for D_d vector, representing the diagonal matrix D_d.  If
     *  D_d has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_d_tag_;
    /** Most recent value of delta_d from Set method */
    double delta_d_;

    /** This is the tag of the matrix storing the augmented system.  Since
     *  this object owns this matrix, no changes should happen outside.
     *  However, since it is given away as a smart pointer, someone outside
     *  might change it.  For debugging purposes, we now track its tag as
     *  well.
     */
    TaggedObject::Tag augsys_tag_;
    //@}

    /** @name Information about the number of nonzeros for the
     *  different parts of the augmented system.
     */
    //@{
    /** Number of optimization variables x
     */
    Index nx_;
    /** Number of inequality constraints (same as number of variables s).
     */
    Index nd_;
    /** Number of equality constraints
     */
    Index nc_;

    /** Number of nonzeros in the Hessian part.  This value is set to
     *	0 in the Initialize method, if the upper left block of the
     *  augmented system is only the diagonal matrix.
     */
    Index nnz_w_;
    /** Number of nonzeros in the equality constraint Jacobian part.
     *	This value is set to 0 in the Initialize method, if J_c was passed
     *  as zero (i.e. n_c_ = 0, too).
     */
    Index nnz_j_c_;
    /** Number of nonzeros in the inequality constraint Jacobian part.
     *	This value is set to 0 in the Initialize method, if J_d was passed
     *  as zero (i.e. n_s_ = 0, too).
     */
    Index nnz_j_d_;
    //@}

    /** The resulting augmented matrix.
     *  This matrix is stored as follows:  First we have the diagonal elements
     *  for the upper left block (for D_W and delta_W), then the elements for
     *  the Hessian W, then the Jacobian A, and finally the diagonal elements
     *  for the lower right block (for D_C and delta_C).
     */
    SmartPtr<SymTMatrix> augsystem_;

    /** Flag indicating whether the data structures to store the
     *	augmented system have already been allcated.
     */
    bool initialized_;

    /** Flag indicating whether values for the augmented matrix have
     *  ever been evalutated
     */
    bool have_values_;

    /**@name local functions for convenience */
    //@{
    /** Little function to decide whether a part of the matrix has changed */
    bool ChangeInPart(const TaggedObject* Obj,
                      const TaggedObject::Tag tag);

    /** @name Little helper functions to copy the matrix values into the
     * augmented system matrix
     */
    //@{
    void FillDPart(Index len, Number delta,
                   const Vector* D, Number* vals2fill);
    void FillWPart(Index len, const SymMatrix* W,
                   Number* vals2fill);
    void FillJPart(Index len, const Matrix* J,
                   Number* vals2fill);
    //@}
    /** Tiny helper function to determine the value of the tag that
     * should be stored for the next call.  If the incoming pointer is Null,
     * this returns 0, otherwise, it returns the GetTag value of the Obj
     */
    TaggedObject::Tag NewTag(const TaggedObject* Obj);
  };

} // namespace Ipopt

#endif
