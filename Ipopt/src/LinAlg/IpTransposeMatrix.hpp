// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2008-08-25

#ifndef __IPTRANSPOSEMATRIX_HPP__
#define __IPTRANSPOSEMATRIX_HPP__

#include "IpMatrix.hpp"

namespace Ipopt
{

  /* forward declarations */
  class TransposeMatrixSpace;

  /** Class for Matrices which are the transpose of another matrix.
   *
   */
  class TransposeMatrix : public Matrix
  {
  public:

    /**@name Constructors / Destructors */
    //@{

    /** Constructor, initializing with dimensions of the matrix.
     */
    TransposeMatrix(const TransposeMatrixSpace* owner_space);

    /** Destructor */
    ~TransposeMatrix()
    {}

    SmartPtr<const Matrix> OrigMatrix() const
    {
      return ConstPtr(orig_matrix_);
    }
    //@}

  protected:
    /**@name Methods overloaded from matrix */
    //@{
    virtual void MultVectorImpl(Number alpha, const Vector& x,
                                Number beta, Vector& y) const
    {
      DBG_ASSERT(IsValid(orig_matrix_));
      orig_matrix_->TransMultVector(alpha, x, beta, y);
    }

    virtual void TransMultVectorImpl(Number alpha, const Vector& x,
                                     Number beta, Vector& y) const
    {
      DBG_ASSERT(IsValid(orig_matrix_));
      orig_matrix_->MultVector(alpha, x, beta, y);
    }

    /** Method for determining if all stored numbers are valid (i.e.,
     *  no Inf or Nan). */
    virtual bool HasValidNumbersImpl() const
    {
      DBG_ASSERT(IsValid(orig_matrix_));
      return orig_matrix_->HasValidNumbers();
    }

    virtual void ComputeRowAMaxImpl(Vector& rows_norms, bool init) const
    {
      DBG_ASSERT(IsValid(orig_matrix_));
      orig_matrix_->ComputeColAMax(rows_norms, init);
    }

    virtual void ComputeColAMaxImpl(Vector& rows_norms, bool init) const
    {
      DBG_ASSERT(IsValid(orig_matrix_));
      orig_matrix_->ComputeRowAMax(rows_norms, init);
    }

    virtual void PrintImpl(const Journalist& jnlst,
                           EJournalLevel level,
                           EJournalCategory category,
                           const std::string& name,
                           Index indent,
                           const std::string& prefix) const;
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
    TransposeMatrix();

    /** Copy Constructor */
    TransposeMatrix(const TransposeMatrix&);

    /** Overloaded Equals Operator */
    void operator=(const TransposeMatrix&);
    //@}

    /** Pointer to original matrix */
    SmartPtr<Matrix> orig_matrix_;
  };

  /** This is the matrix space for TransposeMatrix. */
  class TransposeMatrixSpace : public MatrixSpace
  {
  public:
    /** @name Constructors / Destructors */
    //@{
    /** Constructor, given the dimension of the matrix. */
    TransposeMatrixSpace(const MatrixSpace* orig_matrix_space)
        :
        MatrixSpace(orig_matrix_space->NCols(), orig_matrix_space->NRows()),
        orig_matrix_space_(orig_matrix_space)
    {}

    /** Destructor */
    virtual ~TransposeMatrixSpace()
    {}
    //@}

    /** Overloaded MakeNew method for the MatrixSpace base class.
     */
    virtual Matrix* MakeNew() const
    {
      return MakeNewTransposeMatrix();
    }

    /** Method for creating a new matrix of this specific type. */
    TransposeMatrix* MakeNewTransposeMatrix() const
    {
      return new TransposeMatrix(this);
    }

    Matrix* MakeNewOrigMatrix() const
    {
      return orig_matrix_space_->MakeNew();
    }

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
    TransposeMatrixSpace();

    /** Copy Constructor */
    TransposeMatrixSpace(const TransposeMatrixSpace&);

    /** Overloaded Equals Operator */
    void operator=(const TransposeMatrixSpace&);
    //@}

    /** Matrix space of the original matrix */
    SmartPtr<const MatrixSpace> orig_matrix_space_;
  };

} // namespace Ipopt
#endif
