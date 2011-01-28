// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter             IBM    2009-11-05

#ifndef __IPEXPANDEDMULTIVECTORMATRIX_HPP__
#define __IPEXPANDEDMULTIVECTORMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpMatrix.hpp"
#include "IpExpansionMatrix.hpp"

namespace Ipopt
{

  /** forward declarations */
  class ExpandedMultiVectorMatrixSpace;

  /** Class for Matrices with few rows that consists of Vectors,
   *  together with a premultiplied Expansion matrix.  So, the matrix
   *  is V^T*P^T.  If P is NULL, it is assumed to be the identity
   *  matrix.  If a row vector of V is NULL, it is assumed to be all
   *  zero.  This is used to construct the KKT system with low-rank
   *  Hessian approximation.
   */
  class ExpandedMultiVectorMatrix : public Matrix
  {
  public:

    /**@name Constructors / Destructors */
    //@{

    /** Constructor, taking the owner_space.
     */
    ExpandedMultiVectorMatrix(const ExpandedMultiVectorMatrixSpace* owner_space);

    /** Destructor */
    virtual ~ExpandedMultiVectorMatrix()
    {}
    //@}

    SmartPtr<ExpandedMultiVectorMatrix> MakeNewExpandedMultiVectorMatrix() const;

    /** Set a particular Vector at a given row position, replacing
     *  another vector if there has been one.  */
    void SetVector(Index i, SmartPtr<const Vector> vec);

    /** Get a Vector in a particular row as a const Vector */
    inline SmartPtr<const Vector> GetVector(Index i) const
    {
      DBG_ASSERT(i < NRows());
      return vecs_[i];
    }

    /** Vector space for the rows */
    SmartPtr<const VectorSpace> RowVectorSpace() const;

    /** Return the ExpandedMultiVectorMatrixSpace */
    SmartPtr<const ExpandedMultiVectorMatrixSpace>
    ExpandedMultiVectorMatrixOwnerSpace() const;

    /** Return the Expansion matrix.  If NULL, there is no expansion,
     *  the vector is used as is. */
    SmartPtr<const ExpansionMatrix> GetExpansionMatrix() const;

  protected:
    /**@name Overloaded methods from Matrix base class */
    //@{
    virtual void MultVectorImpl(Number alpha, const Vector &x, Number beta,
                                Vector &y) const;

    virtual void TransMultVectorImpl(Number alpha, const Vector& x,
                                     Number beta, Vector& y) const;

    /** Method for determining if all stored numbers are valid (i.e.,
     *  no Inf or Nan). */
    virtual bool HasValidNumbersImpl() const;

    virtual void ComputeRowAMaxImpl(Vector& rows_norms, bool init) const;

    virtual void ComputeColAMaxImpl(Vector& cols_norms, bool init) const;

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
    ExpandedMultiVectorMatrix();

    /** Copy Constructor */
    ExpandedMultiVectorMatrix(const ExpandedMultiVectorMatrix&);

    /** Overloaded Equals Operator */
    void operator=(const ExpandedMultiVectorMatrix&);
    //@}

    const ExpandedMultiVectorMatrixSpace* owner_space_;

    /** space for storing the const Vector's */
    std::vector<SmartPtr<const Vector> > vecs_;

  };

  /** This is the matrix space for ExpandedMultiVectorMatrix.
   */
  class ExpandedMultiVectorMatrixSpace : public MatrixSpace
  {
  public:
    /** @name Constructors / Destructors */
    //@{
    /** Constructor, given the number of rows (i.e., Vectors to be
     *  stored) and given the VectorSpace for the Vectors.
     */
    ExpandedMultiVectorMatrixSpace(Index nrows,
                                   const VectorSpace& vec_space,
                                   SmartPtr<const ExpansionMatrix> exp_matrix);
    /** Destructor */
    virtual ~ExpandedMultiVectorMatrixSpace()
    {}
    //@}

    /** Method for creating a new matrix of this specific type. */
    ExpandedMultiVectorMatrix* MakeNewExpandedMultiVectorMatrix() const
    {
      return new ExpandedMultiVectorMatrix(this);
    }

    /** Overloaded MakeNew method for the MatrixSpace base class.
     */
    virtual Matrix* MakeNew() const
    {
      return MakeNewExpandedMultiVectorMatrix();
    }

    /** Accessor method for the VectorSpace for the rows */
    SmartPtr<const VectorSpace> RowVectorSpace() const
    {
      return vec_space_;
    }

    SmartPtr<const ExpansionMatrix> GetExpansionMatrix() const
    {
      return exp_matrix_;
    }

  private:
    SmartPtr<const VectorSpace> vec_space_;

    SmartPtr<const ExpansionMatrix> exp_matrix_;
  };

  inline
  SmartPtr<ExpandedMultiVectorMatrix> ExpandedMultiVectorMatrix::MakeNewExpandedMultiVectorMatrix() const
  {
    return owner_space_->MakeNewExpandedMultiVectorMatrix();
  }

  inline
  SmartPtr<const VectorSpace> ExpandedMultiVectorMatrix::RowVectorSpace() const
  {
    return owner_space_->RowVectorSpace();
  }

  inline
  SmartPtr<const ExpansionMatrix> ExpandedMultiVectorMatrix::GetExpansionMatrix() const
  {
    return owner_space_->GetExpansionMatrix();
  }

  inline
  SmartPtr<const ExpandedMultiVectorMatrixSpace>
  ExpandedMultiVectorMatrix::ExpandedMultiVectorMatrixOwnerSpace() const
  {
    return owner_space_;
  }

} // namespace Ipopt
#endif
