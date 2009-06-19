// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-06-17

#ifndef __IPPAREXPANSIONMATRIX_HPP__
#define __IPPAREXPANSIONMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpParVector.hpp"
#include "IpExpansionMatrix.hpp"

namespace Ipopt
{

  /** forward declarations */
  class ParExpansionMatrixSpace;

  /** Class for parallel expansion/projection matrices. Each processor
   *  stores an expansion matrix, and the collection of expansion
   *  matrices define the global expansion matrix, which can be viewed
   *  as mapping a vector to another with the same non-zero elements in
   *  a higher dimensional space.
   * 
   *  It is important to note that the rows are assumed to be numbered
   *  from zero. Secondly, the non-zero rows of the expansion matrix on
   *  also start from zero, i.e., they are given relative to a starting
   *  row for each processor.
   */
  class ParExpansionMatrix : public Matrix
  {
  public:

    /**@name Constructors / Destructors */
    //@{

    /** Constructor, taking the owner_space.
     */
    ParExpansionMatrix(const ParExpansionMatrixSpace* owner_space);

    /** Destructor */
    ~ParExpansionMatrix();
    //@}

    /** @name Additional public methods not in Matrix base class. */
    //@{
    /** Create a new ParExpansionMatrix from same MatrixSpace */
    SmartPtr<ParExpansionMatrix> MakeNewParExpansionMatrix() const;

    /** Obtain pointer to local matrix. 
     */
    inline ExpansionMatrix* LocalMatrix()
    {
      ObjectChanged();
      return GetRawPtr(local_matrix_);
    }

    /** Obtain pointer to const local matrix. 
     */
    inline const ExpansionMatrix* LocalMatrix() const
    {
      return GetRawPtr(local_matrix_);
    }

    /** Offset of row indices for local part */
    Index RowStartPos() const;

    /** Offset of column indices for local part */
    Index ColStartPos() const;
    //@}

  protected:
    /**@name Overloaded methods from Matrix base class*/
    //@{
    virtual void MultVectorImpl(Number alpha, const Vector &x, Number beta,
                                Vector &y) const;

    virtual void TransMultVectorImpl(Number alpha, const Vector& x,
                                     Number beta, Vector& y) const;

    /** X = beta*X + alpha*(Matrix S^{-1} Z).  Specialized implementation.
     */
    virtual void AddMSinvZImpl(Number alpha, const Vector& S, const Vector& Z,
                               Vector& X) const;

    /** X = S^{-1} (r + alpha*Z*M^Td).  Specialized implementation.
     */
    virtual void SinvBlrmZMTdBrImpl(Number alpha, const Vector& S,
                                    const Vector& R, const Vector& Z,
                                    const Vector& D, Vector& X) const;

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
    ParExpansionMatrix();

    /** Copy Constructor */
    ParExpansionMatrix(const ParExpansionMatrix&);

    /** Overloaded Equals Operator */
    void operator=(const ParExpansionMatrix&);
    //@}

    /** local piece of ParExpansionMatrix stored as an ExpansionMatrix */
    SmartPtr<ExpansionMatrix> local_matrix_;

    /** Copy of the owner_space ptr as a ParExpansionMatrixSpace instead
     *  of a MatrixSpace
     */
    const ParExpansionMatrixSpace* owner_space_;
  };

  /** This is the matrix space for ParExpansionMatrix.
   */
  class ParExpansionMatrixSpace : public MatrixSpace
  {
  public:
    /** @name Constructors / Destructors */
    //@{
    /** Constructor, given the list of elements of the local piece of
     *  the large vector to be filtered into the local piece of the
     *  small vector.
     */
    ParExpansionMatrixSpace(SmartPtr<const ParVectorSpace> LargeVectorSpace, 
			    SmartPtr<const ParVectorSpace> SmallVectorSpace, 
			    const Index *ExpPos,
			    const int offset = 0);

    /** Destructor */
    ~ParExpansionMatrixSpace()
    {}
    //@}

    /** Method for creating a new matrix of this specific type. */
    ParExpansionMatrix* MakeNewParExpansionMatrix() const
    {
      return new ParExpansionMatrix(this);
    }

    /** Overloaded MakeNew method for the MatrixSpace base class.
     */
    virtual Matrix* MakeNew() const
    {
      return MakeNewParExpansionMatrix();
    }

    SmartPtr<ExpansionMatrixSpace> LocalSpace() const
    {
      return local_space_;
    }

    SmartPtr<const ParVectorSpace> LargeVectorSpace() const
    {
      return large_vector_space_;
    }

    SmartPtr<const ParVectorSpace> SmallVectorSpace() const
    {
      return small_vector_space_;
    }

  private:
    SmartPtr<ExpansionMatrixSpace> local_space_;
    SmartPtr<const ParVectorSpace> large_vector_space_;
    SmartPtr<const ParVectorSpace> small_vector_space_;
  };

  inline Index ParExpansionMatrix::RowStartPos() const
  {
    return owner_space_->LargeVectorSpace()->StartPos();
  }

  inline Index ParExpansionMatrix::ColStartPos() const
  {
    return owner_space_->SmallVectorSpace()->StartPos();
  }

} // namespace Ipopt
#endif
