// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-06-11

#ifndef __IPPARSYMTMATRIX_HPP__
#define __IPPARSYMTMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpSymTMatrix.hpp"

namespace Ipopt
{

  /* forward declarations */
  class ParSymMatrixSpace;

  /** Class for parallel symmetric matrices. A ParSymMatrix is partitioned
   *  across processors in a different way than a ParGenMatrix. If the
   *  the matrix to be partitioned is an m X m matrix, and there are n
   *  processors, then each processor contains an m X m symmetric matrix
   *  stored in triplet form, and the ParSymMatrix is the sum of the
   *  n different symmetric matrices. Obviously, the most effective use
   *  of this class would be when the symmetric matrix pieces in the different
   *  processors do not have nonzero values in common positions.
   *  Any function to be executed on the ParSymMatrix has to be executed in
   *  each processor, even if the symmetric matrix piece on a processor has
   *  no nonzero values.
   *
   *  It stores the local piece of the matrix as a SymTMatrix, and therefore
   *  the rows and columns of the non-zero values start from 1.
   */
  class ParSymMatrix : public SymMatrix
  {
  public:
    /**@name Constructors / Destructors */
    //@{
    /** Constructor, given the corresponding matrix space.
     */
    ParSymMatrix(const ParSymMatrixSpace* owner_space);

    /** Destructor */
    ~ParSymMatrix();
    //@}

    /** @name Additional public methods not in Matrix base class. */
    //@{
    /** Create a new ParSymMatrix from same MatrixSpace */
    SmartPtr<ParSymMatrix> MakeNewParSymMatrix() const;

    /** Obtain pointer to local matrix. 
     */
    inline SymTMatrix* LocalMatrix()
    {
      ObjectChanged();
      return GetRawPtr(local_matrix_);
    }

    /** Obtain pointer to const local matrix. 
     */
    inline const SymTMatrix* LocalMatrix() const
    {
      return GetRawPtr(local_matrix_);
    }

    /** Rank of the processor handling this vector */
    int Rank() const;

    /** Total number of processors. */
    int NumProc() const;
    //@}

  protected:
    /**@name Methods overloaded from matrix */
    //@{
    /** Since the matrix is symetric, this is the same operation as
     *  MultVector */
    virtual void TransMultVectorImpl(Number alpha, const Vector& x, Number beta,
                                     Vector& y) const
    {
      MultVectorImpl(alpha, x, beta, y);
    }

    virtual void MultVectorImpl(Number alpha, const Vector& x, Number beta,
                                Vector& y) const;

    /** Method for determining if all stored numbers are valid (i.e.,
     *  no Inf or Nan). */
    virtual bool HasValidNumbersImpl() const;

    virtual void ComputeRowAMaxImpl(Vector& rows_norms, bool init) const
    {}

    virtual void ComputeColAMaxImpl(Vector& cols_norms, bool init) const
    {}

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
    ParSymMatrix();

    /** Copy Constructor */
    ParSymMatrix(const ParSymMatrix&);

    /** Overloaded Equals Operator */
    void operator=(const ParSymMatrix&);
    //@}

    /** Copy of the owner_space ptr as a ParSymMatrixSpace instead
     *  of a MatrixSpace
     */
    const ParSymMatrixSpace* owner_space_;

    Index Nonzeros() const;

    const Index* Irows() const;

    const Index* Jcols() const;

    /** local piece of ParSymMatrix stored as a SymTMatrix */
    SmartPtr<SymTMatrix> local_matrix_;
  };

  /** This is the matrix space for a ParSymMatrix with fixed sparsity
   *  structure per processor. The sparsity structure is stored here in
   *  the matrix space.
   */
  class ParSymMatrixSpace : public SymMatrixSpace
  {
  public:
    /** @name Constructors / Destructors */
    //@{
    /** Constructor, has to be given the total ndimension of the matrix space
     *  and the local sparsity structure. Rows and Columns are numbered from 1.
     */
    ParSymMatrixSpace(Index dim, Index nonZeros, const Index* iRows,
		      const Index* jCols);

    /** Destructor */
    ~ParSymMatrixSpace()
    {}
    //@}

    /** Overloaded MakeNew method for the SymMatrixSpace base class.
     */
    virtual SymMatrix* MakeNewSymMatrix() const
    {
      return MakeNewParSymMatrix();
    }

    /** Method for creating a new matrix of this specific type. */
    virtual ParSymMatrix* MakeNewParSymMatrix() const
    {
      return new ParSymMatrix(this);
    }

    SmartPtr<SymTMatrixSpace> getLocalSpace() const
    {
      return local_space_;
    }

    int Rank() const
    {
      return rank_;
    }

    int NumProc() const
    {
      return num_proc_;
    }

    /**@name Methods describing Matrix structure */
    //@{
    /** Number of non-zeros in the sparse matrix */
    Index Nonzeros() const
    {
      return local_space_->Nonzeros();
    }

    /** Row index of each non-zero element */
    const Index* Irows() const
    {
      return local_space_->Irows();
    }

    /** Column index of each non-zero element */
    const Index* Jcols() const
    {
      return local_space_->Jcols();
    }
    //@}

  private:
    int rank_;
    int num_proc_;

    SmartPtr<SymTMatrixSpace> local_space_;
  };

  /* Inline Methods */
  inline
  SmartPtr<ParSymMatrix> ParSymMatrix::MakeNewParSymMatrix() const
  {
    return owner_space_->MakeNewParSymMatrix();
  }

  inline
  int ParSymMatrix::Rank() const
  {
    return owner_space_->Rank();
  }

  inline
  int ParSymMatrix::NumProc() const
  {
    return owner_space_->NumProc();
  }

  inline
  Index ParSymMatrix::Nonzeros() const
  {
    return owner_space_->Nonzeros();
  }

  inline
  const Index* ParSymMatrix::Irows() const
  {
    return owner_space_->Irows();
  }

  inline
  const Index* ParSymMatrix::Jcols() const
  {
    return owner_space_->Jcols();
  }

} // namespace Ipopt
#endif
