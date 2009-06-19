// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-05-29

#ifndef __IPPARGENMATRIX_HPP__
#define __IPPARGENMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpParVector.hpp"
#include "IpGenTMatrix.hpp"
#include <vector>

namespace Ipopt
{

  /* forward declarations */
  class ParGenMatrixSpace;

  /** Class of parallel general matrices. A ParGenMatrix consists of a
   *  matrix whose rows are partitioned across the universe of
   *  processors (assume n in number).  Each processor is responsible
   *  for its piece of the entire matrix, and for manipulating this
   *  piece. All ParGenMatrix functions are assumed to be executed
   *  simultaneously at each of the n processors. This class does not
   *  support the functionality of a matrix being partitioned across a
   *  subset of processors. Thus if the number of rows of the
   *  ParGenMatrix is less than n, some processor (say k) will not
   *  have any rows of the matrix.  Yet a function to be executed on a
   *  ParGenMatrix must also be executed in processor k. As an
   *  example, if the maximum element of a column of a ParGenMatrix is
   *  to be computed, each processor finds the maximum element of its
   *  piece of the column, and then uses MPI calls to find the global
   *  maximum.
   *
   *  The rows of a ParGenMatrix have the same numbering convention as
   *  a GenTMatrix, i.e., the rows are numbered from 1. A ParGenMatrix
   *  stores the piece of the matrix on a processor as a
   *  GenTMatrix. Therefore the collection of row numbers on a
   *  processor need to be passed relative to the starting row
   *  number. For example, if the part of the matrix in a processor
   *  consists of rows 15-20, then the row numbers 15 - 20 should be
   *  specified as row 1, 2, .., 5.
   */
  class ParGenMatrix : public Matrix
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor, given the corresponding ParGenMatrixSpace.
     */
    ParGenMatrix(const ParGenMatrixSpace* owner_space);

    /** Default destructor */
    virtual ~ParGenMatrix();
    //@}

    /** @name Additional public methods not in Matrix base class. */
    //@{
    /** Create a new ParGenMatrix from same MatrixSpace */
    SmartPtr<ParGenMatrix> MakeNewParGenMatrix() const;

    /** Obtain pointer to local matrix. 
     */
    inline GenTMatrix* LocalMatrix()
    {
      ObjectChanged();
      return GetRawPtr(local_matrix_);
    }

    /** Obtain pointer to const local matrix. 
     */
    inline const GenTMatrix* LocalMatrix() const
    {
      return GetRawPtr(local_matrix_);
    }

    /** Rank of the processor handling this vector */
    int Rank() const;

    /** Total number of processors. */
    int NumProc() const;

    /** number of rows of local matrix */
    int LocalNRows() const;

    /** get offset of rows in local matrix */
    int RowStartPos() const;
    //@}

  protected:
    /**@name Overloaded methods from Matrix base class*/
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
     * they will not be implicitly created/called.
     */
    //@{
    /** Default Constructor */
    ParGenMatrix();

    /** Copy Constructor */
    ParGenMatrix(const ParGenMatrix&);

    /** Overloaded Equals Operator */
    void operator=(const ParGenMatrix&);
    //@}

    /** Copy of the owner_space ptr as a ParGenMatrixSpace instead
     *  of a MatrixSpace
     */
    const ParGenMatrixSpace* owner_space_;

    /** local piece of ParGenMatrix stored as a GenTMatrix */
    SmartPtr<GenTMatrix> local_matrix_;
  };

  /** This matrix space is the matrix space for ParGenMatrix.
   *  Before a ParGenMatrix can be created, all information relevant to
   *  a ParGenMatrix, such as dimension of the space, and the range of
   *  indices on the current processor have to be set in the ParGenMatrixSpace.
   */
  class ParGenMatrixSpace : public MatrixSpace
  {
  public:
    /** @name Constructors/Destructors. */
    //@{
    /** Constructor, has to be given the total ndimension of the matrix space
     *  and the local range of contiguous rows + local matrix space information */
    ParGenMatrixSpace(SmartPtr<const ParVectorSpace> RowVectorSpace, Index nCols,
		      Index nonZeros,
		      const Index* iRows, const Index* jCols);

    /** Destructor */
    ~ParGenMatrixSpace()
    {}
    //@}

    /** Method for creating a new matrix of this specific type. */
    virtual ParGenMatrix* MakeNewParGenMatrix() const
    {
      return new ParGenMatrix(this);
    }

    /** Overloaded MakeNew method for the MatrixSpace base class.
     */
    virtual Matrix* MakeNew() const
    {
      return MakeNewParGenMatrix();
    }

    SmartPtr<const ParVectorSpace> RowVectorSpace() const
    {
      return rowVectorSpace_;
    }

    SmartPtr<GenTMatrixSpace> LocalSpace() const
    {
      return local_space_;
    }

    int RowStartPos() const
    {
      return rowVectorSpace_->StartPos();
    }

    const std::vector<int> & RecvCounts() const
    {
      return rowVectorSpace_->RecvCounts();
    }

    const std::vector<int> & Displs() const
    {
      return rowVectorSpace_->Displs();
    }

    int LocalNRows() const
    {
      return rowVectorSpace_->LocalSize();
    }

    int Rank() const
    {
      return rowVectorSpace_->Rank();
    }

    int NumProc() const
    {
      return rowVectorSpace_->NumProc();
    }

  private:
    SmartPtr<const ParVectorSpace> rowVectorSpace_;
    SmartPtr<GenTMatrixSpace> local_space_;
  };

  inline
  SmartPtr<ParGenMatrix> ParGenMatrix::MakeNewParGenMatrix() const
  {
    return owner_space_->MakeNewParGenMatrix();
  }

  inline
  int ParGenMatrix::Rank() const
  {
    return owner_space_->Rank();
  }

  inline
  int ParGenMatrix::NumProc() const
  {
    return owner_space_->NumProc();
  }

  inline
  int ParGenMatrix::LocalNRows() const
  {
    return owner_space_->LocalNRows();
  }

  inline
  int ParGenMatrix::RowStartPos() const
  {
    return owner_space_->RowStartPos();
  }

} // namespace Ipopt

#endif
