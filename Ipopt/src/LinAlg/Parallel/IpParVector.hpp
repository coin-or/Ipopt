// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-05-12

#ifndef __IPPARVECTOR_HPP__
#define __IPPARVECTOR_HPP__

#include "IpUtils.hpp"
#include "IpDenseVector.hpp"
#include <vector>

namespace Ipopt
{

  /* forward declarations */
  class ParVectorSpace;

  /** Class of parallel vectors. A ParVector consists of a vector which
   *  is partitioned across the universe of processors (assume n in number).
   *  Each processor is responsible for its piece of the entire vector,
   *  and for manipulating this piece. All ParVector functions are assumed
   *  to be executed simultaneously at each of the n processors. This class
   *  does not support the functionality of a vector being partitioned across
   *  a subset of processors. Thus if the dimension of the ParVector is less
   *  than n, some processor (say k) will not have any elements of the vector.
   *  Yet a function to be executed on a ParVector must also be executed in
   *  processor k. As an example, if the maximum element of a ParVector is to
   *  be computed, each processor finds the maximum element of its piece of the
   *  vector, and then uses MPI calls to find the global maximum.
   *  A ParVector is similar to a CompoundVector in that it is implemented as
   *  a block vector. However, each processor can only directly access a single
   *  component or block.
   *
   *  It stores the local piece of the vector in as a DenseVector. An important
   *  point to keep in mind is that rows in a ParVector are numbered from 0.
   */
  class ParVector : public Vector
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor, given the corresponding ParVectorSpace.
     */
    ParVector(const ParVectorSpace* owner_space);

    /** Default destructor */
    virtual ~ParVector();
    //@}

    /** @name Additional public methods not in Vector base class. */
    //@{
    /** Create a new ParVector from same VectorSpace */
    SmartPtr<ParVector> MakeNewParVector() const;

    /** Create a new DenseVector from GlobalVectorSpace */
    SmartPtr<DenseVector> MakeNewGlobalVector() const;

    /** Obtain pointer to local vector. 
     */
    inline DenseVector* LocalVector()
    {
      ObjectChanged();
      return GetRawPtr(local_vector_);
    }

    /** Obtain pointer to const local vector. 
     */
    inline const DenseVector* LocalVector() const
    {
      return GetRawPtr(local_vector_);
    }

    /** Obtain pointer to const global vector. 
     */
    // TODO: Maybe we want to cache this?
    SmartPtr<const DenseVector> GlobalVector() const;

    /** Obtain pointer to global vector.  This is a copy of all data,
	changing it will not change data in the vector.
     */
    SmartPtr<DenseVector> GlobalVectorNonConst() const;

    /* Extract the local part from a global vector into the data of
       this vector. */
    void ExtractLocalVector(const DenseVector& global_vector);

    /* Extract the local part from a global values array into the data of
       this vector. */
    void ExtractLocalValues(const Number* vals);

    /** Rank of the processor handling this vector */
    int Rank() const;

    /** Total number of processors. */
    int NumProc() const;

    /** Offset of local vector */
    int StartPos() const;

    /** Size of local vector */
    int LocalSize() const;
    //@}

  protected:
    /** @name Overloaded methods from Vector base class */
    //@{
    /** Copy the data of the vector x into this vector (DCOPY). */
    virtual void CopyImpl(const Vector& x);

    /** Scales the vector by scalar alpha (DSCAL) */
    virtual void ScalImpl(Number alpha);

    /** Add the multiple alpha of vector x to this vector (DAXPY) */
    virtual void AxpyImpl(Number alpha, const Vector &x);

    /** Computes inner product of vector x with this (DDOT) */
    virtual Number DotImpl(const Vector &x) const;

    /** Computes the 2-norm of this vector (DNRM2) */
    virtual Number Nrm2Impl() const;

    /** Computes the 1-norm of this vector (DASUM) */
    virtual Number AsumImpl() const;

    /** Computes the max-norm of this vector (based on IDAMAX) */
    virtual Number AmaxImpl() const;

    /** Set each element in the vector to the scalar alpha. */
    virtual void SetImpl(Number value);

    /** Element-wise division  \f$y_i \gets y_i/x_i\f$.*/
    virtual void ElementWiseDivideImpl(const Vector& x);

    /** Element-wise multiplication \f$y_i \gets y_i*x_i\f$.*/
    virtual void ElementWiseMultiplyImpl(const Vector& x);

    /** Set entry to max of itself and the corresponding element in x */
    virtual void ElementWiseMaxImpl(const Vector& x);

    /** Set entry to min of itself and the corresponding element in x */
    virtual void ElementWiseMinImpl(const Vector& x);

    /** reciprocates the elements of the vector */
    virtual void ElementWiseReciprocalImpl();

    /** take abs of the elements of the vector */
    virtual void ElementWiseAbsImpl();

    /** take square-root of the elements of the vector */
    virtual void ElementWiseSqrtImpl();

    /** Changes each entry in the vector to its sgn value */
    virtual void ElementWiseSgnImpl();

    /** Add scalar to every component of the vector.*/
    virtual void AddScalarImpl(Number scalar);

    /** Max value in the vector */
    virtual Number MaxImpl() const;

    /** Min value in the vector */
    virtual Number MinImpl() const;

    /** Computes the sum of the lements of vector */
    virtual Number SumImpl() const;

    /** Computes the sum of the logs of the elements of vector */
    virtual Number SumLogsImpl() const;

    /** @name Implemented specialized functions */
    //@{
    /** Add two vectors (a * v1 + b * v2).  Result is stored in this
    vector. */
    void AddTwoVectorsImpl(Number a, const Vector& v1,
                           Number b, const Vector& v2, Number c);

    /** Fraction to the boundary parameter. */
    Number FracToBoundImpl(const Vector& delta, Number tau) const;

    /** Add the quotient of two vectors, y = a * z/s + c * y. */
    void AddVectorQuotientImpl(Number a, const Vector& z, const Vector& s,
                               Number c);
    //@}

    /** @name Output methods */
    //@{
    /* Print the entire vector with padding */
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
    ParVector();

    /** Copy Constructor */
    ParVector(const ParVector&);

    /** Overloaded Equals Operator */
    void operator=(const ParVector&);
    //@}

    /** Copy of the owner_space ptr as a ParVectorSpace instead
     *  of a VectorSpace
     */
    const ParVectorSpace* owner_space_;

    /** local piece of ParVector stored as a DenseVector */
    SmartPtr<DenseVector> local_vector_;
  };

  /** This vectors space is the vector space for ParVector.
   *  Before a ParVector can be created, all information relevant to
   *  a ParVector, such as dimension of the space, and the range of
   *  indices on the current processor have to be set in the ParVectorSpace.
   */
  class ParVectorSpace : public VectorSpace
  {
  public:
    /** @name Constructors/Destructors. */
    //@{
    /** Constructor, has to be given the total dimension of the vector
     *  and the local range of contiguous indices. */
    ParVectorSpace(Index total_dim, Index start_pos, Index size);

    /** Destructor */
    ~ParVectorSpace()
    {}
    //@}

    /** Method for creating a new vector of this specific type. */
    virtual ParVector* MakeNewParVector() const
    {
      return new ParVector(this);
    }

    /** Method for creating a new dense vector of global dimensions. */
    virtual DenseVector* MakeNewGlobalVector() const
    {
      return global_space_->MakeNewDenseVector();
    }

    /** Overloaded MakeNew method for the VectorSpace base class.
     */
    virtual Vector* MakeNew() const
    {
      return MakeNewParVector();
    }

    SmartPtr<DenseVectorSpace> LocalSpace() const
    {
      return local_space_;
    }

    SmartPtr<DenseVectorSpace> GlobalSpace() const
    {
      return global_space_;
    }

    int StartPos() const
    {
      return start_pos_;
    }

    const std::vector<int> & RecvCounts() const
    {
      return recvcounts_;
    }

    const std::vector<int> & Displs() const
    {
      return displs_;
    }

    int LocalSize() const
    {
      return size_;
    }

    int Rank() const
    {
      return rank_;
    }

    int NumProc() const
    {
      return num_proc_;
    }

  private:
    int start_pos_;
    int size_;
    int num_proc_;
    int rank_;
    std::vector<int> recvcounts_;
    std::vector<int> displs_;
    SmartPtr<DenseVectorSpace> local_space_;
    SmartPtr<DenseVectorSpace> global_space_;
  };

  inline
  SmartPtr<ParVector> ParVector::MakeNewParVector() const
  {
    return owner_space_->MakeNewParVector();
  }

  inline
  SmartPtr<DenseVector> ParVector::MakeNewGlobalVector() const
  {
    return owner_space_->MakeNewGlobalVector();
  }

  inline
  int ParVector::Rank() const
  {
    return owner_space_->Rank();
  }

  inline
  int ParVector::NumProc() const
  {
    return owner_space_->NumProc();
  }

  inline
  int ParVector::StartPos() const
  {
    return owner_space_->StartPos();
  }

  inline
  int ParVector::LocalSize() const
  {
    return owner_space_->LocalSize();
  }

} // namespace Ipopt

#endif
