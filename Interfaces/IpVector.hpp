// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPVECTOR_HPP__
#define __IPVECTOR_HPP__

#include "IpTypes.hpp"
#include "IpTaggedObject.hpp"
#include "IpCachedResults.hpp"
#include "IpSmartPtr.hpp"

#ifdef OLD_C_HEADERS
# include <stdio.h>
#else
# include <cstdio>
#endif
#include <string>

namespace Ipopt
{
  /* forward declarations */
  class VectorSpace;

  /** Vector Base Class.
   * This is the base class for all derived vector types.  Those vectors
   * are meant to store entities like iterates, Lagrangian multipliers,
   * constraint values etc.  The implementation of a vector type depends
   * on the computational environment (e.g. just a double array on a shared
   * memory machine, or distributed double arrays for a distributed
   * memory machine.)
   * 
   * Deriving from Vector: This class inherits from tagged object to
   * implement an advanced caching scheme. Because of this, the
   * TaggedObject method ObjectChanged() must be called each time the
   * Vector changes. If you overload the XXXX_Impl protected methods,
   * this taken care of (along with caching if possible) for you. If
   * you have additional methods in your derived class that change the
   * underlying data (vector values), you MUST remember to call
   * ObjectChanged() AFTER making the change!
   */
  class Vector : public TaggedObject
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    /** Constructor.  It has to be given a pointer to the
     *  corresponding VectorSpace.
     */
    Vector(const VectorSpace* owner_space);

    /** Destructor */
    virtual ~Vector();
    //@}

    /** Create new Vector of the same type with uninitialized data */
    Vector* MakeNew() const;

    /**@name Standard BLAS-1 Operations
     *  (derived classes do NOT overload these 
     *  methods, instead, overload the 
     *  protected versions of these methods). */
    //@{
    /** Copy the data of the vector x into this vector (DCOPY). */
    void Copy(const Vector& x);

    /** Scales the vector by scalar alpha (DSCAL) */
    void Scal(Number alpha);

    /** Add the multiple alpha of vector x to this vector (DAXPY) */
    void Axpy(Number alpha, const Vector &x);

    /** Computes inner product of vector x with this (DDOT) */
    Number Dot(const Vector &x) const;

    /** Computes the 2-norm of this vector (DNRM2) */
    Number Nrm2() const;

    /** Computes the 1-norm of this vector (DASUM) */
    Number Asum() const;

    /** Computes the max-norm of this vector (based on IDAMAX) */
    Number Amax() const;
    //@}

    /** @name Additional (Non-BLAS) Vector Methods
     *  (derived classes do NOT overload these 
     *  methods, instead, overload the 
     *  protected versions of these methods). */
    //@{
    /** Set each element in the vector to the scalar alpha. */
    void Set(Number alpha);

    /** Element-wise division  \f$y_i \gets y_i/x_i\f$*/
    void ElementWiseDivide(const Vector& x);

    /** Element-wise multiplication \f$y_i \gets y_i*x_i\f$ */
    void ElementWiseMultiply(const Vector& x);

    /** Element-wise max against entries in x */
    void ElementWiseMax(const Vector& x);

    /** Element-wise min against entries in x */
    void ElementWiseMin(const Vector& x);

    /** Reciprocates the entries in the vector */
    void ElementWiseReciprocal();

    /** Absolute values of the entries in the vector */
    void ElementWiseAbs();

    /** Element-wise square root of the entries in the vector */
    void ElementWiseSqrt();

    /** Replaces the vector values with their sgn values
    ( -1 if x_i < 0, 0 if x_i == 0, and 1 if x_i > 0)
    */
    void ElementWiseSgn();

    /** Add scalar to every vector component */
    void AddScalar(Number scalar);

    /** Returns the maximum value in the vector */
    Number Max() const;

    /** Returns the minimum value in the vector */
    Number Min() const;

    /** Returns the sum of the vector entries */
    Number Sum() const;

    /** Returns the sum of the logs of each vector entry */
    Number SumLogs() const;
    //@}

    /** @name Methods for specialized operations.  A prototype
     *  implementation is provided, but for efficient implementation
     *  those should be specially implemented.
     */
    //@{
    /** Add one vector, y = a * v1 + c * y.  This is automatically
     *  reduced to call AddTwoVectors.  */
    void AddOneVector(Number a, const Vector& v1, Number c);

    /** Add two vectors, y = a * v1 + b * v2 + c * y.  Here, this
	vector is y */
    void AddTwoVectors(Number a, const Vector& v1,
		       Number b, const Vector& v2, Number c);
    //@}

    /** @name Accessor methods */
    //@{
    /** Dimension of the Vector */
    Index Dim() const;

    /** Return the owner VectorSpace*/
    SmartPtr<const VectorSpace> OwnerSpace() const;
    //@}

    /** @name Output methods
     *  (derived classes do NOT overload these 
     *  methods, instead, overload the 
     *  protected versions of these methods). */
    //@{
    /** Print the entire vector */
    void Print(FILE* fp, std::string name,
               Index indent, std::string prefix) const;
    //@}

  protected:
    /** @name implementation methods (derived classes MUST
     *  overload these pure virtual protected methods.)
     */
    //@{
    /** Copy the data of the vector x into this vector (DCOPY). */
    virtual void CopyImpl(const Vector& x)=0;

    /** Scales the vector by scalar alpha (DSCAL) */
    virtual void ScalImpl(Number alpha)=0;

    /** Add the multiple alpha of vector x to this vector (DAXPY) */
    virtual void AxpyImpl(Number alpha, const Vector &x)=0;

    /** Computes inner product of vector x with this (DDOT) */
    virtual Number DotImpl(const Vector &x) const =0;

    /** Computes the 2-norm of this vector (DNRM2) */
    virtual Number Nrm2Impl() const =0;

    /** Computes the 1-norm of this vector (DASUM) */
    virtual Number AsumImpl() const =0;

    /** Computes the max-norm of this vector (based on IDAMAX) */
    virtual Number AmaxImpl() const =0;

    /** Set each element in the vector to the scalar alpha. */
    virtual void SetImpl(Number alpha)=0;

    /** Element-wise division  \f$y_i \gets y_i/x_i\f$*/
    virtual void ElementWiseDivideImpl(const Vector& x)=0;

    /** Element-wise multiplication \f$y_i \gets y_i*x_i\f$ */
    virtual void ElementWiseMultiplyImpl(const Vector& x)=0;

    /** Element-wise max against entries in x */
    virtual void ElementWiseMaxImpl(const Vector& x)=0;

    /** Element-wise min against entries in x */
    virtual void ElementWiseMinImpl(const Vector& x)=0;

    /** Reciprocates the elements of the vector */
    virtual void ElementWiseReciprocalImpl()=0;

    /** Take elementwise absolute values of the elements of the vector */
    virtual void ElementWiseAbsImpl()=0;

    /** Take elementwise square-root of the elements of the vector */
    virtual void ElementWiseSqrtImpl()=0;

    /** Replaces entries with sgn of the entry */
    virtual void ElementWiseSgnImpl()=0;

    /** Add scalar to every component of vector */
    virtual void AddScalarImpl(Number scalar)=0;

    /** Max value in the vector */
    virtual Number MaxImpl() const=0;

    /** Min number in the vector */
    virtual Number MinImpl() const=0;

    /** Sum of entries in the vector */
    virtual Number SumImpl() const=0;

    /** Sum of logs of entries in the vector */
    virtual Number SumLogsImpl() const=0;

    /** Add two vectors (a * v1 + b * v2).  Result is stored in this
	vector. */
    virtual void AddTwoVectorsImpl(Number a, const Vector& v1,
				   Number b, const Vector& v2, Number c);

    /** Print the entire vector */
    virtual void PrintImpl(FILE* fp, std::string name = "Vector",
                           Index indent=0, std::string prefix = "") const =0;
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
    /** Default constructor */
    Vector();

    /** Copy constructor */
    Vector(const Vector&);

    /** Overloaded Equals Operator */
    Vector& operator=(const Vector&);
    //@}

    /** Vector Space */
    const SmartPtr<const VectorSpace> owner_space_;

    /**@name CachedResults data members */
    //@{
    /** Cache for dot products */
    mutable CachedResults<Number> dot_cache_;
    /** Cache for 2-norm */
    mutable CachedResults<Number> nrm2_cache_;
    /** Cache for Asum */
    mutable CachedResults<Number> asum_cache_;
    /** Cache for Amax */
    mutable CachedResults<Number> amax_cache_;
    /** Cache for Sum */
    mutable CachedResults<Number> sum_cache_;
    /** Cache for SumLogs */
    mutable CachedResults<Number> sumlogs_cache_;
    //@}

  };

  /** VectorSpace base class, corresponding to the Vector base class.
   *  For each Vector implementation, a corresponding VectorSpace has
   *  to be implemented.  A VectorSpace is able to create new Vectors
   *  of a specific type.  The VectorSpace should also store
   *  information that is common to all Vectors of that type.  For
   *  example, the dimension of a Vector is stored in the VectorSpace
   *  base class.
   */
  class VectorSpace : public ReferencedObject
  {
  public:
    /** @name Constructors/Destructors */
    //@{
    /** Constructor, given the dimension of all vectors generated by
     *  this VectorSpace.
     */
    VectorSpace(Index dim);

    /** Destructor */
    virtual ~VectorSpace();
    //@}

    /** Pure virtual method for creating a new Vector of the
     *  corresponding type.
     */
    virtual Vector* MakeNew() const=0;

    /** Accessor function for the dimension of the vectors of this type.*/
    Index Dim() const
    {
      return dim_;
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
    /** default constructor */
    VectorSpace();

    /** Copy constructor */
    VectorSpace(const VectorSpace&);

    /** Overloaded Equals Operator */
    VectorSpace& operator=(const VectorSpace&);
    //@}

    /** Dimension of the vectors in this vector space. */
    const Index dim_;
  };

  /* inline methods */
  inline
  Vector* Vector::MakeNew() const
  {
    return owner_space_->MakeNew();
  }

  inline
  void Vector::Copy(const Vector& x)
  {
    CopyImpl(x);
    ObjectChanged();
  }

  inline
  void Vector::Scal(Number alpha)
  {
    if (alpha!=1.) {
      ScalImpl(alpha);
      ObjectChanged();
    }
  }

  inline
  void Vector::Axpy(Number alpha, const Vector &x)
  {
    AxpyImpl(alpha, x);
    ObjectChanged();
  }

  inline
  Number Vector::Dot(const Vector &x) const
  {
    Number retValue = 0.0;
    if (!dot_cache_.GetCachedResult2Dep(retValue, this, &x)) {
      retValue = DotImpl(x);
      dot_cache_.AddCachedResult2Dep(retValue, this, &x);
    }
    return retValue;
  }

  inline
  Number Vector::Nrm2() const
  {
    Number retValue = 0.0;
    if (!nrm2_cache_.GetCachedResult1Dep(retValue, this)) {
      retValue = Nrm2Impl();
      nrm2_cache_.AddCachedResult1Dep(retValue, this);
    }
    return retValue;
  }

  inline
  Number Vector::Asum() const
  {
    Number retValue = 0.0;
    if (!asum_cache_.GetCachedResult1Dep(retValue, this)) {
      retValue = AsumImpl();
      asum_cache_.AddCachedResult1Dep(retValue, this);
    }
    return retValue;
  }

  inline
  Number Vector::Amax() const
  {
    Number retValue = 0.0;
    if (!amax_cache_.GetCachedResult1Dep(retValue, this)) {
      retValue = AmaxImpl();
      amax_cache_.AddCachedResult1Dep(retValue, this);
    }
    return retValue;
  }

  inline
  Number Vector::Sum() const
  {
    Number retValue = 0.0;
    if (!sum_cache_.GetCachedResult1Dep(retValue, this)) {
      retValue = SumImpl();
      sum_cache_.AddCachedResult1Dep(retValue, this);
    }
    return retValue;
  }

  inline
  Number Vector::SumLogs() const
  {
    Number retValue = 0.0;
    if (!sumlogs_cache_.GetCachedResult1Dep(retValue, this)) {
      retValue = SumLogsImpl();
      sumlogs_cache_.AddCachedResult1Dep(retValue, this);
    }
    return retValue;
  }

  inline
  void Vector::ElementWiseSgn()
  {
    ElementWiseSgnImpl();
    ObjectChanged();
  }

  inline
  void Vector::Set(Number alpha)
  {
    SetImpl(alpha);
    ObjectChanged();
  }

  inline
  void Vector::ElementWiseDivide(const Vector& x)
  {
    ElementWiseDivideImpl(x);
    ObjectChanged();
  }

  inline
  void Vector::ElementWiseMultiply(const Vector& x)
  {
    ElementWiseMultiplyImpl(x);
    ObjectChanged();
  }

  inline
  void Vector::ElementWiseReciprocal()
  {
    ElementWiseReciprocalImpl();
    ObjectChanged();
  }

  inline
  void Vector::ElementWiseMax(const Vector& x)
  {
    ElementWiseMaxImpl(x);
    ObjectChanged();
  }

  inline
  void Vector::ElementWiseMin(const Vector& x)
  {
    ElementWiseMinImpl(x);
    ObjectChanged();
  }

  inline
  void Vector::ElementWiseAbs()
  {
    ElementWiseAbsImpl();
    ObjectChanged();
  }

  inline
  void Vector::ElementWiseSqrt()
  {
    ElementWiseSqrtImpl();
    ObjectChanged();
  }

  inline
  void Vector::AddScalar(Number scalar)
  {
    AddScalarImpl(scalar);
    ObjectChanged();
  }

  inline
  Number Vector::Max() const
  {
    return MaxImpl();
  }

  inline
  Number Vector::Min() const
  {
    return MinImpl();
  }

  inline
  void Vector::AddOneVector(Number a, const Vector& v1, Number c)
  {
    AddTwoVectorsImpl(a, v1, 0., v1, c);
  }

  inline
  void Vector::AddTwoVectors(Number a, const Vector& v1,
			     Number b, const Vector& v2, Number c)
  {
    AddTwoVectorsImpl(a, v1, b, v2, c);
    ObjectChanged();
  }

  inline
  Index Vector::Dim() const
  {
    return owner_space_->Dim();
  }

  inline
  void Vector::Print(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    PrintImpl(fp, name, indent, prefix);
  }

  inline
  SmartPtr<const VectorSpace> Vector::OwnerSpace() const
  {
    return owner_space_;
  }

} // namespace Ipopt


#endif
