// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPCOMPOUNDVECTOR_HPP__
#define __IPCOMPOUNDVECTOR_HPP__

#include "IpUtils.hpp"
#include "IpVector.hpp"
#include <vector>

namespace Ipopt
{

/* forward declarations */
class CompoundVectorSpace;

/** Class of Vectors consisting of other vectors.
 *
 *  This vector is a
 *  vector that consists of zero, one or more Vector's which are
 *  stacked on each others: \f$ x_{\rm compound} =
 *  \left(\begin{array}{c}x_0\\\dots\\x_{{\rm
 *  ncomps} - 1}\end{array}\right)\f$.  The individual components can be
 *  associated to different VectorSpaces.  The individual components
 *  can also be const and non-const Vectors.
 */
class IPOPTLIB_EXPORT CompoundVector: public Vector
{
public:
   /**@name Constructors/Destructors */
   ///@{
   /** Constructor, given the corresponding CompoundVectorSpace.
    *
    *  Before this constructor can be called, all components of the
    *  CompoundVectorSpace have to be set, so that the constructors
    *  for the individual components can be called.  If the flag
    *  create_new is true, then the individual components of the new
    *  CompoundVector are initialized with the MakeNew methods of
    *  each VectorSpace (and are non-const).  Otherwise, the
    *  individual components can later be set using the SetComp and
    *  SetCompNonConst method.
    */
   CompoundVector(
      const CompoundVectorSpace* owner_space,
      bool                       create_new
   );

   /** Default destructor */
   virtual ~CompoundVector() { }
   ///@}

   /** Method for setting the pointer for a component that is a const Vector */
   void SetComp(
      Index         icomp,
      const Vector& vec
   );

   /** Method for setting the pointer for a component that is a non-const Vector */
   void SetCompNonConst(
      Index   icomp,
      Vector& vec
   );

   /** Number of components of this compound vector */
   inline Index NComps() const;

   /** Check if a particular component is const or not */
   bool IsCompConst(
      Index i
   ) const
   {
      DBG_ASSERT(i > 0 && i < NComps());
      DBG_ASSERT(IsValid(comps_[i]) || IsValid(const_comps_[i]));
      if( IsValid(const_comps_[i]) )
      {
         return true;
      }
      return false;
   }

   /** Check if a particular component is null or not */
   bool IsCompNull(
      Index i
   ) const
   {
      DBG_ASSERT(i >= 0 && i < NComps());
      if( IsValid(comps_[i]) || IsValid(const_comps_[i]) )
      {
         return false;
      }
      return true;
   }

   /** Return a particular component (const version) */
   SmartPtr<const Vector> GetComp(
      Index i
   ) const
   {
      return ConstComp(i);
   }

   /** Return a particular component (non-const version).
    *
    *  Note that calling this method with mark the CompoundVector as changed.
    *  Therefore, only use this method if you are intending to change
    *  the Vector that you receive.
    */
   SmartPtr<Vector> GetCompNonConst(
      Index i
   )
   {
      ObjectChanged();
      return Comp(i);
   }

protected:
   /** @name Overloaded methods from Vector base class */
   ///@{
   virtual void CopyImpl(
      const Vector& x
   );

   virtual void ScalImpl(
      Number alpha
   );

   virtual void AxpyImpl(
      Number        alpha,
      const Vector& x
   );

   virtual Number DotImpl(
      const Vector& x
   ) const;

   virtual Number Nrm2Impl() const;

   virtual Number AsumImpl() const;

   virtual Number AmaxImpl() const;

   virtual void SetImpl(
      Number value
   );

   virtual void ElementWiseDivideImpl(
      const Vector& x
   );

   virtual void ElementWiseMultiplyImpl(
      const Vector& x
   );

   virtual void ElementWiseSelectImpl(
      const Vector& x
   );

   virtual void ElementWiseMaxImpl(
      const Vector& x
   );

   virtual void ElementWiseMinImpl(
      const Vector& x
   );

   virtual void ElementWiseReciprocalImpl();

   virtual void ElementWiseAbsImpl();

   virtual void ElementWiseSqrtImpl();

   virtual void ElementWiseSgnImpl();

   virtual void AddScalarImpl(
      Number scalar
   );

   virtual Number MaxImpl() const;

   virtual Number MinImpl() const;

   virtual Number SumImpl() const;

   virtual Number SumLogsImpl() const;
   ///@}

   /** @name Implemented specialized functions */
   ///@{
   void AddTwoVectorsImpl(
      Number        a,
      const Vector& v1,
      Number        b,
      const Vector& v2,
      Number        c
   );

   Number FracToBoundImpl(
      const Vector& delta,
      Number        tau
   ) const;

   void AddVectorQuotientImpl(
      Number        a,
      const Vector& z,
      const Vector& s,
      Number        c
   );
   ///@}

   /** Method for determining if all stored numbers are valid (i.e., no Inf or Nan). */
   virtual bool HasValidNumbersImpl() const;

   /** @name Output methods */
   ///@{
   virtual void PrintImpl(
      const Journalist&  jnlst,
      EJournalLevel      level,
      EJournalCategory   category,
      const std::string& name,
      Index              indent,
      const std::string& prefix
   ) const;
   ///@}

private:
   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called.
    */
   ///@{
   /** Default Constructor */
   CompoundVector();

   /** Copy Constructor */
   CompoundVector(
      const CompoundVector&
   );

   /** Default Assignment Operator */
   void operator=(
      const CompoundVector&
   );
   ///@}

   /** Components of the compound vector.
    *
    * The components are stored by SmartPtrs in a std::vector
    */
   std::vector<SmartPtr<Vector> > comps_;
   std::vector<SmartPtr<const Vector> > const_comps_;

   const CompoundVectorSpace* owner_space_;

   bool vectors_valid_;

   bool VectorsValid();

   inline const Vector* ConstComp(
      Index i
   ) const;

   inline Vector* Comp(
      Index i
   );
};

/** This vectors space is the vector space for CompoundVector.
 *
 *  Before a CompoundVector can be created, all components of this
 *  CompoundVectorSpace have to be set.  When calling the constructor,
 *  the number of component has to be specified.  The individual
 *  VectorSpaces can be set with the SetComp method.
 */
class IPOPTLIB_EXPORT CompoundVectorSpace: public VectorSpace
{
public:
   /** @name Constructors/Destructors. */
   ///@{
   /** Constructor, has to be given the number of components and the
    *  total dimension of all components combined.
    */
   CompoundVectorSpace(
      Index ncomp_spaces,
      Index total_dim
   );

   /** Destructor */
   ~CompoundVectorSpace()
   { }
   ///@}

   /** Method for setting the individual component VectorSpaces */
   virtual void SetCompSpace(
      Index              icomp,    /**< Number of the component to be set */
      const VectorSpace& vec_space /**< VectorSpace for component icomp */
   );

   /** Method for obtaining an individual component VectorSpace */
   SmartPtr<const VectorSpace> GetCompSpace(
      Index icomp
   ) const;

   /** Accessor method to obtain the number of components */
   Index NCompSpaces() const
   {
      return ncomp_spaces_;
   }

   /** Method for creating a new vector of this specific type. */
   virtual CompoundVector* MakeNewCompoundVector(
      bool create_new = true
   ) const
   {
      return new CompoundVector(this, create_new);
   }

   virtual Vector* MakeNew() const
   {
      return MakeNewCompoundVector();
   }

private:
   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called. */
   ///@{
   /** Default constructor */
   CompoundVectorSpace();

   /** Copy Constructor */
   CompoundVectorSpace(
      const CompoundVectorSpace&
   );

   /** Default Assignment Operator */
   CompoundVectorSpace& operator=(
      const CompoundVectorSpace&
   );
   ///@}

   /** Number of components */
   const Index ncomp_spaces_;

   /** std::vector of vector spaces for the components */
   std::vector<SmartPtr<const VectorSpace> > comp_spaces_;
};

/* inline methods */
inline Index CompoundVector::NComps() const
{
   return owner_space_->NCompSpaces();
}

inline const Vector* CompoundVector::ConstComp(
   Index i
) const
{
   DBG_ASSERT(i < NComps());
   DBG_ASSERT(IsValid(comps_[i]) || IsValid(const_comps_[i]));
   if( IsValid(comps_[i]) )
   {
      return GetRawPtr(comps_[i]);
   }
   else if( IsValid(const_comps_[i]) )
   {
      return GetRawPtr(const_comps_[i]);
   }

   DBG_ASSERT(false && "shouldn't be here");
   return NULL;
}

inline Vector* CompoundVector::Comp(
   Index i
)
{
   DBG_ASSERT(i < NComps());
   DBG_ASSERT(IsValid(comps_[i]));
   return GetRawPtr(comps_[i]);
}

} // namespace Ipopt

#endif
