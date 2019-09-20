// Copyright (C) 2004, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPDIAGMATRIX_HPP__
#define __IPDIAGMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpSymMatrix.hpp"

namespace Ipopt
{

/** Class for diagonal matrices.
 *
 *  The diagonal is stored as a Vector.
 */
class IPOPTLIB_EXPORT DiagMatrix: public SymMatrix
{
public:

   /**@name Constructors / Destructors */
   //@{
   /** Constructor, given the corresponding matrix space. */
   DiagMatrix(
      const SymMatrixSpace* owner_space
   );

   /** Destructor */
   ~DiagMatrix();
   //@}

   /** Set the diagonal elements (as a Vector). */
   void SetDiag(
      const Vector& diag
   )
   {
      diag_ = &diag;
   }

   /** Get the diagonal elements. */
   SmartPtr<const Vector> GetDiag() const
   {
      return diag_;
   }

protected:
   /**@name Methods overloaded from matrix */
   //@{
   virtual void MultVectorImpl(
      Number        alpha,
      const Vector& x,
      Number        beta,
      Vector&       y
   ) const;

   virtual bool HasValidNumbersImpl() const;

   virtual void ComputeRowAMaxImpl(
      Vector& rows_norms,
      bool    init
   ) const;

   virtual void PrintImpl(
      const Journalist&  jnlst,
      EJournalLevel      level,
      EJournalCategory   category,
      const std::string& name,
      Index              indent,
      const std::string& prefix
   ) const;
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
   DiagMatrix();

   /** Copy Constructor */
   DiagMatrix(
      const DiagMatrix&
   );

   /** Default Assignment Operator */
   void operator=(
      const DiagMatrix&
   );
   //@}

   /** Vector storing the diagonal elements */
   SmartPtr<const Vector> diag_;
};

/** This is the matrix space for DiagMatrix. */
class IPOPTLIB_EXPORT DiagMatrixSpace: public SymMatrixSpace
{
public:
   /** @name Constructors / Destructors */
   //@{
   /** Constructor, given the dimension of the matrix. */
   DiagMatrixSpace(
      Index dim
   )
      : SymMatrixSpace(dim)
   { }

   /** Destructor */
   virtual ~DiagMatrixSpace()
   { }
   //@}

   virtual SymMatrix* MakeNewSymMatrix() const
   {
      return MakeNewDiagMatrix();
   }

   /** Method for creating a new matrix of this specific type. */
   DiagMatrix* MakeNewDiagMatrix() const
   {
      return new DiagMatrix(this);
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
   DiagMatrixSpace();

   /** Copy Constructor */
   DiagMatrixSpace(
      const DiagMatrixSpace&
   );

   /** Default Assignment Operator */
   void operator=(
      const DiagMatrixSpace&
   );
   //@}

};

} // namespace Ipopt
#endif
