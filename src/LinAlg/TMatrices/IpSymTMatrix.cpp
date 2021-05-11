// Copyright (C) 2004, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpSymTMatrix.hpp"
#include "IpDenseVector.hpp"
#include "IpBlas.hpp"

#include <cmath>

namespace Ipopt
{

SymTMatrix::SymTMatrix(
   const SymTMatrixSpace* owner_space
)
   : SymMatrix(owner_space),
     owner_space_(owner_space),
     values_(NULL),
     initialized_(false)
{
   values_ = owner_space_->AllocateInternalStorage();

   if( Nonzeros() == 0 )
   {
      initialized_ = true; // I guess ?!? what does this mean ?!?
   }
}

SymTMatrix::~SymTMatrix()
{
   owner_space_->FreeInternalStorage(values_);
}

void SymTMatrix::SetValues(
   const Number* Values
)
{
   IpBlasCopy(Nonzeros(), Values, 1, values_, 1);
   initialized_ = true;
   ObjectChanged();
}

void SymTMatrix::MultVectorImpl(
   Number        alpha,
   const Vector& x,
   Number        beta,
   Vector&       y
) const
{
   //  A few sanity checks
   DBG_ASSERT(Dim() == x.Dim());
   DBG_ASSERT(Dim() == y.Dim());

   // Take care of the y part of the addition
   DBG_ASSERT(initialized_);
   if( beta != 0.0 )
   {
      y.Scal(beta);
   }
   else
   {
      y.Set(0.0);  // In case y hasn't been initialized yet
   }

   // See if we can understand the data
   const DenseVector* dense_x = static_cast<const DenseVector*>(&x);
   DBG_ASSERT(dynamic_cast<const DenseVector*>(&x));
   DenseVector* dense_y = static_cast<DenseVector*>(&y);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&y));

   if( dense_x && dense_y )
   {
      const Index* irn = Irows();
      const Index* jcn = Jcols();
      const Number* val = values_;
      Number* yvals = dense_y->Values();

      if( dense_x->IsHomogeneous() )
      {
         Number as = alpha * dense_x->Scalar();
         for( Index i = 0; i < Nonzeros(); i++ )
         {
            yvals[*irn - 1] += as * (*val);
            if( *irn != *jcn )
            {
               // this is not a diagonal element
               yvals[*jcn - 1] += as * (*val);
            }
            val++;
            irn++;
            jcn++;
         }
      }
      else
      {
         const Number* xvals = dense_x->Values();
         for( Index i = 0; i < Nonzeros(); i++ )
         {
            yvals[*irn - 1] += alpha * (*val) * xvals[*jcn - 1];
            if( *irn != *jcn )
            {
               // this is not a diagonal element
               yvals[*jcn - 1] += alpha * (*val) * xvals[*irn - 1];
            }
            val++;
            irn++;
            jcn++;
         }
      }
   }
}

Number* SymTMatrix::Values()
{
   // cannot check for initialized values here, in case this pointer is
   // requested for setting the first values
   //DBG_ASSERT(initialized_);

   // Here we assume that every time someone requests this direct raw
   // pointer, the data is going to change and the Tag for this
   // vector has to be updated.
   ObjectChanged();
   initialized_ = true;
   return values_;
}

const Number* SymTMatrix::Values() const
{
   DBG_ASSERT(initialized_);
   return values_;
}

void SymTMatrix::FillStruct(
   Index* Irn,
   Index* Jcn
) const
{
   DBG_ASSERT(initialized_);
   for( Index i = 0; i < Nonzeros(); i++ )
   {
      Irn[i] = Irows()[i];
      Jcn[i] = Jcols()[i];
   }
}

void SymTMatrix::FillValues(
   Number* Values
) const
{
   DBG_ASSERT(initialized_);
   IpBlasCopy(Nonzeros(), values_, 1, Values, 1);
}

bool SymTMatrix::HasValidNumbersImpl() const
{
   DBG_ASSERT(initialized_);
   Number sum = IpBlasAsum(Nonzeros(), values_, 1);
   return IsFiniteNumber(sum);
}

void SymTMatrix::ComputeRowAMaxImpl(
   Vector& rows_norms,
   bool    /*init*/
) const
{
   DBG_ASSERT(initialized_);

   if( NRows() == 0 )
   {
      return;
   }

   DenseVector* dense_vec = static_cast<DenseVector*>(&rows_norms);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&rows_norms));

   const Index* irn = Irows();
   const Index* jcn = Jcols();
   const Number* val = values_;
   Number* vec_vals = dense_vec->Values();
   DBG_ASSERT(vec_vals != NULL);

   const Number zero = 0.;
   IpBlasCopy(NRows(), &zero, 0, vec_vals, 1);

   vec_vals--; // to deal with 1-based indexing in irn and jcn (I believe)
   for( Index i = 0; i < Nonzeros(); i++ )
   {
      const Number f = std::abs(*val);
      vec_vals[*irn] = Max(vec_vals[*irn], f);
      vec_vals[*jcn] = Max(vec_vals[*jcn], f);
      val++;
      irn++;
      jcn++;
   }
}

void SymTMatrix::PrintImpl(
   const Journalist&  jnlst,
   EJournalLevel      level,
   EJournalCategory   category,
   const std::string& name,
   Index              indent,
   const std::string& prefix
) const
{
   jnlst.Printf(level, category,
                "\n");
   jnlst.PrintfIndented(level, category, indent,
                        "%sSymTMatrix \"%s\" of dimension %" IPOPT_INDEX_FORMAT " with %" IPOPT_INDEX_FORMAT " nonzero elements:\n", prefix.c_str(), name.c_str(), Dim(), Nonzeros());
   if( initialized_ )
   {
      for( Index i = 0; i < Nonzeros(); i++ )
      {
         jnlst.PrintfIndented(level, category, indent,
                              "%s%s[%5" IPOPT_INDEX_FORMAT ",%5" IPOPT_INDEX_FORMAT "]=%23.16e  (%" IPOPT_INDEX_FORMAT ")\n", prefix.c_str(), name.c_str(), Irows()[i], Jcols()[i], values_[i], i);
      }
   }
   else
   {
      jnlst.PrintfIndented(level, category, indent,
                           "%sUninitialized!\n", prefix.c_str());
   }
}

SymTMatrixSpace::SymTMatrixSpace(
   Index        dim,
   Index        nonZeros,
   const Index* iRows,
   const Index* jCols
)
   : SymMatrixSpace(dim),
     nonZeros_(nonZeros),
     iRows_(NULL),
     jCols_(NULL)
{
   iRows_ = new Index[nonZeros];
   jCols_ = new Index[nonZeros];
   for( Index i = 0; i < nonZeros; i++ )
   {
      iRows_[i] = iRows[i];
      jCols_[i] = jCols[i];
   }
}

SymTMatrixSpace::~SymTMatrixSpace()
{
   delete[] iRows_;
   delete[] jCols_;
}

Number* SymTMatrixSpace::AllocateInternalStorage() const
{
   return new Number[Nonzeros()];
}

void SymTMatrixSpace::FreeInternalStorage(
   Number* values) const
{
   delete[] values;
}

} // namespace Ipopt
