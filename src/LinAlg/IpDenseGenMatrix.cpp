// Copyright (C) 2005, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter             IBM    2005-12-24

#include "IpDenseGenMatrix.hpp"
#include "IpBlas.hpp"
#include "IpLapack.hpp"

#include <cmath>

namespace Ipopt
{

#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
#endif

DenseGenMatrix::DenseGenMatrix(
   const DenseGenMatrixSpace* owner_space
)
   : Matrix(owner_space),
     owner_space_(owner_space),
     values_(new Number[NCols() * NRows()]),
     initialized_(false),
     factorization_(NONE),
     pivot_(NULL)
{
}

DenseGenMatrix::~DenseGenMatrix()
{
   DBG_START_METH("DenseGenMatrix::~DenseGenMatrix()", dbg_verbosity);
   delete[] values_;
   delete[] pivot_;
}

void DenseGenMatrix::ScaleColumns(
   const DenseVector& scal_vec
)
{
   DBG_ASSERT(scal_vec.Dim() == NCols());
   DBG_ASSERT(initialized_);

   const Number* scal_values = scal_vec.Values();

   for( Index j = 0; j < NCols(); j++ )
   {
      IpBlasScal(NRows(), scal_values[j], &values_[j * NRows()], 1);
   }
   ObjectChanged();
}

void DenseGenMatrix::Copy(
   const DenseGenMatrix& M
)
{
   DBG_ASSERT(NCols() == M.NCols());
   DBG_ASSERT(NRows() == M.NRows());

   IpBlasCopy(NCols() * NRows(), M.Values(), 1, values_, 1);
   initialized_ = true;
   ObjectChanged();
}

void DenseGenMatrix::FillIdentity(
   Number factor /*=1.*/
)
{
   DBG_ASSERT(NCols() == NRows());

   const Number zero = 0.;
   IpBlasCopy(NCols() * NRows(), &zero, 0, values_, 1);

   if( factor != 0. )
   {
      for( Index i = 0; i < NRows(); i++ )
      {
         values_[i + i * NRows()] = factor;
      }
   }
   ObjectChanged();
   initialized_ = true;
}

void DenseGenMatrix::AddMatrixProduct(
   Number                alpha,
   const DenseGenMatrix& A,
   bool                  transA,
   const DenseGenMatrix& B,
   bool                  transB,
   Number                beta
)
{
   Index m = NRows();
   DBG_ASSERT((transA && A.NCols() == m) || (!transA && A.NRows() == m));
   Index n = NCols();
   DBG_ASSERT((transB && B.NRows() == n) || (!transB && B.NCols() == n));
   Index k;
   if( transA )
   {
      k = A.NRows();
   }
   else
   {
      k = A.NCols();
   }
   DBG_ASSERT((transB && B.NCols() == k) || (!transB && B.NRows() == k));
   DBG_ASSERT(beta == 0. || initialized_);

   IpBlasGemm(transA, transB, m, n, k, alpha, A.Values(), A.NRows(), B.Values(), B.NRows(), beta, values_, NRows());
   initialized_ = true;
   ObjectChanged();
}

void DenseGenMatrix::HighRankUpdateTranspose(
   Number                   alpha,
   const MultiVectorMatrix& V1,
   const MultiVectorMatrix& V2,
   Number                   beta
)
{
   DBG_ASSERT(NRows() == V1.NCols());
   DBG_ASSERT(NCols() == V2.NCols());
   DBG_ASSERT(beta == 0. || initialized_);

   if( beta == 0. )
   {
      for( Index j = 0; j < NCols(); j++ )
      {
         for( Index i = 0; i < NRows(); i++ )
         {
            values_[i + j * NRows()] = alpha * V1.GetVector(i)->Dot(*V2.GetVector(j));
         }
      }
   }
   else
   {
      for( Index j = 0; j < NCols(); j++ )
      {
         for( Index i = 0; i < NRows(); i++ )
         {
            values_[i + j * NRows()] = alpha * V1.GetVector(i)->Dot(*V2.GetVector(j)) + beta * values_[i + j * NRows()];
         }
      }
   }
   initialized_ = true;
   ObjectChanged();
}

bool DenseGenMatrix::ComputeCholeskyFactor(
   const DenseSymMatrix& M
)
{
   Index dim = M.Dim();
   DBG_ASSERT(dim == NCols());
   DBG_ASSERT(dim == NRows());

   ObjectChanged();

   // First we copy the content of the symmetric matrix into J
   const Number* Mvalues = M.Values();
   for( Index j = 0; j < dim; j++ )
   {
      for( Index i = j; i < dim; i++ )
      {
         values_[i + j * dim] = Mvalues[i + j * dim];
      }
   }

   // Now call the lapack subroutine to perform the factorization
   Index info;
   IpLapackPotrf(dim, values_, dim, info);

   DBG_ASSERT(info >= 0);
   if( info != 0 )
   {
      initialized_ = false;
      return false;
   }

   // We set all strictly upper values to zero
   // ToDo: This might not be necessary?!?
   for( Index j = 1; j < dim; j++ )
   {
      for( Index i = 0; i < j; i++ )
      {
         values_[i + j * dim] = 0.;
      }
   }

   factorization_ = CHOL;
   initialized_ = true;
   return true;
}

bool DenseGenMatrix::ComputeEigenVectors(
   const DenseSymMatrix& M,
   DenseVector&          Evalues
)
{
   Index dim = M.Dim();
   DBG_ASSERT(Evalues.Dim() == dim);
   DBG_ASSERT(NRows() == dim);
   DBG_ASSERT(NCols() == dim);

   // First we copy the content of the matrix into Q
   const Number* Mvalues = M.Values();
   for( Index j = 0; j < dim; j++ )
   {
      for( Index i = j; i < dim; i++ )
      {
         values_[i + j * dim] = Mvalues[i + j * dim];
      }
   }

   bool compute_eigenvectors = true;
   Number* Evals = Evalues.Values();
   Index info;
   IpLapackSyev(compute_eigenvectors, dim, values_, dim, Evals, info);

   initialized_ = (info == 0);
   ObjectChanged();
   return (info == 0);
}

void DenseGenMatrix::CholeskyBackSolveMatrix(
   bool            trans,
   Number          alpha,
   // cppcheck-suppress constParameter  // cannot be const since it is modified in Blas
   DenseGenMatrix& B
) const
{
   DBG_ASSERT(NRows() == NCols());
   DBG_ASSERT(B.NRows() == NRows());
   DBG_ASSERT(initialized_);

   IpBlasTrsm(trans, NRows(), B.NCols(), alpha, values_, NRows(), B.Values(), B.NRows());
}

void DenseGenMatrix::CholeskySolveVector(
   // cppcheck-suppress constParameter  // cannot be const since it is modified in Blas
   DenseVector& b
) const
{
   DBG_ASSERT(NRows() == NCols());
   DBG_ASSERT(b.Dim() == NRows());
   DBG_ASSERT(initialized_);
   DBG_ASSERT(factorization_ == CHOL);

   IpLapackPotrs(NRows(), 1, values_, NRows(), b.Values(), b.Dim());
}

void DenseGenMatrix::CholeskySolveMatrix(
   // cppcheck-suppress constParameter  // cannot be const since it is modified in Blas
   DenseGenMatrix& B
) const
{
   DBG_ASSERT(NRows() == NCols());
   DBG_ASSERT(B.NRows() == NRows());
   DBG_ASSERT(initialized_);
   DBG_ASSERT(factorization_ == CHOL);

   IpLapackPotrs(NRows(), B.NCols(), values_, NRows(), B.Values(), B.NRows());
}

bool DenseGenMatrix::ComputeLUFactorInPlace()
{
   Index dim = NRows();
   DBG_ASSERT(dim == NCols());

   DBG_ASSERT(factorization_ == NONE);

   ObjectChanged();

   // create pivot space
   delete[] pivot_;
   pivot_ = NULL; // set to NULL so that destructor will not try to
   // delete again if the new in following line fails
   pivot_ = new Index[dim];

   // call the lapack subroutine for the factorization (dgetrf )
   Index info;
   IpLapackGetrf(dim, values_, pivot_, dim, info);

   DBG_ASSERT(info >= 0);
   if( info != 0 )
   {
      delete[] pivot_;
      pivot_ = NULL;
      initialized_ = false;
      return false;
   }
   else
   {
      initialized_ = true;
   }

   factorization_ = LU;
   return true;
}

void DenseGenMatrix::LUSolveMatrix(
   // cppcheck-suppress constParameter  // cannot be const since it is modified in Blas
   DenseGenMatrix& B
) const
{
   DBG_ASSERT(NRows() == NCols());
   DBG_ASSERT(B.NRows() == NRows());
   DBG_ASSERT(initialized_);
   DBG_ASSERT(factorization_ == LU);

   IpLapackGetrs(NRows(), B.NCols(), values_, NRows(), pivot_, B.Values(), B.NRows());
}

void DenseGenMatrix::LUSolveVector(
   // cppcheck-suppress constParameter  // cannot be const since it is modified in Blas
   DenseVector& b
) const
{
   DBG_ASSERT(NRows() == NCols());
   DBG_ASSERT(b.Dim() == NRows());
   DBG_ASSERT(initialized_);
   DBG_ASSERT(factorization_ == LU);

   IpLapackGetrs(NRows(), 1, values_, NRows(), pivot_, b.Values(), b.Dim());
}

void DenseGenMatrix::MultVectorImpl(
   Number        alpha,
   const Vector& x,
   Number        beta,
   Vector&       y
) const
{
   //  A few sanity checks
   DBG_ASSERT(NCols() == x.Dim());
   DBG_ASSERT(NRows() == y.Dim());
   DBG_ASSERT(initialized_);

   // See if we can understand the data
   const DenseVector* dense_x = static_cast<const DenseVector*>(&x);
   DBG_ASSERT(dynamic_cast<const DenseVector*>(&x));
   DenseVector* dense_y = static_cast<DenseVector*>(&y);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&y));

   bool trans = false;
   IpBlasGemv(trans, NRows(), NCols(), alpha, values_, NRows(), dense_x->Values(), 1, beta, dense_y->Values(), 1);
}

void DenseGenMatrix::TransMultVectorImpl(
   Number        alpha,
   const Vector& x,
   Number        beta,
   Vector&       y
) const
{
   //  A few sanity checks
   DBG_ASSERT(NCols() == y.Dim());
   DBG_ASSERT(NRows() == x.Dim());
   DBG_ASSERT(initialized_);

   // See if we can understand the data
   const DenseVector* dense_x = static_cast<const DenseVector*>(&x);
   DBG_ASSERT(dynamic_cast<const DenseVector*>(&x));
   DenseVector* dense_y = static_cast<DenseVector*>(&y);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&y));

   bool trans = true;
   IpBlasGemv(trans, NRows(), NCols(), alpha, values_, NRows(), dense_x->Values(), 1, beta, dense_y->Values(), 1);
}

void DenseGenMatrix::ComputeRowAMaxImpl(
   Vector& rows_norms,
   bool    /*init*/
) const
{
   //  A few sanity checks
   DBG_ASSERT(initialized_);

   DenseVector* dense_vec = static_cast<DenseVector*>(&rows_norms);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&rows_norms));
   Number* vec_vals = dense_vec->Values();

   const Number* vals = values_;
   for( Index irow = 0; irow < NRows(); irow++ )
   {
      for( Index jcol = 0; jcol < NCols(); jcol++ )
      {
         vec_vals[irow] = Max(vec_vals[irow], std::abs(*vals));
         vals++;
      }
   }
}

void DenseGenMatrix::ComputeColAMaxImpl(
   Vector& cols_norms,
   bool    /*init*/
) const
{
   //  A few sanity checks
   DBG_ASSERT(initialized_);

   DenseVector* dense_vec = static_cast<DenseVector*>(&cols_norms);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&cols_norms));
   Number* vec_vals = dense_vec->Values();

   const Number* vals = values_;
   for( Index jcol = 0; jcol < NCols(); jcol++ )
   {
      Index i = IpBlasIamax(NRows(), vals, 1);
      vec_vals[jcol] = Max(vec_vals[jcol], std::abs(vals[i]));
      vals += NRows();
   }
}

bool DenseGenMatrix::HasValidNumbersImpl() const
{
   DBG_ASSERT(initialized_);
   Number sum = IpBlasAsum(NRows() * NCols(), values_, 1);
   return IsFiniteNumber(sum);
}

void DenseGenMatrix::PrintImpl(
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
                        "%sDenseGenMatrix \"%s\" with %" IPOPT_INDEX_FORMAT " rows and %" IPOPT_INDEX_FORMAT " columns:\n", prefix.c_str(), name.c_str(), NRows(), NCols());

   if( initialized_ )
   {
      for( Index j = 0; j < NCols(); j++ )
      {
         for( Index i = 0; i < NRows(); i++ )
         {
            jnlst.PrintfIndented(level, category, indent,
                                 "%s%s[%5" IPOPT_INDEX_FORMAT ",%5" IPOPT_INDEX_FORMAT "]=%23.16e\n", prefix.c_str(), name.c_str(), i, j, values_[i + NRows() * j]);
         }
      }
   }
   else
   {
      jnlst.PrintfIndented(level, category, indent,
                           "The matrix has not yet been initialized!\n");
   }
}

DenseGenMatrixSpace::DenseGenMatrixSpace(
   Index nRows,
   Index nCols
)
   : MatrixSpace(nRows, nCols)
{ }

} // namespace Ipopt
