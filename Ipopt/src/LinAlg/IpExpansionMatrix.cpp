// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpExpansionMatrix.hpp"
#include "IpDenseVector.hpp"

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  ExpansionMatrix::ExpansionMatrix(const ExpansionMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      owner_space_(owner_space)
  {}

  ExpansionMatrix::~ExpansionMatrix()
  {}

  void ExpansionMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                       Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NCols()==x.Dim());
    DBG_ASSERT(NRows()==y.Dim());

    // Take care of the y part of the addition
    if ( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    // See if we can understand the data
    const DenseVector* dense_x = static_cast<const DenseVector*>(&x);
    DBG_ASSERT(dynamic_cast<const DenseVector*>(&x));
    DenseVector* dense_y = static_cast<DenseVector*>(&y);
    DBG_ASSERT(dynamic_cast<DenseVector*>(&y));

    const Index* exp_pos = ExpandedPosIndices();

    if (dense_x && dense_y) {
      Number* yvals=dense_y->Values();
      if (dense_x->IsHomogeneous()) {
        Number val = alpha * dense_x->Scalar();
        if (val != 0.) {
          for (Index i=0; i<NCols(); i++) {
            yvals[exp_pos[i]] += val;
          }
        }
      }
      else {
        const Number* xvals=dense_x->Values();
        if (alpha == 1.) {
          for (Index i=0; i<NCols(); i++) {
            yvals[exp_pos[i]] += xvals[i];
          }
        }
        else if (alpha == -1.) {
          for (Index i=0; i<NCols(); i++) {
            yvals[exp_pos[i]] -= xvals[i];
          }
        }
        else {
          for (Index i=0; i<NCols(); i++) {
            yvals[exp_pos[i]] += alpha * xvals[i];
          }
        }
      }
    }
  }

  void ExpansionMatrix::TransMultVectorImpl(Number alpha, const Vector &x,
      Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NCols()==y.Dim());
    DBG_ASSERT(NRows()==x.Dim());

    // Take care of the y part of the addition
    if ( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    // See if we can understand the data
    const DenseVector* dense_x = static_cast<const DenseVector*>(&x);
    DBG_ASSERT(dynamic_cast<const DenseVector*>(&x));
    DenseVector* dense_y = static_cast<DenseVector*>(&y);
    DBG_ASSERT(dynamic_cast<DenseVector*>(&y));

    const Index* exp_pos = ExpandedPosIndices();

    if (dense_x && dense_y) {
      Number* yvals=dense_y->Values();
      if (dense_x->IsHomogeneous()) {
        Number val = alpha * dense_x->Scalar();
        if (val != 0.) {
          for (Index i=0; i<NCols(); i++) {
            yvals[i] += val;
          }
        }
      }
      else {
        const Number* xvals=dense_x->Values();
        if (alpha == 1.) {
          for (Index i=0; i<NCols(); i++) {
            yvals[i] += xvals[exp_pos[i]];
          }
        }
        else if (alpha == -1.) {
          for (Index i=0; i<NCols(); i++) {
            yvals[i] -= xvals[exp_pos[i]];
          }
        }
        else {
          for (Index i=0; i<NCols(); i++) {
            yvals[i] += alpha * xvals[exp_pos[i]];
          }
        }
      }
    }
  }

  // Specialized method (overloaded from IpMatrix)
  void ExpansionMatrix::AddMSinvZImpl(Number alpha, const Vector& S,
                                      const Vector& Z, Vector& X) const
  {
    DBG_ASSERT(NCols()==S.Dim());
    DBG_ASSERT(NCols()==Z.Dim());
    DBG_ASSERT(NRows()==X.Dim());

    const DenseVector* dense_S = static_cast<const DenseVector*>(&S);
    DBG_ASSERT(dynamic_cast<const DenseVector*>(&S));
    const DenseVector* dense_Z = static_cast<const DenseVector*>(&Z);
    DBG_ASSERT(dynamic_cast<const DenseVector*>(&Z));
    DenseVector* dense_X = static_cast<DenseVector*>(&X);
    DBG_ASSERT(dynamic_cast<DenseVector*>(&X));

    // if vector S is homogeneous type, call the default implementation
    // ToDo: find out how often the default implementation is called and
    // if we should implement specialized method for the homogenenous
    // case
    if (dense_S->IsHomogeneous()) {
      DBG_ASSERT(false && "dense_S is homogeneous - implement specialized option?");
      Matrix::AddMSinvZImpl(alpha, S, Z, X);
      return;
    }

    const Index* exp_pos = ExpandedPosIndices();
    const Number* vals_S = dense_S->Values();
    Number* vals_X = dense_X->Values();

    if (dense_Z->IsHomogeneous()) {
      Number val = alpha*dense_Z->Scalar();
      if (val != 0.) {
        for (Index i=0; i<NCols(); i++) {
          vals_X[exp_pos[i]] += val/vals_S[i];
        }
      }
    }
    else {
      const Number* vals_Z = dense_Z->Values();
      if (alpha==1.) {
        for (Index i=0; i<NCols(); i++) {
          vals_X[exp_pos[i]] += vals_Z[i]/vals_S[i];
        }
      }
      else if (alpha==-1.) {
        for (Index i=0; i<NCols(); i++) {
          vals_X[exp_pos[i]] -= vals_Z[i]/vals_S[i];
        }
      }
      else {
        for (Index i=0; i<NCols(); i++) {
          vals_X[exp_pos[i]] += alpha*vals_Z[i]/vals_S[i];
        }
      }
    }
  }

  void ExpansionMatrix::SinvBlrmZMTdBrImpl(Number alpha, const Vector& S,
      const Vector& R, const Vector& Z,
      const Vector& D, Vector& X) const
  {
    DBG_START_METH("ExpansionMatrix::SinvBlrmZMTdBrImpl",
                   dbg_verbosity);

    DBG_ASSERT(NCols()==S.Dim());
    DBG_ASSERT(NCols()==R.Dim());
    DBG_ASSERT(NCols()==Z.Dim());
    DBG_ASSERT(NRows()==D.Dim());
    DBG_ASSERT(NCols()==X.Dim());

    const DenseVector* dense_S = static_cast<const DenseVector*>(&S);
    DBG_ASSERT(dynamic_cast<const DenseVector*>(&S));
    const DenseVector* dense_R = static_cast<const DenseVector*>(&R);
    DBG_ASSERT(dynamic_cast<const DenseVector*>(&R));
    const DenseVector* dense_Z = static_cast<const DenseVector*>(&Z);
    DBG_ASSERT(dynamic_cast<const DenseVector*>(&Z));
    const DenseVector* dense_D = static_cast<const DenseVector*>(&D);
    DBG_ASSERT(dynamic_cast<const DenseVector*>(&D));
    DenseVector* dense_X = static_cast<DenseVector*>(&X);
    DBG_ASSERT(dynamic_cast<DenseVector*>(&X));

    // if the vectors S or D are of the homogeneous type, revert to the
    // default implementation
    // ToDo: find out how often the default implementation is called and
    // if we should implement specialized method for the homogenenous
    // case

    if (dense_S->IsHomogeneous() ||
        dense_D->IsHomogeneous()) {
      DBG_ASSERT(false && "dense_S or dense_D is homogeneous - implement specialized option?");
      Matrix::SinvBlrmZMTdBrImpl(alpha, S, R, Z, D, X);
      return;
    }

    const Index* exp_pos = ExpandedPosIndices();
    const Number* vals_S = dense_S->Values();
    const Number* vals_D = dense_D->Values();
    Number* vals_X = dense_X->Values();

    if (dense_R->IsHomogeneous()) {
      Number scalar_R = dense_R->Scalar();
      if (dense_Z->IsHomogeneous()) {
        Number val = alpha*dense_Z->Scalar();
        if (val == 0.) {
          for (Index i=0; i<NCols(); i++) {
            // ToDo could treat val == 0 extra
            vals_X[i] = (scalar_R)/vals_S[i];
          }
        }
        else {
          for (Index i=0; i<NCols(); i++) {
            // ToDo could treat val == 0 extra
            vals_X[i] = (scalar_R + val*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
      }
      else {
        const Number* vals_Z = dense_Z->Values();
        if (alpha==1.) {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (scalar_R + vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
        else if (alpha==-1.) {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (scalar_R - vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
        else {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (scalar_R + alpha*vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
      }
    }
    else {
      const Number* vals_R = dense_R->Values();
      if (dense_Z->IsHomogeneous()) {
        Number val = alpha*dense_Z->Scalar();
        for (Index i=0; i<NCols(); i++) {
          vals_X[i] = (vals_R[i] + val*vals_D[exp_pos[i]])/vals_S[i];
        }
      }
      else {
        const Number* vals_Z = dense_Z->Values();
        if (alpha==1.) {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (vals_R[i] + vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
        else if (alpha==-1.) {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (vals_R[i] - vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
        else {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (vals_R[i] + alpha*vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
      }
    }
  }

  void ExpansionMatrix::ComputeRowAMaxImpl(Vector& rows_norms, bool init) const
  {
    DenseVector* dense_vec = static_cast<DenseVector*>(&rows_norms);
    DBG_ASSERT(dynamic_cast<DenseVector*>(&rows_norms));
    Number* vec_vals=dense_vec->Values();

    const Index* exp_pos = ExpandedPosIndices();

    for (Index i=0; i<NCols(); i++) {
      vec_vals[exp_pos[i]] = Max(vec_vals[exp_pos[i]], 1.);
    }
  }

  void ExpansionMatrix::ComputeColAMaxImpl(Vector& cols_norms, bool init) const
  {
    if (init) {
      cols_norms.Set(1.);
    }
    else {
      SmartPtr<Vector> v = cols_norms.MakeNew();
      v->Set(1.);
      cols_norms.ElementWiseMax(*v);
    }
  }

  void ExpansionMatrix::PrintImplOffset(const Journalist& jnlst,
                                        EJournalLevel level,
                                        EJournalCategory category,
                                        const std::string& name,
                                        Index indent,
                                        const std::string& prefix,
                                        Index row_offset,
                                        Index col_offset) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sExpansionMatrix \"%s\" with %d rows and %d columns:\n",
                         prefix.c_str(), name.c_str(), NRows(), NCols());

    const Index* exp_pos = ExpandedPosIndices();

    for (Index i=0; i<NCols(); i++) {
      jnlst.PrintfIndented(level, category, indent,
                           "%s%s[%5d,%5d]=%23.16e  (%d)\n",
                           prefix.c_str(), name.c_str(), exp_pos[i]+row_offset,
                           i+col_offset, 1., i);
    }
  }

  ExpansionMatrixSpace::ExpansionMatrixSpace(Index NLargeVec,
      Index NSmallVec,
      const Index *ExpPos,
      const int offset /*= 0*/)
      :
      MatrixSpace(NLargeVec, NSmallVec),
      expanded_pos_(NULL),
      compressed_pos_(NULL)
  {
    if (NCols()>0) {
      expanded_pos_  = new Index[NCols()];
    }
    if (NRows()>0) {
      compressed_pos_ = new Index[NRows()];
    }
    for (Index j=0; j<NRows(); j++) {
      compressed_pos_[j] = -1;
    }
    for (Index i=0; i<NCols(); i++) {
      //ToDo decide for offset
      DBG_ASSERT(ExpPos[i]-offset<NRows() && ExpPos[i]-offset>=0);
      expanded_pos_[i]=ExpPos[i]-offset;
      compressed_pos_[ExpPos[i]-offset] = i;
    }
  }

} // namespace Ipopt
