// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpExpansionMatrix.hpp"
#include "IpDenseVector.hpp"

namespace Ipopt
{

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
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    // See if we can understand the data
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x); /* ToDo: Implement others */
    DenseVector* dense_y = dynamic_cast<DenseVector*>(&y);
    DBG_ASSERT(dense_y); /* ToDo: Implement others */

    const Index* exp_pos = ExpandedPosIndices();

    if (dense_x && dense_y) {
      const Number* xvals=dense_x->Values();
      Number* yvals=dense_y->Values();
      for(Index i=0; i<NCols(); i++) {
        yvals[exp_pos[i]] += alpha * xvals[i];
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
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    // See if we can understand the data
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x); /* ToDo: Implement others */
    DenseVector* dense_y = dynamic_cast<DenseVector*>(&y);
    DBG_ASSERT(dense_y); /* ToDo: Implement others */

    const Index* exp_pos = ExpandedPosIndices();

    if (dense_x && dense_y) {
      const Number* xvals=dense_x->Values();
      Number* yvals=dense_y->Values();
      for(Index i=0; i<NCols(); i++) {
        yvals[i] += alpha * xvals[exp_pos[i]];
      }
    }
  }

  void ExpansionMatrix::PrintImpl(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    fprintf(fp, "\n");
    for (Index ind=0; ind<indent; ind++) {
      fprintf(fp, " ");
    }
    fprintf(fp, "%sExpansionMatrix \"%s\" with %d nonzero elements:\n",
            prefix.c_str(), name.c_str(), NCols());

    const Index* exp_pos = ExpandedPosIndices();

    for (Index i=0; i<NCols(); i++) {
      for (Index ind=0; ind<indent; ind++) {
        fprintf(fp, " ");
      }
      fprintf(fp, "%s%s[%5d,%5d]=%23.16e  (%d)\n", prefix.c_str(), name.c_str(), exp_pos[i]+1,
              i+1, 1., i);
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
    expanded_pos_  = new Index[NCols()];
    compressed_pos_ = new Index[NRows()];
    for (Index j=0; j<NRows(); j++) {
      compressed_pos_[j] = -1;
    }
    for(Index i=0; i<NCols(); i++) {
      //ToDo decide for offset
      DBG_ASSERT(ExpPos[i]-offset<NRows() && ExpPos[i]-offset>=0);
      expanded_pos_[i]=ExpPos[i]-offset;
      compressed_pos_[ExpPos[i]-offset] = i;
    }
  }

} // namespace Ipopt
