// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpCompoundMatrix.hpp"
#include "IpCompoundVector.hpp"

namespace Ipopt
{

  CompoundMatrix::CompoundMatrix(const CompoundMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      owner_space_(owner_space),
      matrices_valid_(false)
  {
    for (Index irow=0; irow<NComps_Rows(); irow++) {
      std::vector<SmartPtr<Matrix> > row(NComps_Cols());
      std::vector<SmartPtr<const Matrix> > const_row(NComps_Cols());
      const_comps_.push_back(const_row);
      comps_.push_back(row);
    }
  }


  CompoundMatrix::~CompoundMatrix()
  {}

  void CompoundMatrix::SetComp(Index irow, Index jcol, const Matrix& matrix)
  {
    DBG_ASSERT(irow<NComps_Rows());
    DBG_ASSERT(jcol<NComps_Cols());
    DBG_ASSERT(owner_space_->GetCompSpace(irow, jcol)->IsMatrixFromSpace(matrix));

    comps_[irow][jcol] = NULL;
    const_comps_[irow][jcol] = &matrix;
    ObjectChanged();
  }

  void CompoundMatrix::SetCompNonConst(Index irow, Index jcol, Matrix& matrix)
  {
    DBG_ASSERT(irow < NComps_Rows());
    DBG_ASSERT(jcol < NComps_Cols());
    DBG_ASSERT(owner_space_->GetCompSpace(irow, jcol)->IsMatrixFromSpace(matrix));

    const_comps_[irow][jcol] = NULL;
    comps_[irow][jcol] = &matrix;
    ObjectChanged();
  }

  void CompoundMatrix::CreateBlockFromSpace(Index irow, Index jcol)
  {
    DBG_ASSERT(irow < NComps_Rows());
    DBG_ASSERT(jcol < NComps_Cols());
    DBG_ASSERT(IsValid(owner_space_->GetCompSpace(irow, jcol)));
    SetCompNonConst(irow, jcol, *owner_space_->GetCompSpace(irow,jcol)->MakeNew());
  }

  void CompoundMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                      Number beta, Vector &y) const
  {
    if (!matrices_valid_) {
      matrices_valid_ = MatricesValid();
    }
    DBG_ASSERT(matrices_valid_);

    // The vectors are assumed to be compound Vectors as well
    const CompoundVector* comp_x = dynamic_cast<const CompoundVector*>(&x);
    CompoundVector* comp_y = dynamic_cast<CompoundVector*>(&y);

    //  A few sanity checks
    if (comp_x) {
      DBG_ASSERT(NComps_Cols()==comp_x->NComps());
    }
    else {
      DBG_ASSERT(NComps_Cols() == 1);
    }

    if (comp_y) {
      DBG_ASSERT(NComps_Rows()==comp_y->NComps());
    }
    else {
      DBG_ASSERT(NComps_Rows() == 1);
    }

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    for( Index irow = 0; irow < NComps_Rows(); irow++ ) {
      SmartPtr<Vector> y_i;
      if (comp_y) {
        y_i = comp_y->GetCompNonConst(irow);
      }
      else {
        y_i = &y;
      }
      DBG_ASSERT(IsValid(y_i));

      for( Index jcol = 0; jcol < NComps_Cols(); jcol++ ) {
        SmartPtr<const Vector> x_j;
        if (comp_x) {
          x_j = comp_x->GetComp(jcol);
        }
        else if (NComps_Cols() == 1) {
          x_j = &x;
        }
        DBG_ASSERT(IsValid(x_j));

        if (ConstComp(irow,jcol)) {
          ConstComp(irow, jcol)->MultVector(alpha, *x_j,
                                            1., *y_i);
        }
      }
    }
  }

  void CompoundMatrix::TransMultVectorImpl(Number alpha, const Vector &x,
      Number beta, Vector &y) const
  {
    if (!matrices_valid_) {
      matrices_valid_ = MatricesValid();
    }
    DBG_ASSERT(matrices_valid_);

    // The vectors are assumed to be compound Vectors as well
    const CompoundVector* comp_x = dynamic_cast<const CompoundVector*>(&x);
    CompoundVector* comp_y = dynamic_cast<CompoundVector*>(&y);

    //  A few sanity checks
    if (comp_y) {
      DBG_ASSERT(NComps_Cols()==comp_y->NComps());
    }
    else {
      DBG_ASSERT(NComps_Cols() == 1);
    }

    if (comp_x) {
      DBG_ASSERT(NComps_Rows()==comp_x->NComps());
    }
    else {
      DBG_ASSERT(NComps_Rows() == 1);
    }

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    for( Index irow = 0; irow < NComps_Cols(); irow++ ) {
      SmartPtr<Vector> y_i;
      if (comp_y) {
        y_i = comp_y->GetCompNonConst(irow);
      }
      else {
        y_i = &y;
      }
      DBG_ASSERT(IsValid(y_i));

      for( Index jcol = 0; jcol < NComps_Rows(); jcol++ ) {
        SmartPtr<const Vector> x_j;
        if (comp_x) {
          x_j = comp_x->GetComp(jcol);
        }
        else {
          x_j = &x;
        }
        DBG_ASSERT(IsValid(x_j));

        if (ConstComp(jcol, irow)) {
          ConstComp(jcol, irow)->TransMultVector(alpha, *x_j,
                                                 1., *y_i);
        }
      }
    }
  }

  // Specialized method (overloaded from IpMatrix)
  void CompoundMatrix::AddMSinvZImpl(Number alpha, const Vector& S,
				     const Vector& Z, Vector& X) const
  {
    // The vectors are assumed to be compound Vectors as well (unless they
    // are assumed to consist of only one component
    const CompoundVector* comp_S = dynamic_cast<const CompoundVector*>(&S);
    const CompoundVector* comp_Z = dynamic_cast<const CompoundVector*>(&Z);
    CompoundVector* comp_X = dynamic_cast<CompoundVector*>(&X);

    //  A few sanity checks for sizes
    if (comp_S) {
      DBG_ASSERT(NComps_Cols()==comp_S->NComps());
    }
    else {
      DBG_ASSERT(NComps_Cols() == 1);
    }
    if (comp_Z) {
      DBG_ASSERT(NComps_Cols()==comp_Z->NComps());
    }
    else {
      DBG_ASSERT(NComps_Cols() == 1);
    }
    if (comp_X) {
      DBG_ASSERT(NComps_Rows()==comp_X->NComps());
    }
    else {
      DBG_ASSERT(NComps_Rows() == 1);
    }

    for( Index irow = 0; irow < NComps_Cols(); irow++ ) {
      SmartPtr<Vector> X_i;
      if (comp_X) {
        X_i = comp_X->GetCompNonConst(irow);
      }
      else {
        X_i = &X;
      }
      DBG_ASSERT(IsValid(X_i));

      for( Index jcol = 0; jcol < NComps_Rows(); jcol++ ) {
        SmartPtr<const Vector> S_j;
        if (comp_S) {
          S_j = comp_S->GetComp(jcol);
        }
        else {
          S_j = &S;
        }
        DBG_ASSERT(IsValid(S_j));
        SmartPtr<const Vector> Z_j;
        if (comp_Z) {
          Z_j = comp_Z->GetComp(jcol);
        }
        else {
          Z_j = &Z;
        }
        DBG_ASSERT(IsValid(Z_j));

        if (ConstComp(jcol, irow)) {
          ConstComp(jcol, irow)->AddMSinvZ(alpha, *S_j, *Z_j, *X_i);
        }
      }
    }
  }

  // Specialized method (overloaded from IpMatrix)
  void CompoundMatrix::SinvBlrmZMTdBrImpl(Number alpha, const Vector& S,
					  const Vector& R, const Vector& Z,
					  const Vector& D, Vector& X) const
  {
    // First check if the matrix is indeed such that we can use the
    // special methods from the component spaces (this only works if
    // we have exactly one submatrix per column)

    // ToDo: Do this test only once?
    bool fast_SinvBlrmZMTdBr = true;
    for (Index jcol=0; jcol < NComps_Cols(); jcol++ ) {
      Index nblocks = 0;
      for (Index irow=0; irow < NComps_Rows(); irow++ ) {
	if (ConstComp(irow, jcol)) {
	  nblocks++;
	}
      }
      if (nblocks!=1) {
	fast_SinvBlrmZMTdBr = false;
	break;
      }
    }

    if (!fast_SinvBlrmZMTdBr) {
      // Use the standard replacement implementation
      Matrix::SinvBlrmZMTdBrImpl(alpha, S, R, Z, D, X);
      DBG_ASSERT(false && "Found a matrix where we can't use the fast SinvBlrmZMTdBr implementation in CompoundMatrix");
      return;
    }

    // The vectors are assumed to be compound Vectors as well (unless they
    // are assumed to consist of only one component
    const CompoundVector* comp_S = dynamic_cast<const CompoundVector*>(&S);
    const CompoundVector* comp_R = dynamic_cast<const CompoundVector*>(&R);
    const CompoundVector* comp_Z = dynamic_cast<const CompoundVector*>(&Z);
    const CompoundVector* comp_D = dynamic_cast<const CompoundVector*>(&D);
    CompoundVector* comp_X = dynamic_cast<CompoundVector*>(&X);

    //  A few sanity checks for sizes
    if (comp_S) {
      DBG_ASSERT(NComps_Cols()==comp_S->NComps());
    }
    else {
      DBG_ASSERT(NComps_Cols() == 1);
    }
    if (comp_Z) {
      DBG_ASSERT(NComps_Cols()==comp_Z->NComps());
    }
    else {
      DBG_ASSERT(NComps_Cols() == 1);
    }
    if (comp_R) {
      DBG_ASSERT(NComps_Cols()==comp_R->NComps());
    }
    else {
      DBG_ASSERT(NComps_Cols() == 1);
    }
    if (comp_D) {
      DBG_ASSERT(NComps_Rows()==comp_D->NComps());
    }
    else {
      DBG_ASSERT(NComps_Rows() == 1);
    }
    if (comp_X) {
      DBG_ASSERT(NComps_Cols()==comp_X->NComps());
    }
    else {
      DBG_ASSERT(NComps_Cols() == 1);
    }

    for (Index irow=0; irow<NComps_Cols(); irow++ ) {
      // Find the entry in the matrix for that row
      Index jcol_found=NComps_Rows();
      for (Index jcol=0; jcol<NComps_Rows(); jcol++ ) {
	if (ConstComp(irow, jcol)) {
	  jcol_found = jcol;
	  break;
	}
      }
      DBG_ASSERT(jcol_found<NComps_Rows());

      SmartPtr<const Vector> S_i;
      if (comp_S) {
	S_i = comp_S->GetComp(irow);
      }
      else {
	S_i = &S;
      }
      DBG_ASSERT(IsValid(S_i));
      SmartPtr<const Vector> Z_i;
      if (comp_Z) {
	Z_i = comp_Z->GetComp(irow);
      }
      else {
	Z_i = &Z;
      }
      DBG_ASSERT(IsValid(Z_i));
      SmartPtr<const Vector> R_i;
      if (comp_R) {
	R_i = comp_R->GetComp(irow);
      }
      else {
	R_i = &R;
      }
      DBG_ASSERT(IsValid(R_i));
      SmartPtr<const Vector> D_i;
      if (comp_D) {
	D_i = comp_D->GetComp(jcol_found);
      }
      else {
	D_i = &D;
      }
      DBG_ASSERT(IsValid(D_i));
      SmartPtr<Vector> X_i;
      if (comp_X) {
	X_i = comp_X->GetCompNonConst(irow);
      }
      else {
	X_i = &X;
      }
      DBG_ASSERT(IsValid(X_i));

      ConstComp(jcol_found,irow)->SinvBlrmZMTdBr(alpha, *S_i, *R_i, *Z_i,
						 *D_i, *X_i);
    }
  }

  void CompoundMatrix::PrintImpl(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    fprintf(fp, "\n");
    for (Index ind=0; ind<indent; ind++) {
      fprintf(fp, " ");
    }
    fprintf(fp, "%sCompoundMatrix \"%s\" with %d row and %d columns components:\n",
            prefix.c_str(), name.c_str(), NComps_Rows(), NComps_Cols());
    for (Index irow = 0; irow < NComps_Rows(); irow++ ) {
      for (Index jcol = 0; jcol < NComps_Cols(); jcol++ ) {
        for (Index ind=0; ind<indent; ind++) {
          fprintf(fp, " ");
        }
        fprintf(fp, "%sComponent for row %d and column %d:\n", prefix.c_str(), irow, jcol);
        if (ConstComp(irow, jcol)) {
          char buffer[256];
          sprintf(buffer, "%s[%d][%d]", name.c_str(), irow, jcol);
          std::string term_name = buffer;
          ConstComp(irow, jcol)->Print(fp, term_name, indent+1, prefix);
        }
        else {
          for (Index ind=0; ind<indent; ind++) {
            fprintf(fp, " ");
          }
          fprintf(fp, "%sComponent has not been set.\n", prefix.c_str());
        }
      }
    }
  }

  bool CompoundMatrix::MatricesValid() const
  {
    // Check to make sure we have matrices everywhere the space has matrices
    // We already check that the matrix agrees with the block space
    // in the SetComp methods
    bool retValue = true;
    for (Index i=0; i<NComps_Rows(); i++) {
      for (Index j=0; j<NComps_Cols(); j++) {
        if ( (!ConstComp(i,j) && IsValid(owner_space_->GetCompSpace(i,j)))
             || (ConstComp(i,j) && IsNull(owner_space_->GetCompSpace(i,j))) ) {
          retValue = false;
          break;
        }
      }
    }
    return retValue;
  }


  CompoundMatrixSpace::CompoundMatrixSpace(Index ncomps_rows, Index ncomps_cols, Index total_nRows, Index total_nCols)
      :
      MatrixSpace(total_nRows, total_nCols),
      ncomps_rows_(ncomps_rows),
      ncomps_cols_(ncomps_cols),
      block_rows_(ncomps_rows, -1),
      block_cols_(ncomps_cols, -1),
      dimensions_set_(false)
  {
    DBG_START_METH("CompoundMatrixSpace::CompoundMatrixSpace", 0);
    std::vector<SmartPtr<const MatrixSpace> > row(ncomps_cols_);
    std::vector< bool > allocate_row(ncomps_cols_, false);
    for (Index i=0; i<ncomps_rows_; i++) {
      DBG_PRINT((1, "block_rows_[%d] = %d\n", i, block_rows_[i]));
      comp_spaces_.push_back(row);
      allocate_block_.push_back(allocate_row);
    }
  }

  void CompoundMatrixSpace::SetBlockCols(Index jcol, Index ncols)
  {
    DBG_ASSERT(!dimensions_set_ && "for now, if the dimensions have all been set, they cannot be changed");
    DBG_ASSERT(jcol < ncomps_cols_);
    DBG_ASSERT(block_cols_[jcol] == -1 && "This dimension has already been set - sanity check");
    block_cols_[jcol] = ncols;
  }

  void CompoundMatrixSpace::SetBlockRows(Index irow, Index nrows)
  {
    DBG_ASSERT(!dimensions_set_ && "for now, if the dimensions have all been set, they cannot be changed");
    DBG_ASSERT(irow < ncomps_rows_);
    DBG_ASSERT(block_rows_[irow] == -1 && "This dimension has already been set - sanity check");
    block_rows_[irow] = nrows;
  }

  Index CompoundMatrixSpace::GetBlockRows(Index irow) const
  {
    DBG_ASSERT(dimensions_set_);
    DBG_ASSERT(irow < ncomps_rows_);
    return block_rows_[irow];
  }

  Index CompoundMatrixSpace::GetBlockCols(Index jcol) const
  {
    DBG_ASSERT(dimensions_set_);
    DBG_ASSERT(jcol < ncomps_cols_);
    return block_cols_[jcol];
  }

  void CompoundMatrixSpace::SetCompSpace(Index irow, Index jcol,
                                         const MatrixSpace& mat_space,
                                         bool auto_allocate /*=false*/)
  {
    if (!dimensions_set_) {
      dimensions_set_ = DimensionsSet();
    }
    DBG_ASSERT(dimensions_set_);
    DBG_ASSERT(irow<ncomps_rows_);
    DBG_ASSERT(jcol<ncomps_cols_);
    DBG_ASSERT(IsNull(comp_spaces_[irow][jcol]));
    DBG_ASSERT(block_cols_[jcol] != -1 && block_cols_[jcol] == mat_space.NCols());
    DBG_ASSERT(block_rows_[irow] != -1 && block_rows_[irow] == mat_space.NRows());

    comp_spaces_[irow][jcol] = &mat_space;
    allocate_block_[irow][jcol] = auto_allocate;
  }

  CompoundMatrix* CompoundMatrixSpace::MakeNewCompoundMatrix() const
  {
    if (!dimensions_set_) {
      dimensions_set_ = DimensionsSet();
    }
    DBG_ASSERT(dimensions_set_);

    CompoundMatrix* mat = new CompoundMatrix(this);
    for(Index i=0; i<ncomps_rows_; i++) {
      for (Index j=0; j<ncomps_cols_; j++) {
        if (allocate_block_[i][j]) {
          mat->SetCompNonConst(i, j, *GetCompSpace(i, j)->MakeNew());
        }
      }
    }

    return mat;
  }

  bool CompoundMatrixSpace::DimensionsSet() const
  {
    DBG_START_METH("CompoundMatrixSpace::DimensionsSet", 0);
    Index total_nrows = 0;
    Index total_ncols = 0;
    bool valid = true;
    for (Index i=0; i<ncomps_rows_; i++) {
      if (block_rows_[i] == -1) {
        valid = false;
        break;
      }
      total_nrows += block_rows_[i];
    }
    if (valid) {
      for (Index j=0; j<ncomps_cols_; j++) {
        if (block_cols_[j] == -1) {
          valid = false;
          break;
        }
        total_ncols += block_cols_[j];
      }
    }

    if (valid) {
      DBG_ASSERT(total_nrows == NRows() && total_ncols == NCols());
    }

    return valid;
  }

} // namespace Ipopt
