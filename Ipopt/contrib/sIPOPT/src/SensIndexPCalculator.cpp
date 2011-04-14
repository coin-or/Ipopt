// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-27

#include "SensIndexPCalculator.hpp"
#include "SensIndexSchurData.hpp"
#include "IpDenseVector.hpp"
#include "IpDenseGenMatrix.hpp"
#include "IpBlas.hpp"
#include <vector>

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  IndexPCalculator::IndexPCalculator(SmartPtr<SensBacksolver> backsolver,
				     SmartPtr<SchurData> A_data)
    :
    PCalculator(backsolver, A_data),
    nrows_(0),
    ncols_(A_data->GetNRowsAdded())
  {
    DBG_START_METH("IndexPCalculator::IndexPCalculator", dbg_verbosity);
  }

  IndexPCalculator::~IndexPCalculator()
  {
    DBG_START_METH("IndexPCalculator::~IndexPCalculator", dbg_verbosity);
  }

  bool IndexPCalculator::InitializeImpl(const OptionsList& options,
					const std::string& prefix)
  {
    DBG_START_METH("IndexPCalculator::InitializeImpl", dbg_verbosity);

    SmartPtr<const IteratesVector> iv = IpData().curr();
    nrows_ = 0;
    for (Index i=0; i<iv->NComps(); ++i) {
      nrows_+=iv->GetComp(i)->Dim();
    }
    data_A()->Print(Jnlst(),J_VECTOR,J_USER1,"PCalc SchurData");

    return true;
  }

  bool IndexPCalculator::ComputeP()
  {
    DBG_START_METH("IndexPCalculator::ComputeP", dbg_verbosity);
    bool retval = true;

    // 1. check whether all columns needed by data_A() are in map cols_ - we suppose data_A is IndexSchurData
    const std::vector<Index>* p2col_idx = dynamic_cast<const IndexSchurData*>(GetRawPtr(data_A()))->GetColIndices();
    Index col;
    Number* col_values = NULL;
    Index curr_dim, curr_schur_row=0;
    SmartPtr<const DenseVector> comp_vec;
    const Number* comp_values;
    std::map< Index, SmartPtr<PColumn> >::iterator find_it;
    SmartPtr<IteratesVector> col_vec = IpData().curr()->MakeNewIteratesVector();
    SmartPtr<IteratesVector> sol_vec = col_vec->MakeNewIteratesVector();
    for (std::vector<Index>::const_iterator col_it=p2col_idx->begin(); col_it!=p2col_idx->end(); ++col_it){
      col = *col_it;

      find_it = cols_.find(col);
      if (find_it==cols_.end()) {
	// column is in data_A but not in P-matrix ->create
	data_A()->GetRow(curr_schur_row, *col_vec);
	retval = Solver()->Solve(sol_vec, ConstPtr(col_vec));
	DBG_ASSERT(retval);

	/* This part is for displaying norm2(I_z*K^(-1)*I_1) */
	DBG_PRINT((dbg_verbosity,"\ncurr_schur_row=%d, ",curr_schur_row));
	DBG_PRINT((dbg_verbosity,"norm2(z)=%23.16e\n",sol_vec->x()->Nrm2()));
	/* end displaying norm2 */

	DBG_ASSERT(col_values== NULL);
	col_values = new Number[nrows_];
	curr_dim = 0;
	for (Index j=0; j<sol_vec->NComps(); ++j) {
	  comp_vec = dynamic_cast<const DenseVector*>(GetRawPtr(sol_vec->GetComp(j)));
	  comp_values = comp_vec->Values();
	  IpBlasDcopy(comp_vec->Dim(), comp_values, 1, col_values+curr_dim,1);
	  curr_dim += comp_vec->Dim();
	}
	cols_[col] = new PColumn(nrows_, col_values);
	col_values = NULL;
      }
      curr_schur_row++;
    }

    return retval;
  }

  bool IndexPCalculator::GetSchurMatrix(const SmartPtr<const SchurData>& B, SmartPtr<Matrix>& S)
  {
    DBG_START_METH("IndexPCalculator::GetSchurMatrix", dbg_verbosity);
    bool retval = true;

    Number* S_values;
    if (!IsValid(S)) {
      if ( B==data_A() ) {
	SmartPtr<DenseSymMatrixSpace> S_sym_space = new DenseSymMatrixSpace(B->GetNRowsAdded());
	SmartPtr<DenseSymMatrix> dS = new DenseSymMatrix(GetRawPtr(S_sym_space));
	S_values = dS->Values();
	S = GetRawPtr(dS);
      }
      else {
	SmartPtr<DenseGenMatrixSpace> S_sym_space = new DenseGenMatrixSpace(B->GetNRowsAdded(), B->GetNRowsAdded());
	SmartPtr<DenseGenMatrix> dS = new DenseGenMatrix(GetRawPtr(S_sym_space));
	S_values = dS->Values();
	S = GetRawPtr(dS);
      }
    }
    else {
      // Try DenseGenMatrix - if NULL, try DenseSymMatrix
      SmartPtr<DenseGenMatrix> dS_gen = dynamic_cast<DenseGenMatrix*>(GetRawPtr(S));
      if (!IsValid(dS_gen)) {
	SmartPtr<DenseSymMatrix> dS_sym = dynamic_cast<DenseSymMatrix*>(GetRawPtr(S));
	S_values = dS_sym->Values();
      }
      else {
	S_values = dS_gen->Values();
      }
    }
    /*
      DenseGenMatrix* dS = static_cast<DenseGenMatrix*>(&S);
      DBG_ASSERT(dynamic_cast<const DenseGenMatrix*>(&S));
    */
    // Check whether data_A was changed from the outside
    if (ncols_!=data_A()->GetNRowsAdded()) {
      ncols_ = data_A()->GetNRowsAdded();
      ComputeP();
    }
    /*
      DBG_ASSERT(dS->NRows()==dS->NCols());
      DBG_ASSERT(dS->NRows()==data_A()->GetNRowsAdded());
    */
    std::vector<Index> indices;
    std::vector<Number> factors;

    // Compute S = B^T*P from indices, factors and P
    const std::vector<Index>* data_A_idx = dynamic_cast<const IndexSchurData*>(GetRawPtr(data_A()))->GetColIndices();
    const std::vector<Index>* data_B_idx = dynamic_cast<const IndexSchurData*>(GetRawPtr(B))->GetColIndices();
    Index col_count = 0;
    for (std::vector<Index>::const_iterator a_it=data_A_idx->begin(); a_it!=data_A_idx->end(); ++a_it) {
      cols_[*a_it]->GetSchurMatrixRows(data_B_idx, S_values+col_count*ncols_);
      col_count++;
    }
    return retval;
  }

  void IndexPCalculator::PrintImpl(const Journalist& jnlst,
				   EJournalLevel level,
				   EJournalCategory category,
				   const std::string& name,
				   Index indent,
				   const std::string& prefix) const
  {
    DBG_START_METH("IndexPCalculator::PrintImpl", dbg_verbosity);

    const Number* col_val;
    jnlst.PrintfIndented(level, category, indent,
                         "%sIndexPCalculator \"%s\" with %d rows and %d columns:\n",
                         prefix.c_str(), name.c_str(), nrows_, ncols_ );
    Index col_counter = 0;
    for (std::map< Index, SmartPtr<PColumn> >::const_iterator j=cols_.begin(); j!=cols_.end(); ++j) {
      col_val = j->second->Values();
      for (Index i=0; i<nrows_; ++i) {
	jnlst.PrintfIndented(level, category, indent,
			     "%s%s[%5d,%5d]=%23.16e\n",
			     prefix.c_str(), name.c_str(), i, col_counter, col_val[i]);
      }
      col_counter++;
    }
  }

  PColumn::PColumn(Index nrows, Number* values)
    :
    nrows_(nrows),
    val_(values)
  {
    DBG_START_METH("PColumn::PColumn", dbg_verbosity);
  }

  PColumn::~PColumn()
  {
    DBG_START_METH("PColumn::~PColumn", dbg_verbosity);
    delete[] val_;
  }

  void PColumn::GetSchurMatrixRows(const std::vector<Index>* row_idx_B, Number* S_col) const
  {
    DBG_START_METH("PColumn::GetSchurMatrixRows", dbg_verbosity);

    for (Index i=0; i<row_idx_B->size(); ++i) {
      S_col[i] = -val_[(*row_idx_B)[i]];
    }
  }

  const Number* PColumn::Values() const
  {
    DBG_START_METH("PColumn::Values", dbg_verbosity);
    return val_;
  }
}
