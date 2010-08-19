// Copyright 2009, 2010 Hans Pirnay
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Date   : 2009-05-06

#include "AsStdPCalculator.hpp"
#include "IpDenseVector.hpp"
#include "IpDenseGenMatrix.hpp"
#include "IpBlas.hpp"
#include <vector>

namespace Ipopt 
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  StdPCalculator::StdPCalculator(SmartPtr<AsBacksolver> backsolver,
				 SmartPtr<SchurData> A_data)
    :
    PCalculator(backsolver, A_data),
    nrows_(0),
    ncols_(A_data->GetNRowsAdded()),
    P_(NULL)
  {
    DBG_START_METH("StdPCalculator::StdPCalculator", dbg_verbosity);
  }

  StdPCalculator::~StdPCalculator()
  {
    DBG_START_METH("StdPCalculator::~StdPCalculator", dbg_verbosity);
    delete[] P_;
  }

  /** Overloaded from AlgorithmStrategyObject */
  bool StdPCalculator::InitializeImpl(const OptionsList& options,
		      const std::string& prefix)
  {
    DBG_START_METH("StdPCalculator::InitializeImpl", dbg_verbosity);

    SmartPtr<const IteratesVector> iv = IpData().curr(); 

    nrows_ = 0;

    for (Index i=0; i<iv->NComps(); ++i) {
      nrows_+=iv->GetComp(i)->Dim();
    }
    data_A()->Print(Jnlst(),J_VECTOR,J_USER1,"PCalc SchurData");

    return true;
  }


  /** Function to start the computation of  P from E_0 and KKT*/
  bool StdPCalculator::ComputeP()
  {
    DBG_START_METH("StdPCalculator::ComputeP", dbg_verbosity);
    bool retval;

    DBG_ASSERT(P_==NULL);
    P_ = new Number[nrows_*ncols_];
    // For debugging purposes, set all elements of P_ to a recognizable value
    for (Index i=0; i<nrows_*ncols_; ++i) {
      P_[i] = -4;
    }
    /* Declarations */
    SmartPtr<IteratesVector> col_vec = IpData().curr()->MakeNewIteratesVector();
    SmartPtr<IteratesVector> sol_vec   = col_vec->MakeNewIteratesVector();
    SmartPtr<const DenseVector> comp_vec;
    const Number* comp_values;
    Index curr_dim; // keeps track of at which row in the matrix we are (in copying the compound vector)

    // fill up the matrix by looping through each column an solving for it
    for (Index col=0; col<ncols_; ++col) {
      data_A()->GetRow(col, *col_vec);
      retval = Solver()->Solve(sol_vec, ConstPtr(col_vec));
      col_vec->Print(Jnlst(),J_VECTOR,J_USER1,"col_vec PCalc");
      sol_vec->Print(Jnlst(),J_VECTOR,J_USER1,"sol_vec PCalc");
      DBG_ASSERT(retval);
      
      curr_dim = 0;
      for (Index j=0; j<sol_vec->NComps(); ++j) {
	comp_vec = dynamic_cast<const DenseVector*>(GetRawPtr(sol_vec->GetComp(j)));
	comp_values = comp_vec->Values();
	IpBlasDcopy(comp_vec->Dim(), comp_values, 1, P_+curr_dim+col*nrows_,1);
      	curr_dim += comp_vec->Dim();
      }
    }
    return retval;
  }
    
  /** Function to extract a SchurMatrix corresponding to $B K^{-1} A^T$ */
  bool StdPCalculator::GetSchurMatrix(const SmartPtr<const SchurData>& B, SmartPtr<Matrix>& S)
  {
    DBG_START_METH("StdPCalculator::GetSchurMatrix", dbg_verbosity);
    bool retval;

    Number* S_values;
    if (!IsValid(S)) {
      // Create new matrix. If B==A, the matrix is symmetric.
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
  
    DBG_ASSERT(dS->NRows()==dS->NCols());
    DBG_ASSERT(dS->NRows()==ncols_);
    */
    std::vector<Index> indices;
    std::vector<Number> factors;

    // Compute S = B^T*P from indices, factors and P
    Index n_ind;
    for (Index row=0; row<ncols_; ++row) {
      indices.clear();
      factors.clear();
      B->GetMultiplyingVectors(row, indices, factors); 
      n_ind = indices.size();
      DBG_ASSERT(n_ind==factors.size());
      for (Index col=0; col<ncols_; ++col) {
	S_values[row+col*ncols_] = 0;
	for (Index i=0; i<n_ind; ++i) {
	  S_values[row+col*ncols_] -= P_[indices[i]+col*nrows_]*factors[i]; // minus from formula S = -B^T*S*A comes into play here.
	}
      }
    }

    return retval;
  }

  void StdPCalculator::PrintImpl(const Journalist& jnlst,
		 EJournalLevel level,
		 EJournalCategory category,
		 const std::string& name,
		 Index indent,
		 const std::string& prefix) const
  {
    DBG_START_METH("StdPCalculator::PrintImpl", dbg_verbosity);

    jnlst.PrintfIndented(level, category, indent,
                         "%sPCalculator \"%s\" with %d rows and %d columns:\n",
                         prefix.c_str(), name.c_str(), nrows_, ncols_ );
    for (Index j=0; j<ncols_; ++j) {
      for (Index i=0; i<nrows_; ++i) {
	jnlst.PrintfIndented(level, category, indent,
			     "%s%s[%5d,%5d]=%23.16e\n",
			     prefix.c_str(), name.c_str(), i, j, P_[i+j*nrows_]);
      }
    }
  }


}
