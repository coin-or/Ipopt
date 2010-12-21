// Copyright 2009, 2010 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-08

#include "AsStdSchurData.hpp"
#include "IpIteratesVector.hpp"
#include "IpDenseVector.hpp"
#include "IpBlas.hpp"
#include "AsNmpcUtils.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  StdSchurData::StdSchurData()
  {
  }

  StdSchurData::~StdSchurData()
  {
  }

  SmartPtr<SchurData> StdSchurData::MakeNewSchurDataCopy() const
  {
    DBG_START_METH("StdSchurData::MakeNewSchurDataCopy", dbg_verbosity);
    
    SmartPtr<SchurData> retval = new StdSchurData();
   
    return retval;   
  }

  void StdSchurData::SetData_Flag(Index dim, const Index* flags, Number v)
  {
    DBG_START_METH("StdSchurData::SetData_FullIndex", dbg_verbosity);
    DBG_ASSERT(!Is_Initialized());
    DBG_ASSERT(idx1_.empty());
    Index rowcounter = -1;
    for (Index i=0; i<dim; ++i) {
      DBG_ASSERT(flags[i]==1 || flags[i]==0);
      if (flags[i]) {
	idx1_.push_back(++rowcounter);
	idx2_.push_back(i);
	val_.push_back(v);
      }
    }
    Set_NRows(rowcounter+1);
    Set_Initialized();
  }

  void StdSchurData::SetData_Flag(Index dim, const Index* flags, const Number* values)
  {
    DBG_START_METH("StdSchurData::SetData_Full", dbg_verbosity);
    DBG_ASSERT(!Is_Initialized());
    DBG_ASSERT(idx1_.empty());
    Index rowcounter = -1;
    for (Index i=0; i<dim; ++i) {
      DBG_ASSERT(flags[i]==1 || flags[i]==0);
      if (flags[i]) {
	idx1_.push_back(++rowcounter);
	idx2_.push_back(i);
	val_.push_back(values[i]);
      }
    }
    Set_NRows(rowcounter+1);
    Set_Initialized();
  }

  Index StdSchurData::SetData_Index(Index dim, const Index* index, Number v)
  {
    DBG_START_METH("StdSchurData::SetData_Index", dbg_verbosity);
    DBG_ASSERT(!Is_Initialized());
    DBG_ASSERT(idx1_.empty());

    Index n_ind = AsIndexMax(dim, index, 1);
    std::vector<Index> sortvec(n_ind,-1);
    // fill up sortlist
    for (Index i=0; i<dim; ++i) {
      if (index[i]>0) {
	DBG_ASSERT(sortvec[index[i]-1]==-1);
	if ( sortvec[index[i]-1]!=-1 ) {
	  return index[i];
	}
	sortvec[index[i]-1] = i;
      }
    }
    // initialize vectors and write sortlist into vector
    idx1_.resize(n_ind,0);
    idx2_.resize(n_ind,0);
    val_.resize(n_ind,0);
    for (Index i=0; i<n_ind; ++i) {
      DBG_ASSERT(sortvec[i]>-1);
      idx1_[i] = i;
      idx2_[i] = sortvec[i];
      val_[i]  = v;
    }
    Set_NRows(idx1_.size());
    Set_Initialized();
    return 0;
  }

  void StdSchurData::SetData_List(const std::vector<Index>& list, Number v)
  {
    DBG_START_METH("StdSchurData::SetData_List", dbg_verbosity);
    bool not_implemented_yet = false;
    DBG_ASSERT(not_implemented_yet);
  }

  void StdSchurData::GetRow(Index row, IteratesVector& v) const
  {
    DBG_START_METH("StdSchurData::GetRow", dbg_verbosity);
    DBG_ASSERT(Is_Initialized());
    DBG_ASSERT(row<GetNRowsAdded());

    // retrieve structure of IteratesVector - this should probably be cached or sth.
    Index n_comps = v.NComps();

    Index* v_lens = GetVectorLengths(v);

    // set vector v to 0
    v.Set(0.0);
    
    // go through all indices saved in data and find the ones that have the correct row
    Index col, vec_idx;
    for (Index i=0; i<idx2_.size(); ++i) {
      if (idx1_[i]==row) {
	col = idx2_[i];
	
	// find vector in CompoundVector that corresponds to the given col.
	vec_idx = 0;
	while(!(col<v_lens[vec_idx])) { 
	  vec_idx++;
	}
	
	/* I know this is a very long line. It kinda assumes that only few components
	 * will be changed. If that isn't the case, one should cache vectors etc., but
	 * I guess that's for a different implementation of SchurData. */
	dynamic_cast<DenseVector*>(GetRawPtr(v.GetCompNonConst(vec_idx)))->Values()[col+v.GetComp(vec_idx)->Dim()-v_lens[vec_idx]] = val_[i];
      }
    }
    
    delete[] v_lens;
  }

  void StdSchurData::GetMultiplyingVectors(Index row, std::vector<Index>& indices, std::vector<Number>& factors) const
  {
    DBG_START_METH("StdSchurData::GetMultiplyingVectors", dbg_verbosity);

    DBG_ASSERT(indices.size()==0);
    DBG_ASSERT(factors.size()==0);

    for (Index i=0; i<idx1_.size(); ++i) {
      if (idx1_[i] == row) {
	indices.push_back(idx2_[i]);
	factors.push_back(val_[i]);
      }
    }
  }

  void StdSchurData::Multiply(const IteratesVector& v, Vector& u) const
  {
    DBG_START_METH("StdSchurData::Multiply", dbg_verbosity);

    DBG_ASSERT(idx1_.size()==idx2_.size());
    DBG_ASSERT(idx1_.size()==val_.size());
    
    DenseVector* du = static_cast<DenseVector*>(&u);
    du->Set(0.0);
    Number* u_val = du->Values();

    Index* v_lens = GetVectorLengths(v);

    Index v_row, vec_idx;
    for (Index i=0; i<idx1_.size(); ++i) {
      v_row = idx2_[i];
      
      // find vector in CompoundVector that corresponds to the given col in matrix/row in v.
      vec_idx = -1;
      while(!(v_row<v_lens[++vec_idx])) { }

      SmartPtr<const DenseVector> d_ptr = dynamic_cast<const DenseVector*>(GetRawPtr(v.GetComp(vec_idx)));
      if (!d_ptr->IsHomogeneous()) {
	u_val[i] += 
	  val_[i]*d_ptr->Values()[v_row+v.GetComp(vec_idx)->Dim()-v_lens[vec_idx]];
      }
      else {
	u_val[i] += val_[i]*d_ptr->Scalar();
      }
    }
  }

  void StdSchurData::TransMultiply(const Vector& u, IteratesVector& v) const
  {
    DBG_START_METH("StdSchurData::TransMultiply", dbg_verbosity);
    DBG_ASSERT(u.Dim()==GetNRowsAdded());

    const DenseVector* du = static_cast<const DenseVector*>(&u);

    // Get total number of elements of v
    Index ncols = 0;
    for (Index i=0; i<v.NComps(); ++i) {
      ncols += v.GetComp(i)->Dim();
    }

    // Create space in which v_vals will be saved
    Number* v_vals = new Number[ncols];

    const Number* u_vals = du->Values();

    // set v to zero
    for (Index i=0; i<ncols; ++i) {
      v_vals[i] = 0;
    }
    
    // perform v_vals <- A^T*u
    Index row, col;
    Number val;
    for (Index i=0; i<idx1_.size(); ++i) {
      row = idx1_[i];
      col = idx2_[i];
      val = val_[i];

      v_vals[col] += val*u_vals[row];
    }

    // save v_vals in v
    Index v_idx = 0, curr_dim;
    Number* curr_val;
    for (Index i=0; i<v.NComps(); ++i) {
      curr_dim = v.GetCompNonConst(i)->Dim();
      curr_val = dynamic_cast<DenseVector*>(GetRawPtr(v.GetCompNonConst(i)))->Values();
      IpBlasDcopy(curr_dim, v_vals+v_idx, 1, curr_val, 1);

      v_idx +=curr_dim;
    }
  }


  Index* StdSchurData::GetVectorLengths(const IteratesVector& v) const 
  {
    DBG_START_METH("StdSchurData::GetVectorLengths", dbg_verbosity);
    // retrieve structure of IteratesVector - this should probably be cached or sth.
    Index n_comps = v.NComps();
    Index* v_lens = new Index[n_comps];

    // v_lens[i] holds the maximum number up to which component i belongs in there.
    v_lens[0] = v.GetComp(0)->Dim();
    for (Index i=1; i<n_comps; ++i) {
      v_lens[i] = v_lens[i-1]+ v.GetComp(i)->Dim();
    }
    return v_lens;
  }

  void StdSchurData::PrintImpl(const Journalist& jnlst,
			       EJournalLevel level,
			       EJournalCategory category,
			       const std::string& name,
			       Index indent,
			       const std::string& prefix) const
  {
    DBG_START_METH("StdSchurData::PrintImpl", dbg_verbosity);
    
    jnlst.PrintfIndented(level, category, indent,
                         "%sSchurData \"%s\" with %d rows:\n",
                         prefix.c_str(), name.c_str(), GetNRowsAdded());
    if (Is_Initialized()) {
      for (Index i=0; i<idx1_.size(); i++) {
	jnlst.PrintfIndented(level, category, indent,
			     "%s%s[%5d,%5d]=%23.16e\n",
			     prefix.c_str(), name.c_str(), idx1_[i], idx2_[i], val_[i]);
      }
    }
    else {
      jnlst.PrintfIndented(level, category, indent,
                           "%sUninitialized!\n",
                           prefix.c_str());
    }
  }

  void StdSchurData::AddData_List(std::vector<Index> cols, std::vector<Index>& delta_u_sort, Index& new_du_size, Index v)
  {
    DBG_START_METH("StdSchurData::AddData_List", dbg_verbosity);

    //printf("NOOOOOOOOOOOO!, cannot use StdSchurData!!!!");
  }

}
