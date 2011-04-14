// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-27

#include "SensIndexSchurData.hpp"
#include "IpIteratesVector.hpp"
#include "IpDenseVector.hpp"
#include "IpBlas.hpp"
#include "SensUtils.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  IndexSchurData::IndexSchurData()
  {
    DBG_START_METH("IndexSchurData::IndexSchurData", dbg_verbosity);
  }

  IndexSchurData::IndexSchurData(const std::vector<Index> idx, const std::vector<Index> val)
  {
    DBG_START_METH("IndexSchurData::IndexSchurData(vector,vector)", dbg_verbosity);
    idx_ = idx;
    val_ = val;

    Set_NRows((Index)idx_.size());
    Set_Initialized();
  }


  IndexSchurData::~IndexSchurData()
  {
    DBG_START_METH("IndexSchurData::~IndexSchurData", dbg_verbosity);
  }

  SmartPtr<SchurData> IndexSchurData::MakeNewSchurDataCopy() const
  {
    DBG_START_METH("IndexSchurData::MakeNewSchurDataCopy", dbg_verbosity);
    SmartPtr<SchurData> retval = new IndexSchurData(idx_, val_);
    return retval;
  }

  void IndexSchurData::SetData_Flag(Index dim, const Index* flags, Number v)
  {
    DBG_START_METH("IndexSchurData::SetData_Flag", dbg_verbosity);
    DBG_ASSERT(idx_.size()==0);
    DBG_ASSERT(!Is_Initialized());
    DBG_ASSERT(v!=0);

    Index w;
    (v>0) ? w=1 : w=-1;

    for (Index i=0; i<dim; ++i) {
      DBG_ASSERT(flags[i]==1 || flags[i]==0);
      DBG_ASSERT(v!=0);
      if (flags[i]) {
	idx_.push_back(i);
	val_.push_back(w);
      }
    }
    Set_Initialized();
    Set_NRows((Index)val_.size());
  }

  void IndexSchurData::SetData_Flag(Index dim, const Index* flags, const Number* values)
  {
    DBG_START_METH("InexSchurData::SetData_Flag", dbg_verbosity);
    DBG_ASSERT(idx_.size()==0);
    DBG_ASSERT(!Is_Initialized());

    for (Index i=0; i<dim; ++i) {
      DBG_ASSERT(flags[i]==1 || flags[i]==0);
      DBG_ASSERT(values[i]!=0);
      if (flags[i]) {
	idx_.push_back(i);
	(values[i]>0) ?	val_.push_back(1) : val_.push_back(-1);
      }
    }
    Set_Initialized();
    Set_NRows((Index)val_.size());
  }

  Index IndexSchurData::SetData_Index(Index dim, const Index* index, Number v)
  {
    DBG_START_METH("IndexSchurData::SetData_Index", dbg_verbosity);
    DBG_ASSERT(idx_.empty());
    DBG_ASSERT(!Is_Initialized());

    Index w;
    (v>0) ? w=1 : w=-1;
    DBG_PRINT((dbg_verbosity, "Schurdata::w=%d\n", w));
    Index n_ind = AsIndexMax(dim, index, 1);
    std::vector<Index> sortvec(n_ind,-1);
    // fill up sortlist
    for (Index i=0; i<dim; ++i) {
      if (index[i]>0) {
	DBG_ASSERT(sortvec[index[i]-1]==-1); // THIS SHOULD THROW AN EXCEPTION! (OR SWITCH TO FLAG?)
	if ( sortvec[index[i]-1]!=-1 ) {
	  return index[i];
	}
	sortvec[index[i]-1] = i;
      }
    }

    idx_.resize(n_ind,0);
    val_.resize(n_ind,0);
    for (Index i=0; i<n_ind; ++i) {
      DBG_ASSERT(sortvec[i]>-1);
      idx_[i] = sortvec[i];
      val_[i]  = w;
    }

    Set_Initialized();
    Set_NRows((Index)val_.size());
    return 0;
  }

  void IndexSchurData::SetData_List(const std::vector<Index>& list, Number v)
  {
    DBG_START_METH("IndexSchurData::SetData_List", dbg_verbosity);

    DBG_ASSERT(!Is_Initialized());
    DBG_ASSERT(idx_.empty());
    DBG_ASSERT(v!=0);

    Index w;
    (v>0) ? w=1 : w=-1;

    val_.resize(list.size(), w);
    idx_ = list;

    Set_Initialized();
  }

  void IndexSchurData::GetRow(Index row, IteratesVector& v) const
  {
    DBG_START_METH("IndexSchurData::GetRow", dbg_verbosity);
    DBG_ASSERT(Is_Initialized());
    DBG_ASSERT(row<GetNRowsAdded());

    // retrieve structure of IteratesVector - this should probably be cached or sth.
    //Index n_comps = v.NComps();
    Index* v_lens = GetVectorLengths(v);
    // set vector v to 0
    v.Set(0.0);

    // find the vector and index in iteratesvector to which idx_[row] corresponds
    Index col = idx_[row];

    Index vec_idx = 0;
    while(!(col<v_lens[vec_idx])) {
      vec_idx++;
    }

    dynamic_cast<DenseVector*>(GetRawPtr(v.GetCompNonConst(vec_idx)))->Values()[col+v.GetComp(vec_idx)->Dim()-v_lens[vec_idx]] = (Number)val_[row];

    delete[] v_lens;
  }

  void IndexSchurData::GetMultiplyingVectors(Index row, std::vector<Index>& indices, std::vector<Number>& factors) const
  {
    DBG_START_METH("IndexSchurData::GetMultiplyingVectors", dbg_verbosity);

    DBG_ASSERT(indices.size()==0);
    DBG_ASSERT(factors.size()==0);

    indices.push_back(idx_[row]);
    factors.push_back(val_[row]);
  }

  void IndexSchurData::Multiply(const IteratesVector& v, Vector& u) const
  {
    DBG_START_METH("IndexSchurData::Multiply", dbg_verbosity);

    // this is awful.
    DenseVector* du = static_cast<DenseVector*>(&u);
    du->Set(0.0);
    Number* u_val = du->Values();

    Index* v_lens = GetVectorLengths(v);

    Index v_row, vec_idx;
    for (unsigned int i=0; i<idx_.size(); ++i) {
      v_row = idx_[i];

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

    delete [] v_lens;
  }

  void IndexSchurData::TransMultiply(const Vector& u, IteratesVector& v) const
  {
    DBG_START_METH("IndexSchurData::TransMultiply", dbg_verbosity);
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
    for (unsigned int i=0; i<idx_.size(); ++i) {
      row = i;
      col = idx_[i];
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

    delete [] v_vals;
  }

  Index* IndexSchurData::GetVectorLengths(const IteratesVector& v) const
  {
    DBG_START_METH("IndexSchurData::GetVectorLengths", dbg_verbosity);
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

  void IndexSchurData::PrintImpl(const Journalist& jnlst,
				 EJournalLevel level,
				 EJournalCategory category,
				 const std::string& name,
				 Index indent,
				 const std::string& prefix) const
  {
    DBG_START_METH("IndexSchurData::PrintImpl", dbg_verbosity);

    jnlst.PrintfIndented(level, category, indent,
                         "%sIndexSchurData \"%s\" with %d rows:\n",
                         prefix.c_str(), name.c_str(), GetNRowsAdded());
    if (Is_Initialized()) {
      for (unsigned int i=0; i<idx_.size(); i++) {
	jnlst.PrintfIndented(level, category, indent,
			     "%s%s[%5d,%5d]=%d\n",
			     prefix.c_str(), name.c_str(), i, idx_[i], val_[i]);
      }
    }
    else {
      jnlst.PrintfIndented(level, category, indent,
                           "%sUninitialized!\n",
                           prefix.c_str());
    }
  }

  void IndexSchurData::AddData_Flag(Index dim, Index* flags, std::vector<Index>& delta_u_sort, Index v)
  {
    DBG_START_METH("IndexSchurData::AddData_Flag", dbg_verbosity);

    Index sortcounter = (Index)idx_.size();
    bool oldindex;
    for (Index i=0; i<dim; ++i) {
      if (flags[i]) {
	oldindex = false;
	for (unsigned int j=0; j<idx_.size(); ++j) {
	  if (i==idx_[j]) {
	    delta_u_sort.push_back(j);
	    val_[j] = v;
	    oldindex = true;
	    break;
	  }
	}
	if (!oldindex) {
	  delta_u_sort.push_back(sortcounter++);
	  idx_.push_back(i);
	  val_.push_back(v);
	}
      }
    }
  }

  void IndexSchurData::AddData_List(std::vector<Index> cols, std::vector<Index>& delta_u_sort, Index& new_du_size, Index v)
  {
    DBG_START_METH("IndexSchurData::AddData_List", dbg_verbosity);

    new_du_size = (Index)idx_.size();
    bool oldindex;
    for (unsigned int i=0; i<cols.size(); ++i) {
      oldindex = false;
      for (unsigned int j=0; j<idx_.size(); ++j) {
	if (cols[i]==idx_[j]) {
	  delta_u_sort.push_back(j);
	  val_[j] = v;
	  oldindex = true;
	  break;
	}
      }
      if (!oldindex) {
	delta_u_sort.push_back(new_du_size++);
	idx_.push_back(cols[i]);
	val_.push_back(v);
      }
    }
    Set_NRows((Index)idx_.size());
    if(!Is_Initialized()) {
      Set_Initialized();
    }
  }

  Index IndexSchurData::GetNRowsAdded() const
  {
    DBG_START_METH("IndexSchurData::GetNRowsAdded", dbg_verbosity);

    return (Index)idx_.size();
  }

  const std::vector<Index>* IndexSchurData::GetColIndices() const
  {
    DBG_START_METH("IndexSchurData::GetColIndices", dbg_verbosity);
    return &idx_;
  }
}
