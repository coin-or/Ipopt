// Copyright 2009, 2010 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-08

#ifndef __ASSTDSCHURDATA_HPP__
#define __ASSTDSCHURDATA_HPP__

#include "AsSchurData.hpp"

namespace Ipopt
{

  class StdSchurData : public SchurData
  {
    /** This class is the standard implementation of the SchurData interface.
     *  It uses std::vectors to save the data in lists and is build with 
     *  the goals of understandability, changeability and flexibility in 
     *  mind rather than performance. 
     *  Use this class only for academic purposes. If you can exploit a 
     *  special structure in your SchurData (e.g. all indices are ones),
     *  please do in a different implementation of the SchurData interface.*/

  public:

    StdSchurData();
    
    virtual ~StdSchurData();

    virtual SmartPtr<SchurData> MakeNewSchurDataCopy() const;

    virtual void SetData_Flag(Index dim, const Index* flags, Number v=1.0);

    virtual void SetData_Flag(Index dim, const Index* flags, const Number* values);

    virtual Index SetData_Index(Index dim, const Index* index, Number v=1.0);

    virtual void SetData_List(const std::vector<Index>& list, Number v=1.0);

    virtual void AddData_List(std::vector<Index> cols, std::vector<Index>& delta_u_sort, Index& new_du_size, Index v);

    virtual void GetRow(Index i, IteratesVector& v) const;

    virtual void GetMultiplyingVectors(Index i, std::vector<Index>& indices, std::vector<Number>& factors) const;

    virtual void Multiply(const IteratesVector& v, Vector& u) const;

    virtual void TransMultiply(const Vector& u, IteratesVector& v) const;

    virtual void PrintImpl(const Journalist& jnlst,
			   EJournalLevel level,
			   EJournalCategory category,
			   const std::string& name,
			   Index indent,
			   const std::string& prefix) const;
  private:

    /** returns a vector that holds the accumulated length of each vector component:
     *  v_len[0] = v.GetComp(0)->Dim()
     *  v_len[i] = sum(k=0..i, v.GetComp(k)->Dim()) */
    Index* GetVectorLengths(const IteratesVector& v) const;

    /** Row (maximal in 1:nrows_) */
    std::vector<Index> idx1_;
    /** Column (maximal in sum(i, curr()->GetComp(i)->Dim()) */
    std::vector<Index> idx2_;
    /** value at (row,col) */
    std::vector<Number> val_;
  };

}

#endif
