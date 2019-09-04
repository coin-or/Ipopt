// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-08

#ifndef __ASINDEXSCHURDATA_HPP__
#define __ASINDEXSCHURDATA_HPP__

#include "SensSchurData.hpp"

namespace Ipopt
{

  class IndexSchurData : public SchurData
  {
    /** This class is the implementation aimed at applications where
     *  only SchurData matrices with entries 1 or -1 appear. */

  public:

    IndexSchurData();

    IndexSchurData(const std::vector<Index> idx, const std::vector<Index> val);

    virtual ~IndexSchurData();

    virtual SmartPtr<SchurData> MakeNewSchurDataCopy() const;

    virtual Index GetNRowsAdded() const;

    virtual void SetData_Flag(Index dim, const Index* flags, Number v=1.0);

    virtual void SetData_Flag(Index dim, const Index* flags, const Number* values);

    virtual Index SetData_Index(Index dim, const Index* index, Number v=1.0);

    virtual void SetData_List(const std::vector<Index>& list, Number v=1.0);

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

    /** Functions specific to IndexSchurData */

    /** This function is for adding data to a SchurData object. It takes a set of column-indices
     *  a value v and adds indices accordingly. If the column is already set in the data,
     *  it stays at the same place, otherwise the new indices are added at the bottom,
     *  in the order specified by the indices. The vector delta_u_sort returns the actual
     *  sorting so that the user knows how to place the new values inside the elongated
     *  delta_u vector. These places are in C++ index style, so they correspond exactly
     *  to the indices used for the C++-array of the delta_u DenseVector*/
    void AddData_Flag(Index dim, Index* flags, std::vector<Index>& delta_u_sort, Index v);

    void AddData_List(std::vector<Index> cols, std::vector<Index>& delta_u_sort, Index& new_du_size, Index v);

    const std::vector<Index>* GetColIndices() const;

  private:

    /** returns a vector that holds the accumulated length of each vector component:
     *  v_len[0] = v.GetComp(0)->Dim()
     *  v_len[i] = sum(k=0..i, v.GetComp(k)->Dim()) */
    Index* GetVectorLengths(const IteratesVector& v) const;

    std::vector<Index> idx_;
    std::vector<Index> val_;
  };

}

#endif
