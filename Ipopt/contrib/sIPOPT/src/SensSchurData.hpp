// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-06

#ifndef __ASSCHURDATA_HPP__
#define __ASSCHURDATA_HPP__

#include "IpVector.hpp"
#include "IpIteratesVector.hpp"
#include <vector>

namespace Ipopt
{


  class SchurData : public ReferencedObject
  {
    /** This interface serves as a reference point for multiple classes
     *  that need to use SchurData (PCalculator, SchurDriver). It
     *  declares as little as possible, so that SchurData implementations
     *  can be very special and fast.
     *
     *  I have not decided yet if there are certain ways I want to impose
     *  that SchurData can be set. I will figure this out as soon as I
     *  write the upstream classes that need to do that
     *
     *  Nomenclature in this program is based on Victor Zavalas thesis. */

  public:

    SchurData()  : initialized_(false), nrows_(0)
    {}

    virtual ~SchurData()
    {
    }

    virtual SmartPtr<SchurData> MakeNewSchurDataCopy() const =0;

    /** Functions to set the Schurdata. At least one must be overloaded */

    /** Set Data to one for given indices. Size of vector is ipopt_x_<full_x_ */
    virtual void SetData_Flag(Index dim, const Index* flags, Number v=1.0)=0;

    /** Set Data to corresponing Number */
    virtual void SetData_Flag(Index dim, const Index* flags, const Number* values)=0;

    virtual Index SetData_Index(Index dim, const Index* flags, Number v=1.0)=0;

    virtual void SetData_List(const std::vector<Index>& list, Number v=1.0) =0;

    virtual void AddData_List(std::vector<Index> cols, std::vector<Index>& delta_u_sort, Index& new_du_size, Index v)=0;

    /** Returns number of rows/columns in schur matrix */
    virtual Index GetNRowsAdded() const
    {
      return nrows_;
    }

    virtual bool Is_Initialized() const
    {
      return initialized_;
    }


    /** Returns the i-th column vector of the matrix */
    virtual void GetRow(Index i, IteratesVector& v) const = 0;

    /** Returns two vectors that are needed for matrix-vector
     *  multiplication of B and P.
     *  The index is the row, the first vector are the indices
     *  of non-zero components, in this row of B,
     *  the second vector gives the numbers in B(row,indices) */
    virtual void GetMultiplyingVectors(Index row, std::vector<Index>& indices, std::vector<Number>& factors) const =0;

    /** Computes B*v with B in R(mxn) */
    virtual void Multiply(const IteratesVector& v, Vector& u) const =0;

    /** Computes A*u with A in R(nxm), KKT in R(n,n) */
    virtual void TransMultiply(const Vector& u, IteratesVector& v) const =0;

    virtual void PrintImpl(const Journalist& jnlst,
			   EJournalLevel level,
			   EJournalCategory category,
			   const std::string& name,
			   Index indent,
			   const std::string& prefix) const =0;

    void Print(const Journalist& jnlst,
	       EJournalLevel level,
	       EJournalCategory category,
	       const std::string& name,
	       Index indent=0,
	       const std::string& prefix="") const
    {
      if (jnlst.ProduceOutput(level, category)) {
	PrintImpl(jnlst, level, category, name, indent, prefix);
      }
    }

    void Print(SmartPtr<const Journalist> jnlst,
	       EJournalLevel level,
	       EJournalCategory category,
	       const std::string& name,
	       Index indent,
	       const std::string& prefix) const
    {
      if (IsValid(jnlst) && jnlst->ProduceOutput(level, category)) {
	PrintImpl(*jnlst, level, category, name, indent, prefix);
      }
    }

  protected:

    virtual void Set_Initialized()
    {
      initialized_ = true;
    }

    virtual void Set_NRows(Index nrows)
    {
      nrows_ = nrows;
    }

  private:

    /** Makes sure that data is not set twice accidentially */
    bool initialized_;

    /** Number of columns/rows of corresponding Schur Matrix*/
    Index nrows_;

  };

}

#endif
