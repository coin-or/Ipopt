// Copyright (C) 2005, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-13

#include "IpTripletToCSRConverter.hpp"
#include <vector>
#include <algorithm>

#ifdef HAVE_CSTDDEF
# include <cstddef>
#else
# ifdef HAVE_STDDEF_H
#  include <stddef.h>
# else
#  error "don't have header file for stddef"
# endif
#endif

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  TripletToCSRConverter::
  TripletToCSRConverter(Index offset, ETriFull hf /*= Triangular_Format*/)
      :
      offset_(offset),
      hf_(hf),
      ia_(NULL),
      ja_(NULL),
      dim_(0),
      nonzeros_triplet_(0),
      nonzeros_compressed_(0),
      initialized_(false),
      ipos_first_(NULL),
      ipos_double_triplet_(NULL),
      ipos_double_compressed_(NULL)
  {
    DBG_ASSERT(offset==0|| offset==1);
  }

  TripletToCSRConverter::~TripletToCSRConverter()
  {
    delete[] ia_;
    delete[] ja_;
    delete[] ipos_first_;
    delete[] ipos_double_triplet_;
    delete[] ipos_double_compressed_;
  }

  Index TripletToCSRConverter::InitializeConverter(Index dim, Index nonzeros,
      const Index* airn,
      const Index* ajcn)
  {
    DBG_START_METH("TSymLinearSolver::InitializeStructure",
                   dbg_verbosity);

    DBG_ASSERT(dim>0);
    DBG_ASSERT(nonzeros>0);

    delete[] ia_;
    delete[] ja_;
    delete[] ipos_first_;
    delete[] ipos_double_triplet_;
    delete[] ipos_double_compressed_;

    dim_ = dim;
    nonzeros_triplet_ = nonzeros;

    // Create a list with all triplet entries
    std::vector<TripletEntry> entry_list(nonzeros);
    std::vector<TripletEntry>::iterator list_iterator = entry_list.begin();
    for (Index i=0; i<nonzeros; i++) {
      list_iterator->Set(airn[i], ajcn[i], i);
      list_iterator++;
    }
    DBG_ASSERT(list_iterator == entry_list.end());

    if (DBG_VERBOSITY()>=2) {
      for (Index i=0; i<nonzeros; i++) {
        DBG_PRINT((2, "airn[%5d] = %5d acjn[%5d] = %5d\n", i, airn[i], i, ajcn[i]));
      }
    }

    // sort the list
    std::sort(entry_list.begin(), entry_list.end());

    // Now got through the list and compute ipos_ arrays and the
    // number of elements in the compressed format
    Index* ja_tmp = new Index[nonzeros];    // overestimate memory requirement
    Index* rc_tmp = NULL;
    if (hf_ == Full_Format) {
      rc_tmp = new Index[dim_+1];
    }
    ia_ = new Index[dim_+1];
    Index* ipos_first_tmp = new Index[nonzeros];  // overestimate memory requirement
    Index* ipos_double_triplet_tmp = new Index[nonzeros];  // overestimate memory requirement
    Index* ipos_double_compressed_tmp = new Index[nonzeros];  // overestimate memory requirement

    Index nonzeros_compressed_full = 0;
    nonzeros_compressed_ = 0;
    Index cur_row = 1;

    if (hf_ == Full_Format) {
      // Zero row counts
      for (Index i=0; i<dim_+1; i++) {
        rc_tmp[i] = 0;
      }
    }

    // Take care of possible empty rows
    list_iterator = entry_list.begin();
    while (cur_row < list_iterator->IRow()) {
      ia_[cur_row-1] = 0;
      cur_row++;
    }
    ia_[cur_row-1] = 0;
    ja_tmp[0] = list_iterator->JCol();
    ipos_first_tmp[0] = list_iterator->PosTriplet();
    if (hf_ == Full_Format) {
      // Count in both lower and upper triangles. Count diagonal only once.
      nonzeros_compressed_full++;
      rc_tmp[cur_row-1]++;
      if (cur_row!=list_iterator->JCol()) {
        nonzeros_compressed_full++;
        rc_tmp[list_iterator->JCol()-1]++;
      }
    }

    list_iterator++;
    Index idouble = 0;
    Index idouble_full = 0;
    while (list_iterator != entry_list.end()) {
      Index irow = list_iterator->IRow();
      Index jcol = list_iterator->JCol();
      if (cur_row == irow && ja_tmp[nonzeros_compressed_] == jcol) {
        // This element appears repeatedly, add to the double list
        ipos_double_triplet_tmp[idouble] = list_iterator->PosTriplet();
        ipos_double_compressed_tmp[idouble] = nonzeros_compressed_;
        idouble++;
        idouble_full++;
        if (hf_==Full_Format && irow!=jcol) {
          idouble_full++;
        }
      }
      else {
        // This is a new element
        if (hf_==Full_Format) {
          // Count in both lower and upper triangles. Count diagonal only once.
          nonzeros_compressed_full++;
          rc_tmp[jcol-1]++;
          if (irow!=jcol) {
            nonzeros_compressed_full++;
            rc_tmp[irow-1]++;
          }
        }
        nonzeros_compressed_++;
        ja_tmp[nonzeros_compressed_] = jcol;
        ipos_first_tmp[nonzeros_compressed_] = list_iterator->PosTriplet();
        if (cur_row != irow) {
          // this is in a new row

          ia_[cur_row] = nonzeros_compressed_;
          cur_row++;
        }
      }

      list_iterator++;
    }
    nonzeros_compressed_++;
    for (Index i=cur_row; i<=dim_; i++) {
      ia_[i] = nonzeros_compressed_;
    }
    DBG_ASSERT(idouble == nonzeros_triplet_-nonzeros_compressed_);

    // Now copy the ja_tmp array to the (shorter) final one and make
    // sure that the correct offset is used
    if (hf_==Triangular_Format) {
      ja_ = new Index[nonzeros_compressed_];
      if (offset_==0) {
        for (Index i=0; i<nonzeros_compressed_; i++) {
          ja_[i] = ja_tmp[i] - 1;
        }
      }
      else {
        for (Index i=0; i<nonzeros_compressed_; i++) {
          ja_[i] = ja_tmp[i];
        }
        for (Index i=0; i<=dim_; i++) {
          ia_[i] = ia_[i] + 1;
        }
      }
      delete[] ja_tmp;

      // Reallocate memory for the "first" array
      ipos_first_ = new Index[nonzeros_compressed_];
      for (Index i=0; i<nonzeros_compressed_; i++) {
        ipos_first_[i] = ipos_first_tmp[i];
      }
      delete[] ipos_first_tmp;

      // Reallocate memory for the "double" arrays
      ipos_double_triplet_ = new Index[idouble];
      ipos_double_compressed_ = new Index[idouble];
      for (Index i=0; i<idouble; i++) {
        ipos_double_triplet_[i] = ipos_double_triplet_tmp[i];
        ipos_double_compressed_[i] = ipos_double_compressed_tmp[i];
      }
      delete[] ipos_double_triplet_tmp;
      delete[] ipos_double_compressed_tmp;
      num_doubles_ = nonzeros_triplet_ - nonzeros_compressed_;
    }
    else { // hf_==Full_Format

      // Setup ia_tmp to contain insert position for column i as ia_tmp[i+1]
      Index *ia_tmp = new Index[dim_+1];
      ia_tmp[0] = 0;
      ia_tmp[1] = 0;
      for (Index i=1; i<dim_; i++) {
        ia_tmp[i+1] = ia_tmp[i] + rc_tmp[i-1];
      }
      delete[] rc_tmp;

      // Loop over elements of matrix, copying them and duplicating as required
      ja_ = new Index[nonzeros_compressed_full];
      ipos_first_ = new Index[nonzeros_compressed_full];
      ipos_double_triplet_ = new Index[idouble_full];
      ipos_double_compressed_ = new Index[idouble_full];
      Index jd1=0; // Entry into ipos_double_compressed_tmp
      Index jd2=0; // Entry into ipos_double_compressed_
      for (Index i=0; i<dim_; i++) {
        for (Index j=ia_[i]; j<ia_[i+1]; j++) {
          Index jrow = ja_tmp[j]-1;
          ja_[ia_tmp[i+1]] = jrow + offset_;
          ipos_first_[ia_tmp[i+1]] = ipos_first_tmp[j];
          while (jd1<idouble && j==ipos_double_compressed_tmp[jd1]) {
            ipos_double_triplet_[jd2] = ipos_double_triplet_tmp[jd1];
            ipos_double_compressed_[jd2] = ia_tmp[i+1];
            jd2++;
            if (jrow!=i) {
              ipos_double_triplet_[jd2] = ipos_double_triplet_tmp[jd1];
              ipos_double_compressed_[jd2] = ia_tmp[jrow+1];
              jd2++;
            }
            jd1++;
          }
          ia_tmp[i+1]++;
          if (jrow!=i) {
            ja_[ia_tmp[jrow+1]] = i + offset_;
            ipos_first_[ia_tmp[jrow+1]] = ipos_first_tmp[j];
            ia_tmp[jrow+1]++;
          }
        }
      }
      delete[] ja_tmp;
      delete[] ipos_first_tmp;
      delete[] ipos_double_triplet_tmp;
      delete[] ipos_double_compressed_tmp;

      // Copy ia_tmp to ia_ with offset_ adjustment
      for (Index i=0; i<dim_+1; i++) {
        ia_[i] = ia_tmp[i] + offset_;
      }
      delete[] ia_tmp;

      // Set nonzeros_compressed_ to correct size
      nonzeros_compressed_ = nonzeros_compressed_full;
      num_doubles_ = idouble_full;
    }

    initialized_ = true;

    if (DBG_VERBOSITY()>=2) {
      for (Index i=0; i<=dim_; i++) {
        DBG_PRINT((2, "ia[%5d] = %5d\n", i, ia_[i]));
      }
      for (Index i=0; i<nonzeros_compressed_; i++) {
        DBG_PRINT((2, "ja[%5d] = %5d ipos_first[%5d] = %5d\n", i, ja_[i], i, ipos_first_[i]));
      }
      for (Index i=0; i<nonzeros_triplet_-nonzeros_compressed_; i++) {
        DBG_PRINT((2, "ipos_double_triplet[%5d] = %5d ipos_double_compressed[%5d] = %5d\n", i, ipos_double_triplet_[i], i, ipos_double_compressed_[i]));
      }
    }

    return nonzeros_compressed_;
  }

  void TripletToCSRConverter::ConvertValues(Index nonzeros_triplet,
      const Number* a_triplet,
      Index nonzeros_compressed,
      Number* a_compressed)
  {
    DBG_START_METH("TSymLinearSolver::ConvertValues",
                   dbg_verbosity);

    DBG_ASSERT(initialized_);

    DBG_ASSERT(nonzeros_triplet_==nonzeros_triplet);
    DBG_ASSERT(nonzeros_compressed_==nonzeros_compressed);

    for (Index i=0; i<nonzeros_compressed_; i++) {
      a_compressed[i] = a_triplet[ipos_first_[i]];
    }
    for (Index i=0; i<num_doubles_; i++) {
      a_compressed[ipos_double_compressed_[i]] +=
        a_triplet[ipos_double_triplet_[i]];
    }

    if (DBG_VERBOSITY()>=2) {
      for (Index i=0; i<nonzeros_triplet; i++) {
        DBG_PRINT((2, "atriplet[%5d] = %24.16e\n", i, a_triplet[i]));
      }
      for (Index i=0; i<nonzeros_compressed; i++) {
        DBG_PRINT((2, "acompre[%5d] = %24.16e\n", i, a_compressed[i]));
      }
    }
  }

} // namespace Ipopt
