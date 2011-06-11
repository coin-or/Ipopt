// Copyright (C) 2009, Jonathan Hogg <jdh41.at.cantab.net>
// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors: Jonathan Hogg                           2009-07-29 
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#ifdef COIN_HAS_HSL
#include "CoinHslConfig.h"
#endif

#ifdef COINHSL_HAS_MA77

#include "IpMa77SolverInterface.hpp"
#include <iostream>
using namespace std;

extern "C" {
   /*
    * Easier to just have our own definition than include the full metis.h
    */
   extern void METIS_NodeND(int *n, int *xadj, int *adjncy,
      int *numflag, int *options, int *perm, int *iperm);

   /*
    * Initialise control parameters
    */
   extern void F77_FUNC (ma77_iface_set_control, MA77_IFACE_SET_CONTROL) (
         const ipfint* print_level,
         const ipfint* bits,
         const ipfint* buffer_lpage,
         const ipfint* buffer_npage,
         const ipfint* file_size,
         const ipfint* maxstore,
         const ipfint* nemin,
         const double* small,
         const double* static_,
         const double* u);
   /* 
    * Initialise solver, load pattern, call analyse
    */
   extern void F77_FUNC (ma77_iface_initstructure, MA77_IFACE_INITSTRUCTURE) (
         const ipfint* ndim,  /* order of matrix */
         const ipfint* nnz,   /* number of non-zeroes */
         const ipfint* ptr,   /* column starts */
         const ipfint* row,   /* row numbers */
         ipfint* order,       /* permutation */
         const ipfint* icntl, /* integer controls */
         const double* rcntl, /* float controls */
         ipfint* info);       /* return code from ma77 routines */
   /*
    * Load reals, call factorise
    */
   extern void F77_FUNC (ma77_iface_factor,MA77_IFACE_FACTOR) (
         const ipfint* ndim,     /* order of matrix */
         const ipfint* ptr,      /* column starts */
         const double* val,      /* matrix entries */
         const ipfint* icntl,    /* integer controls */
         const double* rcntl,    /* float controls */
         ipfint* numneg,         /* returns number of negative pivots */
         ipfint* info);          /* return code from ma77 routines */
   /*
    * Perform a solve with pre-calculated factors
    */
   extern void F77_FUNC (ma77_iface_solve, MA77_IFACE_SOLVE) (
         const ipfint* nrhs,     /* number of right-hand sides */
         const ipfint *lx,       /* leading edge of rhs */
         double* rhs,            /* values of rhs */
         const ipfint* icntl,    /* integer controls */
         const double* rcntl,    /* float controls */
         ipfint* info);          /* return code from ma77 routines */
   /*
    * Cleanup memory
    */
   extern void F77_FUNC (ma77_iface_finalise, MA77_IFACE_FINALISE) (
         const ipfint* icntl,    /* integer controls */
         const double* rcntl);   /* float controls */
}

namespace Ipopt
{

Ma77SolverInterface::~Ma77SolverInterface()
{
   delete [] val_;

   F77_FUNC (ma77_iface_finalise, MA77_IFACE_FINALISE)
      (icntl_, rcntl_);
}

void Ma77SolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
{
   roptions->AddIntegerOption(
      "ma77_print_level",
      "Debug printing level for the linear solver MA77",
      0,
      "Meep");
      /*
      "<0 no printing.\n"
      "0  Error and warning messages only.\n"
      "=1 Limited diagnostic printing.\n"
      ">1 Additional diagnostic printing.");
      */
   roptions->AddLowerBoundedIntegerOption(
      "ma77_buffer_lpage",
      "Number of scalars per MA77 buffer page",
      1, 4096,
      "Number of scalars per an in-core buffer in the out-of-core solver "
      "MA77. Must be at most ma77_file_size.");
   roptions->AddLowerBoundedIntegerOption(
      "ma77_buffer_npage",
      "Number of pages that make up MA77 buffer",
      1, 1600,
      "Number of pages of size buffer_lpage that exist in-core for the "
      "out-of-core solver MA77.");
   roptions->AddLowerBoundedIntegerOption(
      "ma77_file_size",
      "Target size of each temporary file for MA77, scalars per type",
      1, 2097152,
      "MA77 uses many temporary files, this option controls the size of "
      "each one. It is measured in the number of entries (int or double), "
      "NOT bytes.");
   roptions->AddLowerBoundedIntegerOption(
      "ma77_maxstore",
      "Maximum storage size for MA77 in-core mode",
      0, 0,
      "If greater than zero, the maximum size of factors stored in core "
      "before out-of-core mode is invoked.");
   roptions->AddLowerBoundedIntegerOption(
      "ma77_nemin",
      "Node Amalgamation parameter",
      1, 8,
      "Two nodes in elimination tree are merged if result has fewer than "
      "ma77_nemin variables.");
   roptions->AddLowerBoundedNumberOption(
      "ma77_small",
      "Zero Pivot Threshold",
      0.0, false, 1e-308,
      "Any pivot less than ma77_small is treated as zero.");
   roptions->AddLowerBoundedNumberOption(
      "ma77_static",
      "Static Pivoting Threshold",
      0.0, false, 0.0,
      "See MA77 documentation. Either ma77_static=0.0 or "
      "ma77_static>ma77_small. ma77_static=0.0 disables static pivoting.");
   roptions->AddBoundedNumberOption(
      "ma77_u",
      "Pivoting Threshold",
      0.0, false, 0.5, false, 0.01,
      "See MA77 documentation.");
}
  
bool Ma77SolverInterface::InitializeImpl(const OptionsList& options,
   const std::string& prefix)
{
   int bits=32;

   options.GetIntegerValue("ma77_print_level", ma77_print_level_, prefix);
   options.GetIntegerValue("ma77_buffer_lpage", ma77_buffer_lpage_, prefix);
   options.GetIntegerValue("ma77_buffer_npage", ma77_buffer_npage_, prefix);
   options.GetIntegerValue("ma77_file_size", ma77_file_size_, prefix);
   options.GetIntegerValue("ma77_maxstore", ma77_maxstore_, prefix);
   options.GetIntegerValue("ma77_nemin", ma77_nemin_, prefix);
   options.GetNumericValue("ma77_small", ma77_small_, prefix);
   options.GetNumericValue("ma77_static", ma77_static_, prefix);
   options.GetNumericValue("ma77_u", ma77_u_, prefix);

   icntl_[0] = ma77_print_level_;
   icntl_[1] = bits;
   icntl_[2] = ma77_buffer_lpage_;
   icntl_[3] = ma77_buffer_npage_;
   icntl_[4] = ma77_file_size_;
   icntl_[5] = ma77_maxstore_;
   icntl_[6] = ma77_nemin_;

   rcntl_[0] = ma77_small_;
   rcntl_[1] = ma77_static_;
   rcntl_[2] = ma77_u_;

   /*F77_FUNC (ma77_iface_set_control, MA77_IFACE_SET_CONTROL)
      (&ma77_print_level_, &bits, &ma77_buffer_lpage_,
       &ma77_buffer_npage_, &ma77_file_size_, &ma77_maxstore_, &ma77_nemin_,
       &ma77_small_, &ma77_static_, &ma77_u_);*/

   return true; // All is well
}

/*  Method for initializing internal stuctures.  Here, ndim gives
 *  the number of rows and columns of the matrix, nonzeros give
 *  the number of nonzero elements, and ia and ja give the
 *  positions of the nonzero elements, given in the matrix format
 *  determined by MatrixFormat.
 */
ESymSolverStatus Ma77SolverInterface::InitializeStructure(Index dim, 
   Index nonzeros, const Index* ia, const Index* ja)
{
   int info;

   // Store size for later use
   ndim_ = dim;

   // Determine an ordering
   Index *perm = new Index[dim];
   MetisOrder(dim, ia, ja, perm);
   //for(int i=0; i<dim; i++) perm[i] = i+1;

   // Setup memory for values
   if(val_!=NULL) delete[] val_;
   val_ = new double[nonzeros];

   // perform analyse
   F77_FUNC (ma77_iface_initstructure, MA77_IFACE_INITSTRUCTURE)
      (&dim, &nonzeros, ia, ja, perm, icntl_, rcntl_, &info);
   delete[] perm; // Done with order

   if(info>=0) {
      return SYMSOLVER_SUCCESS;
   } else {
      return SYMSOLVER_FATAL_ERROR;
   }
}

/*  Solve operation for multiple right hand sides.  Solves the
 *  linear system A * x = b with multiple right hand sides, where
 *  A is the symmtric indefinite matrix.  Here, ia and ja give the
 *  positions of the values (in the required matrix data format).
 *  The actual values of the matrix will have been given to this
 *  object by copying them into the array provided by
 *  GetValuesArrayPtr. ia and ja are identical to the ones given
 *  to InitializeStructure.  The flag new_matrix is set to true,
 *  if the values of the matrix has changed, and a refactorzation
 *  is required.
 *
 *  The return code is SYMSOLV_SUCCESS if the factorization and
 *  solves were successful, SYMSOLV_SINGULAR if the linear system
 *  is singular, and SYMSOLV_WRONG_INERTIA if check_NegEVals is
 *  true and the number of negative eigenvalues in the matrix does
 *  not match numberOfNegEVals.  If SYMSOLV_CALL_AGAIN is
 *  returned, then the calling function will request the pointer
 *  for the array for storing a again (with GetValuesPtr), write
 *  the values of the nonzero elements into it, and call this
 *  MultiSolve method again with the same right-hand sides.  (This
 *  can be done, for example, if the linear solver realized it
 *  does not have sufficient memory and needs to redo the
 *  factorization; e.g., for MA27.)
 *
 *  The number of right-hand sides is given by nrhs, the values of
 *  the right-hand sides are given in rhs_vals (one full right-hand
 *  side stored immediately after the other), and solutions are
 *  to be returned in the same array.
 *
 *  check_NegEVals will not be chosen true, if ProvidesInertia()
 *  returns false.
 */
ESymSolverStatus Ma77SolverInterface::MultiSolve(bool new_matrix,
   const Index* ia, const Index* ja, Index nrhs, double* rhs_vals,
   bool check_NegEVals, Index numberOfNegEVals)
{
   int info;

   if(new_matrix)
   {
      F77_FUNC (ma77_iface_factor, MA77_IFACE_FACTOR)
         (&ndim_, ia, val_, icntl_, rcntl_, &numneg_, &info);
      if(info<0) return SYMSOLVER_FATAL_ERROR;

      if(check_NegEVals && numneg_!=numberOfNegEVals)
         return SYMSOLVER_WRONG_INERTIA;
   }

   F77_FUNC (ma77_iface_solve, MA77_IFACE_SOLVE) 
      (&nrhs, &ndim_, rhs_vals, icntl_, rcntl_, &info);

   return SYMSOLVER_SUCCESS;
}

/*
 * Call metis_NodeND to perform ordering on the graph, return it in perm
 */
void Ma77SolverInterface::MetisOrder(const int ndim, const Index *ptr, 
   const Index *row, Index *perm)
{
   int options[8];
   options[0] = 0; // Defaults
   int numflag = 1;
   int ndim_nc = ndim;

   Index *ptr_tmp = new Index[ndim+1];
   Index *row_tmp = new Index[ptr[ndim]-1];
   ptr_tmp[0] = 1;
   for(int i=0; i<ndim; i++)
   {
      ptr_tmp[i+1] = ptr_tmp[i];
      for(int j=ptr[i]-1; j<ptr[i+1]-1; j++)
      {
         if(i==row[j]-1) continue; // Skip diagonals
         row_tmp[ptr_tmp[i+1]-1] = row[j];
         ptr_tmp[i+1]++;
      }
   }

   // Note that MeTiS's iperm is our perm and vice-versa
   Index *iperm = new Index[ndim];
   METIS_NodeND(&ndim_nc, ptr_tmp, row_tmp, &numflag, options, iperm, perm);
   delete[] iperm;
   delete[] row_tmp;
   delete[] ptr_tmp;
}

} // namespace Ipopt

#endif /* COINHSL_HAS_MA77 */
