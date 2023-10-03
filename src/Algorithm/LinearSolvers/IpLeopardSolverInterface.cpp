// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Morten Lysgaard               2023-09-19
//               original version (based on MA57TSolverInterface.hpp)

#include "IpoptConfig.h"
#include "IpLeopardSolverInterface.hpp"

#include <cmath>
#include <iostream>

namespace Ipopt {
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
#endif

LeopardSolverInterface::LeopardSolverInterface(
) :
    dim_(0),
    nonzeros_(0),
    initialized_(false),
    pivtol_changed_(false),
    refactorize_(false) {
  DBG_START_METH("LeopardSolverInterface::LeopardSolverInterface()", dbg_verbosity);
}

LeopardSolverInterface::~LeopardSolverInterface() {
  DBG_START_METH("LeopardSolverInterface::~LeopardSolverInterface()", dbg_verbosity);
}

void LeopardSolverInterface::RegisterOptions(
    SmartPtr<RegisteredOptions> roptions
) {
  // TODO
}

bool LeopardSolverInterface::InitializeImpl(
    const OptionsList &options,
    const std::string &prefix
) {
  return true;
}

ESymSolverStatus LeopardSolverInterface::MultiSolve(
    bool new_matrix,
    const Index *ia,
    const Index *ja,
    Index nrhs,
    Number *rhs_vals,
    bool check_NegEVals,
    Index numberOfNegEVals
) {
  DBG_START_METH("LeopardSolverInterface::MultiSolve", dbg_verbosity);
   DBG_ASSERT(!check_NegEVals || ProvidesInertia());
   DBG_ASSERT(initialized_);

   // check if a factorization has to be done
   if( new_matrix )
   {
      // perform the factorization
      ESymSolverStatus retval;
      retval = Factorization(ia, ja, check_NegEVals, numberOfNegEVals);
      if( retval != SYMSOLVER_SUCCESS )
      {
         DBG_PRINT((1, "FACTORIZATION FAILED!\n"));
         return retval;  // Matrix singular or error occurred
      }
   }

   // do the solve
  return Backsolve(nrhs, rhs_vals);
}

Number *LeopardSolverInterface::GetValuesArrayPtr() {
  DBG_START_METH("LeopardSolverInterface::GetValuesArrayPtr", dbg_verbosity);
  DBG_ASSERT(initialized_);
  DBG_ASSERT(a_);

  return a_;
}

/** Initialize the local copy of the positions of the nonzero elements */
ESymSolverStatus LeopardSolverInterface::InitializeStructure(
    Index dim,
    Index nonzeros,
    const Index *ia,
    const Index *ja
) {
  DBG_START_METH("LeopardSolverInterface::InitializeStructure", dbg_verbosity);
  dim_ = dim;
  nonzeros_ = nonzeros;

  // Make space for storing the matrix elements
  if(!a_){
    delete[] a_;
  }
  a_ = NULL;
  a_ = new Number[nonzeros_];

  // Do the symbolic facotrization
  ESymSolverStatus retval = SymbolicFactorization(ia, ja);
  if( retval != SYMSOLVER_SUCCESS )
  {
     return retval;
  }

  initialized_ = true;

  return retval;
}

ESymSolverStatus LeopardSolverInterface::SymbolicFactorization(
    const Index *airn,
    const Index *ajcn
) {
  DBG_START_METH("LeopardSolverInterface::SymbolicFactorization", dbg_verbosity);

#ifdef IPOPT_INT64
  ij_ = leopard_ij_new_csr_array_i64(dim_, nonzeros_, airn, ajcn, true);
#else
  ij_ = leopard_ij_new_csr_array_i32(dim_, nonzeros_, airn, ajcn, true);
#endif
  if(ij_.code != LeopardReturnCode_Ok){
    printf("LeopardSolverInterface::SymbolicFactorization::ERROR: %s\n", leopard_explain_return_code(ij_.code));
    return SYMSOLVER_FATAL_ERROR;
  }
  DBG_ASSERT(ij_.object!=nullptr);

  ordering_ = leopard_ordering_new_amd_default(ij_.object);
  if(ordering_.code != LeopardReturnCode_Ok){
    printf("LeopardSolverInterface::SymbolicFactorization::ERROR: %s\n", leopard_explain_return_code(ordering_.code));
    return SYMSOLVER_FATAL_ERROR;
  }
  DBG_ASSERT(ordering_.object!=nullptr);

  ordered_ij_ = leopard_ordered_ij_new(ij_.object, ordering_.object);
  if(ordered_ij_.code != LeopardReturnCode_Ok) {
    printf("LeopardSolverInterface::SymbolicFactorization::ERROR: %s\n", leopard_explain_return_code(ordered_ij_.code));
    return SYMSOLVER_FATAL_ERROR;
  }
  DBG_ASSERT(ordered_ij_.object!=nullptr);

  assembly_tree_ = leopard_assembly_tree_new(ordered_ij_.object, true);
  if(assembly_tree_.code != LeopardReturnCode_Ok)
  {
    printf("LeopardSolverInterface::SymbolicFactorization::ERROR: %s\n", leopard_explain_return_code(assembly_tree_.code));
    return SYMSOLVER_FATAL_ERROR;
  }
  DBG_ASSERT(assembly_tree_.object!=nullptr);

  return SYMSOLVER_SUCCESS;
}

ESymSolverStatus LeopardSolverInterface::Factorization(
    const Index * /*ia*/,
    const Index * /*ja*/,
    bool check_NegEVals,
    Index numberOfNegEVals
) {
  DBG_START_METH("LeopardSolverInterface::Factorization", dbg_verbosity);

  DBG_ASSERT(ordered_ij_.object != nullptr);
  DBG_ASSERT(ordered_ij_.code == OrderedIjReturnCode::Ok);
  DBG_ASSERT(ordering_.object != nullptr);
  DBG_ASSERT(ordering_.code == OrderingReturnCode::Ok);
  DBG_ASSERT(assembly_tree_.object != nullptr);
  DBG_ASSERT(assembly_tree_.code == AssemblyTreeReturnCode::Ok);

  // Do factorization
#ifdef IPOPT_SINGLE
  ldl_factorization_ = leopard_ldl_factorize_f32(ordered_ij_.object, a_, ordering_.object, assembly_tree_.object);
#else
  ldl_factorization_ = leopard_ldl_factorize_f64(ordered_ij_.object, a_, ordering_.object, assembly_tree_.object);
#endif
  if(ldl_factorization_.code!=LeopardReturnCode_Ok)
  {
    printf("LeopardSolverInterface::Factorization::ERROR: %s\n", leopard_explain_return_code(ldl_factorization_.code));
    return SYMSOLVER_FATAL_ERROR;
  }
  DBG_ASSERT(ldl.object != nullptr);

  // Compute matrix intertia
#ifdef IPOPT_SINGLE
  eigenvalue_sign_count_ = leopard_ldl_eigenvalue_sign_count_f32(ldl_factorization_.object);
#else
  eigenvalue_sign_count_ = leopard_ldl_eigenvalue_sign_count_f64(ldl_factorization_.object);
#endif

  // Check if we have correct intertia
  // TODO: Check if this is corrent way to do it
  if(check_NegEVals && eigenvalue_sign_count_.negative!=numberOfNegEVals){
    return SYMSOLVER_WRONG_INERTIA;
  }

  return SYMSOLVER_SUCCESS;
}

ESymSolverStatus LeopardSolverInterface::Backsolve(
    Index nrhs,
    Number *rhs_vals
) {
  DBG_START_METH("LeopardSolverInterface::Backsolve", dbg_verbosity);
  DBG_ASSERT(ldl_factorization_.object != nullptr);
  DBG_ASSERT(ldl_factorization_.code == LDLFactorizationReturnCode::Ok);

  // TODO
  for (int i=0; i<nrhs; i++){
    Number *cur_rhs = rhs_vals+i*dim_;
#ifdef IPOPT_SINGLE
    LeopardReturnCode retval = leopard_ldl_solve_f32(ldl_factorization_.object, cur_rhs);
#else
    LeopardReturnCode retval = leopard_ldl_solve_f64(ldl_factorization_.object, cur_rhs);
#endif
    if(retval!=LeopardReturnCode_Ok)
    {
      return SYMSOLVER_FATAL_ERROR;
    }
  }
  return SYMSOLVER_SUCCESS;
}

Index LeopardSolverInterface::NumberOfNegEVals() const {
  DBG_START_METH("LeopardSolverInterface::NumberOfNegEVals", dbg_verbosity);
  DBG_ASSERT(ProvidesInertia());
  DBG_ASSERT(initialized_);
  return eigenvalue_sign_count_.negative;
}

bool LeopardSolverInterface::IncreaseQuality() {
  DBG_START_METH("LeopardSolverInterface::IncreaseQuality", dbg_verbosity);
  // TODO
  if (pivtol_ == pivtolmax_) {
    return false;
  }
  pivtol_changed_ = true;

  Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                 "Increasing pivot tolerance for Leopard from %7.2e ", pivtol_);
  pivtol_ = Min(pivtolmax_, std::pow(pivtol_, Number(0.75)));
  Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                 "to %7.2e.\n", pivtol_);
  return true;
}

} // namespace Ipopt
