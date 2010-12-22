// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter     IBM    2007-03-01

#include "IpGenAugSystemSolver.hpp"
#include "IpTripletHelper.hpp"
#include "IpDenseVector.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  GenAugSystemSolver::GenAugSystemSolver(GenKKTSolverInterface& SolverInterface)
      :
      AugSystemSolver(),
      solver_interface_(&SolverInterface),
      w_tag_(0),
      d_x_tag_(0),
      delta_x_(0.),
      d_s_tag_(0),
      delta_s_(0.),
      j_c_tag_(0),
      d_c_tag_(0),
      delta_c_(0.),
      j_d_tag_(0),
      d_d_tag_(0),
      delta_d_(0.),
      dx_vals_copy_(NULL),
      ds_vals_copy_(NULL),
      dc_vals_copy_(NULL),
      dd_vals_copy_(NULL)
  {
    DBG_START_METH("GenAugSystemSolver::GenAugSystemSolver()",dbg_verbosity);
    DBG_ASSERT(IsValid(solver_interface_));
  }

  GenAugSystemSolver::~GenAugSystemSolver()
  {
    DBG_START_METH("GenAugSystemSolver::~GenAugSystemSolver()",dbg_verbosity);
    delete [] dx_vals_copy_;
    delete [] ds_vals_copy_;
    delete [] dc_vals_copy_;
    delete [] dd_vals_copy_;
  }


  bool GenAugSystemSolver::InitializeImpl(const OptionsList& options,
                                          const std::string& prefix)
  {
    // This option is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);

    if (!warm_start_same_structure_) {
      delete [] dx_vals_copy_;
      delete [] ds_vals_copy_;
      delete [] dc_vals_copy_;
      delete [] dd_vals_copy_;
    }

    return solver_interface_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                         options, prefix);
  }

  ESymSolverStatus GenAugSystemSolver::MultiSolve(
    const SymMatrix* W,
    double W_factor,
    const Vector* D_x,
    double delta_x,
    const Vector* D_s,
    double delta_s,
    const Matrix* J_c,
    const Vector* D_c,
    double delta_c,
    const Matrix* J_d,
    const Vector* D_d,
    double delta_d,
    std::vector<SmartPtr<const Vector> >& rhs_xV,
    std::vector<SmartPtr<const Vector> >& rhs_sV,
    std::vector<SmartPtr<const Vector> >& rhs_cV,
    std::vector<SmartPtr<const Vector> >& rhs_dV,
    std::vector<SmartPtr<Vector> >& sol_xV,
    std::vector<SmartPtr<Vector> >& sol_sV,
    std::vector<SmartPtr<Vector> >& sol_cV,
    std::vector<SmartPtr<Vector> >& sol_dV,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {
    DBG_START_METH("GenAugSystemSolver::MultiSolve",dbg_verbosity);
    DBG_ASSERT(J_c && J_d && "Currently, you MUST specify J_c and J_d in the augmented system");

    DBG_ASSERT(W_factor == 0.0 || W_factor == 1.0);

    const Index nrhs = (Index)rhs_xV.size();
    DBG_ASSERT(nrhs>0);
    DBG_ASSERT(nrhs==(Index)rhs_sV.size());
    DBG_ASSERT(nrhs==(Index)rhs_cV.size());
    DBG_ASSERT(nrhs==(Index)rhs_dV.size());
    DBG_ASSERT(nrhs==(Index)sol_xV.size());
    DBG_ASSERT(nrhs==(Index)sol_sV.size());
    DBG_ASSERT(nrhs==(Index)sol_cV.size());
    DBG_ASSERT(nrhs==(Index)sol_dV.size());

    // Check if the input data has changed:
    bool new_matrix =
      AugmentedSystemChanged(W, W_factor, D_x, delta_x, D_s, delta_s,
                             *J_c, D_c, delta_c, *J_d, D_d, delta_d);

    // Get the individual arrays to be given to the
    // GenKKTSolverInterface
    const Index n_x = rhs_xV[0]->Dim();
    const Index n_c = rhs_cV[0]->Dim();
    const Index n_d = rhs_dV[0]->Dim();
    const Number* dx_vals=NULL;
    if (D_x) {
      const DenseVector* dD_x = dynamic_cast<const DenseVector*> (D_x);
      if (dD_x && !dD_x->IsHomogeneous()) {
        dx_vals = dD_x->Values();
      }
      else if (D_x->GetTag() != d_x_tag_) {
        delete [] dx_vals_copy_;
        dx_vals_copy_ = new Number[n_x];
        TripletHelper::FillValuesFromVector(n_x, *D_x, dx_vals_copy_);
        dx_vals = dx_vals_copy_;
      }
    }
    const Number* ds_vals=NULL;
    if (D_s) {
      const DenseVector* dD_s = dynamic_cast<const DenseVector*> (D_s);
      if (dD_s && !dD_s->IsHomogeneous()) {
        ds_vals = dD_s->Values();
      }
      else if (D_s->GetTag() != d_s_tag_) {
        delete [] ds_vals_copy_;
        ds_vals_copy_ = new Number[n_d];
        TripletHelper::FillValuesFromVector(n_d, *D_s, ds_vals_copy_);
        ds_vals = ds_vals_copy_;
      }
    }
    const Number* dc_vals=NULL;
    if (D_c) {
      const DenseVector* dD_c = dynamic_cast<const DenseVector*> (D_c);
      if (dD_c && !dD_c->IsHomogeneous()) {
        dc_vals = dD_c->Values();
      }
      else if (D_c->GetTag() != d_c_tag_) {
        delete [] dc_vals_copy_;
        dc_vals_copy_ = new Number[n_c];
        TripletHelper::FillValuesFromVector(n_c, *D_c, dc_vals_copy_);
        dc_vals = dc_vals_copy_;
      }
    }
    const Number* dd_vals=NULL;
    if (D_d) {
      const DenseVector* dD_d = dynamic_cast<const DenseVector*> (D_d);
      if (dD_d && !dD_d->IsHomogeneous()) {
        dd_vals = dD_d->Values();
      }
      else if (D_d->GetTag() != d_d_tag_) {
        delete [] dd_vals_copy_;
        dd_vals_copy_ = new Number[n_d];
        TripletHelper::FillValuesFromVector(n_d, *D_d, dd_vals_copy_);
        dd_vals = dd_vals_copy_;
      }
    }

    const Index dim = n_x+n_d+n_c+n_d;
    Number* rhssol = new Number[nrhs*dim];
    for (Index irhs=0; irhs<nrhs; irhs++) {
      // TODO: make order an option
      TripletHelper::FillValuesFromVector(n_x, *rhs_xV[irhs],
                                          &rhssol[irhs*dim]);
      TripletHelper::FillValuesFromVector(n_c, *rhs_cV[irhs],
                                          &rhssol[irhs*dim+n_x]);
      TripletHelper::FillValuesFromVector(n_d, *rhs_dV[irhs],
                                          &rhssol[irhs*dim+n_x+n_c]);
      TripletHelper::FillValuesFromVector(n_d, *rhs_sV[irhs],
                                          &rhssol[irhs*dim+n_x+n_c+n_d]);
    }

    bool done = false;
    ESymSolverStatus retval = SYMSOLVER_FATAL_ERROR;
    const SymMatrix* Wgive = NULL;
    if (W && W_factor==1.0) {
      Wgive = W;
    }
    while (!done) {
      retval = solver_interface_->MultiSolve(new_matrix, n_x, n_c, n_d,
                                             Wgive, J_c, J_d,
                                             dx_vals, ds_vals,
                                             dc_vals, dd_vals,
                                             delta_x, delta_s,
                                             delta_c, delta_d,
                                             nrhs, rhssol,
                                             check_NegEVals, numberOfNegEVals);
      if (retval==SYMSOLVER_CALL_AGAIN) {
        DBG_PRINT((1, "Solver interface asks to be called again.  Don't really se why...?\n"));
      }
      else {
        done = true;
      }
    }

    // Copy the values back into the vectors
    if (retval==SYMSOLVER_SUCCESS) {
      for (Index irhs=0; irhs<nrhs; irhs++) {
        TripletHelper::PutValuesInVector(n_x, &rhssol[irhs*dim],
                                         *sol_xV[irhs]);
        TripletHelper::PutValuesInVector(n_c, &rhssol[irhs*dim+n_x],
                                         *sol_cV[irhs]);
        TripletHelper::PutValuesInVector(n_d, &rhssol[irhs*dim+n_x+n_c],
                                         *sol_dV[irhs]);
        TripletHelper::PutValuesInVector(n_d, &rhssol[irhs*dim+n_x+n_c+n_d],
                                         *sol_sV[irhs]);
      }
    }
    else if (retval==SYMSOLVER_FATAL_ERROR) {
      delete [] rhssol;
      THROW_EXCEPTION(FATAL_ERROR_IN_LINEAR_SOLVER,"A fatal error occured in the linear solver.");
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Factorization failed with retval = %d\n", retval);
    }

    delete [] rhssol;

    return retval;
  }

  void GenAugSystemSolver::UpdateTags(const SymMatrix* W,
                                      double W_factor,
                                      const Vector* D_x,
                                      double delta_x,
                                      const Vector* D_s,
                                      double delta_s,
                                      const Matrix& J_c,
                                      const Vector* D_c,
                                      double delta_c,
                                      const Matrix& J_d,
                                      const Vector* D_d,
                                      double delta_d)
  {
    if (W) {
      w_tag_ = W->GetTag();
    }
    else {
      w_tag_ = 0;
    }
    w_factor_ = W_factor;

    if (D_x) {
      d_x_tag_ = D_x->GetTag();
    }
    else {
      d_x_tag_ = 0;
    }
    delta_x_ = delta_x;
    if (D_s) {
      d_s_tag_ = D_s->GetTag();
    }
    else {
      d_s_tag_ = 0;
    }
    delta_s_ = delta_s;
    if (D_c) {
      d_c_tag_ = D_c->GetTag();
    }
    else {
      d_c_tag_ = 0;
    }
    delta_c_ = delta_c;
    if (D_d) {
      d_d_tag_ = D_d->GetTag();
    }
    else {
      d_d_tag_ = 0;
    }
    delta_d_ = delta_d;
    j_c_tag_ = J_c.GetTag();
    j_d_tag_ = J_d.GetTag();
  }


  bool GenAugSystemSolver::AugmentedSystemChanged(
    const SymMatrix* W,
    double W_factor,
    const Vector* D_x,
    double delta_x,
    const Vector* D_s,
    double delta_s,
    const Matrix& J_c,
    const Vector* D_c,
    double delta_c,
    const Matrix& J_d,
    const Vector* D_d,
    double delta_d)
  {
    DBG_START_METH("GenAugSystemSolver::AugmentedSystemRequiresChange",dbg_verbosity);

#if COIN_IPOPT_VERBOSITY > 0

    bool Wtest = (W && W->GetTag() != w_tag_);
    bool iWtest = (!W && w_tag_ != 0);
    bool wfactor_test = (W_factor != w_factor_);
    bool D_xtest = (D_x && D_x->GetTag() != d_x_tag_);
    bool iD_xtest = (!D_x && d_x_tag_ != 0);
    bool delta_xtest = (delta_x != delta_x_);
    bool D_stest = (D_s && D_s->GetTag() != d_s_tag_);
    bool iD_stest = (!D_s && d_s_tag_ != 0);
    bool delta_stest = (delta_s != delta_s_);
    bool J_ctest = (J_c.GetTag() != j_c_tag_);
    bool D_ctest = (D_c && D_c->GetTag() != d_c_tag_);
    bool iD_ctest = (!D_c && d_c_tag_ != 0);
    bool delta_ctest = (delta_c != delta_c_);
    bool J_dtest = (J_d.GetTag() != j_d_tag_);
    bool D_dtest = (D_d && D_d->GetTag() != d_d_tag_);
    bool iD_dtest = (!D_d && d_d_tag_ != 0);
    bool delta_dtest = (delta_d != delta_d_);
#endif

    DBG_PRINT((2,"Wtest = %d\n", Wtest));
    DBG_PRINT((2,"iWtest = %d\n", iWtest));
    DBG_PRINT((2,"wfactor_test = %d\n", wfactor_test));
    DBG_PRINT((2,"D_xtest = %d\n", D_xtest));
    DBG_PRINT((2,"iD_xtest = %d\n", iD_xtest));
    DBG_PRINT((2,"delta_xtest = %d\n", delta_xtest));
    DBG_PRINT((2,"D_stest = %d\n", D_stest));
    DBG_PRINT((2,"iD_stest = %d\n", iD_stest));
    DBG_PRINT((2,"delta_stest = %d\n", delta_stest));
    DBG_PRINT((2,"J_ctest = %d\n", J_ctest));
    DBG_PRINT((2,"D_ctest = %d\n", D_ctest));
    DBG_PRINT((2,"iD_ctest = %d\n", iD_ctest));
    DBG_PRINT((2,"delta_ctest = %d\n", delta_ctest));
    DBG_PRINT((2,"J_dtest = %d\n", J_dtest));
    DBG_PRINT((2,"D_dtest = %d\n", D_dtest));
    DBG_PRINT((2,"iD_dtest = %d\n", iD_dtest));
    DBG_PRINT((2,"delta_dtest = %d\n", delta_dtest));

    if ( (W && W->GetTag() != w_tag_)
         || (!W && w_tag_ != 0)
         || (W_factor != w_factor_)
         || (D_x && D_x->GetTag() != d_x_tag_)
         || (!D_x && d_x_tag_ != 0)
         || (delta_x != delta_x_)
         || (D_s && D_s->GetTag() != d_s_tag_)
         || (!D_s && d_s_tag_ != 0)
         || (delta_s != delta_s_)
         || (J_c.GetTag() != j_c_tag_)
         || (D_c && D_c->GetTag() != d_c_tag_)
         || (!D_c && d_c_tag_ != 0)
         || (delta_c != delta_c_)
         || (J_d.GetTag() != j_d_tag_)
         || (D_d && D_d->GetTag() != d_d_tag_)
         || (!D_d && d_d_tag_ != 0)
         || (delta_d != delta_d_) ) {
      return true;
    }

    return false;
  }

  Index GenAugSystemSolver::NumberOfNegEVals() const
  {
    return solver_interface_->NumberOfNegEVals();
  }

  bool GenAugSystemSolver::ProvidesInertia() const
  {
    return solver_interface_->ProvidesInertia();
  }

  bool GenAugSystemSolver::IncreaseQuality()
  {
    return solver_interface_->IncreaseQuality();
  }

} // namespace Ipopt
