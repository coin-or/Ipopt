// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpOrigIpoptNLP.cpp 321 2005-06-20 21:53:55Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpGradientScaling.hpp"
#include "IpTripletHelper.hpp"

namespace Ipopt
{

  DefineIpoptType(GradientScaling);

  void GradientScaling::RegisterOptions(SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->AddLowerBoundedNumberOption("scaling_max_gradient", "maximum gradient after scaling",
                                          1, true, 1000.0,
                                          "This is the gradient scaling cut-off. If the maximum"
                                          " gradient is above this value, then gradient based scaling"
                                          " will be performed. Scaling parameters will scale the maximum"
                                          " gradient back to this value. Note: This option is only used"
                                          " if scaling_method = gradient_based");
  }

  bool GradientScaling::Initialize(const Journalist& jnlst,
                                   const OptionsList& options,
                                   const std::string& prefix)
  {
    options.GetNumericValue("scaling_max_gradient", scaling_max_gradient_, prefix);
    return StandardScalingBase::Initialize(jnlst, options, prefix);
  }


  void GradientScaling::DetermineScalingParametersImpl(
    const SmartPtr<const VectorSpace> x_space,
    const SmartPtr<const VectorSpace> c_space,
    const SmartPtr<const VectorSpace> d_space,
    const SmartPtr<const MatrixSpace> jac_c_space,
    const SmartPtr<const MatrixSpace> jac_d_space,
    const SmartPtr<const SymMatrixSpace> h_space,
    Number& df, Vector& dx,
    Vector& dc, Vector& dd)
  {
    DBG_ASSERT(IsValid(nlp_));

    SmartPtr<Vector> x = x_space->MakeNew();
    nlp_->GetStartingPoint(GetRawPtr(x), true,
                           NULL, false,
                           NULL, false,
                           NULL, false,
                           NULL, false,
                           NULL, false,
                           NULL, false);

    //
    // Calculate grad_f scaling
    //
    SmartPtr<Vector> grad_f = x_space->MakeNew();
    nlp_->Eval_grad_f(*x, *grad_f);
    double max_grad_f = grad_f->Amax();
    df = 1.0;
    if (max_grad_f > scaling_max_gradient_) {
      df = scaling_max_gradient_ / max_grad_f;
    }

    //
    // calculate x scaling
    //
    dx.Set(1.0);

    //
    // Calculate c scaling
    //
    SmartPtr<Matrix> jac_c = jac_c_space->MakeNew();
    nlp_->Eval_jac_c(*x, *jac_c);
    Index nnz = TripletHelper::GetNumberEntries(*jac_c);
    Index* irow = new Index[nnz];
    Index* jcol = new Index[nnz];
    Number* values = new Number[nnz];
    TripletHelper::FillRowCol(nnz, *jac_c, irow, jcol);
    TripletHelper::FillValues(nnz, *jac_c, values);
    Number* c_scaling = new Number[jac_c->NRows()];

    for (Index r=0; r<jac_c->NRows(); r++) {
      c_scaling[r] = 0;
    }

    // put the max of each row into c_scaling...
    for (Index i=0; i<nnz; i++) {
      if (values[i] > scaling_max_gradient_
          && values[i] > c_scaling[irow[i]-1]) {
        c_scaling[irow[i]-1] = values[i];
      }
    }

    // now compute the scaling factors for each row
    for (Index r=0; r<jac_c->NRows(); r++) {
      Number scaling = 1.0;
      if (c_scaling[r] > scaling_max_gradient_) {
        scaling = scaling_max_gradient_/c_scaling[r];
      }
      c_scaling[r] = scaling;
    }

    TripletHelper::PutValuesInVector(jac_c->NRows(), c_scaling, dc);
    delete [] irow;
    irow = NULL;
    delete [] jcol;
    jcol = NULL;
    delete [] values;
    values = NULL;
    delete [] c_scaling;
    c_scaling = NULL;

    //
    // Calculate d scaling
    //
    SmartPtr<Matrix> jac_d = jac_d_space->MakeNew();
    nlp_->Eval_jac_d(*x, *jac_d);
    nnz = TripletHelper::GetNumberEntries(*jac_d);
    irow = new Index[nnz];
    jcol = new Index[nnz];
    values = new Number[nnz];
    TripletHelper::FillRowCol(nnz, *jac_d, irow, jcol);
    TripletHelper::FillValues(nnz, *jac_d, values);
    Number* d_scaling = new Number[jac_d->NRows()];

    for (Index r=0; r<jac_d->NRows(); r++) {
      d_scaling[r] = 0;
    }

    // put the max of each row into c_scaling...
    for (Index i=0; i<nnz; i++) {
      if (values[i] > scaling_max_gradient_
          && values[i] > d_scaling[irow[i]-1]) {
        d_scaling[irow[i]-1] = values[i];
      }
    }

    // now compute the scaling factors for each row
    for (Index r=0; r<jac_d->NRows(); r++) {
      Number scaling = 1.0;
      if (d_scaling[r] > scaling_max_gradient_) {
        scaling = scaling_max_gradient_/d_scaling[r];
      }
      d_scaling[r] = scaling;
    }

    TripletHelper::PutValuesInVector(jac_d->NRows(), d_scaling, dd);
    delete [] irow;
    irow = NULL;
    delete [] jcol;
    jcol = NULL;
    delete [] values;
    values = NULL;
    delete [] d_scaling;
    d_scaling = NULL;
  }

} // namespace Ipopt
