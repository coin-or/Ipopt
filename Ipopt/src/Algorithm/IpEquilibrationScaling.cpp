// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2007-05-21

#include "IpoptConfig.h"
#include "IpEquilibrationScaling.hpp"
#include "IpTripletHelper.hpp"

#ifdef COIN_HAS_HSL
#include "CoinHslConfig.h"
#endif

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

/** Prototypes for MA27's Fortran subroutines */
extern "C"
{
  // here we assume that float corresponds to Fortran's single
  // precision
  void F77_FUNC(mc19ad,MC19AD)(const ipfint *N, const ipfint *NZ,
                               const double* A, const ipfint *IRN,
                               const ipfint* ICN, float* R, float* C,
                               float* W);
}

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  void EquilibrationScaling::
  RegisterOptions(const SmartPtr<RegisteredOptions>& roptions)
  {}

  bool EquilibrationScaling::
  InitializeImpl(const OptionsList& options,
                 const std::string& prefix)
  {
    options.GetNumericValue("point_perturbation_radius",
                            point_perturbation_radius_, prefix);
    return StandardScalingBase::InitializeImpl(options, prefix);
  }


  void EquilibrationScaling::DetermineScalingParametersImpl(
    const SmartPtr<const VectorSpace> x_space,
    const SmartPtr<const VectorSpace> c_space,
    const SmartPtr<const VectorSpace> d_space,
    const SmartPtr<const MatrixSpace> jac_c_space,
    const SmartPtr<const MatrixSpace> jac_d_space,
    const SmartPtr<const SymMatrixSpace> h_space,
    const Matrix& Px_L, const Vector& x_L,
    const Matrix& Px_U, const Vector& x_U,
    Number& df,
    SmartPtr<Vector>& dx,
    SmartPtr<Vector>& dc,
    SmartPtr<Vector>& dd)
  {
    DBG_ASSERT(IsValid(nlp_));

    SmartPtr<Vector> x0 = x_space->MakeNew();
    if (!nlp_->GetStartingPoint(GetRawPtr(x0), true,
                                NULL, false,
                                NULL, false,
                                NULL, false,
                                NULL, false)) {
      THROW_EXCEPTION(FAILED_INITIALIZATION,
                      "Error getting initial point from NLP in EquilibrationScaling.\n");
    }

    // We store the added absolute values of the Jacobian and
    // objective function gradient in an array of sufficient size

    SmartPtr<Matrix> jac_c = jac_c_space->MakeNew();
    SmartPtr<Matrix> jac_d = jac_d_space->MakeNew();
    SmartPtr<Vector> grad_f = x_space->MakeNew();
    const Index nnz_jac_c = TripletHelper::GetNumberEntries(*jac_c);
    const Index nnz_jac_d = TripletHelper::GetNumberEntries(*jac_d);
    const Index nc = jac_c_space->NRows();
    const Index nd = jac_d_space->NRows();
    const Index nx = x_space->Dim();
    Number* avrg_values = new Number[nnz_jac_c+nnz_jac_d+nx];
    Number* val_buffer = new Number[Max(nnz_jac_c,nnz_jac_d,nx)];

    SmartPtr<PointPerturber> perturber =
      new PointPerturber(*x0, point_perturbation_radius_,
                         Px_L, x_L, Px_U, x_U);

    const Index num_evals = 4;
    const Index max_num_eval_errors = 10;
    Index num_eval_errors = 0;
    for (Index ieval=0; ieval<num_evals; ieval++) {
      // Compute obj gradient and Jacobian at random perturbation point
      bool done = false;
      while (!done) {
        SmartPtr<Vector> xpert = perturber->MakeNewPerturbedPoint();
        done = (nlp_->Eval_grad_f(*xpert, *grad_f) &&
                nlp_->Eval_jac_c(*xpert, *jac_c) &&
                nlp_->Eval_jac_d(*xpert, *jac_d));
        if (!done) {
          Jnlst().Printf(J_WARNING, J_INITIALIZATION,
                         "Error evaluating first derivatives as at perturbed point for equilibration-based scaling.\n");
          num_eval_errors++;
        }
        if (num_eval_errors>max_num_eval_errors) {
          delete [] val_buffer;
          delete [] avrg_values;
          THROW_EXCEPTION(FAILED_INITIALIZATION,
                          "Too many evaluation failures during equilibiration-based scaling.");
        }
      }
      // Get the numbers out of the matrices and vectors, and add it
      // to avrg_values
      TripletHelper::FillValues(nnz_jac_c, *jac_c, val_buffer);
      if (ieval==0) {
        for (Index i=0; i<nnz_jac_c; i++) {
          avrg_values[i] = fabs(val_buffer[i]);
        }
      }
      else {
        for (Index i=0; i<nnz_jac_c; i++) {
          avrg_values[i] += fabs(val_buffer[i]);
        }
      }
      TripletHelper::FillValues(nnz_jac_d, *jac_d, val_buffer);
      if (ieval==0) {
        for (Index i=0; i<nnz_jac_d; i++) {
          avrg_values[nnz_jac_c+i] = fabs(val_buffer[i]);
        }
      }
      else {
        for (Index i=0; i<nnz_jac_d; i++) {
          avrg_values[nnz_jac_c+i] += fabs(val_buffer[i]);
        }
      }
      TripletHelper::FillValuesFromVector(nx, *grad_f, val_buffer);
      if (ieval==0) {
        for (Index i=0; i<nx; i++) {
          avrg_values[nnz_jac_c+nnz_jac_d+i] = fabs(val_buffer[i]);
        }
      }
      else {
        for (Index i=0; i<nx; i++) {
          avrg_values[nnz_jac_c+nnz_jac_d+i] += fabs(val_buffer[i]);
        }
      }
    }
    delete [] val_buffer;
    for (Index i=0; i<nnz_jac_c+nnz_jac_d+nx; i++) {
      avrg_values[i] /= (Number)num_evals;
    }

    // Get the sparsity structure
    ipfint* AIRN = new ipfint[nnz_jac_c+nnz_jac_d+nx];
    ipfint* AJCN = new ipfint[nnz_jac_c+nnz_jac_d+nx];
    if (sizeof(ipfint)==sizeof(Index)) {
      TripletHelper::FillRowCol(nnz_jac_c, *jac_c, &AIRN[0], &AJCN[0]);
      TripletHelper::FillRowCol(nnz_jac_d, *jac_d, &AIRN[nnz_jac_c], &AJCN[nnz_jac_c], nc);
    }
    else {
      THROW_EXCEPTION(INTERNAL_ABORT,
                      "Need to implement missing code in EquilibriationScaling.");
    }
    // sort out the zero entries in objective function gradient
    Index nnz_grad_f = 0;
    const Index idx = nnz_jac_c+nnz_jac_d;
    for (Index i=0; i<nx; i++) {
      if (avrg_values[idx+i] != 0.) {
        AIRN[idx+nnz_grad_f] = nc+nd+1;
        AJCN[idx+nnz_grad_f] = i+1;
        avrg_values[idx+nnz_grad_f] = avrg_values[idx+i];
        nnz_grad_f++;
      }
    }

    // Now call MC19 to compute the scaling factors
    const ipfint N = Max(nc+nd+1,nx);
    float* R = new float[N];
    float* C = new float[N];
    float* W = new float[5*N];
#if defined(COINHSL_HAS_MC19) || defined(HAVE_LINEARSOLVERLOADER)
    const ipfint NZ = nnz_jac_c+nnz_jac_d+nnz_grad_f;
    //F77_FUNC(mc19ad,MC19AD)(&N, &NZ, avrg_values, AIRN, AJCN, R, C, W);
    F77_FUNC(mc19ad,MC19AD)(&N, &NZ, avrg_values, AJCN, AIRN, C, R, W);
#else

    THROW_EXCEPTION(OPTION_INVALID,
                    "Currently cannot do equilibration-based NLP scaling if MC19 is not available.");
#endif

    delete[] W;

    delete [] avrg_values;
    delete [] AIRN;
    delete [] AJCN;

    // Correct the scaling values
    Number* row_scale = new Number[nc+nd+1];
    Number* col_scale = new Number[nx];
    for (Index i=0; i<nc+nd+1; i++) {
      row_scale[i] = exp((Number)R[i]);
    }
    for (Index i=0; i<nx; i++) {
      col_scale[i] = exp((Number)C[i]);
    }
    delete [] R;
    delete [] C;

    // get the scaling factors
    df = row_scale[nc+nd];
    dc = c_space->MakeNew();
    TripletHelper::PutValuesInVector(nc, &row_scale[0], *dc);
    dd = d_space->MakeNew();
    TripletHelper::PutValuesInVector(nd, &row_scale[nc], *dd);
    dx = x_space->MakeNew();
    TripletHelper::PutValuesInVector(nx, col_scale, *dx);

    delete [] row_scale;
    delete [] col_scale;
  }

  PointPerturber::PointPerturber(const Vector& x0,
                                 Number random_pert_radius,
                                 const Matrix& Px_L, const Vector& x_L,
                                 const Matrix& Px_U, const Vector& x_U)
  {
    DBG_START_METH("PointPerturber::PointPerturber", dbg_verbosity);
    const Number very_large = 1e300;
    // First we compute full-space lower and upper bounds
    SmartPtr<Vector> full_x_L = x0.MakeNew();
    full_x_L->Set(-very_large);
    SmartPtr<Vector> tmp = x_L.MakeNew();
    tmp->Set(very_large);
    Px_L.MultVector(1., *tmp, 1., *full_x_L);
    DBG_PRINT_VECTOR(1,"full_x_L1", *full_x_L);
    Px_L.MultVector(1., x_L, 1., *full_x_L);
    DBG_PRINT_VECTOR(1,"full_x_L2", *full_x_L);

    SmartPtr<Vector> full_x_U = x0.MakeNew();
    full_x_U->Set(very_large);
    tmp = x_U.MakeNew();
    tmp->Set(-very_large);
    Px_U.MultVector(1., *tmp, 1., *full_x_U);
    DBG_PRINT_VECTOR(1,"full_x_U1", *full_x_U);
    Px_U.MultVector(1., x_U, 1., *full_x_U);
    DBG_PRINT_VECTOR(1,"full_x_U2", *full_x_U);

    pert_dir_ = full_x_U->MakeNew();
    pert_dir_->AddTwoVectors(.5, *full_x_U, -.5, *full_x_L, 0.);
    tmp = full_x_U->MakeNew();
    tmp->Set(random_pert_radius);
    pert_dir_->ElementWiseMin(*tmp);
    DBG_PRINT_VECTOR(1,"pert_dir", *pert_dir_);
    ref_point_ = x0.MakeNewCopy();
    DBG_PRINT_VECTOR(1,"ref_point1", *ref_point_);
    full_x_U->AddOneVector(-1., *pert_dir_, 1.);
    ref_point_->ElementWiseMin(*full_x_U);
    DBG_PRINT_VECTOR(1,"ref_point2", *ref_point_);
    full_x_L->AddOneVector(1., *pert_dir_, 1.);
    ref_point_->ElementWiseMax(*full_x_L);
    DBG_PRINT_VECTOR(1,"ref_point3", *ref_point_);
  }

  SmartPtr<Vector> PointPerturber::
  MakeNewPerturbedPoint() const
  {
    const Index nx = ref_point_->Dim();
    Number* vals = new Number[nx];
    TripletHelper::FillValuesFromVector(nx, *ref_point_, vals);
    Number* pert_vals = new Number[nx];
    TripletHelper::FillValuesFromVector(nx, *pert_dir_, pert_vals);

    for (Index i=0; i<nx; i++) {
      Number random = IpRandom01();
      vals[i] += 2.*(random-0.5)*pert_vals[i];
    }
    delete [] pert_vals;

    SmartPtr<Vector> retval = ref_point_->MakeNew();
    TripletHelper::PutValuesInVector(nx, vals, *retval);

    delete [] vals;

    return retval;
  }

} // namespace Ipopt
