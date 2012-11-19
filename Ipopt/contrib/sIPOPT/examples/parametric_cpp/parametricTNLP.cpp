// Copyright 2010, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2010-30-04


#include "parametricTNLP.hpp"
#include "IpDenseVector.hpp"
#include "IpIpoptData.hpp"
#include <stdio.h>

using namespace Ipopt;


/* Constructor */
ParametricTNLP::ParametricTNLP() :
  nominal_eta1_(5.0),
  nominal_eta2_(1.0),
  eta_1_perturbed_value_(4.5),
  eta_2_perturbed_value_(1.0)
{
}

ParametricTNLP::~ParametricTNLP()
{
}

bool ParametricTNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
				  Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // x1, x2, x3, eta1, eta2
  n = 5;

  // 2 constraints + 2 parametric initial value constraints
  m = 4;

  nnz_jac_g = 10;

  nnz_h_lag = 5;

  index_style = FORTRAN_STYLE;

  return true;
}

bool ParametricTNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
				     Index m, Number* g_l, Number* g_u)
{
  for (Index k=0; k<3; ++k) {
    x_l[k] = 0.0;
    x_u[k] = 1.0e19;
  }
  x_l[3] = -1.0e19;
  x_u[3] = 1.0e19;
  x_l[4] = -1.0e19;
  x_u[4] = 1.0e19;

  g_l[0] = 0.0;
  g_u[0] = 0.0;
  g_l[1] = 0.0;
  g_u[1] = 0.0;

  // initial value constraints
  g_l[2] = nominal_eta1_;
  g_u[2] = nominal_eta1_;
  g_l[3] = nominal_eta2_;
  g_u[3] = nominal_eta2_;

  return true;
}

bool ParametricTNLP::get_starting_point(Index n, bool init_x, Number* x,
					bool init_z, Number* z_L, Number* z_U,
					Index m, bool init_lambda,
					Number* lambda)
{
  x[0] = 0.15;
  x[1] = 0.15;
  x[2] = 0.0;
  x[3] = 0.0;
  x[4] = 0.0;

  return true;
}

bool ParametricTNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  obj_value = 0;
  for (Index k=0; k<3; ++k) {
    obj_value += x[k]*x[k];
  }
  return true;
}

bool ParametricTNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  grad_f[0] = 2*x[0];
  grad_f[1] = 2*x[1];
  grad_f[2] = 2*x[2];
  grad_f[3] = 0.0;
  grad_f[4] = 0.0;
  return true;
}

bool ParametricTNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  Number x1, x2, x3, eta1, eta2;
  x1 = x[0];
  x2 = x[1];
  x3 = x[2];
  eta1 = x[3];
  eta2 = x[4];
  g[0] = 6*x1+3*x2+2*x3-eta1;
  g[1] = eta2*x1+x2-x3-1;
  g[2] = eta1;
  g[3] = eta2;
  return true;
}

bool ParametricTNLP::eval_jac_g(Index n, const Number* x, bool new_x,
				Index m, Index nele_jac, Index* iRow, Index *jCol,
				Number* values)
{
  if (values == NULL) {
    iRow[0] = 1; // dg1/dx1
    jCol[0] = 1;
    iRow[1] = 1; // dg1/dx2
    jCol[1] = 2;
    iRow[2] = 1; // dg1/dx3
    jCol[2] = 3;
    iRow[3] = 1; // dg1/deta1
    jCol[3] = 4;
    iRow[4] = 2; // dg2/dx1
    jCol[4] = 1;
    iRow[5] = 2; // dg2/dx2
    jCol[5] = 2;
    iRow[6] = 2; // dg2/dx3
    jCol[6] = 3;
    iRow[7] = 2; // dg2/deta2
    jCol[7] = 5;
    iRow[8] = 3;
    jCol[8] = 4;
    iRow[9] = 4;
    jCol[9] = 5;
  }
  else {
    values[0] = 6.0;
    values[1] = 3.0;
    values[2] = 2.0;
    values[3] = -1.0;
    values[4] = x[4];
    values[5] = 1.0;
    values[6] = -1.0;
    values[7] = x[0];
    values[8] = 1.0;
    values[9] = 1.0;
  }
  return true;
}

bool ParametricTNLP::eval_h(Index n, const Number* x, bool new_x,
			    Number obj_factor, Index m, const Number* lambda,
			    bool new_lambda, Index nele_hess, Index* iRow,
			    Index* jCol, Number* values)
{
  if (values == NULL) {
    iRow[0] = 1;
    jCol[0] = 1;

    iRow[1] = 2;
    jCol[1] = 2;

    iRow[2] = 3;
    jCol[2] = 3;

    iRow[3] = 1;
    jCol[3] = 5;

    iRow[4] = 5;
    jCol[4] = 1;
  }
  else {
    values[0] = 2.0*obj_factor;
    values[1] = 2.0*obj_factor;
    values[2] = 2.0*obj_factor;
    values[3] = 0.5*lambda[1];
    values[4] = 0.5*lambda[1];
  }
  return true;
}

bool ParametricTNLP::get_var_con_metadata(Index n,
					  StringMetaDataMapType& var_string_md,
					  IntegerMetaDataMapType& var_integer_md,
					  NumericMetaDataMapType& var_numeric_md,
					  Index m,
					  StringMetaDataMapType& con_string_md,
					  IntegerMetaDataMapType& con_integer_md,
					  NumericMetaDataMapType& con_numeric_md)
{
  /* In this function, the indices for the parametric computations are set.
   * To keep track of the parameters, each parameter gets an index from 1 to n_parameters.
   * In this case, [1] eta_1, [2] eta_2.
   * The following metadata vectors are important:
   */


  /* 1. sens_init_constr: in this list, the constraints that set the initial
   *    values for the parameters are indicated.
   *    For parameter 1 (eta_1) this is constraint 3 (e.g. C++ index 2), which is
   *    the constraint   eta_1 = eta_1_nominal;
   *    For parameter 2 (eta_2) this is constraint 4 (e.g. C++ index 3).
   */
  std::vector<Index> sens_init_constr(m,0);
  sens_init_constr[2] = 1;
  sens_init_constr[3] = 2;
  con_integer_md["sens_init_constr"] = sens_init_constr;

  /* 2. sens_state_1: in this index list, the parameters are indicated:
   *    Here: [1] eta_1, [2] eta_2
   */
  std::vector<Index> sens_state_1(n,0);
  sens_state_1[3] = 1;
  sens_state_1[4] = 2;
  var_integer_md["sens_state_1"] = sens_state_1;

  /* 3. sens_state_values_1: In this list of Numbers (=doubles), the perturbed
   *    values for the parameters are set.
   */
  std::vector<Number> sens_state_value_1(n,0);
  sens_state_value_1[3] = eta_1_perturbed_value_;
  sens_state_value_1[4] = eta_2_perturbed_value_;
  var_numeric_md["sens_state_value_1"] = sens_state_value_1;

  return true;
}

void ParametricTNLP::finalize_solution(SolverReturn status,
				       Index n, const Number* x, const Number* z_L, const Number* z_U,
				       Index m, const Number* g, const Number* lambda,
				       Number obj_value,
				       const IpoptData* ip_data,
				       IpoptCalculatedQuantities* ip_cq)
{
  // Check whether sIPOPT Algorithm aborted internally
  //  bool sens_internal_abort;
  //options_->GetBoolValue("sens_internal_abort", sens_internal_abort, "");

  // Get access to the metadata, where the solutions are stored. The metadata is part of the DenseVectorSpace.
  SmartPtr<const DenseVectorSpace> x_owner_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(ip_data->curr()->x()->OwnerSpace()));

  if (!IsValid(x_owner_space)) {
    printf("Error IsValid(x_owner_space) failed\n");
    return;
  }
  std::string state;
  std::vector<Number> sens_sol_vec;
  state = "sens_sol_state_1";
  sens_sol_vec = x_owner_space->GetNumericMetaData(state.c_str());

  // Print the solution vector
  printf("\n"
	 "                Nominal                    Perturbed\n");
  for (Index k=0; k<(Index) sens_sol_vec.size(); ++k) {
    printf("x[%3d]   % .23f   % .23f\n", k, x[k], sens_sol_vec[k]);
  }
}

void ParametricTNLP::finalize_metadata(Index n,
				       const StringMetaDataMapType& var_string_md,
				       const IntegerMetaDataMapType& var_integer_md,
				       const NumericMetaDataMapType& var_numeric_md,
				       Index m,
				       const StringMetaDataMapType& con_string_md,
				       const IntegerMetaDataMapType& con_integer_md,
				       const NumericMetaDataMapType& con_numeric_md)
{
  // bound multipliers for lower and upper bounds
  printf("\nDual bound multipliers:\n");
  NumericMetaDataMapType::const_iterator z_L_solution = var_numeric_md.find("sens_sol_state_1_z_L");
  NumericMetaDataMapType::const_iterator z_U_solution = var_numeric_md.find("sens_sol_state_1_z_U");
  if (z_L_solution!=var_numeric_md.end() && z_U_solution!=var_numeric_md.end()) {
    for (Index k=0; k<n; ++k) {
      printf("z_L[%d] = %f      z_U[%d] = %f\n", k, z_L_solution->second[k], k, z_U_solution->second[k]);
    }
  }

  // constraint mutlipliers
  printf("\nConstraint multipliers:\n");
  NumericMetaDataMapType::const_iterator lambda_solution = con_numeric_md.find("sens_sol_state_1");
  if (lambda_solution!=con_numeric_md.end()) {
    for (Index k=0; k<m; ++k) {
      printf("lambda[%d] = %f\n", k, lambda_solution->second[k]);
    }
  }
}
