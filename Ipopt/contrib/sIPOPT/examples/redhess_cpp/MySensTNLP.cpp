// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2010-10-04


#include "MySensTNLP.hpp"
#include "IpDenseVector.hpp"

#include <cassert>

using namespace Ipopt;

/* Constructor. */
MySensTNLP::MySensTNLP()
{}

MySensTNLP::~MySensTNLP()
{}

bool MySensTNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
			     Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem has 3 variables
  n = 3;

  // one equality constraint
  m = 1;

  nnz_jac_g = 3;

  nnz_h_lag = 3;

  index_style = FORTRAN_STYLE;

  return true;
}

bool MySensTNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
				 Index m, Number* g_l, Number* g_u)
{
  assert(n==3);
  assert(m==1);

  for (Index k=0; k<3; k++) {
    x_l[k] = -1.0e19;
    x_u[k] = +1.0e19;
  }

  g_l[0] = 0.0;
  g_u[0] = 0.0;

  return true;
}

bool MySensTNLP::get_starting_point(Index n, bool init_x, Number* x,
				   bool init_z, Number* z_L, Number* z_U,
				   Index m, bool init_lambda,
				   Number* lambda)
{

  x[0] = 25;
  x[1] = 0;
  x[2] = 0;

  return true;
}

bool MySensTNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  // return the value of the objective function
  Number x1 = x[0];
  Number x2 = x[1];
  Number x3 = x[2];
  obj_value = (x1-1)*(x1-1) + (x2-2)*(x2-2) + (x3-3)*(x3-3);

  return true;
}

bool MySensTNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // return the gradient of the objective function grad_{x} f(x)

  Number x1 = x[0];
  Number x2 = x[1];
  Number x3 = x[2];

  grad_f[0] = 2*(x1-1);
  grad_f[1] = 2*(x2-2);
  grad_f[2] = 2*(x3-3);

  return true;
}

bool MySensTNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  // return the value of the constraints: g(x)
  Number x1 = x[0];
  Number x2 = x[1];
  Number x3 = x[2];

  g[0] = x1+2*x2+3*x3;

  return true;
}

bool MySensTNLP::eval_jac_g(Index n, const Number* x, bool new_x,
			   Index m, Index nele_jac, Index* iRow, Index *jCol,
			   Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian of the constraints

    // element at 1,1: grad_{x1} g_{1}(x)
    iRow[0] = 1;
    jCol[0] = 1;

    // element at 1,2: grad_{x2} g_{1}(x)
    iRow[1] = 1;
    jCol[1] = 2;

    // element at 1,3: grad_{x3} g_{1}(x)
    iRow[2] = 1;
    jCol[2] = 3;
  }
  else {
    // return the values of the jacobian of the constraints

    // element at 1,1: grad_{x1} g_{1}(x)
    values[0] = 1.0;

    // element at 1,2: grad_{x1} g_{1}(x)
    values[1] = 2.0;

    values[2] = 3.0;
  }

  return true;
}

bool MySensTNLP::eval_h(Index n, const Number* x, bool new_x,
		       Number obj_factor, Index m, const Number* lambda,
		       bool new_lambda, Index nele_hess, Index* iRow,
		       Index* jCol, Number* values)
{
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // element at 1,1: grad^2_{x1,x1} L(x,lambda)
    iRow[0] = 1;
    jCol[0] = 1;

    // element at 2,2: grad^2_{x2,x2} L(x,lambda)
    iRow[1] = 2;
    jCol[1] = 2;

    iRow[2] = 3;
    jCol[2] = 3;

    // Note: off-diagonal elements are zero for this problem
  }
  else {
    // return the values

    // element at 1,1: grad^2_{x1,x1} L(x,lambda)
    values[0] = 2.0;

    // element at 2,2: grad^2_{x2,x2} L(x,lambda)
    values[1] = 2.0;

    values[2] = 2.0;
    // Note: off-diagonal elements are zero for this problem
  }

  return true;
}


bool MySensTNLP::get_var_con_metadata(Index n,
				      StringMetaDataMapType& var_string_md,
				      IntegerMetaDataMapType& var_integer_md,
                                      NumericMetaDataMapType& var_numeric_md,
                                      Index m,
                                      StringMetaDataMapType& con_string_md,
                                      IntegerMetaDataMapType& con_integer_md,
                                      NumericMetaDataMapType& con_numeric_md)
{
  std::vector<Index> red_hess_idx(3,0);
  red_hess_idx[1] = 1;
  red_hess_idx[2] = 2;
  
  var_integer_md["red_hessian"] = red_hess_idx;

  return true;
}


void MySensTNLP::finalize_solution(SolverReturn status,
				  Index n, const Number* x, const Number* z_L, const Number* z_U,
				  Index m, const Number* g, const Number* lambda,
				  Number obj_value,
				  const IpoptData* ip_data,
				  IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution. Since the solution is displayed to the console,
  // we currently do nothing here.
}
