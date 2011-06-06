// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "MittelmannDistCntrlDiri.hpp"

#ifdef HAVE_CASSERT
# include <cassert>
#else
# ifdef HAVE_ASSERT_H
#  include <assert.h>
# else
#  error "don't have header file for assert"
# endif
#endif

using namespace Ipopt;

/* Constructor. */
MittelmannDistCntrlDiriBase::MittelmannDistCntrlDiriBase()
    :
    y_d_(NULL)
{}

MittelmannDistCntrlDiriBase::~MittelmannDistCntrlDiriBase()
{
  delete [] y_d_;
}

void
MittelmannDistCntrlDiriBase::SetBaseParameters(Index N, Number alpha, Number lb_y,
    Number ub_y, Number lb_u, Number ub_u,
    Number u_init)
{
  N_ = N;
  h_ = 1./(N+1);
  hh_ = h_*h_;
  lb_y_ = lb_y;
  ub_y_ = ub_y;
  lb_u_ = lb_u;
  ub_u_ = ub_u;
  u_init_ = u_init;
  alpha_ = alpha;

  // Initialize the target profile variables
  delete [] y_d_;
  y_d_ = new Number[(N_+2)*(N_+2)];
  for (Index j=0; j<= N_+1; j++) {
    for (Index i=0; i<= N_+1; i++) {
      y_d_[y_index(i,j)] = y_d_cont(x1_grid(i),x2_grid(j));
    }
  }
}

bool MittelmannDistCntrlDiriBase::get_nlp_info(
  Index& n, Index& m, Index& nnz_jac_g,
  Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // We for each of the N_+2 times N_+2 mesh points we have the value
  // of the functions y, and for each N_ tiems N_ interior mesh points
  // we have values for u
  n = (N_+2)*(N_+2) + N_*N_;

  // For each of the N_ times N_ interior mesh points we have the
  // discretized PDE.
  m = N_*N_;

  // y(i,j), y(i-1,j), y(i+1,j), y(i,j-1), y(i,j+1), u(i,j) for each
  // of the N_*N_ discretized PDEs
  nnz_jac_g = 6*N_*N_;

  // diagonal entry for each y(i,j) in the interior
  nnz_h_lag = N_*N_;
  if (alpha_>0.) {
    // and one entry for u(i,j) in the interior if alpha is not zero
    nnz_h_lag += N_*N_;
  }

  // We use the C indexing style for row/col entries (corresponding to
  // the C notation, starting at 0)
  index_style = C_STYLE;

  return true;
}

bool
MittelmannDistCntrlDiriBase::get_bounds_info(Index n, Number* x_l, Number* x_u,
    Index m, Number* g_l, Number* g_u)
{
  // Set overall bounds on the variables
  for (Index i=0; i<=N_+1; i++) {
    for (Index j=0; j<=N_+1; j++) {
      Index iy = y_index(i,j);
      x_l[iy] = lb_y_;
      x_u[iy] = ub_y_;
    }
  }
  for (Index i=1; i<=N_; i++) {
    for (Index j=1; j<=N_; j++) {
      Index iu = u_index(i,j);
      x_l[iu] = lb_u_;
      x_u[iu] = ub_u_;
    }
  }

  // Define the boundary condition on y as bounds
  for (Index i=0; i<=N_+1; i++) {
    x_l[y_index(i,0)] = 0.;
    x_u[y_index(i,0)] = 0.;
  }
  for (Index i=0; i<=N_+1; i++) {
    x_l[y_index(0,i)] = 0.;
    x_u[y_index(0,i)] = 0.;
  }
  for (Index i=0; i<=N_+1; i++) {
    x_l[y_index(i,N_+1)] = 0.;
    x_u[y_index(i,N_+1)] = 0.;
  }
  for (Index i=0; i<=N_+1; i++) {
    x_l[y_index(N_+1,i)] = 0.;
    x_u[y_index(N_+1,i)] = 0.;
  }

  // all discretized PDE constraints have right hand side zero
  for (Index i=0; i<m; i++) {
    g_l[i] = 0.;
    g_u[i] = 0.;
  }

  return true;
}

bool
MittelmannDistCntrlDiriBase::get_starting_point(Index n, bool init_x, Number* x,
    bool init_z, Number* z_L, Number* z_U,
    Index m, bool init_lambda,
    Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // set all y's to the perfect match with y_d
  for (Index i=0; i<= N_+1; i++) {
    for (Index j=0; j<= N_+1; j++) {
      x[y_index(i,j)] = y_d_[y_index(i,j)];
      //      x[y_index(i,j)] += h_*x1_grid(i) + 2*h_*x2_grid(j);
    }
  }

  // Set the initial (constant) value for the u's
  for (Index i=1; i<= N_; i++) {
    for (Index j=1; j<= N_; j++) {
      x[u_index(i,j)] = u_init_;
      //      x[u_index(i,j)] -= h_*x1_grid(i) + 2*h_*x2_grid(j);
    }
  }

  return true;
}

bool
MittelmannDistCntrlDiriBase::get_scaling_parameters(Number& obj_scaling,
    bool& use_x_scaling, Index n, Number* x_scaling,
    bool& use_g_scaling, Index m, Number* g_scaling)
{
  obj_scaling = 1./hh_;
  use_x_scaling = false;
  use_g_scaling = false;
  return true;
}

bool
MittelmannDistCntrlDiriBase::eval_f(Index n, const Number* x,
                                    bool new_x, Number& obj_value)
{
  // return the value of the objective function
  obj_value = 0.;
  for (Index i=1; i<=N_; i++) {
    for (Index j=1; j<= N_; j++) {
      Index iy = y_index(i,j);
      Number tmp = x[iy] - y_d_[iy];
      obj_value += tmp*tmp;
    }
  }
  obj_value *= hh_/2.;

  if (alpha_>0.) {
    Number usum = 0.;
    for (Index i=1; i<=N_; i++) {
      for (Index j=1; j<= N_; j++) {
        Index iu = u_index(i,j);
        usum += x[iu]*x[iu];
      }
    }
    obj_value += alpha_*hh_/2.*usum;
  }

  return true;
}

bool
MittelmannDistCntrlDiriBase::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // return the gradient of the objective function grad_{x} f(x)

  // The values are zero for variables on the boundary
  for (Index i=0; i<= N_+1; i++) {
    grad_f[y_index(i,0)] = 0.;
  }
  for (Index i=0; i<= N_+1; i++) {
    grad_f[y_index(i,N_+1)] = 0.;
  }
  for (Index j=1; j<= N_; j++) {
    grad_f[y_index(0,j)] = 0.;
  }
  for (Index j=1; j<= N_; j++) {
    grad_f[y_index(N_+1,j)] = 0.;
  }

  // now let's take care of the nonzero values
  for (Index i=1; i<=N_; i++) {
    for (Index j=1; j<= N_; j++) {
      Index iy = y_index(i,j);
      grad_f[iy] = hh_*(x[iy] - y_d_[iy]);
    }
  }

  if (alpha_>0.) {
    for (Index i=1; i<=N_; i++) {
      for (Index j=1; j<= N_; j++) {
        Index iu = u_index(i,j);
        grad_f[iu] = alpha_*hh_*x[iu];
      }
    }
  }
  else {
    for (Index i=1; i<=N_; i++) {
      for (Index j=1; j<=N_; j++) {
        Index iu = u_index(i,j);
        grad_f[iu] = 0.;
      }
    }
  }

  return true;
}

bool MittelmannDistCntrlDiriBase::eval_g(Index n, const Number* x, bool new_x,
    Index m, Number* g)
{
  // return the value of the constraints: g(x)

  // compute the discretized PDE for each interior grid point
  for (Index i=1; i<=N_; i++) {
    for (Index j=1; j<=N_; j++) {
      Number val;

      // Start with the discretized Laplacian operator
      val = 4.* x[y_index(i,j)]
            - x[y_index(i-1,j)] - x[y_index(i+1,j)]
            - x[y_index(i,j-1)] - x[y_index(i,j+1)];

      // Add the forcing term (including the step size here)
      val += hh_*d_cont(x1_grid(i), x2_grid(j),
                        x[y_index(i,j)], x[u_index(i,j)]);
      g[pde_index(i,j)] = val;
    }
  }

  return true;
}

bool MittelmannDistCntrlDiriBase::eval_jac_g(Index n, const Number* x, bool new_x,
    Index m, Index nele_jac, Index* iRow, Index *jCol,
    Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian of the constraints

    Index ijac = 0;
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
        Index ig = pde_index(i,j);

        // y(i,j)
        iRow[ijac] = ig;
        jCol[ijac] = y_index(i,j);
        ijac++;

        // y(i-1,j)
        iRow[ijac] = ig;
        jCol[ijac] = y_index(i-1,j);
        ijac++;

        // y(i+1,j)
        iRow[ijac] = ig;
        jCol[ijac] = y_index(i+1,j);
        ijac++;

        // y(i,j-1)
        iRow[ijac] = ig;
        jCol[ijac] = y_index(i,j-1);
        ijac++;

        // y(i,j+1)
        iRow[ijac] = ig;
        jCol[ijac] = y_index(i,j+1);
        ijac++;

        // u(i,j)
        iRow[ijac] = ig;
        jCol[ijac] = u_index(i,j);
        ijac++;
      }
    }

    DBG_ASSERT(ijac==nele_jac);
  }
  else {
    // return the values of the jacobian of the constraints
    Index ijac = 0;
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
        // y(i,j)
        values[ijac] = 4. + hh_*d_cont_dy(x1_grid(i), x2_grid(j),
                                          x[y_index(i,j)], x[u_index(i,j)]);
        ijac++;

        // y(i-1,j)
        values[ijac] = -1.;
        ijac++;

        // y(i+1,j)
        values[ijac] = -1.;
        ijac++;

        // y(1,j-1)
        values[ijac] = -1.;
        ijac++;

        // y(1,j+1)
        values[ijac] = -1.;
        ijac++;

        // y(i,j)
        values[ijac] = hh_*d_cont_du(x1_grid(i), x2_grid(j),
                                     x[y_index(i,j)], x[u_index(i,j)]);
        ijac++;
      }
    }

    DBG_ASSERT(ijac==nele_jac);
  }

  return true;
}

bool
MittelmannDistCntrlDiriBase::eval_h(Index n, const Number* x, bool new_x,
                                    Number obj_factor, Index m,
                                    const Number* lambda,
                                    bool new_lambda, Index nele_hess, Index* iRow,
                                    Index* jCol, Number* values)
{
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    Index ihes=0;
    // First the diagonal entries for y(i,j)
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
        iRow[ihes] = y_index(i,j);
        jCol[ihes] = y_index(i,j);
        ihes++;
      }
    }

    if (alpha_>0.) {
      // Now the diagonal entries for u(i,j)
      for (Index i=1; i<= N_; i++) {
        for (Index j=1; j<= N_; j++) {
          iRow[ihes] = u_index(i,j);
          jCol[ihes] = u_index(i,j);
          ihes++;
        }
      }
    }

    DBG_ASSERT(ihes==nele_hess);
  }
  else {
    // return the values

    Index ihes=0;
    // First the diagonal entries for y(i,j)
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {

        // Contribution from the objective function
        values[ihes] = obj_factor*hh_;

        // Contribution from the PDE constraint
        values[ihes] += lambda[pde_index(i,j)] * hh_ *
                        d_cont_dydy(x1_grid(i), x2_grid(j),
                                    x[y_index(i,j)], x[u_index(i,j)]);

        ihes++;
      }
    }

    // Now the diagonal entries for u(i,j)
    if (alpha_>0.) {
      for (Index i=1; i<= N_; i++) {
        for (Index j=1; j<= N_; j++) {
          // Contribution from the objective function
          values[ihes] = obj_factor*hh_*alpha_;
          ihes++;
        }
      }
    }

  }

  return true;
}

void
MittelmannDistCntrlDiriBase::finalize_solution(SolverReturn status,
    Index n, const Number* x, const Number* z_L, const Number* z_U,
    Index m, const Number* g, const Number* lambda, Number obj_value,
    const IpoptData* ip_data,
    IpoptCalculatedQuantities* ip_cq)
{
  /*
  FILE* fp = fopen("solution.txt", "w+");

  for (Index i=0; i<=N_+1; i++) {
    for (Index j=0; j<=N_+1; j++) {
      fprintf(fp, "y[%6d,%6d] = %15.8e\n", i, j, x[y_index(i,j)]);
    }
  }
  for (Index i=1; i<=N_; i++) {
    for (Index j=1; j<=N_; j++) {
      fprintf(fp, "u[%6d,%6d] = %15.8e\n", i, j ,x[u_index(i,j)]);
    }
  }

  fclose(fp);
  */
}
