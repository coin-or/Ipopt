// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2005-10-18
//                     based on MyNLP.cpp

#include "MittelmannBndryCntrlNeum.hpp"

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
MittelmannBndryCntrlNeumBase::MittelmannBndryCntrlNeumBase()
    :
    y_d_(NULL)
{}

MittelmannBndryCntrlNeumBase::~MittelmannBndryCntrlNeumBase()
{
  delete [] y_d_;
}

void
MittelmannBndryCntrlNeumBase::SetBaseParameters(Index N, Number alpha,
    Number lb_y, Number ub_y,
    Number lb_u, Number ub_u,
    Number u_init)
{
  N_ = N;
  h_ = 1./(N+1);
  hh_ = h_*h_;
  alpha_ = alpha;
  lb_y_ = lb_y;
  ub_y_ = ub_y;
  lb_u_ = lb_u;
  ub_u_ = ub_u;
  u_init_ = u_init;

  // Initialize the target profile variables
  delete [] y_d_;
  y_d_ = new Number[(N_+2)*(N_+2)];
  for (Index j=0; j<= N_+1; j++) {
    for (Index i=0; i<= N_+1; i++) {
      y_d_[y_index(i,j)] = y_d_cont(x1_grid(i),x2_grid(j));
    }
  }
}

bool MittelmannBndryCntrlNeumBase::get_nlp_info(
  Index& n, Index& m, Index& nnz_jac_g,
  Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // We for each of the N_+2 times N_+2 mesh points we have the value
  // of the functions y, and for each 4*N_ boundary mesh points we
  // have values for u
  n = (N_+2)*(N_+2) + 4*N_;

  // For each of the N_ times N_ interior mesh points we have the
  // discretized PDE, and we have one constriant for each boundary
  // point (except for the corners)
  m = N_*N_ + 4*N_;

  // y(i,j), y(i-1,j), y(i+1,j), y(i,j-1), y(i,j+1) for each of the
  // N_*N_ discretized PDEs, and for the Neumann boundary conditions
  // we have entries for two y's and one u
  nnz_jac_g = 5*N_*N_ + 3*4*N_;

  // diagonal entry for each dydy, dudu, dydu in the interior
  nnz_h_lag = N_*N_;
  if (!b_cont_dydy_alwayszero()) {
    nnz_h_lag += 4*N_;
  }
  if (alpha_!=0.) {
    nnz_h_lag += 4*N_;
  }

  // We use the C indexing style for row/col entries (corresponding to
  // the C notation, starting at 0)
  index_style = C_STYLE;

  return true;
}

bool
MittelmannBndryCntrlNeumBase::get_bounds_info(Index n, Number* x_l, Number* x_u,
    Index m, Number* g_l, Number* g_u)
{
  // Set overall bounds on the y variables
  for (Index i=0; i<=N_+1; i++) {
    for (Index j=0; j<=N_+1; j++) {
      Index iy = y_index(i,j);
      x_l[iy] = lb_y_;
      x_u[iy] = ub_y_;
    }
  }

  // Set overall bounds on the u variables
  for (Index j=1; j<=N_; j++) {
    Index iu = u0j_index(j);
    x_l[iu] = lb_u_;
    x_u[iu] = ub_u_;
  }
  for (Index j=1; j<=N_; j++) {
    Index iu = u1j_index(j);
    x_l[iu] = lb_u_;
    x_u[iu] = ub_u_;
  }
  for (Index i=1; i<=N_; i++) {
    Index iu = ui0_index(i);
    x_l[iu] = lb_u_;
    x_u[iu] = ub_u_;
  }
  for (Index i=1; i<=N_; i++) {
    Index iu = ui1_index(i);
    x_l[iu] = lb_u_;
    x_u[iu] = ub_u_;
  }

  // There is no information for the y's at the corner points, so just
  // take those variables out
  x_l[y_index(0,0)] = x_u[y_index(0,0)] = 0.;
  x_l[y_index(0,N_+1)] = x_u[y_index(0,N_+1)] = 0.;
  x_l[y_index(N_+1,0)] = x_u[y_index(N_+1,0)] = 0.;
  x_l[y_index(N_+1,N_+1)] = x_u[y_index(N_+1,N_+1)] = 0.;

  // all discretized PDE constraints have right hand side zero
  for (Index i=0; i<m; i++) {
    g_l[i] = 0.;
    g_u[i] = 0.;
  }

  return true;
}

bool
MittelmannBndryCntrlNeumBase::get_starting_point(Index n, bool init_x, Number* x,
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
      //x[y_index(i,j)] += h_*x1_grid(i) + 2*h_*x2_grid(j);
    }
  }

  // Set the initial (constant) value for the u's
  for (Index j=1; j<= N_; j++) {
    x[u0j_index(j)] = u_init_;
    x[u1j_index(j)] = u_init_;
  }
  for (Index i=1; i<= N_; i++) {
    x[ui0_index(i)] = u_init_;
    x[ui1_index(i)] = u_init_;
  }

  return true;
}

bool
MittelmannBndryCntrlNeumBase::get_scaling_parameters(Number& obj_scaling,
    bool& use_x_scaling,
    Index n, Number* x_scaling,
    bool& use_g_scaling,
    Index m, Number* g_scaling)
{
  obj_scaling = 1./hh_;
  use_x_scaling = false;
  use_g_scaling = false;
  return true;
}

bool
MittelmannBndryCntrlNeumBase::eval_f(Index n, const Number* x,
                                     bool new_x, Number& obj_value)
{
  // return the value of the objective function
  obj_value = 0.;

  // First the integration of y-td over the interior
  for (Index i=1; i<=N_; i++) {
    for (Index j=1; j<= N_; j++) {
      Index iy = y_index(i,j);
      Number tmp = x[iy] - y_d_[iy];
      obj_value += tmp*tmp;
    }
  }
  obj_value *= hh_/2.;

  // Now the integration of u over the boundary
  if (alpha_!=0.) {
    Number usum = 0.;
    for (Index j=1; j<=N_; j++) {
      Index iu = u0j_index(j);
      usum += x[iu]*x[iu];
    }
    for (Index j=1; j<=N_; j++) {
      Index iu = u1j_index(j);
      usum += x[iu]*x[iu];
    }
    for (Index i=1; i<=N_; i++) {
      Index iu = ui0_index(i);
      usum += x[iu]*x[iu];
    }
    for (Index i=1; i<=N_; i++) {
      Index iu = ui1_index(i);
      usum += x[iu]*x[iu];
    }
    obj_value += alpha_*h_/2.*usum;
  }

  return true;
}

bool
MittelmannBndryCntrlNeumBase::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // return the gradient of the objective function grad_{x} f(x)

  // now let's take care of the nonzero values coming from the
  // integrant over the interior
  for (Index i=1; i<=N_; i++) {
    for (Index j=1; j<= N_; j++) {
      Index iy = y_index(i,j);
      grad_f[iy] = hh_*(x[iy] - y_d_[iy]);
    }
  }

  // The values for variables on the boundary
  if (alpha_!=0.) {
    for (Index j=1; j<= N_; j++) {
      Index iu = u0j_index(j);
      grad_f[iu] = alpha_*h_*x[iu];
    }
    for (Index j=1; j<= N_; j++) {
      Index iu = u1j_index(j);
      grad_f[iu] = alpha_*h_*x[iu];
    }
    for (Index i=1; i<= N_; i++) {
      Index iu = ui0_index(i);
      grad_f[iu] = alpha_*h_*x[iu];
    }
    for (Index i=1; i<= N_; i++) {
      Index iu = ui1_index(i);
      grad_f[iu] = alpha_*h_*x[iu];
    }
  }
  else {
    for (Index j=1; j<= N_; j++) {
      Index iu = u0j_index(j);
      grad_f[iu] = 0.;
    }
    for (Index j=1; j<= N_; j++) {
      Index iu = u1j_index(j);
      grad_f[iu] = 0.;
    }
    for (Index i=1; i<= N_; i++) {
      Index iu = ui0_index(i);
      grad_f[iu] = 0.;
    }
    for (Index i=1; i<= N_; i++) {
      Index iu = ui1_index(i);
      grad_f[iu] = 0.;
    }
  }

  // The values are zero for y variables on the boundary
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

  return true;
}

bool MittelmannBndryCntrlNeumBase::eval_g(Index n, const Number* x, bool new_x,
    Index m, Number* g)
{
  // return the value of the constraints: g(x)

  // compute the discretized PDE for each interior grid point
  Index ig=0;
  for (Index i=1; i<=N_; i++) {
    for (Index j=1; j<=N_; j++) {
      Number val;

      // Start with the discretized Laplacian operator
      val = 4.* x[y_index(i,j)]
            - x[y_index(i-1,j)] - x[y_index(i+1,j)]
            - x[y_index(i,j-1)] - x[y_index(i,j+1)];

      // Add the forcing term (including the step size here)
      val += hh_*d_cont(x1_grid(i), x2_grid(j), x[y_index(i,j)]);
      g[ig] = val;
      ig++;
    }
  }

  // set up the Neumann boundary conditions
  for (Index j=1; j<= N_; j++) {
    g[ig] = x[y_index(0,j)] - x[y_index(1,j)]
            - h_*b_cont(x1_grid(0), x2_grid(j), x[y_index(0,j)], x[u0j_index(j)]);
    ig++;
  }
  for (Index j=1; j<= N_; j++) {
    g[ig] = x[y_index(N_+1,j)] - x[y_index(N_,j)]
            - h_*b_cont(x1_grid(N_+1), x2_grid(j), x[y_index(N_+1,j)], x[u1j_index(j)]);
    ig++;
  }
  for (Index i=1; i<= N_; i++) {
    g[ig] = x[y_index(i,0)] - x[y_index(i,1)]
            - h_*b_cont(x1_grid(i), x2_grid(0), x[y_index(i,0)], x[ui0_index(i)]);
    ig++;
  }
  for (Index i=1; i<= N_; i++) {
    g[ig] = x[y_index(i,N_+1)] - x[y_index(i,N_)]
            - h_*b_cont(x1_grid(i), x2_grid(N_+1), x[y_index(i,N_+1)], x[ui1_index(i)]);
    ig++;
  }

  DBG_ASSERT(ig==m);

  return true;
}

bool MittelmannBndryCntrlNeumBase::eval_jac_g(Index n, const Number* x, bool new_x,
    Index m, Index nele_jac, Index* iRow, Index *jCol,
    Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian of the constraints

    // distretized PDEs
    Index ijac = 0;
    Index ig = 0;
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
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
        ig++;
      }
    }

    // set up the Neumann boundary conditions
    for (Index j=1; j<=N_; j++) {
      iRow[ijac] = ig;
      jCol[ijac] = y_index(0,j);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = y_index(1,j);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = u0j_index(j);
      ijac++;
      ig++;
    }
    for (Index j=1; j<=N_; j++) {
      iRow[ijac] = ig;
      jCol[ijac] = y_index(N_,j);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = y_index(N_+1,j);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = u1j_index(j);
      ijac++;
      ig++;
    }
    for (Index i=1; i<=N_; i++) {
      iRow[ijac] = ig;
      jCol[ijac] = y_index(i,0);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = y_index(i,1);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = ui0_index(i);
      ijac++;
      ig++;
    }
    for (Index i=1; i<=N_; i++) {
      iRow[ijac] = ig;
      jCol[ijac] = y_index(i,N_);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = y_index(i,N_+1);
      ijac++;
      iRow[ijac] = ig;
      jCol[ijac] = ui1_index(i);
      ijac++;
      ig++;
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
                                          x[y_index(i,j)]);
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
      }
    }

    for (Index j=1; j<=N_; j++) {
      values[ijac] = 1. -
                     h_*b_cont_dy(x1_grid(0), x2_grid(j), x[y_index(0,j)], x[u0j_index(j)]);
      ijac++;
      values[ijac] = -1.;
      ijac++;
      values[ijac] = -h_*b_cont_du(x1_grid(0), x2_grid(j), x[y_index(0,j)], x[u0j_index(j)]);
      ijac++;
    }
    for (Index j=1; j<=N_; j++) {
      values[ijac] = -1.;
      ijac++;
      values[ijac] = 1. -
                     h_*b_cont_dy(x1_grid(N_+1), x2_grid(j), x[y_index(N_+1,j)], x[u1j_index(j)]);
      ijac++;
      values[ijac] = -h_*b_cont_du(x1_grid(N_+1), x2_grid(j), x[y_index(N_+1,j)], x[u1j_index(j)]);
      ijac++;
    }
    for (Index i=1; i<=N_; i++) {
      values[ijac] = 1. -
                     h_*b_cont_dy(x1_grid(i), x2_grid(0), x[y_index(i,0)], x[ui0_index(i)]);
      ijac++;
      values[ijac] = -1.;
      ijac++;
      values[ijac] = -h_*b_cont_du(x1_grid(i), x2_grid(0), x[y_index(i,0)], x[ui0_index(i)]);
      ijac++;
    }
    for (Index i=1; i<=N_; i++) {
      values[ijac] = -1.;
      ijac++;
      values[ijac] = 1. -
                     h_*b_cont_dy(x1_grid(i), x2_grid(N_+1), x[y_index(i,N_+1)], x[ui1_index(i)]);
      ijac++;
      values[ijac] = -h_*b_cont_du(x1_grid(i), x2_grid(N_+1), x[y_index(i,N_+1)], x[ui1_index(i)]);
      ijac++;
    }

    DBG_ASSERT(ijac==nele_jac);
  }

  return true;
}

bool
MittelmannBndryCntrlNeumBase::eval_h(Index n, const Number* x, bool new_x,
                                     Number obj_factor, Index m,
                                     const Number* lambda,
                                     bool new_lambda, Index nele_hess, Index* iRow,
                                     Index* jCol, Number* values)
{
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    Index ihes=0;
    // First the diagonal entries for dydy in the interior
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
        iRow[ihes] = y_index(i,j);
        jCol[ihes] = y_index(i,j);
        ihes++;
      }
    }

    // Now, if necessary, the dydy entries on the boundary
    if (!b_cont_dydy_alwayszero()) {
      // Now the diagonal entries for dudu
      for (Index j=1; j<= N_; j++) {
        iRow[ihes] = y_index(0,j);
        jCol[ihes] = y_index(0,j);
        ihes++;
      }
      for (Index j=1; j<= N_; j++) {
        iRow[ihes] = y_index(N_+1,j);
        jCol[ihes] = y_index(N_+1,j);
        ihes++;
      }
      for (Index i=1; i<= N_; i++) {
        iRow[ihes] = y_index(i,0);
        jCol[ihes] = y_index(i,0);
        ihes++;
      }
      for (Index i=1; i<= N_; i++) {
        iRow[ihes] = y_index(i,N_+1);
        jCol[ihes] = y_index(i,N_+1);
        ihes++;
      }
    }

    if (alpha_!=0.) {
      // Now the diagonal entries for dudu
      for (Index j=1; j<=N_; j++) {
        Index iu = u0j_index(j);
        iRow[ihes] = iu;
        jCol[ihes] = iu;
        ihes++;
      }
      for (Index j=1; j<=N_; j++) {
        Index iu = u1j_index(j);
        iRow[ihes] = iu;
        jCol[ihes] = iu;
        ihes++;
      }
      for (Index i=1; i<= N_; i++) {
        Index iu = ui0_index(i);
        iRow[ihes] = iu;
        jCol[ihes] = iu;
        ihes++;
      }
      for (Index i=1; i<= N_; i++) {
        Index iu = ui1_index(i);
        iRow[ihes] = iu;
        jCol[ihes] = iu;
        ihes++;
      }
    }

    DBG_ASSERT(ihes==nele_hess);
  }
  else {
    // return the values

    Index ihes=0;
    Index ihes_store;
    // First the diagonal entries for dydy
    ihes_store = ihes;
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
        // Contribution from the objective function
        values[ihes] = obj_factor*hh_;

        ihes++;
      }
    }
    // If we have something from the discretized PDEs, add this now
    if (!d_cont_dydy_alwayszero()) {
      Index ig = 0;
      ihes = ihes_store;
      for (Index i=1; i<=N_; i++) {
        for (Index j=1; j<=N_; j++) {
          values[ihes] += lambda[ig]*hh_*d_cont_dydy(x1_grid(i), x2_grid(j), x[y_index(i,j)]);
          ihes++;
          ig++;
        }
      }
    }

    // Now include the elements for dydy on the boudary if there are any
    if (!b_cont_dydy_alwayszero()) {
      Index ig = N_*N_;
      // Now the diagonal entries for dudu
      for (Index j=1; j<= N_; j++) {
        values[ihes] = -lambda[ig]*h_*b_cont_dydy(x1_grid(0), x2_grid(j), x[y_index(0,j)], x[u0j_index(j)]);
        ig++;
        ihes++;
      }
      for (Index j=1; j<= N_; j++) {
        values[ihes] = -lambda[ig]*h_*b_cont_dydy(x1_grid(N_+1), x2_grid(j), x[y_index(N_+1,j)], x[u1j_index(j)]);
        ig++;
        ihes++;
      }
      for (Index i=1; i<= N_; i++) {
        values[ihes] = -lambda[ig]*h_*b_cont_dydy(x1_grid(i), x2_grid(0), x[y_index(i,0)], x[ui0_index(i)]);
        ig++;
        ihes++;
      }
      for (Index i=1; i<= N_; i++) {
        values[ihes] = -lambda[ig]*h_*b_cont_dydy(x1_grid(i), x2_grid(N_+1), x[y_index(i,N_+1)], x[ui1_index(i)]);
        ig++;
        ihes++;
      }
    }

    // Finally, we take care of the dudu entries
    if (alpha_!=0.) {
      // Now the diagonal entries for u at the boundary
      for (Index i=1; i<=N_; i++) {
        values[ihes] = obj_factor*h_*alpha_;
        ihes++;
      }
      for (Index i=1; i<=N_; i++) {
        values[ihes] = obj_factor*h_*alpha_;
        ihes++;
      }
      for (Index j=1; j<=N_; j++) {
        values[ihes] = obj_factor*h_*alpha_;
        ihes++;
      }
      for (Index i=1; i<=N_; i++) {
        values[ihes] = obj_factor*h_*alpha_;
        ihes++;
      }
    }

  }

  return true;
}

void
MittelmannBndryCntrlNeumBase::finalize_solution(SolverReturn status,
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
  for (Index j=1; j<=N_; j++) {
    fprintf(fp, "u[%6d,%6d] = %15.8e\n", 0, j, x[u0j_index(j)]);
  }
  for (Index j=1; j<=N_; j++) {
    fprintf(fp, "u[%6d,%6d] = %15.8e\n", N_+1, j, x[u1j_index(j)]);
  }
  for (Index i=1; i<=N_; i++) {
    fprintf(fp, "u[%6d,%6d] = %15.8e\n", i, 0, x[ui0_index(i)]);
  }
  for (Index i=1; i<=N_; i++) {
    fprintf(fp, "u[%6d,%6d] = %15.8e\n", i, N_+1, x[ui1_index(i)]);
  }

  fclose(fp);
  */
}
