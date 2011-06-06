// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter             IBM    2005-10-18
//           Olaf Schenk   (Univ. of Basel)      2007-08-01
//              modified MittelmannBndryCntrlDiri.cpp for 3-dim problem

#include "MittelmannBndryCntrlDiri3D.hpp"

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
MittelmannBndryCntrlDiriBase3D::MittelmannBndryCntrlDiriBase3D()
    :
    y_d_(NULL)
{}

MittelmannBndryCntrlDiriBase3D::~MittelmannBndryCntrlDiriBase3D()
{
  delete [] y_d_;
}

void
MittelmannBndryCntrlDiriBase3D::SetBaseParameters(Index N, Number alpha, Number lb_y,
    Number ub_y, Number lb_u, Number ub_u,
    Number d_const, Number B, Number C)
{
  N_ = N;
  h_ = 1./(N+1);
  hh_ = h_*h_;
  hhh_ = hh_*h_;
  lb_y_ = lb_y;
  ub_y_ = ub_y;
  lb_u_ = lb_u;
  ub_u_ = ub_u;
  d_const_ = d_const;
  alpha_ = alpha;
  B_ = B;
  C_ = C;

  PenA_ = 1.5 - 1.125*C_/B_;
  PenB_ = 1.75*C_/pow(B_,3) - 1.5/(B_*B_);
  PenC_ = 0.5/pow(B_,4) - 0.625*C_/pow(B_,5);

  // Initialize the target profile variables
  delete [] y_d_;
  y_d_ = new Number[(N_+2)*(N_+2)*(N_+2)];
  for (Index k=0; k<= N_+1; k++) {
    for (Index j=0; j<= N_+1; j++) {
      for (Index i=0; i<= N_+1; i++) {
        y_d_[y_index(i,j,k)] = y_d_cont(x1_grid(i),x2_grid(j),x3_grid(k));
      }
    }
  }
}

bool MittelmannBndryCntrlDiriBase3D::get_nlp_info(
  Index& n, Index& m, Index& nnz_jac_g,
  Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // We for each of the N_+2 times N_+2 times N_+2  mesh points we have
  // the value of the functions y, including the control parameters on
  // the boundary
  n = (N_+2)*(N_+2)*(N_+2);

  // For each of the N_ times N_ times N_ interior mesh points we have the
  // discretized PDE.
  m = N_*N_*N_;

  // y(i,j,k), y(i-1,j,k), y(i+1,j,k), y(i,j-1,k), y(i,j+1,k),
  // y(i-1,j,k-1), y(i,j,k+1)
  // of the N_*N_*N_ discretized PDEs
  nnz_jac_g = 7*N_*N_*N_;

  // diagonal entry for each y(i,j) in the interior
  nnz_h_lag = N_*N_*N_;
  if (alpha_>0.) {
    // and one entry for u(i,j) in the bundary if alpha is not zero
    nnz_h_lag += 6*N_*N_;
  }

  // We use the C indexing style for row/col entries (corresponding to
  // the C notation, starting at 0)
  index_style = C_STYLE;

  return true;
}

bool
MittelmannBndryCntrlDiriBase3D::get_bounds_info(Index n, Number* x_l, Number* x_u,
    Index m, Number* g_l, Number* g_u)
{
  // Set overall bounds on the y variables
  for (Index i=0; i<=N_+1; i++) {
    for (Index j=0; j<=N_+1; j++) {
      for (Index k=0; k<=N_+1; k++) {
        Index iy = y_index(i,j,k);
        x_l[iy] = lb_y_;
        x_u[iy] = ub_y_;
      }
    }
  }
  // Set the overall 3D bounds on the control variables
  for (Index i=0; i<=N_+1; i++) {
    for (Index j=0; j<=N_+1; j++) {
      Index iu = y_index(i,j,0);
      x_l[iu] = lb_u_;
      x_u[iu] = ub_u_;
    }
  }
  for (Index i=0; i<=N_+1; i++) {
    for (Index j=0; j<=N_+1; j++) {
      Index iu = y_index(i,j,N_+1);
      x_l[iu] = lb_u_;
      x_u[iu] = ub_u_;
    }
  }
  for (Index i=0; i<=N_+1; i++) {
    for (Index k=0; k<=N_+1; k++) {
      Index iu = y_index(i,0,k);
      x_l[iu] = lb_u_;
      x_u[iu] = ub_u_;
    }
  }
  for (Index i=0; i<=N_+1; i++) {
    for (Index k=0; k<=N_+1; k++) {
      Index iu = y_index(i,N_+1,k);
      x_l[iu] = lb_u_;
      x_u[iu] = ub_u_;
    }
  }
  for (Index j=0; j<=N_+1; j++) {
    for (Index k=0; k<=N_+1; k++) {
      Index iu = y_index(0,j,k);
      x_l[iu] = lb_u_;
      x_u[iu] = ub_u_;
    }
  }
  for (Index j=0; j<=N_+1; j++) {
    for (Index k=0; k<=N_+1; k++) {
      Index iu = y_index(N_+1,j,k);
      x_l[iu] = lb_u_;
      x_u[iu] = ub_u_;
    }
  }

  // The values of y on the corners doens't appear anywhere, so we fix
  // them to zero
  for (Index j=0; j<=N_+1; j++) {
    x_l[y_index(0,j,0)] = x_u[y_index(0,j,0)] = 0.;
    x_l[y_index(0,j,N_+1)] = x_u[y_index(0,j,N_+1)] = 0.;
    x_l[y_index(N_+1,j,0)] = x_u[y_index(N_+1,j,0)] = 0.;
    x_l[y_index(N_+1,j,N_+1)] = x_u[y_index(N_+1,j,N_+1)] = 0.;
  }
  for (Index k=0; k<=N_+1; k++) {
    x_l[y_index(0,0,k)] = x_u[y_index(0,0,k)] = 0.;
    x_l[y_index(0,N_+1,k)] = x_u[y_index(0,N_+1,k)] = 0.;
    x_l[y_index(N_+1,0,k)] = x_u[y_index(N_+1,0,k)] = 0.;
    x_l[y_index(N_+1,N_+1,k)] = x_u[y_index(N_+1,N_+1,k)] = 0.;
  }
  for (Index i=0; i<=N_+1; i++) {
    x_l[y_index(i,0,0)] = x_u[y_index(i,0,0)] = 0.;
    x_l[y_index(i,0,N_+1)] = x_u[y_index(i,0,N_+1)] = 0.;
    x_l[y_index(i,N_+1,0)] = x_u[y_index(i,N_+1,0)] = 0.;
    x_l[y_index(i,N_+1,N_+1)] = x_u[y_index(i,N_+1,N_+1)] = 0.;
  }


  // all discretized PDE constraints have right hand side equal to
  // minus the constant value of the function d
  for (Index i=0; i<m; i++) {
    g_l[i] = -hh_*d_const_;
    g_u[i] = -hh_*d_const_;
  }

  return true;
}

bool
MittelmannBndryCntrlDiriBase3D::get_starting_point(Index n, bool init_x, Number* x,
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
      for (Index k=0; k<= N_+1; k++) {
        x[y_index(i,j,k)] = y_d_[y_index(i,j,k)]; // 0 in AMPL model
      }
    }
  }

  Number umid = (ub_u_ + lb_u_)/2.;
  for (Index i=0; i<=N_+1; i++) {
    for (Index j=0; j<=N_+1; j++) {
      Index iu = y_index(i,j,0);
      x[iu] = umid;
    }
  }
  for (Index i=0; i<=N_+1; i++) {
    for (Index j=0; j<=N_+1; j++) {
      Index iu = y_index(i,j,N_+1);
      x[iu] = umid;
    }
  }
  for (Index i=0; i<=N_+1; i++) {
    for (Index k=0; k<=N_+1; k++) {
      Index iu = y_index(i,0,k);
      x[iu] = umid;
    }
  }
  for (Index i=0; i<=N_+1; i++) {
    for (Index k=0; k<=N_+1; k++) {
      Index iu = y_index(i,N_+1,k);
      x[iu] = umid;
    }
  }
  for (Index k=0; k<=N_+1; k++) {
    for (Index j=0; j<=N_+1; j++) {
      Index iu = y_index(0,j,k);
      x[iu] = umid;
    }
  }
  for (Index k=0; k<=N_+1; k++) {
    for (Index j=0; j<=N_+1; j++) {
      Index iu = y_index(N_+1,j,k);
      x[iu] = umid;
    }
  }


  return true;
}

bool
MittelmannBndryCntrlDiriBase3D::get_scaling_parameters(Number& obj_scaling,
    bool& use_x_scaling, Index n, Number* x_scaling,
    bool& use_g_scaling, Index m, Number* g_scaling)
{
  obj_scaling = 1./hhh_;
  use_x_scaling = false;
  use_g_scaling = false;
  return true;
}

bool
MittelmannBndryCntrlDiriBase3D::eval_f(Index n, const Number* x,
                                       bool new_x, Number& obj_value)
{
  // return the value of the objective function
  obj_value = 0.;

  // First the integration of y-td over the interior
  for (Index i=1; i<=N_; i++) {
    for (Index j=1; j<= N_; j++) {
      for (Index k=1; k<= N_; k++) {
        Index iy = y_index(i,j,k);
        Number tmp = x[iy] - y_d_[iy];
        obj_value += PenObj(tmp);
        //obj_value += tmp*tmp;
      }
    }
  }
  obj_value *= hhh_;

  // Now the integration of u over the boundary
  if (alpha_>0.) {
    Number usum = 0.;

    for (Index i=1; i<=N_; i++) {
      for (Index j=1; j<=N_; j++) {
        Index iu = y_index(i,j,0);
        usum += x[iu]*x[iu];
      }
    }
    for (Index i=1; i<=N_; i++) {
      for (Index j=1; j<=N_; j++) {
        Index iu = y_index(i,j,N_+1);
        usum += x[iu]*x[iu];
      }
    }

    for (Index k=1; k<=N_; k++) {
      for (Index j=1; j<=N_; j++) {
        Index iu = y_index(0,j,k);
        usum += x[iu]*x[iu];
      }
    }
    for (Index k=1; k<=N_; k++) {
      for (Index j=1; j<=N_; j++) {
        Index iu = y_index(N_+1,j,k);
        usum += x[iu]*x[iu];
      }
    }
    for (Index i=1; i<=N_; i++) {
      for (Index k=1; k<=N_; k++) {
        Index iu = y_index(i,0,k);
        usum += x[iu]*x[iu];
      }
    }
    for (Index i=1; i<=N_; i++) {
      for (Index k=1; k<=N_; k++) {
        Index iu = y_index(i,N_+1,k);
        usum += x[iu]*x[iu];
      }
    }

    obj_value += alpha_*hh_*0.5*usum;
  }

  return true;
}

bool
MittelmannBndryCntrlDiriBase3D::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // return the gradient of the objective function grad_{x} f(x)

  // now let's take care of the nonzero values coming from the
  // integrant over the interior
  for (Index i=1; i<=N_; i++) {
    for (Index j=1; j<= N_; j++) {
      for (Index k=1; k<= N_; k++) {
        Index iy = y_index(i,j,k);
        grad_f[iy] = hhh_*PenObj_1(x[iy] - y_d_[iy]);
      }
    }
  }

  // The values for variables on the boundary
  if (alpha_>0.) {
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
        Index iu = y_index(i,j,0);
        grad_f[iu] = alpha_*hh_*x[iu];
      }
    }
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
        Index iu = y_index(i,j,N_+1);
        grad_f[iu] = alpha_*hh_*x[iu];
      }
    }
    for (Index k=1; k<= N_; k++) {
      for (Index j=1; j<= N_; j++) {
        Index iu = y_index(0,j,k);
        grad_f[iu] = alpha_*hh_*x[iu];
      }
    }
    for (Index k=1; k<= N_; k++) {
      for (Index j=1; j<= N_; j++) {
        Index iu = y_index(N_+1,j,k);
        grad_f[iu] = alpha_*hh_*x[iu];
      }
    }
    for (Index i=1; i<= N_; i++) {
      for (Index k=1; k<= N_; k++) {
        Index iu = y_index(i,0,k);
        grad_f[iu] = alpha_*hh_*x[iu];
      }
    }
    for (Index i=1; i<= N_; i++) {
      for (Index k=1; k<= N_; k++) {
        Index iu = y_index(i,N_+1,k);
        grad_f[iu] = alpha_*hh_*x[iu];
      }
    }

  }
  else {
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
        grad_f[y_index(i,j,0)] = 0.;
      }
    }
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
        grad_f[y_index(i,j,N_+1)] = 0.;
      }
    }

    for (Index k=1; k<= N_; k++) {
      for (Index j=1; j<= N_; j++) {
        grad_f[y_index(0,j,k)] = 0.;
      }
    }
    for (Index k=1; k<= N_; k++) {
      for (Index j=1; j<= N_; j++) {
        grad_f[y_index(N_+1,j,k)] = 0.;
      }
    }

    for (Index i=1; i<= N_; i++) {
      for (Index k=1; k<= N_; k++) {
        grad_f[y_index(i,0,k)] = 0.;
      }
    }
    for (Index i=1; i<= N_; i++) {
      for (Index k=1; k<= N_; k++) {
        grad_f[y_index(i,N_+1,k)] = 0.;
      }
    }
  }

  // Nothing on the corner points
  for (Index i=0; i<= N_+1; i++) {
    grad_f[y_index(i,0   ,0   )] = 0.;
    grad_f[y_index(i,0   ,N_+1)] = 0.;
    grad_f[y_index(i,N_+1,0   )] = 0.;
    grad_f[y_index(i,N_+1,N_+1)] = 0.;
  }
  for (Index j=0; j<= N_+1; j++) {
    grad_f[y_index(0,   j,0   )] = 0.;
    grad_f[y_index(N_+1,j,0   )] = 0.;
    grad_f[y_index(0,   j,N_+1)] = 0.;
    grad_f[y_index(N_+1,j,N_+1)] = 0.;
  }
  for (Index k=0; k<= N_+1; k++) {
    grad_f[y_index(0,   0,   k)] = 0.;
    grad_f[y_index(N_+1,0,   k)] = 0.;
    grad_f[y_index(0,   N_+1,k)] = 0.;
    grad_f[y_index(N_+1,N_+1,k)] = 0.;
  }


  return true;
}

bool MittelmannBndryCntrlDiriBase3D::eval_g(Index n, const Number* x, bool new_x,
    Index m, Number* g)
{
  // return the value of the constraints: g(x)

  // compute the discretized PDE for each interior grid point
  Index ig = 0;
  for (Index i=1; i<=N_; i++) {
    for (Index j=1; j<=N_; j++) {
      for (Index k=1; k<=N_; k++) {
        Number val;

        // Start with the discretized Laplacian operator
        val = 6.* x[y_index(i,j,k)]
              - x[y_index(i-1,j,k)] - x[y_index(i+1,j,k)]
              - x[y_index(i,j-1,k)] - x[y_index(i,j+1,k)]
              - x[y_index(i,j,k-1)] - x[y_index(i,j,k+1)];

        g[ig] = val;
        ig++;
      }
    }
  }

  DBG_ASSERT(ig==m);

  return true;
}

bool MittelmannBndryCntrlDiriBase3D::eval_jac_g(Index n, const Number* x, bool new_x,
    Index m, Index nele_jac, Index* iRow, Index *jCol,
    Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian of the constraints

    Index ijac = 0;
    Index ig = 0;
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
        for (Index k=1; k<= N_; k++) {

          // y(i,j,k)
          iRow[ijac] = ig;
          jCol[ijac] = y_index(i,j,k);
          ijac++;

          // y(i-1,j,k)
          iRow[ijac] = ig;
          jCol[ijac] = y_index(i-1,j,k);
          ijac++;

          // y(i+1,j,k)
          iRow[ijac] = ig;
          jCol[ijac] = y_index(i+1,j,k);
          ijac++;

          // y(i,j-1,k)
          iRow[ijac] = ig;
          jCol[ijac] = y_index(i,j-1,k);
          ijac++;

          // y(i,j+1,k)
          iRow[ijac] = ig;
          jCol[ijac] = y_index(i,j+1,k);
          ijac++;

          // y(i,j,k-1)
          iRow[ijac] = ig;
          jCol[ijac] = y_index(i,j,k-1);
          ijac++;

          // y(i,j,k+1)
          iRow[ijac] = ig;
          jCol[ijac] = y_index(i,j,k+1);
          ijac++;

          ig++;
        }
      }
    }

    DBG_ASSERT(ijac==nele_jac);
  }
  else {
    // return the values of the jacobian of the constraints
    Index ijac = 0;
    for (Index i=1; i<= N_; i++) {
      for (Index j=1; j<= N_; j++) {
        for (Index k=1; k<= N_; k++) {
          // y(i,j,k)
          values[ijac] = 6.;
          ijac++;

          // y(i-1,j,k)
          values[ijac] = -1.;
          ijac++;

          // y(i+1,j,k)
          values[ijac] = -1.;
          ijac++;

          // y(1,j-1,k)
          values[ijac] = -1.;
          ijac++;

          // y(1,j+1,k)
          values[ijac] = -1.;
          ijac++;

          // y(1,j,k-1)
          values[ijac] = -1.;
          ijac++;

          // y(1,j,k+1)
          values[ijac] = -1.;
          ijac++;
        }
      }
    }

    DBG_ASSERT(ijac==nele_jac);
  }

  return true;
}

bool
MittelmannBndryCntrlDiriBase3D::eval_h(Index n, const Number* x, bool new_x,
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
        for (Index k=1; k<= N_; k++) {
          iRow[ihes] = y_index(i,j,k);
          jCol[ihes] = y_index(i,j,k);
          ihes++;
        }
      }
    }

    if (alpha_>0.) {
      // Now the diagonal entries for u at the boundary
      for (Index i=1; i<=N_; i++) {
        for (Index j=1; j<=N_; j++) {
          Index iu = y_index(i,j,0);
          iRow[ihes] = iu;
          jCol[ihes] = iu;
          ihes++;
        }
      }
      for (Index i=1; i<=N_; i++) {
        for (Index j=1; j<=N_; j++) {
          Index iu = y_index(i,j,N_+1);
          iRow[ihes] = iu;
          jCol[ihes] = iu;
          ihes++;
        }
      }
      for (Index k=1; k<=N_; k++) {
        for (Index j=1; j<=N_; j++) {
          Index iu = y_index(0,j,k);
          iRow[ihes] = iu;
          jCol[ihes] = iu;
          ihes++;
        }
      }
      for (Index k=1; k<=N_; k++) {
        for (Index j=1; j<=N_; j++) {
          Index iu = y_index(N_+1,j,k);
          iRow[ihes] = iu;
          jCol[ihes] = iu;
          ihes++;
        }
      }
      for (Index i=1; i<=N_; i++) {
        for (Index k=1; k<=N_; k++) {
          Index iu = y_index(i,0,k);
          iRow[ihes] = iu;
          jCol[ihes] = iu;
          ihes++;
        }
      }
      for (Index i=1; i<=N_; i++) {
        for (Index k=1; k<=N_; k++) {
          Index iu = y_index(i,N_+1,k);
          iRow[ihes] = iu;
          jCol[ihes] = iu;
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
        for (Index k=1; k<= N_; k++) {
          // Contribution from the objective function
          Index iy = y_index(i,j,k);
          values[ihes] = obj_factor*hhh_*PenObj_2(x[iy] - y_d_[iy]);

          ihes++;
        }
      }
    }

    // Now the diagonal entries for u(i,j)
    if (alpha_>0.) {
      // Now the diagonal entries for u at the boundary
      for (Index i=1; i<=N_; i++) {
        for (Index j=1; j<=N_; j++) {
          values[ihes] = obj_factor*hh_*alpha_;
          ihes++;
        }
      }
      for (Index i=1; i<=N_; i++) {
        for (Index j=1; j<=N_; j++) {
          values[ihes] = obj_factor*hh_*alpha_;
          ihes++;
        }
      }
      for (Index i=1; i<=N_; i++) {
        for (Index j=1; j<=N_; j++) {
          values[ihes] = obj_factor*hh_*alpha_;
          ihes++;
        }
      }
      for (Index i=1; i<=N_; i++) {
        for (Index j=1; j<=N_; j++) {
          values[ihes] = obj_factor*hh_*alpha_;
          ihes++;
        }
      }

      for (Index i=1; i<=N_; i++) {
        for (Index j=1; j<=N_; j++) {
          values[ihes] = obj_factor*hh_*alpha_;
          ihes++;
        }
      }
      for (Index i=1; i<=N_; i++) {
        for (Index j=1; j<=N_; j++) {
          values[ihes] = obj_factor*hh_*alpha_;
          ihes++;
        }
      }
    }
  }

  return true;
}

void
MittelmannBndryCntrlDiriBase3D::finalize_solution(SolverReturn status,
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

  fclose(fp);
  */
}
