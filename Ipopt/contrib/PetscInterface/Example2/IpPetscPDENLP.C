// Copyright (C) 2010 International Business Machines.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: IpPetscPDENLP.C 1887 2011-01-28 03:17:57Z andreasw $
//
// Authors:  Johannes Huber, Andreas Waechter     IBM        2010-09-03
#include "mpi.h"
#include "IpPetscPDENLP.hpp"
#include "IpJournalist.hpp"

#include "petscvec.h"
#include "petscmat.h"

#include <iostream>

int GetProcID();

#define PRINT_LEVEL 10
//#define MY_MY_DBG_PRINT(s) {std::cout << GetProcID() << ": " << s << std::endl;}
#define MY_DBG_PRINT(s) {}
#define START_FUNCTION { if(PRINT_LEVEL>10) {std::cout << "[" << GetProcID() << "]" << __FILE__ << ":" << __LINE__ <<":" << "called " << __func__ << std::endl;} }
#define END_FUNCTION { if(PRINT_LEVEL>10) {std::cout << "[" << GetProcID() << "]" << __FILE__ << ":" << __LINE__ <<":" << "finished " << __func__ << std::endl;} }
#define CHECKPOINT { std::cout << GetProcID() << ", " << __func__ << ", Line " << __LINE__ << std::endl; }

using namespace Ipopt;

using namespace std;

PetscPDENLP::PetscPDENLP(PetscPDEBase& PetscPDE,
                             Journalist& jnlst)
    :
    PetscPDE_(&PetscPDE),
    jnlst_(&jnlst)
{}

bool
PetscPDENLP::
get_nlp_info(Index num_proc, Index proc_id,
             Index& n, Index& n_first, Index& n_last,
             Index& m, Index& m_first, Index& m_last,
             Index& nnz_jac_g_part, Index& nnz_h_lag_part,
             IndexStyleEnum& index_style)
{
  START_FUNCTION
  int m_pde;
  int m_aux;
  int n_state;
  int n_control;
  Vec petsc_vec;
  PetscInt tmp_end;
  // Distribution of constraints
  {
    PetscPDE_->getAuxConstrVector(petsc_vec);
    VecGetOwnershipRange(petsc_vec,&aux_first_,&tmp_end);
    aux_last_ = tmp_end-1;
    VecGetSize(petsc_vec,&m_aux);

    PetscPDE_->getPDEResVector(petsc_vec);
    VecGetOwnershipRange(petsc_vec,&pde_first_,&tmp_end);
    pde_last_ = tmp_end-1;
    VecGetSize(petsc_vec,&m_pde);
    
    int m_local = aux_last_-aux_first_+1 + pde_last_-pde_first_+1;
    int m_glob_end;
    MPI_Scan(&m_local, &m_glob_end, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    m_first_ = m_first = m_glob_end - m_local;
    m_last_ = m_last = m_glob_end-1;
    m = m_pde+m_aux;
  }

  // Distribution of variables (state & controls)
  {
    // states distribution
    PetscPDE_->getStateVector(petsc_vec);
    VecGetOwnershipRange(petsc_vec, &state_first_, &tmp_end);
    state_last_ = tmp_end-1;
    VecGetSize(petsc_vec,&n_state);
    
    // control distribution
    PetscPDE_->getControlVector(petsc_vec);
    VecGetOwnershipRange(petsc_vec, &control_first_, &tmp_end);
    control_last_ = tmp_end-1;
    VecGetSize(petsc_vec,&n_control);
    // Get the local offset in large x vector (ParNLP distribution)
    int n_local = state_last_-state_first_+1 + control_last_-control_first_+1;
    int n_glob_end;
    MPI_Scan(&n_local, &n_glob_end, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    n_first_ = n_first = n_glob_end - n_local;
    n_last_ = n_last = n_glob_end-1;
    n = n_state+n_control;
  }

  PetscPDE_->m_CalcType = PetscPDEBase::StructureOnly;
  Mat jac_state;
  Mat jac_control;
  PetscPDE_->calcPDE_jacobians(jac_state,jac_control);

  MatInfo info;
  MatGetInfo(jac_state,MAT_LOCAL,&info);
  nnz_jac_g_part = info.nz_used;
  MatGetInfo(jac_control,MAT_LOCAL,&info);
  nnz_jac_g_part += info.nz_used;

  Mat jac_aux_state;
  Mat jac_aux_control;
  PetscPDE_->calcAux_jacobians(jac_aux_state,jac_aux_control);
  MatGetInfo(jac_aux_state,MAT_LOCAL,&info);
  nnz_jac_g_part += info.nz_used;
  MatGetInfo(jac_aux_control,MAT_LOCAL,&info);
  nnz_jac_g_part += info.nz_used;

  nnz_h_lag_part = 0;

  Mat hss;
  Mat hcc;
  Mat hcs;
  Vec lambda_pde;
  Vec lambda_aux;

  VecCreateSeq(MPI_COMM_WORLD,m_pde,&lambda_pde);
  VecSet(lambda_pde,1.0);
  VecCreateSeq(MPI_COMM_WORLD,m_aux,&lambda_aux);
  VecSet(lambda_aux,1.0);
  PetscPDE_->calc_hessians(1.0,lambda_pde,lambda_aux,hcc,hcs,hss);
  MatGetInfo(hcc,MAT_LOCAL,&info);
  nnz_h_lag_part += info.nz_used;
  MatGetInfo(hcs,MAT_LOCAL,&info);
  nnz_h_lag_part += info.nz_used;
  MatGetInfo(hss,MAT_LOCAL,&info);
  nnz_h_lag_part += info.nz_used;
  
  VecDestroy(lambda_pde);
  VecDestroy(lambda_aux);

  index_style = C_STYLE;
  PetscPDE_->m_CalcType = PetscPDEBase::Values;
  END_FUNCTION
  return true;
}

void PetscPDENLP::optim_var_libMesh2local(Vec& state, Vec& control, Number* plocalDest)
{ // creates [local_state local_control]
  START_FUNCTION
  int iVal = 0;
  int iOffset = 0;
  PetscScalar *Vals;
  VecGetArray(state,&Vals);
  for ( iVal=0; iVal<=(state_last_-state_first_); ++iVal, ++iOffset)
    *(plocalDest+iOffset) = Vals[iVal];
  VecRestoreArray(state,&Vals);
  VecGetArray(control,&Vals);
  for ( iVal=0; iVal<=(control_last_-control_first_); ++iVal, ++iOffset)
    *(plocalDest+iOffset) = Vals[iVal];
  VecRestoreArray(control,&Vals);
  END_FUNCTION
}

void PetscPDENLP::optim_var_global2libMesh(const Number* pglobal, Vec& state, Vec& control)
{ // [state_proc0 contrl_proc0 state_proc1 contrl_proc1 .... state_procN contrl_procNn] -> local_state, local_cntrl
  START_FUNCTION
  int iVal = 0;
  int iOffset = n_first_;
  PetscScalar *Vals;
  VecGetArray(state,&Vals);
  for ( iVal=0; iVal<=(state_last_-state_first_); ++iVal, ++iOffset)
    Vals[iVal] = *(pglobal+iOffset);
  VecRestoreArray(state,&Vals);
  VecGetArray(control,&Vals);
  for ( iVal=0; iVal<=(control_last_-control_first_); ++iVal, ++iOffset)
    Vals[iVal] = *(pglobal+iOffset);
  VecRestoreArray(control,&Vals);
  assert(iOffset == n_last_+1);
  END_FUNCTION
}

bool PetscPDENLP::get_bounds_info(Index num_proc, Index proc_id,
                                    Index n, Index n_first, Index n_last,
                                    Number* x_l_part, Number* x_u_part,
                                    Index m, Index m_first, Index m_last,
                                    Number* g_l_part, Number* g_u_part)
{
  START_FUNCTION
  Vec state_l, state_u;
  Vec control_l, control_u;
  Vec aux_constr_l, aux_constr_u;
  
  {
    Vec tmp;
    PetscPDE_->getStateVector(tmp);
    VecDuplicate(tmp,&state_l);
    VecDuplicate(tmp,&state_u);
    PetscPDE_->getControlVector(tmp);
    VecDuplicate(tmp,&control_l);
    VecDuplicate(tmp,&control_u);
    PetscPDE_->getAuxConstrVector(tmp);
    VecDuplicate(tmp,&aux_constr_l);
    VecDuplicate(tmp,&aux_constr_u);
  }

  PetscPDE_->get_bounds(state_l,state_u,control_l,control_u,aux_constr_l,aux_constr_u);
  optim_var_libMesh2local(state_l, control_l, x_l_part);
  optim_var_libMesh2local(state_u, control_u, x_u_part);
  int iVal=0;
  for (iVal=0;iVal<=(pde_last_-pde_first_);++iVal) {
    g_l_part[iVal] = 0.0;
    g_u_part[iVal] = 0.0;
  }
  int iAuxVal=0;
  PetscScalar *Vals_l(NULL),*Vals_u(NULL);
  VecGetArray(aux_constr_l,&Vals_l);
  VecGetArray(aux_constr_u,&Vals_u);
  for (iAuxVal=0;iAuxVal<=(aux_last_-aux_first_);++iVal,++iAuxVal) {
    g_l_part[iVal] = Vals_l[iAuxVal];
    g_u_part[iVal] = Vals_u[iAuxVal];
  }
  VecRestoreArray(aux_constr_l,&Vals_l);
  VecRestoreArray(aux_constr_u,&Vals_u);
  
  VecDestroy(state_l);
  VecDestroy(state_u);
  VecDestroy(control_l);
  VecDestroy(control_u);
  VecDestroy(aux_constr_l);
  VecDestroy(aux_constr_u);
  END_FUNCTION
  return true;
}


bool
PetscPDENLP::get_starting_point(Index num_proc, Index proc_id,
                                  Index n, Index n_first, Index n_last,
                                  bool init_x, Number* x_part,
                                  bool init_z, Number* z_L_part, Number* z_U_part,
                                  Index m, Index m_first, Index m_last,
                                  bool init_lambda, Number* lambda_part)
{
  START_FUNCTION
  assert(!init_z && !init_lambda);

  Vec state, control, pde_residual, aux_residual;
  PetscPDE_->getStateVector(state);
  PetscPDE_->getControlVector(control);
  PetscPDE_->calcPDE_residual(pde_residual);
  PetscPDE_->calcPDE_residual(aux_residual);

  Vec state_lb_mults, state_ub_mults;
  Vec control_lb_mults, control_ub_mults;
  Vec pde_residual_mults, aux_constr_mults;

  VecDuplicate(state,&state_lb_mults);
  VecDuplicate(state,&state_ub_mults);
  VecDuplicate(control,&control_lb_mults);
  VecDuplicate(control,&control_ub_mults);
  VecDuplicate(pde_residual,&pde_residual_mults);
  VecDuplicate(aux_residual,&aux_constr_mults);

  PetscPDE_->get_starting_point(state, control, init_x, state_lb_mults,
				  state_ub_mults, control_lb_mults,
				  control_ub_mults, init_lambda,
				  pde_residual_mults, aux_constr_mults);
  optim_var_libMesh2local(state, control, x_part);
  if (init_z) {
    optim_var_libMesh2local(state_lb_mults, control_lb_mults, z_L_part);
    optim_var_libMesh2local(state_ub_mults, control_ub_mults, z_U_part);
  }
  if (init_lambda) {
    PetscScalar *Vals = NULL;
    VecGetArray(pde_residual_mults,&Vals);
    int iVal;
    for (iVal=0; iVal<=(pde_last_-pde_first_); ++iVal) {
      lambda_part[iVal] = Vals[iVal];
    }
    VecRestoreArray(pde_residual_mults,&Vals);

    VecGetArray(aux_constr_mults,&Vals);
    for (int iAuxVal=0;iAuxVal<=(aux_last_-aux_first_);++iAuxVal,++iVal) {
      lambda_part[iVal] = Vals[iAuxVal];
    }
    VecRestoreArray(aux_constr_mults,&Vals);
  }

  VecDestroy(state_lb_mults);
  VecDestroy(state_ub_mults);
  VecDestroy(control_lb_mults);
  VecDestroy(control_ub_mults);
  VecDestroy(pde_residual_mults);
  VecDestroy(aux_constr_mults);

  END_FUNCTION
  return true;
}

void PetscPDENLP::update_x(const Number* x)
{
  START_FUNCTION
  Vec state, control;
  PetscPDE_->getStateVector(state);
  PetscPDE_->getControlVector(control);
  optim_var_global2libMesh(x, state, control);
  // Need to make sure that local values are up-to-date
  PetscPDE_->update();
  END_FUNCTION
}

bool
PetscPDENLP::eval_f(Index num_proc, Index proc_id,
                      Index n, Index n_first, Index n_last,
                      const Number* x, bool new_x, Number& obj_value)
{
  START_FUNCTION
  if (new_x) {
    update_x(x);
  }
  PetscPDE_->calc_objective_part(obj_value);
  END_FUNCTION
  return true;
}


bool
PetscPDENLP::eval_grad_f(Index num_proc, Index proc_id,
                           Index n,  Index n_first, Index n_last,
                           const Number* x, bool new_x,
                           Number* grad_f_part)
{
  START_FUNCTION
  if (new_x) {
    update_x(x);
  }
  Vec grad_state, grad_control;
  PetscPDE_->calc_objective_gradient(grad_state, grad_control);
  optim_var_libMesh2local(grad_state, grad_control, grad_f_part);
  END_FUNCTION
  return true;
}

bool
PetscPDENLP::eval_g(Index num_proc, Index proc_id,
                      Index n, const Number* x, bool new_x,
                      Index m, Index m_first, Index m_last,
                      Number* g_part)
{
  // g = [Aux_Proc1 PDE_Proc_1 Aux_Proc2 PDE_Proc_2 .... Aux_ProcN PDE_ProcN]
  START_FUNCTION
  if (new_x) {
    update_x(x);
  }

  // PDE resuduals
  Vec v;
  PetscPDE_->calcPDE_residual(v);
  PetscScalar *Vals = NULL;
  VecGetArray(v,&Vals);
  int iVal;
  for (iVal=0; iVal<=(pde_last_-pde_first_); ++iVal) {
    g_part[iVal] = Vals[iVal];
  }
  VecRestoreArray(v,&Vals);

  // auxiliary constraints
  PetscPDE_->calcAux_constr(v);
  VecGetArray(v,&Vals);
  for (int iAuxVal=0;iAuxVal<=(aux_last_-aux_first_);++iAuxVal,++iVal) {
    g_part[iVal] = Vals[iAuxVal];
  }
  VecRestoreArray(v,&Vals);
  END_FUNCTION
  return true;
}

bool
PetscPDENLP::eval_jac_g(Index num_proc, Index proc_id,
                          Index n, const Number* x, bool new_x,
                          Index m, Index m_first, Index m_last,
                          Index nele_jac_part, Index* iRow_part,
                          Index *jCol_part, Number* values_part)
{
  START_FUNCTION
  if (new_x) {
    update_x(x);
  }
  if (values_part==NULL) {
    PetscPDE_->m_CalcType = PetscPDEBase::StructureOnly;
    assert(iRow_part && jCol_part);
  }
  else {
    PetscPDE_->m_CalcType = PetscPDEBase::Values;
    assert(!iRow_part && !jCol_part);
  }

  Mat mat_state, mat_control;
  Mat mat_aux_state, mat_aux_control;

  PetscPDE_->calcPDE_jacobians(mat_state, mat_control);
  PetscPDE_->calcAux_jacobians(mat_aux_state, mat_aux_control);

  int iOffset = 0;
  if (values_part==NULL) {
    int *p_state_first = new int[num_proc];
    int *p_state_last = new int[num_proc];
    int *p_control_first = new int[num_proc];
    int *p_control_last = new int[num_proc];
    int *p_optvar_first = new int[num_proc];
    int *p_loc_control_offset = new int[num_proc];
    MPI_Allgather(&state_first_,1,MPI_INT,p_state_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&state_last_ ,1,MPI_INT,p_state_last ,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&control_first_,1,MPI_INT,p_control_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&control_last_ ,1,MPI_INT,p_control_last ,1,MPI_INT,MPI_COMM_WORLD);

    MY_DBG_PRINT( "p_state_first[0]:"<<p_state_first[0] );
    MY_DBG_PRINT( "p_state_last[0]:"<<p_state_last[0] );
    MY_DBG_PRINT( "p_control_first[0]:"<<p_control_first[0] );
    MY_DBG_PRINT( "p_control_last[0]:"<<p_control_last[0] );

    int i_var(0);
    for (int i_proc=0;i_proc<num_proc;++i_proc) {
      p_optvar_first[i_proc] = i_var;
      i_var += (p_state_last[i_proc]-p_state_first[i_proc]+1) +
               (p_control_last[i_proc]-p_control_first[i_proc]+1);
      p_loc_control_offset[i_proc] = p_state_last[i_proc] - p_state_first[i_proc] + 1;
    }

    int prev_proc(0);
    // PDE
    PetscInt irow;
    for (irow = pde_first_; irow<=pde_last_; ++irow) {
      PetscInt ncols;
      const PetscInt* cols;
      // state
      MatGetRow(mat_state, irow, &ncols, &cols, PETSC_NULL);
      for (int i=0; i<ncols; ++i) {
        iRow_part[iOffset] = (irow-pde_first_);
        unsigned int i_var_glob = GlobStateControlIdx2OptVarIdx(cols[i],num_proc,p_state_first,p_state_last,
                                  p_optvar_first,NULL,prev_proc);
        jCol_part[iOffset] = i_var_glob;//n_first_state + (cols[i]-state_first_);
        iOffset++;
      }
      MatRestoreRow(mat_state,irow,&ncols,&cols,PETSC_NULL);
      // control
      MatGetRow(mat_control, irow, &ncols, &cols, PETSC_NULL);
      for (int i=0; i<ncols; ++i) {
        iRow_part[iOffset] = (irow-pde_first_);
        unsigned int i_var_glob = GlobStateControlIdx2OptVarIdx(cols[i],num_proc,
                                  p_control_first,p_control_last,p_optvar_first,p_loc_control_offset,prev_proc);
        jCol_part[iOffset] = i_var_glob;//n_first_control + (cols[i]-control_first_);
        iOffset++;
      }
      MatRestoreRow(mat_control,irow,&ncols,&cols,PETSC_NULL);
    }
    // auxiliary constraints
    PetscInt msz,nsz;
    MatGetSize(mat_aux_state,&msz,&nsz);
    const int m_loc_aux_offset = pde_last_-pde_first_+1;
    for (irow=aux_first_; irow<=aux_last_; ++irow) {
      PetscInt ncols;
      const PetscInt* cols;
      // state
      MatGetRow(mat_aux_state, irow, &ncols, &cols, PETSC_NULL);
      for (int i=0; i<ncols; ++i) {
        iRow_part[iOffset] = m_loc_aux_offset + (irow-aux_first_);
        unsigned int i_var_glob = GlobStateControlIdx2OptVarIdx(cols[i],num_proc,p_state_first,p_state_last,
                                  p_optvar_first,NULL,prev_proc);
        jCol_part[iOffset] = i_var_glob;//n_first_state + (cols[i]-state_first_);
        iOffset++;
      }
      MatRestoreRow(mat_aux_state,irow,&ncols,&cols,PETSC_NULL);
      // control
      MatGetRow(mat_aux_control, irow, &ncols, &cols, PETSC_NULL);
      for (int i=0; i<ncols; ++i) {
        iRow_part[iOffset] = m_loc_aux_offset + (irow-aux_first_);
        unsigned int i_var_glob = GlobStateControlIdx2OptVarIdx(cols[i],num_proc,
                                  p_control_first,p_control_last,p_optvar_first,p_loc_control_offset,prev_proc);
        jCol_part[iOffset] = i_var_glob;//n_first_control + (cols[i]-control_first_);
        iOffset++;
      }
      MatRestoreRow(mat_aux_control,irow,&ncols,&cols,PETSC_NULL);
    }
    delete [] p_state_first;
    delete [] p_state_last;
    delete [] p_control_first;
    delete [] p_control_last;
    delete [] p_optvar_first;
    delete [] p_loc_control_offset;
  }
  else {
    // PDE constraints
    PetscInt irow;
    for (irow = pde_first_; irow<=pde_last_; ++irow) {
      PetscInt ncols;
      const PetscScalar* vals;
      // state
      MatGetRow(mat_state, irow, &ncols, PETSC_NULL, &vals);
      for (int i=0; i<ncols; ++i) {
        values_part[iOffset] = vals[i];
        iOffset++;
      }
      MatRestoreRow(mat_state,irow,&ncols,PETSC_NULL,&vals);
      // control
      MatGetRow(mat_control, irow, &ncols, PETSC_NULL, &vals);
      for (int i=0; i<ncols; ++i) {
        values_part[iOffset] = vals[i];
        iOffset++;
      }
      MatRestoreRow(mat_control,irow,&ncols,PETSC_NULL,&vals);
    }

    // auxiliary constraints
    for (irow=aux_first_; irow<=aux_last_; ++irow) {
      PetscInt ncols;
      const PetscScalar* vals;
      // state
      MatGetRow(mat_aux_state, irow, &ncols, PETSC_NULL, &vals);
      for (int i=0; i<ncols; ++i) {
        values_part[iOffset] = vals[i];
        iOffset++;
      }
      MatRestoreRow(mat_aux_state,irow,&ncols,PETSC_NULL,&vals);
      // control
      MatGetRow(mat_aux_control, irow, &ncols, PETSC_NULL, &vals);
      for (int i=0; i<ncols; ++i) {
        values_part[iOffset] = vals[i];
        iOffset++;
      }
      MatRestoreRow(mat_aux_control,irow,&ncols,PETSC_NULL,&vals);
    }
  }
  if (nele_jac_part != iOffset) {
    MY_DBG_PRINT( "n_ele_jac_part=" << nele_jac_part << "!=" << iOffset << "=iOffset" );
    assert(nele_jac_part == iOffset);
  }

  END_FUNCTION
  return true;
}

int PetscPDENLP::GlobStateControlIdx2OptVarIdx(int i_var, int num_of_procs,
    const int* proc_first_arr, const int* proc_last_arr,
    const int* proc_var_first, const int* proc_loc_offset, int & prev_proc)
{
  START_FUNCTION
  // start with previous proc, in hope to accelerate search, if distribution is contigious
  int rv;
  for (int i_proc=prev_proc;i_proc<num_of_procs; ++i_proc) {
    if ( (proc_first_arr[i_proc]<=i_var)
         && (proc_last_arr[i_proc]>=i_var) ) {
      prev_proc = i_proc;
      rv = proc_var_first[i_proc] + (i_var-proc_first_arr[i_proc]);
      if (proc_loc_offset)
        rv += proc_loc_offset[i_proc];
      END_FUNCTION
      return rv;
    }
  }

  for (int i_proc=0;i_proc<prev_proc; ++i_proc) {
    if ( (proc_first_arr[i_proc]<=i_var)
         && (proc_last_arr[i_proc]>=i_var) ) {
      prev_proc = i_proc;
      rv = proc_var_first[i_proc] + (i_var-proc_first_arr[i_proc]);
      if (proc_loc_offset)
        rv += proc_loc_offset[i_proc];
      END_FUNCTION
      return rv;
    }
  }

  int sz;
  MPI_Comm_size(MPI_COMM_WORLD,&sz);
  MY_DBG_PRINT( "PetscPDENLP::GlobStateControlIdx2OptVarIdx(" << i_var << ")=? (no proc found)" );
  std::cout << "proc_first_arr:";
  for (int i=0;i<sz;i++) std::cout << ", " << proc_first_arr[i];
  std::cout << std::endl;
  std::cout << "proc_last_arr:";
  for (int i=0;i<sz;i++) std::cout << ", " << proc_last_arr[i];
  std::cout << std::endl;

  assert(false);
  return -1;
}

bool
PetscPDENLP::eval_h(Index num_proc, Index proc_id,
                      Index n, Index n_first, Index n_last,
                      const Number* x, bool new_x, Number obj_factor,
                      Index m, Index m_first, Index m_last,
                      const Number* lambda,
                      bool new_lambda, Index nele_hess_part,
                      Index* iRow_part, Index* jCol_part,
                      Number* values_part)
{
  START_FUNCTION
  if (new_x) {
    update_x(x);
  }
  if (values_part==NULL) {
    PetscPDE_->m_CalcType = PetscPDEBase::StructureOnly;
    assert(iRow_part && jCol_part);
  }
  else {
    PetscPDE_->m_CalcType = PetscPDEBase::Values;
    assert(!iRow_part && !jCol_part);
  }

  Vec lambda_pde, lambda_aux;
  Vec tmp;
  int m_pde, m_aux;
  PetscPDE_->getPDEResVector(tmp);
  VecGetSize(tmp,&m_pde);
  PetscPDE_->getAuxConstrVector(tmp);
  VecGetSize(tmp,&m_aux);
  
  VecCreateSeq(MPI_COMM_WORLD,m_pde,&lambda_pde);
  VecCreateSeq(MPI_COMM_WORLD,m_aux,&lambda_aux);
  if (NULL==values_part) {
  VecSet(lambda_pde,1.0);
  VecSet(lambda_aux,1.0);
  }
  else {
    int *p_pde_first = new int[num_proc];
    int *p_pde_last = new int[num_proc];
    int *p_aux_first = new int[num_proc];
    int *p_aux_last = new int[num_proc];
    MPI_Allgather(&pde_first_,1,MPI_INT,p_pde_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&pde_last_ ,1,MPI_INT,p_pde_last ,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&aux_first_,1,MPI_INT,p_aux_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&aux_last_ ,1,MPI_INT,p_aux_last ,1,MPI_INT,MPI_COMM_WORLD);

    int iVal(0);
    int ni, *Idx;
    for (int i_proc=0;i_proc<num_proc;i_proc++) {
      ni = p_pde_last[i_proc]-p_pde_first[i_proc]+1;
      Idx = new int[ni];
      for (int i_constr=0;i_constr<ni;++i_constr) {
        Idx[i_constr] = iVal+i_constr;
      }
      VecSetValues(lambda_pde,ni,Idx,lambda+iVal,INSERT_VALUES);
      iVal += ni;
      delete [] Idx;
      
      ni = p_aux_last[i_proc]-p_aux_first[i_proc]+1;
      Idx = new int[ni];
      for (int i_constr=0;i_constr<ni;++i_constr) {
        Idx[i_constr] = iVal + i_constr;
      }
      VecSetValues(lambda_aux,ni,Idx,lambda+iVal,INSERT_VALUES);
      iVal += ni;
      delete [] Idx;
    }
    delete[] p_pde_first;
    delete[] p_pde_last;
    delete[] p_aux_first;
    delete[] p_aux_last;
  }

  Mat mat_hcc, mat_hcs, mat_hss;
  PetscPDE_->calc_hessians(obj_factor, lambda_pde, lambda_aux, mat_hcc, mat_hcs, mat_hss);
  PetscInt n_control, n_state;
  {
    PetscInt tmp;
    MatGetSize(mat_hcc,&n_control,&tmp);
    MatGetSize(mat_hss,&n_state,&tmp);
  }
  PetscInt ncols;
  const PetscInt *cols;
  const PetscScalar* vals;
  unsigned long iOffset(0);
  if (NULL==values_part) {
    int *p_state_first = new int[num_proc];
    int *p_state_last = new int[num_proc];
    int *p_control_first = new int[num_proc];
    int *p_control_last = new int[num_proc];
    int *p_optvar_first = new int[num_proc];
    int *p_loc_control_offset = new int[num_proc];
    MPI_Allgather(&state_first_,1,MPI_INT,p_state_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&state_last_ ,1,MPI_INT,p_state_last ,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&control_first_,1,MPI_INT,p_control_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&control_last_ ,1,MPI_INT,p_control_last ,1,MPI_INT,MPI_COMM_WORLD);
    int i_var(0);
    for (int i_proc=0;i_proc<num_proc;++i_proc) {
      p_optvar_first[i_proc] = i_var;
      i_var += (p_state_last[i_proc]-p_state_first[i_proc]+1) +
               (p_control_last[i_proc]-p_control_first[i_proc]+1);
      p_loc_control_offset[i_proc] = p_state_last[i_proc] - p_state_first[i_proc] + 1;
    }

    int prev_proc_row(0);
    int prev_proc_col(0);
    // state_state
    for (unsigned int i_row=0;i_row<(unsigned int)n_state;++i_row) {
      MatGetRow(mat_hss, i_row, &ncols, &cols, PETSC_NULL);
      unsigned int i_var_row = GlobStateControlIdx2OptVarIdx(i_row,num_proc,
                               p_state_first,p_state_last,p_optvar_first,NULL,prev_proc_row);
      for (unsigned int i_col=0; i_col<(unsigned int)ncols; ++i_col) {
        unsigned int i_var_col = GlobStateControlIdx2OptVarIdx(cols[i_col],num_proc,
                                 p_state_first,p_state_last,p_optvar_first,NULL,prev_proc_col);
        iRow_part[iOffset] = i_var_row;
        jCol_part[iOffset] = i_var_col;
        ++iOffset;
      }
      MatRestoreRow(mat_hss, i_row, &ncols, &cols, PETSC_NULL);
    }

    // control control
    prev_proc_row=0;
    prev_proc_col=0;
    for (unsigned int i_row=0;i_row<(unsigned int)n_control;++i_row) {
      MatGetRow(mat_hcc, i_row, &ncols, &cols, PETSC_NULL);
      unsigned int i_var_row = GlobStateControlIdx2OptVarIdx(i_row,num_proc,
                                 p_control_first,p_control_last,p_optvar_first,p_loc_control_offset,prev_proc_row);
      for (unsigned int i_col=0; i_col<(unsigned int)ncols; ++i_col) {
        unsigned int i_var_col = GlobStateControlIdx2OptVarIdx(cols[i_col],num_proc,p_control_first,p_control_last,
                                 p_optvar_first,p_loc_control_offset,prev_proc_col);
        iRow_part[iOffset] = i_var_row;
        jCol_part[iOffset] = i_var_col;
        ++iOffset;
      }
      MatRestoreRow(mat_hcc, i_row, &ncols, &cols, PETSC_NULL);
    }

    // control state
    prev_proc_row=0;
    prev_proc_col=0;
    for (unsigned int i_row=0;i_row<(unsigned int)n_control;++i_row) {
      MatGetRow(mat_hcs, i_row, &ncols, &cols, PETSC_NULL);
      unsigned int i_var_row = GlobStateControlIdx2OptVarIdx(i_row,num_proc,
                               p_control_first,p_control_last,p_optvar_first,p_loc_control_offset,prev_proc_row);
      for (unsigned int i_col=0; i_col<(unsigned int)ncols; ++i_col) {
        unsigned int i_var_col = GlobStateControlIdx2OptVarIdx(cols[i_col],num_proc,p_state_first,p_state_last,
                                 p_optvar_first,NULL,prev_proc_col);
        iRow_part[iOffset] = i_var_row;
        jCol_part[iOffset] = i_var_col;
        ++iOffset;
      }
      MatRestoreRow(mat_hcs, i_row, &ncols, &cols, PETSC_NULL);
    }

    assert(iOffset==nele_hess_part);

    delete[] p_state_first;
    delete[] p_state_last;
    delete[] p_control_first;
    delete[] p_control_last;
    delete[] p_optvar_first;
    delete[] p_loc_control_offset;
  }
  else {
    // state_state
    for (unsigned int i_row=0;i_row<(unsigned int)n_state;++i_row) {
      MatGetRow(mat_hss, i_row, &ncols, PETSC_NULL, &vals);
      for (unsigned int i_col=0; i_col<(unsigned int)ncols; ++i_col) {
        values_part[iOffset] = vals[i_col];
        ++iOffset;
      }
      MatRestoreRow(mat_hss, i_row, &ncols, PETSC_NULL, &vals);
    }

    // control control
    for (unsigned int i_row=0;i_row<(unsigned int)n_control;++i_row) {
      MatGetRow(mat_hcc, i_row, &ncols, PETSC_NULL, &vals);
      for (unsigned int i_col=0; i_col<(unsigned int)ncols; ++i_col) {
        values_part[iOffset] = vals[i_col];
        ++iOffset;
      }
      MatRestoreRow(mat_hcc, i_row, &ncols, PETSC_NULL, &vals);
    }

    // control state
    for (unsigned int i_row=0;i_row<(unsigned int)n_control;++i_row) {
      MatGetRow(mat_hcs, i_row, &ncols, PETSC_NULL, &vals);
      for (unsigned int i_col=0; i_col<(unsigned int)ncols; ++i_col) {
        values_part[iOffset] = vals[i_col];
        ++iOffset;
      }
      MatRestoreRow(mat_hcs, i_row, &ncols, PETSC_NULL, &vals);
    }
    assert(iOffset==nele_hess_part);
  }

  VecDestroy(lambda_pde);
  VecDestroy(lambda_aux);

  END_FUNCTION
  return true;
}

void PetscPDENLP::finalize_solution(SolverReturn status,
                                      Index n, const Number* x, const Number* z_L, const Number* z_U,
                                      Index m, const Number* g, const Number* lambda,
                                      Number obj_value,
                                      const IpoptData* ip_data,
                                      IpoptCalculatedQuantities* ip_cq)
{
  START_FUNCTION
  update_x(x);
  PetscPDE_->Write2File("Sol");

  Vec state, control, pde_residual, aux_residual;
  PetscPDE_->getStateVector(state);
  PetscPDE_->getControlVector(control);
  PetscPDE_->calcPDE_residual(pde_residual);
  PetscPDE_->calcPDE_residual(aux_residual);

  Vec state_lb_mults, state_ub_mults;
  Vec control_lb_mults, control_ub_mults;
  Vec pde_residual_mults, aux_constr_mults;

  VecDuplicate(state,&state_lb_mults);
  VecDuplicate(state,&state_ub_mults);
  VecDuplicate(control,&control_lb_mults);
  VecDuplicate(control,&control_ub_mults);
  VecDuplicate(pde_residual,&pde_residual_mults);
  VecDuplicate(aux_residual,&aux_constr_mults);

  optim_var_global2libMesh(z_L, state_lb_mults, control_lb_mults);
  optim_var_global2libMesh(z_U, state_ub_mults, control_ub_mults);

  {
    PetscScalar *Vals = NULL;
    VecGetArray(pde_residual_mults,&Vals);
    int iVal;
    int iLam = m_first_;
      for (iVal=0; iVal<=(pde_last_-pde_first_); ++iVal, ++iLam) {
      Vals[iVal] = lambda[iLam];
    }
    VecRestoreArray(pde_residual_mults,&Vals);
  }
  {
    PetscScalar *Vals = NULL;
    VecGetArray(aux_constr_mults,&Vals);
    int iVal;
    int iLam = m_first_ + (pde_last_-pde_first_)+1;
    for (iVal=0; iVal<=(aux_last_-aux_first_); ++iVal, ++iLam) {
      Vals[iVal] = lambda[iLam];
    }
    VecRestoreArray(aux_constr_mults,&Vals);
  }

  VecAssemblyBegin(state);
  VecAssemblyBegin(state);
  VecAssemblyBegin(control);
  VecAssemblyBegin(control);
  VecAssemblyBegin(pde_residual);
  VecAssemblyBegin(aux_residual);

  VecAssemblyEnd(state);
  VecAssemblyEnd(state);
  VecAssemblyEnd(control);
  VecAssemblyEnd(control);
  VecAssemblyEnd(pde_residual);
  VecAssemblyEnd(aux_residual);
  
  PetscPDE_->set_finalize_vectors(state_lb_mults, state_ub_mults,
				    control_lb_mults, control_ub_mults,
				    pde_residual_mults, aux_constr_mults);

  VecDestroy(state_lb_mults);
  VecDestroy(state_ub_mults);
  VecDestroy(control_lb_mults);
  VecDestroy(control_ub_mults);
  VecDestroy(pde_residual_mults);
  VecDestroy(aux_constr_mults);

  END_FUNCTION
}
