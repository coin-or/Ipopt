// Copyright (C) 2010 International Business Machines.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: IpPetscPDE.hpp 1938 2011-03-18 17:05:04Z huber $
//
// Authors:  Johannes Huber, Andreas Waechter     IBM    2010-09-03

#ifndef __IPLIBMESHPDE_HPP__
#define __IPLIBMESHPDE_HPP__

#include "IpReferenced.hpp"
#include "petscvec.h"
#include "petscmat.h"
#include <string>

//#define SIMULATION_CONVERGENCE

/** Base class for problem description of a PDE constrainted
    optimization problem with libMesh. */
class PetscPDEBase
{
public:
  enum CalculationModeType { StructureOnly, Values };
  CalculationModeType m_CalcType;

  /**@name Constructors/Destructors */
  //@{
  PetscPDEBase() : m_CalcType(Values) {}

  /** Default destructor */
  virtual ~PetscPDEBase() {}
  //@}

  /** receive libmesh vector that stores the current values of the
    *  controls. */
  virtual void getControlVector(Vec& Control)=0;

  /** receive libmesh vector that stores the current values of the
    *  state, i.e. which is saved as solution vector in the libmesh system. */
  virtual void getStateVector(Vec& State)=0;
  virtual void getAuxConstrVector(Vec& AuxConstr)=0;
  virtual void getPDEResVector(Vec& PdeResVec)=0;

  /** Update internal data after change of values of state and
  control variables. */
  virtual void update()=0;
  /** Calculate one part of the objective function. All parts of the
  processors are sumed up by Ipopt, so the intension is, to
  compute only the local contribution to the objective  */
  virtual void calc_objective_part(double& Val)=0;
  /** Calculate the local part of the objective gradient, i.e.the
      part that consists of the derivatives w.r.t. local optimization
      variables, uses lm_control_vec_ as control and lm_sys->solution
      as state */
  virtual void calc_objective_gradient(Vec& grad_state,
                                       Vec& grad_control)=0;
  /** Calculate the local part of the residual (e.g. for a linerar
      PDE typically written as A*x-b) of the constraining PDE,
      i.e.the part that consists of all derivatives of the local
      constraints, (i.e. the mesh), uses lm_control_vec_ as control
      and lm_sys->solution as state. Saves the residual value to
      residual_ and returns a handle to it */
  virtual void calcPDE_residual(Vec& residual)=0;
  /** Calculate the local part of the jacobians of the constraining
      PDE, i.e.the part that consists of the derivatives w.r.t to all
      control ans state variables of the local constraints, (i.e. the
      local part of the mesh), uses lm_control_vec_ as control and
      lm_sys_->solution as state. Saves the state jacobian value to
      lm_sys_->matrix and the control jacobian in jac_control_ and
      returns a handle to them */
  virtual void calcPDE_jacobians(Mat& jac_state,
                                 Mat& jac_control)=0;
  /** Calculate the local part of the hessian of the constraining PDE,
        i.e. the ???
        , uses lm_control_vec_ as control
        and lm_sys_->solution as state. Saves the intermediate results
        in hess_control_control_, hess_control_state_ and
        hess_state_state_ and returns a handle to them */
  virtual void calc_hessians(double sigma, Vec& lambda_loc_pde, Vec& lambda_loc_aux, Mat& Hcc, Mat& Hcs, Mat& Hss)=0;
  virtual void calcAux_constr(Vec& constr)=0;
  virtual void calcAux_jacobians(Mat& jac_state,
                                 Mat& jac_control)=0;
  virtual void Write2File( const std::string& pre_filename)=0;
  virtual void get_bounds(Vec& state_l,
                          Vec& state_u,
                          Vec& control_l,
                          Vec& control_u,
                          Vec& aux_constr_l,
                          Vec& aux_constr_u)=0;
  virtual void get_starting_point(Vec& state,
                                  Vec& control,
				                          bool init_z,
				                          Vec& state_lb_mults,
				                          Vec& state_ub_mults,
				                          Vec& control_lb_mults,
				                          Vec& control_ub_mults,
				                          bool init_lambda,
				                          Vec& pde_residual_mults,
				                          Vec& aux_constr_mults)=0;
  virtual void set_finalize_vectors(Vec& lm_state_lb_mults,
				                            Vec& lm_state_ub_mults,
				                            Vec& lm_control_lb_mults,
				                            Vec& lm_control_ub_mults,
				                            Vec& lm_pde_residual_mults,
				                            Vec& lm_aux_constr_mults)=0;
private:
  /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and 
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called. */
  //@{
  /** Default Constructor */
  //PetscPDEBase();
  /** Copy Constructor */
  PetscPDEBase(const PetscPDEBase&) {}
  /** Overloaded Equals Operator */
  void operator=(const PetscPDEBase&) {}
  //@}
};
#endif
