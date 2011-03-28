// Copyright (C) 2010 University Basel
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Johannes Huber

#ifndef __IPLIBMESHPDE_TEMPRADIATION_HPP__
#define __IPLIBMESHPDE_TEMPRADIATION_HPP__

#include "IpPetscPDE.hpp"
#include <boundary_mesh.h>
#include "enum_order.h"
#include <dof_map.h>
#include <petsc_matrix.h>
#include <petsc_vector.h>

/** Base class for problem description of a PDE constrainted
    optimization problem with libMesh. */
class PetscPDETempRadiation : public PetscPDEBase
{
public:
  enum SolutionModeType { Simulate, Optimize };
  SolutionModeType m_SolutionMode;

  /**@name Constructors/Destructors */
  //@{
  PetscPDETempRadiation();

  /** Default destructor */
  virtual ~PetscPDETempRadiation();
  //@}

  /** receive Petsc vector that stores the current values of the
    *  controls. */
  virtual void getControlVector(Vec& rv)
  {
    rv = m_Control->vec();
  }
  /** receive Petsc vector that stores the current values of the
    *  state, i.e. which is saved as solution vector in the Petsc system. */
  virtual void getStateVector(Vec &rv)
  {
    rv = m_State->vec();
  }

  virtual void getAuxConstrVector(Vec& rv)
  {
    rv = m_AuxConstr->vec();
  }

  virtual void getPDEResVector(Vec& rv)
  {
    rv = m_PDEConstr->vec();
  }

  /** Update internal data after change of values of state and
  control variables. */
  virtual void update() {}
  /** Calculate one part of the objective function. All parts of the
  processors are sumed up by Ipopt, so the intension is, to
  compute only the local contribution to the objective  */
  virtual void calc_objective_part(double& Val);
  /** Calculate the local part of the objective gradient, i.e.the
      part that consists of the derivatives w.r.t. local optimization
      variables, uses lm_control_vec_ as control and lm_sys->solution
      as state */
  virtual void calc_objective_gradient(Vec& grad_state,
                                       Vec& grad_control);
  /** Calculate the local part of the residual (e.g. for a linerar
      PDE typically written as A*x-b) of the constraining PDE,
      i.e.the part that consists of all derivatives of the local
      constraints, (i.e. the mesh), uses lm_control_vec_ as control
      and lm_sys->solution as state. Saves the residual value to
      residual_ and returns a handle to it */
  virtual void calcPDE_residual(Vec& residual);
  /** Calculate the local part of the jacobians of the constraining
      PDE, i.e.the part that consists of the derivatives w.r.t to all
      control ans state variables of the local constraints, (i.e. the
      local part of the mesh), uses lm_control_vec_ as control and
      lm_sys_->solution as state. Saves the state jacobian value to
      lm_sys_->matrix and the control jacobian in jac_control_ and
      returns a handle to them */
  void calcPDE_jacobians(Mat& jac_state, Mat& jac_control);
  /** Calculate the local part of the hessian of the constraining PDE,
        i.e. the ???
        , uses lm_control_vec_ as control
        and lm_sys_->solution as state. Saves the intermediate results
        in hess_control_control_, hess_control_state_ and
        hess_state_state_ and returns a handle to them */
  virtual void calc_hessians(double sigma, Vec& lambda_loc_pde, Vec& lambda_loc_aux, Mat& Hcc, Mat& Hcs, Mat& Hss);
  virtual void calcAux_constr(Vec& constr);
  virtual void calcAux_jacobians(Mat& jac_state,
                                 Mat& jac_control);
  virtual void get_bounds(Vec& state_l,
                          Vec& state_u,
                          Vec& control_l,
                          Vec& control_u,
                          Vec& aux_constr_l,
                          Vec& aux_constr_u);
  virtual void get_starting_point(Vec& state,
                                  Vec& control,
				                          bool init_z,
				                          Vec& state_lb_mults,
				                          Vec& state_ub_mults,
				                          Vec& control_lb_mults,
				                          Vec& control_ub_mults,
				                          bool init_lambda,
				                          Vec& pde_residual_mults,
				                          Vec& aux_constr_mults);
  virtual void set_finalize_vectors(Vec& lm_state_lb_mults,
				                            Vec& lm_state_ub_mults,
				                            Vec& lm_control_lb_mults,
				                            Vec& lm_control_ub_mults,
				                            Vec& lm_pde_residual_mults,
				                            Vec& lm_aux_constr_mults);
  void Init(const std::string filename);
  virtual void Write2File( const std::string& pre_filename);
  double GetOuterMaxTemp() {return m_OuterMaxTemp;}
protected:
  AutoPtr<libMesh::PetscVector<libMesh::Number> > m_Control;
  AutoPtr<libMesh::PetscVector<libMesh::Number> > m_State;
  AutoPtr<libMesh::PetscVector<libMesh::Number> > m_PDEConstr;
  AutoPtr<libMesh::PetscVector<libMesh::Number> > m_AuxConstr;
  AutoPtr<libMesh::PetscVector<libMesh::Number> > m_GradObjControl;
  AutoPtr<libMesh::PetscVector<libMesh::Number> > m_GradObjState;
  AutoPtr<libMesh::PetscMatrix<libMesh::Number> > m_JacPDEState;
  AutoPtr<libMesh::PetscMatrix<libMesh::Number> > m_JacPDEControl;
  AutoPtr<libMesh::PetscMatrix<libMesh::Number> > m_JacAuxState;
  AutoPtr<libMesh::PetscMatrix<libMesh::Number> > m_JacAuxControl;
  AutoPtr<libMesh::PetscMatrix<libMesh::Number> > m_HessControlControl;
  AutoPtr<libMesh::PetscMatrix<libMesh::Number> > m_HessControlState;
  AutoPtr<libMesh::PetscMatrix<libMesh::Number> > m_HessStateState;
  int m_Dim;
  libMesh::Mesh m_StateMesh;
  //libMesh::BoundaryMesh m_ControlMesh;
  bool IsInOmegaInner(const libMesh::Point& pt);
  bool IsInOmegaOuter(const libMesh::Point& pt);
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
  PetscPDETempRadiation(const PetscPDETempRadiation& org) {}

  /** Overloaded Equals Operator */
  void operator=(const PetscPDETempRadiation&) {}
  //@}
  void clear_math_obj();
  void DetroySelfOwnedLibMeshPetscMatrix(AutoPtr<PetscMatrix<Number> >& matrix);
  void DetroySelfOwnedLibMeshPetscVector(AutoPtr<PetscVector<Number> >& vector);
  
  libMeshEnums::Order m_FEOrder;
  AutoPtr<DofMap> m_StateDofMap;
  //AutoPtr<DofMap> m_ControlDofMap;
  std::map<unsigned int, unsigned int> m_ControlNodeIDToControlDOF;
  double m_PDEConstrScale;
  double m_Alpha;           // Boundary Condition parameter
  double m_Beta;            // Tikhonov reularization parameter
  double m_OmegaInner[6];   // border of hot inner region [x0Min, x1Min, x2Min, x0Max, x1Max, x2Max]
  double m_OmegaOuter[6];   // border of cold outer region R^n/[x0Min, x1Min, x2Min, x0Max, x1Max, x2Max]
  double m_OuterMaxTemp;    // upper bound of state variable in OmegaOuter
};
#endif
