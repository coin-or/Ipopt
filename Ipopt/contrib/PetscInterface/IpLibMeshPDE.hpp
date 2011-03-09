// Copyright (C) 2010 International Business Machines.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Johannes Huber, Andreas Waechter     IBM    2010-09-03

#ifndef __IPLIBMESHPDE_HPP__
#define __IPLIBMESHPDE_HPP__

#include "IpReferenced.hpp"
#include "libmesh.h"
#include "implicit_system.h"
#include "numeric_vector.h"

#include "ProblemGeometry.h"

#include "petsc_vector.h"
#include "equation_systems.h"
#include "implicit_system.h"
#include "dense_vector.h"
#include "enum_order.h"
#include <string>

//#define SIMULATION_CONVERGENCE

/** Base class for problem description of a PDE constrainted
    optimization problem with libMesh. */
class LibMeshPDEBase
{
public:
  /**@name Constructors/Destructors */
  //@{
  LibMeshPDEBase();

  /** Default destructor */
  virtual ~LibMeshPDEBase();
  //@}

  /** receive libmesh EquationSystem.  It is assumed that this
    *  system has been initialized with the mesh, and the system
    *  matrix and solution vector are initialized. */
  libMesh::ImplicitSystem& getImplicitSystem()
  {
    return lm_eqn_sys_->get_system<ImplicitSystem>("PDE");
  }

  /** receive libmesh vector that stores the current values of the
    *  controls. */
  libMesh::NumericVector<libMesh::Number>& getControlVector()
  {
    return *lm_control_vec_;
  }
  /** receive libmesh vector that stores the current values of the
    *  state, i.e. which is saved as solution vector in the libmesh system. */
  libMesh::NumericVector<libMesh::Number>& getStateVector()
  {
    return *(getImplicitSystem().solution);
  }

  libMesh::NumericVector<libMesh::Number>& getAuxConstrVector()
  {
    return *lm_aux_constr_vec_;
  }

  libMesh::NumericVector<libMesh::Number>& getPDEResVector()
  {
    return *lm_pde_residual_vec_;
  }

  /** Update internal data after change of values of state and
  control variables. */
  void update()
  {
    ConvertControl2PGData();
    getImplicitSystem().update();
  }

  /** Calculate one part of the objective function. All parts of the
  processors are sumed up by Ipopt, so the intension is, to
  compute only the local contribution to the objective  */
  virtual void calc_objective_part(Number& Val);

  /** Calculate the local part of the objective gradient, i.e.the
      part that consists of the derivatives w.r.t. local optimization
      variables, uses lm_control_vec_ as control and lm_sys->solution
      as state */
  virtual void calc_objective_gradient(libMesh::NumericVector<libMesh::Number>& grad_state,
                                       libMesh::NumericVector<libMesh::Number>& grad_control);

  /** Calculate the local part of the residual (e.g. for a linerar
      PDE typically written as A*x-b) of the constraining PDE,
      i.e.the part that consists of all derivatives of the local
      constraints, (i.e. the mesh), uses lm_control_vec_ as control
      and lm_sys->solution as state. Saves the residual value to
      residual_ and returns a handle to it */
  virtual void calcPDE_residual(libMesh::NumericVector<libMesh::Number>*& residual);

  /** Calculate the local part of the jacobians of the constraining
      PDE, i.e.the part that consists of the derivatives w.r.t to all
      control ans state variables of the local constraints, (i.e. the
      local part of the mesh), uses lm_control_vec_ as control and
      lm_sys_->solution as state. Saves the state jacobian value to
      lm_sys_->matrix and the control jacobian in jac_control_ and
      returns a handle to them */
  virtual void calcPDE_jacobians(libMesh::SparseMatrix<libMesh::Number>*& jac_state,
                                 libMesh::SparseMatrix<libMesh::Number>*& jac_control);

  /** Calculate the local part of the hessian of the constraining PDE,
        i.e. the ???
        , uses lm_control_vec_ as control
        and lm_sys_->solution as state. Saves the intermediate results
        in hess_control_control_, hess_control_state_ and
        hess_state_state_ and returns a handle to them */
  virtual void calc_hessians(Number sigma, libMesh::DenseVector<Number>& lambda_loc_pde, libMesh::DenseVector<Number>& lambda_loc_aux, libMesh::SparseMatrix<Number>*& Hcc, libMesh::SparseMatrix<Number>*& Hcs,libMesh::SparseMatrix<Number>*& Hss);

  virtual void calcAux_constr(libMesh::NumericVector<libMesh::Number>*& constr);
  virtual void calcAux_jacobians(libMesh::SparseMatrix<libMesh::Number>*& jac_state,
                                 libMesh::SparseMatrix<libMesh::Number>*& jac_control);

  virtual void Write2File( const std::string& pre_filename);

  enum CalculationModeType
  {
    StructureOnly, Values
  } ;
  CalculationModeType calc_type_;

  virtual void get_bounds(libMesh::NumericVector<libMesh::Number>& state_l,
                          libMesh::NumericVector<libMesh::Number>& state_u,
                          libMesh::NumericVector<libMesh::Number>& control_l,
                          libMesh::NumericVector<libMesh::Number>& control_u,
                          libMesh::NumericVector<libMesh::Number>& aux_constr_l,
                          libMesh::NumericVector<libMesh::Number>& aux_constr_u);

  virtual void get_starting_point(libMesh::NumericVector<libMesh::Number>& state,
                                  libMesh::NumericVector<libMesh::Number>& control,
				  bool init_z,
				  libMesh::NumericVector<libMesh::Number>* state_lb_mults,
				  libMesh::NumericVector<libMesh::Number>* state_ub_mults,
				  libMesh::NumericVector<libMesh::Number>* control_lb_mults,
				  libMesh::NumericVector<libMesh::Number>* control_ub_mults,
				  bool init_lambda,
				  libMesh::NumericVector<libMesh::Number>* pde_residual_mults,
				  libMesh::NumericVector<libMesh::Number>* aux_constr_mults);

  virtual void InitProblemData(std::istream& is);
  virtual void reinit();

  virtual void get_finalize_vectors(libMesh::NumericVector<libMesh::Number>*& lm_state_lb_mults,
				    libMesh::NumericVector<libMesh::Number>*& lm_state_ub_mults,
				    libMesh::NumericVector<libMesh::Number>*& lm_control_lb_mults,
				    libMesh::NumericVector<libMesh::Number>*& lm_control_ub_mults,
				    libMesh::NumericVector<libMesh::Number>*& lm_pde_residual_mults,
				    libMesh::NumericVector<libMesh::Number>*& lm_aux_constr_mults);

  void RefineMesh(int iter);
  bool simulation_mode_;
  double CalcL2Diff(libMesh::NumericVector<libMesh::Number>* jac_state);
  void GetLocalIneqIdx(int* low, int* high);
protected:
  libMesh::SparseMatrix<libMesh::Number>* jac_control_;
  libMesh::SparseMatrix<libMesh::Number>* jac_aux_state_;
  libMesh::SparseMatrix<libMesh::Number>* jac_aux_control_;
  libMesh::SparseMatrix<libMesh::Number>* hess_control_control_;
  libMesh::SparseMatrix<libMesh::Number>* hess_control_state_;
  libMesh::SparseMatrix<libMesh::Number>* hess_state_state_;

  virtual void InitAuxConstr(int *plocal, int *pglobal, std::list<Number>* pFactList);

  virtual void calcPDE_jacobian_state(libMesh::SparseMatrix<libMesh::Number>*& jac_state);
  virtual void calcPDE_jacobian_control(libMesh::SparseMatrix<libMesh::Number>*& jac_control);
  virtual void calcAux_jacobian_state(libMesh::SparseMatrix<libMesh::Number>*& jac_state);
  virtual void calcAux_jacobian_control(libMesh::SparseMatrix<libMesh::Number>*& jac_control);

  static void assemble_Phi_PDE(EquationSystems& es, const std::string& system_name);

  void WriteNodes(const MeshBase& mesh, std::ostream& os);
  void WriteNodeFile(const MeshBase& mesh, const std::string& Filename);
  void WriteElems(const MeshBase& mesh, std::ostream& os);
  void WriteEleFile(const MeshBase& mesh, const std::string& Filename);
  void WritePotentialCSV(const std::string& Filename);
  void WriteAirflowCSVs(const std::string& VolumeFilename, const std::string& SurfFilename);
  void WriteAirflowTKVs(const std::string& VolumeFilename, const std::string& SurfFilename);

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
  //LibMeshPDEBase();

  /** Copy Constructor */
  LibMeshPDEBase(const LibMeshPDEBase&);

  /** Overloaded Equals Operator */
  void operator=(const LibMeshPDEBase&);
  //@}

  int GetPinNodeDof(bool* bLocal=NULL);
  int GetPinConstrIdx();
  void GetControlIdx(int* low, int* high);
  int GetMassConservationConstrIdx();
  /** Store the libMesh system */
  libMesh::EquationSystems* lm_eqn_sys_;

  /** Current values of control vector */
  libMesh::NumericVector<libMesh::Number>* lm_control_vec_;
  libMesh::NumericVector<libMesh::Number>* lm_pde_residual_vec_;
  libMesh::NumericVector<libMesh::Number>* lm_aux_constr_vec_;
  libMesh::NumericVector<libMesh::Number>* lm_aux_constr_vec_low_bd_;

  /** Multipliers at solution */
  AutoPtr<libMesh::NumericVector<libMesh::Number> > lm_control_lb_mults_;
  AutoPtr<libMesh::NumericVector<libMesh::Number> > lm_control_ub_mults_;
  AutoPtr<libMesh::NumericVector<libMesh::Number> > lm_aux_constr_mults_;

  libMesh::Number min_airflow;
  unsigned int first_aux_constr_;
  ProblemGeometry PG_;
  libMesh::Mesh mesh_;
  std::set<int> AuxConstrBoundMarkerList_;  // vector of all boundary markers, for which inequality constraints hold

  void ConvertControl2PGData();
  void DetroySelfOwnedLibMeshPetscMatrix(SparseMatrix<Number>*& matrix);
  void DetroySelfOwnedLibMeshPetscVector(NumericVector<Number>*& vector);
  void clear_math_obj();
  const libMeshEnums::Order lm_Num_quadrature_order_;
  int pin_down_node_;
  int pin_down_constr_;
  int mass_conservation_constr_;
};
#endif
