// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-16

#include "SensAlgorithm.hpp"
#include "SensUtils.hpp"
#include "IpSmartPtr.hpp"

#include "IpVector.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  SensAlgorithm::SensAlgorithm(std::vector< SmartPtr<SchurDriver> >& driver_vec,
			       SmartPtr<SensitivityStepCalculator> sens_step_calc,
			       SmartPtr<Measurement> measurement,
			       Index n_sens_steps)
    :
    driver_vec_(driver_vec),
    sens_step_calc_(sens_step_calc),
    measurement_(measurement),
    n_sens_steps_(n_sens_steps), // why doesn't he get this from the options?
    Sensitivity_X_(NULL),
    Sensitivity_L_(NULL),
    Sensitivity_Z_L_(NULL),
    Sensitivity_Z_U_(NULL)
  {
    DBG_START_METH("SensAlgorithm::SensAlgorithm", dbg_verbosity);
    DBG_ASSERT(n_sens_steps<=driver_vec.size());
  }

  SensAlgorithm::~SensAlgorithm()
  {
    DBG_START_METH("SensAlgorithm::~SensAlgorithm", dbg_verbosity);
    if (NULL != Sensitivity_X_) delete [] Sensitivity_X_ ;
    if (NULL != Sensitivity_L_) delete [] Sensitivity_L_ ;
    if (NULL != Sensitivity_Z_U_) delete [] Sensitivity_Z_U_ ;
    if (NULL != Sensitivity_Z_L_) delete [] Sensitivity_Z_L_ ;
  }

  bool SensAlgorithm::InitializeImpl(const OptionsList& options,
				     const std::string& prefix)
  {
    // initialize values for variable sizes, and allocate memory for sensitivity vectors
    nx_ = dynamic_cast<const DenseVector*>( GetRawPtr( IpData().curr()->x() ) )->Dim() ;
    nceq_ = dynamic_cast<const DenseVector*>( GetRawPtr( IpData().curr()->y_c() ) )->Dim() ;
    ncineq_ = dynamic_cast<const DenseVector*>( GetRawPtr( IpData().curr()->y_d() ) )->Dim() ;
    nzl_ = dynamic_cast<const DenseVector*>( GetRawPtr( IpData().curr()->z_L() ) )->Dim() ;
    nzu_ = dynamic_cast<const DenseVector*>( GetRawPtr( IpData().curr()->z_U() ) )->Dim() ;    
    nl_ = nceq_ + ncineq_ ;
    
    Sensitivity_X_ = new Number[nx_] ;
    if (NULL == Sensitivity_X_) return false ;
    Sensitivity_L_ = new Number[nceq_+ncineq_] ;
    if (NULL == Sensitivity_L_) return false ;
    Sensitivity_Z_L_ = new Number[nzl_] ;
    if (NULL == Sensitivity_Z_L_) return false ;
    Sensitivity_Z_U_ = new Number[nzu_] ;
    if (NULL == Sensitivity_Z_U_) return false ;
    return true;
  }

  /** Main loop: Wait for new measurement, Get new step, maybe deal with
   *  bounds,  see to it that everything happens in the required
   *  timeframe. */
  SensAlgorithmExitStatus SensAlgorithm::Run()
  {
    DBG_START_METH("SensAlgorithm::Run", dbg_verbosity);

    SensAlgorithmExitStatus retval = SOLVE_SUCCESS;

    /* Loop through all steps */
    SmartPtr<IteratesVector> sol = IpData().curr()->MakeNewIteratesVector();
    SmartPtr<DenseVector> delta_u;
    SmartPtr<const Vector> unscaled_x;
    SmartPtr<const Vector> unscaled_yc;
    
    SmartPtr<IteratesVector> trialcopy;
    for (Index step_i=0; step_i<n_sens_steps_; ++step_i) {
      sens_step_calc_->SetSchurDriver(driver_vec_[step_i]);
      delta_u = measurement_->GetMeasurement(step_i+1);
      delta_u->Print(Jnlst(),J_VECTOR,J_USER1,"delta_u");
      sens_step_calc_->Step(*delta_u, *sol);
      SmartPtr<IteratesVector> saved_sol = sol->MakeNewIteratesVectorCopy();
      saved_sol->Print(Jnlst(),J_VECTOR,J_USER1,"sol_vec");

      // unscale solution...
      UnScaleIteratesVector(&saved_sol) ;

      // update variables
      measurement_->SetSolution(step_i+1, saved_sol);

      // get sensitivity vector
      GetSensitivities() ;
      
    }

    return retval;
  }

  void SensAlgorithm::GetSensitivities(void) {

    /*
      Extract sensitivity vector for each vector type
      primal, lagrange, and bound multipliers(zl,zu)
    */
    SmartPtr<IteratesVector> SV = sens_step_calc_->GetSensitivityVector() ;    
    UnScaleIteratesVector(&SV) ;


    const Number* X_ = dynamic_cast<const DenseVector*>( GetRawPtr( (*SV).x() ) )->Values();

    for (int i = 0; i < nx_; ++i) { 
      //printf(" ds/dp(X)[%3d] = %.14g\n", i+1, X_[i]);
      Sensitivity_X_[i] = X_[i] ;
    }

    const Number* Z_L_ = dynamic_cast<const DenseVector*>( GetRawPtr( (*SV).z_L() ) )->Values();
    for (int i = 0; i < nzl_; ++i) { 
      //printf(" ds/dp(X)[%3d] = %.14g\n", i+1, X_[i]);
      Sensitivity_Z_L_[i] = Z_L_[i] ;
    }

    const Number* Z_U_ = dynamic_cast<const DenseVector*>( GetRawPtr( (*SV).z_U() ) )->Values();
    for (int i = 0; i < nzu_; ++i) { 
      //printf(" ds/dp(X)[%3d] = %.14g\n", i+1, X_[i]);
      Sensitivity_Z_U_[i] = Z_U_[i] ;
    }

    const Number* LE_  = dynamic_cast<const DenseVector*>( GetRawPtr( (*SV).y_c() ) )->Values();
    for (int i = 0; i < nceq_; ++i) {
      //printf(" ds/dp(LE)[%3d] = %.14g\n", i+1, LE_[i]);
      Sensitivity_L_[i] = LE_[i] ;
    }

    const Number* LIE_ = dynamic_cast<const DenseVector*>( GetRawPtr( (*SV).y_d() ) )->Values();
      for (int i = 0; i < ncineq_; ++i) {
	//printf(" ds/dp(LIE)[%3d] = %.14g\n", i+1, LIE_[i]);
	Sensitivity_L_[i+nceq_] = LIE_[i] ;
      }

  }

  void SensAlgorithm::UnScaleIteratesVector(SmartPtr<IteratesVector> *V) {

    // unscale the iterates vector
    // pretty much a copy from IpOrigIpopt::finalize_solution
    
    SmartPtr<const Vector> unscaled_x;
    unscaled_x = IpNLP().NLP_scaling()->unapply_vector_scaling_x((*V)->x());
    DBG_ASSERT(IsValid(unscaled_x));
    (*V)->Set_x(*unscaled_x);
    unscaled_x = NULL ;

    SmartPtr<const Matrix> Px_L = IpNLP().Px_L();
    SmartPtr<const Matrix> Px_U = IpNLP().Px_U();
    SmartPtr<const VectorSpace> x_space = IpNLP().x_space();

    SmartPtr<const Vector> y_c = (*V)->y_c();
    SmartPtr<const Vector> y_d = (*V)->y_d();
    
    SmartPtr<const Vector> z_L = (*V)->z_L();
    SmartPtr<const Vector> z_U = (*V)->z_U();

    
    // unscale y_c
    SmartPtr<const Vector> unscaled_yc;
    SmartPtr<const Vector> unscaled_yd;
    SmartPtr<const Vector> unscaled_z_L;
    SmartPtr<const Vector> unscaled_z_U;


    Number obj_unscale_factor = IpNLP().NLP_scaling()->unapply_obj_scaling(1.);
    if (obj_unscale_factor!=1.) {
      
      SmartPtr<Vector> tmp = IpNLP().NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*Px_L, z_L, *x_space);
      tmp->Scal(obj_unscale_factor);
      unscaled_z_L = ConstPtr(tmp);    
      
      tmp = IpNLP().NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*Px_U, z_U, *x_space);
      tmp->Scal(obj_unscale_factor);
      unscaled_z_U = ConstPtr(tmp);

      tmp = IpNLP().NLP_scaling()->apply_vector_scaling_c_NonConst(y_c);
      tmp->Scal(obj_unscale_factor);
      unscaled_yc = ConstPtr(tmp);
      
      tmp = IpNLP().NLP_scaling()->apply_vector_scaling_d_NonConst(y_d);
      tmp->Scal(obj_unscale_factor);
      unscaled_yd = ConstPtr(tmp);
      
    }
    else {
      
      unscaled_z_L = IpNLP().NLP_scaling()->apply_vector_scaling_x_LU(*Px_L, z_L, *x_space);
      unscaled_z_U = IpNLP().NLP_scaling()->apply_vector_scaling_x_LU(*Px_U, z_U, *x_space);
      unscaled_yc = IpNLP().NLP_scaling()->apply_vector_scaling_c(y_c);
      unscaled_yd = IpNLP().NLP_scaling()->apply_vector_scaling_d(y_d);
      
    }

    (*V)->Set_z_U(*unscaled_z_U);
    (*V)->Set_z_L(*unscaled_z_L);
    (*V)->Set_y_c(*unscaled_yc);
    (*V)->Set_y_d(*unscaled_yd);
    
  }

}
