/* 
 * Copyright (C) 2007 VRTech Industrial Technologies - www.vrtech.com.br.
 * All Rights Reserved.
 * This code is published under the Common Public License.
 * 
 * $Id$
 * Authors: Rafael de Pelegrini Soares
 *
 *
 * Copyright (C) 2007 Tong Kewei, BeiHang University - www.buaa.edu.cn.
 * All Rights Reserved.
 * This code is published under the Common Public License.
 * 
 * $Id$
 * Authors: Tong Kewei, E-mail:tongkewei@126.com
 * Beihang University, website: www.buaa.edu.cn
 * Beijing,China.
 * 2007-11-11
 */

#include <jni.h>
#include <ipopt/IpTNLP.hpp>
#include "IpIpoptApplication.hpp"
#include "org_coinor_Ipopt.h"
#include "ipsmartptr.hpp"


//#include"string_igno_case.h"


using namespace std;
using namespace Ipopt;




/**
 * Main structure for Ipopt JNI implementation.
 * 
 * All functions will receive a pointer to this structure as
 * an integer argument (the address in memory of the structure).
 */

struct ipopt_jni {
	
	JNIEnv *env;
	jobject solver;
	
	jint n;
	jint m;
	jint nele_jac;
	jint nele_hess;

	jint index_style;

   // some cached arguments
	jdoubleArray mult_gj;
	
	// the callback arguments
	jdoubleArray xj;
	jdoubleArray fj;
	jdoubleArray grad_fj;
	jdoubleArray gj;
	jdoubleArray jac_gj;
	jdoubleArray hessj;

	jdoubleArray mult_x_Lj;
	jdoubleArray mult_x_Uj;

	jboolean using_scaling_parameters;
	jboolean using_LBFGS;


	SmartPtr<TNLP> problem;//IpoptProblem problem;
	
	SmartPtr<IpoptApplication> application; 

};



class Jipopt:public TNLP
{
public:
  /**  constructor */
  Jipopt(ipopt_jni* ipopt);

  /** default destructor */
  virtual ~Jipopt();

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
				 const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq);
  //@}



  /** overload this method to return scaling parameters. This is
     *  only called if the options are set to retrieve user scaling.
     *  There, use_x_scaling (or use_g_scaling) should get set to true
     *  only if the variables (or constraints) are to be scaled.  This
     *  method should return true only if the scaling parameters could
     *  be provided.
     */
    virtual bool get_scaling_parameters(Number& obj_scaling,
                                        bool& use_x_scaling, Index n,
                                        Number* x_scaling,
                                        bool& use_g_scaling, Index m,
                                        Number* g_scaling);
   


  /** @name Methods for quasi-Newton approximation.  If the second
     *  derivatives are approximated by Ipopt, it is better to do this
     *  only in the space of nonlinear variables.  The following
     *  methods are call by Ipopt if the quasi-Newton approximation is
     *  selected.  If -1 is returned as number of nonlinear variables,
     *  Ipopt assumes that all variables are nonlinear.  Otherwise, it
     *  calls get_list_of_nonlinear_variables with an array into which
     *  the indices of the nonlinear variables should be written - the
     *  array has the lengths num_nonlin_vars, which is identical with
     *  the return value of get_number_of_nonlinear_variables().  It
     *  is assumed that the indices are counted starting with 1 in the
     *  FORTRAN_STYLE, and 0 for the C_STYLE. */
    //@{
    virtual Index get_number_of_nonlinear_variables();
    

    virtual bool get_list_of_nonlinear_variables(Index num_nonlin_vars,Index* pos_nonlin_vars);

    ipopt_jni* ipopt;//

  	// the callback methods	
	//jmethodID get_nlp_info_;
	jmethodID get_bounds_info_;
	jmethodID get_starting_point_;
	jmethodID eval_f_;
	jmethodID eval_grad_f_;
	jmethodID eval_g_;
	jmethodID eval_jac_g_;
	jmethodID eval_h_;

	jmethodID get_scaling_parameters_;
	jmethodID get_number_of_nonlinear_variables_;
	jmethodID get_list_of_nonlinear_variables_;
	

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  //  MyNLP();
  Jipopt(const Jipopt&);
  Jipopt& operator=(const Jipopt&);
  
  //@}
};


//
Jipopt::Jipopt( ipopt_jni* ipopt){
	Jipopt::ipopt=ipopt;

	JNIEnv *env=ipopt->env;
	jobject solver=ipopt->solver;

	// the solver class
	jclass solverCls = env->GetObjectClass(solver);

    // get the methods
	get_bounds_info_=env->GetMethodID(solverCls,"get_bounds_info","(I[D[DI[D[D)Z");
	get_starting_point_=env->GetMethodID(solverCls,"get_starting_point","(IZ[DZ[D[DIZ[D)Z");
	eval_f_ = env->GetMethodID(solverCls, "eval_f", "(I[DZ[D)Z");
	eval_grad_f_ = env->GetMethodID(solverCls, "eval_grad_f", "(I[DZ[D)Z");
	eval_g_ = env->GetMethodID(solverCls, "eval_g", "(I[DZI[D)Z");
	eval_jac_g_ = env->GetMethodID(solverCls, "eval_jac_g", "(I[DZII[I[I[D)Z");
	eval_h_ = env->GetMethodID(solverCls, "eval_h", "(I[DZDI[DZI[I[I[D)Z");
	
	get_scaling_parameters_=env->GetMethodID(solverCls,"get_scaling_parameters","([DI[DI[D[Z)Z");
	get_number_of_nonlinear_variables_=env->GetMethodID(solverCls,"get_number_of_nonlinear_variables","()I");
	get_list_of_nonlinear_variables_=env->GetMethodID(solverCls,"get_list_of_nonlinear_variables","(I[I)Z");
	

	if(get_starting_point_==0||eval_f_==0 || eval_grad_f_==0 || eval_g_==0 || eval_jac_g_==0 ||    //
		eval_h_==0||get_number_of_nonlinear_variables_==0||get_list_of_nonlinear_variables_==0){
        std::cerr << "Expected callback methods missing on JIpopt.java" << std::endl;
		
		return;
	}

}

Jipopt::~Jipopt(){
}

bool Jipopt::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
						  Index& nnz_h_lag, IndexStyleEnum& index_style){

	n=Jipopt::ipopt->n;
	m=Jipopt::ipopt->m;
	nnz_jac_g=Jipopt::ipopt->nele_jac;
	nnz_h_lag=Jipopt::ipopt->nele_hess;
	
	index_style=(IndexStyleEnum)Jipopt::ipopt->index_style;
	  

	return true;

};


bool Jipopt::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u){

	
	JNIEnv *env=ipopt->env;
	jobject solver=ipopt->solver;
	jdoubleArray x_lj;
	jdoubleArray x_uj;
	jdoubleArray g_lj;
	jdoubleArray g_uj;

	if(x_l!=NULL||x_u!=NULL){
		x_lj=env->NewDoubleArray(n);
		x_uj=env->NewDoubleArray(n);
		
	}
	if(g_l!=NULL||g_u!=NULL){
		g_lj=env->NewDoubleArray(m);
		g_uj=env->NewDoubleArray(m);
		
	}

	if(!env->CallBooleanMethod(solver,get_bounds_info_,n,x_lj,x_uj,m,g_lj,g_uj)){
		
		return false;
	}
	
	// Copy from Java to native value 
	if(x_l!=NULL||x_u!=NULL){
	env->GetDoubleArrayRegion(x_lj, 0, n, x_l);//
	env->GetDoubleArrayRegion(x_uj, 0, n, x_u);//
	}
	if(g_l!=NULL||g_u!=NULL){
	env->GetDoubleArrayRegion(g_lj, 0, m, g_l);//
	env->GetDoubleArrayRegion(g_uj, 0, m, g_u);//
	}
		
	return true;	
	
}

bool Jipopt::get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
									Number* lambda){

										
	// cast back our structure
	JNIEnv *env=ipopt->env;
	jobject solver=ipopt->solver;

	jdoubleArray xj=NULL;
	jdoubleArray z_lj=NULL;
	jdoubleArray z_uj=NULL;
	jdoubleArray lambdaj=NULL;

	if(init_x){
		xj=Jipopt::ipopt->xj;
		
	}
	if(init_z){
		z_lj=ipopt->mult_x_Lj;
		z_uj=ipopt->mult_x_Uj;

	}
	if(init_lambda){
		lambdaj=Jipopt::ipopt->mult_gj;
	}

	if(!env->CallBooleanMethod(solver,get_starting_point_,
		n,init_x,xj,
		init_z,z_lj,z_uj,
		m,init_lambda,lambdaj)){			
		return false;
	}


	/* Copy from Java to native value */
	if(init_x){
		env->GetDoubleArrayRegion(xj, 0, n, x);//
		
	}
	
	if(init_z){
		
		env->GetDoubleArrayRegion(z_lj,0,n,z_L);
		env->GetDoubleArrayRegion(z_uj,0,n,z_U);
		
	}
	if(init_lambda){		
		env->GetDoubleArrayRegion(lambdaj,0,m,lambda);
	}
   
	return true;

}

bool Jipopt::eval_f(Index n, const Number* x, bool new_x,
                        Number& obj_value){

	
	

	// cast back our structure
	JNIEnv *env=ipopt->env;
	jobject solver=ipopt->solver;

	
	jdoubleArray xj = ipopt->xj;
	jdoubleArray fj = ipopt->fj;

	if(new_x){
	/* Copy the native double x to the Java double array xj */ 
	env->SetDoubleArrayRegion(xj, 0, n, x);
	}
    
	/* Call the java method */
	jboolean new_xj = new_x;
    if(!env->CallBooleanMethod(solver, eval_f_, n, xj, new_xj, fj))
    	return false;
    
	/* Copy from Java to native value */
    env->GetDoubleArrayRegion(fj, 0, 1, &obj_value);//should be a pointer
	
	
    return true;
}

bool Jipopt::eval_grad_f(Index n, const Number* x, bool new_x,
                             Number* grad_f){

	
	// cast back our structure
	JNIEnv *env=ipopt->env;
	jobject solver=ipopt->solver;
	

	jdoubleArray xj ;
	jdoubleArray grad_fj;
	
	
	xj=ipopt->xj;
	grad_fj=ipopt->grad_fj;

	if(new_x){
	/* Copy the native double x to the Java double array xj */ 
	env->SetDoubleArrayRegion(xj, 0, n, x);	
	}

	/* Call the java method */
	jboolean new_xj = new_x;
	if(!env->CallBooleanMethod(solver, eval_grad_f_, n, xj, new_xj, grad_fj)){    	
		return false;
	}	
	
	env->GetDoubleArrayRegion(grad_fj, 0, n, grad_f);
	
	return true;
}

bool Jipopt::eval_g(Index n, const Number* x, bool new_x,
                        Index m, Number* g){
	
	// cast back our structure
	JNIEnv *env=ipopt->env;
	jobject solver=ipopt->solver;

	
	jdoubleArray xj = ipopt->xj;
	jdoubleArray gj = ipopt->gj; 

	if(new_x){
		/* Copy the native double x to the Java double array xj */ 
		env->SetDoubleArrayRegion(xj, 0, n, x);
	}
	
	
	/* Call the java method */
	jboolean new_xj = new_x;
    jboolean ret = env->CallBooleanMethod(solver, eval_g_, n, xj, new_xj, m, gj);	
    if(!ret)
    	return false;
	/* Copy from Java to native value */
    env->GetDoubleArrayRegion(gj, 0, m, g);
	
	return true;
}

bool Jipopt::eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow,
                            Index *jCol, Number* jac_g){
	
	JNIEnv *env=ipopt->env;
	jobject solver=ipopt->solver;
 
	

 	jdoubleArray xj = NULL;
	jdoubleArray jac_gj = NULL;
	
	xj = ipopt->xj;
	jac_gj = ipopt->jac_gj;

  	
  	if( new_x &&x!=NULL){
		// Copy the native double x to the Java double array xj  
    	env->SetDoubleArrayRegion(xj, 0, n, x);
  	}
  	
  	/// Create the index arrays if needed
  	jintArray iRowj = NULL;
  	jintArray jColj = NULL;
  	if(iRow != NULL && jCol != NULL){
  		iRowj = env->NewIntArray(nele_jac);
  		jColj = env->NewIntArray(nele_jac);
  	}
	
	/* Call the java method */
	jboolean new_xj = new_x;
    if(!env->CallBooleanMethod(solver, eval_jac_g_, n, xj, new_xj,
		m, nele_jac, iRowj, jColj, jac_g == NULL ? NULL : jac_gj)){
    	
		return false;
	}
	
	/* Copy from Java to native value */
	if(jac_g != NULL){
    	env->GetDoubleArrayRegion(jac_gj, 0, nele_jac, jac_g);
		
	}
	if(iRow != NULL && jCol != NULL){    	
		jint *iRow_jint = env->GetIntArrayElements(iRowj, 0);
		jint *jCol_jint = env->GetIntArrayElements(jColj, 0);
		for(int i=0; i<nele_jac; ++i){
			iRow[i] = iRow_jint[i];
			jCol[i] = jCol_jint[i];
		}		
    }

	return true;

}

bool Jipopt::eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess,
                        Index* iRow, Index* jCol, Number* hess)
{
	
	// cast back our structure
	JNIEnv *env=ipopt->env;
	jobject solver=ipopt->solver;
 
	//NULL is necessary here.
	jdoubleArray xj = ipopt->xj;;
	jdoubleArray hessj = ipopt->hessj;
	jdoubleArray mult_gj = ipopt->mult_gj;;

  	

  	if(  new_x && x!=NULL){		
		/* Copy the native double x to the Java double array xj */ 
    	env->SetDoubleArrayRegion(xj, 0, n, x);
  	}
  	if( new_lambda && lambda!=NULL){		
		/* Copy the native double lambda to the Java double array lambdaj */ 
    	env->SetDoubleArrayRegion(mult_gj, 0, m, lambda);//multi_gj <==> lambdaj
  	}
  	
  	/// Create the index arrays if needed
  	jintArray iRowj = NULL;
  	jintArray jColj = NULL;
  	if(iRow != NULL && jCol != NULL){
  		iRowj = env->NewIntArray(nele_hess);
  		jColj = env->NewIntArray(nele_hess);
  	}


	/* Call the java method */
	jboolean new_xj = new_x;
	jboolean new_lambdaj = new_lambda;
    if(!env->CallBooleanMethod(solver, eval_h_, n, xj, new_xj,
    	obj_factor, m, mult_gj, new_lambdaj,
		nele_hess, iRowj, jColj, hess == NULL ? NULL : hessj)){			
    	return false;
	}
   
	/* Copy from Java to native value */
	if(hess != NULL){
    	env->GetDoubleArrayRegion(hessj, 0, nele_hess, hess);
		
	}
    if(iRow != NULL && jCol != NULL){
			
 		jint *iRow_jint = env->GetIntArrayElements(iRowj, 0);
		jint *jCol_jint = env->GetIntArrayElements(jColj, 0);
		
		
		for(int i=0; i<nele_hess; ++i){
			iRow[i] = iRow_jint[i];
			jCol[i] = jCol_jint[i];

			
		}
		
	}
	
	return true;
}

void Jipopt::finalize_solution(SolverReturn status, Index n, const Number *x, 
							   const Number *z_L, const Number *z_U, Index m, 
							   const Number *g, const Number *lambda, Number obj_value, 
							   const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq)
{
	/*
	//you can add anything you like, it is not need for me!

	JNIEnv *env=ipopt->env;
	jobject solver=ipopt->solver;

	
	*/
	

}


/** overload this method to return scaling parameters. This is
     *  only called if the options are set to retrieve user scaling.
     *  There, use_x_scaling (or use_g_scaling) should get set to true
     *  only if the variables (or constraints) are to be scaled.  This
     *  method should return true only if the scaling parameters could
     *  be provided.
     */
bool Jipopt::get_scaling_parameters(Number& obj_scaling,
                                        bool& use_x_scaling, Index n,
                                        Number* x_scaling,
                                        bool& use_g_scaling, Index m,
                                        Number* g_scaling)
    {
		

		if(Jipopt::ipopt->using_scaling_parameters){
			JNIEnv *env=ipopt->env;
			jobject solver=ipopt->solver;

			jdoubleArray obj_scaling_j=env->NewDoubleArray(1);
			jdoubleArray x_scaling_j=env->NewDoubleArray(n);
			jdoubleArray g_scaling_j=env->NewDoubleArray(m);

			jbooleanArray use_x_g_scaling_j=env->NewBooleanArray(2);			

			env->CallBooleanMethod(solver,get_scaling_parameters_,
				obj_scaling_j,
				n,x_scaling_j,
				m,g_scaling_j,
				use_x_g_scaling_j);

			jboolean* use_x_g_scaling=env->GetBooleanArrayElements(use_x_g_scaling_j,0);
			/* Copy from Java to native value */
			env->GetDoubleArrayRegion(obj_scaling_j, 0, 1, &obj_scaling);
			/* Copy from Java to native value */
			if(use_x_g_scaling[0]){			
				env->GetDoubleArrayRegion(x_scaling_j, 0, n, x_scaling);
				use_x_scaling=true;
			}else{
				use_x_scaling=false;
			}
			/* Copy from Java to native value */
			if(use_x_g_scaling[1]){			
				env->GetDoubleArrayRegion(g_scaling_j, 0, n, g_scaling);
				use_g_scaling=true;
			}else{
				use_g_scaling=false;
			}
			env->ReleaseBooleanArrayElements(use_x_g_scaling_j,use_x_g_scaling,0);
			return true;
		}else{
			return false;
		}
    }



Index Jipopt::get_number_of_nonlinear_variables(){
	if(Jipopt::ipopt->using_LBFGS){
		JNIEnv *env=ipopt->env;
		jobject solver=ipopt->solver;

		return env->CallIntMethod(solver,get_number_of_nonlinear_variables_);

	}else{
		return -1;
	}
}
    

bool Jipopt::get_list_of_nonlinear_variables(Index num_nonlin_vars,Index* pos_nonlin_vars){
	if(Jipopt::ipopt->using_LBFGS){
		JNIEnv *env=ipopt->env;
		jobject solver=ipopt->solver;

		jintArray pos_nonlin_vars_j=env->NewIntArray(num_nonlin_vars);
		
		jboolean ret=env->CallBooleanMethod(solver,get_list_of_nonlinear_variables_,
			num_nonlin_vars,pos_nonlin_vars_j);
		if(!ret){
			return false;
		}

		if(pos_nonlin_vars!=NULL){	
			//env->GetIntArrayRegion(pos_nonlin_vars_j,0,num_nonlin_vars,pos_nonlin_vars);
			jint *pos_nonlin_vars_jp = env->GetIntArrayElements(pos_nonlin_vars_j, 0);
			for(int i=0; i<num_nonlin_vars; ++i){
				pos_nonlin_vars[i]=pos_nonlin_vars_jp[i];
			}			
			//env->ReleaseIntArrayElements(pos_nonlin_vars_j,pos_nonlin_vars,0);
			
		}		
		return true;
	}else{
		return false;
	}

}



#ifdef __cplusplus
extern "C" {
#endif


JNIEXPORT jlong JNICALL Java_org_coinor_Ipopt_CreateIpoptProblem 
(JNIEnv *env, jobject obj_this, 
 jint n,  jint m,
 jint nele_jac, jint nele_hess,
 jint index_style)
{
	 ipopt_jni *ipopt=new ipopt_jni;

	ipopt->env=env;
	ipopt->solver=obj_this;	

	
	ipopt->n=n;
	ipopt->m=m;
	ipopt->nele_jac=nele_jac;
	ipopt->nele_hess=nele_hess;

	ipopt->index_style=index_style;
	
	ipopt->using_LBFGS=false;
	ipopt->using_scaling_parameters=false;

	/* create the IpoptProblem */
	Jipopt* problem=new Jipopt(ipopt);
	

	if(problem == NULL){
		
		delete problem;
		return 0;
	}
	
	ipopt->problem=problem;
	ipopt->application=new IpoptApplication();	
	

	if(ipopt == NULL){		
		delete ipopt;
		return 0;
	}
	   	
	return (jlong)ipopt;
}



JNIEXPORT jint JNICALL Java_org_coinor_Ipopt_OptimizeTNLP
  (JNIEnv *env, jobject obj_this, jlong pipopt, jstring outfilename ,
  jdoubleArray xj,
jdoubleArray gj,
jdoubleArray obj_valj,
jdoubleArray mult_gj,
jdoubleArray mult_x_Lj,
jdoubleArray mult_x_Uj,

jdoubleArray callback_grad_f,
jdoubleArray callback_jac_g,
jdoubleArray callback_hess)
{
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;
	// cast back our structure
	ipopt->env=env;
	ipopt->solver=obj_this;
	

	ipopt->xj = xj;
    ipopt->fj = obj_valj;
	ipopt->gj = gj;
    ipopt->grad_fj = callback_grad_f;
	ipopt->jac_gj = callback_jac_g;
	ipopt->hessj = callback_hess;
	ipopt->mult_gj = mult_gj;
	ipopt->mult_x_Lj=mult_x_Lj;
	ipopt->mult_x_Uj=mult_x_Uj;
	

	const char *pparameterName = env->GetStringUTFChars(outfilename, 0);
	string outfile=pparameterName;

	 //  (use a SmartPtr, not raw)	
	ApplicationReturnStatus status;
	
	if(outfilename){
		status=ipopt->application->Initialize(outfile);	
	}else{
		status=ipopt->application->Initialize();
	}	
	if(outfilename){
		env->ReleaseStringUTFChars(outfilename, pparameterName);
	}
	
	if (status != Solve_Succeeded) {
		printf("\n\n*** Error during initialization!\n");
		return (int) status;
	}

	
	

	/* solve the problem */
	status = ipopt->application->OptimizeTNLP(ipopt->problem);

  if (status == Solve_Succeeded) {
    printf("\n\n*** The problem solved!\n");
  }
  else {
    printf("\n\n*** The problem FAILED!\n");
  }
  

  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.

  return (int) status;
  
}



JNIEXPORT void JNICALL Java_org_coinor_Ipopt_FreeIpoptProblem
(JNIEnv *env, 
jobject obj_this, 
jlong pipopt){
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;	
	ipopt->env=env;
	ipopt->solver=obj_this;

	if(ipopt!=NULL){	
		delete ipopt;
	}
}



JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_AddIpoptIntOption
(JNIEnv * env, jobject obj_this, jlong pipopt, jstring jparname, jint jparvalue){
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;
	ipopt->env=env;
	ipopt->solver=obj_this;

	
	const char *pparameterName = env->GetStringUTFChars(jparname, 0);
	string parameterName=pparameterName;

	// Try to apply the integer option
	jboolean ret=ipopt->application->Options()->SetIntegerValue(parameterName,jparvalue);

	env->ReleaseStringUTFChars(jparname, pparameterName);
	
	return ret;
}

JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_AddIpoptNumOption
(JNIEnv * env, jobject obj_this, jlong pipopt, jstring jparname, jdouble jparvalue){  
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;
	ipopt->env=env;
	ipopt->solver=obj_this;
	
	const char *pparameterName = env->GetStringUTFChars(jparname, 0);
	string parameterName=pparameterName;

	// Try to set the real option
	jboolean ret = ipopt->application->Options()->SetNumericValue(parameterName,jparvalue);
	
	env->ReleaseStringUTFChars(jparname, pparameterName);
	
	return ret;
}

JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_AddIpoptStrOption
(JNIEnv * env, jobject obj_this, jlong pipopt, jstring jparname, jstring jparvalue){
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;
	ipopt->env=env;
	ipopt->solver=obj_this;
	
	const char *pparameterName = env->GetStringUTFChars(jparname, NULL);
	string parameterName=pparameterName;
	const char *pparameterValue = env->GetStringUTFChars(jparvalue, NULL);
	string parameterValue=pparameterValue;


	//parameterValue has been changed to LowerCase in Java!
	if(parameterName=="hessian_approximation"&&(parameterValue=="limited-memory")){
		ipopt->using_LBFGS=true;		
	}else if(parameterName=="nlp_scaling_method"&&(parameterValue=="user-scaling")){
		ipopt->using_scaling_parameters=true;		
	}else{		
	}

	// Try to apply the string option
	jboolean ret = ipopt->application->Options()->SetStringValue(parameterName,parameterValue);


	env->ReleaseStringUTFChars(jparname, pparameterName);
	env->ReleaseStringUTFChars(jparname, pparameterValue);
	
	return ret;
}

#ifdef __cplusplus
}
#endif
