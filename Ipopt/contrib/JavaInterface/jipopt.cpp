/* 
 * Copyright (C) 2007 VRTech Industrial Technologies - www.vrtech.com.br.
 * All Rights Reserved.
 * This code is published under the Common Public License.
 * 
 * $Id$
 * Authors: Rafael de Pelegrini Soares, Edson C. do Valle
 */

#include <jni.h>
#include <ipopt/IpStdCInterface.h>
#include "org_coinor_Ipopt.h"

#include <iostream>

/**
 * Main structure for Ipopt JNI implementation.
 * 
 * All functions will receive a pointer to this structure as
 * an integer argument (the address in memory of the structure).
 */
struct ipopt_jni {
	JNIEnv *env;
	
	jobject solver;

	// the callback methods	
	jmethodID eval_f;
	jmethodID eval_grad_f;
	jmethodID eval_g;
	jmethodID eval_jac_g;
	jmethodID eval_h;

	// some cached arguments
	jdoubleArray mult_gj;
	
	// the callback arguments
	jdoubleArray xj;
	jdoubleArray fj;
	jdoubleArray grad_fj;
	jdoubleArray gj;
	jdoubleArray jac_gj;
	jdoubleArray hessj;

	// equation and nnz sizes
	jint n;
	jint m;
	jint nele_jac;
	jint nele_hess;
	
	// ipopt problem
	IpoptProblem problem;
};

#ifdef __cplusplus
extern "C" {
#endif

/* Ipopt Functions (these functions will call back to java) */
Bool eval_f(Index n, Number* x, Bool new_x,
	Number* obj_value, UserDataPtr user_data);
Bool eval_grad_f(Index n, Number* x, Bool new_x,
		Number* grad_f, UserDataPtr user_data);
Bool eval_g(Index n, Number* x, Bool new_x,
	Index m, Number* g, UserDataPtr user_data);
Bool eval_jac_g(Index n, Number *x, Bool new_x,
		Index m, Index nele_jac,
		Index *iRow, Index *jCol, Number *jac_g,
		UserDataPtr user_data);
Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
	Index m, Number *lambda, Bool new_lambda,
	Index nele_hess, Index *iRow, Index *jCol,
	Number *values, UserDataPtr user_data);

JNIEXPORT jlong JNICALL Java_org_coinor_Ipopt_CreateIpoptProblem 
(JNIEnv *env, 
jobject obj_this, 
jint n,
jdoubleArray ylb, 
jdoubleArray yub,
jint m,
jdoubleArray glb, 
jdoubleArray gub,
jint nele_jac,
jint nele_hess,
jint index_style){

    // create the ipopt_jni structure
	ipopt_jni *ipopt = new ipopt_jni;
	ipopt->env = env;
	ipopt->solver = obj_this;
	ipopt->n = n;
	ipopt->m = m;
	ipopt->nele_jac = nele_jac;
	ipopt->nele_hess = nele_hess;
	
	// the solver class
	jclass solverCls = env->GetObjectClass(ipopt->solver);
    // get the methods
	ipopt->eval_f = env->GetMethodID(solverCls, "eval_f", "(I[DZ[D)Z");
	ipopt->eval_grad_f = env->GetMethodID(solverCls, "eval_grad_f", "(I[DZ[D)Z");
	ipopt->eval_g = env->GetMethodID(solverCls, "eval_g", "(I[DZI[D)Z");
	ipopt->eval_jac_g = env->GetMethodID(solverCls, "eval_jac_g", "(I[DZII[I[I[D)Z");
	ipopt->eval_h = env->GetMethodID(solverCls, "eval_h", "(I[DZDI[DZI[I[I[D)Z");
	if(ipopt->eval_f==0 || ipopt->eval_grad_f==0 || ipopt->eval_g==0 || ipopt->eval_jac_g==0 ||
		ipopt->eval_h==0){
        std::cerr << "Expected callback methods missing on Ipopt.java" << std::endl;
		delete ipopt;
		return 0;
	}
	
	jdouble  *x_L, *x_U;
	x_L = env->GetDoubleArrayElements(ylb, 0);
	x_U = env->GetDoubleArrayElements(yub, 0);

	/* set the number of constraints and allocate space for the bounds */
	Number* g_L = env->GetDoubleArrayElements(glb, 0); /* lower bounds on g */
	Number* g_U = env->GetDoubleArrayElements(gub, 0); /* upper bounds on g */

	/* create the IpoptProblem */
	ipopt->problem = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
			index_style, &eval_f, &eval_g, &eval_grad_f,
			&eval_jac_g, &eval_h);
			
	if(ipopt->problem == NULL){
        std::cerr << "CreateIpoptProblem returned NULL" << std::endl;
		delete ipopt;
		return 0;
	}

	return (jlong)ipopt;
}

JNIEXPORT jint JNICALL Java_org_coinor_Ipopt_IpoptSolve
(JNIEnv *env,
jobject obj_this,

jlong pipopt,
jdoubleArray xj,
jdoubleArray gj,
jdoubleArray obj_valj,
jdoubleArray mult_gj,
jdoubleArray mult_x_Lj,
jdoubleArray mult_x_Uj,

jdoubleArray callback_grad_f,
jdoubleArray callback_jac_g,
jdoubleArray callback_hess){

	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;
	ipopt->solver = obj_this;
	
    if(xj==NULL || gj==NULL || obj_valj==NULL || mult_gj==NULL || mult_x_Lj==NULL || mult_x_Uj==NULL ||
    	callback_grad_f==NULL || callback_jac_g==NULL || callback_hess==NULL){
    	return 0; // out of memory exception already thrown
    }
    // store some of the arguments
    ipopt->xj = xj;
    ipopt->fj = obj_valj;
	ipopt->gj = gj;
    ipopt->grad_fj = callback_grad_f;
	ipopt->jac_gj = callback_jac_g;
	ipopt->hessj = callback_hess;
	ipopt->mult_gj = mult_gj;
	
	Number *x = env->GetDoubleArrayElements(xj, 0);
	Number *g = env->GetDoubleArrayElements(gj, 0);
	Number *obj_val = env->GetDoubleArrayElements(obj_valj, 0);
	Number *mult_g = env->GetDoubleArrayElements(mult_gj, 0);
	Number *mult_x_L = env->GetDoubleArrayElements(mult_x_Lj, 0);
	Number *mult_x_U = env->GetDoubleArrayElements(mult_x_Uj, 0);
	
	/* solve the problem */
	int status = IpoptSolve(ipopt->problem, x, g, obj_val, mult_g, mult_x_L, mult_x_U, (void*)ipopt);

	env->ReleaseDoubleArrayElements(xj, x, 0);
	env->ReleaseDoubleArrayElements(gj, g, 0);
	env->ReleaseDoubleArrayElements(obj_valj, obj_val, 0);
	env->ReleaseDoubleArrayElements(mult_gj, mult_g, 0);
	env->ReleaseDoubleArrayElements(mult_x_Lj, mult_x_L, 0);
	env->ReleaseDoubleArrayElements(mult_x_Uj, mult_x_U, 0);
	
	return status;
}


JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_SetIpoptProblemScaling
(JNIEnv *env,
jobject obj_this,

jlong pipopt,
jdouble obj_scaling,
jdoubleArray x_scalingj,
jdoubleArray g_scalingj)
{
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;
	ipopt->solver = obj_this;
	
	Number *x_scaling = x_scalingj ? env->GetDoubleArrayElements(x_scalingj, 0) : NULL;
	Number *g_scaling = g_scalingj ? env->GetDoubleArrayElements(g_scalingj, 0) : NULL;
	
	/* set the problem scaling */
	jboolean status = SetIpoptProblemScaling(ipopt->problem, obj_scaling, x_scaling, g_scaling);

	return status;
}


JNIEXPORT void JNICALL Java_org_coinor_Ipopt_FreeIpoptProblem
(JNIEnv *env, 
jobject obj_this, 
jlong pipopt){
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;
	
	if(ipopt->problem){
		FreeIpoptProblem(ipopt->problem);
		ipopt->problem = NULL;
	}
	delete ipopt;
}

/* IPOPT Function Implementations */
Bool eval_f(Index n, Number* x, Bool new_x,
            Number* obj_value, UserDataPtr user_data)
{
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)user_data;
	JNIEnv *env = ipopt->env;

	jdoubleArray xj = ipopt->xj;
	jdoubleArray fj = ipopt->fj;

  	if(new_x){
		/* Copy the native double x to the Java double array xj */ 
    	env->SetDoubleArrayRegion(xj, 0, n, x);
  	}
    
	/* Call the java method */
	jboolean new_xj = new_x;
    if(!env->CallBooleanMethod(ipopt->solver, ipopt->eval_f, n, xj, new_xj, fj))
    	return FALSE;
    
	/* Copy from Java to native value */
    env->GetDoubleArrayRegion(fj, 0, 1, obj_value);

    return TRUE;
}

Bool eval_grad_f(Index n, Number* x, Bool new_x,
	Number* grad_f, UserDataPtr user_data)
{
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)user_data;
	JNIEnv *env = ipopt->env;
	
	jdoubleArray xj = ipopt->xj;
	jdoubleArray grad_fj = ipopt->grad_fj;
  	
  	if(new_x){
		/* Copy the native double x to the Java double array xj */ 
    	env->SetDoubleArrayRegion(xj, 0, n, x);
  	}
	
	/* Call the java method */
	jboolean new_xj = new_x;
    if(!env->CallBooleanMethod(ipopt->solver, ipopt->eval_grad_f, n, xj, new_xj, grad_fj))
    	return FALSE;

	/* Copy from Java to native value */
    env->GetDoubleArrayRegion(grad_fj, 0, n, grad_f);
	
	return TRUE;
}

Bool eval_g(Index n, Number* x, Bool new_x,
	Index m, Number* g, UserDataPtr user_data)
{
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)user_data;
	JNIEnv *env = ipopt->env;

	jdoubleArray xj = ipopt->xj;
	jdoubleArray gj = ipopt->gj;
  	
  	if(new_x){
		/* Copy the native double x to the Java double array xj */ 
    	env->SetDoubleArrayRegion(xj, 0, n, x);
  	}
	
	/* Call the java method */
	jboolean new_xj = new_x;
    jboolean ret = env->CallBooleanMethod(ipopt->solver, ipopt->eval_g, n, xj, new_xj, m, gj);
    if(!ret)
    	return ret;
	/* Copy from Java to native value */
    env->GetDoubleArrayRegion(gj, 0, ipopt->m, g);
    
	return TRUE;
}

Bool eval_jac_g(Index n, Number *x, Bool new_x,
	Index m, Index nele_jac,
	Index *iRow, Index *jCol, Number *jac_g,
	UserDataPtr user_data)
{
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)user_data;
	JNIEnv *env = ipopt->env;
 
 	jdoubleArray xj = NULL;
	jdoubleArray jac_gj = NULL;
	xj = ipopt->xj;
	jac_gj = ipopt->jac_gj;
  	
  	if(new_x && x!=NULL){
		/* Copy the native double x to the Java double array xj */ 
    	env->SetDoubleArrayRegion(xj, 0, n, x);
  	}
  	
  	/// Create the index arrays if needed
  	jintArray iRowj = NULL;
  	jintArray jColj = NULL;
  	if(iRow != NULL && jCol != NULL){
  		iRowj = env->NewIntArray(ipopt->nele_jac);
  		jColj = env->NewIntArray(ipopt->nele_jac);
  	}
	
	/* Call the java method */
	jboolean new_xj = new_x;
    if(!env->CallBooleanMethod(ipopt->solver, ipopt->eval_jac_g, n, xj, new_xj,
    	m, nele_jac, iRowj, jColj, jac_g == NULL ? NULL : jac_gj))
    	return FALSE;
	/* Copy from Java to native value */
    if(jac_g != NULL)
    	env->GetDoubleArrayRegion(jac_gj, 0, ipopt->nele_jac, jac_g);
    if(iRow != NULL && jCol != NULL){
    	// env->GetIntArrayRegion(iRowj, 0, ipopt->nele_jac, iRow);
    	// env->GetIntArrayRegion(jColj, 0, ipopt->nele_jac, jCol);
		jint *iRow_jint = env->GetIntArrayElements(iRowj, 0);
		jint *jCol_jint = env->GetIntArrayElements(jColj, 0);
		for(int i=0; i<nele_jac; ++i){
			iRow[i] = iRow_jint[i];
			jCol[i] = jCol_jint[i];
		}
    }

	return TRUE;
}

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
	Index m, Number *lambda, Bool new_lambda,
	Index nele_hess, Index *iRow, Index *jCol,
	Number *hess, UserDataPtr user_data)
{
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)user_data;
	JNIEnv *env = ipopt->env;
 
 	jdoubleArray xj = ipopt->xj;
	jdoubleArray hessj = ipopt->hessj;
	jdoubleArray mult_gj = ipopt->mult_gj;
  	
  	if(new_x && x!=NULL){
		/* Copy the native double x to the Java double array xj */ 
    	env->SetDoubleArrayRegion(xj, 0, n, x);
  	}
  	if(new_lambda && lambda!=NULL){
		/* Copy the native double lambda to the Java double array lambdaj */ 
    	env->SetDoubleArrayRegion(mult_gj, 0, m, lambda);
  	}
  	
  	/// Create the index arrays if needed
  	jintArray iRowj = NULL;
  	jintArray jColj = NULL;
  	if(iRow != NULL && jCol != NULL){
  		iRowj = env->NewIntArray(ipopt->nele_hess);
  		jColj = env->NewIntArray(ipopt->nele_hess);
  	}

	/* Call the java method */
	jboolean new_xj = new_x;
	jboolean new_lambdaj = new_lambda;
    if(!env->CallBooleanMethod(ipopt->solver, ipopt->eval_h, n, xj, new_xj,
    	obj_factor, m, mult_gj, new_lambdaj,
		nele_hess, iRowj, jColj, hess == NULL ? NULL : hessj))
    	return FALSE;
    
	/* Copy from Java to native value */
    if(hess != NULL)
    	env->GetDoubleArrayRegion(hessj, 0, ipopt->nele_hess, hess);
    if(iRow != NULL && jCol != NULL){
    	// env->GetIntArrayRegion(iRowj, 0, ipopt->nele_hess, iRow);
    	// env->GetIntArrayRegion(jColj, 0, ipopt->nele_hess, jCol);
 		jint *iRow_jint = env->GetIntArrayElements(iRowj, 0);
		jint *jCol_jint = env->GetIntArrayElements(jColj, 0);
		for(int i=0; i<nele_hess; ++i){
			iRow[i] = iRow_jint[i];
			jCol[i] = jCol_jint[i];
		}
   }

	return TRUE;
}

JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_OpenIpoptOutputFile
(JNIEnv * env, jobject obj_this, jlong pipopt, jstring file_namej, jint print_level){
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;
	
	const char *file_name = env->GetStringUTFChars(file_namej, NULL);

	// Try to apply the string option
	jboolean ret = OpenIpoptOutputFile(ipopt->problem, (char*)file_name, print_level);		
	
	env->ReleaseStringUTFChars(file_namej, file_name);
	
	return ret;
}

JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_AddIpoptIntOption
(JNIEnv * env, jobject obj_this, jlong pipopt, jstring jparname, jint jparvalue){
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;
	
	const char *pparameterName = env->GetStringUTFChars(jparname, 0);

	// Try to apply the integer option
	jboolean ret = AddIpoptIntOption(ipopt->problem, (char*)pparameterName, jparvalue);		

	env->ReleaseStringUTFChars(jparname, pparameterName);
	
	return ret;
}

JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_AddIpoptNumOption
(JNIEnv * env, jobject obj_this, jlong pipopt, jstring jparname, jdouble jparvalue){  
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;
	
	const char *pparameterName = env->GetStringUTFChars(jparname, 0);

	// Try to set the real option
	jboolean ret = AddIpoptNumOption(ipopt->problem, (char*)pparameterName, jparvalue);
	env->ReleaseStringUTFChars(jparname, pparameterName);
	
	return ret;
}

JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_AddIpoptStrOption
(JNIEnv * env, jobject obj_this, jlong pipopt, jstring jparname, jstring jparvalue){
	// cast back our structure
	ipopt_jni *ipopt = (ipopt_jni *)pipopt;
	
	const char *pparameterName = env->GetStringUTFChars(jparname, NULL);
	const char *pparameterValue = env->GetStringUTFChars(jparvalue, NULL);

	// Try to apply the string option
	jboolean ret = AddIpoptStrOption(ipopt->problem, (char*)pparameterName, (char*)pparameterValue);		
	
	env->ReleaseStringUTFChars(jparname, pparameterName);
	env->ReleaseStringUTFChars(jparname, pparameterValue);
	
	return ret;
}

#ifdef __cplusplus
}
#endif
