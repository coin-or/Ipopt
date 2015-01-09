/** Copyright (C) 2007 VRTech Industrial Technologies - www.vrtech.com.br.
 * Copyright (C) 2007 Tong Kewei, BeiHang University - www.buaa.edu.cn.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * $Id$
 */

#include <cassert>
#include <jni.h>
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "org_coinor_Ipopt.h"

using namespace std;
using namespace Ipopt;

/** Main structure for Ipopt JNI implementation.
 * 
 * All functions will receive a pointer to this structure as
 * an integer argument (the address in memory of the structure).
 */
class Jipopt : public TNLP
{
public:
   /**  constructor */
   Jipopt(JNIEnv *env, jobject solver, jint n, jint m, jint nele_jac, jint nele_hess, jint index_style);

   /** default destructor */
   virtual ~Jipopt() { }

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
   //@}

public:
   // The JNI Environment
   JNIEnv *env;
   jobject solver;

   jint n;
   jint m;
   jint nele_jac;
   jint nele_hess;
   jint index_style;

   // some cached arguments
   jdoubleArray mult_gj;
   jdoubleArray mult_x_Lj;
   jdoubleArray mult_x_Uj;

   // the callback arguments
   jdoubleArray xj;
   jdoubleArray fj;
   jdoubleArray grad_fj;
   jdoubleArray gj;
   jdoubleArray jac_gj;
   jdoubleArray hessj;

   jboolean using_scaling_parameters;
   jboolean using_LBFGS;

   SmartPtr<IpoptApplication> application;

   // the callback methods
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
   Jipopt(const Jipopt&);
   Jipopt& operator=(const Jipopt&);
};

Jipopt::Jipopt(JNIEnv *env_, jobject solver_, jint n_, jint m_, jint nele_jac_, jint nele_hess_, jint index_style_)
: env(env_), solver(solver_), n(n_), m(m_), nele_jac(nele_jac_), nele_hess(nele_hess_), index_style(index_style_),
  mult_gj(NULL), mult_x_Lj(NULL), mult_x_Uj(NULL),
  xj(NULL), fj(NULL), grad_fj(NULL), gj(NULL), jac_gj(NULL), hessj(NULL),
  using_scaling_parameters(false), using_LBFGS(false)
{
   application = new IpoptApplication();
   application->RethrowNonIpoptException(false);

   // the solver class
   jclass solverCls = env->GetObjectClass(solver);

   // get the methods
   get_bounds_info_ = env->GetMethodID(solverCls,"get_bounds_info","(I[D[DI[D[D)Z");
   get_starting_point_ = env->GetMethodID(solverCls,"get_starting_point","(IZ[DZ[D[DIZ[D)Z");
   eval_f_ = env->GetMethodID(solverCls, "eval_f", "(I[DZ[D)Z");
   eval_grad_f_ = env->GetMethodID(solverCls, "eval_grad_f", "(I[DZ[D)Z");
   eval_g_ = env->GetMethodID(solverCls, "eval_g", "(I[DZI[D)Z");
   eval_jac_g_ = env->GetMethodID(solverCls, "eval_jac_g", "(I[DZII[I[I[D)Z");
   eval_h_ = env->GetMethodID(solverCls, "eval_h", "(I[DZDI[DZI[I[I[D)Z");
   get_scaling_parameters_ = env->GetMethodID(solverCls,"get_scaling_parameters","([DI[DI[D[Z)Z");
   get_number_of_nonlinear_variables_ = env->GetMethodID(solverCls,"get_number_of_nonlinear_variables","()I");
   get_list_of_nonlinear_variables_ = env->GetMethodID(solverCls,"get_list_of_nonlinear_variables","(I[I)Z");

   if( get_bounds_info_ == 0 || get_starting_point_ == 0 ||
      eval_f_ == 0 || eval_grad_f_ == 0 || eval_g_ == 0 || eval_jac_g_ == 0 || eval_h_ == 0 ||
      get_scaling_parameters_ == 0 || get_number_of_nonlinear_variables_ == 0 || get_list_of_nonlinear_variables_ == 0 )
      std::cerr << "Expected callback methods missing on JIpopt.java" << std::endl;
}

bool Jipopt::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
   Index& nnz_h_lag, IndexStyleEnum& index_style)
{
   n = this->n;
   m = this->m;
   nnz_jac_g = this->nele_jac;
   nnz_h_lag = this->nele_hess;

   index_style = (IndexStyleEnum)this->index_style;

   return true;
}

bool Jipopt::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u)
{
   jdoubleArray x_lj = NULL;
   jdoubleArray x_uj = NULL;
   jdoubleArray g_lj = NULL;
   jdoubleArray g_uj = NULL;

   assert(x_l != NULL);
   assert(x_u != NULL);
   assert(g_l != NULL);
   assert(g_u != NULL);

   x_lj = env->NewDoubleArray(n);
   x_uj = env->NewDoubleArray(n);
   g_lj = env->NewDoubleArray(m);
   g_uj = env->NewDoubleArray(m);

   if( !env->CallBooleanMethod(solver, get_bounds_info_, n, x_lj, x_uj, m, g_lj, g_uj) )
      return false;

   // Copy from Java to native value
   env->GetDoubleArrayRegion(x_lj, 0, n, x_l);
   env->GetDoubleArrayRegion(x_uj, 0, n, x_u);
   env->GetDoubleArrayRegion(g_lj, 0, m, g_l);
   env->GetDoubleArrayRegion(g_uj, 0, m, g_u);

   return true;
}

bool Jipopt::get_starting_point(Index n, bool init_x, Number* x,
   bool init_z, Number* z_L, Number* z_U,
   Index m, bool init_lambda, Number* lambda)
{
   jdoubleArray xj      = this->xj;
   jdoubleArray z_lj    = this->mult_x_Lj;
   jdoubleArray z_uj    = this->mult_x_Uj;
   jdoubleArray lambdaj = this->mult_gj;

   if( !env->CallBooleanMethod(solver, get_starting_point_, n, init_x, xj, init_z, z_lj, z_uj, m, init_lambda, lambdaj) )
      return false;

   /* Copy from Java to native value */
   if( init_x )
      env->GetDoubleArrayRegion(xj, 0, n, x);

   if( init_z )
   {
      env->GetDoubleArrayRegion(z_lj, 0, n, z_L);
      env->GetDoubleArrayRegion(z_uj, 0, n, z_U);
   }

   if( init_lambda )
      env->GetDoubleArrayRegion(lambdaj, 0, m, lambda);

   return true;
}

bool Jipopt::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
   /* Copy the native double x to the Java double array xj, if new values */
   if( new_x )
      env->SetDoubleArrayRegion(xj, 0, n, const_cast<Number*>(x));

   /* Call the java method */
   jboolean new_xj = new_x;
   if( !env->CallBooleanMethod(solver, eval_f_, n, xj, new_xj, fj) )
      return false;

   /* Copy from Java to native value */
   env->GetDoubleArrayRegion(fj, 0, 1, &obj_value);

   return true;
}

bool Jipopt::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
   /* Copy the native double x to the Java double array xj, if new values */
   if( new_x )
      env->SetDoubleArrayRegion(xj, 0, n, const_cast<Number*>(x));

   /* Call the java method */
   jboolean new_xj = new_x;
   if( !env->CallBooleanMethod(solver, eval_grad_f_, n, xj, new_xj, grad_fj) )
      return false;

   env->GetDoubleArrayRegion(grad_fj, 0, n, grad_f);

   return true;
}

bool Jipopt::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
   /* Copy the native double x to the Java double array xj, if new values */
   if( new_x )
      env->SetDoubleArrayRegion(xj, 0, n, const_cast<Number*>(x));

   /* Call the java method */
   jboolean new_xj = new_x;
   if( !env->CallBooleanMethod(solver, eval_g_, n, xj, new_xj, m, gj) )
      return false;

   /* Copy from Java to native value */
   env->GetDoubleArrayRegion(gj, 0, m, g);

   return true;
}

bool Jipopt::eval_jac_g(Index n, const Number* x, bool new_x,
   Index m, Index nele_jac, Index* iRow,
   Index *jCol, Number* jac_g)
{
   // Copy the native double x to the Java double array xj, if new values
   if( new_x && x != NULL )
      env->SetDoubleArrayRegion(xj, 0, n, const_cast<Number*>(x));

   /// Create the index arrays if needed
   jintArray iRowj = NULL;
   jintArray jColj = NULL;
   if( iRow != NULL && jCol != NULL )
   {
      iRowj = env->NewIntArray(nele_jac);
      jColj = env->NewIntArray(nele_jac);
   }

   /* Call the java method */
   jboolean new_xj = new_x;
   if( !env->CallBooleanMethod(solver, eval_jac_g_, n, xj, new_xj, m, nele_jac, iRowj, jColj, jac_g == NULL ? NULL : jac_gj) )
      return false;

   /* Copy from Java to native value */
   if( jac_g != NULL )
      env->GetDoubleArrayRegion(jac_gj, 0, nele_jac, jac_g);

   if( iRow != NULL && jCol != NULL )
   {
      if( sizeof(jint) == sizeof(Index) )
      {
         env->GetIntArrayRegion(iRowj, 0, nele_jac, reinterpret_cast<jint*>(iRow));
         env->GetIntArrayRegion(jColj, 0, nele_jac, reinterpret_cast<jint*>(jCol));
      }
      else
      {
         jint* tmp = new jint[nele_jac];

         env->GetIntArrayRegion(iRowj, 0, nele_jac, tmp);
         for( int i = 0; i < nele_jac; ++i )
            iRow[i] = (Index)tmp[i];

         env->GetIntArrayRegion(jColj, 0, nele_jac, tmp);
         for( int i = 0; i < nele_jac; ++i )
            jCol[i] = (Index)tmp[i];

         delete[] tmp;
      }
   }

   return true;
}

bool Jipopt::eval_h(Index n, const Number* x, bool new_x,
   Number obj_factor, Index m, const Number* lambda,
   bool new_lambda, Index nele_hess,
   Index* iRow, Index* jCol, Number* hess)
{
   /* Copy the native double x to the Java double array xj, if new values */
   if( new_x && x != NULL )
      env->SetDoubleArrayRegion(xj, 0, n, const_cast<Number*>(x));

   /* Copy the native double lambda to the Java double array lambdaj, if new values */
   if( new_lambda && lambda != NULL )
      env->SetDoubleArrayRegion(mult_gj, 0, m, const_cast<Number*>(lambda));

   /* Create the index arrays if needed */
   jintArray iRowj = NULL;
   jintArray jColj = NULL;
   if( iRow != NULL && jCol != NULL )
   {
      iRowj = env->NewIntArray(nele_hess);
      jColj = env->NewIntArray(nele_hess);
   }

   /* Call the java method */
   jboolean new_xj = new_x;
   jboolean new_lambdaj = new_lambda;
   if( !env->CallBooleanMethod(solver, eval_h_, n, xj, new_xj, obj_factor, m, mult_gj, new_lambdaj, nele_hess, iRowj, jColj, hess == NULL ? NULL : hessj) )
      return false;

   /* Copy from Java to native value */
   if( hess != NULL )
      env->GetDoubleArrayRegion(hessj, 0, nele_hess, hess);

   if( iRow != NULL && jCol != NULL )
   {
      if( sizeof(jint) == sizeof(Index) )
      {
         env->GetIntArrayRegion(iRowj, 0, nele_hess, reinterpret_cast<jint*>(iRow));
         env->GetIntArrayRegion(jColj, 0, nele_hess, reinterpret_cast<jint*>(jCol));
      }
      else
      {
         jint* tmp = new jint[nele_hess];

         env->GetIntArrayRegion(iRowj, 0, nele_hess, tmp);
         for( int i = 0; i < nele_hess; ++i )
            iRow[i] = (Index)tmp[i];

         env->GetIntArrayRegion(jColj, 0, nele_hess, tmp);
         for( int i = 0; i < nele_hess; ++i )
            jCol[i] = (Index)tmp[i];

         delete[] tmp;
      }
   }

   return true;
}

void Jipopt::finalize_solution(SolverReturn status, Index n, const Number *x, 
   const Number *z_L, const Number *z_U, Index m,
   const Number *g, const Number *lambda, Number obj_value,
   const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq)
{
   /* Copy the native arrays to Java double arrays */

   if( x != NULL )
      env->SetDoubleArrayRegion(xj, 0, n, const_cast<Number*>(x));

   if( z_L != NULL )
      env->SetDoubleArrayRegion(mult_x_Lj, 0, n, const_cast<Number*>(z_L));

   if( z_U != NULL )
      env->SetDoubleArrayRegion(mult_x_Uj, 0, n, const_cast<Number*>(z_U));

   if( g != NULL )
      env->SetDoubleArrayRegion(gj, 0, m, const_cast<Number*>(g));

   if( lambda != NULL )
      env->SetDoubleArrayRegion(mult_gj, 0, m, const_cast<Number*>(lambda));

   env->GetDoubleArrayRegion(fj, 0, 1, &obj_value);
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
   if( using_scaling_parameters )
   {
      jdoubleArray obj_scaling_j = env->NewDoubleArray(1);
      jdoubleArray x_scaling_j   = env->NewDoubleArray(n);
      jdoubleArray g_scaling_j   = env->NewDoubleArray(m);

      jbooleanArray use_x_g_scaling_j = env->NewBooleanArray(2);

      env->CallBooleanMethod(solver, get_scaling_parameters_,
         obj_scaling_j,
         n, x_scaling_j,
         m, g_scaling_j,
         use_x_g_scaling_j);

      jboolean* use_x_g_scaling = env->GetBooleanArrayElements(use_x_g_scaling_j,0);

      /* Copy from Java to native value */
      env->GetDoubleArrayRegion(obj_scaling_j, 0, 1, &obj_scaling);

      /* Copy from Java to native value */
      if( use_x_g_scaling[0] )
      {
         env->GetDoubleArrayRegion(x_scaling_j, 0, n, x_scaling);
         use_x_scaling = true;
      }
      else
      {
         use_x_scaling = false;
      }

      /* Copy from Java to native value */
      if( use_x_g_scaling[1] )
      {
         env->GetDoubleArrayRegion(g_scaling_j, 0, n, g_scaling);
         use_g_scaling = true;
      }
      else
      {
         use_g_scaling=false;
      }

      env->ReleaseBooleanArrayElements(use_x_g_scaling_j, use_x_g_scaling, 0);

      return true;
   }

   return false;
}

Index Jipopt::get_number_of_nonlinear_variables()
{
   if( using_LBFGS )
      return env->CallIntMethod(solver, get_number_of_nonlinear_variables_);

   return -1;
}

bool Jipopt::get_list_of_nonlinear_variables(Index num_nonlin_vars,Index* pos_nonlin_vars)
{
   if( using_LBFGS )
   {
      jintArray pos_nonlin_vars_j = env->NewIntArray(num_nonlin_vars);

      if( !env->CallBooleanMethod(solver, get_list_of_nonlinear_variables_, num_nonlin_vars, pos_nonlin_vars_j) )
         return false;

      if( pos_nonlin_vars != NULL )
      {
         if( sizeof(jint) == sizeof(Index) )
         {
            env->GetIntArrayRegion(pos_nonlin_vars_j, 0, num_nonlin_vars, reinterpret_cast<jint*>(pos_nonlin_vars));
         }
         else
         {
            jint* tmp = new jint[num_nonlin_vars];

            env->GetIntArrayRegion(pos_nonlin_vars_j, 0, num_nonlin_vars, tmp);
            for( int i = 0; i < num_nonlin_vars; ++i )
               pos_nonlin_vars[i] = (Index)tmp[i];

            delete[] tmp;
         }
      }

      return true;
   }

   return false;
}

extern "C" {

JNIEXPORT jlong JNICALL Java_org_coinor_Ipopt_CreateIpoptProblem(
   JNIEnv *env, jobject obj_this,
   jint n,  jint m,
   jint nele_jac, jint nele_hess,
   jint index_style)
{
   /* create the smart pointer to the Ipopt problem */
   SmartPtr<Jipopt>* pproblem = new SmartPtr<Jipopt>;

   /* create the IpoptProblem */
   *pproblem = new Jipopt(env, obj_this, n, m, nele_jac, nele_hess, index_style);

   /* return the smart pointer to our class */
   return (jlong)pproblem;
}

JNIEXPORT jint JNICALL Java_org_coinor_Ipopt_OptimizeTNLP(
   JNIEnv *env, jobject obj_this, jlong pipopt,
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
   Jipopt* problem = GetRawPtr(*(SmartPtr<Jipopt>*)pipopt);

   problem->env = env;
   problem->solver = obj_this;

   problem->xj = xj;
   problem->gj = gj;
   problem->fj = obj_valj;
   problem->mult_gj = mult_gj;
   problem->mult_x_Lj = mult_x_Lj;
   problem->mult_x_Uj = mult_x_Uj;

   problem->grad_fj = callback_grad_f;
   problem->jac_gj = callback_jac_g;
   problem->hessj = callback_hess;

   //  (use a SmartPtr, not raw)
   ApplicationReturnStatus status;

   status = problem->application->Initialize();

   if( status != Solve_Succeeded )
   {
      printf("\n\n*** Error during initialization!\n");
      return (int) status;
   }

   /* solve the problem */
   status = problem->application->OptimizeTNLP(problem);

   return (jint) status;
}

JNIEXPORT void JNICALL Java_org_coinor_Ipopt_FreeIpoptProblem(JNIEnv *env,  jobject obj_this, jlong pipopt)
{
   SmartPtr<Jipopt>* pproblem = (SmartPtr<Jipopt>*)pipopt;

   if( pproblem != NULL && IsValid(*pproblem) )
   {
      /* if OptimizeTNLP has been called, the application holds a SmartPtr to our problem class
       * to resolve this circular dependency we first free the application explicitly
       */
      (*pproblem)->application = NULL;

      /* now free our JIpopt itself */
      *pproblem = NULL;
   }
}

JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_AddIpoptIntOption(
   JNIEnv * env, jobject obj_this, jlong pipopt, jstring jparname, jint jparvalue)
{
   Jipopt* problem = GetRawPtr(*(SmartPtr<Jipopt>*)pipopt);

   const char* pparameterName = env->GetStringUTFChars(jparname, 0);
   string parameterName = pparameterName;

   // Try to apply the integer option
   jboolean ret = problem->application->Options()->SetIntegerValue(parameterName, jparvalue);

   env->ReleaseStringUTFChars(jparname, pparameterName);

   return ret;
}

JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_AddIpoptNumOption(
   JNIEnv * env, jobject obj_this, jlong pipopt, jstring jparname, jdouble jparvalue)
{
   Jipopt* problem = GetRawPtr(*(SmartPtr<Jipopt>*)pipopt);

   const char* pparameterName = env->GetStringUTFChars(jparname, 0);
   string parameterName=pparameterName;

   // Try to set the real option
   jboolean ret = problem->application->Options()->SetNumericValue(parameterName,jparvalue);

   env->ReleaseStringUTFChars(jparname, pparameterName);

   return ret;
}

JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_AddIpoptStrOption(
   JNIEnv * env, jobject obj_this, jlong pipopt, jstring jparname, jstring jparvalue)
{
   Jipopt* problem = GetRawPtr(*(SmartPtr<Jipopt>*)pipopt);

   const char* pparameterName = env->GetStringUTFChars(jparname, NULL);
   string parameterName = pparameterName;
   const char* pparameterValue = env->GetStringUTFChars(jparvalue, NULL);
   string parameterValue = pparameterValue;

   // parameterValue has been changed to LowerCase in Java!
   if( parameterName == "hessian_approximation" && (parameterValue == "limited-memory") )
      problem->using_LBFGS = true;
   else if( parameterName == "nlp_scaling_method" && (parameterValue == "user-scaling") )
      problem->using_scaling_parameters = true;

   // Try to apply the string option
   jboolean ret = problem->application->Options()->SetStringValue(parameterName, parameterValue);

   env->ReleaseStringUTFChars(jparname, pparameterName);
   env->ReleaseStringUTFChars(jparname, pparameterValue);

   return ret;
}

} // extern "C"
