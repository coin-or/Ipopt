/* Copyright (C) 2007 VRTech Industrial Technologies - www.vrtech.com.br.
 * Copyright (C) 2007 Tong Kewei, BeiHang University - www.buaa.edu.cn.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 */

#include <cassert>
#include <jni.h>
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "org_coinor_Ipopt.h"

using namespace std;
using namespace Ipopt;

#ifdef IPOPT_SINGLE
typedef jfloat jnumber;
typedef jfloatArray jnumberArray;
#define NewNumberArray NewFloatArray
#define GetNumberArrayRegion GetFloatArrayRegion
#define SetNumberArrayRegion SetFloatArrayRegion
#else
typedef jdouble jnumber;
typedef jdoubleArray jnumberArray;
#define NewNumberArray NewDoubleArray
#define GetNumberArrayRegion GetDoubleArrayRegion
#define SetNumberArrayRegion SetDoubleArrayRegion
#endif

/** Main structure for Ipopt JNI implementation.
 *
 * All functions will receive a pointer to this structure as
 * an integer argument (the address in memory of the structure).
 */
class Jipopt: public TNLP
{
public:
   /** constructor */
   Jipopt(
      JNIEnv* env,        /**< JNI environment */
      jobject solver,     /**< Java Ipopt solver object */
      jint    n,          /**< number of variables */
      jint    m,          /**< number of constraints */
      jint    nele_jac,   /**< number of nonzero elements in Jacobian */
      jint    nele_hess,  /**< number of nonzero elements in Hessian */
      jint    index_style /**< whether to use C or Fortran index style (0 for C, 1 for Fortran) */
   );

   /** default destructor */
   virtual ~Jipopt()
   {
   }

   /**@name Overloaded from TNLP */
   ///@{
   /* Method to return some info about the NLP */
   virtual bool get_nlp_info(
      Index&          n,
      Index&          m,
      Index&          nnz_jac_g,
      Index&          nnz_h_lag,
      IndexStyleEnum& index_style
   );

   /* Method to return the bounds for my problem */
   virtual bool get_bounds_info(
      Index   n,
      Number* x_l,
      Number* x_u,
      Index   m,
      Number* g_l,
      Number* g_u
   );

   /* Method to return the starting point for the algorithm */
   virtual bool get_starting_point(
      Index   n,
      bool    init_x,
      Number* x,
      bool    init_z,
      Number* z_L,
      Number* z_U,
      Index   m,
      bool    init_lambda,
      Number* lambda
   );

   /* Method to return the objective value */
   virtual bool eval_f(
      Index         n,
      const Number* x,
      bool          new_x,
      Number&       obj_value
   );

   /* Method to return the gradient of the objective */
   virtual bool eval_grad_f(
      Index         n,
      const Number* x,
      bool          new_x,
      Number*       grad_f
   );

   /* Method to return the constraint activities */
   virtual bool eval_g(
      Index         n,
      const Number* x,
      bool          new_x,
      Index         m,
      Number*       g
   );

   /* Method to return:
    *  1) The structure of the Jacobian (if "values" is NULL)
    *  2) The values of the Jacobian (if "values" is not NULL)
    */
   virtual bool eval_jac_g(
      Index         n,
      const Number* x,
      bool          new_x,
      Index         m,
      Index         nele_jac,
      Index*        iRow,
      Index*        jCol,
      Number*       values
   );

   /* Method to return:
    *  1) The structure of the Hessian of the Lagrangian (if "values" is NULL)
    *  2) The values of the Hessian of the Lagrangian (if "values" is not NULL)
    */
   virtual bool eval_h(
      Index         n,
      const Number* x,
      bool          new_x,
      Number        obj_factor,
      Index         m,
      const Number* lambda,
      bool          new_lambda,
      Index         nele_hess,
      Index*        iRow,
      Index*        jCol,
      Number*       values
   );

   /* This method is called when the algorithm is complete so the TNLP can store/write the solution */
   virtual void finalize_solution(
      SolverReturn               status,
      Index                      n,
      const Number*              x,
      const Number*              z_L,
      const Number*              z_U,
      Index                      m,
      const Number*              g,
      const Number*              lambda,
      Number                     obj_value,
      const IpoptData*           ip_data,
      IpoptCalculatedQuantities* ip_cq
   );

   /* Intermediate Callback method for the user. */
   virtual bool intermediate_callback(
      AlgorithmMode              mode,
      Index                      iter,
      Number                     obj_value,
      Number                     inf_pr,
      Number                     inf_du,
      Number                     mu,
      Number                     d_norm,
      Number                     regularization_size,
      Number                     alpha_du,
      Number                     alpha_pr,
      Index                      ls_trials,
      const IpoptData*           ip_data,
      IpoptCalculatedQuantities* ip_cq
   );

   /* Method to return scaling parameters. */
   virtual bool get_scaling_parameters(
      Number&  obj_scaling,
      bool&    use_x_scaling,
      Index    n,
      Number*  x_scaling,
      bool&    use_g_scaling,
      Index    m,
      Number*  g_scaling
   );

   /* Method for quasi-Newton approximation. */
   virtual Index get_number_of_nonlinear_variables();

   virtual bool get_list_of_nonlinear_variables(
      Index  num_nonlin_vars,
      Index* pos_nonlin_vars
   );
   ///@}

public:
   /// The JNI Environment
   JNIEnv* env;
   jobject solver;

   jint n;
   jint m;
   jint nele_jac;
   jint nele_hess;
   jint index_style;

   // some cached arguments
   jnumberArray mult_gj;
   jnumberArray mult_x_Lj;
   jnumberArray mult_x_Uj;

   // the callback arguments
   jnumberArray xj;
   jnumberArray fj;
   jnumberArray grad_fj;
   jnumberArray gj;
   jnumberArray jac_gj;
   jnumberArray hessj;

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

   jmethodID intermediate_callback_;
   jmethodID get_scaling_parameters_;
   jmethodID get_number_of_nonlinear_variables_;
   jmethodID get_list_of_nonlinear_variables_;

private:
   Jipopt(const Jipopt&);
   Jipopt& operator=(const Jipopt&);
};

Jipopt::Jipopt(
   JNIEnv* env_,
   jobject solver_,
   jint    n_,
   jint    m_,
   jint    nele_jac_,
   jint    nele_hess_,
   jint    index_style_
)
   : env(env_), solver(solver_), n(n_), m(m_), nele_jac(nele_jac_), nele_hess(nele_hess_), index_style(index_style_),
     mult_gj(NULL), mult_x_Lj(NULL), mult_x_Uj(NULL), xj(NULL), fj(NULL), grad_fj(NULL),
     gj(NULL), jac_gj(NULL), hessj(NULL),
     using_scaling_parameters(false), using_LBFGS(false),
     application(new IpoptApplication())
{
   application->RethrowNonIpoptException(false);

   // the solver class
   jclass solverCls = env->GetObjectClass(solver);

   // get the methods
#ifndef IPOPT_SINGLE
   get_bounds_info_    = env->GetMethodID(solverCls, "get_bounds_info", "(I[D[DI[D[D)Z");
   get_starting_point_ = env->GetMethodID(solverCls, "get_starting_point", "(IZ[DZ[D[DIZ[D)Z");
   eval_f_             = env->GetMethodID(solverCls, "eval_f", "(I[DZ[D)Z");
   eval_grad_f_        = env->GetMethodID(solverCls, "eval_grad_f", "(I[DZ[D)Z");
   eval_g_             = env->GetMethodID(solverCls, "eval_g", "(I[DZI[D)Z");
   eval_jac_g_         = env->GetMethodID(solverCls, "eval_jac_g", "(I[DZII[I[I[D)Z");
   eval_h_             = env->GetMethodID(solverCls, "eval_h", "(I[DZDI[DZI[I[I[D)Z");
   get_scaling_parameters_ = env->GetMethodID(solverCls, "get_scaling_parameters", "([DI[DI[D[Z)Z");
   intermediate_callback_ = env->GetMethodID(solverCls, "intermediate_callback", "(IIDDDDDDDDIJJ)Z");
#else
   get_bounds_info_    = env->GetMethodID(solverCls, "get_bounds_info", "(I[F[FI[F[F)Z");
   get_starting_point_ = env->GetMethodID(solverCls, "get_starting_point", "(IZ[FZ[F[FIZ[F)Z");
   eval_f_             = env->GetMethodID(solverCls, "eval_f", "(I[FZ[F)Z");
   eval_grad_f_        = env->GetMethodID(solverCls, "eval_grad_f", "(I[FZ[F)Z");
   eval_g_             = env->GetMethodID(solverCls, "eval_g", "(I[FZI[F)Z");
   eval_jac_g_         = env->GetMethodID(solverCls, "eval_jac_g", "(I[FZII[I[I[F)Z");
   eval_h_             = env->GetMethodID(solverCls, "eval_h", "(I[FZFI[FZI[I[I[F)Z");
   get_scaling_parameters_ = env->GetMethodID(solverCls, "get_scaling_parameters", "([FI[FI[F[Z)Z");
   intermediate_callback_ = env->GetMethodID(solverCls, "intermediate_callback", "(IIFFFFFFFFIJJ)Z");
#endif
   get_number_of_nonlinear_variables_ = env->GetMethodID(solverCls, "get_number_of_nonlinear_variables", "()I");
   get_list_of_nonlinear_variables_   = env->GetMethodID(solverCls, "get_list_of_nonlinear_variables", "(I[I)Z");

   if( get_bounds_info_ == 0 || get_starting_point_ == 0 || eval_f_ == 0
       || eval_grad_f_ == 0 || eval_g_ == 0 || eval_jac_g_ == 0 || eval_h_ == 0
       || get_scaling_parameters_ == 0 || get_number_of_nonlinear_variables_ == 0
       || get_list_of_nonlinear_variables_ == 0 )
   {
      std::cerr << "Expected callback methods missing on JIpopt.java" << std::endl;
   }

   assert(get_bounds_info_    != 0);
   assert(get_starting_point_ != 0);
   assert(eval_f_      != 0);
   assert(eval_grad_f_ != 0);
   assert(eval_g_      != 0);
   assert(eval_jac_g_  != 0);
   assert(eval_h_      != 0);
   assert(get_scaling_parameters_ != 0);
   assert(get_number_of_nonlinear_variables_ != 0);
   assert(get_list_of_nonlinear_variables_   != 0);
}

bool Jipopt::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style)
{
   n = this->n;
   m = this->m;
   nnz_jac_g = this->nele_jac;
   nnz_h_lag = this->nele_hess;

   index_style = (IndexStyleEnum) this->index_style;

   return true;
}

bool Jipopt::get_bounds_info(
   Index   n,
   Number* x_l,
   Number* x_u,
   Index   m,
   Number* g_l,
   Number* g_u)
{
   jnumberArray x_lj = NULL;
   jnumberArray x_uj = NULL;
   jnumberArray g_lj = NULL;
   jnumberArray g_uj = NULL;

   assert(x_l != NULL);
   assert(x_u != NULL);
   assert(g_l != NULL);
   assert(g_u != NULL);

   x_lj = env->NewNumberArray(n);
   x_uj = env->NewNumberArray(n);
   g_lj = env->NewNumberArray(m);
   g_uj = env->NewNumberArray(m);

   if( !env->CallBooleanMethod(solver, get_bounds_info_, n, x_lj, x_uj, m, g_lj, g_uj) )
   {
      return false;
   }

   // Copy from Java to native value
   env->GetNumberArrayRegion(x_lj, 0, n, x_l);
   env->GetNumberArrayRegion(x_uj, 0, n, x_u);
   env->GetNumberArrayRegion(g_lj, 0, m, g_l);
   env->GetNumberArrayRegion(g_uj, 0, m, g_u);

   return true;
}

bool Jipopt::get_starting_point(
   Index   n,
   bool    init_x,
   Number* x,
   bool    init_z,
   Number* z_L,
   Number* z_U,
   Index   m,
   bool    init_lambda,
   Number* lambda)
{
   jnumberArray xj      = this->xj;
   jnumberArray z_lj    = this->mult_x_Lj;
   jnumberArray z_uj    = this->mult_x_Uj;
   jnumberArray lambdaj = this->mult_gj;

   if( !env->CallBooleanMethod(solver, get_starting_point_, n, init_x, xj, init_z, z_lj, z_uj, m, init_lambda, lambdaj) )
   {
      return false;
   }

   /* Copy from Java to native value */
   if( init_x )
   {
      env->GetNumberArrayRegion(xj, 0, n, x);
   }

   if( init_z )
   {
      env->GetNumberArrayRegion(z_lj, 0, n, z_L);
      env->GetNumberArrayRegion(z_uj, 0, n, z_U);
   }

   if( init_lambda )
   {
      env->GetNumberArrayRegion(lambdaj, 0, m, lambda);
   }

   return true;
}

bool Jipopt::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value)
{
   /* Copy the native double x to the Java double array xj, if new values */
   if( new_x )
   {
      env->SetNumberArrayRegion(xj, 0, n, const_cast<Number*>(x));
   }

   /* Call the java method */
   jboolean new_xj = new_x;
   if( !env->CallBooleanMethod(solver, eval_f_, n, xj, new_xj, fj) )
   {
      return false;
   }

   /* Copy from Java to native value */
   env->GetNumberArrayRegion(fj, 0, 1, &obj_value);

   return true;
}

bool Jipopt::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f)
{
   /* Copy the native double x to the Java double array xj, if new values */
   if( new_x )
   {
      env->SetNumberArrayRegion(xj, 0, n, const_cast<Number*>(x));
   }

   /* Call the java method */
   jboolean new_xj = new_x;
   if( !env->CallBooleanMethod(solver, eval_grad_f_, n, xj, new_xj, grad_fj) )
   {
      return false;
   }

   env->GetNumberArrayRegion(grad_fj, 0, n, grad_f);

   return true;
}

bool Jipopt::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g)
{
   /* Copy the native double x to the Java double array xj, if new values */
   if( new_x )
   {
      env->SetNumberArrayRegion(xj, 0, n, const_cast<Number*>(x));
   }

   /* Call the java method */
   jboolean new_xj = new_x;
   if( !env->CallBooleanMethod(solver, eval_g_, n, xj, new_xj, m, gj) )
   {
      return false;
   }

   /* Copy from Java to native value */
   env->GetNumberArrayRegion(gj, 0, m, g);

   return true;
}

bool Jipopt::eval_jac_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Index         nele_jac,
   Index*        iRow,
   Index*        jCol,
   Number*       jac_g)
{
   // Copy the native double x to the Java double array xj, if new values
   if( new_x && x != NULL )
   {
      env->SetNumberArrayRegion(xj, 0, n, const_cast<Number*>(x));
   }

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
   {
      return false;
   }

   /* Copy from Java to native value */
   if( jac_g != NULL )
   {
      env->GetNumberArrayRegion(jac_gj, 0, nele_jac, jac_g);
   }

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
         for( Index i = 0; i < nele_jac; ++i )
         {
            iRow[i] = (Index) tmp[i];
         }

         env->GetIntArrayRegion(jColj, 0, nele_jac, tmp);
         for( Index i = 0; i < nele_jac; ++i )
         {
            jCol[i] = (Index) tmp[i];
         }

         delete[] tmp;
      }
   }

   return true;
}

bool Jipopt::eval_h(
   Index         n,
   const Number* x,
   bool          new_x,
   Number        obj_factor,
   Index         m,
   const Number* lambda,
   bool          new_lambda,
   Index         nele_hess,
   Index*        iRow,
   Index*        jCol,
   Number*       hess)
{
   /* Copy the native double x to the Java double array xj, if new values */
   if( new_x && x != NULL )
   {
      env->SetNumberArrayRegion(xj, 0, n, const_cast<Number*>(x));
   }

   /* Copy the native double lambda to the Java double array lambdaj, if new values */
   if( new_lambda && lambda != NULL )
   {
      env->SetNumberArrayRegion(mult_gj, 0, m, const_cast<Number*>(lambda));
   }

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
   {
      return false;
   }

   /* Copy from Java to native value */
   if( hess != NULL )
   {
      env->GetNumberArrayRegion(hessj, 0, nele_hess, hess);
   }

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
         for( Index i = 0; i < nele_hess; ++i )
         {
            iRow[i] = (Index) tmp[i];
         }

         env->GetIntArrayRegion(jColj, 0, nele_hess, tmp);
         for( Index i = 0; i < nele_hess; ++i )
         {
            jCol[i] = (Index) tmp[i];
         }

         delete[] tmp;
      }
   }

   return true;
}

void Jipopt::finalize_solution(
   SolverReturn               /*status*/,
   Index                      n,
   const Number*              x,
   const Number*              z_L,
   const Number*              z_U,
   Index                      m,
   const Number*              g,
   const Number*              lambda,
   Number                     obj_value,
   const IpoptData*           /*ip_data*/,
   IpoptCalculatedQuantities* /*ip_cq*/
)
{
   /* Copy the native arrays to Java double arrays */

   if( x != NULL )
   {
      env->SetNumberArrayRegion(xj, 0, n, const_cast<Number*>(x));
   }

   if( z_L != NULL )
   {
      env->SetNumberArrayRegion(mult_x_Lj, 0, n, const_cast<Number*>(z_L));
   }

   if( z_U != NULL )
   {
      env->SetNumberArrayRegion(mult_x_Uj, 0, n, const_cast<Number*>(z_U));
   }

   if( g != NULL )
   {
      env->SetNumberArrayRegion(gj, 0, m, const_cast<Number*>(g));
   }

   if( lambda != NULL )
   {
      env->SetNumberArrayRegion(mult_gj, 0, m, const_cast<Number*>(lambda));
   }

   env->GetNumberArrayRegion(fj, 0, 1, &obj_value);
}

/* Intermediate Callback method for the user. */
bool Jipopt::intermediate_callback(
   AlgorithmMode              mode,
   Index                      iter,
   Number                     obj_value,
   Number                     inf_pr,
   Number                     inf_du,
   Number                     mu,
   Number                     d_norm,
   Number                     regularization_size,
   Number                     alpha_du,
   Number                     alpha_pr,
   Index                      ls_trials,
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq
)
{
   return env->CallBooleanMethod(solver, intermediate_callback_, (int)mode, iter, obj_value, inf_pr, inf_du, mu, d_norm, regularization_size, alpha_du, alpha_pr, ls_trials, ip_data, ip_cq);
}

bool Jipopt::get_scaling_parameters(
   Number& obj_scaling,
   bool&   use_x_scaling,
   Index   n,
   Number* x_scaling,
   bool&   use_g_scaling,
   Index   m,
   Number* g_scaling)
{
   if( !using_scaling_parameters )
   {
      return false;
   }

   jnumberArray obj_scaling_j = env->NewNumberArray(1);
   jnumberArray x_scaling_j   = env->NewNumberArray(n);
   jnumberArray g_scaling_j   = env->NewNumberArray(m);

   jbooleanArray use_x_g_scaling_j = env->NewBooleanArray(2);

   env->CallBooleanMethod(solver, get_scaling_parameters_, obj_scaling_j, n, x_scaling_j, m, g_scaling_j, use_x_g_scaling_j);

   jboolean* use_x_g_scaling = env->GetBooleanArrayElements(use_x_g_scaling_j, 0);

   /* Copy from Java to native value */
   env->GetNumberArrayRegion(obj_scaling_j, 0, 1, &obj_scaling);

   /* Copy from Java to native value */
   if( use_x_g_scaling[0] )
   {
      env->GetNumberArrayRegion(x_scaling_j, 0, n, x_scaling);
      use_x_scaling = true;
   }
   else
   {
      use_x_scaling = false;
   }

   /* Copy from Java to native value */
   if( use_x_g_scaling[1] )
   {
      env->GetNumberArrayRegion(g_scaling_j, 0, n, g_scaling);
      use_g_scaling = true;
   }
   else
   {
      use_g_scaling = false;
   }

   env->ReleaseBooleanArrayElements(use_x_g_scaling_j, use_x_g_scaling, 0);

   return true;
}

Index Jipopt::get_number_of_nonlinear_variables()
{
   if( using_LBFGS )
   {
      return env->CallIntMethod(solver, get_number_of_nonlinear_variables_);
   }

   return -1;
}

bool Jipopt::get_list_of_nonlinear_variables(
   Index  num_nonlin_vars,
   Index* pos_nonlin_vars)
{
   if( !using_LBFGS )
   {
      return false;
   }

   jintArray pos_nonlin_vars_j = env->NewIntArray(num_nonlin_vars);

   if( !env->CallBooleanMethod(solver, get_list_of_nonlinear_variables_, num_nonlin_vars, pos_nonlin_vars_j) )
   {
      return false;
   }

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
         for( Index i = 0; i < num_nonlin_vars; ++i )
         {
            pos_nonlin_vars[i] = (Index) tmp[i];
         }

         delete[] tmp;
      }
   }

   return true;
}

extern "C"
{

   JNIEXPORT jlong JNICALL Java_org_coinor_Ipopt_CreateIpoptProblem(
      JNIEnv* env,
      jobject obj_this,
      jint    n,
      jint    m,
      jint    nele_jac,
      jint    nele_hess,
      jint    index_style)
   {
      /* create the smart pointer to the Ipopt problem */
      SmartPtr<Jipopt>* pproblem = new SmartPtr<Jipopt>;

      /* create the IpoptProblem */
      *pproblem = new Jipopt(env, obj_this, n, m, nele_jac, nele_hess, index_style);

      /* return the smart pointer to our class */
      return (jlong) pproblem;
   }

   JNIEXPORT jint JNICALL Java_org_coinor_Ipopt_OptimizeTNLP(
      JNIEnv*      env,
      jobject      obj_this,
      jlong        pipopt,
      jnumberArray xj,
      jnumberArray gj,
      jnumberArray obj_valj,
      jnumberArray mult_gj,
      jnumberArray mult_x_Lj,
      jnumberArray mult_x_Uj,
      jnumberArray callback_grad_f,
      jnumberArray callback_jac_g,
      jnumberArray callback_hess)
   {
      Jipopt* problem = GetRawPtr(*(SmartPtr<Jipopt>*) pipopt);

      problem->env       = env;
      problem->solver    = obj_this;

      problem->xj        = xj;
      problem->gj        = gj;
      problem->fj        = obj_valj;
      problem->mult_gj   = mult_gj;
      problem->mult_x_Lj = mult_x_Lj;
      problem->mult_x_Uj = mult_x_Uj;

      problem->grad_fj   = callback_grad_f;
      problem->jac_gj    = callback_jac_g;
      problem->hessj     = callback_hess;

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

   JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_GetCurrIterate(
      JNIEnv*      env,
      jobject      /* obj_this */,
      jlong        pipopt,
      jlong        jip_data,
      jlong        jip_cq,
      jboolean     jscaled,
      jint         jn,
      jnumberArray jx,
      jnumberArray jz_L,
      jnumberArray jz_U,
      jint         jm,
      jnumberArray jg,
      jnumberArray jlambda)
   {
      Jipopt* problem = GetRawPtr(*(SmartPtr<Jipopt>*) pipopt);
      IpoptData* ip_data = (IpoptData*)jip_data;
      IpoptCalculatedQuantities* ip_cq = (IpoptCalculatedQuantities*)jip_cq;

      Index n = jn;
      Index m = jm;

      Number* x = NULL;
      Number* z_L = NULL;
      Number* z_U = NULL;
      Number* g = NULL;
      Number* lambda = NULL;

      if( jx != NULL )
      {
         x = new Number[n];
      }
      if( jz_L != NULL )
      {
         z_L = new Number[n];
      }
      if( jz_U != NULL )
      {
         z_U = new Number[n];
      }
      if( jg != NULL )
      {
         g = new Number[m];
      }
      if( jlambda != NULL )
      {
         lambda = new Number[m];
      }

      bool ok = problem->get_curr_iterate(ip_data, ip_cq, jscaled, n, x, z_L, z_U, m, g, lambda);
      if( ok )
      {
         if( jx != NULL )
         {
            env->SetNumberArrayRegion(jx, 0, n, const_cast<Number*>(x));
         }
         if( jz_L != NULL )
         {
            env->SetNumberArrayRegion(jz_L, 0, n, const_cast<Number*>(z_L));
         }
         if( jz_U != NULL )
         {
            env->SetNumberArrayRegion(jz_U, 0, n, const_cast<Number*>(z_U));
         }
         if( jg != NULL )
         {
            env->SetNumberArrayRegion(jg, 0, m, const_cast<Number*>(g));
         }
         if( jlambda != NULL )
         {
            env->SetNumberArrayRegion(jlambda, 0, m, const_cast<Number*>(lambda));
         }
      }

      delete[] lambda;
      delete[] g;
      delete[] z_U;
      delete[] z_L;
      delete[] x;

      return ok;
   }

   JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_GetCurrViolations(
      JNIEnv*      env,
      jobject      /* obj_this */,
      jlong        pipopt,
      jlong        jip_data,
      jlong        jip_cq,
      jboolean     jscaled,
      jint         jn,
      jnumberArray jx_L_violation,
      jnumberArray jx_U_violation,
      jnumberArray jcompl_x_L,
      jnumberArray jcompl_x_U,
      jnumberArray jgrad_lag_x,
      jint         jm,
      jnumberArray jnlp_constraint_violation,
      jnumberArray jcompl_g)
   {
      Jipopt* problem = GetRawPtr(*(SmartPtr<Jipopt>*) pipopt);
      IpoptData* ip_data = (IpoptData*)jip_data;
      IpoptCalculatedQuantities* ip_cq = (IpoptCalculatedQuantities*)jip_cq;

      Index n = jn;
      Index m = jm;

      Number* x_L_violation = NULL;
      Number* x_U_violation = NULL;
      Number* compl_x_L = NULL;
      Number* compl_x_U = NULL;
      Number* grad_lag_x = NULL;
      Number* nlp_constraint_violation = NULL;
      Number* compl_g = NULL;

      if( jx_L_violation != NULL )
      {
         x_L_violation = new Number[n];
      }
      if( jx_U_violation != NULL )
      {
         x_U_violation = new Number[n];
      }
      if( jcompl_x_L != NULL )
      {
         compl_x_L = new Number[n];
      }
      if( jcompl_x_U != NULL )
      {
         compl_x_U = new Number[n];
      }
      if( jgrad_lag_x != NULL )
      {
         grad_lag_x = new Number[n];
      }
      if( jnlp_constraint_violation != NULL )
      {
         nlp_constraint_violation = new Number[m];
      }
      if( jcompl_g != NULL )
      {
         compl_g = new Number[m];
      }

      bool ok = problem->get_curr_violations(ip_data, ip_cq, jscaled, n, x_L_violation, x_U_violation, compl_x_L, compl_x_U, grad_lag_x, m, nlp_constraint_violation, compl_g);
      if( ok )
      {
         if( jx_L_violation != NULL )
         {
            env->SetNumberArrayRegion(jx_L_violation, 0, n, const_cast<Number*>(x_L_violation));
         }
         if( jx_U_violation != NULL )
         {
            env->SetNumberArrayRegion(jx_U_violation, 0, n, const_cast<Number*>(x_U_violation));
         }
         if( jcompl_x_L != NULL )
         {
            env->SetNumberArrayRegion(jcompl_x_L, 0, n, const_cast<Number*>(compl_x_L));
         }
         if( jcompl_x_U != NULL )
         {
            env->SetNumberArrayRegion(jcompl_x_U, 0, n, const_cast<Number*>(compl_x_U));
         }
         if( jgrad_lag_x != NULL )
         {
            env->SetNumberArrayRegion(jgrad_lag_x, 0, n, const_cast<Number*>(grad_lag_x));
         }
         if( jnlp_constraint_violation != NULL )
         {
            env->SetNumberArrayRegion(jnlp_constraint_violation, 0, m, const_cast<Number*>(nlp_constraint_violation));
         }
         if( jcompl_g != NULL )
         {
            env->SetNumberArrayRegion(jcompl_g, 0, m, const_cast<Number*>(compl_g));
         }
      }

      delete[] compl_g;
      delete[] nlp_constraint_violation;
      delete[] grad_lag_x;
      delete[] compl_x_U;
      delete[] compl_x_L;
      delete[] x_U_violation;
      delete[] x_L_violation;

      return ok;
   }

   JNIEXPORT void JNICALL Java_org_coinor_Ipopt_FreeIpoptProblem(
      JNIEnv* /*env*/,
      jobject /*obj_this*/,
      jlong   pipopt
   )
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
      JNIEnv*  env,
      jobject  /*obj_this*/,
      jlong    pipopt,
      jstring  jparname,
      jint     jparvalue
   )
   {
      Jipopt* problem = GetRawPtr(*(SmartPtr<Jipopt>*) pipopt);

      const char* pparameterName = env->GetStringUTFChars(jparname, 0);
      string parameterName = pparameterName;

      // Try to apply the integer option
      jboolean ret = problem->application->Options()->SetIntegerValue(parameterName, jparvalue);

      env->ReleaseStringUTFChars(jparname, pparameterName);

      return ret;
   }

   JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_AddIpoptNumOption(
      JNIEnv* env,
      jobject /*obj_this*/,
      jlong   pipopt,
      jstring jparname,
      jnumber jparvalue
   )
   {
      Jipopt* problem = GetRawPtr(*(SmartPtr<Jipopt>*) pipopt);

      const char* pparameterName = env->GetStringUTFChars(jparname, 0);
      string parameterName = pparameterName;

      // Try to set the real option
      jboolean ret = problem->application->Options()->SetNumericValue(parameterName, jparvalue);

      env->ReleaseStringUTFChars(jparname, pparameterName);

      return ret;
   }

   JNIEXPORT jboolean JNICALL Java_org_coinor_Ipopt_AddIpoptStrOption(
      JNIEnv* env,
      jobject /*obj_this*/,
      jlong   pipopt,
      jstring jparname,
      jstring jparvalue
   )
   {
      Jipopt* problem = GetRawPtr(*(SmartPtr<Jipopt>*) pipopt);

      const char* pparameterName = env->GetStringUTFChars(jparname, NULL);
      string parameterName = pparameterName;
      const char* pparameterValue = env->GetStringUTFChars(jparvalue, NULL);
      string parameterValue = pparameterValue;

      // parameterValue has been changed to LowerCase in Java!
      if( parameterName == "hessian_approximation" && parameterValue == "limited-memory" )
      {
         problem->using_LBFGS = true;
      }
      else if( parameterName == "nlp_scaling_method" && parameterValue == "user-scaling" )
      {
         problem->using_scaling_parameters = true;
      }

      // Try to apply the string option
      jboolean ret = problem->application->Options()->SetStringValue(parameterName, parameterValue);

      env->ReleaseStringUTFChars(jparname, pparameterName);
      env->ReleaseStringUTFChars(jparname, pparameterValue);

      return ret;
   }

} // extern "C"
