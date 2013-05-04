/** Copyright (C) 2007 VRTech Industrial Technologies - www.vrtech.com.br.
 * Copyright (C) 2007 Tong Kewei, Beihang University, - www.buaa.edu.cn.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 * 
 * $Id$
 */

package org.coinor;

import java.io.File;

/** A Java Native Interface for the Ipopt optimization solver.
 * <p>
 * Ipopt is a solver for large scale nonlinear optimization problems (NLP).
 * <p>
 * The Java Native Interface (JNI) is a programming framework that allows
 * Java code running in the Java Virtual Machine (JVM) to call and be
 * called by native applications (programs specific to a hardware and
 * operating system platform) and libraries written in other languages,
 * such as C and C++.
 * <p>
 * This class is a JNI hook around the C++ interface of Ipopt, as a consequence
 * it will need a nativelly compiled DLL to run.
 * For more details about Ipopt
 * <a href="https://projects.coin-or.org/Ipopt">click here</a>.
 * <p>
 * The user should subclass this class and implement the abstract methods.
 * At some point before solving the problem the
 * {@link #create(int, int, int, int, int)}
 * function should be called.
 * For simple cases you can call this function in the constructor of your class.
 * <p>
 * Once the problem was created, {@link #OptimizeNLP()} will solve the problem.
 * Objects of this class can be reused to solve different problems, in other words,
 * {@link #create(int, int, int, int, int)}
 * and {@link #OptimizeNLP()} can be called multiple times.
 * <p>
 * Programmers should, for efficiency, call {@link #dispose()} when finished using a
 * Ipopt object, otherwise the nativelly allocated memory will be disposed of only
 * when the JVM call {@link #finalize()} on it.
 * 
 * @author Rafael de Pelegrini Soares
 * @author Edson C. do Valle
 * @author Tong Kewei, BeiHang University
 */
public abstract class Ipopt
{
	/** Native function should not be used directly */
	private native boolean AddIpoptIntOption(long ipopt, String keyword, int val);

	/** Native function should not be used directly */
	private native boolean AddIpoptNumOption(long ipopt, String keyword, double val);

	/** Native function should not be used directly */
	private native boolean AddIpoptStrOption(long ipopt, String keyword, String val);

	/** Native function should not be used directly */
	private native long CreateIpoptProblem(int n,int m, 
			int nele_jac, int nele_hess, int index_style);

	/** Native function should not be used directly */
	private native void FreeIpoptProblem(long ipopt);

	/** Native function should not be used directly */
	private native int OptimizeTNLP(long ipopt,
			double x[], double g[],
			double obj_val[], double mult_g[], double mult_x_L[], double mult_x_U[],
			double callback_grad_f[], double callback_jac_g[], double callback_hess[]);


	/** The default DLL name of the native implementation (without any platform dependent prefixes or suffixes) */
	public static final String DLLNAME = "jipopt";
	/** The relative path where the native DLL is found */
	public static final String DLLPATH = "lib";

	/** Use C index style for iRow and jCol vectors */
	public final static int C_STYLE = 0;
	/** Use FORTRAN index style for iRow and jCol vectors */
	public final static int FORTRAN_STYLE = 1;

	/** The possible Ipopt status return codes: should be kept in sync with Ipopt return codes */
	public final static int SOLVE_SUCCEEDED = 0;
	public final static int ACCEPTABLE_LEVEL = 1;
	public final static int INFEASIBLE_PROBLEM = 2;
	public final static int SEARCH_DIRECTION_TOO_SMALL = 3;
	public final static int DIVERGING_ITERATES = 4;
	public final static int USER_REQUESTED_STOP = 5;
	public final static int ITERATION_EXCEEDED = -1;
	public final static int RESTORATION_FAILED = -2;
	public final static int ERROR_IN_STEP_COMPUTATION = -3;
	public final static int CPUTIME_EXCEEDED = -4;
	public final static int NOT_ENOUGH_DEGREES_OF_FRE = -10;
	public final static int INVALID_PROBLEM_DEFINITION = -11;
	public final static int INVALID_OPTION = -12;
	public final static int INVALID_NUMBER_DETECTED = -13;
	public final static int UNRECOVERABLE_EXCEPTION = -100;
	public final static int NON_IPOPT_EXCEPTION = -101;
	public final static int INSUFFICIENT_MEMORY = -102;
	public final static int INTERNAL_ERROR = -199;

	/** Pointer to the native optimization object */
	private long ipopt;

	/// Callback arguments
	private double callback_grad_f[];
	private double callback_jac_g[];
	private double callback_hess[];

	/** Final value of variable values */
	private double x[];
	
	/** Final value of objective function */
	private double obj_val[] = {0.0};

	/** Values of constraint at final point */
	private double g[];

	/** Final multipliers for lower variable bounds */
	private double mult_x_L[];

	/** Final multipliers for upper variable bounds */
	private double mult_x_U[];

	/** Final multipliers for constraints */
	private double mult_g[];

	/** Status returned by the solver*/
	private int status = INVALID_PROBLEM_DEFINITION;

	/** Creates a new NLP Solver using {@value #DLLPATH} as path and {@value #DLLNAME} as the DLL name.
	 * 
	 * @see #Ipopt(String, String)
	 */
	public Ipopt()
	{
		this(DLLPATH, DLLNAME);
	}

	/** Creates a NLP Solver for the given DLL file.
	 * The given file must implement the native interface required by this class.
	 * 
	 * @param path the path where the DLL is found.
	 * @param DLL the name of the DLL (without the extension or any platform dependent prefix).
	 * 
	 * @see #Ipopt()
	 */
	public Ipopt(String path, String DLL)
	{
		// Loads the library
		File file = new File(path, System.mapLibraryName(DLL));
		System.load(file.getAbsolutePath());
	}

	/** Callback function for the variable bounds and constraint sides. */  
	abstract protected boolean get_bounds_info(int n, double[] x_l, double[] x_u,
			int m, double[] g_l, double[] g_u);
	
	/** Callback function for retrieving a starting point. */
	abstract protected boolean get_starting_point(int n, boolean init_x, double[] x,
			boolean init_z, double[] z_L, double[] z_U,
			int m, boolean init_lambda, double[] lambda);
	
	/** Callback function for the objective function. */
	abstract protected boolean eval_f(int n, double[] x, boolean new_x, double[] obj_value);
	
	/** Callback function for the objective function gradient. */
	abstract protected boolean eval_grad_f(int n, double[] x, boolean new_x, double[] grad_f);
	
	/** Callback function for the constraints. */
	abstract protected boolean eval_g(int n, double[] x, boolean new_x, int m, double[] g);
	
	/** Callback function for the constraints Jacobian. */
	abstract protected boolean eval_jac_g(int n, double[] x, boolean new_x,
			int m, int nele_jac, int[] iRow, int[] jCol, double[] values);
	
	/** Callback function for the hessian. */
	abstract protected boolean eval_h(int n, double[] x, boolean new_x, double obj_factor,
			int m, double[] lambda, boolean new_lambda,
			int nele_hess, int[] iRow, int[] jCol,
			double[] values);

	/** Dispose of the natively allocated memory.
	 * Programmers should, for efficiency, call the dispose method when finished
	 * using a Ipopt object.
	 * <p>
	 * An JIpopt object can be reused to solve different problems by calling again
	 * {@link #create(int, int, int, int, int)}.
	 * In this case, you should call the dispose method only when you
	 * finished with the object and it is not needed anymore.
	 */
	public void dispose()
	{
		// dispose the native implementation
		if( ipopt != 0 )
		{
			FreeIpoptProblem(ipopt);
			ipopt = 0;
		}
	}
	
	protected void finalize() throws Throwable
	{
		dispose();
	}

	/** Create a new problem. the use is the same as get_nlp_info, change the name for clarity in java.
	 * 
	 * @param n the number of variables in the problem.
	 * @param m the number of constraints in the problem.
	 * @param nele_jac the number of nonzero entries in the Jacobian.
	 * @param nele_hess the number of nonzero entries in the Hessian.
	 * @param index_style the numbering style used for row/col entries in the sparse matrix format (C_STYLE or FORTRAN_STYLE).
	 *
	 * @return true on success, otherwise false
	 */
	public boolean create(int n, int m,  int nele_jac, int nele_hess, int index_style)
	{
		// delete any previously created native memory
		dispose();

		x = new double[n];
		// allocate the callback arguments
		callback_grad_f = new double[n];
		callback_jac_g = new double[nele_jac];
		callback_hess = new double[nele_hess];

		// the multiplier
		mult_x_U = new double[n];
		mult_x_L = new double[n];
		g = new double[m];
		mult_g = new double[m];

		//	Create the optimization problem and return a pointer to it
		ipopt = CreateIpoptProblem(n, m,  nele_jac, nele_hess, index_style);

		//System.out.println("Finish Java Obj");
		return ipopt == 0 ? false : true;
	}

	/** Function for setting an integer option.
	 * <p>
	 * For a list of valid keywords check the Ipopt documentation.
	 * 
	 * @param keyword the option keyword
	 * @param val the value
	 * @return false if the option could not be set (e.g., if keyword is unknown)
	 */
	public boolean setIntegerOption(String keyword, int val)
	{
		return ipopt == 0 ? false : AddIpoptIntOption(ipopt, keyword, val);
	}

	/** Function for setting a number option.
	 * <p>
	 * For a list of valid keywords check the Ipopt documentation.
	 * 
	 * @param keyword the option keyword
	 * @param val the value
	 * @return false if the option could not be set (e.g., if keyword is unknown)
	 */
	public boolean setNumericOption(String keyword, double val)
	{
		return ipopt == 0 ? false : AddIpoptNumOption(ipopt, keyword, val);
	}

	/** Function for setting a string option.
	 * <p>
	 * For a list of valid keywords check the Ipopt documentation.
	 * 
	 * @param keyword the option keyword
	 * @param val the value
	 * @return false if the option could not be set (e.g., if keyword is unknown)
	 */
	public boolean setStringOption(String keyword, String val)
	{
		return ipopt == 0 ? false : AddIpoptStrOption(ipopt, keyword, val.toLowerCase());
	}

	/** This function actually solve the problem.
	 * <p>
	 * The solve status returned is one of the constant fields of this class,
	 * e.g. SOLVE_SUCCEEDED. For more details about the valid solve status
	 * check the Ipopt documentation.
	 * 
	 * @return the solve status
	 * 
	 * @see #getStatus()
	 */
	public int OptimizeNLP()
	{
		this.status = this.OptimizeTNLP(ipopt,
				x, g, obj_val, mult_g, mult_x_L, mult_x_U,
				callback_grad_f, callback_jac_g, callback_hess);
		
		return this.status;
	}

	/** Gives primal variable values at final point.
	 * @return the primal variable values at the final point.
	 */
	public double[] getVariableValues()
	{
		return x;
	}

	/** Gives objective function value at final point.
	 * @return the final value of the objective function.
	 */
	public double getObjectiveValue()
	{
		return obj_val[0];
	}

	/** Gives Ipopt status of last OptimizeNLP call.
	 * @return the status of the solver.
	 * 
	 * @see #OptimizeNLP()
	 */
	public int getStatus()
	{
		return status;
	}

	/** Gives constraint function values at final point.
	 * @return Returns the final values for the constraints functions. 
	 */
	public double[] getConstraintValues()
	{
		return g;
	}

	/** Gives constraint dual multipliers in final point.
	 * @return Returns the final multipliers for the constraints. 
	 */
	public double[] getConstraintMultipliers()
	{
		return mult_g;
	}

	/** Gives dual multipliers for variable lower bounds in final point.
	 * @return Returns the final multipliers for the variable lower bounds.
	 */
	public double[] getLowerBoundMultipliers()
	{
		return mult_x_L;
	}

	/** Gives dual multipliers for variable upper bounds in final point.
	 * @return Returns the final multipliers for the variable upper bounds.
	 */
	public double[] getUpperBoundMultipliers()
	{
		return mult_x_U;
	}

	/** If you using_scaling_parameters = true, please overload this method, 
	 *
	 * To instruct IPOPT to use scaling values for variables, the first element of use_x_g_scaling should be set.
	 * To instruct IPOPT to use scaling values for constraints, the second element of use_x_g_scaling should be set.
	 *  
	 * @param obj_scaling  double[1] to store a scaling factor for the objective (negative value leads to maximizing the objective function) 
	 * @param n  the number of variables in the problem
	 * @param x_scaling  array to store the scaling factors for the variables
	 * @param m  the number of constraints in the problem
	 * @param g_scaling  array to store the scaling factors for the constraints
	 * @param use_x_g_scaling boolean[2] to store whether scaling factors for variables (1st entry) and constraints (2nd entry) should be used
	 *
	 * @return true on success, otherwise false
	 */
	public boolean get_scaling_parameters(double[] obj_scaling,
			int n, double[] x_scaling,
			int m, double[] g_scaling,
			boolean[] use_x_g_scaling)
	{
		return false;
	}

	/** When LBFGS hessian approximation is used, this method should be overloaded.
	 *
	 * @return number of nonlinear variables, a negative value indicates that all variables are negative
	 */
	public int get_number_of_nonlinear_variables()
	{
		return -1;
	}

	/** When LBFGS hessian approximation is used, this method should be overloaded.
	 *
	 * @param num_nonlin_vars number of nonlinear variables and length of pos_nonlin_vars array
	 * @param pos_nonlin_vars the indices of all nonlinear variables
	 *
	 * @return true on success, otherwise false
	 */
	public boolean get_list_of_nonlinear_variables(int num_nonlin_vars,
			int[] pos_nonlin_vars)
	{
		return false;
	}
}
