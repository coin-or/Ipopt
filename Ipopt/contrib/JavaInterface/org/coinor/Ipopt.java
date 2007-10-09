/* 
 * Copyright (C) 2007 VRTech Industrial Technologies - www.vrtech.com.br.
 * All Rights Reserved.
 * This code is published under the Common Public License.
 * 
 * $Id$
 * Authors: Rafael de Pelegrini Soares
 */

package org.coinor;



import java.io.File;


/**
 * A Java Native Interface for the Ipopt optimization solver.
 * <p>
 * Ipopt is a solver for large scale nonlinear optimization problems (NLP) with 
 * non-linear contraints.
 * <p>
 * The Java Native Interface (JNI) is a programming framework that allows
 * Java code running in the Java Virtual Machine (JVM) to call and be
 * called by native applications (programs specific to a hardware and
 * operating system platform) and libraries written in other languages,
 * such as C and C++.
 * <p>
 * This class is a JNI hook around the C interface of Ipopt, as a consequence
 * it will need a nativelly compiled DLL to run.
 * For more details about Ipopt
 * <a href="https://projects.coin-or.org/Ipopt">click here</a>.
 * <p>
 * The user should subclass this class and implement the abstract methods.
 * At some point before solving the problem the
 * {@link #create(int, double[], double[], int, double[], double[], int, int, int)}
 * function should be called.
 * For simple cases you can call this function in the constructor of your class.
 * <p>
 * Once the problem was created, {@link #solve(double[])} will solve the problem.
 * Objects of this class can be reused to solve different problems, in other words,
 * {@link #create(int, double[], double[], int, double[], double[], int, int, int)}
 * and {@link #solve(double[])} can be called multiple times.
 * <p>
 * Programmers should, for efficiency, call {@link #dispose()} when finished using a
 * Ipopt object, otherwise the nativelly allocated memory will be disposed of only
 * when the JVM call {@link #finalize()} on it.
 * 
 * @author Rafael de Pelegrini Soares, Edson C. do Valle
 * 
 */
public abstract class Ipopt {

	/** Native function should not be used directly */
	private native boolean AddIpoptIntOption(long ipopt, String keyword, int val);

	/** Native function should not be used directly */
	private native boolean AddIpoptNumOption(long ipopt, String keyword, double val);

	/** Native function should not be used directly */
	private native boolean AddIpoptStrOption(long ipopt, String keyword, String val);

	/** Native function should not be used directly */
	private native long CreateIpoptProblem(int n, double ylb[], double yub[],
			int m, double glb[], double gub[],
			int nele_jac, int nele_hess, int index_style);
	
	/** Native function should not be used directly */
	private native boolean OpenIpoptOutputFile(long ipopt, String file_name,
            int print_level);
	
	/** Native function should not be used directly */
	private native void FreeIpoptProblem(long ipopt);

	/** Native function should not be used directly */
	private native int IpoptSolve(long ipopt, double x[], double g[],
			double obj_val[], double mult_g[], double mult_x_L[], double mult_x_U[],
			
			double callback_grad_f[], double callback_jac_g[], double callback_hess[]);
	
	/** Native function should not be used directly */
	private native boolean SetIpoptProblemScaling(long ipopt, double obj_scaling,
			double x_scaling[], double g_scaling[]);
	
	/** The default DLL name of the native implementation (without any platform dependent
	 * prefixes or sufixes) */
	public static final String DLLNAME = "jipopt";
	/** The relative path where the native DLL is found */
	public static final String DLLPATH = "lib";
	
	/** Use C index style for iRow and jCol vectors */
	public final static int C_STYLE = 0;
	/** Use FORTRAN index style for iRow and jCol vectors */
	public final static int FORTRAN_STYLE = 1;
	
	/* The possible return codes: should be kept in sync with Ipopt return codes */
	public final static int SOLVE_SUCCEEDED = 0;
	public final static int ACCEPTABLE_LEVEL = 1;
	public final static int INFEASIBLE_PROBLEM = 2;
	public final static int SEARCH_DIRECTION_TOO_SMALL = 3;
	public final static int DIVERGING_ITERATES = 4;
	public final static int USER_REQUESTED_STOP = 5;
	public final static int ITERATION_EXCEEDED = -1;
	public final static int RESTORATION_FAILED = -2;
	public final static int ERROR_IN_STEP_COMPUTATION = -3;
	public final static int NOT_ENOUGH_DEGREES_OF_FRE = -10;
	public final static int INVALID_PROBLEM_DEFINITION = -11;
	public final static int INVALID_OPTION = -12;
	public final static int INVALID_NUMBER_DETECTED = -13;
	public final static int UNRECOVERABLE_EXCEPTION = -100;
	public final static int NON_IPOPT_EXCEPTION = -101;
	public final static int INSUFFICIENT_MEMORY = -102;
	public final static int INTERNAL_ERROR = -199;

	/* The possible parameter names: should be kept in sync with Ipopt parameters */
	public final static String KEY_TOL = "tol";
	public final static String KEY_COMPL_INF_TOL = "compl_inf_tol"; 
	public final static String KEY_DUAL_INF_TOL = "dual_inf_tol";
	public final static String KEY_CONSTR_VIOL_TOL = "constr_viol_tol";
	public final static String KEY_ACCEPTABLE_TOL = "acceptable_tol"; 
	public final static String KEY_ACCEPTABLE_COMPL_INF_TOL = "acceptable_compl_inf_tol";
	public final static String KEY_ACCEPTABLE_CONSTR_VIOL_TOL = "acceptable_constr_viol_tol";
	public final static String KEY_ACCEPTABLE_DUAL_INF_TOL= "acceptable_dual_inf_tol"; 
	public final static String KEY_BARRIER_TOL_FACTOR = "barrier_tol_factor";
	public final static String KEY_OBJ_SCALING_FACTOR = "obj_scaling_factor";
	public final static String KEY_BOUND_RELAX_FACTOR = "bound_relax_factor"; 
	public final static String KEY_MAX_ITER = "max_iter";
	public final static String KEY_LIMITED_MEMORY_MAX_HISTORY = "limited_memory_max_history";
	public final static String KEY_FILE_PRINT_LEVEL = "file_print_level";
	public final static String KEY_PRINT_LEVEL = "print_level";
	public final static String KEY_MU_STRATEGY = "mu_strategy";
	public final static String KEY_OUTPUT_FILE = "output_file";
	public final static String KEY_DERIVATIVE_TEST_TOL = "derivative_test_tol";
	public final static String KEY_DERIVATIVE_TEST = "derivative_test";
	public final static String KEY_DERIVATIVE_TEST_PRINT_ALL = "derivative_test_print_all";
	public final static String KEY_PRINT_USER_OPTIONS = "print_user_options";
	public final static String KEY_LINEAR_SOLVER = "linear_solver";
	
	/** The hessian approximation, set to "limited-memory" if no hessian is available */
	public final static String KEY_HESSIAN_APPROXIMATION = "hessian_approximation";
	
	/** Pointer to the native optimization object */
	private long ipopt;

	/// Callback arguments
	private double callback_grad_f[];
	private double callback_jac_g[];
	private double callback_hess[];

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

	/**Status returned by the solver*/
	private int status = INVALID_PROBLEM_DEFINITION;

	/**
	 * Creates a new NLP Solver using {@value #DLLPATH} as path and {@value #DLLNAME}
	 * as the DLL name.
	 * 
	 * @see #Ipopt(String, String)
	 */
	public Ipopt(){
		this(DLLPATH, DLLNAME);
	}
	
	/**
	 * Creates a NLP Solver for the given DLL file.
	 * The given file must implement the native interface required by this class.
	 * 
	 * @param path the path where the DLL is found.
	 * @param DLL the name of the DLL (without the extension or any platform dependent prefix).
	 * 
	 * @see #Ipopt()
	 */
	public Ipopt(String path, String DLL){
		// Loads the library
		File file = new File(path, System.mapLibraryName(DLL));
		System.load(file.getAbsolutePath());
	}

	/**
	 * Disposes of the natively allocated memory.
	 * Programmers should, for efficiency, call the dispose method when finished
	 * using a Ipopt object.
	 * <p>
	 * An Ipopt object can be reused to solve different problems by calling again
	 * {@link #create(int, double[], double[], int, double[], double[], int, int, int)}.
	 * In this case, you should call the dispose method only when you
	 * finished with the object and it is not needed anymore.
	 * 
	 * @see #finalize()
	 */
	public void dispose(){
		// dispose the native implementation
		if(ipopt!=0){
			FreeIpoptProblem(ipopt);
			ipopt = 0;
		}
	}

	/**
	 * Disposes of the object once it is no longer referenced (automatically
	 * called by the JVM).
	 * 
	 * @see #dispose()
	 */
	public void finalize(){
		dispose();
	}
	
	/**
	 * Create the NLP problem to be solved.
	 * <p>
	 * This function should be called before {@link #solve(double[])}
	 * and before any of the option configuration methods,
	 * e.g, {@link #addIntOption(String, int)}.
	 * <p>
	 * In simple cases this function in the constructor of the
	 * deriving class.
	 *
	 * @param n number of variables
	 * @param x_L lower bound vector
	 * @param x_U upper bound vector
	 * @param m number of constraints
	 * @param g_L lower bounds on constraints
	 * @param g_U upper bounds on constraints
	 * @param nele_jac number of non-zero elements in constraint Jacobian
	 * @param nele_hess number of non-zero elements in Hessian of Lagrangian
	 * @param index_style indexing style for iRow & jCol ({@link #C_STYLE} or {@link #FORTRAN_STYLE}) 
	 * @return false if the problem could not be created, otherwise true
	 */
	public boolean create(int n, double x_L[], double x_U[],
			int m, double g_L[], double g_U[],
			int nele_jac, int nele_hess, int index_style)
	{
		// delete any previously created native memory
		dispose();

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
		ipopt = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, index_style);
		return ipopt == 0 ? false : true;
	}
	
	/**
	 * Function for adding an integer option.
	 * <p>
	 * The valid keywords are public static members of this class, with names
	 * beginning with <code>KEY_</code>, e.g, {@link #KEY_TOL}.
	 * For more details about the valid options check the Ipopt documentation.
	 * 
	 * @param keyword the option keyword
	 * @param val the value
	 * @return false if the option could not be set (e.g., if keyword is unknown)
	 */
	public boolean addIntOption(String keyword, int val){
		return ipopt==0 ? false : AddIpoptIntOption(ipopt, keyword, val);
	}

	/**
	 * Function for adding a number option.
	 * 
	 * @param keyword the option keyword
	 * @param val the value
	 * @return false if the option could not be set (e.g., if keyword is unknown)
	 * 
	 * @see #addIntOption(String, int)
	 */
	public boolean addNumOption(String keyword, double val){
		return ipopt==0 ? false : AddIpoptNumOption(ipopt, keyword, val);
	}

	/**
	 * Function for adding a string option.
	 * 
	 * @param keyword the option keyword
	 * @param val the value
	 * @return false if the option could not be set (e.g., if keyword is unknown)
	 * 
	 * @see #addIntOption(String, int)
	 */
	public boolean addStrOption(String keyword, String val){
		return ipopt==0 ? false : AddIpoptStrOption(ipopt, keyword, val);
	}
	
	/**
	 * Set the scaling for the optimization problem.
	 * 
	 * @param obj_scaling objective function scaling
	 * @param x_scaling variables scaling
	 * @param g_scaling constraints scaling
	 */
	public void setProblemScaling(double obj_scaling,
			double x_scaling[], double g_scaling[]){
		if(ipopt!=0)
			SetIpoptProblemScaling(ipopt, obj_scaling, x_scaling, g_scaling);
	}
	
	/**
	 * Function for opening an output file for a given name with given printlevel.
	 * @param file_name the file name
	 * @param print_level the print level
	 * @return false, if there was a problem opening the file
	 */
	public boolean openOutputFile(String file_name, int print_level){
		return ipopt==0 ? false : OpenIpoptOutputFile(ipopt, file_name, print_level);
	}

	/**
	 * This function actually solve the problem.
	 * <p>
	 * The solve status returned is one of the constant fields of this class,
	 * e.g. SOLVE_SUCCEEDED. For more details about the valid solve status
	 * check the Ipopt documentation or the <code>ReturnCodes_inc.h<\code>
	 * which is installed in the Ipopt include directory.
	 * 
	 * @param x the start point, and the solution when returns
	 * @return the solve status
	 * 
	 * @see #getStatus()
	 */
	public int solve(double[] x) {
		if(ipopt != 0)
			status = IpoptSolve(ipopt, x, g, obj_val, mult_g, mult_x_L, mult_x_U,
				callback_grad_f, callback_jac_g, callback_hess);
		return status;
	}

	/** Callback function for the objective function. */
	abstract protected boolean eval_f(int n, double []x, boolean new_x, double []obj_value);
	/** Callback function for the objective function gradient */
	abstract protected boolean eval_grad_f(int n, double []x, boolean new_x, double []grad_f);
	/** Callback function for the constraints */
	abstract protected boolean eval_g(int n, double []x, boolean new_x, int m, double []g);
	/** Callback function for the constraints Jacobian */
	abstract protected boolean eval_jac_g(int n, double []x, boolean new_x,
			int m, int nele_jac, int []iRow, int []jCol, double []values);
	/** Callback function for the hessian */
	abstract protected boolean eval_h(int n, double[] x, boolean new_x, double obj_factor,
            int m, double []lambda, boolean new_lambda,
            int nele_hess, int[] iRow, int []jCol,
            double []values);

	/**
	 * @return the final value of the objective function.
	 */
	public double getObjVal() {
		return obj_val[0];
	}

	/**
	 * @return the status of the solver.
	 * 
	 * @see #solve(double[])
	 */
	public int getStatus(){
		return status;
	}

	/**
	 * @return Returns the final multipliers for constraints.
	 */
	public double[] getMultConstraints() {
		return mult_g;
	}

	/**
	 * @return Returns the final multipliers for upper variable bounds.
	 */
	public double[] getMultUpperBounds() {
		return mult_x_U;
	}

	/**
	 * @return Returns the final multipliers for lower variable bounds.
	 */
	public double[] getMultLowerBounds() {
		return mult_x_L;
	}
}
