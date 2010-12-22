/* 
 * Copyright (C) 2007 VRTech Industrial Technologies - www.vrtech.com.br.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 * 
 * $Id$
 * Authors: Rafael de Pelegrini Soares 
 *
 *
 * Copyright (C) 2007 Tong Kewei, Beihang University, - www.buaa.edu.cn.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 * 
 * $Id$
 *
 * I changed from Rafael de Pelegrini Soares's original code.
 * His codes are originally derived form C version of IPOPT,which has limited functions. 
 * I derived my codes from C++ version of IPOPT, which is much more powerful.  
 * I also fix a bug in Rafael de Pelegrini Soares's code on function setProblemScaling,
 * In his original code the function setProblemScaling has no use to change problem.
 * I added some useful functions in JIpopt, such as get_scaling_parameters or get_number_of_nonlinear_variables
 * or get_list_of_nonlinear_variables. You can add any more functions as you like. Follow my structure it is 
 * very easy.
 *
 * If you have problem or you need me to add another functions, please contact me.
 *
 * Authors: Tong Kewei, E-mail:tongkewei@126.com
 * Beihang University, website: www.buaa.edu.cn
 * Beijing,China.
 * 2007-11-11
 */


package org.coinor.examples.scalable;


/**
 * Implementation of Example 5.1 from "Sparse and Parially Separable
 * Test Problems for Unconstrained and Equality Constrained
 * Optimization" by L. Luksan and J. Vlcek.
 * <p>
 * This code was based on Ipopt file with same name.
 * 
 * @author Rafael de Pelegrini Soares, Tong Kewei
 */
public class LuksanVlcek1 extends Scalable {

	/**
	 * Constructor.
	 * Here, gl and gu are the bounds for the constraints.
	 * The original formulation is obtained by setting
	 * gl and gu to zero. Using gl lower than gu allows the obtain a problem
	 * formulation with inequality constraints.
	 */
	public LuksanVlcek1(String name, double gl, double gu) {
		super(name, gl, gu);
	}

	@Override
	public boolean initialize(int n) {
		if (n<=2) {
			System.out.print("N needs to be at least 3.\n");
			return false;
		}

		// The problem described in LuksanVlcek1.hpp has 4 variables, x[0] through x[3]
		this.n = n;

		m = n-2;

		nnz_jac_g = m*3;

		nnz_h_lag = n + n-1;

		// use the C style numbering of matrix indices (starting at 0)
		index_style = C_STYLE;

		// none of the variables have bounds
		x_l = new double[n];
		x_u = new double[n];
		for (int i=0; i<n; i++) {
			x_l[i] = -1e20;
			x_u[i] =  1e20;
		}

		// Set the bounds for the constraints
		g_l = new double[m];
		g_u = new double[m];
		for (int i=0; i<m; i++) {
			g_l[i] = gl;
			g_u[i] = gu;
		}

		// set the starting point
		x = new double[n];
		for (int i=0; i<n/2; i++) {
			x[2*i] = -1.2;
			x[2*i+1] = 1.;
		}
		if (n != 2*(n/2)) {
			x[n-1] = -1.2;
		}


		return true;
	}
      
     // Callback function for the objective function.  
     protected boolean get_bounds_info(int n, double[] x_l, double[] x_u,
                            int m, double[] g_l, double[] g_u){
         // none of the variables have bounds
		//x_l = new double[n];
		//x_u = new double[n];
		for (int i=0; i<n; i++) {
			x_l[i] = -1e20;
			x_u[i] =  1e20;
		}

		// Set the bounds for the constraints
		//g_l = new double[m];
		//g_u = new double[m];
		for (int i=0; i<m; i++) {
			g_l[i] = gl;
			g_u[i] = gu;
		}
                
                return true;
     }
    
     /** Callback function for the objective function. */
     protected boolean get_starting_point(int n, boolean init_x, double[] x,
                               boolean init_z, double[] z_L, double[] z_U,
                               int m, boolean init_lambda,double[] lambda){
         for (int i=0; i<n/2; i++) {
			x[2*i] = -1.2;
			x[2*i+1] = 1.;
		}
		if (n != 2*(n/2)) {
			x[n-1] = -1.2;
		}
         
         return true;

     }
        
        
	@Override
	protected boolean eval_f(int n, double[] x, boolean new_x,
			double[] obj_value) {
		obj_value[0] = 0.;
		for (int i=0; i<n-1; i++) {
			double a1 = x[i]*x[i]-x[i+1];
			double a2 = x[i] - 1.;
			obj_value[0] += 100.*a1*a1 + a2*a2;
		}

		return true;
	}

	@Override
	protected boolean eval_g(int n, double[] x, boolean new_x, int m, double[] g) {
		for (int i=0; i<n-2; i++) {
			g[i] = 3.*Math.pow(x[i+1],3.) + 2.*x[i+2] - 5.
			+ Math.sin(x[i+1]-x[i+2])*Math.sin(x[i+1]+x[i+2]) + 4.*x[i+1]
			                                                         - x[i]*Math.exp(x[i]-x[i+1]) - 3.;
		}

		return true;
	}

	@Override
	protected boolean eval_grad_f(int n, double[] x, boolean new_x,
			double[] grad_f) {
		grad_f[0] = 0.;
		for (int i=0; i<n-1; i++) {
			grad_f[i] += 400.*x[i]*(x[i]*x[i]-x[i+1]) + 2.*(x[i]-1.);
			grad_f[i+1] = -200.*(x[i]*x[i]-x[i+1]);
		}

		return true;
	}

	@Override
	protected boolean eval_h(int n, double[] x, boolean new_x,
			double obj_factor, int m, double[] lambda, boolean new_lambda,
			int nele_hess, int[] iRow, int[] jCol, double[] values) {
		if (values == null) {
			int ihes=0;
			for (int i=0; i<n; i++) {
				iRow[ihes] = i;
				jCol[ihes] = i;
				ihes++;
				if (i<n-1) {
					iRow[ihes] = i;
					jCol[ihes] = i+1;
					ihes++;
				}
			}
			assert ihes == nele_hess;
		}
		else {                    
			int ihes=0;
			for (int i=0; i<n; i++) {
				// x[i],x[i]
				if (i<n-1) {
					values[ihes] = obj_factor*(2.+400.*(3.*x[i]*x[i]-x[i+1]));
					if (i<n-2) {
						values[ihes] -= lambda[i]*(2.+x[i])*Math.exp(x[i]-x[i+1]);
					}
				}
				else {
					values[ihes] = 0.;
				}
				if (i>0) {
					// x[i+1]x[i+1]
					values[ihes] += obj_factor*200.;
					if (i<n-1) {
						values[ihes] += lambda[i-1]*(18.*x[i]
						                                   - 2.*Math.sin(x[i]-x[i+1])*Math.sin(x[i]+x[i+1])
						                                   + 2.*Math.cos(x[i]-x[i+1])*Math.cos(x[i]+x[i+1])
						                                   - x[i-1]*Math.exp(x[i-1]-x[i]));
					}
				}
				if (i>1) {
					// x[i+2]x[i+2]
					values[ihes] +=
						lambda[i-2]*(- 2.*Math.sin(x[i-1]-x[i])*Math.sin(x[i-1]+x[i])
								- 2.*Math.cos(x[i-1]-x[i])*Math.cos(x[i-1]+x[i]));
				}
				ihes++;

				if (i<n-1) {
					// x[i],x[i+1]
					values[ihes] = obj_factor*(-400.*x[i]);
					if (i<n-2) {
						values[ihes] += lambda[i]*(1.+x[i])*Math.exp(x[i]-x[i+1]);
					}
					/*
			        if (i>0) {
			        // x[i+1],x[i+2]
			        values[ihes] +=
			        lambda[i-1]*(  sin(x[i]-x[i+1])*sin(x[i]+x[i+1])
			           + cos(x[i]-x[i+1])*cos(x[i]+x[i+1])
			           - cos(x[i]-x[i+1])*cos(x[i]+x[i+1])
			           - sin(x[i]-x[i+1])*sin(x[i]+x[i+1])
			        );
			        }
					 */
					ihes++;
				}
			}
			assert ihes == nele_hess;
		}

		return true;
	}

	@Override
	protected boolean eval_jac_g(int n, double[] x, boolean new_x, int m,
			int nele_jac, int[] iRow, int[] jCol, double[] values) {
		if (values == null) {
			// return the structure of the jacobian

			int ijac=0;
			for (int i=0; i<n-2; i++) {
				iRow[ijac] = i;
				jCol[ijac] = i;
				ijac++;
				iRow[ijac] = i;
				jCol[ijac] = i+1;
				ijac++;
				iRow[ijac] = i;
				jCol[ijac] = i+2;
				ijac++;
			}
		}
		else {
			// return the values of the jacobian of the constraints

			int ijac=0;
			for (int i=0; i<n-2; i++) {
				// x[i]
				values[ijac] = -(1.+x[i])*Math.exp(x[i]-x[i+1]);
				ijac++;
				// x[i+1]
				values[ijac] = 9.*x[i+1]*x[i+1]
				                           + Math.cos(x[i+1]-x[i+2])*Math.sin(x[i+1]+x[i+2])
				                           + Math.sin(x[i+1]-x[i+2])*Math.cos(x[i+1]+x[i+2])
				                           + 4. + x[i]*Math.exp(x[i]-x[i+1]);
				ijac++;
				// x[i+2]
				values[ijac] = 2.
				- Math.cos(x[i+1]-x[i+2])*Math.sin(x[i+1]+x[i+2])
				+ Math.sin(x[i+1]-x[i+2])*Math.cos(x[i+1]+x[i+2]);
				ijac++;
			}
		}

		return true;
	}

}
