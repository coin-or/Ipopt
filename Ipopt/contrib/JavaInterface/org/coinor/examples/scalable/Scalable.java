/* 
 * Copyright (C) 2007 VRTech Industrial Technologies - www.vrtech.com.br.
 * All Rights Reserved.
 * This code is published under the Common Public License.
 * 
 * $Id$
 * Authors: Rafael de Pelegrini Soares
 *
 * 
 * Copyright (C) 2007 Tong Kewei, Beihang University, - www.buaa.edu.cn.
 * All Rights Reserved.
 * This code is published under the Common Public License.
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

import org.coinor.Ipopt;
/**
 * Abstract class for the scalable problems.
 * 
 * Implementations should derive from this class and
 * implement {@link #initialize(int)} where all problem size,
 * bounds, and inital guess should be initialized.
 * Besides the initialization, the abstract functions for evaluation
 * of objective, bounds, etc need to be implemented.
 * 
 * @author Rafael de Pelegrini Soares, Tong Kewei
 */
public abstract class Scalable   extends Ipopt {
	
	// Problem sizes
	int n, m, nnz_jac_g, nnz_h_lag;

	// The bounds
	double x_l[], x_u[];
	double g_l[], g_u[];
	
	// the index style
	int index_style;
	
	// The initial guess and solution
	double x[];

	private String name;

	protected double gl;

	protected double gu;
	
	/**
	 * @param name
	 * @param gl the constraint lower bound value for all elements of g
	 * @param gu the constraint upper bound value for all elements of g
	 */
	public Scalable(String name, double gl, double gu) {
		this.name = name;
		this.gl = gl;
		this.gu = gu;
	}
	
	public String toString(){
		return name;
	}
	
	/**
	 * In this function all problem sizes, bounds and initial guess
	 * should be initialized.
	 * 
	 * @param n the problem size
	 * 
	 * @return true if the given size is valid for this problem
	 */
	abstract public boolean initialize(int n);
	
	/**
	 * Creates the problem based on the already computed problem sizes and bounds.
	 */
	public void create(){
            
            super.create(n,  m, nnz_jac_g, nnz_h_lag, index_style);
	}
	
	public double[] getInitialGuess(){
		return x;
	}
        
        //add this to show result
        public void print(double[]x,String str){
            System.out.println(str);
            for(int i=0;i<x.length;i++){
                System.out.println(x[i]+"\t");
            }
            System.out.println();
        }
}


