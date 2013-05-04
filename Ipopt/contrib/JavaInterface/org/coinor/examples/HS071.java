/** Copyright (C) 2007 VRTech Industrial Technologies - www.vrtech.com.br.
 * Copyright (C) 2007 Tong Kewei, Beihang University, - www.buaa.edu.cn.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 * 
 * $Id$
 */

package org.coinor.examples;
import org.coinor.Ipopt;

/** Java example for interfacing with IPOPT.
 * 
 * HS071 implements a Java example of problem 71 of the
 * Hock-Schittkowsky test suite.
 * <p>
 * The optimal solution is
 * x = (1.00000000, 4.74299963, 3.82114998, 1.37940829).
 * <p>
 * This code was based on same problem of the Ipopt distribution.
 *
 * @author Rafael de Pelegrini Soares, Tong Kewei
 */
public class HS071 extends Ipopt
{
    // Problem sizes
    int n;
    int m;
    int nele_jac;
    int nele_hess;
    
    int count_bounds = 0;
    int dcount_start = 0;

    /** Creates a new instance of HS071cpp */
    public HS071()
    {
        /* Number of nonzeros in the Jacobian of the constraints */
        nele_jac = 8;
        
        /* Number of nonzeros in the Hessian of the Lagrangian (lower or
         * upper triangual part only) */
        nele_hess = 10;
        
        /* Number of variables */
        n = 4;
        
        /* Number of constraints */
        m = 2;
        
        /* Index style for the irow/jcol elements */
        int index_style = Ipopt.C_STYLE;
        
        /* create the IpoptProblem */
        create(n, m, nele_jac, nele_hess, index_style);
    }
    
    /** Callback function for variable bounds and constraint sides. */
    protected boolean get_bounds_info(int n, double[] x_L, double[] x_U,
            int m, double[] g_L, double[] g_U)
    {
        assert n == this.n;
        assert m == this.m;
        
    	/* set the values of the variable bounds */
        for( int i = 0; i < n; ++i )
        {
            x_L[i] = 1.0;
            x_U[i] = 5.0;
        }
        
        /* set the values of the constraint bounds */
        g_L[0] = 25.0;  g_U[0] = 2e19;
        g_L[1] = 40.0;  g_U[1] = 40.0;
        
        return true;
    }
    
    /** Callback function for the starting point. */
    protected boolean get_starting_point(int n, boolean init_x, double[] x,
            boolean init_z, double[] z_L, double[] z_U,
            int m, boolean init_lambda,double[] lambda)
    {
    	assert init_z == false;
    	assert init_lambda = false;
    	
    	if( init_x )
    	{
    		x[0] = 1.0;
    		x[1] = 5.0;
    		x[2] = 5.0;
    		x[3] = 1.0;
    	}

        /*
        x[0] = 0.9999999923240762;
        x[1] = 4.742999641809297;
        x[2] = 3.8211499817883072;
        x[3] = 1.3794082897556983;

        z_L[0] = 1.0878712258676539e+00;
        z_L[1] = 0;
        z_L[2] = 0;	
        z_L[3] = 0;
        
        z_U[0] = 0;
        z_U[1] = 0;
        z_U[2] = 0;
        z_U[3] = 0;
        
        lambda[0] = -0.552293195627571;
        lambda[1] = 0.16146777361782;
        */
        
        return true;
    }
    
    protected boolean eval_f(int n, double[] x, boolean new_x, double[] obj_value)
    {
        assert n == this.n;
        
        obj_value[0] = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
        
        return true;
    }
    
    protected boolean eval_grad_f(int n, double[] x, boolean new_x, double[] grad_f)
    {
        assert n == this.n;
        
        grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
        grad_f[1] = x[0] * x[3];
        grad_f[2] = x[0] * x[3] + 1;
        grad_f[3] = x[0] * (x[0] + x[1] + x[2]);        
       
        return true;
    }
    
    protected boolean eval_g(int n, double[] x, boolean new_x, int m, double[] g)
    {
        assert n == this.n;
        assert m == this.m;
        
        g[0] = x[0] * x[1] * x[2] * x[3];
        g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
        
        return true;
    }
    
    protected boolean eval_jac_g(int n, double[] x, boolean new_x,
            int m, int nele_jac, int[] iRow, int[] jCol, double[] values)
    {
        assert n == this.n;
        assert m == this.m;
        
        if( values == null )
        {
            /* return the structure of the jacobian */
            
            /* this particular jacobian is dense */
            iRow[0] = 0;  jCol[0] = 0;
            iRow[1] = 0;  jCol[1] = 1;
            iRow[2] = 0;  jCol[2] = 2;
            iRow[3] = 0;  jCol[3] = 3;
            iRow[4] = 1;  jCol[4] = 0;
            iRow[5] = 1;  jCol[5] = 1;
            iRow[6] = 1;  jCol[6] = 2;
            iRow[7] = 1;  jCol[7] = 3;
        }
        else
        {
            /* return the values of the jacobian of the constraints */
            
            values[0] = x[1]*x[2]*x[3]; /* 0,0 */
            values[1] = x[0]*x[2]*x[3]; /* 0,1 */
            values[2] = x[0]*x[1]*x[3]; /* 0,2 */
            values[3] = x[0]*x[1]*x[2]; /* 0,3 */
            
            values[4] = 2*x[0];         /* 1,0 */
            values[5] = 2*x[1];         /* 1,1 */
            values[6] = 2*x[2];         /* 1,2 */
            values[7] = 2*x[3];         /* 1,3 */
        }

        return true;
    }
    
    protected boolean eval_h(int n, double[] x, boolean new_x, double obj_factor, int m, double[] lambda, boolean new_lambda, int nele_hess, int[] iRow, int[] jCol, double[] values)
    {
        assert n == this.n;
        assert m == this.m;
        
        int idx = 0; /* nonzero element counter */
        int row = 0; /* row counter for loop */
        int col = 0; /* col counter for loop */
        if( values == null )
        {
        	/* return the structure. This is a symmetric matrix, fill the lower left triangle only. */

        	/* the hessian for this problem is actually dense */
            idx = 0;
            for( row = 0; row < n; ++row )
            {
                for( col = 0; col <= row; ++col)
                {
                    iRow[idx] = row;
                    jCol[idx] = col;
                    idx++;
                }
            }
            
            assert idx == nele_hess;
            assert nele_hess == this.nele_hess;
        }
        else
        {
        	/* return the values. This is a symmetric matrix, fill the lower left triangle only */

            /* fill the objective portion */
            values[0] = obj_factor * (2*x[3]);               /* 0,0 */
            values[1] = obj_factor * (x[3]);                 /* 1,0 */
            values[2] = 0.0;                                 /* 1,1 */
            values[3] = obj_factor * (x[3]);                 /* 2,0 */
            values[4] = 0.0;                                 /* 2,1 */
            values[5] = 0.0;                                 /* 2,2 */
            values[6] = obj_factor * (2*x[0] + x[1] + x[2]); /* 3,0 */
            values[7] = obj_factor * (x[0]);                 /* 3,1 */
            values[8] = obj_factor * (x[0]);                 /* 3,2 */
            values[9] = 0.0;                                 /* 3,3 */
            
            /* add the portion for the first constraint */
            values[1] += lambda[0] * (x[2] * x[3]);          /* 1,0 */
            values[3] += lambda[0] * (x[1] * x[3]);          /* 2,0 */
            values[4] += lambda[0] * (x[0] * x[3]);          /* 2,1 */
            values[6] += lambda[0] * (x[1] * x[2]);          /* 3,0 */
            values[7] += lambda[0] * (x[0] * x[2]);          /* 3,1 */
            values[8] += lambda[0] * (x[0] * x[1]);          /* 3,2 */
            
            /* add the portion for the second constraint */
            values[0] += lambda[1] * 2.0;                    /* 0,0 */
            values[2] += lambda[1] * 2.0;                    /* 1,1 */
            values[5] += lambda[1] * 2.0;                    /* 2,2 */
            values[9] += lambda[1] * 2.0;                    /* 3,3 */
        }
        
        return true;
    }
    
    public boolean get_scaling_parameters(double[] obj_scaling,
            int n, double[] x_scaling,
            int m, double[] g_scaling,
            boolean[] use_x_g_scaling)
    {
        return false;
    }


    public void print(double[] x, String str)
    {
        System.out.println(str);
        for( int i = 0; i < x.length; ++i )
            System.out.println(x[i]);
        System.out.println();
    }

    /**
     * Main function for running this example.
     */
    public static void main(String []args)
    {
        //Create the problem
        HS071 hs071 = new HS071();
        
        //Set some options
        //hs071.setNumericOption("tol",1E-7);
        //hs071.setStringOption("nlp_scaling_method","user-scaling");
        //hs071.setStringOption("print_options_documentation","yes");
        //hs071.setStringOption("warm_start_init_point","yes");
        //hs071.setNumericOption("warm_start_bound_push",1e-9);
        //hs071.setNumericOption("warm_start_bound_frac",1e-9);
        //hs071.setNumericOption("warm_start_slack_bound_frac",1e-9);
        //hs071.setNumericOption("warm_start_slack_bound_push",1e-9);
        //hs071.setNumericOption("warm_start_mult_bound_push",1e-9);
        
        //Solve the problem
        int status = hs071.OptimizeNLP();        
        
        //Print the solution
        if( status == SOLVE_SUCCEEDED )
        	System.out.println("\n\n*** The problem solved!");
        else
        	System.out.println("\n\n*** The problem was not solved successfully!");
        
        double obj = hs071.getObjectiveValue();
        System.out.println("\nObjective Value = " + obj + "\n");

        double x[] = hs071.getVariableValues();
        hs071.print(x, "Primal Variable Values:");
        
        double constraints[] = hs071.getConstraintValues();
        hs071.print(constraints, "Constraint Values:");

        double MLB[] = hs071.getLowerBoundMultipliers();
        hs071.print(MLB, "Dual Multipliers for Variable Lower Bounds:");
        
        double MUB[] = hs071.getUpperBoundMultipliers();
        hs071.print(MUB, "Dual Multipliers for Variable Upper Bounds:");
        
        double lam[] = hs071.getConstraintMultipliers();
        hs071.print(lam, "Dual Multipliers for Constraints:");
    }
}
