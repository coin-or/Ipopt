/*
 *
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
import java.util.HashMap;

import org.coinor.Ipopt;
/**
 * Class for running several different Scalable problems.
 *
 * @author Rafael de Pelegrini Soares, Tong Kewei
 */
public class SolveProblem{
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        
        HashMap<String, Scalable> list = new HashMap<String, Scalable>();
        
        // adding all problems here
        list.put("LukVlE1", new LuksanVlcek1("LukVlE1", 0,0));//E means equal
        list.put("LukVlI1", new LuksanVlcek1("LukVlI1", -1,0));//I means inequal
        
        
        /*
        if(args.length < 2){
            System.out.println("Usage: ProblemName N\n");
            System.out.println("  - N is a positive parameter determining problem size");
            System.out.println("  - ProblemName is one of:");
            // list all problems
            for(Scalable s : list.values()){
                System.out.println("       " + s);
            }
            return;
        }
        */
        //I changed this in order to facilitate to execute on IDE such as NetBeans/Eclipse
        args=new String[2];
        args[0]="LukVlI1";//"LukVlE1";//"LukVlI1""
        args[1]="20000";
         
        
        String problem = args[0];
        int n = Integer.parseInt(args[1]);
        
        System.out.println("Solving problem " + problem + " for N=" + n);
        Scalable p = list.get(problem);
        if(p == null){
            System.out.println("Problem not found!");
            return;
        }
        if(!p.initialize(n))
            return;
        p.create();
        
        //p.addStrOption(Ipopt.KEY_HESSIAN_APPROXIMATION,"limited-memory");//LBFGS
        
        double x[] = p.getInitialGuess();
        p.OptimizeNLP();//.solve(x);
        
        switch (p.getStatus()){
            case Ipopt.SOLVE_SUCCEEDED:
            case Ipopt.ACCEPTABLE_LEVEL:
                System.out.println("Solution found.");
                break;
            default:
                System.out.println("** Could not solve problem " + problem + " for N=" + n
                        + " , status:" + p.getStatus());
        }
        
    }
    
    
    
}
