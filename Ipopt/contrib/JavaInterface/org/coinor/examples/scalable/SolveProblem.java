/** Copyright (C) 2007 VRTech Industrial Technologies - www.vrtech.com.br.
 * Copyright (C) 2007 Tong Kewei, Beihang University, - www.buaa.edu.cn.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 * 
 * $Id$
 */

package org.coinor.examples.scalable;

import java.util.HashMap;
import org.coinor.Ipopt;

/**
 * Class for running several different Scalable problems.
 *
 * @author Rafael de Pelegrini Soares, Tong Kewei
 */
public class SolveProblem
{
    public static void main(String[] args)
    {
        HashMap<String, Scalable> list = new HashMap<String, Scalable>();
        
        // adding all problems here
        list.put("LukVlE1", new LuksanVlcek1("LukVlE1",  0.0, 0.0));//E means equal
        list.put("LukVlI1", new LuksanVlcek1("LukVlI1", -1.0, 0.0));//I means inequal
        
        if( args.length < 2 )
        {
            System.out.println("Usage: ProblemName N\n");
            System.out.println("  - N is a positive parameter determining problem size");
            System.out.println("  - ProblemName is one of:");
            // list all problems
            for( Scalable s : list.values() )
                System.out.println("       " + s);
            
            return;
        }

        String problem = args[0];
        int n = Integer.parseInt(args[1]);
        
        System.out.println("Solving problem " + problem + " for N=" + n);
        
        Scalable p = list.get(problem);
        if( p == null )
        {
            System.out.println("Problem not found!");
            return;
        }
        
        if( !p.initialize(n) )
            return;
        
        p.create();
        
        p.OptimizeNLP();
        
        switch( p.getStatus() )
        {
            case Ipopt.SOLVE_SUCCEEDED:
            case Ipopt.ACCEPTABLE_LEVEL:
                System.out.println("Solution found.");
                break;
            default:
                System.out.println("** Could not solve problem " + problem + " for N=" + n + ", status: " + p.getStatus());
        }
    }    
}
