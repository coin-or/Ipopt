# Copyright (C) 2009, International Business Machines
#
# This file is part of the Ipopt open source package, published under
# the Eclipse Public License
#
# Author:  Andreas Waechter         IBM       2009-04-13

# This is an example of bad modeling.
#
# In this reformulation, the product of the positive variables is
# replaced by the sum of the log terms.  This makes the first set of
# constraints convex, and it also makes the Hessian sparse (it is
# dense in the original formulation)

param n := 50;

var x{1..n} >= 0, := 1;
var p{1..n} >= 0, <= 1, :=0.1;

minimize obj: 
  sum{i in 1..n} x[i]
  ;

# subject to constr1: prod{i in 1..n} p[i] >= 0.1;
subject to constr1: sum{i in 1..n} log(p[i]) >= log(0.1);

subject to constr2{i in 1..n}: x[i]/p[i] >= i/(10*n);
