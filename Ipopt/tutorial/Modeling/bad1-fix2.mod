# Copyright (C) 2009, International Business Machines
#
# This file is part of the Ipopt open source package, published under
# the Eclipse Public License
#
# Author:  Andreas Waechter         IBM       2009-04-13

# This is an example of bad modeling.
#
# In this second step of the reformulation, the second set of
# nonlinear constraints is replaced by some linear constraints.  In
# this way, we ended up with a nice convex opitmization problem, even
# though the original one was very nonconvex, and Ipopt, in some
# cases, was not even able to find a feasible point!

param n := 50;

var x{1..n} >= 0, := 1;
var p{1..n} >= 0, <= 1, :=0.1;

minimize obj: 
  sum{i in 1..n} x[i]
  ;

# subject to constr1: prod{i in 1..n} p[i] >= 0.1;
subject to constr1: sum{i in 1..n} log(p[i]) >= log(0.1);

# subject to constr2{i in 1..n}: x[i]/p[i] >= i/(10*n);
subject to constr2{i in 1..n}: x[i] >= i/(10*n)*p[i];
