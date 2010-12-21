# Copyright (C) 2009, International Business Machines
#
# This file is part of the Ipopt open source package, published under
# the Eclipse Public License
#
# Author:  Andreas Waechter         IBM       2009-04-13

# This is an example of bad modeling.  The formulation can be improved
# in several way, as shown in the *-fix.mod files.

param n := 50;

var x{1..n} >= 0, := 1;
var p{1..n} >= 0, <= 1, :=0.1;

minimize obj: 
  sum{i in 1..n} x[i]
  ;

subject to constr1: prod{i in 1..n} p[i] >= 0.1;

subject to constr2{i in 1..n}: x[i]/p[i] >= i/(10*n);
