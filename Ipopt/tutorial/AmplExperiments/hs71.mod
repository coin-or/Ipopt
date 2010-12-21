# Copyright (C) 2009, International Business Machines
#
# This file is part of the Ipopt open source package, published under
# the Eclipse Public License
#
# Author:  Andreas Waechter         IBM       2009-04-03

# This is a model of Example 71 from
#
# Hock, W, and Schittkowski, K,
# Test Examples for Nonlinear Programming Codes,
# Lecture Notes in Economics and Mathematical Systems.
# Springer Verlag, 1981.

##############################################################################

# Definition of the variables with bounds
var x {i in 1..4}, >= 1, <= 5;

# objective function
minimize obj:
  x[1]*x[4]*(x[1] + x[2] + x[3]) + x[3];

# and the constraints
subject to c1:
  x[1]*x[2]*x[3]*x[4] >= 25;

subject to c2:
  x[1]^2+x[2]^2+x[3]^2+x[4]^2 = 40;

# Now we set the starting point:

let x[1] := 1;
let x[2] := 5;
let x[3] := 5;
let x[4] := 1;

