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


# Definition of the variables with starting point and bounds
var x{1..2} >= 0;

# objective function
minimize obj:
  -x[1]^2 + x[2];

# here the constraint.  It is infeasible!
subject to c1:
  x[1]^2 <= -1;
