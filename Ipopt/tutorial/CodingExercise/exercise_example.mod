# Copyright (C) 2009, International Business Machines
#
# This file is part of the Ipopt open source package, published under
# the Eclipse Public License
#
# Author:  Andreas Waechter         IBM       2009-04-02
#
# This is the AMPL formulation of the coding example problem from the
# Ipopt tutorial

# Number of variables (this is a scalable formulation)
param n := 8;

# Definition of the variables with bounds
var x {1..n} <= 0, >= -1.5, := -0.5;

# The objective function....
minimize obj: 
  sum{i in 1..n} (x[i]-1)^2;

# ... and the constraints
subject to constr {i in 2..n-1}:
  (x[i]^2+1.5*x[i]-i/n)*cos(x[i+1]) - x[i-1] = 0;
