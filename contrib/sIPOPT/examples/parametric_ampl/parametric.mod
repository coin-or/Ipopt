# Copyright 2009, 2011 Hans Pirnay
# All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# Date   : 2010-10-04

# variables
var x1>=0 := 0.15;
var x2>=0 := 0.15;
var x3>=0 := 0.0 ;

# parameters
var eta1;
var eta2;

# model
const1: 6*x1+3*x2+2*x3-eta1=0;
const2: eta2*x1+x2-x3-1=0;

# initial constraints for parameters
consteta1: eta1=nominal_eta1;
consteta2: eta2=nominal_eta2;

# objective
minimize cost: x1^2+x2^2+x3^2;
