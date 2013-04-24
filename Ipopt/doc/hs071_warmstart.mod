# tell ampl to use the ipopt executable as a solver
# make sure ipopt is in the path!
option solver ipopt;

# declare the variables and their bounds,
# set notation could be used, but this is straightforward
var x1 >= 1, <= 5;
var x2 >= 1, <= 5;
var x3 >= 1, <= 5;
var x4 >= 1, <= 5;

# specify the objective function
minimize obj:
x1 * x4 * (x1 + x2 + x3) + x3;

# specify the constraints
s.t.
inequality:
x1 * x2 * x3 * x4 >= 25;

equality:
x1^2 + x2^2 + x3^2 +x4^2 = 40;

# specify the starting point
let x1 := 1;
let x2 := 5;
let x3 := 5;
let x4 := 1;

# solve the problem
solve;

# ipopt returns the bound multipliers through
# the suffixes ipopt_zL_out and ipopt_zU_out
# the user does not need to specify them
# print the solution and variable bounds multipliers
display x1;
display x2;
display x3;
display x4;
display x1.ipopt_zL_out;
display x1.ipopt_zU_out;
display x2.ipopt_zL_out;
display x2.ipopt_zU_out;
display x3.ipopt_zL_out;
display x3.ipopt_zU_out;
display x4.ipopt_zL_out;
display x4.ipopt_zU_out;

# define initial conditions for bound multipliers
# to solve new problem

suffix ipopt_zL_in, IN;
suffix ipopt_zU_in, IN;

let x1.ipopt_zL_in := x1.ipopt_zL_out;
let x1.ipopt_zU_in := x1.ipopt_zU_out;
let x2.ipopt_zL_in := x2.ipopt_zL_out;
let x2.ipopt_zU_in := x2.ipopt_zU_out;
let x3.ipopt_zL_in := x3.ipopt_zL_out;
let x3.ipopt_zU_in := x3.ipopt_zU_out;
let x4.ipopt_zL_in := x4.ipopt_zL_out;
let x4.ipopt_zU_in := x4.ipopt_zU_out;
    
# set options for warm-start 
option ipopt_options "warm_start_init_point yes  warm_start_bound_push 1e-6 warm_start_mult_bound_push 1e-6  mu_init 1e-6";

# solve the problem again
solve;
