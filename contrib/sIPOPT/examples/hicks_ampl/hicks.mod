# ===================================================================
# dynamic optimization formulation of the hicks-ray reactor
# model declaration
# victor m zavala  march 2006
# adapted for asNMPC by Hans Pirnay 2009, 2011
# ===================================================================

# define indexes and general variables

param nfe >= 1 integer       ;
param ncp >= 1 integer       ;

# define mathematical model parameters

param time      ;  
param jj        ;
param cf        ;
param alpha     ;
param tf        ;
param k10       ;
param tc        ;
param n         ;
param alpha1    ;
param alpha2    ;
param alpha3    ;
param c_des     ;
param t_des     ;
param u_des     ;
param c_init    ;
param t_init    ;
param u_init    ;
param r1        ;
param r2        ;
param r3        ;
param theta     ;
param point     ;
param slopec    ;
param slopet    ;
param slopeu    ;

# define dimensions for all indexed variables

set fe := 1..nfe ;  # number of finite elements
set cp := 1..ncp ;  # number of collocation points

param a{cp,cp} ;    # collocation matrix
param h{fe}    ;    # finite element length

# define the decision variables

var c {fe,cp}  >= 0           ;  
var t {fe,cp}  >= 0           ;  
var u {fe,cp}  >= 0           ;  

# auxiliary equations

param yc := tc/(jj*cf)         ;
param yf := tf/(jj*cf)         ;

# states first order derivatives
var cdot{i in fe, j in cp} = (1-c[i,j])/theta-k10*exp(-n/t[i,j])*c[i,j]                           ;
var tdot{i in fe, j in cp} = (yf-t[i,j])/theta+k10*exp(-n/t[i,j])*c[i,j]-alpha*u[i,j]*(t[i,j]-yc) ;

#---------------------------------
# This is specific to the asNMPC code: 
# The initial constraints have to be defined as variables.
# They have to be set explictly with initial constraints.
# These constraints need to be identified by the 
# sens_init_constr suffix.
var c_init_var;
var t_init_var;

c_init_constr: c_init_var = c_init;
t_init_constr: t_init_var = t_init; #0.7293;
#---------------------------------

# collocation equations
fecolc{i in fe diff{1},j in cp}: c[i,j] = c[i-1,ncp]+time*h[i]*sum{k in cp} a[k,j]*cdot[i,k];
fecolt{i in fe diff{1},j in cp}: t[i,j] = t[i-1,ncp]+time*h[i]*sum{k in cp} a[k,j]*tdot[i,k];

fecolc0{i in 1..1,j in cp}:       c[i,j] = c_init_var+time*h[i]*sum{k in cp} a[k,j]*cdot[i,k];
fecolt0{i in 1..1,j in cp}:       t[i,j] = t_init_var+time*h[i]*sum{k in cp} a[k,j]*tdot[i,k];

# objective function...

minimize cost: sum{i in 2..nfe} (h[i]*sum{j in cp} ((alpha1*(c[i,j]-c_des)^2+ alpha2*(t[i,j]-t_des)^2+alpha3*(u[i,j]-u_des)^2 )*a[j,ncp])) + h[1]*sum{j in cp} ((alpha1*((c_init_var+time*h[1]*sum{k in cp} a[k,j]*cdot[1,k]) - c_des)^2 + alpha2*((t_init_var+time*h[1]*sum{k in cp} a[k,j]*tdot[1,k])-t_des)^2 + alpha3*(u[1,j]-u_des)^2)*a[j,ncp]);

#-- end of the hicks.mod file --
