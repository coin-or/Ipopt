g0 0 0 0	# toy problem - the 'g' indicates the file here is in ascii (other numbers (0,0,0,0) ignored)
4 3 1 0 0	# num vars, num constraints, num objectives, 0,0
1 0	        # number nonlinear constraints, num nonlinear objectives
0 0	        # num nonlinear network constraints, num linear network constraints, 
2 0 0	        # num nonlinear vars in constraints, objectives, both
0 0 0 1	        # num linear network variables; functions; arith, flags
1 1 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
7 3	        # number nonzeros in Jacobian, number nonzeros in linear part of objectives
0 0	        # 0,0
0 0 0 0 0	# 0,0,0,0,0    // Note that the variables are in order y[1],y[2],x,z
C0	#           Here starts the nonlinear part of constraint-0
o0	#             = add
o5	#             = raise-to-power
o0	#             = add
n-0.5   #             = constant -0.5
v0	#             = the continuous variable y[1]
n2      #             = constant 2
o5	#             = raise-to-power
o0	#             = add
n-0.5   #             = constant -0.5
v1	#             = the continuous variable y[2]
n2      #             = constant 2     ....Here ends the nonlinear part of constraint-0                   
C1	#           Here starts the nonlinear part of constraint-1
n0      #             Here ends the nonlinear part of constraint-1                              
C2	#           Here starts the nonlinear part of constraint-2
n0      #             Here ends the nonlinear part of constraint-2                               
O0 0	#           Here starts the nonlinear part of objective(-0)
n0      #             Here ends the nonlinear part of objective(-0)
r	#           Here starts the right-hand side section
1 0.25  #             constraint-0 body <= 0.25
1 0     #             constraint-1 body <=0
1 2     #             constraint-2 body <= 2
b	#           Here starts the bounds on variables section
2 0     #             var-0 >= 0
2 0     #             var-1 >= 0
0 0 1   #             0 <= var-2 <= 1
0 0 5   #             0 <= var-3 <= 5
k3	#           INFO FOR JACOBIAN COLUMN COUNTS. THIS RECORD IS ALWAYS 'k' FOLLOWED BY NUMBER OF VARIABLES MINUS 1
2       #             the number of nonzeros in the column2 0..0 of the jacobian
4       #             the number of nonzeros in the columns 0..1 of the jacobian
6       #             the number of nonzeros in the columns 0..2 of the jacobian 
J0 2    #           Two coefficients for the row-0 of the (explicitly) linear part of the constraints
0 0     #              var-0, linear coefficient 0
1 0     #              var-1, linear coefficient 0
J1 2    #           Two coefficients for the row-1 of the linear part of the constraints
0 -1    #             var-0, linear coefficient -1
2 1     #             var-2, linear coefficient 1
J2 3    #           Three coefficients for the row-2 of the linear part of the constraints
1 1     #             var-1, linear coefficient 1
2 1     #             var-2, linear coefficient 1
3 1     #             var-3, linear coefficient 1
G0 3    #           Three coefficients for the linear part of objective(-0)
0 -1    #             var-0, linear coefficient -1
1 -1    #             var-1, linear coefficient -1
2 -1    #             var-2, linear coefficient -1
