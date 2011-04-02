# ===================================================================
# dynamic optimization formulation of the hicks-ray reactor
# create results file for matlab
# victor m zavala  march 2006
# Adapted for sIPOPT by Hans Pirnay 2009, 2011
# ===================================================================

printf "clc\n" >graph.m;
printf "clear all" >graph.m;
printf "\n" >graph.m;

#display variables and parameters

printf "c=[" >graph.m;
printf "%5.5f ",c_init>>graph.m;
printf {j in fe,k in cp}:
"%5.5f ",c[j,k]>>graph.m;
printf "]';\n" >>graph.m;

printf "t=[" >graph.m;
printf "%5.5f ",t_init>>graph.m;
printf {j in fe,k in cp}:
"%5.5f ",t[j,k]>>graph.m;
printf "]';\n" >>graph.m;

printf "u=[" >graph.m;
printf "%5.5f ",u_init>>graph.m;
printf {j in fe,k in cp}:
"%5.5f ",u[j,k]>>graph.m;
printf "]';\n" >>graph.m;

printf "c1=[" >graph.m;
printf "%5.5f ",c_init>>graph.m;
printf {j in fe,k in cp}:
"%5.5f ",c[j,k].sens_sol_state_1>>graph.m;
printf "]';\n" >>graph.m;

printf "t1=[" >graph.m;
printf "%5.5f ",t_init>>graph.m;
printf {j in fe,k in cp}:
"%5.5f ",t[j,k].sens_sol_state_1>>graph.m;
printf "]';\n" >>graph.m;

printf "u1=[" >graph.m;
printf "%5.5f ",u_init>>graph.m;
printf {j in fe,k in cp}:
"%5.5f ",u[j,k].sens_sol_state_1>>graph.m;
printf "]';\n" >>graph.m;

printf "nfe= %5i ",nfe>>graph.m; printf ";\n">>graph.m;
printf "ncp= %5i ",ncp>>graph.m; printf ";\n">>graph.m;
printf "ncp= ncp-1">>graph.m; printf ";\n">>graph.m;
printf "tt= %5.5f ",time>>graph.m; printf ";\n">>graph.m;
printf "r1= %5.5f ",r1>>graph.m; printf ";\n">>graph.m;
printf "r2= %5.5f ",r2>>graph.m; printf ";\n">>graph.m;
printf "r3= %5.5f ",r3>>graph.m; printf ";\n">>graph.m;
printf "roots(1)= r1 ">>graph.m; printf ";\n">>graph.m;
printf "roots(2)= r2 ">>graph.m; printf ";\n">>graph.m;
printf "roots(3)= r3 ">>graph.m; printf ";\n">>graph.m;
printf "h= tt/nfe ">>graph.m; printf ";\n">>graph.m;
printf "index= 0 ">>graph.m; printf ";\n">>graph.m;
printf "nelement= 0 ">>graph.m; printf ";\n">>graph.m;
printf "for i=1:nfe,">>graph.m; printf "\n">>graph.m;
printf "nelement=nelement+1">>graph.m; printf ";\n">>graph.m;
printf "for j=1:ncp+1,">>graph.m; printf "\n">>graph.m;
printf "index=index+1">>graph.m; printf ";\n">>graph.m;
printf "if i == 1">>graph.m; printf "\n">>graph.m;
printf "time(index) = h*roots(j)">>graph.m; printf ";\n">>graph.m;
printf "else">>graph.m; printf "\n">>graph.m;
printf "time(index) = (nelement-1)*h+h*roots(j)">>graph.m; printf ";\n">>graph.m;
printf "end">>graph.m; printf "\n">>graph.m;
printf "end">>graph.m; printf "\n">>graph.m;
printf "end">>graph.m; printf "\n">>graph.m;
printf "time= [0,time]'">>graph.m; printf ";\n">>graph.m;

# plot the results

printf "subplot(3,1,1)\n" >>graph.m;
printf "plot(time,c,'k')\nhold on\n" >>graph.m;
printf "plot(time,c1,'k')\n" >>graph.m;
printf "xlabel ('time')">>graph.m; printf ";\n">>graph.m;
printf "ylabel ('concentration')">>graph.m; printf ";\n">>graph.m;
printf "subplot(3,1,2)\n" >>graph.m;
printf "plot(time,[t,t1],'k')\n" >>graph.m;
printf "xlabel ('time')">>graph.m; printf ";\n">>graph.m;
printf "ylabel ('temperature')">>graph.m; printf ";\n">>graph.m;
printf "subplot(3,1,3)\n" >>graph.m;
printf "plot(time,[u,u1],'k')\n" >>graph.m;
printf "xlabel ('time')">>graph.m; printf ";\n">>graph.m;
printf "ylabel ('cooling w flowrate')">>graph.m; printf ";\n">>graph.m;

printf "print -deps hicks.eps\n">>graph.m;

close graph.m;
printf "output written to matlab m-file: graph.m ...\n";

#-- end of the graphshicks.inc file -

