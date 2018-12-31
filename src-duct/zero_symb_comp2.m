syms lo_R;
syms A_L A_R;
syms U1 U2 U3;
syms u_L u_R;
syms gama;
RHL = 0.5*U2^2/A_L^2/((2*U1-A_R*lo_R)/A_L)^2+gama/(gama-1)*((2*U3-0.5*U2^2/(2*U1-A_R*lo_R)-0.5*U2^2/A_R/lo_R)*(gama-1)/(A_L/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama+A_R)/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama)/((2*U1-A_R*lo_R)/A_L)-0.5*U2^2/A_R^2/lo_R^2-gama/(gama-1)*((2*U3-0.5*U2^2/(2*U1-A_R*lo_R)-0.5*U2^2/A_R/lo_R)*(gama-1)/(A_L/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama+A_R))/lo_R;
dfun=diff(RHL,'lo_R')
