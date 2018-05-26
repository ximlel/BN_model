RHL(1) = 2*U1-A_L*lo_L-A_R*lo_R;
RHL(2) = 2*U3-0.5*U2^2/A_L/lo_L-A_L*p_L/(gama-1)-0.5*U2^2/A_R/lo_R-A_R*p_R/(gama-1);
RHL(3) = 0.5*U2^2/A_L^2/lo_L^2+gama/(gama-1)*p_L/lo_L-0.5*U2^2/A_R^2/lo_R^2-gama/(gama-1)*p_R/lo_R;
RHL(4) = p_L/lo_L^gama-p_R/lo_R^gama;
lo_L = (2*U1-A_R*lo_R)/A_L;
p_R  = (2*U3-0.5*U2^2/(2*U1-A_R*lo_R)-0.5*U2^2/A_R/lo_R)*(gama-1)/(A_L/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama+A_R);
p_L  = (2*U3-0.5*U2^2/(2*U1-A_R*lo_R)-0.5*U2^2/A_R/lo_R)*(gama-1)/(A_L/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama+A_R)/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama;
