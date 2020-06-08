syms U1 U2 U3 U4 U5 U6;
syms lo_gR p_gR;
syms u_s lo_s;
syms gama_s gama_g;
syms area_L area_R phi_gL phi_gR phi_sL phi_sR phi_s;
%fun(7) = 0.5*(u_gL-u_s)^2+gama_g/(gama_g-1)*p_gL/lo_gL-0.5*(u_gR-u_s)^2-gama_g/(gama_g-1)*p_gR/lo_gR;
%fun(8) = p_gL/lo_gL^gama_g-p_gR/lo_gR^gama_g;
fun(1) = 0.5*((U2-(U1+U4)*u_s)/(U1-area_R*phi_gR*lo_gR)*area_L)^2+gama_g/(gama_g-1)*(((U3-0.5*area_R*phi_gR*lo_gR*((U2-(U1+U4)*u_s)/phi_gR/lo_gR+u_s)^2-0.5*area_L*phi_gL*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*((U2-(U1+U4)*u_s)/(U1-area_R*phi_gR*lo_gR)*area_L+u_s)^2)*(gama_g-1)-area_R*phi_gR*p_gR)/area_L/phi_gL)/((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)-0.5*((U2-(U1+U4)*u_s)/phi_gR/lo_gR)^2-gama_g/(gama_g-1)*p_gR/lo_gR;
fun(2) = (((U3-0.5*area_R*phi_gR*lo_gR*((U2-(U1+U4)*u_s)/phi_gR/lo_gR+u_s)^2-0.5*area_L*phi_gL*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*((U2-(U1+U4)*u_s)/(U1-area_R*phi_gR*lo_gR)*area_L+u_s)^2)*(gama_g-1)-area_R*phi_gR*p_gR)/area_L/phi_gL)/((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)^gama_g-p_gR/lo_gR^gama_g;
dfun=[diff(fun,'lo_gR');diff(fun,'p_gR')]
