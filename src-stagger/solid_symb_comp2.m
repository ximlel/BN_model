syms U1 U2 U3 U4 U5 U6;
syms lo_gR u_gR p_gR p_sR;
syms u_s lo_s;
syms gama_s gama_g;
syms area_L area_R phi_gL phi_gR phi_sL phi_sR phi_s;
fun(1) = U3-area_L*phi_gL*(0.5*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR))^2+(p_gR/lo_gR^gama_g*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)^gama_g)/(gama_g-1))-area_R*phi_gR*(0.5*lo_gR*u_gR^2+p_gR/(gama_g-1));
fun(2) = phi_gL*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*(((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR))-u_s)-phi_gR*lo_gR*(u_gR-u_s);
fun(3) = phi_gL*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*(((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR))-u_s)^2+phi_gL*(p_gR/lo_gR^gama_g*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)^gama_g)+phi_sL*(((U6-0.5*phi_s*lo_s*u_s^2)*(gama_s-1)-area_R*phi_sR*p_sR)/area_L/phi_sL)-phi_gR*lo_gR*(u_gR-u_s)^2-phi_gR*p_gR-phi_sR*p_sR;
fun(4) = 0.5*(((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR))-u_s)^2+gama_g/(gama_g-1)*(p_gR/lo_gR^gama_g*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)^gama_g)/((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)-0.5*(u_gR-u_s)^2-gama_g/(gama_g-1)*p_gR/lo_gR;
dfun=[diff(fun,'lo_gR');diff(fun,'p_gR');diff(fun,'u_gR');diff(fun,'p_sR')]
