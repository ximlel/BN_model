syms U1 U2 U3 U4 U5 U6;
syms lo_gL u_gL p_gL lo_gR u_gR p_gR p_sL p_sR;
syms u_s lo_s;
syms gama_s gama_g;
syms area_L area_R phi_gL phi_gR phi_sL phi_sR;
fun(1) = U1-area_L*phi_gL*lo_gL-area_R*phi_gR*lo_gR;
fun(2) = U2-area_L*phi_gL*lo_gL*u_gL-area_R*phi_gR*lo_gR*u_gR;
fun(3) = U3-area_L*phi_gL*(0.5*lo_gL*u_gL^2+p_gL/(gama_g-1))-area_R*phi_gR*(0.5*lo_gR*u_gR^2+p_gR/(gama_g-1));
fun(4) = U6-area_L*phi_sL*(0.5*lo_s*u_s^2  +p_sL/(gama_s-1))-area_R*phi_sR*(0.5*lo_s*u_s^2  +p_sR/(gama_s-1));
fun(5) = phi_gL*lo_gL*(u_gL-u_s)-phi_gR*lo_gR*(u_gR-u_s);
fun(6) = phi_gL*lo_gL*(u_gL-u_s)^2+phi_gL*p_gL+phi_sL*p_sL-phi_gR*lo_gR*(u_gR-u_s)^2-phi_gR*p_gR-phi_sR*p_sR;
fun(7) = 0.5*(u_gL-u_s)^2+gama_g/(gama_g-1)*p_gL/lo_gL-0.5*(u_gR-u_s)^2-gama_g/(gama_g-1)*p_gR/lo_gR;
fun(8) = p_gL/lo_gL^gama_g-p_gR/lo_gR^gama_g;
dfun=[diff(fun,'lo_gL');diff(fun,'lo_gR');diff(fun,'p_gL');diff(fun,'p_gR');diff(fun,'u_gL');diff(fun,'u_gR');diff(fun,'p_sL');diff(fun,'p_sR')]
