fun(1) = U1-area_L*phi_gL*lo_gL-area_R*phi_gR*lo_gR;
fun(2) = U2-area_L*phi_gL*lo_gL*u_gL-area_R*phi_gR*lo_gR*u_gR;
fun(4) = U6-area_L*phi_sL*(0.5*lo_s*u_s^2  +p_sL/(gama_s-1))-area_R*phi_sR*(0.5*lo_s*u_s^2  +p_sR/(gama_s-1));
fun(8) = p_gL/lo_gL^gama_g-p_gR/lo_gR^gama_g;
lo_gL = ((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL);
u_gL  = ((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR));
p_gL  = (p_gR/lo_gR^gama_g*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)^gama_g);
p_sL  = (((U6-0.5*phi_s*lo_s*u_s^2)*(gama_s-1)-area_R*phi_sR*p_sR)/area_L/phi_sL);
