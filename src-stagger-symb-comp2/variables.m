fun(1) = U1-area_L*phi_gL*lo_gL-area_R*phi_gR*lo_gR;
fun(2) = U2-area_L*phi_gL*lo_gL*u_gL-area_R*phi_gR*lo_gR*u_gR;
fun(3) = U3-area_L*phi_gL*(0.5*lo_gL*u_gL^2+p_gL/(gama_g-1))-area_R*phi_gR*(0.5*lo_gR*u_gR^2+p_gR/(gama_g-1));
fun(4) = U6-area_L*phi_sL*(0.5*lo_s*u_s^2  +p_sL/(gama_s-1))-area_R*phi_sR*(0.5*lo_s*u_s^2  +p_sR/(gama_s-1));
fun(4) = U6-area_L*phi_sL*(0.5*lo_s*u_s^2  +p_sL/(gama_s-1))-area_R*phi_sR*(0.5*lo_s*u_s^2  +p_sR/(gama_s-1));
fun(5) = phi_gL*lo_gL*(u_gL-u_s)-phi_gR*lo_gR*(u_gR-u_s);
fun(6) = phi_gL*lo_gL*(u_gL-u_s)^2+phi_gL*p_gL+phi_sL*p_sL-phi_gR*lo_gR*(u_gR-u_s)^2-phi_gR*p_gR-phi_sR*p_sR;

Q = (U2-(U1+U4)*u_s);
u_gR = (Q/phi_gR/lo_gR+u_s);
lo_gL  = ((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL);
%u_gL  = Q/phi_gL/lo_gL+u_s;
u_gL = (Q/(U1-area_R*phi_gR*lo_gR)*area_L+u_s);
%p_gL  = (((U3-0.5*area_R*phi_gR*lo_gR*u_gR^2-0.5*area_L*phi_gL*lo_gL*u_gL^2)*(gama_g-1)-area_R*phi_gR*p_gR)/area_L/phi_gL);
p_gL  = (((U3-0.5*area_R*phi_gR*lo_gR*(Q/phi_gR/lo_gR+u_s)^2-0.5*area_L*phi_gL*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*(Q/(U1-area_R*phi_gR*lo_gR)*area_L+u_s)^2)*(gama_g-1)-area_R*phi_gR*p_gR)/area_L/phi_gL);
