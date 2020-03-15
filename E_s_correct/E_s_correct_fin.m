%comput Roe mean
function [E_sum_I,E_sum_II]=E_s_correct_fin(S_sM,lo_g1,u_g1,p_g1,lo_s1,p_s1,phi_s1,lo_g2,u_g2,p_g2,lo_s2,p_s2,phi_s2,ratio_t_x)
phi_g1 = 1.0-phi_s1;
phi_g2 = 1.0-phi_s2;
global gama_s gama_g p0;

P_M = phi_s1*p_s1+phi_g1*p_g1+phi_g1*lo_g1*(u_g1-S_sM)^2;
P_M = phi_s2*p_s2+phi_g2*p_g2+phi_g2*lo_g2*(u_g2-S_sM)^2;

Q_M = phi_g1*lo_g1*(u_g1-S_sM);
Q_M = phi_g2*lo_g2*(u_g2-S_sM);
H_M = gama_g/(gama_g-1.0)*p_g1/lo_g1+0.5*(u_g1-S_sM)^2;
H_M = gama_g/(gama_g-1.0)*p_g2/lo_g2+0.5*(u_g2-S_sM)^2;
eta_M = p_g1/lo_g1^gama_g;
eta_M = p_g2/lo_g2^gama_g;

KL = ratio_t_x*S_sM;
KR = 1.0-KL;
E_s1 = (p_s1+gama_s*p0)/(gama_s-1)+0.5*lo_s1*S_sM^2;
E_s2 = (p_s2+gama_s*p0)/(gama_s-1)+0.5*lo_s2*S_sM^2;
E_g1 = p_g1/(gama_g-1)+0.5*lo_g1*u_g1^2;
E_g2 = p_g2/(gama_g-1)+0.5*lo_g2*u_g2^2;
E_sum_II =  KL*(phi_g1*E_g1+phi_s1*E_s1) + KR*(phi_g2*E_g2+phi_s2*E_s2);

A_sM = KL*phi_s1*lo_s1+KR*phi_s2*lo_s2;
A_gM = KL*phi_g1*lo_g1+KR*phi_g2*lo_g2;
%M = KL*(phi_s1*lo_s1*u_s1+phi_g1*lo_g1*u_g1)+KR*(phi_s2*lo_s2*u_s2+phi_g2*lo_g2*u_g2);

u_gM = Q_M/A_gM + S_sM;
lo_gM = ((H_M-(u_gM-S_sM)^2/2)/eta_M*(gama_g-1)/gama_g)^(1/(gama_g-1));
p_gM = lo_gM^gama_g*eta_M;
phi_gM = A_gM/lo_gM;
phi_sM = 1.0-phi_gM;
lo_sM = A_sM/phi_sM;
p_sM = (P_M-Q_M*(u_gM-S_sM)-phi_gM*p_gM)/phi_sM;
E_sum_I = phi_gM*(p_gM/(gama_g-1)+0.5*lo_gM*u_gM^2)+phi_sM*((p_sM+gama_s*p0)/(gama_s-1)+0.5*lo_sM*S_sM^2);

end
