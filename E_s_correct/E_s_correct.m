%comput Roe mean
function [phi_s_out,dE_g]=E_s_correct(S_sM,lo_g1,u_g1,p_g1,lo_s1,p_s1,phi_s1,lo_g2,u_g2,p_g2,lo_s2,p_s2,phi_s2,ratio_t_x)
phi_g1 = 1.0-phi_s1;
phi_g2 = 1.0-phi_s2;

P_M = phi_s1*p_s1+phi_g1*p_g1+phi_g1*lo_g1*(u_g1-S_sM)^2
P_M = phi_s2*p_s2+phi_g2*p_g2+phi_g2*lo_g2*(u_g2-S_sM)^2

global gama_s gama_g p0;
KL = ratio_t_x*S_sM;
KR = 1.0-KL;
phi_sM = KL*phi_s1+KR*phi_s2;
lo_sM  = (KL*phi_s1*lo_s1+KR*phi_s2*lo_s2)/phi_sM;
phi_gM = 1.0-phi_sM;
lo_gM  = (KL*phi_g1*lo_g1+KR*phi_g2*lo_g2)/phi_gM;
u_gM   = (KL*phi_g1*lo_g1*u_g1+KR*phi_g2*lo_g2*u_g2)/lo_gM/phi_gM;
E_s1 = (p_s1+gama_s*p0)/(gama_s-1)+0.5*lo_s1*S_sM^2;
E_s2 = (p_s2+gama_s*p0)/(gama_s-1)+0.5*lo_s2*S_sM^2;
E_g1 = p_g1/(gama_g-1)+0.5*lo_g1*u_g1^2;
E_g2 = p_g2/(gama_g-1)+0.5*lo_g2*u_g2^2;
E_gsum = phi_g2*E_g2 + KL*(phi_s1*p_s1-phi_s2*p_s2)+ratio_t_x*(phi_g1*u_g1*(p_g1+E_g1)-phi_g2*u_g2*(p_g2+E_g2));
E_sum = phi_g2*E_g2+phi_s2*E_s2 + KL*(phi_s1*(p_s1+E_s1)-phi_s2*(p_s2+E_s2))+ratio_t_x*(phi_g1*u_g1*(p_g1+E_g1)-phi_g2*u_g2*(p_g2+E_g2));
e_sum = E_sum - 0.5*phi_sM*lo_sM*S_sM^2 - 0.5*phi_gM*lo_gM*u_gM^2;

P_sum = (gama_s-1)*e_sum + phi_gM*lo_gM*(u_gM-S_sM)^2

%Riemann invariants
HL = 0.5*(u_g1-S_sM)^2+gama_g/(gama_g-1)*p_g1/lo_g1;

AP = (HL-0.5*(u_gM-S_sM)^2)*(gama_g-1)*phi_gM*lo_gM/gama_g;
(gama_s-1)*(e_sum-AP/(gama_g-1))
P_M - phi_gM*lo_gM*(u_gM-S_sM)^2 - AP

E_gM = (HL-0.5*(u_gM-S_sM)^2)/gama_g*phi_gM*lo_gM + 0.5*phi_gM*lo_gM*u_gM^2; %
dE_g = E_gM-E_gsum;
phi_p_gM = (HL-0.5*(u_gM-S_sM)^2)/gama_g*(gama_g-1)*phi_gM*lo_gM;
P_M = phi_s1*p_s1+phi_g1*p_g1+phi_g1*lo_g1*(u_g1-S_sM)^2-phi_gM*lo_gM*(u_gM-S_sM)^2;
phi_p_sM = P_M - phi_p_gM;
phi_p_s0 = (e_sum - phi_p_gM/(gama_g-1))*(gama_s-1) - phi_p_sM;
phi_s_out = phi_p_s0/p0/gama_s;
end
