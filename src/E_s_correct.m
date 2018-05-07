%comput Roe mean
function [phi_sM, dE_g]=E_s_correct(S_sM,lo_g1,u_g1,p_g1,lo_s1,p_s1,phi_s1,lo_g2,u_g2,p_g2,lo_s2,p_s2,phi_s2,ratio_t_x)
global gama_s gama_g p0;
phi_g1 = 1.0-phi_s1;
phi_g2 = 1.0-phi_s2;
KL = ratio_t_x*S_sM;
KR = 1.0-KL;
phi_sM = KL*phi_s1+KR*phi_s2;
lo_sM  = (KL*phi_s1*lo_s1+KR*phi_s2*lo_s2)/phi_sM;
phi_gM = 1.0-phi_sM;
lo_gM  = (KL*phi_g1*lo_g1+KR*phi_g2*lo_g2)/phi_gM;
u_gM   = (KL*phi_g1*lo_g1*u_g1+KR*phi_g2*lo_g2*u_g2)/lo_gM/phi_gM;
HL = 0.5*(u_g1-u_s1)^2+gama_g/(gama_g-1)*p_g1/rho_g1;
HR = 0.5*(u_g2-u_s2)^2+gama_g/(gama_g-1)*p_g2/rho_g2;
HM = HL*KL+HR*KR;
E_gM = (HM-0.5*(u_gM-S_sM)^2)/gama_g*phi_gM*lo_gM + 0.5*phi_gM*lo_gM*u_gM^2;
E_s1 = (p_s1+gama_s*p0)/(gama_s-1)+0.5*lo_s1*S_sM^2;
E_s2 = (p_s2+gama_s*p0)/(gama_s-1)+0.5*lo_s2*S_sM^2;
E_g1 = p_g1/(gama_g-1)+0.5*lo_g1*u_g1^2;
E_g2 = p_g2/(gama_g-1)+0.5*lo_g2*u_g2^2;
E_gsum = phi_g2*E_g2 + KL*(phi_s1*p_s1-phi_s2*p_s2)+ratio_t_x*(phi_g1*u_g1*(p_g1+E_g1)-phi_g2*u_g2*(p_g2+E_g2);
E_sum = phi_g2*E_g2+phi_s2*E_s2 + KL*(phi_s1*(p_s1+E_s1)-phi_s2*(p_s2+E_s2))+ratio_t_x*(phi_g1*u_g1*(p_g1+E_g1)-phi_g2*u_g2*(p_g2+E_g2));
e_sum = E_sum - 0.5*phi_sM*lo_sM*S_sM^2 - 0.5*phi_gM*lo_gM*u_gM^2;
p_sum = phi_s1*p_s1+phi_g1*p_g1+phi_g1*lo_g1*(u_g1-S_sM)^2-phi_gM*lo_gM*(u_gM-S_sM)^2;
p_sM = (e_sum/phi_sM - p_sum/(gama_g-1)/phi_sM - gama_s*p0/(gama_s-1))/(1/(gama_s-1)-1/(gama_g-1));
E_sM = KL*phi_s1*E_s1+KR*phi_s2*E_s2;
dE_s = E_sM-phi_sM*(p_sM+gama_s*p0)/(gama_s-1)-0.5*phi_sM*lo_sM*S_sM^2;
end
