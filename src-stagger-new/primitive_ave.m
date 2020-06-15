%compute flux
function [lo_gL_n,u_gL_n,p_gL_n,lo_sL_n,u_sL_n,p_sL_n,lo_gR_n,u_gR_n,p_gR_n,lo_sR_n,u_sR_n,p_sR_n]=primitive_ave(U,phi_sL,phi_sR)
global gama_s gama_g;
U1=U(1);
U2=U(2);
U3=U(3);
U4=U(4);
U5=U(5);
U6=U(6);
phi_s = 0.5*phi_sL+0.5*phi_sR;
phi_g = 1-phi_s;
lo_g = U1/phi_g;
u_g  = U2/U1;
p_g  = (U3/phi_g - 0.5*lo_g*u_g*2)*(gama_g-1);
lo_s = U4/phi_s;
u_s  = U5/U4;
p_s  = (U6/phi_s - 0.5*lo_s*u_s^2)*(gama_s-1);
[lo_gL_n,u_gL_n,p_gL_n,p_sL_n]=Riemann_inv(phi_s,lo_g,u_g,p_g,u_s,p_s,phi_sL);
[lo_gR_n,u_gR_n,p_gR_n,p_sR_n]=Riemann_inv(phi_s,lo_g,u_g,p_g,u_s,p_s,phi_sR);
lo_sL_n = lo_s;
lo_sR_n = lo_s;
u_sL_n = u_s;
u_sR_n = u_s;
end

