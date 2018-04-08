%compute flux
function [lo_g u_g p_g phi_g lo_s u_s p_s phi_s]=primitive_comp(U)
global gama_s gama_g p0;
phi_s = U(7);
phi_g = 1.0-phi_s;
lo_g  = U(1)/phi_g;
u_g   = U(2)/U(1);
p_g   = (U(3)/phi_g - 0.5*lo_g*u_g^2)*(gama_g-1);
lo_s  = U(4)/phi_s;
u_s   = U(5)/U(4);
p_s   = (U(6)/phi_s - 0.5*lo_s*u_s^2)*(gama_s-1)-gama_s*p0;
end
