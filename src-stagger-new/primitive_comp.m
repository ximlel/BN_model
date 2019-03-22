%compute flux
function [lo_gL,u_gL,p_gL,lo_sL,u_sL,p_sL,lo_gR,u_gR,p_gR,lo_sR,u_sR,p_sR]=primitive_comp(U,phi_sL,phi_sR,area_L,area_R)
phi_gL=1-phi_sL;
phi_gR=1-phi_sR;
global gama_s gama_g;
global ep;
U1=U(1);
U2=U(2);
U3=U(3);
U4=U(4);
U5=U(5);
U6=U(6);
phi_s = area_L*phi_sL+area_R*phi_sR;
phi_g = 1-phi_s;
lo_s  = U4/phi_s;
u_s   = U5/U4;
lo_gR = U1/phi_g;
u_gR  = U2/U1;
p_gR  = (U3/phi_g - 0.5*lo_gR*u_gR^2)*(gama_g-1);
p_sR  = (U6/phi_s - 0.5*lo_s*u_s^2)*(gama_s-1);

% lo_gR = U1/phi_gR;
% u_gR  = U2/U1;
% p_gR  = (U3/phi_gR - 0.5*lo_gR*u_gR^2)*(gama_s-1);
% p_sR  = (U6/phi_sR - 0.5*lo_s*u_s^2)*(gama_s-1);
% if abs(phi_sL-phi_sR)>ep
%     lo_gL   =1;
%     u_gL    =2;
%     p_gL    =1;
%     lo_s    =2;
%     u_s     =0.3;
%     p_sL    =5;
%     lo_gR   =0.1941934235006083;
%     u_gR    =2.801188129642115;
%     p_gR    =0.1008157360849781;
%     p_sR    =12.85675006887399;
% end
fun  = zeros(1,4);
it_max = 500;
k = 0; err2 = 1e50;
while (k<it_max && err2>ep && abs(phi_sL-phi_sR)>ep)
    fun(1) = U3-area_L*phi_gL*(0.5*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR))^2+(p_gR/lo_gR^gama_g*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)^gama_g)/(gama_g-1))-area_R*phi_gR*(0.5*lo_gR*u_gR^2+p_gR/(gama_g-1));
    fun(2) = phi_gL*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*(((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR))-u_s)-phi_gR*lo_gR*(u_gR-u_s);
    fun(3) = phi_gL*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)*(((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR))-u_s)^2+phi_gL*(p_gR/lo_gR^gama_g*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)^gama_g)+phi_sL*(((U6-0.5*phi_s*lo_s*u_s^2)*(gama_s-1)-area_R*phi_sR*p_sR)/area_L/phi_sL)-phi_gR*lo_gR*(u_gR-u_s)^2-phi_gR*p_gR-phi_sR*p_sR;
    fun(4) = 0.5*(((U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR))-u_s)^2+gama_g/(gama_g-1)*(p_gR/lo_gR^gama_g*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)^gama_g)/((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)-0.5*(u_gR-u_s)^2-gama_g/(gama_g-1)*p_gR/lo_gR;
    dfun=[
	 area_L*phi_gL*((gama_g*p_gR*((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL))^gama_g)/(lo_gR^(gama_g + 1)*(gama_g - 1)) - (area_R*phi_gR*(U2 - area_R*lo_gR*phi_gR*u_gR)^2)/(2*area_L*phi_gL*(U1 - area_R*lo_gR*phi_gR)^2) + (area_R*phi_gR*u_gR*(U2 - area_R*lo_gR*phi_gR*u_gR))/(area_L*phi_gL*(U1 - area_R*lo_gR*phi_gR)) + (area_R*gama_g*p_gR*phi_gR*((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL))^(gama_g - 1))/(area_L*lo_gR^gama_g*phi_gL*(gama_g - 1))) - (area_R*phi_gR*u_gR^2)/2, (area_R*phi_gR*(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR)))/area_L - ((U1 - area_R*lo_gR*phi_gR)*((area_R*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR) - (area_R*phi_gR*(U2 - area_R*lo_gR*phi_gR*u_gR))/(U1 - area_R*lo_gR*phi_gR)^2))/area_L - phi_gR*(u_gR - u_s), (2*(U1 - area_R*lo_gR*phi_gR)*((area_R*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR) - (area_R*phi_gR*(U2 - area_R*lo_gR*phi_gR*u_gR))/(U1 - area_R*lo_gR*phi_gR)^2)*(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR)))/area_L - phi_gR*(u_gR - u_s)^2 - (area_R*phi_gR*(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR))^2)/area_L - (gama_g*p_gR*phi_gL*((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL))^gama_g)/lo_gR^(gama_g + 1) - (area_R*gama_g*p_gR*phi_gR*((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL))^(gama_g - 1))/(area_L*lo_gR^gama_g), ((area_R*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR) - (area_R*phi_gR*(U2 - area_R*lo_gR*phi_gR*u_gR))/(U1 - area_R*lo_gR*phi_gR)^2)*(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR)) + (gama_g*p_gR)/(lo_gR^2*(gama_g - 1)) - (area_L*gama_g^2*p_gR*phi_gL*((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL))^gama_g)/(lo_gR^(gama_g + 1)*(U1 - area_R*lo_gR*phi_gR)*(gama_g - 1)) - (area_R*gama_g^2*p_gR*phi_gR*((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL))^(gama_g - 1))/(lo_gR^gama_g*(U1 - area_R*lo_gR*phi_gR)*(gama_g - 1)) + (area_L*area_R*gama_g*p_gR*phi_gL*phi_gR*((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL))^gama_g)/(lo_gR^gama_g*(U1 - area_R*lo_gR*phi_gR)^2*(gama_g - 1));
	                                                                                                                                                                                                                                                                                                                                                      - (area_R*phi_gR)/(gama_g - 1) - (area_L*phi_gL*((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL))^gama_g)/(lo_gR^gama_g*(gama_g - 1)),                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          (phi_gL*((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL))^gama_g)/lo_gR^gama_g - phi_gR,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (area_L*gama_g*phi_gL*((U1 - area_R*lo_gR*phi_gR)/(area_L*phi_gL))^gama_g)/(lo_gR^gama_g*(U1 - area_R*lo_gR*phi_gR)*(gama_g - 1)) - gama_g/(lo_gR*(gama_g - 1));
	                                                                                                                                                                                                                                                                                                                                                                           (area_R*lo_gR*phi_gR*(U2 - area_R*lo_gR*phi_gR*u_gR))/(U1 - area_R*lo_gR*phi_gR) - area_R*lo_gR*phi_gR*u_gR,                                                                                                                                                                                                                                          - lo_gR*phi_gR - (area_R*lo_gR*phi_gR)/area_L,                                                                                                                                                                                                                                                                                                                                                                                                                                           (2*area_R*lo_gR*phi_gR*(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR)))/area_L - lo_gR*phi_gR*(2*u_gR - 2*u_s),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      u_s - u_gR + (area_R*lo_gR*phi_gR*(u_s - (U2 - area_R*lo_gR*phi_gR*u_gR)/(U1 - area_R*lo_gR*phi_gR)))/(U1 - area_R*lo_gR*phi_gR);
	                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     0,                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           - phi_sR - (area_R*phi_sR)/area_L,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     0];
    [x_star, err2] = NewtonRapshon(fun,dfun,[lo_gR p_gR u_gR p_sR],ep);
    lo_gR=max(real(x_star(1)),ep);
    lo_gR=min(lo_gR,U1/area_R/phi_gR-ep);
    p_gR =max(real(x_star(2)),ep);
    u_gR =real(x_star(3));
    % u_gR =max(real(x_star(3)),ep);
    p_sR =max(real(x_star(4)),ep);
    p_sR =min(p_sR,(U6-0.5*phi_s*lo_s*u_s^2)*(gama_s-1)/area_R/phi_sR-ep);
    k=k+1;
end
if k>=it_max
    err2
end
lo_gL= (U1-area_R*phi_gR*lo_gR)/area_L/phi_gL;
u_gL = (U2-area_R*phi_gR*lo_gR*u_gR)/(U1-area_R*phi_gR*lo_gR);
p_gL = p_gR/lo_gR^gama_g*((U1-area_R*phi_gR*lo_gR)/area_L/phi_gL)^gama_g;
p_sL = ((U6-0.5*phi_s*lo_s*u_s^2)*(gama_s-1)-area_R*phi_sR*p_sR)/area_L/phi_sL;
lo_sL= lo_s;
lo_sR= lo_s;
% lo_sL= U_lo_sL/phi_sL;
% lo_sR= U_lo_sR/phi_sR;
u_sL = u_s;
u_sR = u_s;
end
