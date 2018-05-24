%compute flux
function [lo_gL,u_gL,p_gL,lo_sL,u_sL,p_sL,lo_gR,u_gR,p_gR,lo_sR,u_sR,p_sR]=primitive_comp(U,phi_sL,phi_sR,area_L,area_R);
phi_gL=1-phi_sL;
phi_gR=1-phi_sR;
global gama_s gama_g;
U1=U(1);
U2=U(2);
U3=U(3);
U4=U(4);
U5=U(5);
U6=U(6);
phi_s = area_L*phi_sL+area_R*phi_sR;
phi_g = 1-phi_s;
lo_gL = U1/phi_g;
u_gL  = U2/U1;
p_gL  = (U3/phi_g - 0.5*lo_gL*u_gL^2)*(gama_g-1);
lo_gR = lo_gL;
u_gR  = u_gL;
p_gR  = p_gL;
lo_s  = U1/phi_s;
u_s   = U5/U4;
p_sL  = (U3/phi_s - 0.5*lo_s*u_s^2)*(gama_s-1);
p_sR  = p_gL;
global ep;
fun  = zeros(1,8);
it_max = 50;
%it_max = 500;
k = 0; err = 1e50;
while (k<it_max && err>ep && abs(phi_sL-phi_sR)>ep)
    fun(1) = U1-area_L*phi_gL*lo_gL-area_R*phi_gR*lo_gR;
    fun(2) = U2-area_L*phi_gL*lo_gL*u_gL-area_R*phi_gR*lo_gR*u_gR;
    fun(3) = U3-area_L*phi_gL*(0.5*lo_gL*u_gL^2+p_gL/(gama_g-1))-area_R*phi_gR*(0.5*lo_gR*u_gR^2+p_gR/(gama_g-1));
    fun(4) = U6-area_L*phi_sL*(0.5*lo_s*u_s^2  +p_sL/(gama_s-1))-area_R*phi_sR*(0.5*lo_s*u_s^2  +p_sR/(gama_s-1));
    fun(5) = phi_gL*lo_gL*(u_gL-u_s)-phi_gR*lo_gR*(u_gR-u_s);
    fun(6) = phi_gL*lo_gL*(u_gL-u_s)^2-phi_gL*p_gL-phi_sL*p_sL-phi_gR*lo_gR*(u_gR-u_s)^2+phi_gR*p_gR+phi_sR*p_sR;
    fun(7) = 0.5*(u_gL-u_s)^2+gama_g/(gama_g-1)*p_gL/lo_gL-0.5*(u_gR-u_s)-gama_g/(gama_g-1)*p_gR/lo_gR;
    fun(8) = p_gL/lo_gL^gama_g-p_gR/lo_gR^gama_g;
    dfun=[
    -area_L*phi_gL,  -area_L*phi_gL*u_gL,     -(area_L*phi_gL*u_gL^2)/2,                             0,  phi_gL*(u_gL - u_s),          phi_gL*(u_gL - u_s)^2, -(gama_g*p_gL)/(lo_gL^2*(gama_g - 1)), -(gama_g*p_gL)/lo_gL^(gama_g + 1);
    -area_R*phi_gR,  -area_R*phi_gR*u_gR,     -(area_R*phi_gR*u_gR^2)/2,                             0, -phi_gR*(u_gR - u_s),         -phi_gR*(u_gR - u_s)^2,  (gama_g*p_gR)/(lo_gR^2*(gama_g - 1)),  (gama_g*p_gR)/lo_gR^(gama_g + 1);
                 0,                    0, -(area_L*phi_gL)/(gama_g - 1),                             0,                    0,                        -phi_gL,           gama_g/(lo_gL*(gama_g - 1)),                    1/lo_gL^gama_g;
                 0,                    0, -(area_R*phi_gR)/(gama_g - 1),                             0,                    0,                         phi_gR,          -gama_g/(lo_gR*(gama_g - 1)),                   -1/lo_gR^gama_g;
                 0, -area_L*lo_gL*phi_gL,     -area_L*lo_gL*phi_gL*u_gL,                             0,         lo_gL*phi_gL,  lo_gL*phi_gL*(2*u_gL - 2*u_s),                            u_gL - u_s,                                 0;
                 0, -area_R*lo_gR*phi_gR,     -area_R*lo_gR*phi_gR*u_gR,                             0,        -lo_gR*phi_gR, -lo_gR*phi_gR*(2*u_gR - 2*u_s),                                  -1/2,                                 0;
                 0,                    0,                             0, -(area_L*phi_sL)/(gama_s - 1),                    0,                        -phi_sL,                                     0,                                 0;
                 0,                    0,                             0, -(area_R*phi_sR)/(gama_s - 1),                    0,                         phi_sR,                                     0,                                 0];
    [x_star, err] = NewtonRapshon(fun,dfun,[lo_gL lo_gR p_gL p_gR u_gL u_gR p_sL p_sR],ep);
    lo_gL=max(x_star(1),ep);
    lo_gR=max(x_star(2),ep);
    p_gL =max(x_star(3),ep);
    p_gR =max(x_star(4),ep);
    u_gL =max(x_star(5),ep);
    u_gR =max(x_star(6),ep);
    p_sL =max(x_star(7),ep);
    p_sR =max(x_star(8),ep);
    k=k+1;
end
if k>=it_max
    err
end
lo_sL=lo_s;
lo_sR=lo_s;
u_sL =u_s;
u_sR =u_s;
end
