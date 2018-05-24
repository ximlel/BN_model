%compute flux
function [lo_L,u_L,p_L,lo_R,u_R,p_R]=primitive_comp(U,A_L,A_R)
global gama;
U1=U(1);
U2=U(2);
U3=U(3);
u=U2/U1;
% lo_L  = U1/A_L;
% p_L   = (U3/A_L - 0.5*lo_L*u^2)*(gama-1);
% lo_R  = U1/A_R;
% p_R   = (U3/A_R - 0.5*lo_R*u^2)*(gama-1);
A=0.5*(A_L+A_R);
lo_L  = U1/A;
p_L   = (U3/A - 0.5*lo_L*u^2)*(gama-1);
lo_R  = lo_L;
p_R   = p_L;
% if abs(A_L-A_R)>ep
%     lo_L  =151.13;
%     p_L   =2.4836e8;
%     lo_R  =95.199;
%     p_R   =1.4067e8;
% end
global ep;
fun  = zeros(1,4);
it_max = 50;
%it_max = 500;
k = 0; err = 1e50;
while (k<it_max && err>ep && abs(A_L-A_R)>ep)
    fun(1) = 2*U1-A_L*lo_L-A_R*lo_R;
    fun(2) = 2*U3-0.5*U2^2/A_L/lo_L-A_L*p_L/(gama-1)-0.5*U2^2/A_R/lo_R-A_R*p_R/(gama-1);
    fun(3) = 0.5*U2^2/A_L^2/lo_L^2+gama/(gama-1)*p_L/lo_L-0.5*U2^2/A_R^2/lo_R^2-gama/(gama-1)*p_R/lo_R;
    fun(4) = p_L/lo_L^gama-p_R/lo_R^gama;
    dfun=[
     -A_L, U2^2/(2*A_L*lo_L^2), - U2^2/(A_L^2*lo_L^3) - (gama*p_L)/(lo_L^2*(gama - 1)), -(gama*p_L)/lo_L^(gama + 1);
     -A_R, U2^2/(2*A_R*lo_R^2),   U2^2/(A_R^2*lo_R^3) + (gama*p_R)/(lo_R^2*(gama - 1)),  (gama*p_R)/lo_R^(gama + 1);
        0,     -A_L/(gama - 1),                                 gama/(lo_L*(gama - 1)),                 1/lo_L^gama;
        0,     -A_R/(gama - 1),                                -gama/(lo_R*(gama - 1)),                -1/lo_R^gama];
    [x_star, err] = NewtonRapshon(fun,dfun,[lo_L lo_R p_L p_R],ep);
    lo_L=max(x_star(1),ep);
    lo_R=max(x_star(2),ep);
    p_L=max(x_star(3),ep);
    p_R=max(x_star(4),ep);
    k=k+1;
end
if k>=it_max
    err
end
u_L=U2/A_L/lo_L;
u_R=U2/A_R/lo_R;
end
