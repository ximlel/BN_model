%compute flux
function [lo_L,u_L,p_L,lo_R,u_R,p_R]=primitive_comp(U,A_L,A_R)
global gama;
U1=U(1);
U2=U(2);
U3=U(3);
% lo_R  = U1/A_R;
A=0.5*(A_L+A_R);
lo_R = U1/A;
% if abs(A_L-A_R)>ep
%     lo_R  =95.199;
% end
global ep;
%it_max = 50;
it_max = 500;
k = 0; err = 1e50;
while (k<it_max && err>ep && abs(A_L-A_R)>ep)
    fun = 0.5*U2^2/A_L^2/((2*U1-A_R*lo_R)/A_L)^2+gama/(gama-1)*((2*U3-0.5*U2^2/(2*U1-A_R*lo_R)-0.5*U2^2/A_R/lo_R)*(gama-1)/(A_L/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama+A_R)/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama)/((2*U1-A_R*lo_R)/A_L)-0.5*U2^2/A_R^2/lo_R^2-gama/(gama-1)*((2*U3-0.5*U2^2/(2*U1-A_R*lo_R)-0.5*U2^2/A_R/lo_R)*(gama-1)/(A_L/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama+A_R))/lo_R;
    dfun= U2^2/(A_R^2*lo_R^3) + (A_R*U2^2)/(2*U1 - A_R*lo_R)^3 - (gama*(U2^2/(2*A_R*lo_R^2) - (A_R*U2^2)/(2*(2*U1 - A_R*lo_R)^2)))/(lo_R*(A_R + (A_L*((2*U1 - A_R*lo_R)/A_L)^gama)/lo_R^gama)) - (gama*(U2^2/(2*(2*U1 - A_R*lo_R)) - 2*U3 + U2^2/(2*A_R*lo_R)))/(lo_R^2*(A_R + (A_L*((2*U1 - A_R*lo_R)/A_L)^gama)/lo_R^gama)) + (gama*((A_L*gama*((2*U1 - A_R*lo_R)/A_L)^gama)/lo_R^(gama + 1) + (A_R*gama*((2*U1 - A_R*lo_R)/A_L)^(gama - 1))/lo_R^gama)*(U2^2/(2*(2*U1 - A_R*lo_R)) - 2*U3 + U2^2/(2*A_R*lo_R)))/(lo_R*(A_R + (A_L*((2*U1 - A_R*lo_R)/A_L)^gama)/lo_R^gama)^2) + (A_L*gama^2*((2*U1 - A_R*lo_R)/A_L)^gama*(U2^2/(2*(2*U1 - A_R*lo_R)) - 2*U3 + U2^2/(2*A_R*lo_R)))/(lo_R^(gama + 1)*(2*U1 - A_R*lo_R)*(A_R + (A_L*((2*U1 - A_R*lo_R)/A_L)^gama)/lo_R^gama)) + (A_R*gama^2*((2*U1 - A_R*lo_R)/A_L)^(gama - 1)*(U2^2/(2*(2*U1 - A_R*lo_R)) - 2*U3 + U2^2/(2*A_R*lo_R)))/(lo_R^gama*(2*U1 - A_R*lo_R)*(A_R + (A_L*((2*U1 - A_R*lo_R)/A_L)^gama)/lo_R^gama)) + (A_L*gama*((2*U1 - A_R*lo_R)/A_L)^gama*(U2^2/(2*A_R*lo_R^2) - (A_R*U2^2)/(2*(2*U1 - A_R*lo_R)^2)))/(lo_R^gama*(2*U1 - A_R*lo_R)*(A_R + (A_L*((2*U1 - A_R*lo_R)/A_L)^gama)/lo_R^gama)) - (A_L*gama*((2*U1 - A_R*lo_R)/A_L)^gama*((A_L*gama*((2*U1 - A_R*lo_R)/A_L)^gama)/lo_R^(gama + 1) + (A_R*gama*((2*U1 - A_R*lo_R)/A_L)^(gama - 1))/lo_R^gama)*(U2^2/(2*(2*U1 - A_R*lo_R)) - 2*U3 + U2^2/(2*A_R*lo_R)))/(lo_R^gama*(2*U1 - A_R*lo_R)*(A_R + (A_L*((2*U1 - A_R*lo_R)/A_L)^gama)/lo_R^gama)^2) - (A_L*A_R*gama*((2*U1 - A_R*lo_R)/A_L)^gama*(U2^2/(2*(2*U1 - A_R*lo_R)) - 2*U3 + U2^2/(2*A_R*lo_R)))/(lo_R^gama*(2*U1 - A_R*lo_R)^2*(A_R + (A_L*((2*U1 - A_R*lo_R)/A_L)^gama)/lo_R^gama));
    [x_star, err] = NewtonRapshon(fun,dfun,lo_R,ep);
    lo_R=max(x_star,ep);
    lo_R=min(lo_R,2*U1/A_R);
    k=k+1;
end
if k>=it_max
    err
end
lo_L = (2*U1-A_R*lo_R)/A_L;
p_R  = (2*U3-0.5*U2^2/(2*U1-A_R*lo_R)-0.5*U2^2/A_R/lo_R)*(gama-1)/(A_L/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama+A_R);
p_L  = (2*U3-0.5*U2^2/(2*U1-A_R*lo_R)-0.5*U2^2/A_R/lo_R)*(gama-1)/(A_L/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama+A_R)/lo_R^gama*((2*U1-A_R*lo_R)/A_L)^gama;
if abs(A_L-A_R)>ep
    lo_R
    p_R
    lo_L
    p_L
end
u_L=U2/A_L/lo_L;
u_R=U2/A_R/lo_R;
end
