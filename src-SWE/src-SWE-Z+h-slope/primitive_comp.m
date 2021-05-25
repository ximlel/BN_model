%compute flux
function [h_L,u_L,h_R,u_R,H_t]=primitive_comp(U,Z_L,Z_R)
global g;
global ep;
h=U(1);
q=U(2);
h_R = h;
it_max = 500;
k = 0; err = 1e50;
while (k<it_max && err>ep && abs(Z_L-Z_R)>ep)
    fun = 2*(h-h_R)*(1-q^2*h/g/(2*h-h_R)^2/h_R^2)+Z_L-Z_R;
    dfun= (2*h*q^2)/(g*h_R^2*(2*h - h_R)^2) - ((2*h*q^2)/(g*h_R^2*(2*h - h_R)^3) - (2*h*q^2)/(g*h_R^3*(2*h - h_R)^2))*(2*h - 2*h_R) - 2;
    [x_star, err] = NewtonRapshon(fun,dfun,h_R,ep);
    h_R=max([x_star,ep]);
    h_R=min([2*h-ep,h_R]);
    k=k+1;
end
if k>=it_max
    err
end
h_L = 2*h-h_R;
u_L = q/h_L;
u_R = q/h_R;
H_t = 0.5*((h_L + 0.5*q^2/g/h_L^2 + Z_L)+(h_R + 0.5*q^2/g/h_R^2 + Z_R));
end
