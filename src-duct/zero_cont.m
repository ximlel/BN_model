%comput solid contact
function [S_gL,S_gM,S_gR,S_sL,S_sM,S_sR,lo_g_srL,lo_g_srR,p_g_srL,p_g_srR,u_g_srL,u_g_srR,lo_s_srL,lo_s_srR,p_s_srL,p_s_srR,u_s_srL,u_s_srR,phi_sL,phi_sR] = zero_cont(lo_gL,lo_gR,p_gL,p_gR,u_gL,u_gR,lo_sL,lo_sR,p_sL,p_sR,u_sL,u_sR,phi_sL,phi_sR)
phi_gL = 1.0-phi_sL;
phi_gR = 1.0-phi_sR;
global ep;
global gama_g gama_s p0;
%wave speed estimate
S_sL = u_sL - sqrt(gama_s*(p_sL+p0)/lo_sL);
S_sR = u_sR + sqrt(gama_s*(p_sR+p0)/lo_sR);
S_gL = u_gL - sqrt(gama_g*p_gL/lo_gL);
S_gR = u_gR + sqrt(gama_g*p_gR/lo_gR);
S_s_sr = (p_sR-p_sL+lo_sL*u_sL*(S_sL-u_sL)-lo_sR*u_sR*(S_sR-u_sR))/(lo_sL*(S_sL-u_sL)-lo_sR*(S_sR-u_sR));
S_g_sr = (p_gR-p_gL+lo_gL*u_gL*(S_gL-u_gL)-lo_gR*u_gR*(S_gR-u_gR))/(lo_gL*(S_gL-u_gL)-lo_gR*(S_gR-u_gR));
p_s_srL = max(p_sL+lo_sL*(S_sL-u_sL)*(S_s_sr-u_sL),ep);
p_s_srR = max(p_sR+lo_sR*(S_sR-u_sR)*(S_s_sr-u_sR),ep);
p_g_srL = max(p_gL+lo_gL*(S_gL-u_gL)*(S_g_sr-u_gL),ep);
p_g_srR = max(p_gR+lo_gR*(S_gR-u_gR)*(S_g_sr-u_gR),ep);
if p_s_srL <= p_sL
    S_sL = u_sL - sqrt(gama_s*(p_sL+p0)/lo_sL);
else
    S_sL = u_sL - sqrt(gama_s*(p_sL+p0)/lo_sL)*sqrt(1+(gama_s+1)/2/gama_s*((p_s_srL+p0)/(p_sL+p0)-1));
end
if p_s_srR <= p_sR
    S_sR = u_sR + sqrt(gama_s*(p_sR+p0)/lo_sR);
else
    S_sR = u_sR + sqrt(gama_s*(p_sR+p0)/lo_sR)*sqrt(1+(gama_s+1)/2/gama_s*((p_s_srR+p0)/(p_sR+p0)-1));
end 
if p_g_srL <= p_gL
    S_gL = u_gL - sqrt(gama_g*p_gL/lo_gL);
else
    S_gL = u_gL - sqrt(gama_g*p_gL/lo_gL)*sqrt(1+(gama_g+1)/2/gama_g*(p_g_srL/p_gL-1));
end
if p_g_srR <= p_gR
    S_gR = u_gR + sqrt(gama_g*p_gR/lo_sR);
else
    S_gR = u_gR + sqrt(gama_g*p_gR/lo_gR)*sqrt(1+(gama_g+1)/2/gama_g*(p_g_srR/p_gR-1));
end
S_s_sr = (p_sR-p_sL+lo_sL*u_sL*(S_sL-u_sL)-lo_sR*u_sR*(S_sR-u_sR))/(lo_sL*(S_sL-u_sL)-lo_sR*(S_sR-u_sR));
S_g_sr = (p_gR-p_gL+lo_gL*u_gL*(S_gL-u_gL)-lo_gR*u_gR*(S_gR-u_gR))/(lo_gL*(S_gL-u_gL)-lo_gR*(S_gR-u_gR));
p_s_srL = max(p_sL+lo_sL*(S_sL-u_sL)*(S_s_sr-u_sL),ep);
p_s_srR = max(p_sR+lo_sR*(S_sR-u_sR)*(S_s_sr-u_sR),ep);
p_g_srL = max(p_gL+lo_gL*(S_gL-u_gL)*(S_g_sr-u_gL),ep);
p_g_srR = max(p_gR+lo_gR*(S_gR-u_gR)*(S_g_sr-u_gR),ep);

%Newton-Rapshon iteration for GNL
it_max = 500; 
k = 0; err = 1e50; repe = 0;
fun_L  = zeros(1,4);
fun_R  = zeros(1,4);
%while (abs(phi_sL-phi_sR)>ep && k<it_max && err>ep)
while (k<it_max)
if p_s_srL <= p_sL
    S_sL = u_sL - sqrt(gama_s*(p_sL+p0)/lo_sL);
else
    S_sL = u_sL - sqrt(gama_s*(p_sL+p0)/lo_sL)*sqrt(1+(gama_s+1)/2/gama_s*((p_s_srL+p0)/(p_sL+p0)-1));
end
if p_s_srR <= p_sR
    S_sR = u_sR + sqrt(gama_s*(p_sR+p0)/lo_sR);
else
    S_sR = u_sR + sqrt(gama_s*(p_sR+p0)/lo_sR)*sqrt(1+(gama_s+1)/2/gama_s*((p_s_srR+p0)/(p_sR+p0)-1));
end 
if p_g_srL <= p_gL
    S_gL = u_gL - sqrt(gama_g*p_gL/lo_gL);
else
    S_gL = u_gL - sqrt(gama_g*p_gL/lo_gL)*sqrt(1+(gama_g+1)/2/gama_g*(p_g_srL/p_gL-1));
end
if p_g_srR <= p_gR
    S_gR = u_gR + sqrt(gama_g*p_gR/lo_sR);
else
    S_gR = u_gR + sqrt(gama_g*p_gR/lo_gR)*sqrt(1+(gama_g+1)/2/gama_g*(p_g_srR/p_gR-1));
end
S_s_sr=u_sL + (p_s_srL-p_sL)/lo_sL/(S_sL-u_sL);
S_g_sr=u_gL + (p_g_srL-p_gL)/lo_gL/(S_gL-u_gL);
%S_s_sr=u_sR + (p_s_srR-p_sR)/lo_sR/(S_sR-u_sR);
%S_g_sr=u_gR + (p_g_srR-p_gR)/lo_gR/(S_gR-u_gR);

u_s_srL  = u_sL + (p_s_srL-p_sL)/lo_sL/(S_sL-u_sL);
u_s_srR  = u_sR + (p_s_srR-p_sR)/lo_sR/(S_sR-u_sR);
u_g_srL  = u_gL + (p_g_srL-p_gL)/lo_gL/(S_gL-u_gL);
u_g_srR  = u_gR + (p_g_srR-p_gR)/lo_gR/(S_gR-u_gR);
lo_s_srL = lo_sL^2*(S_sL-u_sL)^2 / (lo_sL*(S_sL-u_sL)^2-p_s_srL+p_sL);
lo_s_srR = lo_sR^2*(S_sR-u_sR)^2 / (lo_sR*(S_sR-u_sR)^2-p_s_srR+p_sR);
lo_g_srL = lo_gL^2*(S_gL-u_gL)^2 / (lo_gL*(S_gL-u_gL)^2-p_g_srL+p_gL);
lo_g_srR = lo_gR^2*(S_gR-u_gR)^2 / (lo_gR*(S_gR-u_gR)^2-p_g_srR+p_gR);
fun_L(1) = u_s_srR-u_s_srL;
fun_L(2) = phi_gR*(p_g_srR/p_g_srL)^(1/gama_g)*(u_g_srR-u_s_srR)-phi_gL*(u_g_srL-u_s_srL);
fun_L(3) = phi_gL*lo_g_srL*(u_g_srL-u_s_srL)*(u_g_srR-u_g_srL)+phi_sR*p_s_srR+phi_gR*p_g_srR-phi_sL*p_s_srL-phi_gL*p_g_srL;
fun_L(4) = gama_g/(gama_g-1)*p_g_srR/lo_g_srL*(p_g_srL/p_g_srR)^(1/gama_g)+0.5*(u_g_srR-u_s_srR)^2-gama_g/(gama_g-1)*p_g_srL/lo_g_srL-0.5*(u_g_srL-u_s_srL)^2;
fun_R(1) = u_s_srR-u_s_srL;
fun_R(2) = phi_gR*(u_g_srR-u_s_srR)-phi_gL*(p_g_srL/p_g_srR)^(1/gama_g)*(u_g_srL-u_s_srL);
fun_R(3) = phi_gR*lo_g_srR*(u_g_srR-u_s_srR)*(u_g_srR-u_g_srL)+phi_sR*p_s_srR+phi_gR*p_g_srR-phi_sL*p_s_srL-phi_gL*p_g_srL;
fun_R(4) = gama_g/(gama_g-1)*p_g_srR/lo_g_srR+0.5*(u_g_srR-u_s_srR)^2-gama_g/(gama_g-1)*p_g_srL/lo_g_srR*(p_g_srR/p_g_srL)^(1/gama_g)-0.5*(u_g_srL-u_s_srL)^2;
dfun_L = [
 -1/(lo_sL*(S_sL - u_sL)),                                                                                                                                                                                             phi_gL/(lo_sL*(S_sL - u_sL)),                                                                                                                                                                                                                                                                                                                                                                                                                               (lo_gL^2*phi_gL*(S_gL - u_gL)^2*(u_gL - u_gR - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR))))/(lo_sL*(S_sL - u_sL)*(p_gL - p_g_srL + lo_gL*(S_gL - u_gL)^2)) - phi_sL,                                                                                                                                                                                                                                                                                                                                                                  (u_gL - u_sL - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_sL - p_s_srL)/(lo_sL*(S_sL - u_sL)))/(lo_sL*(S_sL - u_sL));
  1/(lo_sR*(S_sR - u_sR)),                                                                                                                                                             -(phi_gR*(p_g_srR/p_g_srL)^(1/gama_g))/(lo_sR*(S_sR - u_sR)),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 phi_sR,                                                                                                                                                                                                                                                                                                                                                                 -(u_gR - u_sR - (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR)) + (p_sR - p_s_srR)/(lo_sR*(S_sR - u_sR)))/(lo_sR*(S_sR - u_sR));
                        0,                    - phi_gL/(lo_gL*(S_gL - u_gL)) - (p_g_srR*phi_gR*(p_g_srR/p_g_srL)^(1/gama_g - 1)*(u_gR - u_sR - (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR)) + (p_sR - p_s_srR)/(lo_sR*(S_sR - u_sR))))/(gama_g*p_g_srL^2), - phi_gL - (lo_gL*phi_gL*(S_gL - u_gL)*(u_gL - u_gR - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR))))/(p_gL - p_g_srL + lo_gL*(S_gL - u_gL)^2) - (lo_gL*phi_gL*(S_gL - u_gL)*(u_gL - u_sL - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_sL - p_s_srL)/(lo_sL*(S_sL - u_sL))))/(p_gL - p_g_srL + lo_gL*(S_gL - u_gL)^2) - (lo_gL^2*phi_gL*(S_gL - u_gL)^2*(u_gL - u_gR - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR)))*(u_gL - u_sL - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_sL - p_s_srL)/(lo_sL*(S_sL - u_sL))))/(p_gL - p_g_srL + lo_gL*(S_gL - u_gL)^2)^2, (gama_g*p_g_srL)/(lo_gL^2*(S_gL - u_gL)^2*(gama_g - 1)) - (u_gL - u_sL - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_sL - p_s_srL)/(lo_sL*(S_sL - u_sL)))/(lo_gL*(S_gL - u_gL)) + ((p_g_srL/p_g_srR)^(1/gama_g - 1)*(p_gL - p_g_srL + lo_gL*(S_gL - u_gL)^2))/(lo_gL^2*(S_gL - u_gL)^2*(gama_g - 1)) - (gama_g*(p_gL - p_g_srL + lo_gL*(S_gL - u_gL)^2))/(lo_gL^2*(S_gL - u_gL)^2*(gama_g - 1)) - (gama_g*p_g_srR*(p_g_srL/p_g_srR)^(1/gama_g))/(lo_gL^2*(S_gL - u_gL)^2*(gama_g - 1));
                        0, (phi_gR*(p_g_srR/p_g_srL)^(1/gama_g))/(lo_gR*(S_gR - u_gR)) + (phi_gR*(p_g_srR/p_g_srL)^(1/gama_g - 1)*(u_gR - u_sR - (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR)) + (p_sR - p_s_srR)/(lo_sR*(S_sR - u_sR))))/(gama_g*p_g_srL),                                                                                                                                                                                                                                                                                                                                                                                                                               phi_gR + (lo_gL^2*phi_gL*(S_gL - u_gL)^2*(u_gL - u_sL - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_sL - p_s_srL)/(lo_sL*(S_sL - u_sL))))/(lo_gR*(S_gR - u_gR)*(p_gL - p_g_srL + lo_gL*(S_gL - u_gL)^2)),                                                                                                     (u_gR - u_sR - (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR)) + (p_sR - p_s_srR)/(lo_sR*(S_sR - u_sR)))/(lo_gR*(S_gR - u_gR)) + (gama_g*(p_g_srL/p_g_srR)^(1/gama_g)*(p_gL - p_g_srL + lo_gL*(S_gL - u_gL)^2))/(lo_gL^2*(S_gL - u_gL)^2*(gama_g - 1)) - (p_g_srL*(p_g_srL/p_g_srR)^(1/gama_g - 1)*(p_gL - p_g_srL + lo_gL*(S_gL - u_gL)^2))/(lo_gL^2*p_g_srR*(S_gL - u_gL)^2*(gama_g - 1))];


dfun_R = [
 -1/(lo_sL*(S_sL - u_sL)),                                                                                                                                                                (phi_gL*(p_g_srL/p_g_srR)^(1/gama_g))/(lo_sL*(S_sL - u_sL)),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              -phi_sL,                                                                                                                                                                                                                                                                                                                                                                  (u_gL - u_sL - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_sL - p_s_srL)/(lo_sL*(S_sL - u_sL)))/(lo_sL*(S_sL - u_sL));
  1/(lo_sR*(S_sR - u_sR)),                                                                                                                                                                                              -phi_gR/(lo_sR*(S_sR - u_sR)),                                                                                                                                                                                                                                                                                                                                                                                                                             phi_sR + (lo_gR^2*phi_gR*(S_gR - u_gR)^2*(u_gL - u_gR - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR))))/(lo_sR*(S_sR - u_sR)*(p_gR - p_g_srR + lo_gR*(S_gR - u_gR)^2)),                                                                                                                                                                                                                                                                                                                                                                 -(u_gR - u_sR - (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR)) + (p_sR - p_s_srR)/(lo_sR*(S_sR - u_sR)))/(lo_sR*(S_sR - u_sR));
                        0, - (phi_gL*(p_g_srL/p_g_srR)^(1/gama_g))/(lo_gL*(S_gL - u_gL)) - (phi_gL*(p_g_srL/p_g_srR)^(1/gama_g - 1)*(u_gL - u_sL - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_sL - p_s_srL)/(lo_sL*(S_sL - u_sL))))/(gama_g*p_g_srR),                                                                                                                                                                                                                                                                                                                                                                                                                           - phi_gL - (lo_gR^2*phi_gR*(S_gR - u_gR)^2*(u_gR - u_sR - (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR)) + (p_sR - p_s_srR)/(lo_sR*(S_sR - u_sR))))/(lo_gL*(S_gL - u_gL)*(p_gR - p_g_srR + lo_gR*(S_gR - u_gR)^2)),                                                                                                     (p_g_srR*(p_g_srR/p_g_srL)^(1/gama_g - 1)*(p_gR - p_g_srR + lo_gR*(S_gR - u_gR)^2))/(lo_gR^2*p_g_srL*(S_gR - u_gR)^2*(gama_g - 1)) - (gama_g*(p_g_srR/p_g_srL)^(1/gama_g)*(p_gR - p_g_srR + lo_gR*(S_gR - u_gR)^2))/(lo_gR^2*(S_gR - u_gR)^2*(gama_g - 1)) - (u_gL - u_sL - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_sL - p_s_srL)/(lo_sL*(S_sL - u_sL)))/(lo_gL*(S_gL - u_gL));
                        0,                        phi_gR/(lo_gR*(S_gR - u_gR)) + (p_g_srL*phi_gL*(p_g_srL/p_g_srR)^(1/gama_g - 1)*(u_gL - u_sL - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_sL - p_s_srL)/(lo_sL*(S_sL - u_sL))))/(gama_g*p_g_srR^2), phi_gR - (lo_gR*phi_gR*(S_gR - u_gR)*(u_gL - u_gR - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR))))/(p_gR - p_g_srR + lo_gR*(S_gR - u_gR)^2) + (lo_gR*phi_gR*(S_gR - u_gR)*(u_gR - u_sR - (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR)) + (p_sR - p_s_srR)/(lo_sR*(S_sR - u_sR))))/(p_gR - p_g_srR + lo_gR*(S_gR - u_gR)^2) - (lo_gR^2*phi_gR*(S_gR - u_gR)^2*(u_gL - u_gR - (p_gL - p_g_srL)/(lo_gL*(S_gL - u_gL)) + (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR)))*(u_gR - u_sR - (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR)) + (p_sR - p_s_srR)/(lo_sR*(S_sR - u_sR))))/(p_gR - p_g_srR + lo_gR*(S_gR - u_gR)^2)^2, (u_gR - u_sR - (p_gR - p_g_srR)/(lo_gR*(S_gR - u_gR)) + (p_sR - p_s_srR)/(lo_sR*(S_sR - u_sR)))/(lo_gR*(S_gR - u_gR)) - (gama_g*p_g_srR)/(lo_gR^2*(S_gR - u_gR)^2*(gama_g - 1)) - ((p_g_srR/p_g_srL)^(1/gama_g - 1)*(p_gR - p_g_srR + lo_gR*(S_gR - u_gR)^2))/(lo_gR^2*(S_gR - u_gR)^2*(gama_g - 1)) + (gama_g*(p_gR - p_g_srR + lo_gR*(S_gR - u_gR)^2))/(lo_gR^2*(S_gR - u_gR)^2*(gama_g - 1)) + (gama_g*p_g_srL*(p_g_srR/p_g_srL)^(1/gama_g))/(lo_gR^2*(S_gR - u_gR)^2*(gama_g - 1))];

% if S_s_sr < S_gL || S_s_sr > S_gR
%     error('Solid contact is not in the gas!');
% else
if S_g_sr >= S_s_sr
    if norm(fun_L,inf) <= ep
        break;
    end
    [x_star, err] = NewtonRapshon(fun_L,dfun_L,[p_s_srL p_s_srR p_g_srL p_g_srR],ep);
else
    if norm(fun_R,inf) <= ep
        break;
    end
    [x_star, err] = NewtonRapshon(fun_R,dfun_R,[p_s_srL p_s_srR p_g_srL p_g_srR],ep);
end
p_s_srL=max(x_star(1),ep);
p_s_srR=max(x_star(2),ep);
p_g_srL=max(x_star(3),ep);
p_g_srR=max(x_star(4),ep);
k=k+1;
if k>=it_max && repe == 0
    p_s_srL = p_sL;
    p_s_srR = p_sR;
    p_g_srL = p_gL;
    p_g_srR = p_gR;
    repe = 1;
    k = 0;
end
end
if k>=it_max 
    err
end

%compute lo, u, p, S
u_s_srL  = u_sL + (p_s_srL-p_sL)/lo_sL/(S_sL-u_sL);
u_s_srR  = u_sR + (p_s_srR-p_sR)/lo_sR/(S_sR-u_sR);
u_g_srL  = u_gL + (p_g_srL-p_gL)/lo_gL/(S_gL-u_gL);
u_g_srR  = u_gR + (p_g_srR-p_gR)/lo_gR/(S_gR-u_gR);
lo_s_srL = lo_sL^2*(S_sL-u_sL)^2 / (lo_sL*(S_sL-u_sL)^2-p_s_srL+p_sL);
lo_s_srR = lo_sR^2*(S_sR-u_sR)^2 / (lo_sR*(S_sR-u_sR)^2-p_s_srR+p_sR);
lo_g_srL = lo_gL^2*(S_gL-u_gL)^2 / (lo_gL*(S_gL-u_gL)^2-p_g_srL+p_gL);
lo_g_srR = lo_gR^2*(S_gR-u_gR)^2 / (lo_gR*(S_gR-u_gR)^2-p_g_srR+p_gR);
if u_g_srL >= u_s_srL
    S_gM = u_g_srR;
    S_sM = u_s_srL;
else
    S_gM = u_g_srL;
    S_sM = u_s_srR;
end
end
