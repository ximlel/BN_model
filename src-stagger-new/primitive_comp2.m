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
p_gR  = (U3/phi_g - 0.5*lo_gR*(U2/U1)^2)*(gama_g-1);

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
% Least Squares
x_k = [lo_gR;p_gR];
% Gauss-Newton
fun  = ones(1,2);
dfun = zeros(2,2);
it_N = 20;
it_max_LS =500; it_max = 7*it_N;
it_max1 = it_N; it_max2 = 2*it_N; it_max3 = 3*it_N; it_max4 = 4*it_N; it_max5 = 5*it_N; it_max6 = 6*it_N;
k = 0; l = 0; err2 = 1e10; err_LS = 1e10;
use_LS = 0;
while (k<it_max && err2>ep && abs(phi_sL-phi_sR)>ep)
    if (k==it_max1)
        lo_gR = U1/phi_g;
        p_gR=(U3/phi_gL - 0.5*lo_gR*(U2/U1)^2)*(gama_g-1);
    elseif (k==it_max2)
        lo_gR = U1/phi_g;
        p_gR=(U3/phi_gR - 0.5*lo_gR*(U2/U1)^2)*(gama_g-1);
    elseif (k==it_max3)
        lo_gR = U1/phi_gL;
        p_gR=(U3/phi_gL - 0.5*lo_gR*(U2/U1)^2)*(gama_g-1);
    elseif (k==it_max4)
        lo_gR = U1/phi_gR;
        p_gR=(U3/phi_gR - 0.5*lo_gR*(U2/U1)^2)*(gama_g-1);
    elseif (k==it_max5)
        lo_gR = U1/phi_gR;
        p_gR=(U3/phi_gL - 0.5*lo_gR*(U2/U1)^2)*(gama_g-1);
    elseif (k==it_max6)
        lo_gR = U1/phi_gL;
        p_gR=(U3/phi_gR - 0.5*lo_gR*(U2/U1)^2)*(gama_g-1);
    end
    fun(1) = ((U3 + (-0.1e1) * 0.5e0 * area_R * phi_gR * lo_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + (-0.1e1) * 0.5e0 * (-area_R * lo_gR * phi_gR + U1) * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2) * (gama_g - 1) - area_R * phi_gR * p_gR) / area_L / phi_gL / ((-area_R * lo_gR * phi_gR + U1) / area_L / phi_gL) ^ gama_g - p_gR / lo_gR ^ gama_g;
    fun(2) = 0.5e0 * (-U1 * u_s + U2) ^ 2 * area_L ^ 2 / (-area_R * lo_gR * phi_gR + U1) ^ 2 + gama_g * ((U3 + (-0.1e1) * 0.5e0 * area_R * phi_gR * lo_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + (-0.1e1) * 0.5e0 * (-area_R * lo_gR * phi_gR + U1) * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2) * (gama_g - 1) - area_R * phi_gR * p_gR) / (gama_g - 1) / (-area_R * lo_gR * phi_gR + U1) + (-0.1e1) * 0.5e0 * (-U1 * u_s + U2) ^ 2 / phi_gR ^ 2 / lo_gR ^ 2 - gama_g * p_gR / (gama_g - 1) / lo_gR;
    if norm(fun)<ep
       break; 
    end
    dfun(1,1) = ((-0.1e1) * 0.5e0 * area_R * phi_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + 0.10e1 * area_R * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) * (-U1 * u_s + U2) / lo_gR + 0.5e0 * area_R * phi_gR * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2 + (-0.1e1) * 0.10e1 * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) * (-U1 * u_s + U2) * area_L * area_R * phi_gR / (-area_R * lo_gR * phi_gR + U1)) * (gama_g - 1) / area_L / phi_gL / ((-area_R * lo_gR * phi_gR + U1) / area_L / phi_gL) ^ gama_g + ((U3 + (-0.1e1) * 0.5e0 * area_R * phi_gR * lo_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + (-0.1e1) * 0.5e0 * (-area_R * lo_gR * phi_gR + U1) * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2) * (gama_g - 1) - area_R * phi_gR * p_gR) * gama_g * area_R * phi_gR / area_L / phi_gL / ((-area_R * lo_gR * phi_gR + U1) / area_L / phi_gL) ^ gama_g / (-area_R * lo_gR * phi_gR + U1) + p_gR * gama_g / lo_gR ^ gama_g / lo_gR;
    dfun(1,2) = -area_R * phi_gR / area_L / phi_gL / ((-area_R * lo_gR * phi_gR + U1) / area_L / phi_gL) ^ gama_g - 0.1e1 / lo_gR ^ gama_g;
    dfun(2,1) = 0.10e1 * (-U1 * u_s + U2) ^ 2 * area_L ^ 2 * area_R * phi_gR / (-area_R * lo_gR * phi_gR + U1) ^ 3 + gama_g * ((-0.1e1) * 0.5e0 * area_R * phi_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + 0.10e1 * area_R * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) * (-U1 * u_s + U2) / lo_gR + 0.5e0 * area_R * phi_gR * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2 + (-0.1e1) * 0.10e1 * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) * (-U1 * u_s + U2) * area_L * area_R * phi_gR / (-area_R * lo_gR * phi_gR + U1)) / (-area_R * lo_gR * phi_gR + U1) + gama_g * ((U3 + (-0.1e1) * 0.5e0 * area_R * phi_gR * lo_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + (-0.1e1) * 0.5e0 * (-area_R * lo_gR * phi_gR + U1) * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2) * (gama_g - 0.1e1) - area_R * phi_gR * p_gR) * area_R * phi_gR / (gama_g - 0.1e1) / (-area_R * lo_gR * phi_gR + U1) ^ 2 + 0.10e1 * (-U1 * u_s + U2) ^ 2 / phi_gR ^ 2 / lo_gR ^ 3 + gama_g * p_gR / (gama_g - 0.1e1) / lo_gR ^ 2;
    dfun(2,2) = -gama_g * area_R * phi_gR / (gama_g - 1) / (-area_R * lo_gR * phi_gR + U1) - gama_g / (gama_g - 1) / lo_gR;
    [x_star, err2] = NewtonRapshon(fun,dfun',[lo_gR p_gR],ep);
%     if real(x_star(1)) < 1e-6 || real(x_star(1)) > U1/area_R/phi_gR-1e-6 || real(x_star(2)) < 1e-6
%        k = it_max;
%        break; 
%     end
    lo_gR=max(real(x_star(1)),1e-6);
    lo_gR=min(lo_gR,U1/area_R/phi_gR-1e-6);
    p_gR =max(real(x_star(2)),1e-6);
    k=k+1;
end
if abs(phi_sL-phi_sR)>ep && k>=it_max
    err2
    lo_gR
    p_gR
%     FF11=fun(1)
%     FF12=fun(2)
    
    f=@(lo_gR,p_gR) ((U3 + (-0.1e1) * 0.5e0 * area_R * phi_gR * lo_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + (-0.1e1) * 0.5e0 * (-area_R * lo_gR * phi_gR + U1) * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2) * (gama_g - 1) - area_R * phi_gR * p_gR) / area_L / phi_gL / ((-area_R * lo_gR * phi_gR + U1) / area_L / phi_gL) ^ gama_g - p_gR / lo_gR ^ gama_g;
    h=@(lo_gR,p_gR) 0.5e0 * (-U1 * u_s + U2) ^ 2 * area_L ^ 2 / (-area_R * lo_gR * phi_gR + U1) ^ 2 + gama_g * ((U3 + (-0.1e1) * 0.5e0 * area_R * phi_gR * lo_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + (-0.1e1) * 0.5e0 * (-area_R * lo_gR * phi_gR + U1) * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2) * (gama_g - 1) - area_R * phi_gR * p_gR) / (gama_g - 1) / (-area_R * lo_gR * phi_gR + U1) + (-0.1e1) * 0.5e0 * (-U1 * u_s + U2) ^ 2 / phi_gR ^ 2 / lo_gR ^ 2 - gama_g * p_gR / (gama_g - 1) / lo_gR;
    Dxf=@(lo_gR,p_gR) [((-0.1e1) * 0.5e0 * area_R * phi_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + 0.10e1 * area_R * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) * (-U1 * u_s + U2) / lo_gR + 0.5e0 * area_R * phi_gR * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2 + (-0.1e1) * 0.10e1 * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) * (-U1 * u_s + U2) * area_L * area_R * phi_gR / (-area_R * lo_gR * phi_gR + U1)) * (gama_g - 1) / area_L / phi_gL / ((-area_R * lo_gR * phi_gR + U1) / area_L / phi_gL) ^ gama_g + ((U3 + (-0.1e1) * 0.5e0 * area_R * phi_gR * lo_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + (-0.1e1) * 0.5e0 * (-area_R * lo_gR * phi_gR + U1) * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2) * (gama_g - 1) - area_R * phi_gR * p_gR) * gama_g * area_R * phi_gR / area_L / phi_gL / ((-area_R * lo_gR * phi_gR + U1) / area_L / phi_gL) ^ gama_g / (-area_R * lo_gR * phi_gR + U1) + p_gR * gama_g / lo_gR ^ gama_g / lo_gR;
                           -area_R * phi_gR / area_L / phi_gL / ((-area_R * lo_gR * phi_gR + U1) / area_L / phi_gL) ^ gama_g - 0.1e1 / lo_gR ^ gama_g];
    Dxh=@(lo_gR,p_gR) [0.10e1 * (-U1 * u_s + U2) ^ 2 * area_L ^ 2 * area_R * phi_gR / (-area_R * lo_gR * phi_gR + U1) ^ 3 + gama_g * ((-0.1e1) * 0.5e0 * area_R * phi_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + 0.10e1 * area_R * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) * (-U1 * u_s + U2) / lo_gR + 0.5e0 * area_R * phi_gR * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2 + (-0.1e1) * 0.10e1 * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) * (-U1 * u_s + U2) * area_L * area_R * phi_gR / (-area_R * lo_gR * phi_gR + U1)) / (-area_R * lo_gR * phi_gR + U1) + gama_g * ((U3 + (-0.1e1) * 0.5e0 * area_R * phi_gR * lo_gR * ((-U1 * u_s + U2) / phi_gR / lo_gR + u_s) ^ 2 + (-0.1e1) * 0.5e0 * (-area_R * lo_gR * phi_gR + U1) * ((-U1 * u_s + U2) * area_L / (-area_R * lo_gR * phi_gR + U1) + u_s) ^ 2) * (gama_g - 0.1e1) - area_R * phi_gR * p_gR) * area_R * phi_gR / (gama_g - 0.1e1) / (-area_R * lo_gR * phi_gR + U1) ^ 2 + 0.10e1 * (-U1 * u_s + U2) ^ 2 / phi_gR ^ 2 / lo_gR ^ 3 + gama_g * p_gR / (gama_g - 0.1e1) / lo_gR ^ 2;
                       -gama_g * area_R * phi_gR / (gama_g - 1) / (-area_R * lo_gR * phi_gR + U1) - gama_g / (gama_g - 1) / lo_gR];
    while (l<it_max_LS && err_LS >ep)
        if sqrt(f(x_k(1),x_k(2))^2+h(x_k(1),x_k(2))^2)<ep
            break; 
        end
        Dg_Dg = Dxf(x_k(1),x_k(2))*Dxf(x_k(1),x_k(2))'+Dxh(x_k(1),x_k(2))*Dxh(x_k(1),x_k(2))';
        [L, DMC, P, D] = modchol_ldlt(Dg_Dg,1e-4*norm(Dg_Dg,'fro'));
        LL = inv(P)*L*DMC*L'*inv(P');
        Dg_g = Dxf(x_k(1),x_k(2))*f(x_k(1),x_k(2))+Dxh(x_k(1),x_k(2))*h(x_k(1),x_k(2));
        d_k = -LL\Dg_g;
        err_LS=norm(d_k);
        x_k = x_k+d_k;
        x_k(1)=max(x_k(1),1e-6);
        x_k(1)=min(x_k(1),U1/area_R/phi_gR-1e-6);
        x_k(2)=max(x_k(2),1e-6);
        l=l+1;
    end
%     if l>=it_max_LS
        err_LS
%     end
    use_LS = 1;
end
Q = U2 - U1*u_s;
lo_gL= (U1 - area_R*phi_gR*lo_gR)/area_L/phi_gL;
u_gL = Q/(phi_gL*lo_gL) + u_s;
u_gR = Q/(phi_gR*lo_gR) + u_s;
p_gL = ((U3 - 0.5*area_R*phi_gR*lo_gR*u_gR^2 - 0.5*area_L*phi_gL*lo_gL*u_gL^2)*(gama_g - 1) - area_R*phi_gR*p_gR)/(area_L*phi_gL);
p_sR = (((phi_gL * lo_gL + (-0.1e1) * 0.1e1 * phi_gR * lo_gR) * u_s ^ 2 + ((-0.1e1) * 0.2e1 * lo_gL * phi_gL * u_gL + 0.2e1 * lo_gR * phi_gR * u_gR) * u_s + (lo_gL * u_gL ^ 2 + p_gL) * phi_gL + ((-0.1e1) * 0.1e1 * lo_gR * u_gR ^ 2 + (-0.1e1) * 0.1e1 * p_gR) * phi_gR) * area_L + ((-0.1e1) * 0.5e0 * gama_s * lo_s * phi_s + 0.5e0 * phi_s * lo_s) * u_s ^ 2 + U6 * gama_s + (-0.1e1) * 0.1e1 * U6) / phi_sR / (area_L + area_R);
p_sL = ((U6 - 0.5*phi_s*lo_s*u_s^2)*(gama_s - 1) - area_R*phi_sR*p_sR)/(area_L*phi_sL);
lo_sL= lo_s;
lo_sR= lo_s;
% lo_sL= U_lo_sL/phi_sL;
% lo_sR= U_lo_sR/phi_sR;
u_sL = u_s;
u_sR = u_s;
% if use_LS == 1 && 
if abs(phi_sL-phi_sR)>ep
    if lo_gL<0 || lo_sL<0 || p_gL<0 || p_sL<0
    [lo_gL,u_gL,p_gL,lo_sL,u_sL,p_sL]
    elseif  lo_gR<0 || lo_sR<0 || p_gR<0 || p_sR<0
    [lo_gR,u_gR,p_gR,lo_sR,u_sR,p_sR]
    end
    [lo_gL,u_gL,p_gL,lo_sL,u_sL,p_sL,lo_gR,u_gR,p_gR,lo_sR,u_sR,p_sR]=Riemann_ave(lo_gL,u_gL,p_gL,lo_sL,u_sL,p_sL,phi_sL,lo_gR,u_gR,p_gR,lo_sR,u_sR,p_sR,phi_sR);
end
end


