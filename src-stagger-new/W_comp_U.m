%compute flux
function U=W_comp_U(W,lo_g_start)
global gama_s gama_g p0;
global ep;
u_s = W(1);
eta_g = W(2);
eta_s = W(3);
Q = W(4);
P = W(5);
H = W(6);
phi_s = W(7);
phi_g = 1.0-phi_s;
it_max = 500;
k = 0; err5 = 1e50;
lo_g=lo_g_start;
N=100;
for i=1:N
    lo_g = 0.1+lo_g_start*i/N*3.0;
    lo_g_plot(i) = lo_g;
    fun  = H-0.5*(Q/phi_g)^2/lo_g^2-gama_g/(gama_g-1.0)*eta_g*lo_g^(gama_g-1.0);
    fun_plot(i) = fun;
end
figure
plot(lo_g_plot,fun_plot);

while (k<it_max && err5>ep)
    fun  = H-0.5*(Q/phi_g)^2/lo_g^2-gama_g/(gama_g-1.0)*eta_g*lo_g^(gama_g-1.0);
    if abs(fun) < ep
        break;
    end
    dfun = (Q/phi_g)^2/lo_g^3-gama_g*eta_g*lo_g^(gama_g-2.0);
    if abs(dfun) < ep
        fun
        break;
    end
    % 零特征值修正
%     if abs(fun/dfun) > 0.00001
%         dfun = dfun*abs((fun/dfun)/0.00001);
%     end
    [lo_g, err5] = NewtonRapshon(fun,dfun,lo_g,ep);
    lo_g=max(lo_g,ep);
    k=k+1;
    err5 = abs(fun);
end
if k>=it_max
    err5
end
p_g = lo_g^gama_g*eta_g;
u_g = Q/phi_g/lo_g+u_s;
p_s = (P-Q*(u_g-u_s)-phi_g*p_g)/phi_s;
lo_s = (p_s/eta_s)^(1/gama_s);
E_g=p_g/(gama_g-1)+0.5*lo_g*u_g^2;
E_s=(p_s+gama_s*p0)/(gama_s-1)+0.5*lo_s*u_s^2;
U=[phi_g*lo_g;phi_g*lo_g*u_g;phi_g*E_g;phi_s*lo_s;phi_s*lo_s*u_s;phi_s*E_s;phi_s];
end