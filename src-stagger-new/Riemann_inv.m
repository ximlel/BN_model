%compute flux
function [lo_g_n,u_g_n,p_g_n,p_s_n]=Riemann_inv(phi_s,lo_g,u_g,p_g,u_s,p_s,phi_s_n)
global gama_g;
global ep;
    phi_g=1.0-phi_s;
    eta=p_g/lo_g^gama_g;
    Q=phi_g*lo_g*(u_g-u_s);
    P=phi_g*lo_g*(u_g-u_s)^2+phi_g*p_g+phi_s*p_s;
    H=0.5*(u_g-u_s)^2+gama_g/(gama_g-1.0)*p_g/lo_g;
    phi_g_n=1.0-phi_s_n;
    it_max = 500;
    k = 0; err = 1e50;
    lo_g_n = lo_g;
    while (k<it_max && err>ep)
        fun  = H-0.5*(Q/phi_g_n)^2/lo_g_n^2-gama_g/(gama_g-1.0)*eta*lo_g_n^(gama_g-1.0);
        dfun = (Q/phi_g_n)^2/lo_g_n^3-gama_g*eta*lo_g_n^(gama_g-2.0);
        [lo_g_n, err] = NewtonRapshon(fun,dfun,lo_g_n,ep);
        k=k+1;
    end
    if k>=it_max
        err
    end
    p_g_n = lo_g_n^gama_g*eta;
    u_g_n = Q/phi_g_n/lo_g_n+u_s;
    p_s_n = (P-Q*(u_g_n-u_s)-phi_g_n*p_g_n)/phi_s_n;
end
