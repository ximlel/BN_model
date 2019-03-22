%compute flux
function [lo_g,u_g,p_g,u_s,p_s,lo_s]=RI2U_cal(phi_s,RI,lo_g_start)
global gama_g;
global ep;
    phi_g=1.0-phi_s;
    lo_s=RI(1);
    u_s=RI(2);
    Q=RI(3);
    P=RI(4);
    H=RI(5);
    eta_g=RI(6);
    it_max = 500;
    k = 0; err = 1e50;
    lo_g=lo_g_start;
    while (k<it_max && err>ep)
        fun  = H-0.5*(Q/phi_g)^2/lo_g^2-gama_g/(gama_g-1.0)*eta_g*lo_g^(gama_g-1.0);
        dfun = (Q/phi_g)^2/lo_g^3-gama_g*eta_g*lo_g^(gama_g-2.0);
        [lo_g, err] = NewtonRapshon(fun,dfun,lo_g,ep);
        k=k+1;
    end
    if k>=it_max
        err
    end
    p_g = lo_g^gama_g*eta_g;
    u_g = Q/phi_g/lo_g+u_s;
    p_s = (P-Q*(u_g-u_s)-phi_g*p_g)/phi_s;
end
