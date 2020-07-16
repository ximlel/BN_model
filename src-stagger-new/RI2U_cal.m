%compute flux
function [lo_g,u_g,p_g,u_s,p_s,lo_s]=RI2U_cal(phi_s,RI,lo_g_start)
global gama_g;
global ep;
    phi_g=1.0-phi_s;
    lo_s=RI(1);
    u_s=RI(2);
    P=RI(3);
    Q=RI(4);
    H=RI(5);
    eta_g=RI(6);
    it_max = 500;
    k = 0; err1 = 1e50;
    lo_g=lo_g_start;
    while (k<it_max && err1>ep)
        fun  = H-0.5*(Q/phi_g)^2/lo_g^2-gama_g/(gama_g-1.0)*eta_g*lo_g^(gama_g-1.0);
        if abs(fun) < ep
            break;
        end
        dfun = (Q/phi_g)^2/lo_g^3-gama_g*eta_g*lo_g^(gama_g-2.0);
        if dfun^2 < ep
            dfun = ep/dfun;
            break;
        end
        if abs(dfun) < ep
            fun
            break;
        end
        [lo_g, err1] = NewtonRapshon(fun,dfun,lo_g,ep);
        lo_g = max(lo_g,ep);
        k=k+1;
        err1 = abs(fun);
    end
    if k>=it_max
        err1

% N=100;
% for i=1:N
%     lo_g_c = 0.1+lo_g_start*i/N*3.0;
%     lo_g_plot(i) = lo_g_c;
%     fun  = H-0.5*(Q/phi_g)^2/lo_g_c^2-gama_g/(gama_g-1.0)*eta_g*lo_g_c^(gama_g-1.0);
%     fun_plot(i) = fun;
% end
% figure
% plot(lo_g_plot,fun_plot);
  if err1>1e-2
    k = 0; err2 = 1e50;
    lo_g=lo_g_start;
    while (k<it_max && err2>ep)
        fun  = H-0.5*(Q/phi_g)^2/lo_g^2-gama_g/(gama_g-1.0)*eta_g*lo_g^(gama_g-1.0);
        dfun = (Q/phi_g)^2/lo_g^3-gama_g*eta_g*lo_g^(gama_g-2.0);
        if abs(dfun) < ep
            fun
            break;
        end
        ddfun = -3.0*(Q/phi_g)^2/lo_g^4-gama_g*(gama_g-2.0)*eta_g*lo_g^(gama_g-3.0);
        [lo_g, err2] = NewtonRapshon(fun*dfun,dfun^2+fun*ddfun,lo_g,ep);
        lo_g = max(lo_g,ep);
        k=k+1;
        err2 = abs(dfun);
    end
    err2
  end
    end
    p_g = lo_g^gama_g*eta_g;
    u_g = Q/phi_g/lo_g+u_s;
    p_s = (P-Q*(u_g-u_s)-phi_g*p_g)/phi_s;  
    
% if abs(lo_g_start-6.3311) > 0.001 && abs(lo_g_start-0.4141) > 0.001
%     lo_g,u_g,p_g,u_s,p_s,lo_s
%     phi_s
% end

% if lo_g < 0.01
%     lo_g,u_g,p_g,u_s,p_s,lo_s
%     phi_s
%     N=100;
% for i=1:N
%     lo_g_c = 0.1+lo_g_start*i/N*3.0;
%     lo_g_plot(i) = lo_g_c;
%     fun  = H-0.5*(Q/phi_g)^2/lo_g_c^2-gama_g/(gama_g-1.0)*eta_g*lo_g_c^(gama_g-1.0);
%     fun_plot(i) = fun;
% end
% figure
% plot(lo_g_plot,fun_plot);
% end
end
