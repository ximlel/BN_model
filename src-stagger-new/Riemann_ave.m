%compute flux
function [lo_gL_n,u_gL_n,p_gL_n,lo_sL_n,u_sL_n,p_sL_n,lo_gR_n,u_gR_n,p_gR_n,lo_sR_n,u_sR_n,p_sR_n]=Riemann_ave(lo_gL,u_gL,p_gL,lo_sL,u_sL,p_sL,phi_sL,lo_gR,u_gR,p_gR,lo_sR,u_sR,p_sR,phi_sR)
global gama_g;
global ep;
    phi_gL=1.0-phi_sL;
    phi_gR=1.0-phi_sR;
    lo_s=(phi_sL*lo_sL+phi_sR*lo_sR)/(phi_sL+phi_sR);
    lo_sL_n=lo_s;
    lo_sR_n=lo_s;
    u_s=(phi_sL*lo_sL*u_sL+phi_sR*lo_sR*u_sR)/(phi_sL*lo_sL+phi_sR*lo_sR);
    u_sL_n=u_s;
    u_sR_n=u_s;
    etaL=p_gL/lo_gL^gama_g;
    QL=phi_gL*lo_gL*(u_gL-u_s);
    PL=phi_gL*lo_gL*(u_gL-u_s)^2+phi_gL*p_gL+phi_sL*p_sL;
    HL=0.5*(u_gL-u_s)^2+gama_g/(gama_g-1.0)*p_gL/lo_gL;
    etaR=p_gR/lo_gR^gama_g;
    QR=phi_gR*lo_gR*(u_gR-u_s);
    PR=phi_gR*lo_gR*(u_gR-u_s)^2+phi_gR*p_gR+phi_sR*p_sR;
    HR=0.5*(u_gR-u_s)^2+gama_g/(gama_g-1.0)*p_gR/lo_gR;
    eta=0.5*(etaL+etaR);
    eta=(phi_sL*lo_sL^gama_g*etaL+phi_sR*lo_sR^gama_g*etaR)/(phi_sL*lo_sL^gama_g+phi_sR*lo_sR^gama_g);
    Q=0.5*(QL+QR);
    P=0.5*(PL+PR);
    H=0.5*(HL+HR);
    H=(phi_sL*lo_sL*HL+phi_sR*lo_sR*HR)/(phi_sL*lo_sL+phi_sR*lo_sR);
    it_N =2500;
    it_max = 2*it_N; it_max1 = it_N;
    k = 0; errA1 = 1e50;
    lo_gL_n = lo_gL;
    while (k<it_max && errA1>ep)
        if (k == it_max1)
            lo_gL_n = lo_gR;
        end
        fun  = H-0.5*(Q/phi_gL)^2/lo_gL_n^2-gama_g/(gama_g-1.0)*eta*lo_gL_n^(gama_g-1.0);
        if abs(fun) < ep
            break;
        end
        dfun = (Q/phi_gL)^2/lo_gL_n^3-gama_g*eta*lo_gL_n^(gama_g-2.0);
        if abs(dfun) < ep
            fun
            break;
        end
        % 零特征值修正
        if abs(fun/dfun) > abs(lo_gL-lo_gR)
            dfun = dfun*abs((fun/dfun)/(lo_gL-lo_gR));
        end
        [lo_gL_n, errA1] = NewtonRapshon(fun,dfun,lo_gL_n,ep);
        lo_gL_n=max(lo_gL_n,ep);
        k=k+1;
        errA1 = abs(fun);
    end
    if k>=it_max
        errA1
    end
    p_gL_n = lo_gL_n^gama_g*eta;
    u_gL_n = Q/phi_gL/lo_gL_n+u_s;
    p_sL_n = (P-Q*(u_gL_n-u_s)-phi_gL*p_gL_n)/phi_sL;
    
    k = 0; errA2 = 1e50;
    lo_gR_n = lo_gR;
    while (k<it_max && errA2>ep)
        if (k == it_max1)
            lo_gR_n = lo_gL;
        end
        fun  = H-0.5*(Q/phi_gR)^2/lo_gR_n^2-gama_g/(gama_g-1.0)*eta*lo_gR_n^(gama_g-1.0);
        if abs(fun) < ep
            break;
        end
        dfun = (Q/phi_gR)^2/lo_gR_n^3-gama_g*eta*lo_gR_n^(gama_g-2.0);
        if abs(dfun) < ep
            fun
            break;
        end
        % 零特征值修正
        if abs(fun/dfun) > abs(lo_gL-lo_gR)
            dfun = dfun*abs((fun/dfun)/(lo_gL-lo_gR));
        end
        [lo_gR_n, errA2] = NewtonRapshon(fun,dfun,lo_gR_n,ep);
        lo_gR_n=max(lo_gR_n,ep);
        k=k+1;
        errA2 = abs(fun);
    end
    if k>=it_max
        errA2
    end
    p_gR_n = lo_gR_n^gama_g*eta;
    u_gR_n = Q/phi_gR/lo_gR_n+u_s;
    p_sR_n = (P-Q*(u_gR_n-u_s)-phi_gR*p_gR_n)/phi_sR;
end
