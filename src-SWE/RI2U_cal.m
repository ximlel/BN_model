%compute flux
function [lo,u,p]=RI2U_cal(phi,RI,lo_start)
global g;
global ep;
    Q=RI(1);
    eta=RI(2);
    H=RI(3);
    it_max = 5000;
    k = 0; err1 = 1e50;
    lo=lo_start;
    while (k<it_max && err1>ep)
        fun  = H-0.5*(Q/phi)^2/lo^2-gama/(gama-1.0)*eta*lo^(gama-1.0);
        dfun = (Q/phi)^2/lo^3-gama*eta*lo^(gama-2.0);
        [lo, err1] = NewtonRapshon(fun,dfun,lo,ep);
        k=k+1;
    end
    if k>=it_max
        err1
    end
    p = lo^gama*eta;
    u = Q/phi/lo;
end
