%compute flux
function [h,u]=RI2U_cal(q,H_s,F)
global g;
global ep;
alpha = 0.5*q^2/g/H_s^3;
p=[1, -1, 0, alpha];
if alpha > 0 && alpha < 4/27
R=roots(p);
    if F > 1 % supercritical
        x = R(find(R>0 & R<2/3));
    elseif F < 1 % subcritical
        x = R(find(R>2/3 & R<1));
    end
elseif abs(alpha)<ep % zero velocity
    x = 0;
elseif abs(alpha - 4/27)<ep % critical
    x = 2/3;
else
    alpha
    error('alpha error');
end
h = H_s*x;
u = q/h;
end
