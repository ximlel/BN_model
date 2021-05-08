%compute flux
function [h,u,dh,du]=dRI2dU_cal(q,H_s,dq,dH_t,dZ,F)
global g;
global ep;
[h,u]=RI2U_cal(q,H_s,F);
if abs(F-1) > ep
    du = (u*dH_t -     dq - u*dZ)/(u^2/g - h);
    dh = (h*dH_t - u/g*dq - h*dZ)/(h - u^2/g);
else
    error('F=1');
end
end
