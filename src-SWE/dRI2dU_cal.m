%compute flux
function [h,u,dh,du]=dRI2dU_cal(q,H_s,dq,dH_t,dZ,F)
global g;
global ep;
[h,u]=RI2U_cal(q,H_s,F);
dH_s = dH_t - dZ;
if abs(F-1) > ep
    du = (u*dH_s -     dq)/(u^2/g - h);
    dh = (h*dH_s - u/g*dq)/(h - u^2/g);
else
    error('F=1');
end
end