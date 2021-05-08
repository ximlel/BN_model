%compute flux
function [h_out,u,dh_out,du]=dRI2dU_cal(q,h,dq,dh,dZ,F)
u = q/h;
du = (dq - u*dh)/h;
h_out = h;
dh_out = dh;
end
