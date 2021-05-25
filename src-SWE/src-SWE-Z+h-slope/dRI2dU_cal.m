%compute flux
function [h_out,u,dh,du]=dRI2dU_cal(q,h,dq,dzeta,dZ,F)
u = q/h;
dh = dzeta - dZ;
du = (dq - u*dh)/h;
h_out = h;
end
