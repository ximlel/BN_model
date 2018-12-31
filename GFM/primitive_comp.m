%compute flux
function [lo,u,p]=primitive_comp(U,gama)
lo = U(1);
u  = U(2)/U(1);
p  = (U(3) - 0.5*lo*u^2)*(gama-1);
end
