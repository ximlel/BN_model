%Riemann Solver with Roe scheme
function out_flux=Riemann_solver_HLL(hL,uL,hR,uR)
%state constant
global g;
%comput Roe mean
aL=sqrt(g*hL);
aR=sqrt(g*hR);
%solve averaged eigenvalues
SL=min([uL-aL,uR-aR]);
SR=max([uL+aL,uR+aR]);
UL=[hL;hL*uL];
UR=[hR;hR*uR];
FL=[hL*uL;hL*uL^2+g*hL^2/2];
FR=[hR*uR;hR*uR^2+g*hR^2/2];
%entropy fix
FHLL=(SR*FL-SL*FR+SL*SR*(UR-UL))/(SR-SL);
if SL>=0
    out_flux=FL;
elseif SR<=0
    out_flux=FR;
else
    out_flux=FHLL;
end
