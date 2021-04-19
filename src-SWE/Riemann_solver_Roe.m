%Riemann Solver with Roe scheme
function out_flux=Riemann_solver_Roe(loL,loR,pL,pR,uL,uR,A)
%state constant
global gama;
%comput Roe mean
[SL SR u_wave a_wave H_wave]=Roe_Mean(loL,loR,pL,pR,uL,uR);
%solve averaged eigenvalues
lamda1=u_wave-a_wave;
lamda2=u_wave;
lamda3=u_wave+a_wave;
lamda=[lamda1 lamda2 lamda3];
%entropy fix
global ep;
abs_lamda=zeros(3,1);
for rr=1:3
    if abs(lamda(rr))>=ep
        abs_lamda(rr)=abs(lamda(rr));
    else
        abs_lamda(rr)=(lamda(rr)^2+ep^2)/2/ep;
    end
end
%solve right eigenvetors
K1=[1;u_wave-a_wave;H_wave-u_wave*a_wave];
K2=[1;u_wave;1/2*u_wave^2];
K5=[1;u_wave+a_wave;H_wave+u_wave*a_wave];
K=[K1 K2 K5];
%solve wave strengths
d_u=zeros(3,1);
EL=pL/(gama-1)+1/2*loL*uL^2;
ER=pR/(gama-1)+1/2*loR*uR^2;
UL=[loL;loL*uL;EL];
UR=[loR;loR*uR;ER];
d_u=UR-UL;
alpha2=(gama-1)/a_wave^2*(d_u(1)*(H_wave-u_wave^2)+u_wave*d_u(2)-d_u(3));
alpha1=1/2/a_wave*(d_u(1)*(u_wave+a_wave)-d_u(2)-a_wave*alpha2);
alpha5=d_u(1)-alpha1-alpha2;
alpha=[alpha1;alpha2;alpha5];
%solve flux at i+1/2
FL=[loL*uL;loL*uL^2+pL;(EL+pL)*uL];
FR=[loR*uR;loR*uR^2+pR;(ER+pR)*uR];
out_flux=1/2*(FL+FR)-1/2*K*(alpha.*abs_lamda);
out_flux=A*out_flux;
end
