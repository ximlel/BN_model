%comput Roe mean
function [SL,SR,lo_wave,u_wave,a_wave,H_wave,E_wave,p_wave]=Roe_Mean_g(loL,loR,phiL,phiR,pL,pR,uL,uR)
global gama_g;
u_wave=(sqrt(loL*phiL)*uL+sqrt(loR*phiR)*uR)/(sqrt(loL*phiL)+sqrt(loR*phiR));
% u_wave=(sqrt(loL)*uL+sqrt(loR)*uR)/(sqrt(loL)+sqrt(loR));
EL=pL/(gama_g-1)+0.5*loL*uL^2;
ER=pR/(gama_g-1)+0.5*loR*uR^2;
HL=(EL+pL)/loL;
HR=(ER+pR)/loR;
H_wave=(sqrt(loL*phiL)*HL+sqrt(loR*phiR)*HR)/(sqrt(loL*phiL)+sqrt(loR*phiR));
% H_wave=(sqrt(loL)*HL+sqrt(loR)*HR)/(sqrt(loL)+sqrt(loR));
a_wave=sqrt((gama_g-1)*(H_wave-0.5*u_wave^2));
eta_wave=0.5*(pL/loL^gama_g+pR/loR^gama_g);
lo_wave=(a_wave^2/eta_wave/gama_g)^(1.0/(gama_g-1.0));
% lo_wave= sqrt(loL*loR);
E_wave = lo_wave*(H_wave - a_wave^2/gama_g);
p_wave = lo_wave*a_wave^2/gama_g;
SL=u_wave-a_wave;
SR=u_wave+a_wave;
end
