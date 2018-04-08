%comput Roe mean
function [SL SR lo_wave u_wave a_wave H_wave E_wave p_wave]=Roe_Mean_s(loL,loR,phiL,phiR,pL,pR,uL,uR)
global gama_s p0;
u_wave=(sqrt(loL*phiL)*uL+sqrt(loR*phiR)*uR)/(sqrt(loL*phiL)+sqrt(loR*phiR));
EL=(pL+gama_s*p0)/(gama_s-1)+0.5*loL*uL^2;
ER=(pR+gama_s*p0)/(gama_s-1)+0.5*loR*uR^2;
HL=(EL+pL)/loL;
HR=(ER+pR)/loR;
H_wave=(sqrt(loL*phiL)*HL+sqrt(loR*phiR)*HR)/(sqrt(loL*phiL)+sqrt(loR*phiR));
a_wave=sqrt((gama_s-1)*(H_wave-0.5*u_wave^2));
eta_wave=0.5*((pL+p0)/loL^gama_s+(pR+p0)/loR^gama_s);
lo_wave=(a_wave^2/eta_wave/gama_s)^(1.0/(gama_s-1.0));
%lo_wave=sqrt(loL*loR)
E_wave = lo_wave*(H_wave - a_wave^2/gama_s - p0/lo_wave);
p_wave = lo_wave*a_wave^2/gama_s-p0;
SL=u_wave-a_wave;
SR=u_wave+a_wave;
end
