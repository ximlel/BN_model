%comput Roe mean
function [SL SR u_wave a_wave H_wave]=Roe_Mean(loL,loR,pL,pR,vel_uL,vel_uR)
global gama;
u_wave=(sqrt(loL)*vel_uL+sqrt(loR)*vel_uR)/(sqrt(loL)+sqrt(loR));
EL=pL/(gama-1)+1/2*loL*vel_uL^2;
ER=pR/(gama-1)+1/2*loR*vel_uR^2;
HL=(EL+pL)/loL;
HR=(ER+pR)/loR;
H_wave=(sqrt(loL)*HL+sqrt(loR)*HR)/(sqrt(loL)+sqrt(loR));
a_wave=sqrt((gama-1)*(H_wave-1/2*u_wave^2));
SL=u_wave-a_wave;
SR=u_wave+a_wave;
end