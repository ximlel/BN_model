%Riemann Solver with Roe scheme
function [out_flux_L,out_flux_R]=Riemann_solver_Roe(lo_gL,lo_gR,p_gL,p_gR,u_gL,u_gR,lo_sL,lo_sR,p_sL,p_sR,u_sL,u_sR,phi_sL,phi_sR,ratio_t_x)
%state constant
global gama_g gama_s p0;
phi_gL = 1.0-phi_sL;
phi_gR = 1.0-phi_sR;
%comput Roe mean
[S_gL,S_gR,lo_wave_g,u_wave_g,a_wave_g,H_wave_g,E_wave_g,p_wave_g]=Roe_Mean_g(lo_gL,lo_gR,phi_gL,phi_gR,p_gL,p_gR,u_gL,u_gR);
[S_sL,S_sR,lo_wave_s,u_wave_s,a_wave_s,H_wave_s,E_wave_s,p_wave_s]=Roe_Mean_s(lo_sL,lo_sR,phi_sL,phi_sR,p_sL,p_sR,u_sL,u_sR);
%solve averaged eigenvalues
lamda1=u_wave_g-a_wave_g;
lamda2=u_wave_g;
lamda3=u_wave_g+a_wave_g;
lamda4=u_wave_s-a_wave_s;
lamda5=u_wave_s;
lamda6=u_wave_s+a_wave_s;
lamda7=u_wave_s;
lamda=[lamda1 lamda2 lamda3 lamda4 lamda5 lamda6 lamda7];
%entropy fix
TOL=1e-6;
abs_lamda=zeros(7,1);
for rr=1:7
    if abs(lamda(rr))>=TOL
        abs_lamda(rr)=abs(lamda(rr));
    else
        abs_lamda(rr)=(lamda(rr)^2+TOL^2)/2/TOL;
    end
end
%solve right eigenvetors
K1=[1;      u_wave_g-a_wave_g;H_wave_g-u_wave_g*a_wave_g;0;0;0;0];
K2=[1;      u_wave_g         ;     0.5*u_wave_g^2       ;0;0;0;0];
K3=[1;      u_wave_g+a_wave_g;H_wave_g+u_wave_g*a_wave_g;0;0;0;0];
K4=[0;0;0;1;u_wave_s-a_wave_s;H_wave_s-u_wave_s*a_wave_s;0];
K5=[0;0;0;1;u_wave_s         ;     0.5*u_wave_s^2       ;0];
K6=[0;0;0;1;u_wave_s+a_wave_s;H_wave_s+u_wave_s*a_wave_s;0];
v_wave_rel = u_wave_g-u_wave_s;
K7=[gama_g*p_wave_g; gama_g*p_wave_g*u_wave_s;
lo_wave_g*v_wave_rel^2*H_wave_g-gama_g*p_wave_g*u_wave_g*v_wave_rel-E_wave_g*(v_wave_rel^2-a_wave_g^2);
(v_wave_rel^2-a_wave_g^2)*(gama_s*(p_wave_s+p0)+p_wave_g-p_wave_s)/a_wave_s^2;
(v_wave_rel^2-a_wave_g^2)*(gama_s*(p_wave_s+p0)+p_wave_g-p_wave_s)/a_wave_s^2*u_wave_s;
(v_wave_rel^2-a_wave_g^2)*( E_wave_s+(p_wave_g-p_wave_s)/a_wave_s^2*H_wave_s); v_wave_rel^2-a_wave_g^2];
K=[K1 K2 K3 K4 K5 K6 K7];	
%solve wave strengths
E_gL=p_gL/(gama_g-1)+0.5*lo_gL*u_gL^2;
E_gR=p_gR/(gama_g-1)+0.5*lo_gR*u_gR^2;
E_sL=(p_sL+gama_s*p0)/(gama_s-1)+0.5*lo_sL*u_sL^2;
E_sR=(p_sR+gama_s*p0)/(gama_s-1)+0.5*lo_sR*u_sR^2;
U_gL=[phi_gL*lo_gL;phi_gL*lo_gL*u_gL;phi_gL*E_gL];
U_gR=[phi_gR*lo_gR;phi_gR*lo_gR*u_gR;phi_gR*E_gR];
d_u_g=U_gR-U_gL;
U_sL=[phi_sL*lo_sL;phi_sL*lo_sL*u_sL;phi_sL*E_sL;phi_sL];
U_sR=[phi_sR*lo_sR;phi_sR*lo_sR*u_sR;phi_sR*E_sR;phi_sR];
d_u_s=U_sR-U_sL;
%alpha2=(gama-1)/a_wave_g^2*(d_u_g(1)*(H_wave_g-u_wave_g^2)+u_wave_g*d_u_g(2)-d_u_g(3));
%alpha1=1.0/2.0/a_wave_g*(d_u_g(1)*(u_wave_g+a_wave_g)-d_u_g(2)-a_wave_g*alpha2);
%alpha3=d_u_g(1)-alpha1-alpha2;
alpha1=( phi_gR*p_gR-phi_gL*p_gL-sqrt(phi_gR*lo_gR*phi_gL*lo_gL)*a_wave_g*(u_gR-u_gL))/(2.0*a_wave_g^2);
alpha2=(-phi_gR*p_gR+phi_gL*p_gL+    (phi_gR*lo_gR-phi_gL*lo_gL)*a_wave_g^2)          /(    a_wave_g^2);
alpha3=( phi_gR*p_gR-phi_gL*p_gL+sqrt(phi_gR*lo_gR*phi_gL*lo_gL)*a_wave_g*(u_gR-u_gL))/(2.0*a_wave_g^2);
alpha1=alpha1-d_u_s(4)*p_wave_g*(a_wave_g+(gama_g-1.0)*v_wave_rel)/(2.0*a_wave_g^2*(v_wave_rel-a_wave_g));
alpha2=alpha2+d_u_s(4)*p_wave_g*(gama_g-1.0)/a_wave_g^2;
alpha3=alpha3+d_u_s(4)*p_wave_g*(a_wave_g-(gama_g-1.0)*v_wave_rel)/(2.0*a_wave_g^2*(v_wave_rel+a_wave_g));
%alpha5=(gama-1)/a_wave_s^2*(d_u_s(1)*(H_wave_s-u_wave_s^2)+u_wave_s*d_u_s(2)-d_u_s(3));
%alpha4=1.0/2.0/a_wave_s*(d_u_s(1)*(u_wave_s+a_wave_s)-d_u_s(2)-a_wave_s*alpha5);
%alpha6=d_u_s(1)-alpha4-alpha5;
alpha4=( phi_sR*p_sR-phi_sL*p_sL-sqrt(phi_sR*lo_sR*phi_sL*lo_sL)*a_wave_s*(u_sR-u_sL))/(2.0*a_wave_s^2);
alpha5=(-phi_sR*p_sR+phi_sL*p_sL+    (phi_sR*lo_sR-phi_sL*lo_sL)*a_wave_s^2)          /(    a_wave_s^2);
alpha6=( phi_sR*p_sR-phi_sL*p_sL+sqrt(phi_sR*lo_sR*phi_sL*lo_sL)*a_wave_s*(u_sR-u_sL))/(2.0*a_wave_s^2);
alpha4=alpha4-d_u_s(4)*p_wave_g/(2.0*a_wave_s^2);
alpha5=alpha5-d_u_s(4)*(p_wave_s*(gama_s-1.0)+gama_s*p0)/a_wave_s^2;
alpha6=alpha6-d_u_s(4)*p_wave_g/(2.0*a_wave_s^2);
alpha7=-d_u_s(4)/(a_wave_g^2-v_wave_rel^2);
alpha=[alpha1;alpha2;alpha3;alpha4;alpha5;alpha6;alpha7];
%solve conservative vector at i+1/2
[lamda_sort,idx] = sort(lamda(1:7));
lamda_pos = find(lamda_sort>0.0);
UL=[phi_gL*lo_gL;phi_gL*lo_gL*u_gL;phi_gL*E_gL;phi_sL*lo_sL;phi_sL*lo_sL*u_sL;phi_sL*E_sL;phi_sL];
UR=[phi_gR*lo_gR;phi_gR*lo_gR*u_gR;phi_gR*E_gR;phi_sR*lo_sR;phi_sR*lo_sR*u_sR;phi_sR*E_sR;phi_sR];
if lamda_sort(1)>=0.0
    U0 = UL;
elseif lamda_sort(7)<=0.0
    U0 = UR;
elseif u_wave_s > 0.0
    U0 = UL + K(:,idx(1:(lamda_pos(1)-1)))*alpha(idx(1:(lamda_pos(1)-1)));
else
    U0 = UR - K(:,idx(lamda_pos(1):7))    *alpha(idx(lamda_pos(1):7));
end
[lo_g0,u_g0,p_g0,phi_g0,lo_s0,u_s0,p_s0,phi_s0]=primitive_comp(U0);
U11 = UL + K(:,idx(1:(min(find(idx==5),find(idx==7))-1)))*alpha(idx(1:(min(find(idx==5),find(idx==7))-1)));
[lo_g1,u_g1,p_g1,phi_g1,lo_s1,u_s1,p_s1,phi_s1]=primitive_comp(U11);
U22 = UR - K(:,idx((max(find(idx==5),find(idx==7))+1):7))*alpha(idx((max(find(idx==5),find(idx==7))+1):7));
[lo_g2,u_g2,p_g2,phi_g2,lo_s2,u_s2,p_s2,phi_s2]=primitive_comp(U22);
%solve flux at i+1/2
F0=[phi_g0*lo_g0*u_g0;phi_g0*lo_g0*u_g0^2+phi_g0*p_g0;phi_g0*(gama_g/(gama_g-1.0)*p_g0+0.5*lo_g0*u_g0^2)*u_g0;
    phi_s0*lo_s0*u_s0;phi_s0*lo_s0*u_s0^2+phi_s0*p_s0;phi_s0*(gama_s/(gama_s-1.0)*(p_s0+p0)+0.5*lo_s0*u_s0^2)*u_s0;0.0];
out_flux_L=F0;
out_flux_R=F0;
rat = 1;
if u_wave_s > 0.0
    out_flux_R=out_flux_R-[0;0;0;0;0;0;u_wave_s]*d_u_s(4);
     out_flux_R=out_flux_R-[0;p_g0-p_g1;u_s0*p_g0-u_wave_s*p_g1;0;p_g1-p_g0;u_wave_s*p_g1-u_s0*p_g0;u_s0-u_wave_s]*d_u_s(4);
    out_flux_R=out_flux_R-rat*[0;1;u_wave_s;0;-1;-u_wave_s;0]*(phi_s2*p_s2-phi_s1*p_s1);
    out_flux_R=out_flux_R-(1-rat)*[0;p_gR;u_sR*p_gR;0;-p_gR;-u_sR*p_gR;u_sR-u_wave_s]*d_u_s(4);
    %out_flux_R=out_flux_R-(1-rat)*[0;1;u_wave_s;0;-1;-u_wave_s;0]*(phi_g1*p_g1+phi_g1*lo_g1*(u_g1-u_wave_s)^2-phi_g2*p_g2-phi_g2*lo_g2*(u_g2-u_wave_s)^2);
    %out_flux_R=out_flux_R-(1-rat)*[0;p_g2;u_wave_s*p_g2;0;-p_g2;-u_wave_s*p_g2;0]*d_u_s(4);
    %out_flux_R=out_flux_R+1.0/(gama-1.0)*u_wave_s*[0;0;1;0;0;-1;0]*(p_s2-p_s1)*d_u_s(4);
else
    out_flux_L=out_flux_L+[0;0;0;0;0;0;u_wave_s]*d_u_s(4);
     out_flux_L=out_flux_L+[0;p_g0-p_g2;u_s0*p_g0-u_wave_s*p_g2;0;p_g2-p_g0;u_wave_s*p_g2-u_s0*p_g0;u_s0-u_wave_s]*d_u_s(4);
    out_flux_L=out_flux_L+rat*[0;1;u_wave_s;0;-1;-u_wave_s;0]*(phi_s2*p_s2-phi_s1*p_s1);
    out_flux_L=out_flux_L+(1-rat)*[0;p_gL;u_sL*p_gL;0;-p_gL;-u_sL*p_gL;u_sL-u_wave_s]*d_u_s(4);
    %out_flux_L=out_flux_L+(1-rat)*[0;1;u_wave_s;0;-1;-u_wave_s;0]*(phi_g1*p_g1+phi_g1*lo_g1*(u_g1-u_wave_s)^2-phi_g2*p_g2-phi_g2*lo_g2*(u_g2-u_wave_s)^2);
    %out_flux_L=out_flux_L+(1-rat)*[0;p_g1;u_wave_s*p_g1;0;-p_g1;-u_wave_s*p_g1;0]*d_u_s(4);
    %out_flux_L=out_flux_L-1.0/(gama-1.0)*u_wave_s*[0;0;1;0;0;-1;0]*(p_s2-p_s1)*d_u_s(4);
end
end
