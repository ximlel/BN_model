%GRP Solver with HLLC scheme (subsonic case)
function [out_flux_L,out_flux_R,out_U_int]=GRP_solver_HLLC(lo_gL,lo_gR,p_gL,p_gR,u_gL,u_gR,lo_sL,lo_sR,p_sL,p_sR,u_sL,u_sR,phi_sL,phi_sR,d_UL,d_UR,UL_mid,UR_mid,ratio_t_x,d_t)
%state constant
global gama_s gama_g p0;
phi_gL = 1.0-phi_sL;
phi_gR = 1.0-phi_sR;
%comput wave speed
[S_gL,S_gM,S_gR,S_sL,S_sM,S_sR,lo_g_srL,lo_g_srR,p_g1,p_g2,u_g1,u_g2,lo_s1,lo_s2,p_s1,p_s2,u_s1,u_s2,phi_s1,phi_s2] = solid_cont(lo_gL,lo_gR,p_gL,p_gR,u_gL,u_gR,lo_sL,lo_sR,p_sL,p_sR,u_sL,u_sR,phi_sL,phi_sR);
phi_g1=1.0-phi_s1;
phi_g2=1.0-phi_s2;
%solve conservative vector at i+1/2
E_gL=p_gL/(gama_g-1)+0.5*lo_gL*u_gL^2;
E_gR=p_gR/(gama_g-1)+0.5*lo_gR*u_gR^2;
E_sL=(p_sL+gama_s*p0)/(gama_s-1)+0.5*lo_sL*u_sL^2;
E_sR=(p_sR+gama_s*p0)/(gama_s-1)+0.5*lo_sR*u_sR^2;
UL=[phi_gL*lo_gL;phi_gL*lo_gL*u_gL;phi_gL*E_gL;phi_sL*lo_sL;phi_sL*lo_sL*u_sL;phi_sL*E_sL;phi_sL];
UR=[phi_gR*lo_gR;phi_gR*lo_gR*u_gR;phi_gR*E_gR;phi_sR*lo_sR;phi_sR*lo_sR*u_sR;phi_sR*E_sR;phi_sR];
FL=[phi_gL*lo_gL*u_gL;phi_gL*lo_gL*u_gL^2+phi_gL*p_gL;phi_gL*(E_gL+p_gL)*u_gL;
    phi_sL*lo_sL*u_sL;phi_sL*lo_sL*u_sL^2+phi_sL*p_sL;phi_sL*(E_sL+p_sL)*u_sL;0.0];
FR=[phi_gR*lo_gR*u_gR;phi_gR*lo_gR*u_gR^2+phi_gR*p_gR;phi_gR*(E_gR+p_gR)*u_gR;
    phi_sR*lo_sR*u_sR;phi_sR*lo_sR*u_sR^2+phi_sR*p_sR;phi_sR*(E_sR+p_sR)*u_sR;0.0];
U11=zeros(7,1);
U22=zeros(7,1);
U11(7) = phi_sL;
U22(7) = phi_sR;
U11(4:6) = phi_sL*lo_sL*(S_sL-u_sL)/(S_sL-S_sM)*[1;S_sM;E_sL/lo_sL+(S_sM-u_sL)*(S_sM+p_sL/lo_sL/(S_sL-u_sL))];
U22(4:6) = phi_sR*lo_sR*(S_sR-u_sR)/(S_sR-S_sM)*[1;S_sM;E_sR/lo_sR+(S_sM-u_sR)*(S_sM+p_sR/lo_sR/(S_sR-u_sR))];
U_g_srL  = phi_gL*lo_gL*(S_gL-u_gL)/(S_gL-u_g1)*[1;u_g1;E_gL/lo_gL+(u_g1-u_gL)*(u_g1+p_gL/lo_gL/(S_gL-u_gL))];
U_g_srR  = phi_gR*lo_gR*(S_gR-u_gR)/(S_gR-u_g2)*[1;u_g2;E_gR/lo_gR+(u_g2-u_gR)*(u_g2+p_gR/lo_gR/(S_gR-u_gR))];
if S_gM < S_sM
    lo_g1=lo_g_srR*(p_g1/p_g2)^(1/gama_g);
    E_g1 =p_g1/(gama_g-1)+0.5*lo_g1*u_g1^2;
    lo_g2=lo_g_srR;
    U11(1:3) = phi_g1*[lo_g1;lo_g1*u_g1;E_g1];
    U22(1:3) = U_g_srR;
else
    lo_g1=lo_g_srL;
    lo_g2=lo_g_srL*(p_g2/p_g1)^(1/gama_g);
    E_g2 =p_g2/(gama_g-1)+0.5*lo_g2*u_g2^2;
    U11(1:3) = U_g_srL;
    U22(1:3) = phi_g2*[lo_g2;lo_g2*u_g2;E_g2];
end
%solver GRP at star region
d_U_gL=d_UL(1:3);
d_U_gR=d_UR(1:3);
d_U_sL=d_UL(4:6);
d_U_sR=d_UR(4:6);
d_phi_sL=d_UL(7);
d_phi_sR=d_UR(7);
CA_g_Lsr=(S_gL*eye(4)-quasilinear(lo_g_srL,u_g1,p_g1,p_g1,u_s1,'g','C'))*(S_gL*eye(4)-quasilinear(lo_g_srL,u_g1,p_g1,p_g1,u_s1,'g','A'));
CA_s_Lsr=(S_sL*eye(4)-quasilinear(lo_s1   ,u_s1,p_s1,p_g1,u_s1,'s','C'))*(S_sL*eye(4)-quasilinear(lo_s1   ,u_s1,p_s1,p_g1,u_s1,'s','A'));
CA_g_LL =(S_gL*eye(4)-quasilinear(lo_gL   ,u_gL,p_gL,p_gL,u_sL,'g','C'))*(S_gL*eye(4)-quasilinear(lo_gL   ,u_gL,p_gL,p_gL,u_sL,'g','A'));
CA_s_LL =(S_sL*eye(4)-quasilinear(lo_sL   ,u_sL,p_sL,p_gL,u_sL,'s','C'))*(S_sL*eye(4)-quasilinear(lo_sL   ,u_sL,p_sL,p_gL,u_sL,'s','A'));
CA_g_Rsr=(S_gR*eye(4)-quasilinear(lo_g_srR,u_g2,p_g2,p_g2,u_s2,'g','C'))*(S_gR*eye(4)-quasilinear(lo_g_srR,u_g2,p_g2,p_g2,u_s2,'g','A'));
CA_s_Rsr=(S_sR*eye(4)-quasilinear(lo_s2   ,u_s2,p_s2,p_s2,u_s2,'s','C'))*(S_sR*eye(4)-quasilinear(lo_s2   ,u_s2,p_s2,p_s2,u_s2,'s','A'));
CA_g_RR =(S_gR*eye(4)-quasilinear(lo_gR   ,u_gR,p_gR,p_gR,u_sR,'g','C'))*(S_gR*eye(4)-quasilinear(lo_gR   ,u_gR,p_gR,p_gR,u_sR,'g','A'));
CA_s_RR =(S_sR*eye(4)-quasilinear(lo_sR   ,u_sR,p_sR,p_sR,u_sR,'s','C'))*(S_sR*eye(4)-quasilinear(lo_sR   ,u_sR,p_sR,p_sR,u_sR,'s','A'));
d_U_gsr=least_square(d_U_gL,d_U_gR,d_phi_sL,d_phi_sR,CA_g_LL,CA_g_Lsr,CA_g_RR,CA_g_Rsr);
d_U_ssr=least_square(d_U_sL,d_U_sR,d_phi_sL,d_phi_sR,CA_s_LL,CA_s_Lsr,CA_s_RR,CA_s_Rsr);
%solve flux at i+1/2
U0=zeros(7,1);
F0=zeros(7,1);
d_U0=zeros(7,1);
if S_sL>=0.0
    U0(4:7) = UL(4:7);
    F0(4:6) = FL(4:6);
    d_U0(4:7)=d_UL(4:7);
elseif S_sR<=0.0
    U0(4:7) = UR(4:7);
    F0(4:6) = FR(4:6);
    d_U0(4:7)=d_UR(4:7);
elseif S_sM>0.0
    U0(4:7) = U11(4:7);
    F0(4:6) = FL(4:6)+S_sL*(U0(4:6)-UL(4:6));
    d_U0(4:6)=d_U_ssr;
    d_U0(7) = d_UL(7);
else
    U0(4:7) = U22(4:7);
    F0(4:6) = FR(4:6)+S_sR*(U0(4:6)-UR(4:6));
    d_U0(4:6)=d_U_ssr;
    d_U0(7) = d_UR(7);    
end
if S_gL>=0.0
    U0(1:3) = UL(1:3);
    F0(1:3) = FL(1:3);
    d_U0(1:3)=d_UL(1:3);
elseif S_gR<=0.0
    U0(1:3) = UR(1:3);
    F0(1:3) = FR(1:3);
    d_U0(1:3)=d_UR(1:3);
elseif S_sM<0.0 && S_gM<0.0
    U0(1:3) = U_g_srR;
    F0(1:3) = FR(1:3)+S_gR*(U0(1:3)-UR(1:3));
    d_U0(1:3)=d_U_gsr;
elseif S_sM<0.0
    U0(1:3) = U22(1:3);
    F0(1:3) = FR(1:3)+S_gR*(U_g_srR-UR(1:3))+S_gM*(U0(1:3)-U_g_srR);
    d_U0(1:3)=d_U_gsr;
elseif S_gM<0.0
    U0(1:3) = U11(1:3);
    F0(1:3) = FL(1:3)+S_gL*(U_g_srL-UL(1:3))+S_gM*(U0(1:3)-U_g_srL);
    d_U0(1:3)=d_U_gsr;
else
    U0(1:3) = U_g_srL;
    F0(1:3) = FL(1:3)+S_gL*(U0(1:3)-UL(1:3));
    d_U0(1:3)=d_U_gsr;
end
[lo_g0,u_g0,p_g0,phi_g0,lo_s0,u_s0,p_s0,phi_s0]=primitive_comp(U0);
C_U=zeros(7,7);
C_U(1:4,1:4)=quasilinear(lo_g0,u_g0,p_g0,p_g0,u_s0,'g','C');
C_U([1 5:7],[1 5:7])=quasilinear(lo_s0,u_s0,p_s0,p_g0,u_s0,'s','C');
A_U=zeros(7,7);
A_U(1:4,1:4)=quasilinear(lo_g0,u_g0,p_g0,p_g0,u_s0,'g','A');
A_U([1 5:7],[1 5:7])=quasilinear(lo_s0,u_s0,p_s0,p_g0,u_s0,'s','A');
out_U_int=U0-d_t*C_U*d_U0;
% phi_g0 = 1.0 - phi_s0;
% F0=[phi_g0*lo_g0*u_g0;phi_g0*lo_g0*u_g0^2+phi_g0*p_g0;(U0(3)+phi_g0*p_g0)*u_g0;
%     phi_s0*lo_s0*u_s0;phi_s0*lo_s0*u_s0^2+phi_s0*p_s0;(U0(6)+phi_s0*p_s0)*u_s0;0.0];
out_flux_L=F0-0.5*d_t*A_U*C_U*d_U0;
out_flux_R=F0-0.5*d_t*A_U*C_U*d_U0;
% [lo_g1 u_g1 p_g1 phi_g1 lo_s1 u_s1 p_s1 phi_s1]=primitive_comp(UL);
% [lo_g2 u_g2 p_g2 phi_g2 lo_s2 u_s2 p_s2 phi_s2]=primitive_comp(UR);
% out_flux_L=FL;
% out_flux_R=FL;
% S_sM = u_sL;
d_u_s=phi_sR-phi_sL;
out_U_mid=U0-0.5*d_t*C_U*d_U0;
[lo_g0_mid,u_g0_mid,p_g0_mid,phi_g0_mid,lo_s0_mid,u_s0_mid,p_s0_mid,phi_s0_mid]=primitive_comp(out_U_mid);
F0_mid=[phi_g0_mid*lo_g0_mid*u_g0_mid;phi_g0_mid*lo_g0_mid*u_g0_mid^2+phi_g0_mid*p_g0_mid;(out_U_mid(3)+phi_g0_mid*p_g0_mid)*u_g0_mid;
        phi_s0_mid*lo_s0_mid*u_s0_mid;phi_s0_mid*lo_s0_mid*u_s0_mid^2+phi_s0_mid*p_s0_mid;(out_U_mid(6)+phi_s0_mid*p_s0_mid)*u_s0_mid;0.0];
out_flux_L=F0_mid;
out_flux_R=F0_mid;
[lo_gL_mid,u_gL_mid,p_gL_mid,phi_gL_mid,lo_sL_mid,u_sL_mid,p_sL_mid,phi_sL_mid]=primitive_comp(UL_mid);
[lo_gR_mid,u_gR_mid,p_gR_mid,phi_gR_mid,lo_sR_mid,u_sR_mid,p_sR_mid,phi_sR_mid]=primitive_comp(UR_mid);
out_flux_R=out_flux_R-[0;p_g0_mid;p_g0_mid*u_s0_mid;0;-p_g0_mid;-p_g0_mid*u_s0_mid;u_s0_mid]*(phi_sR_mid-phi_s0_mid);
out_flux_L=out_flux_L+[0;p_g0_mid;p_g0_mid*u_s0_mid;0;-p_g0_mid;-p_g0_mid*u_s0_mid;u_s0_mid]*(phi_s0_mid-phi_sL_mid);

%non-conservative term
rat = 1;
if S_sM > 0.0
    out_flux_R=out_flux_R-[0;0;0;0;0;0;S_sM]*d_u_s;
    %out_flux_R=out_flux_R-[0;p_g0-p_g1;u_s0*p_g0-S_sM*p_g1;0;p_g1-p_g0;S_sM*p_g1-u_s0*p_g0;u_s0-S_sM]*d_u_s;
    out_flux_R=out_flux_R-[0;-p_g1;-S_sM*p_g1;0;p_g1;S_sM*p_g1;-S_sM]*d_u_s;
    out_flux_R=out_flux_R-rat*[0;1;S_sM;0;-1;-S_sM;0]*(phi_s2*p_s2-phi_s1*p_s1);
    out_flux_R=out_flux_R-(1-rat)*[0;p_gR;u_sR*p_gR;0;-p_gR;-u_sR*p_gR;u_sR-S_sM]*d_u_s;
    %dE_s = E_s_correct(S_sM,lo_g1,u_g1,p_g1,lo_s1,p_s1,phi_s1,lo_g2,u_g2,p_g2,lo_s2,p_s2,phi_s2,ratio_t_x);
    %out_flux_R=out_flux_R-[0;0;-dE_s;0;0;dE_s;0]/ratio_t_x;
    %out_flux_R=out_flux_R-(1-rat)*[0;1;S_sM;0;-1;-S_sM;0]*(phi_g1*p_g1+phi_g1*lo_g1*(u_g1-S_sM)^2-phi_g2*p_g2-phi_g2*lo_g2*(u_g2-S_sM)^2);
    %out_flux_R=out_flux_R-(1-rat)*[0;p_g2;S_sM*p_g2;0;-p_g2;-S_sM*p_g2;0]*d_u_s;
    %out_flux_R=out_flux_R+1.0/(gama-1.0)*S_sM*[0;0;1;0;0;-1;0]*(p_s2-p_s1)*d_u_s;
else
    out_flux_L=out_flux_L+[0;0;0;0;0;0;S_sM]*d_u_s;
    %out_flux_L=out_flux_L+[0;p_g0-p_g2;u_s0*p_g0-S_sM*p_g2;0;p_g2-p_g0;S_sM*p_g2-u_s0*p_g0;u_s0-S_sM]*d_u_s;
    out_flux_L=out_flux_L+[0;-p_g2;-S_sM*p_g2;0;p_g2;S_sM*p_g2;-S_sM]*d_u_s;
    out_flux_L=out_flux_L+rat*[0;1;S_sM;0;-1;-S_sM;0]*(phi_s2*p_s2-phi_s1*p_s1);
    out_flux_L=out_flux_L+(1-rat)*[0;p_gL;u_sL*p_gL;0;-p_gL;-u_sL*p_gL;u_sL-S_sM]*d_u_s;
    %dE_s = E_s_correct(-S_sM,lo_g2,-u_g2,p_g2,lo_s2,p_s2,phi_s2,lo_g1,-u_g1,p_g1,lo_s1,p_s1,phi_s1,ratio_t_x);
    %out_flux_L=out_flux_L+[0;0;-dE_s;0;0;dE_s;0]/ratio_t_x;
    %out_flux_L=out_flux_L+(1-rat)*[0;1;S_sM;0;-1;-S_sM;0]*(phi_g1*p_g1+phi_g1*lo_g1*(u_g1-S_sM)^2-phi_g2*p_g2-phi_g2*lo_g2*(u_g2-S_sM)^2);
    %out_flux_L=out_flux_L+(1-rat)*[0;p_g1;S_sM*p_g1;0;-p_g1;-S_sM*p_g1;0]*d_u_s;
    %out_flux_L=out_flux_L-1.0/(gama-1.0)*S_sM*[0;0;1;0;0;-1;0]*(p_s2-p_s1)*d_u_s;
end
end
