%Riemann Solver with Roe scheme
function U_out=Riemann_solver_Roe_invariant(lo_gL,lo_gR,p_gL,p_gR,u_gL,u_gR,lo_sL,lo_sR,p_sL,p_sR,u_sL,u_sR,phi_sL,phi_sR)
%state constant
global gama_g gama_s;
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
eta_sL=p_sL/lo_sL^gama_s;
eta_gL=p_gL/lo_gL^gama_g;
QL=phi_gL*lo_gL*(u_gL-u_sL);
PL=phi_gL*lo_gL*(u_gL-u_sL)^2+phi_gL*p_gL+phi_sL*p_sL;
HL=0.5*(u_gL-u_sL)^2+gama_g/(gama_g-1.0)*p_gL/lo_gL;
eta_sR=p_sR/lo_sR^gama_s;
eta_gR=p_gR/lo_gR^gama_g;
QR=phi_gR*lo_gR*(u_gR-u_sR);
PR=phi_gR*lo_gR*(u_gR-u_sR)^2+phi_gR*p_gR+phi_sR*p_sR;
HR=0.5*(u_gR-u_sR)^2+gama_g/(gama_g-1.0)*p_gR/lo_gR;
eta_wave_g=0.5*(eta_gL+eta_gR);

%solve right eigenvetors
lo_ave_g = sqrt(phi_gR*lo_gR*phi_gL*lo_gL);
lo_ave_s = sqrt(phi_sR*lo_sR*phi_sL*lo_sL);
% p_ave_g = sqrt(phi_gR*phi_gL)*p_wave_g;
% lo_ave_g = lo_wave_g*(phi_gR+phi_gL)/2.0;
% lo_ave_s = lo_wave_s*(phi_sR+phi_sL)/2.0;
% p_ave_g = p_wave_g*(phi_gR+phi_gL)/2.0;
p_ave_g = lo_ave_g*a_wave_g^2/gama_g;

v_wave_rel = u_wave_g-u_wave_s;
K1=[ 0; 0; 0;                                        1;                        v_wave_rel-a_wave_g;                           -a_wave_g/lo_ave_g; 0];
K2=[ 0; 1; 0;-v_wave_rel*p_ave_g/eta_wave_g/a_wave_g^2;-v_wave_rel^2*p_ave_g/eta_wave_g/a_wave_g^2; 1.0/(gama_g-1.0)*p_ave_g/eta_wave_g/lo_ave_g; 0];
K3=[ 0; 0; 0;                                        1;                        v_wave_rel+a_wave_g;                            a_wave_g/lo_ave_g; 0];
K4=[-1; 0; 0;                                 lo_ave_g;    2*lo_ave_g*v_wave_rel+a_wave_s*lo_ave_s;                                   v_wave_rel; 0];
K5=[ 0; 0; 1;                                        0;                                          0;                                            0; 0];
K6=[ 0; 0; 0;                                        0;                                          0;                                            0; 1];
K7=[-1; 0; 0;                                 lo_ave_g;    2*lo_ave_g*v_wave_rel-a_wave_s*lo_ave_s;                                   v_wave_rel; 0];
K =[K1 K2 K3 K4 K5 K6 K7];	
%solve wave strengths
WL=[u_sL; eta_gL; eta_sL; QL; PL; HL; phi_sL];
WR=[u_sR; eta_gR; eta_sR; QR; PR; HR; phi_sR];
dW=WR-WL;

M = p_ave_g/(gama_g-1.0)/eta_wave_g*dW(2)-v_wave_rel*dW(4)+dW(5)-lo_ave_g*dW(6);
alpha1=-lo_ave_g*(v_wave_rel-a_wave_g)/2.0/a_wave_g*dW(1)+p_ave_g*(a_wave_g+(gama_g-1.0)*v_wave_rel)/2.0/(gama_g-1.0)/eta_wave_g/a_wave_g^2*dW(2)+0.5*dW(4)-lo_ave_g/2.0/a_wave_g*dW(6);              
alpha3= lo_ave_g*(v_wave_rel+a_wave_g)/2.0/a_wave_g*dW(1)-p_ave_g*(a_wave_g+(gama_g-1.0)*v_wave_rel)/2.0/(gama_g-1.0)/eta_wave_g/a_wave_g^2*dW(2)+0.5*dW(4)+lo_ave_g/2.0/a_wave_g*dW(6); 
alpha4=-0.5*dW(1)+M/2.0/lo_ave_s/a_wave_s;
alpha7=-0.5*dW(1)-M/2.0/lo_ave_s/a_wave_s;
alpha2=dW(2);
alpha5=dW(3);
alpha6=dW(7);
alpha=[alpha1;alpha2;alpha3;alpha4;alpha5;alpha6;alpha7];
%solve conservative vector at i+1/2
[lamda_sort,idx] = sort(lamda(1:7));
% lamda_pos = find(lamda_sort>0.0);

W11 = WL + K(:,idx(1:(min(find(idx==5),find(idx==7))-1)))*alpha(idx(1:(min(find(idx==5),find(idx==7))-1)));
W22 = WR - K(:,idx((max(find(idx==5),find(idx==7))+1):7))*alpha(idx((max(find(idx==5),find(idx==7))+1):7));
% solve W_out

U11 = W_comp_U(W11,lo_wave_g);
U22 = W_comp_U(W22,lo_wave_g);
U_out = 0.5*(U11 + U22);
end
