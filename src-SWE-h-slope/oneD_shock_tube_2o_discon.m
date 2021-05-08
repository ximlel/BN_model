clear;
clc;
%1D shock_tube by GRP Schemes for shallow water model
%state constant
global g;
g=9.81;
global ep;
ep=1e-6;
x0=0;

%state value
Time=0;
Tend=0.04;
CFL=0.5;
Alpha_GRP=1.0;

%initial condition
% h_L_0  =4.7;
% u_L_0  =0.651;
% h_R_0  =0.8;
% u_R_0  =2.048;
% Z_L_0  =0;
% Z_R_0  =1;
% x_min=-1;
% x_max=2;
% N_0=300;

% h_L_0  =2;
% u_L_0  =0;
% h_R_0  =1;
% u_R_0  =0;
% Z_L_0  =0;
% Z_R_0  =1;
% x_min=-1;
% x_max=2;
% N_0=300;

% h_L_0  =4;
% u_L_0  =-2.894;
% h_R_0  =2.7;
% u_R_0  =-3;
% Z_L_0  =0;
% Z_R_0  =1;
% x_min=-3;
% x_max=1;
% N_0=400;

% x_min=0;
% x_max=25;
% N_0=100;

d_x=(x_max-x_min)/N_0;
N=N_0+1;

% Z_0   = 2;
% q_in  = 4.42;
% h_out = 2;

Z=@(x) Z_L_0*(x <= 0.0) + Z_R_0*(x > 0.0);
%Z=@(x) (0.2-0.05*(x-10)^2)*(x>8 && x<12);

Z_L=zeros(1,N);
Z_R=zeros(1,N);
Z_M=zeros(1,N+1);

U=zeros(2,N);
F=zeros(2,N+1);
h_L_int =zeros(1,N+1);
h_R_int =zeros(1,N+1);
u_L_int =zeros(1,N+1);
u_R_int =zeros(1,N+1);
dh_L_int=zeros(1,N+1);
dh_R_int=zeros(1,N+1);
du_L_int=zeros(1,N+1);
du_R_int=zeros(1,N+1);
h_mid=zeros(1,N+1);
W_int=zeros(2,N+1);
dZ=zeros(1,N-1);
dh=zeros(1,N+1);
%du=zeros(1,N+1);

Fr_L=zeros(1,N);
Fr_R=zeros(1,N);
h_L=zeros(1,N);
h_R=zeros(1,N);
u_L=zeros(1,N);
u_R=zeros(1,N);
h_mL=zeros(1,N);
h_mR=zeros(1,N);
u_mL=zeros(1,N);
u_mR=zeros(1,N);
a_L=zeros(1,N);
a_R=zeros(1,N);
H_t=zeros(1,N);
dq =zeros(1,N);

x=zeros(1,N);
x_M=zeros(1,N+1);
%test begin
for i=1:N+1
    x_M(i) = x_min+(i-0.5)*d_x;
    Z_M(i) = Z(x_M(i));
end
for i=1:N
    x(i)=x_min+(i-1)*d_x;
    Z_L(i) = Z(x(i)-ep);
    Z_R(i) = Z(x(i)+ep);
end
for i=1:N-1
    dZ(i)  = (Z_L(i+1)-Z_R(i))/d_x;
end

U_L_0=[h_L_0;h_L_0*u_L_0];
U_R_0=[h_R_0;h_R_0*u_R_0];
for i=1:N
    if x(i) < -ep
        U(:,i)=U_L_0;
    elseif x(i) > ep
        U(:,i)=U_R_0;
    else
        U(:,i)=0.5*(U_L_0+U_R_0);
    end
end

% for i=1:N
%     U(:,i)=[Z_0-Z(x(i));0];
% end

%Godunov's Method
while Time<Tend && isreal(Time)
    %boundary condition
    %U(2,1)=q_in;
    %U(1,1)=Z_0;
    %U(1,N)=h_out;
    %CFL condition
    for i=1:N
        [h_L(i),u_L(i),h_R(i),u_R(i),H_t(i)]=primitive_comp(U(:,i),Z_L(i),Z_R(i));
        a_L(i)=sqrt(g*h_L(i));
        a_R(i)=sqrt(g*h_R(i));
        if u_L(i) > a_L(i)
            Fr_L(i)=2;
        else
            Fr_L(i)=0;
        end
        if u_R(i) > a_R(i)
            Fr_R(i)=2;
        else
            Fr_L(i)=0;
        end
    end
    hh = U(1,:);
    qq = U(2,:);
    Smax=max(max(abs(u_L)+a_L),max(abs(u_R)+a_R));
    d_t=CFL*d_x/Smax;
    if Time+d_t >= Tend
        d_t = Tend-Time+1e-10;
    end
    %reconstruction (minmod limiter)
    for i=2:N-1
%         if Time < ep
            dh(i) =minmod(Alpha_GRP*(h_R(i)-h_R(i-1))/d_x,(h_L(i+1)-h_R(i-1))/2.0/d_x,Alpha_GRP*(h_L(i+1)-h_L(i))/d_x);
            dh(i) =minmod(Alpha_GRP*(h_L(i)-h_R(i-1))/d_x,dh(i),                      Alpha_GRP*(h_L(i+1)-h_R(i))/d_x);
            dq(i) =minmod(Alpha_GRP*(qq(i) -qq(i-1)) /d_x,(qq(i+1) -qq(i-1)) /2.0/d_x,Alpha_GRP*(qq(i+1) -qq(i)) /d_x);
%         else
%             dh(i) =minmod(Alpha_GRP*(h_R(i)-h_R(i-1))/d_x,(W_int(1,i+1)-W_int(1,i))/1.0/d_x,Alpha_GRP*(h_L(i+1)-h_L(i))/d_x);
%             dh(i) =minmod(Alpha_GRP*(h_L(i)-h_R(i-1))/d_x,dh(i),                            Alpha_GRP*(h_L(i+1)-h_R(i))/d_x);
%             dq(i) =minmod(Alpha_GRP*(qq(i) -qq(i-1)) /d_x,(W_int(2,i+1)-W_int(2,i))/1.0/d_x,Alpha_GRP*(qq(i+1) -qq(i)) /d_x);
%         end
    end
    for i=1:N+1
        if i==1
            h_L_int(i)=h_R(1);
            u_L_int(i)=u_R(1);
            h_R_int(i)=h_R(1);
            u_R_int(i)=u_R(1);           
            dh_L_int(i)=0;
            du_L_int(i)=0;
        elseif i==N+1 
            h_L_int(i)=h_L(N);
            u_L_int(i)=u_L(N);
            h_R_int(i)=h_L(N);
            u_R_int(i)=u_L(N);
            dh_R_int(i)=0;
            du_R_int(i)=0;
        else
            [h_L_int(i),u_L_int(i),dh_L_int(i),du_L_int(i)]=dRI2dU_cal(qq(i-1)+0.5*d_x*dq(i-1),h_R(i-1)+0.5*d_x*dh(i-1),dq(i-1),dh(i-1),dZ(i-1),Fr_R(i-1));
            [h_R_int(i),u_R_int(i),dh_R_int(i),du_R_int(i)]=dRI2dU_cal(qq(i)  -0.5*d_x*dq(i),  h_L(i)  -0.5*d_x*dh(i)  ,dq(i),  dh(i),  dZ(i-1),Fr_L(i));
        end
    end    
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
        if i==1
            dZZ=0.0;
        elseif i==N+1
            dZZ=0.0;
        else
            dZZ=dZ(i-1);            
        end
        [h_mid(:,i),F(:,i),W_int(:,i)]=GRP_solver(h_L_int(i),h_R_int(i),dh_L_int(i),dh_R_int(i),u_L_int(i),u_R_int(i),du_L_int(i),du_R_int(i),Z_M(i),dZZ,dZZ,d_t);
    end    
    for i=1:N
        if i==1 || i==N
            h_mL(i) = h_L(i);
            h_mR(i) = h_R(i);            
            u_mL(i) = u_L(i);
            u_mR(i) = u_R(i); 
        else
            [h_mL(i),u_mL(i),dh_mL,du_mL]=dRI2dU_cal(qq(i),h_L(i),dq(i-1),dh(i-1),dZ(i-1),Fr_L(i));
            [h_mR(i),u_mR(i),dh_mR,du_mR]=dRI2dU_cal(qq(i),h_R(i),dq(i-1),dh(i-1),dZ(i),  Fr_R(i));
            h_mL(i) = h_mL(i) - 0.5*d_t*(h_mL(i)*du_mL+u_mL(i)*dh_mL);
            h_mR(i) = h_mR(i) - 0.5*d_t*(h_mR(i)*du_mR+u_mR(i)*dh_mR);        
            u_mL(i) = u_mL(i) - 0.5*d_t*(dh_mL+u_mL(i)*du_mL/g+dZ(i-1));
            u_mR(i) = u_mR(i) - 0.5*d_t*(dh_mR+u_mR(i)*du_mR/g+dZ(i));
        end
    end
    %compute U in next step
    for i=1:N
        if abs(Z_R(i)-Z_L(i))<ep
            S=-g*0.5*(h_mL(i)+h_mR(i))*(Z_R(i)-Z_L(i));
        else
            S_tmp=(h_mR(i)*u_mR(i)^2+g*h_mR(i)^2/2-h_mL(i)*u_mL(i)^2-g*h_mL(i)^2/2);
%             if (S_tmp/g/(Z_L(i)-Z_R(i))>max(h_mL(i),h_mR(i)))
%                 S=-g*max(h_mL(i),h_mR(i))*(Z_R(i)-Z_L(i));        
%             elseif (S_tmp/g/(Z_L(i)-Z_R(i))<min(h_mL(i),h_mR(i)))
%                 S=-g*min(h_mL(i),h_mR(i))*(Z_R(i)-Z_L(i));
%             else
                S=S_tmp;
%             end
        end
        S = S - 0.5*g*(h_mid(i)+h_mL(i))*(Z_L(i)-Z_M(i)) - 0.5*g*(h_mid(i+1)+h_mR(i))*(Z_M(i+1)-Z_R(i));
        U(:,i)=U(:,i)+d_t/d_x*(F(:,i)-F(:,i+1))+d_t/d_x*[0;S];
    end
    Time = Time+d_t
% if Time > 0.002
%    break;
% end
end
% h = 0.5*(h_L+h_R);
% u = 0.5*(u_L+u_R);
h=h_L;
u=u_L;

N_MAX = 3000;
d_xM=(x_max-x_min)/N_MAX;
W_exact = zeros(N_MAX,2);
% load ./EXACT_SWE1.mat;
% W_exact(:,1)=h_E';;
% W_exact(:,2)=u_E';

%plot
col = '-m';
%col = '+k';
% h=figure(1);
% set(h,'position',[100 100 1150 750]);
subplot(1,2,1);
hold on
%plot(x_min:d_xM:x_max-d_xM,W_exact(:,1),'b','LineWidth',0.4);
plot(x,h+Z_L(1,1:N),col,'MarkerSize',6);%col,'LineWidth',1.0);
%xlabel('Position','FontWeight','bold');
%ylabel('Density','FontWeight','bold');
%ylim([90 160])
%ylim([0 2])
title('h+z');
set(gca,'box','on');
subplot(1,2,2);
hold on
% plot(x_min:d_xM:x_max-d_xM,W_exact(:,2),'b','LineWidth',0.4);
plot(x,u,col,'MarkerSize',6);%col,'LineWidth',1.0);
%xlabel('Position','FontWeight','bold');
%ylabel('Velocity','FontWeight','bold');
%ylim([0 2])
title('u');
set(gca,'box','on');
