clear;
clc;
%1D shock_tube by GRP Schemes for shallow water model
%state constant
global g;
g=9.81;
global ep;
ep=1e-6;
x0=0;
x_min=-1;
x_max=2;
N=300;
d_x=(x_max-x_min)/N;
%state value
Time=0;
Tend=0.5;
Alpha_GRP=0;
CFL=0.5;

%initial condition
h_L_0  =4.7;
u_L_0  =0.651;
h_R_0  =0.8;
u_R_0  =2.048;
Z_L_0  =0;
Z_R_0  =1;

% h_L_0  =2;
% u_L_0  =0;
% h_R_0  =1;
% u_R_0  =0;
% Z_L_0  =0;        
% Z_R_0  =1;


% h_L_0  =4;
% u_L_0  =-2.894;
% h_R_0  =2.7;
% u_R_0  =-3;
% Z_L_0  =0;
% Z_R_0  =1;
% x_min=-3;
% x_max=1;
% N=400;
% d_x=(x_max-x_min)/N;


Z=zeros(1,N+1);
U=zeros(2,N);
F=zeros(2,N+1);
% load ../test/test_new1_pi.mat;
%test begin
U_L_0=[h_L_0;h_L_0*u_L_0];
U_R_0=[h_R_0;h_R_0*u_R_0];

W_int=zeros(2,N+1);
dh =zeros(1,N+1);
du =zeros(1,N+1);
dz =zeros(1,N+1);

for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if i<round(N*(x0-x_min)/(x_max-x_min))
        U(:,i)  =U_L_0;
        Z(i)    =Z_L_0;
    elseif i>round(N*(x0-x_min)/(x_max-x_min))
        U(:,i)  =U_R_0;
        Z(i+1)  =Z_R_0;
    else
        U(:,i)  =0.5*(U_L_0+U_R_0);
        Z(i)    =Z_L_0;
        Z(i+1)  =Z_R_0;
    end
end
%Godunov's Method
while Time<Tend && isreal(Time)
    %CFL condition
    for i=1:N
        [h_L(i),u_L(i),h_R(i),u_R(i)]=primitive_comp(U(:,i),Z(i),Z(i+1));
        a_L(i)=sqrt(g*h_L(i));
        a_R(i)=sqrt(g*h_R(i));
    end
    Smax=max(max(abs(u_L)+a_L),max(abs(u_R)+a_R));
    d_t=CFL*d_x/Smax;
    if Time+d_t >= Tend
        d_t = Tend-Time+1e-10;
    end
    %reconstruction (minmod limiter)
    for i=3:N-1
        if abs(Z(i+1)-Z(i))<ep && abs(Z(i+2)-Z(i+1))<ep && abs(Z(i)-Z(i-1))<ep
            dh(i) =minmod(Alpha_GRP*(h_R(i)-h_R(i-1))/d_x,(h_L(i+1)-h_L(i-1))/2.0/d_x,Alpha_GRP*(h_L(i+1)-h_L(i))/d_x);
            du(i) =minmod(Alpha_GRP*(u_R(i)-u_R(i-1))/d_x,(u_L(i+1)-u_L(i-1))/2.0/d_x,Alpha_GRP*(u_L(i+1)-u_L(i) )/d_x);
%             dh(i) =minmod(Alpha_GRP*(h_R(i)-h_R(i-1))/d_x,(W_int(1,i+1)-W_int(1,i))/1.0/d_x,Alpha_GRP*(h_L(i+1)-h_L(i))/d_x);
%             du(i)  =minmod(Alpha_GRP*(u_R(i)-u_R(i-1))/d_x,(W_int(2,i+1)-W_int(1,i))/1.0/d_x,Alpha_GRP*(u_L(i+1)-u_L(i))/d_x);
        end
    end
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
         if i==1
             F(:,1)=Riemann_solver_HLL(h_L(2),u_L(2),h_L(2),u_L(2));
         elseif i==N+1
             F(:,N+1)=Riemann_solver_HLL(h_R(N),u_R(N),h_R(N),u_R(N));
         else
             F(:,i)=Riemann_solver_HLL(h_R(i-1),u_R(i-1),h_L(i),u_L(i));
             %F(:,i)=Riemann_solver_HLL(h_R(i-1)+0.5*d_x*dh(i-1),u_R(i-1)+0.5*d_x*du(i-1),h_L(i)-0.5*d_x*dh(i),u_L(i)-0.5*d_x*du(i));
             %[F(:,i),W_int(:,i)]=GRP_solver(h_R(i-1)+0.5*d_x*dh(i-1),h_L(i)-0.5*d_x*dh(i),dh(i-1),dh(i),u_R(i-1)+0.5*d_x*du(i-1),u_L(i)-0.5*d_x*du(i),du(i-1),du(i),Z(i),Z(i),0.0,0.0,gama,d_t);
         end
    end
    %compute U in next step
    for i=1:N
        if abs(Z(i+1)-Z(i))<ep
            S=-g*0.5*(h_L(i)+h_R(i))*(Z(i+1)-Z(i));
        else
            S_tmp=(h_R(i)*u_R(i)^2+g*h_R(i)^2/2-h_L(i)*u_L(i)^2-g*h_L(i)^2/2);
            if (S_tmp/g/(Z(i)-Z(i+1))>max(h_L(i),h_R(i)))
                S=-g*max(h_L(i),h_R(i))*(Z(i+1)-Z(i));        
            elseif (S_tmp/g/(Z(i)-Z(i+1))<min(h_L(i),h_R(i)))
                S=-g*min(h_L(i),h_R(i))*(Z(i+1)-Z(i));
            else
                S=S_tmp;
           end
         end
        U(:,i)=U(:,i)+d_t/d_x*(F(:,i)-F(:,i+1))+d_t/d_x*[0;S];
    end
    Time=Time+d_t;
% if Time > 0.001
%    break;
% end
end
h = 0.5*(h_L+h_R);
u = 0.5*(u_L+u_R);

N_MAX = 3000;
d_xM=(x_max-x_min)/N_MAX;
W_exact = zeros(N_MAX,2);
%     load ./EXACT_SWE1.mat;
%     W_exact(:,1)=h_E';;
%     W_exact(:,2)=u_E';

%plot
%col = '-m';
col = '+k';
% h=figure(1);
% set(h,'position',[100 100 1150 750]);
subplot(1,2,1);
hold on
%plot(x_min:d_xM:x_max-d_xM,W_exact(:,1),'b','LineWidth',0.4);
plot(x_min:d_x:x_max-d_x,h+Z(1,1:N),col,'MarkerSize',6);%col,'LineWidth',1.0);
%xlabel('Position','FontWeight','bold');
%ylabel('Density','FontWeight','bold');
%ylim([90 160])
%ylim([0 2])
title('h+z');
set(gca,'box','on');
subplot(1,2,2);
hold on
% plot(x_min:d_xM:x_max-d_xM,W_exact(:,2),'b','LineWidth',0.4);
plot(x_min:d_x:x_max-d_x,u,col,'MarkerSize',6);%col,'LineWidth',1.0);
%xlabel('Position','FontWeight','bold');
%ylabel('Velocity','FontWeight','bold');
%ylim([0 2])
title('u');
set(gca,'box','on');
