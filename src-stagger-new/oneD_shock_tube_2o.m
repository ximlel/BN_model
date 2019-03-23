clear;
clc;
%1D shock_tube by 2-order staggered Schemes for BN model
%state constant
global gama_s gama_g p0;
gama_s=1.4;
gama_g=1.4;
p0=0;
global ep;
ep=1e-11;
x_min=0;
x_max=1;
N=300*1;
d_x=(x_max-x_min)/N;
x0=0.5;
CFL=0.2;
Alpha_GRP=1;
%state value
Time=0;
Tend=0.1;
%Tend=0.15;
Alpha=zeros(1,N+1);
U=zeros(6,N);
F=zeros(6,N+1);
RI=zeros(6,N); %Riemann invariants
d_RI=zeros(6,N);
d_Alpha=zeros(1,N+1);
U_lo_sL=zeros(1,N);
U_lo_sR=zeros(1,N);
%initial condition
% lo_gL_0  =1;
% u_gL_0   =2;
% p_gL_0   =1;
% lo_sL_0  =2;
% u_sL_0   =0.3;
% p_sL_0   =5;
% phi_sL_0 =0.8;
% lo_gR_0  =0.1941934235006083;
% u_gR_0   =2.801188129642115;
% p_gR_0   =0.1008157360849781;
% lo_sR_0  =2;
% u_sR_0   =0.3;
% p_sR_0   =12.85675006887399;
% phi_sR_0 =0.3;
load ../test/test2.mat;
phi_gL_0=1.0-phi_sL_0;
phi_gR_0=1.0-phi_sR_0;
E_gL_0=p_gL_0/(gama_g-1)+0.5*lo_gL_0*u_gL_0^2;
E_sL_0=(p_sL_0+gama_s*p0)/(gama_s-1)+0.5*lo_sL_0*u_sL_0^2;
U_L_0=[phi_gL_0*lo_gL_0;phi_gL_0*lo_gL_0*u_gL_0;phi_gL_0*E_gL_0;phi_sL_0*lo_sL_0;phi_sL_0*lo_sL_0*u_sL_0;phi_sL_0*E_sL_0];
E_gR_0=p_gR_0/(gama_g-1)+0.5*lo_gR_0*u_gR_0^2;
E_sR_0=(p_sR_0+gama_s*p0)/(gama_s-1)+0.5*lo_sR_0*u_sR_0^2;
U_R_0=[phi_gR_0*lo_gR_0;phi_gR_0*lo_gR_0*u_gR_0;phi_gR_0*E_gR_0;phi_sR_0*lo_sR_0;phi_sR_0*lo_sR_0*u_sR_0;phi_sR_0*E_sR_0];

E_gL_1=p_gL_1/(gama_g-1)+0.5*lo_gL_1*u_gL_1^2;
E_sL_1=(p_sL_1+gama_s*p0)/(gama_s-1)+0.5*lo_sL_1*u_sL_1^2;
U_L_1=[phi_gL_0*lo_gL_1;phi_gL_0*lo_gL_1*u_gL_1;phi_gL_0*E_gL_1;phi_sL_0*lo_sL_1;phi_sL_0*lo_sL_1*u_sL_1;phi_sL_0*E_sL_1];
E_gR_1=p_gR_1/(gama_g-1)+0.5*lo_gR_1*u_gR_1^2;
E_sR_1=(p_sR_1+gama_s*p0)/(gama_s-1)+0.5*lo_sR_1*u_sR_1^2;
U_R_1=[phi_gR_0*lo_gR_1;phi_gR_0*lo_gR_1*u_gR_1;phi_gR_0*E_gR_1;phi_sR_0*lo_sR_1;phi_sR_0*lo_sR_1*u_sR_1;phi_sR_0*E_sR_1];
%test begin
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if i<round(N*x0/(x_max-x_min))
        U(:,i) =U_L_0;
        Alpha(i) =phi_sL_0;
    elseif i>round(N*x0/(x_max-x_min))
        U(:,i) =U_R_0;
        Alpha(i+1) =phi_sR_0;
    else
%        U(:,i) =0.5*(U_L_0+U_R_0);
        U(:,i) =0.5*(U_L_1+U_R_1);
        Alpha(i) =phi_sL_0;
        Alpha(i+1) =phi_sR_0;
    end
end
Alpha_next=Alpha;
%Godunov's Method
while Time<Tend && isreal(Time)
    %CFL condition
    for i=1:N
        [lo_gL(i),u_gL(i),p_gL(i),lo_sL(i),u_sL(i),p_sL(i),lo_gR(i),u_gR(i),p_gR(i),lo_sR(i),u_sR(i),p_sR(i)]=primitive_comp(U(:,i),Alpha(i),Alpha(i+1),0.5,0.5);
        RI(:,i)=U2RI_cal(Alpha(i),lo_gL(i),u_gL(i),p_gL(i),u_sL(i),p_sL(i),lo_sL(i));
        a_gL(i)=sqrt(gama_g*p_gL(i)/lo_gL(i));
        a_sL(i)=sqrt(gama_s*(p_sL(i)+p0)/lo_sL(i));
        a_gR(i)=sqrt(gama_g*p_gR(i)/lo_gR(i));
        a_sR(i)=sqrt(gama_s*(p_sR(i)+p0)/lo_sR(i));
    end
    Smax=max([max(abs(u_gL)+a_gL),max(abs(u_sL)+a_sL),max(abs(u_gR)+a_gR),max(abs(u_sR)+a_sR)]);
    d_t=CFL*d_x/Smax;
    if Time+d_t >= Tend
        d_t = Tend-Time+1e-10;
    end
    %reconstruction (minmod limiter)
    for i=2:N
%        d_Alpha(i)=minmod2((Alpha(:,i)-Alpha(:,i-1))/d_x,(Alpha(:,i+1)-Alpha(:,i))/d_x);
         d_Alpha(i)=minmod(Alpha_GRP*(Alpha(:,i)-Alpha(:,i-1))/d_x,(Alpha(:,i+1)-Alpha(:,i-1))/2.0/d_x,Alpha_GRP*(Alpha(:,i+1)-Alpha(:,i))/d_x);
% if Time/d_t < 50
%     d_Alpha(i) = 0;
% end
if Time/d_t < 50
    AG = Time/d_t/50;
end
         d_Alpha(i)=minmod(AG*(Alpha(:,i)-Alpha(:,i-1))/d_x,(Alpha(:,i+1)-Alpha(:,i-1))/2.0/d_x,AG*(Alpha(:,i+1)-Alpha(:,i))/d_x);
    end
    for i=4:N-3
%        d_RI(:,i)=minmod2((RI(:,i)-RI(:,i-1))/d_x,(RI(:,i+1)-RI(:,i))/d_x);
         d_RI(:,i)=minmod(Alpha_GRP*(RI(:,i)-RI(:,i-1))/d_x,(RI(:,i+1)-RI(:,i-1))/2.0/d_x,Alpha_GRP*(RI(:,i+1)-RI(:,i))/d_x);
        if d_Alpha(i-3)>ep||d_Alpha(i-2)>ep||d_Alpha(i-1)>ep||d_Alpha(i)>ep||d_Alpha(i+1)>ep||d_Alpha(i+2)>ep||d_Alpha(i+3)>ep
            d_RI(:,i)=0.0;
        end
    end
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
         if i==1
             Alpha_int=Alpha(1);
             [lo_gL_i,u_gL_i,p_gL_i,u_sL_i,p_sL_i(1),lo_sL_i]=RI2U_cal(Alpha_int,RI(:,1),lo_gL(1));
             [lo_gR_i,u_gR_i,p_gR_i,u_sR_i,p_sR_i(1),lo_sR_i]=RI2U_cal(Alpha_int,RI(:,1),lo_gL(1));
         elseif i==N+1
             Alpha_int=Alpha(N+1);
             [lo_gL_i,u_gL_i,p_gL_i,u_sL_i,p_sL_i(N+1),lo_sL_i]=RI2U_cal(Alpha_int,RI(:,N),lo_gR(N));
             [lo_gR_i,u_gR_i,p_gR_i,u_sR_i,p_sR_i(N+1),lo_sR_i]=RI2U_cal(Alpha_int,RI(:,N),lo_gR(N));
         else
             Alpha_int=Alpha(i);
             [lo_gL_i,u_gL_i,p_gL_i,u_sL_i,p_sL_i(i),lo_sL_i]=RI2U_cal(Alpha_int,RI(:,i-1)+0.5*d_x*d_RI(:,i-1),lo_gR(i-1));
             [lo_gR_i,u_gR_i,p_gR_i,u_sR_i,p_sR_i(i),lo_sR_i]=RI2U_cal(Alpha_int,  RI(:,i)-0.5*d_x*d_RI(:,i),  lo_gL(i));
         end
       F(1:3,i)=Riemann_solver_Exact(lo_gL_i,lo_gR_i,p_gL_i,   p_gR_i,   u_gL_i,u_gR_i,1-Alpha_int,gama_g,0.0);
       F(4:6,i)=Riemann_solver_Exact(lo_sL_i,lo_sR_i,p_sL_i(i),p_sR_i(i),u_sL_i,u_sR_i,  Alpha_int,gama_s,0.0);
    end
    for i=1:N
        if i<N
            rho_s(i) =0.5*(lo_sR(i)+lo_sL(i+1));
            arho_s(i)=0.5*(lo_sR(i)+lo_sL(i+1))*Alpha(i+1);           
        end
        F_rho_s(i) =lo_sL(i)*u_sL(i);
        if u_sL(i)>0
            F_arho_s(i)=  (Alpha(i)+0.5*d_x*d_Alpha(i))  *lo_sL(i)*u_sL(i);
        else
            F_arho_s(i)=(Alpha(i+1)-0.5*d_x*d_Alpha(i+1))*lo_sR(i)*u_sR(i);
        end
    end
    for i=1:N-1
        arho_s(i)=arho_s(i)+d_t/d_x*(F_arho_s(i)-F_arho_s(i+1));
        rho_s(i) = rho_s(i)+d_t/d_x*(F_rho_s(i) -F_rho_s(i+1));
        Alpha_next(i+1)=arho_s(i)/rho_s(i);
    end
    %compute U in next step
    for i=1:N
      if abs(Alpha(i+1)-Alpha(i))<ep
          S=0.5*(p_gL(i)+p_gR(i))*(Alpha(i+1)-Alpha(i));
      else
          S_tmp=Alpha(i+1)*p_sL_i(i+1)-Alpha(i)*p_sR_i(i);
%           if (S_tmp/(Alpha(i+1)-Alpha(i))>max(p_gL(i),p_gR(i)))
%               S=max(p_gL(i),p_gR(i))*(Alpha(i+1)-Alpha(i));
%           elseif (S_tmp/(Alpha(i+1)-Alpha(i))<min(p_gL(i),p_gR(i)))
%               S=min(p_gL(i),p_gR(i))*(Alpha(i+1)-Alpha(i));
%           else
              S=S_tmp;
%           end
      end
        U(:,i)=U(:,i)+d_t/d_x*(F(:,i)-F(:,i+1))+d_t/d_x*[0;-S;-S*u_sL(i);0;S;S*u_sL(i)];
        area_L=0.5+u_sL(i)*d_t/d_x;
        area_R=1.0-area_L;
        [lo_gL(i),u_gL(i),p_gL(i),lo_sL(i),u_sL(i),p_sL(i),lo_gR(i),u_gR(i),p_gR(i),lo_sR(i),u_sR(i),p_sR(i)]=primitive_comp(U(:,i),Alpha(i),Alpha(i+1),area_L,area_R);
        phi_sL=Alpha_next(i);
        phi_sR=Alpha_next(i+1);
        [lo_gL(i),u_gL(i),p_gL(i),p_sL(i)]=Riemann_inv(Alpha(i),  lo_gL(i),u_gL(i),p_gL(i),u_sL(i),p_sL(i),phi_sL);   
        [lo_gR(i),u_gR(i),p_gR(i),p_sR(i)]=Riemann_inv(Alpha(i+1),lo_gR(i),u_gR(i),p_gR(i),u_sR(i),p_sR(i),phi_sR);   
        phi_gL=1.0-phi_sL;
        phi_gR=1.0-phi_sR;
        E_gL=p_gL(i)/(gama_g-1)+0.5*lo_gL(i)*u_gL(i)^2;
        E_sL=(p_sL(i)+gama_s*p0)/(gama_s-1)+0.5*lo_sL(i)*u_sL(i)^2;
        U(:,i)=       0.5*[phi_gL*lo_gL(i);phi_gL*lo_gL(i)*u_gL(i);phi_gL*E_gL;phi_sL*lo_sL(i);phi_sL*lo_sL(i)*u_sL(i);phi_sL*E_sL];
        E_gR=p_gR(i)/(gama_g-1)+0.5*lo_gR(i)*u_gR(i)^2;
        E_sR=(p_sR(i)+gama_s*p0)/(gama_s-1)+0.5*lo_sR(i)*u_sR(i)^2;
        U(:,i)=U(:,i)+0.5*[phi_gR*lo_gR(i);phi_gR*lo_gR(i)*u_gR(i);phi_gR*E_gR;phi_sR*lo_sR(i);phi_sR*lo_sR(i)*u_sR(i);phi_sR*E_sR];
    end
    Alpha=Alpha_next;
    Time=Time+d_t
%     if Time > 100*d_t
%         break;
%     end
end
lo_g = 0.5*(lo_gL+lo_gR);
p_g  = 0.5*(p_gL +p_gR);
u_g  = 0.5*(u_gL +u_gR);
lo_s = 0.5*(lo_sL+lo_sR);
p_s  = 0.5*(p_sL +p_sR);
u_s  = 0.5*(u_sL +u_sR);
phi_sL=Alpha(1,1:N);
phi_gL=1-phi_sL;
phi_sR=Alpha(1,2:N+1);
phi_gR=1-phi_sR;
phi_s=0.5*(phi_sL+phi_sR);
phi_g=1-phi_s;
eta=0.5*(p_gL./lo_gL.^gama_g+p_gR./lo_gR.^gama_g);
Q_inv=0.5*(phi_gL.*lo_gL.*(u_gL-u_sL)+phi_gR.*lo_gR.*(u_gR-u_sR));
P_inv=0.5*(phi_gL.*lo_gL.*(u_gL-u_sL).^2+phi_gL.*p_gL+phi_sL.*p_sL+phi_gR.*lo_gR.*(u_gR-u_sR).^2+phi_gR.*p_gR+phi_sR.*p_sR);
H_inv=0.5*(0.5*(u_gL-u_sL).^2+gama_g/(gama_g-1)*p_gL./lo_gL+0.5*(u_gR-u_sR).^2+gama_g/(gama_g-1)*p_gR./lo_gR);
W_exact = zeros(N,8);
W_exact(:,2)=phi_s';
W_exact(:,3)=lo_s';
W_exact(:,4)=u_s';
W_exact(:,5)=p_s';
W_exact(:,6)=lo_g';
W_exact(:,7)=u_g';
W_exact(:,8)=p_g';
load ../test/test2.exact;
for i=1:N
     W_exact(i,:) = test2(ceil(i/(N/300)),:);
end
%plot
%col = '+k';
%col = 'or';
%col = '*m';
col = '+b';
h1=figure(1);
set(h1,'position',[100 100 800 600]);
subplot(2,2,1);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,3),'k','LineWidth',0.4);
plot(x,lo_s,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Density-solid','FontWeight','bold');
% ylim([min(lo_s)-0.00001 max(lo_s)+0.00001])
title('固体相—密度')
subplot(2,2,2);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,4),'k','LineWidth',0.4);
plot(x,u_s,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Velocity-solid','FontWeight','bold');
% ylim([min(u_s)-0.00001 max(u_s)+0.00001])
title('固体相—速度')
subplot(2,2,3);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,5),'k','LineWidth',0.4);
plot(x,p_s,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Pressure-solid','FontWeight','bold');
title('固体相—压力')
subplot(2,2,4);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,2),'k','LineWidth',0.4);
plot(x,phi_s,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Porosity-solid','FontWeight','bold');
title('固体相—体积分数')
h2=figure(2);
set(h2,'position',[100 100 800 600]);
subplot(2,2,1);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,6),'k','LineWidth',0.4);
plot(x,lo_g,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Density-gas','FontWeight','bold');
title('气体相—密度')
subplot(2,2,2);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,7),'k','LineWidth',0.4);
plot(x,u_g,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Velocity-gas','FontWeight','bold');
title('气体相—速度')
subplot(2,2,3);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,8),'k','LineWidth',0.4);
plot(x,p_g,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Pressure-gas','FontWeight','bold');
title('气体相—压力')
subplot(2,2,4);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,8)./W_exact(:,6).^gama_g,'k','LineWidth',0.4);
plot(x,eta,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Entropy-gas','FontWeight','bold');
% ylim([min(eta)-0.00001 max(eta)+0.00001])
title('气体相—熵')
h3=figure(3)
set(h3,'position',[100 100 800 600]);
subplot(3,1,1);
hold on
plot(x_min:d_x:x_max-d_x,(1-W_exact(:,2)).*W_exact(:,6).*(W_exact(:,7)-W_exact(:,4)),'k','LineWidth',0.4);
plot(x,Q_inv,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Riemann_inv-Q','FontWeight','bold');
% ylim([min(Q_inv)-0.00001 max(Q_inv)+0.00001])
title('黎曼不变量—Q')
subplot(3,1,2);
hold on
plot(x_min:d_x:x_max-d_x,(1-W_exact(:,2)).*W_exact(:,6).*(W_exact(:,7)-W_exact(:,4)).^2+(1-W_exact(:,2)).*W_exact(:,8)+W_exact(:,2).*W_exact(:,5),'k','LineWidth',0.4);
plot(x,P_inv,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Riemann_inv-P','FontWeight','bold');
% ylim([min(P_inv)-0.00001 max(P_inv)+0.00001])
title('黎曼不变量—P')
subplot(3,1,3);
hold on
plot(x_min:d_x:x_max-d_x,0.5*(W_exact(:,7)-W_exact(:,4)).^2+gama_g/(gama_g-1)*W_exact(:,8)./W_exact(:,6),'k','LineWidth',0.4);
plot(x,H_inv,col,'MarkerSize',4);
% ylim([min(H_inv)-0.00001 max(H_inv)+0.00001])
% xlabel('Position','FontWeight','bold');
% ylabel('Riemann_inv-H','FontWeight','bold');
title('黎曼不变量—H')
