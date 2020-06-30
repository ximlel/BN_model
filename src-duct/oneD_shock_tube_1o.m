clear;
clc;
%1D shock_tube by Roe Schemes for duct-flow model
%state constant
global gama;
gama=1.23;
global ep;
ep=1e-6;
x_min=0;
x_max=0.06;
x0=0.02;
%state value
Time=0;
Tend=6.3e-6;
N=110*1;
% gama=1.4;
% x_max=1;
% x0=0.7;
% Tend=0.2;
% N=400*1;
% %it_max = 500;

CFL=0.5;
d_x=(x_max-x_min)/N;
A=zeros(1,N+1);
U=zeros(3,N);
F=zeros(3,N+1);
%initial condition
% lo_L_0  =169.34;
% u_L_0   =0;
% p_L_0   =2.96e8;
% lo_R_0  =0.76278;
% u_R_0   =0;
% p_R_0   =1e5;
% phi_L_0 =1.0;
% phi_R_0 =0.25;

lo_L_0  =151.13;
u_L_0   =212.31;
p_L_0   =2.4836e8;
lo_R_0  =95.199;
u_R_0   =1348.2;
p_R_0   =1.4067e8;
phi_L_0 =1.0;
phi_R_0 =0.25;

% lo_L_0  =1.0555;
% u_L_0   =-1.0651;
% p_L_0   =1.5;
% lo_R_0  =1;
% u_R_0   =-1;
% p_R_0   =1;
% phi_L_0 =1.0;
% phi_R_0 =1.25;

% lo_L_0  =0.6894;
% u_L_0   =-1.6941;
% p_L_0   =1.5;
% lo_R_0  =1;
% u_R_0   =-0.5;
% p_R_0   =1;
% phi_L_0 =1.0;
% phi_R_0 =1.25;
% load ../test/test_new1_pi.mat;
%test begin
E_L_0=p_L_0/(gama-1)+0.5*lo_L_0*u_L_0^2;
U_L_0=[phi_L_0*lo_L_0;phi_L_0*lo_L_0*u_L_0;phi_L_0*E_L_0];
E_R_0=p_R_0/(gama-1)+0.5*lo_R_0*u_R_0^2;
U_R_0=[phi_R_0*lo_R_0;phi_R_0*lo_R_0*u_R_0;phi_R_0*E_R_0];

lo_L_1 = 146.5226;
u_L_1 = 210.3386;
p_L_1 = 2.4775e8;
RI_L=zeros(3,1);
RI_L(1)=phi_L_0*lo_L_0*u_L_0;
RI_L(2)=p_L_0/lo_L_0^gama;
RI_L(3)=p_L_0/lo_L_0*gama/(gama-1.0)+0.5*u_L_0*u_L_0;
% RI_R=zeros(3,1);
% RI_R(1)=phi_R_0*lo_R_0*u_R_0;
% RI_R(2)=p_R_0/lo_R_0^gama;
% RI_R(3)=p_R_0/lo_R_0*gama/(gama-1.0)+0.5*u_R_0*u_R_0;
RI = RI_L;% RI=0.5*(RI_L+RI_R);
%[lo_L_1,u_L_1,p_L_1]=RI2U_cal(phi_L_0,RI,lo_L_0);
[lo_R_1,u_R_1,p_R_1]=RI2U_cal(phi_R_0,RI,lo_R_0);
E_L_1=p_L_1/(gama-1)+0.5*lo_L_1*u_L_1^2;
U_L_1=[phi_L_0*lo_L_1;phi_L_0*lo_L_1*u_L_1;phi_L_0*E_L_1];
E_R_1=p_R_1/(gama-1)+0.5*lo_R_1*u_R_1^2;
U_R_1=[phi_R_0*lo_R_1;phi_R_0*lo_R_1*u_R_1;phi_R_0*E_R_1];

W_int=zeros(4,N+1);
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if i<round(N*x0/(x_max-x_min))
        U(:,i)  =U_L_0;
        A(i)    =phi_L_0;
    elseif i>round(N*x0/(x_max-x_min))
        U(:,i)  =U_R_0;
        A(i+1)  =phi_R_0;
    else
        %U(:,i)  =0.5*(U_L_1+U_R_1);
        U(:,i)  =0.5*(U_L_0+U_R_0);
        A(i)    =phi_L_0;
        A(i+1)  =phi_R_0;
    end
end
%Godunov's Method
while Time<Tend && isreal(Time)
    %CFL condition
    for i=1:N
        [lo_L(i),u_L(i),p_L(i),lo_R(i),u_R(i),p_R(i)]=primitive_comp(U(:,i),A(i),A(i+1));
        a_L(i)=sqrt(gama*p_L(i)/lo_L(i));
        a_R(i)=sqrt(gama*p_R(i)/lo_R(i));
    end
    Smax=max(max(abs(u_L)+a_L),max(abs(u_R)+a_R));
    d_t=CFL*d_x/Smax;
    if Time+d_t >= Tend
        d_t = Tend-Time+1e-10;
    end
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
         if i==1
             F(:,1)=Riemann_solver_Roe(lo_L(1),lo_L(1),p_L(1),p_L(1),u_L(1),u_L(1),A(1));
         elseif i==N+1
             F(:,N+1)=Riemann_solver_Roe(lo_R(N),lo_R(N),p_R(N),p_R(N),u_R(N),u_R(N),A(N+1));
         else
             F(:,i)=Riemann_solver_Roe(lo_R(i-1),lo_L(i),p_R(i-1),p_L(i),u_R(i-1),u_L(i),A(i));
            [F(:,i),W_int(:,i)]=GRP_solver(lo_R(i-1),lo_L(i),0.0,0.0,u_R(i-1),u_L(i),0.0,0.0,p_R(i-1),p_L(i),0.0,0.0,A(i),A(i),0.0,0.0,gama,0.0);
         end
    end
    %compute U in next step
    for i=1:N
        if abs(A(i+1)-A(i))<ep
            S=0.5*(p_L(i)+p_R(i))*(A(i+1)-A(i));
        else
            S_tmp=A(i+1)*(lo_R(i)*u_R(i)^2+p_R(i))-A(i)*(lo_L(i)*u_L(i)^2+p_L(i));
            if (S_tmp/(A(i+1)-A(i))>max(p_L(i),p_R(i)))
                S=max(p_L(i),p_R(i))*(A(i+1)-A(i));
            elseif (S_tmp/(A(i+1)-A(i))<min(p_L(i),p_R(i)))
                S=min(p_L(i),p_R(i))*(A(i+1)-A(i));
            else
                S=S_tmp;
            end
        end
        U(:,i)=U(:,i)+d_t/d_x*(F(:,i)-F(:,i+1))+d_t/d_x*[0;S;0];
    end
    Time=Time+d_t
% if Time > 0.5*d_t
%    break;
% end
end
lo = 0.5*(lo_L+lo_R);
p  = 0.5*(p_L +p_R);
u  = 0.5*(u_L +u_R);
eta= 0.5*(p_L./lo_L.^gama+p_R./lo_R.^gama);
% eta_E=eta;
% p_E=p;
% lo_E=lo;
% u_E=u;
% save('EXACT_duct2.mat','eta_E','p_E','lo_E','u_E');

N_MAX = 1000;
d_xM=(x_max-x_min)/N_MAX;
W_exact = zeros(N_MAX,4);
load ./EXACT_duct2.mat;
W_exact(:,1)=eta_E';%A(1:N)';
W_exact(:,2)=p_E';
W_exact(:,3)=lo_E';
W_exact(:,4)=u_E';
%load ../test/test_new1_pi.exact;
%for i=1:N
%     W_exact(i,:) = test_new1_pi(ceil(i/(N/300)),:);
%end

%plot
col = '-m';
h=figure(1);
set(h,'position',[100 100 1150 750]);
subplot(2,2,1);
hold on
hs(3)=plot(x_min:d_xM:x_max-d_xM,W_exact(:,3),'b','LineWidth',0.4);
hs(1)=plot(x_min:d_x:x_max-d_x,lo,'xr','MarkerSize',6);%col,'LineWidth',1.0);
%xlabel('Position','FontWeight','bold');
%ylabel('Density','FontWeight','bold');
ylim([90 160])
% ylim([0 200])
hs(2)=plot(0,-100,'+k'); 
legend(hs,'Godunov solution','GRP solution','Exact solution');
title('Density');
set(gca,'box','on');
subplot(2,2,3);
hold on
plot(x_min:d_xM:x_max-d_xM,W_exact(:,4),'b','LineWidth',0.4);
plot(x_min:d_x:x_max-d_x,u,'xr','MarkerSize',6);%col,'LineWidth',1.0);
%xlabel('Position','FontWeight','bold');
%ylabel('Velocity','FontWeight','bold');
title('Velocity');
set(gca,'box','on');
subplot(2,2,2);
hold on
plot(x_min:d_xM:x_max-d_xM,W_exact(:,2),'b','LineWidth',0.4);
plot(x_min:d_x:x_max-d_x,p,'xr','MarkerSize',6);%,col,'LineWidth',1.0);
%xlabel('Position','FontWeight','bold');
%ylabel('Pressure','FontWeight','bold');
title('Pressure');
set(gca,'box','on');
subplot(2,2,4);
hold on
plot(x_min:d_xM:x_max-d_xM,W_exact(:,1),'b','LineWidth',0.4);
% plot(x_min:d_x:x_max-d_x,A(1:N),col,'LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,eta,'xr','MarkerSize',6);%,col,'LineWidth',1.0);
%xlabel('Position','FontWeight','bold');
%ylabel('Entropy','FontWeight','bold');
title('Entropy');
set(gca,'box','on');
ylim([5.08*10^5 5.2*10^5])
% ylim([0 12*10^5])
% ylim([min(eta)-0.00001 max(eta)+0.00001])
