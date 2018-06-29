clear;
clc;
%1D shock_tube by GRP for modified ghost fluid method
%state constant
global gama_s gama_g p0;
gama_s=1.4;
gama_g=1.6667;
global ep;
ep=1e-9;
x_min=0;
x_max=1;
N=200*1;
d_x=(x_max-x_min)/N;
x0s=0.2;
x0=0.2;
CFL=0.9;
Alpha=1.0;
%state value
Time=0;
Tend=0.1;
%Tend=0.15;
U=zeros(7,N);
FL=zeros(7,N+1);
FR=zeros(7,N+1);
d_U=zeros(7,N);
U_int=zeros(7,N+1);%U at cell interface
%initial condition
% lo_L_0   =0.9600000000000004;
% u_L_0    =1.083333333333333;
% p_L_0    =2.833333333333334;
% phi_sL_0 =0.4;
% lo_R_0   =1;
% u_R_0    =0;
% p_R_0    =0.303591608441676;
% phi_sR_0 =0.3;
load ./data/test1.mat;
%test begin
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if x(i)<=x0
        lo(i)   =lo_L_0;
        u(i)    =u_L_0;
        p(i)    =p_L_0;
        phi(i)  =x(i)-x_0;
    else
        lo(i)   =lo_R_0;
        u(i)    =u_R_0;
        p(i)    =p_R_0;
        phi(i)  =x(i)-x_0;
    end
    if phi(i)>0.0;
      gama(i)=gama_s;
    else
      gama(i)=gama_g;
    end
    E(i)=p(i)/(gama(i)-1)+0.5*lo(i)*u(i)^2;
    if phi(i)>0.0;
      U(:,i)=[lo(i);lo(i)*u(i);E(i);0;0;0];
    else
      U(:,i)=[0;0;0;lo(i);lo(i)*u(i);E(i)];
    end
end
%Godunov's Method
while Time<Tend && isreal(Time)
    %reconstruction (minmod limiter)
    for i=2:N-1
        d_U(:,i)=minmod(Alpha*(U(:,i)-U(:,i-1))/d_x,(U_int(:,i+1)-U_int(:,i))/d_x,Alpha*(U(:,i+1)-U(:,i))/d_x);
        %d_U(:,i)=minmod(Alpha*(U(:,i)-U(:,i-1))/d_x,(U(:,i+1)-U(:,i-1))/2.0/d_x,Alpha*(U(:,i+1)-U(:,i))/d_x);
        %d_U(:,i)=0;
    end
    %CFL condition
    for i=1:N
        a_g(i)=sqrt(gama_g*p_g(i)/lo_g(i));
        a_s(i)=sqrt(gama_s*(p_s(i)+p0)/lo_s(i));
    end
    Smax=max(max(abs(u_g)+a_g),max(abs(u_s)+a_s));
    d_t=CFL*d_x/Smax;
    if Time+d_t >= Tend
        d_t = Tend-Time+1e-10;
    end
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
        C_U=zeros(7,7);
         if i==1
             [lo_gL,u_gL,p_gL,phi_gL,lo_sL,u_sL,p_sL,phi_sL]=primitive_comp(U(:,1));
             [lo_gR,u_gR,p_gR,phi_gR,lo_sR,u_sR,p_sR,phi_sR]=primitive_comp(U(:,1));
             d_UL=zeros(7,1);
             d_UR=zeros(7,1);
             UL_mid=U(:,1);
             UR_mid=U(:,1);
         elseif i==N+1
             [lo_gL,u_gL,p_gL,phi_gL,lo_sL,u_sL,p_sL,phi_sL]=primitive_comp(U(:,N));
             [lo_gR,u_gR,p_gR,phi_gR,lo_sR,u_sR,p_sR,phi_sR]=primitive_comp(U(:,N));
             d_UL=zeros(7,1);
             d_UR=zeros(7,1);
             UL_mid=U(:,N);
             UR_mid=U(:,N);
         else
             [lo_gL,u_gL,p_gL,phi_gL,lo_sL,u_sL,p_sL,phi_sL]=primitive_comp(U(:,i-1)+0.5*d_x*d_U(:,i-1));
             [lo_gR,u_gR,p_gR,phi_gR,lo_sR,u_sR,p_sR,phi_sR]=primitive_comp(U(:,i)-0.5*d_x*d_U(:,i));
             d_UL=d_U(:,i-1);
             d_UR=d_U(:,i);
             C_U(1:4,1:4)=quasilinear(lo_g(i-1),u_g(i-1),p_g(i-1),p_g(i-1),u_s(i-1),'g','C');
             C_U([1 5:7],[1 5:7])=quasilinear(lo_s(i-1),u_s(i-1),p_s(i-1),p_g(i-1),u_s(i-1),'s','C');
             UL_mid=U(:,i-1)-0.5*d_t*C_U*d_UL;
             C_U(1:4,1:4)=quasilinear(lo_g(i),u_g(i),p_g(i),p_g(i),u_s(i),'g','C');
             C_U([1 5:7],[1 5:7])=quasilinear(lo_s(i),u_s(i),p_s(i),p_g(i),u_s(i),'s','C');
             UR_mid=U(:,i)-0.5*d_t*C_U*d_UR;
         end
       [FL(:,i),FR(:,i),U_int(:,i)]=GRP_solver_HLLC(lo_gL,lo_gR,p_gL,p_gR,u_gL,u_gR,lo_sL,lo_sR,p_sL,p_sR,u_sL,u_sR,phi_sL,phi_sR,d_UL,d_UR,UL_mid,UR_mid,d_t/d_x,d_t);
    end
    %compute U in next step
    for i=1:N
        U(:,i)=U(:,i)+d_t/d_x*(FR(:,i)-FL(:,i+1));
        [lo_g(i),u_g(i),p_g(i),phi_g(i),lo_s(i),u_s(i),p_s(i),phi_s(i)]=primitive_comp(U(:,i));
    end
    Time=Time+d_t;
% if Time > 5*d_t
%     break;
% end
end
W_exact = zeros(N,8);
W_exact(:,2)=phi_s';
W_exact(:,3)=lo_s';
W_exact(:,4)=u_s';
W_exact(:,5)=p_s';
W_exact(:,6)=lo_g';
W_exact(:,7)=u_g';
W_exact(:,8)=p_g';
load ../test/test3.exact;
for i=1:N
     W_exact(i,:) = test3(ceil(i/(N/300)),:);
end

%plot
col = '-r';
figure(1);
subplot(2,2,1);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,3),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,lo_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Density-solid','FontWeight','bold');
ylim([min(lo_s)-0.00001 max(lo_s)+0.00001])
subplot(2,2,2);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,4),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,u_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Velocity-solid','FontWeight','bold');
ylim([min(u_s)-0.00001 max(u_s)+0.00001])
subplot(2,2,3);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,5),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,p_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Pressure-solid','FontWeight','bold');
subplot(2,2,4);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,2),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,phi_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Porosity-solid','FontWeight','bold');
figure(2);
subplot(2,2,1);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,6),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,lo_g,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Density-gas','FontWeight','bold');
subplot(2,2,2);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,7),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,u_g,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Velocity-gas','FontWeight','bold');
subplot(2,2,3);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,8),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,p_g,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Pressure-gas','FontWeight','bold');
subplot(2,2,4);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,8)./W_exact(:,6).^gama_g,'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,p_g./lo_g.^gama_g,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Entropy-gas','FontWeight','bold');
figure(3)
subplot(3,1,1);
hold on
plot(x_min:d_x:x_max-d_x,(1-W_exact(:,2)).*W_exact(:,6).*(W_exact(:,7)-W_exact(:,4)),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,phi_g.*lo_g.*(u_g-u_s),col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Riemann_inv-Q','FontWeight','bold');
ylim([min(phi_g.*lo_g.*(u_g-u_s))-0.00001 max(phi_g.*lo_g.*(u_g-u_s))+0.00001])
subplot(3,1,2);
hold on
plot(x_min:d_x:x_max-d_x,(1-W_exact(:,2)).*W_exact(:,6).*(W_exact(:,7)-W_exact(:,4)).^2+(1-W_exact(:,2)).*W_exact(:,8)+W_exact(:,2).*W_exact(:,5),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,phi_g.*lo_g.*(u_g-u_s).^2+phi_g.*p_g+phi_s.*p_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Riemann_inv-P','FontWeight','bold');
ylim([min(phi_g.*lo_g.*(u_g-u_s).^2+phi_g.*p_g+phi_s.*p_s)-0.00001 max(phi_g.*lo_g.*(u_g-u_s).^2+phi_g.*p_g+phi_s.*p_s)+0.00001])
subplot(3,1,3);
hold on
plot(x_min:d_x:x_max-d_x,0.5*(W_exact(:,7)-W_exact(:,4)).^2+gama_g/(gama_g-1)*W_exact(:,8)./W_exact(:,6),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,0.5*(u_g-u_s).^2+gama_g/(gama_g-1)*p_g./lo_g,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Riemann_inv-H','FontWeight','bold');
