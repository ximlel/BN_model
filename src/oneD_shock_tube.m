clear;
clc;
%1D shock_tube by HLLC-GRP Schemes for BN model
%state constant
global gama_s gama_g p0;
gama_s=1.4;
gama_g=1.4;
p0=0;
global ep;
ep=1e-9;
x_min=0;
x_max=1;
N=300*1;
d_x=(x_max-x_min)/N;
x0=0.5;
CFL=0.4;
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
load ../test/test_new1.mat;
%test begin
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if x(i)<=x0
        lo_g(i)   =lo_gL_0;
        u_g(i)    =u_gL_0;
        p_g(i)    =p_gL_0;
        lo_s(i)   =lo_sL_0;
        u_s(i)    =u_sL_0;
        p_s(i)    =p_sL_0;
        phi_s(i)  =phi_sL_0;
        phi_g(i)  =1.0-phi_s(i);
    else
        lo_g(i)   =lo_gR_0;
        u_g(i)    =u_gR_0;
        p_g(i)    =p_gR_0;
        lo_s(i)   =lo_sR_0;
        u_s(i)    =u_sR_0;
        p_s(i)    =p_sR_0;
        phi_s(i)  =phi_sR_0;
        phi_g(i)  =1.0-phi_s(i);
    end
    E_g(i)=p_g(i)/(gama_g-1)+0.5*lo_g(i)*u_g(i)^2;
    E_s(i)=(p_s(i)+gama_s*p0)/(gama_s-1)+0.5*lo_s(i)*u_s(i)^2;
    U(:,i)=[phi_g(i)*lo_g(i);phi_g(i)*lo_g(i)*u_g(i);phi_g(i)*E_g(i);phi_s(i)*lo_s(i);phi_s(i)*lo_s(i)*u_s(i);phi_s(i)*E_s(i);phi_s(i)];
end
%Godunov's Method
while Time<Tend && isreal(Time)
    %reconstruction (minmod limiter)
    for i=2:N-1
        %d_U(:,i)=minmod(Alpha*(U(:,i)-U(:,i-1))/d_x,(U_int(:,i+1)-U_int(:,i))/d_x,Alpha*(U(:,i+1)-U(:,i))/d_x);
        d_U(:,i)=minmod(Alpha*(U(:,i)-U(:,i-1))/d_x,(U(:,i+1)-U(:,i-1))/2.0/d_x,Alpha*(U(:,i+1)-U(:,i))/d_x);
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
    %Riemann Reoblem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
         if i==1
             [lo_gL u_gL p_gL phi_gL lo_sL u_sL p_sL phi_sL]=primitive_comp(U(:,1));
             [lo_gR u_gR p_gR phi_gR lo_sR u_sR p_sR phi_sR]=primitive_comp(U(:,1));
             d_U_gL=zeros(3,1);
             d_U_gR=zeros(3,1);
             d_U_sL=zeros(3,1);
             d_U_sR=zeros(3,1);
             d_phi_sL=0.0;
             d_phi_sR=0.0;
         elseif i==N+1
             [lo_gL u_gL p_gL phi_gL lo_sL u_sL p_sL phi_sL]=primitive_comp(U(:,N));
             [lo_gR u_gR p_gR phi_gR lo_sR u_sR p_sR phi_sR]=primitive_comp(U(:,N));
             d_U_gL=zeros(3,1);
             d_U_gR=zeros(3,1);
             d_U_sL=zeros(3,1);
             d_U_sR=zeros(3,1);
             d_phi_sL=0.0;
             d_phi_sR=0.0;
         else
             [lo_gL u_gL p_gL phi_gL lo_sL u_sL p_sL phi_sL]=primitive_comp(U(:,i-1)+0.5*d_x*d_U(:,i-1));
             [lo_gR u_gR p_gR phi_gR lo_sR u_sR p_sR phi_sR]=primitive_comp(U(:,i)-0.5*d_x*d_U(:,i));
             d_U_gL=d_U(1:3,i-1);
             d_U_gR=d_U(1:3,i);
             d_U_sL=d_U(4:6,i-1);
             d_U_sR=d_U(4:6,i);
             d_phi_sL=d_U(7,i-1);
             d_phi_sR=d_U(7,i);
         end
       [FL(:,i),FR(:,i),U_int(:,i)]=GRP_solver_HLLC(lo_gL,lo_gR,p_gL,p_gR,u_gL,u_gR,lo_sL,lo_sR,p_sL,p_sR,u_sL,u_sR,phi_sL,phi_sR,d_U_gL,d_U_gR,d_U_sL,d_U_sR,d_phi_sL,d_phi_sR,d_t/d_x);
    end
    %compute U in next step
    for i=1:N
        [lo_g_mid u_g_mid p_g_mid phi_g_mid lo_s_mid u_s_mid p_s_mid phi_s_mid]=primitive_comp(U(:,i)-0.5*d_t*C_U*d_U(:,i));
        U(:,i)=U(:,i)+d_t/d_x*(FR(:,i)-FL(:,i+1))+[0;-p_g_mid;-p_g_mid*u_s_mid;0;p_g_mid;p_g_mid*u_s_mid;-u_s_mid]*d_U(7,i);
        [lo_g(i) u_g(i) p_g(i) phi_g(i) lo_s(i) u_s(i) p_s(i) phi_s(i)]=primitive_comp(U(:,i));
    end
    Time=Time+d_t
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
load ../test/test_new1.exact;
for i=1:N
     W_exact(i,:) = test_new1(ceil(i/(N/300)),:);
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
