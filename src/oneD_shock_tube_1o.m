clear;
clc;
%1D shock_tube by HLLC Schemes for BN model
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
%state value
Time=0;
Tend=0.1;
%Tend=0.15;
U=zeros(7,N);
FL=zeros(7,N+1);
FR=zeros(7,N+1);
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
load ../test/test_new1_pi.mat;
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
         if i==1
             [phi_s_out,FL(:,1),FR(:,1)]=Riemann_solver_HLLC(lo_g(1),lo_g(1),p_g(1),p_g(1),u_g(1),u_g(1),lo_s(1),lo_s(1),p_s(1),p_s(1),u_s(1),u_s(1),phi_s(1),phi_s(1),d_t/d_x);
         elseif i==N+1
             [phi_s_out,FL(:,N+1),         FR(:,N+1)]=Riemann_solver_HLLC(lo_g(N),lo_g(N),p_g(N),p_g(N),u_g(N),u_g(N),lo_s(N),lo_s(N),p_s(N),p_s(N),u_s(N),u_s(N),phi_s(N),phi_s(N),d_t/d_x);
         else
             [phi_s_out,FL(:,i),FR(:,i)]=Riemann_solver_HLLC(lo_g(i-1),lo_g(i),p_g(i-1),p_g(i),u_g(i-1),u_g(i),lo_s(i-1),lo_s(i),p_s(i-1),p_s(i),u_s(i-1),u_s(i),phi_s(i-1),phi_s(i),d_t/d_x);
         end
         if phi_s_out > 0.0 && i<N+1
             U(7,i) = phi_s_out;
         elseif phi_s_out <= 0.0 && i>1
             U(7,i-1) = -phi_s_out;
         end
    end
    %compute U in next step
    for i=1:N
        U(1:6,i)=U(1:6,i)+d_t/d_x*(FR(1:6,i)-FL(1:6,i+1));
%        U(:,i)=U(:,i)+d_t/d_x*(FR(:,i)-FL(:,i+1));
        [lo_g(i),u_g(i),p_g(i),phi_g(i),lo_s(i),u_s(i),p_s(i),phi_s(i)]=primitive_comp(U(:,i));
    end
    Time=Time+d_t
% if Time > d_t
%    break;
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
load ../test/test_new1_pi.exact;
for i=1:N
     W_exact(i,:) = test_new1_pi(ceil(i/(N/300)),:);
end

%plot
col = ':.m';
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
