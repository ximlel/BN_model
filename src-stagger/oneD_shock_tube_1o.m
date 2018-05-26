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
CFL=0.2;
%state value
Time=0;
Tend=0.1;
%Tend=0.15;
Alpha=zeros(1,N+1);
U=zeros(6,N);
F=zeros(6,N+1);
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
%test begin
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if i<round(N*x0/(x_max-x_min))
        U(:,i) =U_L_0;
        U_lo_sL(i) =phi_sL_0*lo_sL_0;
        U_lo_sR(i) =phi_sL_0*lo_sL_0;
        Alpha(i) =phi_sL_0;
    elseif i>round(N*x0/(x_max-x_min))
        U(:,i) =U_R_0;
        U_lo_sL(i) =phi_sR_0*lo_sR_0;
        U_lo_sR(i) =phi_sR_0*lo_sR_0;
        Alpha(i+1) =phi_sR_0;
    else
        U(:,i) =0.5*(U_L_0+U_R_0);
        U_lo_sL(i) =phi_sL_0*lo_sL_0;
        U_lo_sR(i) =phi_sR_0*lo_sR_0;
        Alpha(i) =phi_sL_0;
        Alpha(i+1) =phi_sR_0;
    end
end
%Godunov's Method
while Time<Tend && isreal(Time)
    %CFL condition
    for i=1:N
      if i==1
          x_delta_L=0.5*d_x;
          x_delta_R=0.5*(x(2)-x(1));          
      elseif i==N
          x_delta_L=0.5*(x(N)-x(N-1));
          x_delta_R=0.5*d_x;
      else
          x_delta_L=0.5*(x(i)-x(i-1));
          x_delta_R=0.5*(x(i+1)-x(i));
      end
        [lo_gL(i),u_gL(i),p_gL(i),lo_sL(i),u_sL(i),p_sL(i),lo_gR(i),u_gR(i),p_gR(i),lo_sR(i),u_sR(i),p_sR(i)]=primitive_comp(U(:,i),U_lo_sL(i),U_lo_sR(i),Alpha(i),Alpha(i+1),x_delta_L/(x_delta_L+x_delta_R),x_delta_R/(x_delta_L+x_delta_R));
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
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
        if i==1
            F(1:3,1)=Riemann_solver_Exact(lo_gL(1),lo_gL(1),p_gL(1),p_gL(1),u_gL(1),u_gL(1),1-Alpha(1),gama_g,u_sL(1));
            F(4:6,1)=Riemann_solver_Exact(lo_sL(1),lo_sL(1),p_sL(1),p_sL(1),u_sL(1),u_sL(1),Alpha(1),gama_s,u_sL(1));
        elseif i==N+1
            F(1:3,N+1)=Riemann_solver_Exact(lo_gR(N),lo_gR(N),p_gR(N),p_gR(N),u_gR(N),u_gR(N),1-Alpha(N+1),gama_g,u_sL(N));
            F(4:6,N+1)=Riemann_solver_Exact(lo_sR(N),lo_sR(N),p_sR(N),p_sR(N),u_sR(N),u_sR(N),Alpha(N+1),gama_s,u_sL(N));
        else
            F(1:3,i)=Riemann_solver_Exact(lo_gR(i-1),lo_gL(i),p_gR(i-1),p_gL(i),u_gR(i-1),u_gL(i),1-Alpha(i),gama_g,0.5*(u_sL(i-1)+u_sL(i)));
            F(4:6,i)=Riemann_solver_Exact(lo_sR(i-1),lo_sL(i),p_sR(i-1),p_sL(i),u_sR(i-1),u_sL(i),Alpha(i),gama_s,0.5*(u_sL(i-1)+u_sL(i)));
        end
    end
    %compute U in next step
    for i=1:N
      if abs(Alpha(i+1)-Alpha(i))<ep
          S=0.5*(p_gL(i)+p_gR(i))*(Alpha(i+1)-Alpha(i));
      else
          S_tmp=Alpha(i+1)*p_sR(i)-Alpha(i)*p_sL(i);
          if (S_tmp/(Alpha(i+1)-Alpha(i))>max(p_gL(i),p_gR(i)))
              S=max(p_gL(i),p_gR(i))*(Alpha(i+1)-Alpha(i));
          elseif (S_tmp/(Alpha(i+1)-Alpha(i))<min(p_gL(i),p_gR(i)))
              S=min(p_gL(i),p_gR(i))*(Alpha(i+1)-Alpha(i));
          else
              S=S_tmp;
          end
      end
      if i==1
          x_delta=0.5*(x(2)-x(1)+d_x);
      elseif i==N
          x_delta=0.5*(x(N)-x(N-1)+d_x);
      else
          x_delta=0.5*(x(i+1)-x(i-1));
      end
        U(:,i)=U(:,i)*x_delta+d_t*(F(:,i)-F(:,i+1))+d_t*[0;-S;-S*u_sL(i);0;S;S*u_sL(i)];
        U_lo_sL(i)=U_lo_sL(i)*x_delta+d_t*(F(4,i));
        U_lo_sR(i)=U_lo_sR(i)*x_delta+d_t*(-F(4,i+1));

    end
    for i=1:N
        x(i)=x(i)+u_sL(i)*d_t;
    end
    for i=1:N
      if i==1
          x_delta=0.5*(x(2)-x(1)+d_x);
      elseif i==N
          x_delta=0.5*(x(N)-x(N-1)+d_x);
      else
          x_delta=0.5*(x(i+1)-x(i-1));
      end
        U(:,i)=U(:,i)/x_delta;
        U_lo_sL(i)=U_lo_sL(i)/x_delta;
        U_lo_sR(i)=U_lo_sR(i)/x_delta;
    end
    Time=Time+d_t
% if Time > 0.5*d_t
%    break;
% end
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
col = ':.b';
figure(1);
subplot(2,2,1);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,3),'k','LineWidth',1.0);
plot(x,lo_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Density-solid','FontWeight','bold');
ylim([min(lo_s)-0.00001 max(lo_s)+0.00001])
subplot(2,2,2);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,4),'k','LineWidth',1.0);
plot(x,u_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Velocity-solid','FontWeight','bold');
ylim([min(u_s)-0.00001 max(u_s)+0.00001])
subplot(2,2,3);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,5),'k','LineWidth',1.0);
plot(x,p_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Pressure-solid','FontWeight','bold');
subplot(2,2,4);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,2),'k','LineWidth',1.0);
plot(x,phi_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Porosity-solid','FontWeight','bold');
figure(2);
subplot(2,2,1);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,6),'k','LineWidth',1.0);
plot(x,lo_g,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Density-gas','FontWeight','bold');
subplot(2,2,2);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,7),'k','LineWidth',1.0);
plot(x,u_g,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Velocity-gas','FontWeight','bold');
subplot(2,2,3);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,8),'k','LineWidth',1.0);
plot(x,p_g,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Pressure-gas','FontWeight','bold');
subplot(2,2,4);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,8)./W_exact(:,6).^gama_g,'k','LineWidth',1.0);
plot(x,eta,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Entropy-gas','FontWeight','bold');
ylim([min(eta)-0.00001 max(eta)+0.00001])
figure(3)
subplot(3,1,1);
hold on
plot(x_min:d_x:x_max-d_x,(1-W_exact(:,2)).*W_exact(:,6).*(W_exact(:,7)-W_exact(:,4)),'k','LineWidth',1.0);
plot(x,Q_inv,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Riemann_inv-Q','FontWeight','bold');
ylim([min(Q_inv)-0.00001 max(Q_inv)+0.00001])
subplot(3,1,2);
hold on
plot(x_min:d_x:x_max-d_x,(1-W_exact(:,2)).*W_exact(:,6).*(W_exact(:,7)-W_exact(:,4)).^2+(1-W_exact(:,2)).*W_exact(:,8)+W_exact(:,2).*W_exact(:,5),'k','LineWidth',1.0);
plot(x,P_inv,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Riemann_inv-P','FontWeight','bold');
ylim([min(P_inv)-0.00001 max(P_inv)+0.00001])
subplot(3,1,3);
hold on
plot(x_min:d_x:x_max-d_x,0.5*(W_exact(:,7)-W_exact(:,4)).^2+gama_g/(gama_g-1)*W_exact(:,8)./W_exact(:,6),'k','LineWidth',1.0);
plot(x,H_inv,col,'LineWidth',1.0);
ylim([min(H_inv)-0.00001 max(H_inv)+0.00001])
xlabel('Position','FontWeight','bold');
ylabel('Riemann_inv-H','FontWeight','bold');
