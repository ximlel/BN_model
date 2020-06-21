warning off 
% for IJ=1:100
% save myfile IJ;
clc
clear
% load myfile
% IJ
%1D shock_tube by 1-order staggered Schemes for BN model
%state constant
global gama_s gama_g p0;
gama_s=1.4;
gama_g=1.4;
p0=0;
global ep;
ep=1e-6;
x_min=0;
x_max=1;
N=300*1;
d_x=(x_max-x_min)/N;
x0=0.5;
CFL=0.1;
%state value
Time=0;
Tend=0.1;
%Tend=0.15;
Alpha=zeros(1,N+1);
U=zeros(6,N);
F=zeros(6,N+1);
U_lo_sL=zeros(1,N);
U_lo_sR=zeros(1,N);
W_int_s=zeros(4,N+1);
W_int_g=zeros(4,N+1);
%initial condition
load ../test/test_dele0.mat;
% u_sR_0 = (8.7+0.1*(IJ/100));

phi_gL_0=1.0-phi_sL_0;
phi_gR_0=1.0-phi_sR_0;
E_gL_0=p_gL_0/(gama_g-1)+0.5*lo_gL_0*u_gL_0^2;
E_sL_0=(p_sL_0+gama_s*p0)/(gama_s-1)+0.5*lo_sL_0*u_sL_0^2;
U_L_0=[phi_gL_0*lo_gL_0;phi_gL_0*lo_gL_0*u_gL_0;phi_gL_0*E_gL_0;phi_sL_0*lo_sL_0;phi_sL_0*lo_sL_0*u_sL_0;phi_sL_0*E_sL_0];
E_gR_0=p_gR_0/(gama_g-1)+0.5*lo_gR_0*u_gR_0^2;
E_sR_0=(p_sR_0+gama_s*p0)/(gama_s-1)+0.5*lo_sR_0*u_sR_0^2;
U_R_0=[phi_gR_0*lo_gR_0;phi_gR_0*lo_gR_0*u_gR_0;phi_gR_0*E_gR_0;phi_sR_0*lo_sR_0;phi_sR_0*lo_sR_0*u_sR_0;phi_sR_0*E_sR_0];

%[lo_gL_1,u_gL_1,p_gL_1,lo_sL_1,u_sL_1,p_sL_1,lo_gR_1,u_gR_1,p_gR_1,lo_sR_1,u_sR_1,p_sR_1]=Riemann_ave(lo_gL_0,u_gL_0,p_gL_0,lo_sL_0,u_sL_0,p_sL_0,phi_sL_0,lo_gR_0,u_gR_0,p_gR_0,lo_sR_0,u_sR_0,p_sR_0,phi_sR_0);
%[lo_gL_1,u_gL_1,p_gL_1,lo_sL_1,u_sL_1,p_sL_1,lo_gR_1,u_gR_1,p_gR_1,lo_sR_1,u_sR_1,p_sR_1]=primitive_ave(0.5*(U_L_0+U_R_0),phi_sL_0,phi_sR_0);

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
       U(:,i) =0.5*(U_L_0+U_R_0);
%         U(:,i) =0.5*(U_L_1+U_R_1);
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
        phi_sL=Alpha(i);
        phi_sR=Alpha(i+1);
        if abs(phi_sL-phi_sR)>ep
            phi_gL=1.0-phi_sL;
            phi_gR=1.0-phi_sR;
            E_gL=p_gL(i)/(gama_g-1)+0.5*lo_gL(i)*u_gL(i)^2;
            E_sL=(p_sL(i)+gama_s*p0)/(gama_s-1)+0.5*lo_sL(i)*u_sL(i)^2;
            U(:,i)=       0.5*[phi_gL*lo_gL(i);phi_gL*lo_gL(i)*u_gL(i);phi_gL*E_gL;phi_sL*lo_sL(i);phi_sL*lo_sL(i)*u_sL(i);phi_sL*E_sL];
            E_gR=p_gR(i)/(gama_g-1)+0.5*lo_gR(i)*u_gR(i)^2;
            E_sR=(p_sR(i)+gama_s*p0)/(gama_s-1)+0.5*lo_sR(i)*u_sR(i)^2;
            U(:,i)=U(:,i)+0.5*[phi_gR*lo_gR(i);phi_gR*lo_gR(i)*u_gR(i);phi_gR*E_gR;phi_sR*lo_sR(i);phi_sR*lo_sR(i)*u_sR(i);phi_sR*E_sR];
        end
         
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
            %F(1:3,1)=Riemann_solver_Exact(lo_gL(1),lo_gL(1),p_gL(1),p_gL(1),u_gL(1),u_gL(1),1-Alpha(1),gama_g,0.0);
            %F(4:6,1)=Riemann_solver_Exact(lo_sL(1),lo_sL(1),p_sL(1),p_sL(1),u_sL(1),u_sL(1),Alpha(1),gama_s,0.0);
            [F(1:3,1),W_int_g(:,1)]=GRP_solver(lo_gR(1),lo_gL(1),0.0,0.0,u_gR(1),u_gL(1),0.0,0.0,p_gR(1),p_gL(1),0.0,0.0,1-Alpha(1),1-Alpha(1),0.0,0.0,gama_g,d_t);
            [F(4:6,1),W_int_s(:,1)]=GRP_solver(lo_sR(1),lo_sL(1),0.0,0.0,u_sR(1),u_sL(1),0.0,0.0,p_sR(1),p_sL(1),0.0,0.0,  Alpha(1),  Alpha(1),0.0,0.0,gama_s,d_t);
        elseif i==N+1
            %F(1:3,N+1)=Riemann_solver_Exact(lo_gR(N),lo_gR(N),p_gR(N),p_gR(N),u_gR(N),u_gR(N),1-Alpha(N+1),gama_g,0.0);
            %F(4:6,N+1)=Riemann_solver_Exact(lo_sR(N),lo_sR(N),p_sR(N),p_sR(N),u_sR(N),u_sR(N),Alpha(N+1),gama_s,0.0);
            [F(1:3,N+1),W_int_g(:,N+1)]=GRP_solver(lo_gR(N),lo_gL(N),0.0,0.0,u_gR(N),u_gL(N),0.0,0.0,p_gR(N),p_gL(N),0.0,0.0,1-Alpha(N+1),1-Alpha(N+1),0.0,0.0,gama_g,d_t);
            [F(4:6,N+1),W_int_s(:,N+1)]=GRP_solver(lo_sR(N),lo_sL(N),0.0,0.0,u_sR(N),u_sL(N),0.0,0.0,p_sR(N),p_sL(N),0.0,0.0,  Alpha(N+1),  Alpha(N+1),0.0,0.0,gama_s,d_t);
        else
            %F(1:3,i)=Riemann_solver_Exact(lo_gR(i-1),lo_gL(i),p_gR(i-1),p_gL(i),u_gR(i-1),u_gL(i),1-Alpha(i),gama_g,0.0);
            %F(4:6,i)=Riemann_solver_Exact(lo_sR(i-1),lo_sL(i),p_sR(i-1),p_sL(i),u_sR(i-1),u_sL(i),Alpha(i),gama_s,0.0);
            [F(1:3,i),W_int_g(:,i)]=GRP_solver(lo_gR(i-1),lo_gL(i),0.0,0.0,u_gR(i-1),u_gL(i),0.0,0.0,p_gR(i-1),p_gL(i),0.0,0.0,1-Alpha(i),1-Alpha(i),0.0,0.0,gama_g,d_t);
            [F(4:6,i),W_int_s(:,i)]=GRP_solver(lo_sR(i-1),lo_sL(i),0.0,0.0,u_sR(i-1),u_sL(i),0.0,0.0,p_sR(i-1),p_sL(i),0.0,0.0,  Alpha(i),  Alpha(i),0.0,0.0,gama_s,d_t);
        end
    end
    for i=1:N
        if i<N
            rho_s(i) =0.5*(lo_sR(i)+lo_sL(i+1));
            arho_s(i)=0.5*(lo_sR(i)+lo_sL(i+1))*Alpha(i+1);           
        end
        F_rho_s(i) =lo_sL(i)*u_sL(i);
        if u_sL(i)>0
            F_arho_s(i)=  Alpha(i)*lo_sL(i)*u_sL(i);
        else
            F_arho_s(i)=Alpha(i+1)*lo_sR(i)*u_sR(i);
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
          S_tmp=Alpha(i+1)*p_sR(i)-Alpha(i)*p_sL(i);
          if (S_tmp/(Alpha(i+1)-Alpha(i))>max(p_gL(i),p_gR(i)))
              S=max(p_gL(i),p_gR(i))*(Alpha(i+1)-Alpha(i));
          elseif (S_tmp/(Alpha(i+1)-Alpha(i))<min(p_gL(i),p_gR(i)))
              S=min(p_gL(i),p_gR(i))*(Alpha(i+1)-Alpha(i));
          else
              S=S_tmp;
          end
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
    Time=Time+d_t;
%     if Time > 10*d_t
%         break;
%     end
end
% end
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
% load ../test/test1.exact;
% for i=1:N
%      W_exact(i,:) = test1(ceil(i/(N/300)),:);
% end
%plot
%col = '+k';
%col = 'or';
%col = '*m';
col = 'xr';
POS = [50 50 1200 800];
h1=figure(1);
set(h1,'position',POS);
subplot(2,2,1);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,3),'b','LineWidth',0.4);
plot(x,lo_s,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Density-solid','FontWeight','bold');
% ylim([1.96 2.01])
title('Solid density')
subplot(2,2,2);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,4),'b','LineWidth',0.4);
plot(x,u_s,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Velocity-solid','FontWeight','bold');
% ylim([0.298 0.31])
title('Solid velocity')
subplot(2,2,3);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,5),'b','LineWidth',0.4);
plot(x,p_s,col,'MarkerSize',4);
% ylim([4 14])
% xlabel('Position','FontWeight','bold');
% ylabel('Pressure-solid','FontWeight','bold');
title('Solid pressure')
subplot(2,2,4);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,2),'b','LineWidth',0.4);
plot(x,phi_s,col,'MarkerSize',4);
% ylim([0.3 0.9])
% xlabel('Position','FontWeight','bold');
% ylabel('Porosity-solid','FontWeight','bold');
title('Solid volume fraction')
h2=figure(2);
set(h2,'position',POS);
subplot(2,2,1);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,6),'b','LineWidth',0.4);
plot(x,lo_g,col,'MarkerSize',4);
% ylim([0 1])
% xlabel('Position','FontWeight','bold');
% ylabel('Density-gas','FontWeight','bold');
title('Gas density')
subplot(2,2,2);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,7),'b','LineWidth',0.4);
plot(x,u_g,col,'MarkerSize',4);
% ylim([2 3])
% xlabel('Position','FontWeight','bold');
% ylabel('Velocity-gas','FontWeight','bold');
title('Gas velocity')
subplot(2,2,3);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,8),'b','LineWidth',0.4);
plot(x,p_g,col,'MarkerSize',4);
% ylim([0 1])
% xlabel('Position','FontWeight','bold');
% ylabel('Pressure-gas','FontWeight','bold');
title('Gas pressure')
subplot(2,2,4);
hold on
% plot(x_min:d_x:x_max-d_x,W_exact(:,8)./W_exact(:,6).^gama_g,'b','LineWidth',0.4);
% plot(x,eta,col,'MarkerSize',4);
% % xlabel('Position','FontWeight','bold');
% % ylabel('Entropy-gas','FontWeight','bold');
% ylim([min(eta)-0.00001 max(eta)+0.00001])
% title('Gas entropy')
plot(x_min:d_x:x_max-d_x,1.0-W_exact(:,2),'b','LineWidth',0.4);
plot(x,1.0-phi_s,col,'MarkerSize',4);
% ylim([0.1 0.7])
title('Gas volume fraction')
h3=figure(3)
set(h3,'position',POS);
subplot(3,1,1);
hold on
plot(x_min:d_x:x_max-d_x,(1-W_exact(:,2)).*W_exact(:,6).*(W_exact(:,7)-W_exact(:,4)),'b','LineWidth',0.4);
plot(x,Q_inv,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Riemann_inv-Q','FontWeight','bold');
% ylim([min(Q_inv)-0.00001 max(Q_inv)+0.00001])
title('Riemann invariants Q')
subplot(3,1,2);
hold on
plot(x_min:d_x:x_max-d_x,(1-W_exact(:,2)).*W_exact(:,6).*(W_exact(:,7)-W_exact(:,4)).^2+(1-W_exact(:,2)).*W_exact(:,8)+W_exact(:,2).*W_exact(:,5),'b','LineWidth',0.4);
plot(x,P_inv,col,'MarkerSize',4);
% xlabel('Position','FontWeight','bold');
% ylabel('Riemann_inv-P','FontWeight','bold');
% ylim([min(P_inv)-0.00001 max(P_inv)+0.00001])
title('Riemann invariants P')
subplot(3,1,3);
hold on
plot(x_min:d_x:x_max-d_x,0.5*(W_exact(:,7)-W_exact(:,4)).^2+gama_g/(gama_g-1)*W_exact(:,8)./W_exact(:,6),'b','LineWidth',0.4);
plot(x,H_inv,col,'MarkerSize',4);
% ylim([min(H_inv)-0.00001 max(H_inv)+0.00001])
% xlabel('Position','FontWeight','bold');
% ylabel('Riemann_inv-H','FontWeight','bold');
title('Riemann invariants H')
