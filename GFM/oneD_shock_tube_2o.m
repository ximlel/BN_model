clear;
clc;
%1D shock_tube by GRP for modified ghost fluid method
%state constant
global gama_s gama_g;
global ep;
ep=1e-9;
x_min=0;
x_max=1;
N=200;
d_x=(x_max-x_min)/N;
CFL=0.9;
Alpha=1.0;
%state value
Time=0;
Tend=0.1;
%Tend=0.15;
U_s  =zeros(3,N);
F_s  =zeros(3,N+1);
U_g  =zeros(3,N);
F_g  =zeros(3,N+1);
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
        lo_s(i)=lo_L_0;
        u_s(i) =u_L_0;
        p_s(i) =p_L_0;
    else
        lo_g(i)=lo_R_0;
        u_g(i) =u_R_0;
        p_g(i) =p_R_0;
    end
    phi(i)=x(i)-x_0;
    if phi(i)<0.0;
      E_s(i)=p_s(i)/(gama_s-1)+0.5*lo_s(i)*u_S(i)^2;
      U_s(:,i)=[lo_s(i);lo_s(i)*u_s(i);E_s(i)];
    else
      E_g(i)=p_g(i)/(gama_g-1)+0.5*lo_g(i)*u_g(i)^2;
      U_g(:,i)=[lo_g(i);lo_g(i)*u_g(i);E_g(i)];
    end
end
%Godunov's Method
while Time<Tend && isreal(Time)
    for i=1:N-1
        if (phi(i)*phi(i+1))<0.0
            J=i;
            break;
        end
    end
     if phi(J)<0.0
        [p_g(J),u_g(J),lo_g(J),p_s(J+1),u_s(J+1),lo_s(J+1)]=ghost_cal(lo_s(J-1),u_s(J-1),p_s(J-1),gama_s,lo_g(J+2),u_g(J+2),p_g(J+2),gama_g)
        p_g(J-1)=p_s(J-1);
        u_g(J-1)=u_s(J-1);
        lo_g(J-1)=(p_g(J-1)/p_g(J))^(1/gama_g)*lo_g(J);
        p_g(J-2) =p_s(J-2);
        u_g(J-2) =u_s(J-2);
        lo_g(J-2)=(p_g(J-2)/p_g(J))^(1/gama_g)*lo_g(J);
        p_s(J+2) =p_g(J+2);
        u_s(J+2) =u_g(J+2);
        lo_s(J+2)=(p_s(J+2)/p_s(J+1))^(1/gama_s)*lo_s(J+1);
        p_s(J+3) =p_g(J+3);
        u_s(J+3) =u_g(J+3);
        lo_s(J+3)=(p_s(J+3)/p_s(J+1))^(1/gama_s)*lo_s(J+1);
     end
    for i=(J-2):(J+3)
        E_s(i)=p_s(i)/(gama_s-1)+0.5*lo_s(i)*u_S(i)^2;
        E_g(i)=p_g(i)/(gama_g-1)+0.5*lo_g(i)*u_g(i)^2;
        U_s(:,i)=[lo_s(i);lo_s(i)*u_s(i);E_s(i)];
        U_g(:,i)=[lo_g(i);lo_g(i)*u_g(i);E_g(i)];
    end
    %reconstruction (minmod limiter)
    for i=2:(J+2)
        dlo_s(:,i)=minmod(Alpha*(lo_s(:,i)-lo_s(:,i-1))/d_x,(lo_s(:,i+1)-lo_s(:,i-1))/2.0/d_x,Alpha*(lo_s(:,i+1)-lo_s(:,i))/d_x);
        du_s(:,i) =minmod(Alpha*(u_s(:,i) -u_s(:,i-1) )/d_x,(u_s(:,i+1) -u_s(:,i-1) )/2.0/d_x,Alpha*(u_s(:,i+1) -u_s(:,i) )/d_x);
        dp_s(:,i) =minmod(Alpha*(p_s(:,i) -p_s(:,i-1) )/d_x,(p_s(:,i+1) -p_s(:,i-1) )/2.0/d_x,Alpha*(p_s(:,i+1) -p_s(:,i) )/d_x);
    end
    for i=(J-1):(N-1)
        dlo_g(:,i)=minmod(Alpha*(lo_g(:,i)-lo_g(:,i-1))/d_x,(lo_g(:,i+1)-lo_g(:,i-1))/2.0/d_x,Alpha*(lo_g(:,i+1)-lo_g(:,i))/d_x);
        du_g(:,i) =minmod(Alpha*(u_g(:,i) -u_g(:,i-1) )/d_x,(u_g(:,i+1) -u_g(:,i-1) )/2.0/d_x,Alpha*(u_g(:,i+1) -u_g(:,i) )/d_x);
        dp_g(:,i) =minmod(Alpha*(p_g(:,i) -p_g(:,i-1) )/d_x,(p_g(:,i+1) -p_g(:,i-1) )/2.0/d_x,Alpha*(p_g(:,i+1) -p_g(:,i) )/d_x);
    end
    %CFL condition
    for i=1:J
        S(i)=abs(u_s(i))+sqrt(gama_s*p_s(i)/lo_s(i));
    end
    for i=(J+1):N
        S(i)=abs(u_g(i))+sqrt(gama_g*p_g(i)/lo_g(i));
    end
    d_t=CFL*d_x/max(S);
    if Time+d_t >= Tend
        d_t = Tend-Time+eps;
    end
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
         if i==1
            [F_s(:,1)]=GRP_solver(lo_s(1),lo_s(1),0,0,u_s(1),u_s(1),0,0,p_s(1),p_s(1),0,0,d_t);
         elseif i==N+1
            [F_g(:,N+1)]=GRP_solver(lo_g(N),lo_g(N),0,0,u_g(N),u_g(N),0,0,p_g(N),p_g(1),0,0,d_t);
         else
            if i<=J+2
                [F_s(:,i)]=GRP_solver(lo_s(i-1),lo_s(i),dlo_s(i-1),dlo_s(i),u_s(i-1),u_s(i),du_s(i-1),du_s(i),p_s(i-1),p_s(i),dp_s(i-1),dp_s(i),d_t);
            end
            if i>=J-1
                [F_g(:,i)]=GRP_solver(lo_g(i-1),lo_g(i),dlo_g(i-1),dlo_g(i),u_g(i-1),u_g(i),du_g(i-1),du_g(i),p_g(i-1),p_g(i),dp_g(i-1),dp_g(i),d_t);
            end
         end
    end
    %compute U in next step
    for i=1:(J+1)
        U_s(:,i)=U_s(:,i)+d_t/d_x*(F_s(:,i)-F_s(:,i+1));
        [lo_s(i),u_s(i),p_s(i)]=primitive_comp(U_s(:,i),gama_s);
    end
    for i=J:N
        U_g(:,i)=U_g(:,i)+d_t/d_x*(F_g(:,i)-F_g(:,i+1));
        [lo_g(i),u_g(i),p_g(i)]=primitive_comp(U_g(:,i),gama_g);
    end
    Time=Time+d_t;
% if Time > 5*d_t
%     break;
% end
end
for i=1:J
    lo(i)=lo_s(i);
    u(i) =u_s(i);
    p(i) =p_s(i);
end
for i=(J+1):N
    lo(i)=lo_g(i);
    u(i) =u_g(i);
    p(i) =p_g(i);
end
W_exact = zeros(N,4);
W_exact(:,1)=lo';
W_exact(:,2)=u';
W_exact(:,3)=p';
W_exact(:,4)=phi';
load ../test/test1.exact;
for i=1:N
     W_exact(i,:) = test1(ceil(i/(N/300)),:);
end
%plot
col = '-r';
figure(1);
subplot(2,2,1);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,1),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,lo_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Density-solid','FontWeight','bold');
ylim([min(lo_s)-0.00001 max(lo_s)+0.00001])
subplot(2,2,2);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,2),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,u_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Velocity-solid','FontWeight','bold');
ylim([min(u_s)-0.00001 max(u_s)+0.00001])
subplot(2,2,3);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,3),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,p_s,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Pressure-solid','FontWeight','bold');
subplot(2,2,4);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,4),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,phi,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Porosity-solid','FontWeight','bold');
