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
CFL=0.7;
Alpha=1.9;
%Alpha=0;
%state value
Time=0;
Tend=0.1;
U=zeros(3,N);
F=zeros(3,N+1);
W_int_L=zeros(4,N);
W_int_R=zeros(4,N);
dlo  =zeros(1,N);
du   =zeros(1,N);
dp   =zeros(1,N);
gama =zeros(1,N);
load ./data/test5.mat;
EXACT_LOCAT='./data/exact5.mat';
%test begin
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if x(i)<=x0s
        lo(i)=lo_L_0;
        u(i) =u_L_0;
        p(i) =p_L_0;
        gama(i)=gama_s;
    else
        if x(i)>=x0
          lo(i)=lo_R_0;
          u(i) =u_R_0;
          p(i) =p_R_0;
          gama(i)=gama_g;
        else
          lo(i)=lo_M_0;
          u(i) =u_R_0;
          p(i) =p_R_0;
          gama(i)=gama_s;
        end
    end
    mass(i)=lo(i)*d_x;
    E(i)=p(i)/(gama(i)-1)+0.5*lo(i)*u(i)^2;
    U(:,i)=[lo(i);lo(i)*u(i);E(i)];
end
%Godunov's Method
count=1;
while Time<Tend && isreal(Time)
    %reconstruction (minmod limiter)
    for i=2:(N-1)
        dlo(i)=minmod(Alpha,(lo_s(i)-lo_s(i-1))/d_x,(W_int_R(1,i)-W_int_L(1,i))/d_x,(lo_s(i+1)-lo_s(i))/d_x);
        du(i) =minmod(Alpha,(u_s(i) -u_s(i-1) )/d_x,(W_int_R(2,i)-W_int_L(2,i))/d_x,(u_s(i+1) -u_s(i) )/d_x);
        dp(i) =minmod(Alpha,(p_s(i) -p_s(i-1) )/d_x,(W_int_R(3,i)-W_int_L(3,i))/d_x,(p_s(i+1) -p_s(i) )/d_x);
    end
    %CFL condition
    for i=1:N
        S(i)=abs(u(i))+sqrt(gama(i)*p(i)/lo(i));
    end
    d_t=CFL*d_x/max(S);
    if Time+d_t >= Tend
        d_t = Tend-Time+eps;
    end
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
         if i==1
            [F(:,1),W_int_Rtmp,W_int_L(:,1)]=GRP_solver(lo(1),lo(1),0,0,u(1),u(1),0,0,p(1),p(1),0,0,gama(N),gama(N),d_t);
         elseif i==N+1
            [F(:,N+1),W_int_R(:,N),W_int_Ltmp]=GRP_solver(lo(N),lo(N),0,0,u(N),u(N),0,0,p(N),p(N),0,0,gama(N),gama(N),d_t);
         else
            [F(:,i),W_int_R(:,i-1),W_int_L(:,i)]=GRP_solver(lo(i-1)+0.5*d_x*dlo(i-1),lo(i)-0.5*d_x*dlo(i),dlo(i-1),dlo(i),u(i-1)+0.5*d_x*du(i-1),u(i)-0.5*d_x*du(i),du(i-1),du(i),p(i-1)+0.5*d_x*dp(i-1),p(i)-0.5*d_x*dp(i),dp(i-1),dp(i),gama(i-1),gama(i),d_t);
         end
    end
    %compute U in next step
    x0=x0+u_s(J+1)*d_t;
    for i=1:(J+1)
        U(:,i)=U(:,i)+d_t/d_x*(F(:,i)-F(:,i+1));
        [lo_s(i),u_s(i),p_s(i)]=primitive_comp(U_s(:,i),gama(i));
    end
    count=count+1;
    Time=Time+d_t
% if Time > 1*d_t
%     break;
% end
end
W_exact = zeros(N,4);
W_exact(:,1)=lo';
W_exact(:,2)=u';
W_exact(:,3)=p';
W_exact(:,4)=gama';
load(EXACT_LOCAT);
for i=1:N
     W_exact(i,1) = lo_ex(ceil(i/(N/200)));
     W_exact(i,2) = u_ex(ceil(i/(N/200)));
     W_exact(i,3) = p_ex(ceil(i/(N/200)));
end
% W_exact(:,1)=log(W_exact(:,1));
% W_exact(:,3)=log(W_exact(:,3));
% lo=log(lo);
% p=log(p);

%plot
col = '.b';
figure(1);
subplot(2,2,1);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,1),'k','LineWidth',1.0);
plot(x,lo,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Density','FontWeight','bold');
subplot(2,2,2);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,2),'k','LineWidth',1.0);
plot(x,u,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Velocity','FontWeight','bold');
subplot(2,2,3);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,3),'k','LineWidth',1.0);
plot(x,p,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Pressure','FontWeight','bold');
subplot(2,2,4);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,4),'k','LineWidth',1.0);
plot(x,gamma,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Gamma','FontWeight','bold');
