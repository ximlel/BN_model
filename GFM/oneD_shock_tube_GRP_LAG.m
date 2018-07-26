clear;
clc;
%1D shock_tube by GRP for modified ghost fluid method
%state constant
global gama_s gama_g;
global ep;
ep=1e-9;
x_min=0;
x_max=2;
N=200;
d_x=(x_max-x_min)/N;
CFL=0.5;
Alpha=1.0;
%Alpha=0;
%state value
Time=0;
Tend=0.1;
U=zeros(3,N);
F=zeros(3,N+1);
u_MID=zeros(3,N+1);
W_int_L=zeros(3,N);
W_int_R=zeros(3,N);
W_int_tmp=zeros(3);
dlo  =zeros(1,N);
du   =zeros(1,N);
dp   =zeros(1,N);
gama =zeros(1,N);
dxS  =1e10*ones(1,N);
load ./data/test55.mat;
EXACT_LOCAT='./data/exact55.mat';
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
    E(i)=p(i)/(gama(i)-1)+0.5*lo(i)*u(i)^2;
    U(:,i)=[lo(i);lo(i)*u(i);E(i)]*d_x;
end
for i=1:(N+1)
    x_int(i)=x_min+(i-1)*d_x; 
end
%Godunov's Method
count=1;
while Time<Tend && isreal(Time)
    %reconstruction (minmod limiter)
    for i=2:(N-1)
        dlo(i)=minmod(Alpha,(lo(i)-lo(i-1))/(x(i)-x(i-1)),(W_int_R(1,i)-W_int_L(1,i))/(x_int(i+1)-x_int(i)),(lo(i+1)-lo(i))/(x(i+1)-x(i)));
        du(i) =minmod(Alpha,(u(i) -u(i-1) )/(x(i)-x(i-1)),(W_int_R(2,i)-W_int_L(2,i))/(x_int(i+1)-x_int(i)),(u(i+1) -u(i) )/(x(i+1)-x(i)));
        dp(i) =minmod(Alpha,(p(i) -p(i-1) )/(x(i)-x(i-1)),(W_int_R(3,i)-W_int_L(3,i))/(x_int(i+1)-x_int(i)),(p(i+1) -p(i) )/(x(i+1)-x(i)));
%         dlo(i)=minmod(Alpha,(lo(i)-lo(i-1))/(x(i)-x(i-1)),(lo(i+1)-lo(i-1))/(x(i+1)-x(i-1)),(lo(i+1)-lo(i))/(x(i+1)-x(i)));
%         du(i) =minmod(Alpha,(u(i) -u(i-1) )/(x(i)-x(i-1)),(u(i+1) -u(i-1) )/(x(i+1)-x(i-1)),(u(i+1) -u(i) )/(x(i+1)-x(i)));
%         dp(i) =minmod(Alpha,(p(i) -p(i-1) )/(x(i)-x(i-1)),(p(i+1) -p(i-1) )/(x(i+1)-x(i-1)),(p(i+1) -p(i) )/(x(i+1)-x(i)));
%         dlo(i)=minmod(Alpha,(lo(i)-lo(i-1))/(x(i)-x_int(i)),(W_int_R(1,i)-W_int_L(1,i))/(x_int(i+1)-x_int(i)),(lo(i+1)-lo(i))/(x_int(i+1)-x(i)));
%         du(i) =minmod(Alpha,(u(i) -u(i-1) )/(x(i)-x_int(i)),(W_int_R(2,i)-W_int_L(2,i))/(x_int(i+1)-x_int(i)),(u(i+1) -u(i) )/(x_int(i+1)-x(i)));
%         dp(i) =minmod(Alpha,(p(i) -p(i-1) )/(x(i)-x_int(i)),(W_int_R(3,i)-W_int_L(3,i))/(x_int(i+1)-x_int(i)),(p(i+1) -p(i) )/(x_int(i+1)-x(i)));
    end
    %CFL condition
    for i=2:(N-1)
        dx_min=min([x_int(i+2)-x_int(i+1),x_int(i+1)-x_int(i),x_int(i)-x_int(i-1)]);
        dxS(i)=dx_min/(abs(u(i))+sqrt(gama(i)*p(i)/lo(i)));
    end
    d_t=CFL*min(dxS);
    if Time+d_t >= Tend
        d_t = Tend-Time+eps;
    end
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
         if i==1
            [F(:,1),u_MID(1),W_int_tmp,W_int_L(:,1)]=GRP_solver_LAG(lo(1),lo(1),0,0,u(1),u(1),0,0,p(1),p(1),0,0,gama(N),gama(N),d_t);
         elseif i==N+1
            [F(:,N+1),u_MID(N+1),W_int_R(:,N),W_int_tmp]=GRP_solver_LAG(lo(N),lo(N),0,0,u(N),u(N),0,0,p(N),p(N),0,0,gama(N),gama(N),d_t);
         else
            d_xL=(x_int(i)-x_int(i-1));
            d_xR=(x_int(i+1)-x_int(i));
            [F(:,i),u_MID(i),W_int_R(:,i-1),W_int_L(:,i)]=GRP_solver_LAG(lo(i-1)+0.5*d_xL*dlo(i-1),lo(i)-0.5*d_xR*dlo(i),dlo(i-1),dlo(i),u(i-1)+0.5*d_xL*du(i-1),u(i)-0.5*d_xR*du(i),du(i-1),du(i),p(i-1)+0.5*d_xL*dp(i-1),p(i)-0.5*d_xR*dp(i),dp(i-1),dp(i),gama(i-1),gama(i),d_t);
         end
    end
    %compute U in next step
    for i=1:(N+1)
        x_int(i)=x_int(i)+u_MID(i)*d_t;
    end
    for i=1:N
        x(i)=0.5*(x_int(i)+x_int(i+1)); 
    end
    for i=1:N
        U(:,i)=U(:,i)+d_t*(F(:,i)-F(:,i+1));
        [lo(i),u(i),p(i)]=primitive_comp(U(:,i)/(x_int(i+1)-x_int(i)),gama(i));
    end
    count=count+1;
    Time=Time+d_t
if Time > 100000*d_t
    break;
end
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
W_exact(:,1)=log(W_exact(:,1));
W_exact(:,3)=log(W_exact(:,3));
lo=log(lo);
p=log(p);

%plot
col = '.b';
figure(1);
subplot(2,2,1);
hold on
plot(x_min+0.5*d_x:d_x:x_max-0.5*d_x,W_exact(:,1),'k','LineWidth',1.0);
plot(x,lo,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Density','FontWeight','bold');
subplot(2,2,2);
hold on
plot(x_min+0.5*d_x:d_x:x_max-0.5*d_x,W_exact(:,2),'k','LineWidth',1.0);
plot(x,u,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Velocity','FontWeight','bold');
subplot(2,2,3);
hold on
plot(x_min+0.5*d_x:d_x:x_max-0.5*d_x,W_exact(:,3),'k','LineWidth',1.0);
plot(x,p,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Pressure','FontWeight','bold');
subplot(2,2,4);
hold on
% plot(x_min+0.5*d_x:d_x:x_max-0.5*d_x,W_exact(:,4),'k','LineWidth',1.0);
plot(x,gama,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Gamma','FontWeight','bold');
