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
U_s=zeros(3,N);
F_s=zeros(3,N+1);
U_g=zeros(3,N);
F_g=zeros(3,N+1);
W_int_s=zeros(4,N+1);
W_int_g=zeros(4,N+1);
phi_s=zeros(1,N+1);
phi_g=zeros(1,N+1);
dlo  =zeros(1,N);
du   =zeros(1,N);
dp   =zeros(1,N);
dlo_s=zeros(1,N);
du_s =zeros(1,N);
dp_s =zeros(1,N);
dlo_g=zeros(1,N);
du_g =zeros(1,N);
dp_g =zeros(1,N);
dphi =zeros(1,N);
lo_g =zeros(1,N);
lo_s =zeros(1,N);
u_g =zeros(1,N);
u_s =zeros(1,N);
p_g =zeros(1,N);
p_s =zeros(1,N);
load ./data/test55.mat;
EXACT_LOCAT='./data/exact55.mat';
%test begin
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if x(i)<=x0s
        lo_s(i)=lo_L_0;
        u_s(i) =u_L_0;
        p_s(i) =p_L_0;
    else
        if x(i)>=x0
          lo_g(i)=lo_R_0;
          u_g(i) =u_R_0;
          p_g(i) =p_R_0;
        else
          lo_s(i)=lo_M_0;
          u_s(i) =u_R_0;
          p_s(i) =p_R_0;
        end
    end
    phi(i)=sign(x(i)-x0);
    if phi(i)<0.0
      E_s(i)=p_s(i)/(gama_s-1)+0.5*lo_s(i)*u_s(i)^2;
      U_s(:,i)=[lo_s(i);lo_s(i)*u_s(i);E_s(i)];
    else
      E_g(i)=p_g(i)/(gama_g-1)+0.5*lo_g(i)*u_g(i)^2;
      U_g(:,i)=[lo_g(i);lo_g(i)*u_g(i);E_g(i)];
    end
end
for i=1:N-1
    if (phi(i)*phi(i+1))<0.0
        J=i;
        break;
    end
end
%Godunov's Method
count=1;
while Time<Tend && isreal(Time)
    %reconstruction (minmod limiter)
    for i=2:(J-2)
        dlo_s(i)=minmod(Alpha,(lo_s(i)-lo_s(i-1))/d_x,(lo_s(i+1)-lo_s(i-1))/2.0/d_x,(lo_s(i+1)-lo_s(i))/d_x);
        du_s(i) =minmod(Alpha,(u_s(i) -u_s(i-1) )/d_x,(u_s(i+1) -u_s(i-1) )/2.0/d_x,(u_s(i+1) -u_s(i) )/d_x);
        dp_s(i) =minmod(Alpha,(p_s(i) -p_s(i-1) )/d_x,(p_s(i+1) -p_s(i-1) )/2.0/d_x,(p_s(i+1) -p_s(i) )/d_x);
    end
    for i=(J+3):(N-1)
        dlo_g(i)=minmod(Alpha,(lo_g(i)-lo_g(i-1))/d_x,(lo_g(i+1)-lo_g(i-1))/2.0/d_x,(lo_g(i+1)-lo_g(i))/d_x);
        du_g(i) =minmod(Alpha,(u_g(i) -u_g(i-1) )/d_x,(u_g(i+1) -u_g(i-1) )/2.0/d_x,(u_g(i+1) -u_g(i) )/d_x);
        dp_g(i) =minmod(Alpha,(p_g(i) -p_g(i-1) )/d_x,(p_g(i+1) -p_g(i-1) )/2.0/d_x,(p_g(i+1) -p_g(i) )/d_x);
    end
    i=J-1;
        dlo_s(i)=minmod(Alpha,(lo_s(i)-lo_s(i-1))/d_x,(lo_g(i+3)-lo_s(i-1))/3.0/d_x,(lo_g(i+3)-lo_s(i))/2.0/d_x);
        du_s(i) =minmod(Alpha,(u_s(i) -u_s(i-1) )/d_x,(u_g(i+3) -u_s(i-1) )/3.0/d_x,(u_g(i+3) -u_s(i) )/2.0/d_x);
        dp_s(i) =minmod(Alpha,(p_s(i) -p_s(i-1) )/d_x,(p_g(i+3) -p_s(i-1) )/3.0/d_x,(p_g(i+3) -p_s(i) )/2.0/d_x);
    i=J+2;
        dlo_g(i)=minmod(Alpha,(lo_g(i)-lo_s(i-3))/2.0/d_x,(lo_g(i+1)-lo_s(i-3))/3.0/d_x,(lo_g(i+1)-lo_g(i))/d_x);
        du_g(i) =minmod(Alpha,(u_g(i) -u_s(i-3) )/2.0/d_x,(u_g(i+1) -u_s(i-3) )/3.0/d_x,(u_g(i+1) -u_g(i) )/d_x);
        dp_g(i) =minmod(Alpha,(p_g(i) -p_s(i-3) )/2.0/d_x,(p_g(i+1) -p_s(i-3) )/3.0/d_x,(p_g(i+1) -p_g(i) )/d_x);
        
    [p_g(J+1),u_g(J+1),lo_g(J+1),dp_g(J+1),du_g(J+1),dlo_g(J+1),p_s(J),u_s(J),lo_s(J),dp_s(J),du_s(J),dlo_s(J),u_mat]=ghost_cal_GRP(lo_s(J-1)+d_x*dlo_s(J-1),u_s(J-1)+d_x*du_s(J-1),p_s(J-1)+d_x*dp_s(J-1),dlo_s(J-1),du_s(J-1),dp_s(J-1),gama_s,lo_g(J+2)-d_x*dlo_g(J+2),u_g(J+2)-d_x*du_g(J+2),p_g(J+2)-d_x*dp_g(J+2),dlo_g(J+2),du_g(J+2),dp_g(J+2),gama_g);
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
    [p_sM,u_sM,lo_sM]=value_cal_GRP(lo_s(J),u_s(J),p_s(J),dlo_s(J),du_s(J),dp_s(J),gama_s,0,d_t);
    [p_gM,u_gM,lo_gM]=value_cal_GRP(lo_g(J+1),u_g(J+1),p_g(J+1),dlo_g(J+1),du_g(J+1),dp_g(J+1),gama_g,0,d_t);
    for i=(J-2):(J+1)
        E_g(i)=p_gM/(gama_g-1)+0.5*lo_gM*u_gM^2;
        U_g(:,i)=[lo_gM;lo_gM*u_gM;E_g(i)];
    end
    for i=J:(J+3)
        E_s(i)=p_sM/(gama_s-1)+0.5*lo_sM*u_sM^2;
        U_s(:,i)=[lo_sM;lo_sM*u_sM;E_s(i)];
    end
    for i=2:(N-1)
        dphi(i) =minmod(Alpha,(phi(i) -phi(i-1) )/d_x,(phi(i+1) -phi(i-1) )/2.0/d_x,(phi(i+1) -phi(i) )/d_x);
    end
%     for i=2:J
%         dphi(i) =minmod(Alpha,(phi(i) -phi(i-1) )/d_x,(W_int_s(4,i+1)-W_int_s(4,i))/d_x,(phi(i+1) -phi(i) )/d_x);
%     end
%     for i=(J+1):(N-1)
%         dphi(i) =minmod(Alpha,(phi(i) -phi(i-1) )/d_x,(W_int_g(4,i+1)-W_int_g(4,i))/d_x,(phi(i+1) -phi(i) )/d_x);
%     end
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
         if i==1
            [F_s(:,1),W_int_s(:,1),phi_s(1)]=GRP_solver(lo_s(1),lo_s(1),0,0,u_s(1),u_s(1),0,0,p_s(1),p_s(1),0,0,phi(1),phi(1),0,0,gama_s,d_t);
         elseif i==N+1
            [F_g(:,N+1),W_int_g(:,N+1),phi_g(N+1)]=GRP_solver(lo_g(N),lo_g(N),0,0,u_g(N),u_g(N),0,0,p_g(N),p_g(N),0,0,phi(N),phi(N),0,0,gama_g,d_t);
         else
            if i<=J
                [F_s(:,i),W_int_s(:,i),phi_s(i)]=GRP_solver(lo_s(i-1)+0.5*d_x*dlo_s(i-1),lo_s(i)-0.5*d_x*dlo_s(i),dlo_s(i-1),dlo_s(i),u_s(i-1)+0.5*d_x*du_s(i-1),u_s(i)-0.5*d_x*du_s(i),du_s(i-1),du_s(i),p_s(i-1)+0.5*d_x*dp_s(i-1),p_s(i)-0.5*d_x*dp_s(i),dp_s(i-1),dp_s(i),phi(i-1)+0.5*d_x*dphi(i-1),phi(i)-0.5*d_x*dphi(i),dphi(i-1),dphi(i),gama_s,d_t);
            elseif i>=J+2
                [F_g(:,i),W_int_g(:,i),phi_g(i)]=GRP_solver(lo_g(i-1)+0.5*d_x*dlo_g(i-1),lo_g(i)-0.5*d_x*dlo_g(i),dlo_g(i-1),dlo_g(i),u_g(i-1)+0.5*d_x*du_g(i-1),u_g(i)-0.5*d_x*du_g(i),du_g(i-1),du_g(i),p_g(i-1)+0.5*d_x*dp_g(i-1),p_g(i)-0.5*d_x*dp_g(i),dp_g(i-1),dp_g(i),phi(i-1)+0.5*d_x*dphi(i-1),phi(i)-0.5*d_x*dphi(i),dphi(i-1),dphi(i),gama_g,d_t);
            end
         end
    end
%     if u_g(J)>0
%         phi_s(J+1)=phi(J);
%         phi_g(J+1)=phi(J);
%     else
%         phi_s(J+1)=phi(J+1);
%         phi_g(J+1)=phi(J+1);
%     end
    %compute U in next step
    x0=x0+u_s(J)*d_t+0.5*u_mat*d_t^2;  
    for i=1:(J-1)
        U_s(:,i)=U_s(:,i)+d_t/d_x*(F_s(:,i)-F_s(:,i+1));
        [lo_s(i),u_s(i),p_s(i)]=primitive_comp(U_s(:,i),gama_s);
%           phi(i)  =phi(i)+d_t/d_x*(u_s(i)+0.5*d_t*(-dp_s(i)/lo_s(i)-u_s(i)*du_s(i)))*(phi_s(i)-phi_s(i+1));
%           U_phi(i)=U_phi(i)+d_t/d_x*(phi_s(i)-phi_s(i+1));
%         	phi(i)=U_phi(i)/lo_s(i);
    end
    for i=(J+2):N
        U_g(:,i)=U_g(:,i)+d_t/d_x*(F_g(:,i)-F_g(:,i+1));
        [lo_g(i),u_g(i),p_g(i)]=primitive_comp(U_g(:,i),gama_g);
%            phi(i)  =phi(i)+d_t/d_x*(u_g(i)+0.5*d_t*(-dp_g(i)/lo_g(i)-u_g(i)*du_g(i)))*(phi_g(i)-phi_g(i+1));
%            U_phi(i)=U_phi(i)+d_t/d_x*(phi_g(i)-phi_g(i+1));
%            phi(i)=U_phi(i)/lo_g(i);
    end
%     if u_last>0.0
%         x0=x0+0.5*(u_last+u_s(J+1))*d_t;
%     else
%         x0=x0+0.5*(u_last+u_g(J)  )*d_t;
%     end
    for i=1:N
        phi(i)=x(i)-x0;
    end
    for i=1:N-1
        if (phi(i)*phi(i+1))<0.0
            J=i;
            break;
        end
    end
    count=count+1;
    Time=Time+d_t
% if Time > 4*d_t
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
plot(x_min:d_x:x_max-d_x,lo,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Density','FontWeight','bold');
subplot(2,2,2);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,2),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,u,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Velocity','FontWeight','bold');
subplot(2,2,3);
hold on
plot(x_min:d_x:x_max-d_x,W_exact(:,3),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,p,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Pressure','FontWeight','bold');
subplot(2,2,4);
hold on
%plot(x_min:d_x:x_max-d_x,W_exact(:,4),'k','LineWidth',1.0);
plot(x_min:d_x:x_max-d_x,phi,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('level-set','FontWeight','bold');
