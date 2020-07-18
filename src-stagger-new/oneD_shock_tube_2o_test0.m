clear;
clc;
warning off;
%1D shock_tube by 2-order staggered Schemes for BN model
%state constant
global gama_s gama_g p0;
gama_s=1.4;
gama_g=1.4;
p0=0;
global ep;
ep=1e-9;
x_min=0;
x_max=1;
N=800*1;
d_x=(x_max-x_min)/N;
x0=0.5;
CFL=0.4;
Alpha_G=0.0;
Alpha_GRP=0.0;
%state value
Time=0;
Tend=0.1;
%Tend=0.15;
Alpha=zeros(1,N+1);
U=zeros(6,N);
F=zeros(6,N+1);
F_0=zeros(6,N+1);
RI=zeros(6,N); %Riemann invariants
d_RI=zeros(6,N);
d_Alpha=zeros(1,N+1);
dd_Alpha=zeros(1,N+1);
U_lo_sL=zeros(1,N);
U_lo_sR=zeros(1,N);
%initial condition
lo_g_0  =1.0;
u_g_0   =0.0;
p_g_0   =1.0;
lo_s_0  =1.0;
u_s_0   =@(x) 0.5 + 0.5*tanh(20*x-10);
p_s_0   =1.0;
% phi_s_0 =@(x) 0.5 + 0.4*tanh(20*x-8);
% phi_g_0 =@(x) 0.5 - 0.4*tanh(20*x-8);
phi_s_0 =@(x) 0.5+0*x;
phi_g_0 =@(x) 0.5+0*x;
for i=1:(N+1)
    Alpha(i) = integral(phi_s_0,x_min+(i-1.5)*d_x,x_min+(i-0.5)*d_x)/d_x;
end
E_g_0 = p_g_0/(gama_g-1)+0.5*lo_g_0*u_g_0^2;
E_s_0 = @(x) (p_s_0+gama_s*p0)/(gama_s-1)+0.5*lo_s_0*u_s_0(x).^2;
U_0_1 = @(x) phi_g_0(x)*lo_g_0;
U_0_2 = @(x) phi_g_0(x)*lo_g_0*u_g_0;
U_0_3 = @(x) phi_g_0(x)*E_g_0;
U_0_4 = @(x) phi_s_0(x)*lo_s_0;
U_0_5 = @(x) phi_s_0(x)*lo_s_0.*u_s_0(x);
U_0_6 = @(x) phi_s_0(x).*E_s_0(x);
for i=1:N
    U_1_1 = integral(U_0_1,x_min+(i-1)*d_x,x_min+i*d_x)/d_x;
    U_1_2 = integral(U_0_2,x_min+(i-1)*d_x,x_min+i*d_x)/d_x;
    U_1_3 = integral(U_0_3,x_min+(i-1)*d_x,x_min+i*d_x)/d_x;
    U_1_4 = integral(U_0_4,x_min+(i-1)*d_x,x_min+i*d_x)/d_x;
    U_1_5 = integral(U_0_5,x_min+(i-1)*d_x,x_min+i*d_x)/d_x;
    U_1_6 = integral(U_0_6,x_min+(i-1)*d_x,x_min+i*d_x)/d_x;   
    U(:,i) = [U_1_1;U_1_2;U_1_3;U_1_4;U_1_5;U_1_6];
end

W_int_s=zeros(4,N+1);
W_int_g=zeros(4,N+1);
dlo_s=zeros(1,N+1);
du_s =zeros(1,N+1);
dp_s =zeros(1,N+1);
dlo_g=zeros(1,N+1);
du_g =zeros(1,N+1);
dp_g =zeros(1,N+1);
IDX  =ones(1,N+1);
W_RI_int=zeros(7,N+1);
%test begin
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
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
        RI(:,i)=U2RI_cal(Alpha(i),lo_gL(i),u_gL(i),p_gL(i),u_sL(i),p_sL(i),lo_sL(i));
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
    %reconstruction (minmod limiter)
    for i=2:N
%         d_Alpha(i)=minmod2((Alpha(:,i)-Alpha(:,i-1))/d_x,(Alpha(:,i+1)-Alpha(:,i))/d_x);
        d_Alpha(i)=minmod(Alpha_G*(Alpha(:,i)-Alpha(:,i-1))/d_x,(Alpha(:,i+1)-Alpha(:,i-1))/2.0/d_x,Alpha_G*(Alpha(:,i+1)-Alpha(:,i))/d_x);
%          d_Alpha(i)=(Alpha(:,i+1)-Alpha(:,i-1))/2.0/d_x;
         dd_Alpha(i)= abs((Alpha(:,i+1)-Alpha(:,i-1))/2.0/d_x);
    end
    HN=1;
    for i=HN+1:N-HN
        d_RI(:,i)=minmod2((RI(:,i)-RI(:,i-1))/d_x,(RI(:,i+1)-RI(:,i))/d_x);
        d_RI(:,i)=minmod(Alpha_G*(RI(:,i)-RI(:,i-1))/d_x,(RI(:,i+1)-RI(:,i-1))/2.0/d_x,Alpha_G*(RI(:,i+1)-RI(:,i))/d_x);
        dlo_s(i) =minmod(Alpha_GRP*(lo_sL(i)-lo_sL(i-1))/d_x,(lo_sL(i+1)-lo_sL(i-1))/2.0/d_x,Alpha_GRP*(lo_sL(i+1)-lo_sL(i))/d_x);
        du_s(i)  =minmod(Alpha_GRP*(u_sL(i) -u_sL(i-1) )/d_x,(u_sL(i+1) -u_sL(i-1) )/2.0/d_x,Alpha_GRP*(u_sL(i+1) -u_sL(i) )/d_x);
        dp_s(i)  =minmod(Alpha_GRP*(p_sL(i) -p_sL(i-1) )/d_x,(p_sL(i+1) -p_sL(i-1) )/2.0/d_x,Alpha_GRP*(p_sL(i+1) -p_sL(i) )/d_x);
        dlo_g(i) =minmod(Alpha_GRP*(lo_gL(i)-lo_gL(i-1))/d_x,(lo_gL(i+1)-lo_gL(i-1))/2.0/d_x,Alpha_GRP*(lo_gL(i+1)-lo_gL(i))/d_x);
        du_g(i)  =minmod(Alpha_GRP*(u_gL(i) -u_gL(i-1) )/d_x,(u_gL(i+1) -u_gL(i-1) )/2.0/d_x,Alpha_GRP*(u_gL(i+1) -u_gL(i) )/d_x);
        dp_g(i)  =minmod(Alpha_GRP*(p_gL(i) -p_gL(i-1) )/d_x,(p_gL(i+1) -p_gL(i-1) )/2.0/d_x,Alpha_GRP*(p_gL(i+1) -p_gL(i) )/d_x);
        
%         d_RI(:,i)= (RI(:,i+1)-RI(:,i-1))/2.0/d_x;
%         dlo_s(i) = (lo_sL(i+1)-lo_sL(i-1))/2.0/d_x;
%         du_s(i)  = (u_sL(i+1) -u_sL(i-1) )/2.0/d_x;
%         dp_s(i)  = (p_sL(i+1) -p_sL(i-1) )/2.0/d_x;
%         dlo_g(i) = (lo_gL(i+1)-lo_gL(i-1))/2.0/d_x;
%         du_g(i)  = (u_gL(i+1) -u_gL(i-1) )/2.0/d_x;
%         dp_g(i)  = (p_gL(i+1) -p_gL(i-1) )/2.0/d_x;
%         IDX(i)=0;
%         if max(dd_Alpha(i-HN:i+HN-1)>ep)
            IDX(i)=1;
%         end
    end
    for i=1:N+1
        %flux on the boundary of i-1 and i
      if i>1 && i<N+1
         [F_0(1:3,i),W_int_g(:,i)]=GRP_solver(lo_gL(i-1)+0.5*d_x*dlo_g(i-1),lo_gL(i)-0.5*d_x*dlo_g(i),dlo_g(i-1),dlo_g(i),u_gL(i-1)+0.5*d_x*du_g(i-1),u_gL(i)-0.5*d_x*du_g(i),du_g(i-1),du_g(i),p_gL(i-1)+0.5*d_x*dp_g(i-1),p_gL(i)-0.5*d_x*dp_g(i),dp_g(i-1),dp_g(i),1-Alpha(i-1),1-Alpha(i),0.0,0.0,gama_g,d_t);
         [F_0(4:6,i),W_int_s(:,i)]=GRP_solver(lo_sL(i-1)+0.5*d_x*dlo_s(i-1),lo_sL(i)-0.5*d_x*dlo_s(i),dlo_s(i-1),dlo_s(i),u_sL(i-1)+0.5*d_x*du_s(i-1),u_sL(i)-0.5*d_x*du_s(i),du_s(i-1),du_s(i),p_sL(i-1)+0.5*d_x*dp_s(i-1),p_sL(i)-0.5*d_x*dp_s(i),dp_s(i-1),dp_s(i),  Alpha(i-1),  Alpha(i),0.0,0.0,gama_s,d_t);
      end
         if i==1
             Alpha_int=Alpha(1);
             [lo_gL_i,u_gL_i,p_gL_i,u_sL_i,p_sL_i(1),lo_sL_i]=RI2U_cal(Alpha_int,RI(:,1),lo_gL(1));
             [lo_gR_i,u_gR_i,p_gR_i,u_sR_i,p_sR_i(1),lo_sR_i]=RI2U_cal(Alpha_int,RI(:,1),lo_gL(1));
             d_RI_L=d_RI(:,1);
             d_RI_R=d_RI(:,1);
         elseif i==N+1
             Alpha_int=Alpha(N+1);
             [lo_gL_i,u_gL_i,p_gL_i,u_sL_i,p_sL_i(N+1),lo_sL_i]=RI2U_cal(Alpha_int,RI(:,N),lo_gR(N));
             [lo_gR_i,u_gR_i,p_gR_i,u_sR_i,p_sR_i(N+1),lo_sR_i]=RI2U_cal(Alpha_int,RI(:,N),lo_gR(N));
             d_RI_L=d_RI(:,N);
             d_RI_R=d_RI(:,N);
         else
             Alpha_int=Alpha(i);
             d_RI_L=d_RI(:,i-1);
             d_RI_R=d_RI(:,i);
             [lo_gL_i,u_gL_i,p_gL_i,u_sL_i,p_sL_i(i),lo_sL_i]=RI2U_cal(Alpha_int,RI(:,i-1)+0.5*d_x*d_RI_L,lo_gR(i-1));
             [lo_gR_i,u_gR_i,p_gR_i,u_sR_i,p_sR_i(i),lo_sR_i]=RI2U_cal(Alpha_int,  RI(:,i)-0.5*d_x*d_RI_R,  lo_gL(i));
        end
%        F(1:3,i)=Riemann_solver_Exact(lo_gL_i,lo_gR_i,p_gL_i,   p_gR_i,   u_gL_i,u_gR_i,1-Alpha_int,gama_g,0.0);
%        F(4:6,i)=Riemann_solver_Exact(lo_sL_i,lo_sR_i,p_sL_i(i),p_sR_i(i),u_sL_i,u_sR_i,  Alpha_int,gama_s,0.0);
         [Alpha_mid(i),u_s_mid(i),F(:,i),W_RI_int(:,i)]=GRP_RI_solver(lo_gL_i,lo_gR_i,u_gL_i,u_gR_i,p_gL_i,p_gR_i,lo_sL_i,lo_sR_i,u_sL_i,u_sR_i,p_sL_i(i),p_sR_i(i),Alpha_int,Alpha_int,[d_Alpha(i);d_RI_L],[d_Alpha(i);d_RI_R],d_t);    
    end
    F_0(1:6,1)=F(1:6,1);
    F_0(1:6,N+1)=F(1:6,N+1);
    for i=1:N
        if i<N
            rho_s(i) =0.5*(lo_sR(i)+lo_sL(i+1));
            arho_s(i)=rho_s(i)*Alpha(i+1);           
        end
        [lo_gL_m,u_gL_m,p_gL_m,u_sL_m,p_sL_m,lo_sL_m]=RI2U_cal(Alpha(i)+0.5*d_x*d_Alpha(:,i),    RI(:,i),lo_gL(i));
        [lo_gR_m,u_gR_m,p_gR_m,u_sR_m,p_sR_m,lo_sR_m]=RI2U_cal(Alpha(i+1)-0.5*d_x*d_Alpha(:,i+1),RI(:,i),lo_gR(i));
        [phi_s_mL, lo_s_mL, u_s_mL, p_g_mL(i)]=GRP_RI_solver_mid(lo_gL_m,u_gL_m,p_gL_m,lo_sL_m,u_sL_m,p_sL_m,Alpha(:,i)+0.5*d_x*d_Alpha(:,i),    [d_Alpha(i);  d_RI(:,i)],d_t);
        [phi_s_mR, lo_s_mR, u_s_mR, p_g_mR(i)]=GRP_RI_solver_mid(lo_gR_m,u_gR_m,p_gR_m,lo_sR_m,u_sR_m,p_sR_m,Alpha(:,i+1)-0.5*d_x*d_Alpha(:,i+1),[d_Alpha(i+1);d_RI(:,i)],d_t);
        if (u_sL(i)>0.0)
            phi_s_m(i) = phi_s_mL;
            lo_s_m = lo_s_mL;
            u_s_m(i) = u_s_mL;
        else
            phi_s_m(i) = phi_s_mR;
            lo_s_m = lo_s_mR;
            u_s_m(i) = u_s_mR;
        end
        F_rho_s(i) =lo_s_m*u_s_m(i);
        F_arho_s(i)=phi_s_m(i)*lo_s_m*u_s_m(i);
    end
    for i=1:N-1
        arho_s(i)=arho_s(i)+d_t/d_x*(F_arho_s(i)-F_arho_s(i+1));
        rho_s(i) = rho_s(i)+d_t/d_x*(F_rho_s(i) -F_rho_s(i+1));
        Alpha_next(i+1)=arho_s(i)/rho_s(i);
    end
    
    for i=1:N-1
        Alpha_next(i+1)=W_RI_int(1,i+1);
    end
    
    %compute U in next step
    for i=1:N
%       if abs(Alpha(i+1)-Alpha(i))<ep
          S=0.5*(p_g_mL(i)+p_g_mR(i))*(Alpha(i+1)-Alpha(i));
%       else
%           S_tmp=Alpha_mid(i+1)*p_sL_i(i+1)-Alpha_mid(i)*p_sR_i(i);
%           if (S_tmp/(Alpha(i+1)-Alpha(i))>max(p_g_mL(i),p_g_mR(i)))
%               S=max(p_g_mL(i),p_g_mR(i))*(Alpha(i+1)-Alpha(i));
%           elseif (S_tmp/(Alpha(i+1)-Alpha(i))<min(p_g_mL(i),p_g_mR(i)))
%               S=min(p_g_mL(i),p_g_mR(i))*(Alpha(i+1)-Alpha(i));
%           else
%               S=S_tmp;
%           end
%       end
      if IDX(i) == 0
        U(:,i)=U(:,i)+d_t/d_x*(F_0(:,i)-F_0(:,i+1));
        Alpha_next(i) = Alpha(i);        
      else
        U(:,i)=U(:,i)+d_t/d_x*(F(:,i)-F(:,i+1))+d_t/d_x*[0;-S;-S*u_s_m(i);0;S;S*u_s_m(i)];
%         area_L=0.5+u_s_m(i)*d_t/d_x;
%         area_R=1.0-area_L;
%         Alpha_starL = Alpha(i) - u_s_mid(i)*d_t/d_x/area_L*(Alpha(i)-Alpha_mid(i));
%         Alpha_starR = Alpha(i+1) - u_s_mid(i)*d_t/d_x/area_R*(Alpha(i+1)-Alpha_mid(i+1));
%         [lo_gL(i),u_gL(i),p_gL(i),lo_sL(i),u_sL(i),p_sL(i),lo_gR(i),u_gR(i),p_gR(i),lo_sR(i),u_sR(i),p_sR(i)]=primitive_comp(U(:,i),Alpha_starL,Alpha_starR,area_L,area_R);
%         phi_sL=Alpha_next(i);
%         phi_sR=Alpha_next(i+1);
%         [lo_gL(i),u_gL(i),p_gL(i),p_sL(i)]=Riemann_inv(Alpha(i),  lo_gL(i),u_gL(i),p_gL(i),u_sL(i),p_sL(i),phi_sL);   
%         [lo_gR(i),u_gR(i),p_gR(i),p_sR(i)]=Riemann_inv(Alpha(i+1),lo_gR(i),u_gR(i),p_gR(i),u_sR(i),p_sR(i),phi_sR);   
%         phi_gL=1.0-phi_sL;
%         phi_gR=1.0-phi_sR;
%         E_gL=p_gL(i)/(gama_g-1)+0.5*lo_gL(i)*u_gL(i)^2;
%         E_sL=(p_sL(i)+gama_s*p0)/(gama_s-1)+0.5*lo_sL(i)*u_sL(i)^2;
%         U(:,i)=       0.5*[phi_gL*lo_gL(i);phi_gL*lo_gL(i)*u_gL(i);phi_gL*E_gL;phi_sL*lo_sL(i);phi_sL*lo_sL(i)*u_sL(i);phi_sL*E_sL];
%         E_gR=p_gR(i)/(gama_g-1)+0.5*lo_gR(i)*u_gR(i)^2;
%         E_sR=(p_sR(i)+gama_s*p0)/(gama_s-1)+0.5*lo_sR(i)*u_sR(i)^2;
%         U(:,i)=U(:,i)+0.5*[phi_gR*lo_gR(i);phi_gR*lo_gR(i)*u_gR(i);phi_gR*E_gR;phi_sR*lo_sR(i);phi_sR*lo_sR(i)*u_sR(i);phi_sR*E_sR];
      end
    end
    Alpha=Alpha_next;
    Time=Time+d_t
%     if Time > 100*d_t
%         break;
%     end
end
% U_big = U;
% save('U-exact.mat','U_big');

N_big = 3200;
rat = N_big/N;
load U-exact.mat;
U_exact=zeros(6,N_big);
err_U=zeros(1,N);
for i=1:N
     U_exact(:,i) = sum(U_big(:,(i-1)*rat+1:i*rat),2)/rat;
     err_U(i) = norm(U_exact(:,i)-U(:,i));
end
for i=1:N
     err_U_1(i) = norm(U_exact(1,i)-U(1,i));
     err_U_2(i) = norm(U_exact(2,i)-U(2,i));
     err_U_3(i) = norm(U_exact(3,i)-U(3,i));
     err_U_4(i) = norm(U_exact(4,i)-U(4,i));
     err_U_5(i) = norm(U_exact(5,i)-U(5,i));
     err_U_6(i) = norm(U_exact(6,i)-U(6,i));  
end
E_N = sum(err_U)*d_x
E_N_1 = sum(err_U_1)*d_x
E_N_2 = sum(err_U_2)*d_x
E_N_3 = sum(err_U_3)*d_x
E_N_4 = sum(err_U_4)*d_x
E_N_5 = sum(err_U_5)*d_x
E_N_6 = sum(err_U_6)*d_x
E_N_all = [E_N_1,E_N_2,E_N_3,E_N_4,E_N_5,E_N_6];
