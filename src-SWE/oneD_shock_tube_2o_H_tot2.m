clear;
clc;
%1D shock_tube by GRP Schemes for shallow water model
%state constant
global g;
global ep;
ep=1e-6;
x0=0;

%state value
Time=0;
CFL=0.5;
Alpha_GRP=1.5;

%initial condition
init_discon
%init_con

u_mid  =zeros(1,N+1);
H_t_mid=zeros(1,N+1);
dH_t   =zeros(1,N);

%Godunov's Method
while Time<Tend && isreal(Time)
    %boundary condition
    %U(2,1)=q_in;
    %U(1,1)=Z_0;
    %U(1,N)=h_out;
    %CFL condition
    for i=1:N
        [h_L(i),u_L(i),h_R(i),u_R(i),H_t(i)]=primitive_comp(U(:,i),Z_L(i),Z_R(i));
        a_L(i)=sqrt(g*h_L(i));
        a_R(i)=sqrt(g*h_R(i));
        if u_L(i) > a_L(i)
            Fr_L(i)=2;
        else
            Fr_L(i)=0;
        end
        if u_R(i) > a_R(i)
            Fr_R(i)=2;
        else
            Fr_R(i)=0;
        end
    end
    hh = U(1,:);
    qq = U(2,:);
    Smax=max(max(abs(u_L)+a_L),max(abs(u_R)+a_R));
    d_t=CFL*d_x/Smax;
    if Time+d_t >= Tend
        d_t = Tend-Time+1e-10;
    end
    %reconstruction (minmod limiter)
    for i=2:N-1
%            dh(i) =minmod(Alpha_GRP*(h_R(i)-h_R(i-1))/d_x,(h_L(i+1)-h_L(i-1))/2.0/d_x,      Alpha_GRP*(h_L(i+1)-h_L(i))/d_x);
%            du(i) =minmod(Alpha_GRP*(u_R(i)-u_R(i-1))/d_x,(u_L(i+1)-u_L(i-1))/2.0/d_x,      Alpha_GRP*(u_L(i+1)-u_L(i) )/d_x);
%            dh(i) =minmod(Alpha_GRP*(h_R(i)-h_R(i-1))/d_x,(W_int(1,i+1)-W_int(1,i))/1.0/d_x,Alpha_GRP*(h_L(i+1)-h_L(i))/d_x);
%            du(i) =minmod(Alpha_GRP*(u_R(i)-u_R(i-1))/d_x,(W_int(2,i+1)-W_int(1,i))/1.0/d_x,Alpha_GRP*(u_L(i+1)-u_L(i))/d_x);
        if Time < ep
            dq(i)   =minmod(Alpha_GRP*(qq(i) -qq(i-1)) /d_x,(qq(i+1)-qq(i-1))/2.0/d_x,  Alpha_GRP*(qq(i+1) -qq(i)) /d_x);
            dH_t(i) =minmod(Alpha_GRP*(H_t(i)-H_t(i-1))/d_x,(H_t(i+1)-H_t(i-1))/2.0/d_x,Alpha_GRP*(H_t(i+1)-H_t(i))/d_x);
        else
            dq(i)   =minmod(Alpha_GRP*(qq(i) -qq(i-1)) /d_x,(W_int(2,i+1)-W_int(2,i))/1.0/d_x,Alpha_GRP*(qq(i+1) -qq(i)) /d_x);
            dH_t(i) =minmod(Alpha_GRP*(H_t(i)-H_t(i-1))/d_x,(W_int(4,i+1)-W_int(4,i))/1.0/d_x,Alpha_GRP*(H_t(i+1)-H_t(i))/d_x);
        end
    end
    for i=1:N+1
        if i==1
            h_L_int(i)=h_R(1);
            u_L_int(i)=u_R(1);
            h_R_int(i)=h_R(1);
            u_R_int(i)=u_R(1);           
            dh_L_int(i)=0;
            du_L_int(i)=0;
        elseif i==N+1 
            h_L_int(i)=h_L(N);
            u_L_int(i)=u_L(N);
            h_R_int(i)=h_L(N);
            u_R_int(i)=u_L(N);
            dh_R_int(i)=0;
            du_R_int(i)=0;
        else
            [h_L_int(i),u_L_int(i),dh_L_int(i),du_L_int(i)]=dRI2dU_cal(qq(i-1)+0.5*d_x*dq(i-1),H_t(i-1)+0.5*d_x*dH_t(i-1)-Z_M(i),dq(i-1),dH_t(i-1),dZ(i-1),Fr_R(i-1));
            [h_R_int(i),u_R_int(i),dh_R_int(i),du_R_int(i)]=dRI2dU_cal(qq(i)  -0.5*d_x*dq(i),  H_t(i)  -0.5*d_x*dH_t(i)  -Z_M(i),dq(i),  dH_t(i),  dZ(i-1),Fr_L(i));
        end
    end    
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
        if i==1
            dZZ=0.0;
        elseif i==N+1
            dZZ=0.0;
        else
            dZZ=dZ(i-1);            
        end
        [h_mid(:,i),u_mid(:,i),H_t_mid(:,i),F(:,i),W_int(:,i)]=GRP_solver(h_L_int(i),h_R_int(i),dh_L_int(i),dh_R_int(i),u_L_int(i),u_R_int(i),du_L_int(i),du_R_int(i),Z_M(i),dZZ,dZZ,d_t);
    end
    %compute U in next step
    for i=1:N
        S = -0.5*(h_mid(:,i+1)*u_mid(:,i+1)-h_mid(:,i)*u_mid(:,i))*(u_mid(:,i)+u_mid(:,i+1));
        S = S - 0.5*g*(H_t_mid(:,i+1)-H_t_mid(:,i))*(h_mid(:,i)+h_mid(:,i+1));
        U(1,i)=U(1,i)+d_t/d_x*(F(1,i)-F(1,i+1));
        U(2,i)=U(2,i)+d_t/d_x*S;
    end
    Time = Time+d_t
% if Time > 0.002
%    break;
% end
end

plot_SWE
