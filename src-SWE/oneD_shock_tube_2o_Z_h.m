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
%init_discon
init_con

h_mL =zeros(1,N);
h_mR =zeros(1,N);
u_mL =zeros(1,N);
u_mR =zeros(1,N);
dzeta=zeros(1,N+1);

%Godunov's Method
while Time<Tend && isreal(Time)
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
        if Time < ep
            dzeta(i) =minmod(Alpha_GRP*(h_R(i)+Z_R(i)-h_R(i-1)-Z_R(i-1))/d_x,(h_L(i+1)+Z_L(i+1)-h_R(i-1)-Z_R(i-1))/2.0/d_x,Alpha_GRP*(h_L(i+1)+Z_L(i+1)-h_L(i)-Z_L(i))/d_x);
            dzeta(i) =minmod(Alpha_GRP*(h_L(i)+Z_L(i)-h_R(i-1)-Z_R(i-1))/d_x,dzeta(i),                                     Alpha_GRP*(h_L(i+1)+Z_L(i+1)-h_R(i)-Z_R(i))/d_x);
            dq(i)    =minmod(Alpha_GRP*(qq(i) -qq(i-1)) /d_x,                (qq(i+1) -qq(i-1)) /2.0/d_x,                  Alpha_GRP*(qq(i+1) -qq(i)) /d_x);
        else
            dzeta(i) =minmod(Alpha_GRP*(h_R(i)+Z_R(i)-h_R(i-1)-Z_R(i-1))/d_x,(W_int(3,i+1)-W_int(3,i))/1.0/d_x,Alpha_GRP*(h_L(i+1)+Z_L(i+1)-h_L(i)-Z_L(i))/d_x);
            dzeta(i) =minmod(Alpha_GRP*(h_L(i)+Z_L(i)-h_R(i-1)-Z_R(i-1))/d_x,dzeta(i),                         Alpha_GRP*(h_L(i+1)+Z_L(i+1)-h_R(i)-Z_R(i))/d_x);
            dq(i)    =minmod(Alpha_GRP*(qq(i) -qq(i-1)) /d_x,                (W_int(2,i+1)-W_int(2,i))/1.0/d_x,Alpha_GRP*(qq(i+1) -qq(i)) /d_x);
        end
    end
    if Time >= ep        
        dzeta(1) =minmod(Alpha_GRP*(h_R(1)+Z_R(1)-W_int(3,1))/d_x/2,     (W_int(3,2)-W_int(3,1))/1.0/d_x,  Alpha_GRP*(h_L(2)+Z_L(2)-h_L(1)-Z_L(1))/d_x);
        dzeta(1) =minmod(Alpha_GRP*(h_L(1)+Z_L(1)-W_int(3,1))/d_x/2,     dzeta(1),                         Alpha_GRP*(h_L(2)+Z_L(2)-h_R(1)-Z_R(1))/d_x);
        dq(1)    =minmod(Alpha_GRP*(qq(1) -W_int(2,1)) /d_x/2,           (W_int(2,2)-W_int(2,1))/1.0/d_x,  Alpha_GRP*(qq(2) -qq(1)) /d_x); 
        dzeta(N) =minmod(Alpha_GRP*(h_R(N)+Z_R(N)-h_R(N-1)-Z_R(N-1))/d_x,(W_int(3,N+1)-W_int(3,N))/1.0/d_x,Alpha_GRP*(W_int(3,N+1)-h_L(N)-Z_L(N))/d_x/2);
        dzeta(N) =minmod(Alpha_GRP*(h_L(N)+Z_L(N)-h_R(N-1)-Z_R(N-1))/d_x,dzeta(N),                         Alpha_GRP*(W_int(3,N+1)-h_R(N)-Z_R(N))/d_x/2);
        dq(N)    =minmod(Alpha_GRP*(qq(N) -qq(N-1)) /d_x,                (W_int(2,N+1)-W_int(2,N))/1.0/d_x,Alpha_GRP*(W_int(2,N+1)-qq(N)) /d_x/2); 
    end    
    for i=1:N+1
        if i < N+1
            [h_R_int(i),u_R_int(i),dh_R_int(i),du_R_int(i)]=dqzeta2dU_cal(qq(i)  -0.5*d_x*dq(i),  h_L(i)  +Z_L(i)  -0.5*d_x*dzeta(i)  -Z_M(i),dq(i),  dzeta(i),  dZ(i),Fr_L(i));
        end
        if i > 1
            [h_L_int(i),u_L_int(i),dh_L_int(i),du_L_int(i)]=dqzeta2dU_cal(qq(i-1)+0.5*d_x*dq(i-1),h_R(i-1)+Z_R(i-1)+0.5*d_x*dzeta(i-1)-Z_M(i),dq(i-1),dzeta(i-1),dZ(i),Fr_R(i-1));
        end
    end
    %Riemann problem:compute flux
    for i=1:N+1
        %flux on the boundary of i-1 and i
        if i==1 %boundary conditionh
            [h_mid(:,i),u_mid_0,H_t_mid_0,F(:,i),W_int(:,1)]=GRP_solver_inflow( q_in ,h_R_int(i),dh_R_int(i),u_R_int(i),du_R_int(i),Z_M(i),dZ(i),dZ(i),d_t);
        elseif i==N+1
            [h_mid(:,i),u_mid_0,H_t_mid_0,F(:,i),W_int(:,i)]=GRP_solver_outflow(h_out,h_L_int(i),dh_L_int(i),u_L_int(i),du_L_int(i),Z_M(i),dZ(i),dZ(i),d_t);
        else
            [h_mid(:,i),u_mid_0,H_t_mid_0,F(:,i),W_int(:,i)]=GRP_solver(h_L_int(i),h_R_int(i),dh_L_int(i),dh_R_int(i),u_L_int(i),u_R_int(i),du_L_int(i),du_R_int(i),Z_M(i),dZ(i),dZ(i),d_t);
        end
    end
    for i=1:N
        if i==1 || i==N
            h_mL(i) = h_L(i);
            h_mR(i) = h_R(i);            
            u_mL(i) = u_L(i);
            u_mR(i) = u_R(i); 
        else
            [h_mL(i),u_mL(i),dh_mL,du_mL]=dqzeta2dU_cal(qq(i),h_L(i),dq(i-1),dzeta(i-1),dZ(i),Fr_L(i));
            [h_mR(i),u_mR(i),dh_mR,du_mR]=dqzeta2dU_cal(qq(i),h_R(i),dq(i-1),dzeta(i-1),dZ(i+1),  Fr_R(i));
            h_mL(i) = h_mL(i) - 0.5*d_t*(h_mL(i)*du_mL+u_mL(i)*dh_mL);
            h_mR(i) = h_mR(i) - 0.5*d_t*(h_mR(i)*du_mR+u_mR(i)*dh_mR);        
            u_mL(i) = u_mL(i) - 0.5*d_t*(dh_mL*g+u_mL(i)*du_mL+dZ(i)  *g);
            u_mR(i) = u_mR(i) - 0.5*d_t*(dh_mR*g+u_mR(i)*du_mR+dZ(i+1)*g);
        end
    end
    %compute U in next step
    for i=1:N
        if abs(Z_R(i)-Z_L(i))<ep
            S=-g*0.5*(h_mL(i)+h_mR(i))*(Z_R(i)-Z_L(i));
        else
            S_tmp=(h_mR(i)*u_mR(i)^2+g*h_mR(i)^2/2-h_mL(i)*u_mL(i)^2-g*h_mL(i)^2/2);
            if (S_tmp/g/(Z_L(i)-Z_R(i))>max(h_mL(i),h_mR(i)))
                S=-g*max(h_mL(i),h_mR(i))*(Z_R(i)-Z_L(i));        
            elseif (S_tmp/g/(Z_L(i)-Z_R(i))<min(h_mL(i),h_mR(i)))
                S=-g*min(h_mL(i),h_mR(i))*(Z_R(i)-Z_L(i));
            else
                S=S_tmp;
            end
        end
        S = S - 0.5*g*(h_mid(i)+h_mL(i))*(Z_L(i)-Z_M(i)) - 0.5*g*(h_mid(i+1)+h_mR(i))*(Z_M(i+1)-Z_R(i));
        U(:,i)=U(:,i)+d_t/d_x*(F(:,i)-F(:,i+1))+d_t/d_x*[0;S];
    end
    Time = Time+d_t
%     if Time >= 15
%         break;
%     end
end

plot_SWE