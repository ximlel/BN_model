%discontinuity initial condition

ep = 1e-6;

% Tend=0.5;
% g = 1.0;

% h_L_0  =2;
% u_L_0  =0;
% h_R_0  =1;
% u_R_0  =0;
% Z_L_0  =0;
% Z_R_0  =1;
% x_min = -1;
% x_max = 2;
% N_0 = 300;

% h_L_0  =4.7;
% u_L_0  =0.651;
% h_R_0  =0.8;
% u_R_0  =2.048;
% Z_L_0  =0;
% Z_R_0  =1;
% x_min = -1;
% x_max = 2;
% N_0 = 300;

% h_L_0  =4;
% u_L_0  =-2.894;
% h_R_0  =2.7;
% u_R_0  =-3;
% Z_L_0  =0;
% Z_R_0  =1;
% x_min = -3;
% x_max = 1;
% N_0 = 400;

Tend = 0.5;
g = 9.8;

h_L_0  =4;
u_L_0  =-10;
h_R_0  =1;
u_R_0  =-6;
Z_L_0  =0;
Z_R_0  =1;
x_min = -10;
x_max = 10;
N_0 = 100;

% h_L_0  =4;
% u_L_0  =-10;
% h_R_0  =2;
% u_R_0  =0;
% Z_L_0  =0;
% Z_R_0  =1;
% x_min = -10;
% x_max = 10;
% N_0 = 100;

d_x = (x_max-x_min)/N_0;
N = N_0+1;

Z = @(x) Z_L_0*(x <= 0.0) + Z_R_0*(x > 0.0);

Z_L=zeros(1,N);
Z_R=zeros(1,N);
Z_M=zeros(1,N+1);

U=zeros(2,N);
F=zeros(2,N+1);
h_L_int =zeros(1,N+1);
h_R_int =zeros(1,N+1);
u_L_int =zeros(1,N+1);
u_R_int =zeros(1,N+1);
dh_L_int=zeros(1,N+1);
dh_R_int=zeros(1,N+1);
du_L_int=zeros(1,N+1);
du_R_int=zeros(1,N+1);
h_mid   =zeros(1,N+1);
W_int   =zeros(4,N+1);

Fr_L=zeros(1,N);
Fr_R=zeros(1,N);
h_L =zeros(1,N);
h_R =zeros(1,N);
u_L =zeros(1,N);
u_R =zeros(1,N);
a_L =zeros(1,N);
a_R =zeros(1,N);
H_t =zeros(1,N);
du  =zeros(1,N);
dq  =zeros(1,N);
dZ  =zeros(1,N+1);

x  =zeros(1,N);
x_M=zeros(1,N+1);
%test begin
for i=1:N+1
    x_M(i) = x_min+(i-0.5)*d_x;
    Z_M(i) = Z(x_M(i));
end
for i=1:N
    x(i)   = x_min+(i-1)*d_x;
    Z_L(i) = Z(x(i)-ep);
    Z_R(i) = Z(x(i)+ep);
end
dZ(1) = 0;
for i=2:N
    dZ(i)  = (Z_L(i)-Z_R(i-1))/d_x;
end
dZ(N+1) = 0;

U_L_0=[h_L_0;h_L_0*u_L_0];
U_R_0=[h_R_0;h_R_0*u_R_0];
for i=1:N
    if x(i) < -ep
        U(:,i) = U_L_0;
    elseif x(i) > ep
        U(:,i) = U_R_0;
    else
        U(:,i) = 0.5*(U_L_0+U_R_0);
    end
end