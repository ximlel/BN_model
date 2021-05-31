%continuity initial condition

Tend=300;

Z_0   = 2;
q_in  = 4.42;
h_out = 2;

% Z_0   = 0.66;
% q_in  = 1.53;
% h_out = 0.66;

% Z_0   = 0.33;
% q_in  = 0.18;
% h_out = 0.33;

x_min = 0;
x_max = 25;
N_0 = 100;

Z=@(x) (0.2-0.05*(x-10)^2)*(x>8 && x<12);

% Tend=50;
% 
% Z_0   = 3;
% q_in  = 2;
% h_out = 2;
% 
% x_min = -10;
% x_max = 10;

%Z=@(x) (1.0 + cos(pi*x/8)*(x<4)) * (x>-4);

d_x = (x_max-x_min)/N_0;
N = N_0+1;

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
W_int   =zeros(2,N+1);

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
dZ  =zeros(1,N-1);

x  =zeros(1,N);
x_M=zeros(1,N+1);
%test begin
for i=1:N+1
    x_M(i) = x_min+(i-0.5)*d_x;
    Z_M(i) = Z(x_M(i));
end
for i=1:N
    x(i)   = x_min+(i-1)*d_x;
    Z_L(i) = Z(x(i));
    Z_R(i) = Z(x(i));
end
for i=1:N-1
    dZ(i)  = (Z_L(i+1)-Z_R(i))/d_x;
end

for i=1:N
    U(:,i) = [Z_0-Z(x(i));0];
    %U(:,i) = [Z_0-Z(x(i));2.0/(Z_0-Z(x(i)))];
end