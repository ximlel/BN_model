clear
clc

g = 1.0;
N = 2501;
t = 0.5;
Z = @(x) x>0;

x = linspace(-1, 2, N);
b = Z(x)';
h_L = 4.7;
u_L = 0.651;
h_R = 0.8;
u_R = 2.048;
c_L = sqrt(g*h_L);
c_R = sqrt(g*h_R);
c_1 = 2.051835;
u_1 = 0.883145;
c_2 = sqrt(1.2);
u_2 = 3*sqrt(2.4)-2*sqrt(1.2);
c = zeros(N,1);
u = zeros(N,1);
for i = 1:N
    if x(i) < (u_L-c_L)*t
        c(i) = c_L;
        u(i) = u_L;
    elseif x(i) >= (u_L-c_L)*t && x(i) < (u_1-c_1)*t
        c(i) = 1/3*(u_L+2*c_L-  x(i)/t);
        u(i) = 1/3*(u_L+2*c_L+2*x(i)/t);  
    elseif x(i) >= (u_1-c_1)*t && x(i) < 0
        c(i) = c_1;
        u(i) = u_1;       
    elseif x(i) >= 0           && x(i) < (u_2-c_2)*t
        c(i) = 1/3*(u_2+2*c_2-  x(i)/t);
        u(i) = 1/3*(u_2+2*c_2+2*x(i)/t); 
    elseif x(i) >= (u_2-c_2)*t && x(i) <= (c_R^2*u_R-c_2^2*u_2)/(c_R^2-c_2^2)*t
        c(i) = c_2;
        u(i) = u_2;
    else 
        c(i) = c_R;
        u(i) = u_R;    
    end
end
b_E = b;
h_E = c.^2/g;
u_E = u;
save('RP_1.mat','b_E','h_E','u_E');

x = linspace(-3, 1, N);
b = Z(x)';
h_L = 4;
u_L =-2.894;
h_R = 2.7;
u_R =-3;
c_L = sqrt(g*h_L);
c_R = sqrt(g*h_R);
c_1 = 1.897232;
u_1 =-2.688463;
c_2 = 1.540343;
u_2 =-3.413897;
c = zeros(N,1);
u = zeros(N,1);
for i = 1:N
    if x(i) < (u_L-c_L)*t
        c(i) = c_L;
        u(i) = u_L;
    elseif x(i) >= (u_L-c_L)*t && x(i) < (u_1-c_1)*t
        c(i) = 1/3*(u_L+2*c_L-  x(i)/t);
        u(i) = 1/3*(u_L+2*c_L+2*x(i)/t);  
    elseif x(i) >= (u_1-c_1)*t && x(i) < (c_2^2*u_2-c_1^2*u_1)/(c_2^2-c_1^2)*t
        c(i) = c_1;
        u(i) = u_1;
    elseif x(i) >= (c_2^2*u_2-c_1^2*u_1)/(c_2^2-c_1^2)*t && x(i) <= 0
        c(i) = c_2;
        u(i) = u_2;
    else 
        c(i) = c_R;
        u(i) = u_R;    
    end
end
b_E = b;
h_E = c.^2/g;
u_E = u;
save('RP_2.mat','b_E','h_E','u_E');