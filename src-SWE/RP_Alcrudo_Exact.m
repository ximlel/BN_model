clear
clc

g = 9.8;
N = 2501;
t = 0.5;
Z = @(x) x>=0;

x = linspace(-10, 10, N);
b = Z(x)';
h_L = 4;
u_L = -10;
h_R = 1;
u_R = -6;
xi_1 = -16.26099;
xi_2 = -9.2307;
s = -2.931634;
h = zeros(N,1);
m = zeros(N,1);
for i = 1:N
    if x(i) < xi_1*t
        h(i) = h_L;
        m(i) = h_L*u_L;
    elseif x(i) >= xi_1*t && x(i) < xi_2*t
        h(i) = (-1/3/sqrt(g)*(x(i)/t-xi_1)+2)^2;
        m(i) = (2/3*(x(i)/t-xi_1)-10)*h(i);  
    elseif x(i) >= xi_2*t && x(i) < s*t
        h(i) = 1.566049;
        m(i) = -8.320635;       
    elseif x(i) >= s*t    && x(i) < 0
        h(i) = 0.774464;
        m(i) = h_R*u_R; 
    elseif x(i) >= 0
        h(i) = h_R;
        m(i) = h_R*u_R; 
    end
end
b_E = b;
h_E = h;
u_E = m./h;
save('RP_A1.mat','b_E','h_E','u_E');

x = linspace(-10, 10, N);
b = Z(x)';
h_L = 4;
u_L =-10;
h_R = 2;
u_R = 0;
xi_1 = -16.26099;
xi_2 = -6.91461;
xi_3 =  4.427189;
s = -2.037139;
h = zeros(N,1);
m = zeros(N,1);
for i = 1:N
    if x(i) < xi_1*t
        h(i) = h_L;
        m(i) = h_L*u_L;
    elseif x(i) >= xi_1*t && x(i) < xi_2*t
        h(i) =  (-1/3/sqrt(g)*(x(i)/t-xi_1)+2)^2;
        m(i) =  (2/3*(x(i)/t-xi_1)-10)*h(i);  
    elseif x(i) >= xi_2*t && x(i) < s*t
        h(i) =  1.009629;
        m(i) =  -3.805371;
    elseif x(i) >= s*t && x(i) < 0
        h(i) =  0.429476;
        m(i) =  -2.623519;
    elseif x(i) >= 0 && x(i) <= xi_3*t
        h(i) =  (1/3/sqrt(g)*x(i)/t+0.942809)^2;
        m(i) =  (2/3*x(i)/t-2.951459)*h(i);
    else 
        h(i) = h_R;
        m(i) = h_R*u_R;    
    end
end
b_E = b;
h_E = h;
u_E = m./h;
save('RP_A2.mat','b_E','h_E','u_E');
