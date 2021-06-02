clear
clc

g = 9.81;
N = 2501;
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

b_E = b;
h_E = h;
save('RP_1.mat','b_E','h_E');

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

b_E = b;
h_E = h;
save('RP_2.mat','b_E','h_E');