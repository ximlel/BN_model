clear
clc

N = 81;
x = linspace(2, 10, N);
Z = @(x) (1.0 + cos(pi*x/8).*(x<4)).*(x>-4);
b = Z(x)';

p1 = 2*b*9.8 - 6*9.8 - 1;
p2 = 4*9.8;
p = [ones(N,1) zeros(N,1) p1 p2*ones(N,1)];
h = zeros(N,1);
for i=1:N
    r = 2./roots(p(i,:))
    h(i,1)=max(r);
end
b_E = b;
h_E = h;
save('steady_discon.mat','b_E','h_E');
%plot(x,h+b);