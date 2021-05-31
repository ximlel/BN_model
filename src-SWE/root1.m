N = 2500;
x = linspace(0, 25, N);
Z = @(x) (0.2-0.05*(x-10).^2).*(x>8 & x<12);
b = Z(x)';
p1 = b - 4.42^2/2/9.81/2^2 - 2;
p2 = 4.42^2/2/9.81;
p = [ones(N,1) p1 zeros(N,1) p2*ones(N,1)];
h = zeros(N,1);
for i=1:N
    r = roots(p(i,:));
%     for j=1:length(r)
%         if (imag(r(j)) == 0)
%             h(i,1)=r(j);
%         end
%     end
    h(i,1)=max(r);
end
b_E = b;
h_E = h;
save('steady_1.mat','b_E','h_E');