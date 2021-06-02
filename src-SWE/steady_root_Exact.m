clear
clc

N = 2501;
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

p1 = b - 0.2 - 1.5*(1.53^2/9.81)^(1/3);
p2 = 1.53^2/2/9.81;
p = [ones(N,1) p1 zeros(N,1) p2*ones(N,1)];
h = zeros(N,1);
for i=1:N
    r = roots(p(i,:));
    if x(i) <= 10
        h(i,1) = max(r);
    else
        [rr,pos] = sort(r);
        h(i,1) = rr(2);
    end
end
b_E = b;
h_E = h;
save('steady_2.mat','b_E','h_E');

x_s = 11.65;
pL1 = b - 0.2 - 1.5*(0.18^2/9.81)^(1/3);
pL2 = 0.18^2/2/9.81;
pL = [ones(N,1) pL1 zeros(N,1) pL2*ones(N,1)];
pR1 = b - 0.18^2/2/9.81/0.33^2 - 0.33;
pR2 = 0.18^2/2/9.81;
pR = [ones(N,1) pR1 zeros(N,1) pR2*ones(N,1)];
h = zeros(N,1);
for i=1:N
    rL = roots(pL(i,:));
    rR = roots(pR(i,:));
    if x(i) < 10
        h(i,1) = max(rL);
    elseif x(i) == 10
        h(i,1) = max(real(rL));
    else
        if x(i) <= x_s
            [rrL,pos] = sort(rL);
            h(i,1) = rrL(2);
%             diff = trans_flow_diff(rrL(2),max(rR))
        else
            h(i,1) = max(rR);
        end
    end
end
b_E = b;
h_E = h;
save('steady_3.mat','b_E','h_E');