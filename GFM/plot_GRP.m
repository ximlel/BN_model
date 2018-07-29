x_min=0;
x_max=2;
N=4000;
d_x=(x_max-x_min)/N;
EXACT_LOCAT='./data/exact1.mat';
W_exact = zeros(N,4);
% W_exact(:,1)=lo';
% W_exact(:,2)=u';
% W_exact(:,3)=p';
% W_exact(:,4)=gama';
load(EXACT_LOCAT);
for i=1:N
     W_exact(i,1) = lo_ex(ceil(i/(N/200)));
     W_exact(i,2) = u_ex(ceil(i/(N/200)));
     W_exact(i,3) = p_ex(ceil(i/(N/200)));
end
% W_exact(:,1)=log(W_exact(:,1));
% W_exact(:,3)=log(W_exact(:,3));
% lo=log(lo);
% p=log(p);

%plot
col = '.b';
figure(1);
subplot(2,2,1);
hold on
plot(x_min+0.5*d_x:d_x:x_max-0.5*d_x,W_exact(:,1),'k','LineWidth',1.0);
plot(x,lo,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Density','FontWeight','bold');
subplot(2,2,2);
hold on
plot(x_min+0.5*d_x:d_x:x_max-0.5*d_x,W_exact(:,2),'k','LineWidth',1.0);
plot(x,u,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Velocity','FontWeight','bold');
subplot(2,2,3);
hold on
plot(x_min+0.5*d_x:d_x:x_max-0.5*d_x,W_exact(:,3),'k','LineWidth',1.0);
plot(x,p,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Pressure','FontWeight','bold');
subplot(2,2,4);
hold on
% plot(x_min+0.5*d_x:d_x:x_max-0.5*d_x,W_exact(:,4),'k','LineWidth',1.0);
plot(x,gama,col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Gamma','FontWeight','bold');
hold off