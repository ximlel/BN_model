h = 0.5*(h_L+h_R);
u = 0.5*(u_L+u_R);

N_MAX = 2501;
W_exact = zeros(N_MAX,2);
load ./steady_2.mat;
W_exact(:,1)=h_E'+b_E';
W_exact(:,2)=b_E';

%plot
col = '.m';
% col = '+k';
% h=figure(1);
% set(h,'position',[100 100 1150 750]);
subplot(1,2,1);
hold on
plot(linspace(x_min, x_max, N_MAX),W_exact(:,1),'k','LineWidth',0.4);
plot(x,h+Z_L(1:N),col,'MarkerSize',6);%col,'LineWidth',1.0);
xlim([0,25])
ylim([-0.1 1.1])
plot(x,Z_L(1:N),col,'MarkerSize',6);%col,'LineWidth',1.0);
plot(linspace(x_min, x_max, N_MAX),W_exact(:,2),'k','LineWidth',0.4);
xlabel('Position','FontWeight','bold');
ylabel('h+Z','FontWeight','bold');

title('h+z');
set(gca,'box','on');
subplot(1,2,2);
hold on
plot(x,h.*u,col,'MarkerSize',6);%col,'LineWidth',1.0);
plot(linspace(x_min, x_max, N_MAX),W_exact(:,2),'b','LineWidth',0.4);
xlabel('Position','FontWeight','bold');
ylabel('Velocity','FontWeight','bold');
xlim([0,25])
ylim([1 2])
title('hu');
set(gca,'box','on');
