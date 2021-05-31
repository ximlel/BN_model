h = 0.5*(h_L+h_R);
u = 0.5*(u_L+u_R);

N_MAX = 2500;
d_xM = (x_max-x_min)/N_MAX;
W_exact = zeros(N_MAX,2);
load ./steady_1.mat;
W_exact(:,1)=h_E'+b_E';
W_exact(:,2)=b_E';

%plot
col = '-m';
%col = '+k';
% h=figure(1);
% set(h,'position',[100 100 1150 750]);
subplot(1,2,1);
hold on
plot(x_min:d_xM:x_max-d_xM,W_exact(:,1),'k','LineWidth',0.4);
plot(x,h+Z_L(1:N),col,'MarkerSize',6);%col,'LineWidth',1.0);
plot(x,Z_L(1:N),col,'MarkerSize',6);%col,'LineWidth',1.0);
plot(x_min:d_xM:x_max-d_xM,W_exact(:,2),'k','LineWidth',0.4);
xlabel('Position','FontWeight','bold');
ylabel('h+Z','FontWeight','bold');
ylim([0 2])
title('h+z');
set(gca,'box','on');
subplot(1,2,2);
hold on
% plot(x_min:d_xM:x_max-d_xM,W_exact(:,2),'b','LineWidth',0.4);
plot(x,u,col,'MarkerSize',6);%col,'LineWidth',1.0);
xlabel('Position','FontWeight','bold');
ylabel('Velocity','FontWeight','bold');
%ylim([0 2])
title('u');
set(gca,'box','on');