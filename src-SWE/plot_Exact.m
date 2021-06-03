N_MAX = 2500;
d_xM = (x_max-x_min)/N_MAX;
W_exact = zeros(N_MAX+1,3);
load ./RP_2.mat;
W_exact(:,1)=h_E'+b_E';
W_exact(:,2)=b_E';
W_exact(:,3)=u_E';

%plot
col = '-m';
subplot(1,2,1);
hold on
plot(x_min:d_xM:x_max,W_exact(:,1),'k','LineWidth',0.5);
ylim([-0.1 5])
plot(x_min:d_xM:x_max,W_exact(:,2),'k','LineWidth',0.5);
xlabel('Position','FontWeight','bold');
ylabel('h+Z','FontWeight','bold');

subplot(1,2,2);
hold on
plot(x_min:d_xM:x_max,W_exact(:,3),'k','LineWidth',0.5);
ylim([-3.6 -2.4])
xlabel('Position','FontWeight','bold');
ylabel('Velocity','FontWeight','bold');
