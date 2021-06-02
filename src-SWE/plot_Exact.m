N_MAX = 2500;
d_xM = (x_max-x_min)/N_MAX;
W_exact = zeros(N_MAX+1,2);
load ./steady_3.mat;
W_exact(:,1)=h_E'+b_E';
W_exact(:,2)=b_E';

%plot
col = '-m';
hold on
plot(x_min:d_xM:x_max,W_exact(:,1),'k','LineWidth',0.4);
ylim([-0.1 0.45])
plot(x_min:d_xM:x_max,W_exact(:,2),'k','LineWidth',0.4);
xlabel('Position','FontWeight','bold');
ylabel('h+Z','FontWeight','bold');