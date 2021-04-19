load EXACT_duct2.mat;
Y1=eta_E(7250);
Y2=eta_E(7380);
MAX = eta_E(7500);
MIN = eta_E(6300);
for i=6300:7500
   eta_E(i) = (Y2-Y1)/(7380-7250)*(i-7250)+Y1;
   eta_E(i) = min(MAX,eta_E(i));
   eta_E(i) = max(MIN,eta_E(i));
end
lo_0 = lo_E(3745);
for i=3334:3760
   lo_E(i) = (90-lo_0)/(3334-3760)*(i-3760)+lo_0;
end
p_0 = p_E(3745);
for i=3334:3760
   p_E(i) = (1.35e8-p_0)/(3334-3760)*(i-3760)+p_0;
end
u_0 = u_E(3745);
for i=3334:3760
   u_E(i) = (1400-u_0)/(3334-3760)*(i-3760)+u_0;
end
save('EXACT_duct1-2.mat','eta_E','p_E','lo_E','u_E');
