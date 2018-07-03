%initial condition
%% test1
clear
gama_s=1.6667;
mu_s2 =(gama_s-1)/(gama_s+1);
gama_g=1.4;
p_R_0 =1;
lo_R_0=0.1;
u_R_0 =0;
lo_M_0=1;
p_L_0 =100;
lo_L_0=lo_M_0/(p_R_0+mu_s2*p_L_0)*(p_L_0+mu_s2*p_R_0);
u_L_0 =u_R_0+sqrt((1/lo_L_0-1/lo_M_0)*(p_R_0-p_L_0));
lo_L_M=1.63335093183976;
lo_R_M=0.482920839564190;
u_M   =13.4734944938969;
p_M   =23.8943038792879;
x_0   =0.2;
S_3   = u_R_0+sqrt(gama_g*p_R_0/lo_R_0)*((gama_g+1)/2/gama_g*p_M/p_R_0+(gama_g-1)/2/gama_g);
S_1   = u_L_0-sqrt(gama_s*p_L_0/lo_L_0)*((gama_s+1)/2/gama_s*p_M/p_L_0+(gama_s-1)/2/gama_s);
save test1.mat
