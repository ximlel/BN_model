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
x0s   =0.2;
x0    =0.2;
N_T   =200;
save test1.mat
%% test2
clear
gama_s=1.6667;
mu_s2 =(gama_s-1)/(gama_s+1);
gama_g=1.4;
p_R_0 =1;
lo_R_0=0.1;
u_R_0 =0;
lo_M_0=1;
p_L_0 =1000;
lo_L_0=lo_M_0/(p_R_0+mu_s2*p_L_0)*(p_L_0+mu_s2*p_R_0);
u_L_0 =u_R_0+sqrt((1/lo_L_0-1/lo_M_0)*(p_R_0-p_L_0));
x0s   =0.2;
x0    =0.2;
N_T   =200;
save test2.mat
%% test3
clear
gama_s=1.6667;
mu_s2 =(gama_s-1)/(gama_s+1);
gama_g=1.4;
p_R_0 =1;
lo_R_0=1;
u_R_0 =0;
lo_M_0=0.1;
p_L_0 =100;
lo_L_0=lo_M_0/(p_R_0+mu_s2*p_L_0)*(p_L_0+mu_s2*p_R_0);
u_L_0 =u_R_0+sqrt((1/lo_L_0-1/lo_M_0)*(p_R_0-p_L_0));
x0s   =0.2;
x0    =0.2;
N_T   =350;
save test3.mat
%% test4
clear
gama_s=5/3;
mu_s2 =(gama_s-1)/(gama_s+1);
gama_g=1.2;
p_R_0 =1;
lo_R_0=1;
u_R_0 =0;
lo_M_0=0.82369077;
p_L_0 =100;
lo_L_0=lo_M_0/(p_R_0+mu_s2*p_L_0)*(p_L_0+mu_s2*p_R_0);
u_L_0 =u_R_0+sqrt((1/lo_L_0-1/lo_M_0)*(p_R_0-p_L_0));
x0s   =0.2;
x0    =0.2;
N_T   =200;
save test4.mat