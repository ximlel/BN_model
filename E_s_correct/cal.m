global gama_s gama_g p0;
gama_s=1.4;
gama_g=1.4;
p0=0;

lo_gL   =1;
u_gL    =2;
p_gL    =1;
lo_sL   =2;
u_sL    =0.3;
p_sL    =5;
phi_sL  =0.8;
lo_gR   =0.1941934235006083;
u_gR    =2.801188129642115;
p_gR    =0.1008157360849781;
lo_sR   =2;
u_sR    =0.3;
p_sR    =12.85675006887399;
phi_sR  =0.3;

[E_sum1,E_sum2]=E_s_correct_fin(u_sL,lo_gL,u_gL,p_gL,lo_sL,p_sL,phi_sL,lo_gR,u_gR,p_gR,lo_sR,p_sR,phi_sR,300/2400)