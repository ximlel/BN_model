global gama_s gama_g p0;
gama_s=3;
gama_g=1.4;
p0=0;

lo_gL   =1;
u_gL    =1;
p_gL    =1;
lo_sL   =2;
u_sL    =0.3;
p_sL    =5;
phi_sL  =0.8;
lo_gR   =1.17283230280172;
u_gR    =0.4705273631381315;
p_gR    =1.250058353622465;
lo_sR   =2;
u_sR    =0.3;
p_sR    =11.33028440541645;
phi_sR  =0.3;

[phi_s_out,dE_g]=E_s_correct(u_sL,lo_gL,u_gL,p_gL,lo_sL,p_sL,phi_sL,lo_gR,u_gR,p_gR,lo_sR,p_sR,phi_sR,0.3/u_sL);
