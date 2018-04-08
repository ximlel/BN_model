%test 1
clear
lo_gL_0   =0.2;
u_gL_0    =0.0;
p_gL_0    =0.3;
lo_sL_0   =1.0;
u_sL_0    =0.0;
p_sL_0    =1.0;
phi_sL_0  =0.8;
lo_gR_0   =1.0;
u_gR_0    =0.0;
p_gR_0    =1.0;
lo_sR_0   =1.0;
u_sR_0    =0.0;
p_sR_0    =1.0;
phi_sR_0  =0.3;
save test_toro1.mat

%test 2
clear
gama_s = 3.0;
gama_g = 1.35;
p0 = 3400;
lo_gL_0   =2.0;
u_gL_0    =0.0;
p_gL_0    =3.0;
lo_sL_0   =1900.0;
u_sL_0    =0.0;
p_sL_0    =10.0;
phi_sL_0  =0.2;
lo_gR_0   =1.0;
u_gR_0    =0.0;
p_gR_0    =1.0;
lo_sR_0   =1950;
u_sR_0    =0.0;
p_sR_0    =1000;
phi_sR_0  =0.9;

save test_toro2.mat

%test 3
clear
lo_gL_0   =1.0;
u_gL_0    =0.75;
p_gL_0    =1.0;
lo_sL_0   =1.0;
u_sL_0    =0.75;
p_sL_0    =1.0;
phi_sL_0  =0.8;
lo_gR_0   =0.125;
u_gR_0    =0.0;
p_gR_0    =0.1;
lo_sR_0   =0.125;
u_sR_0    =0.0;
p_sR_0    =0.1;
phi_sR_0  =0.3;
save test_toro3.mat

%test 4
clear
lo_gL_0   =1.0;
u_gL_0    =-2.0;
p_gL_0    =0.4;
lo_sL_0   =1.0;
u_sL_0    =-2.0;
p_sL_0    =0.4;
phi_sL_0  =0.8;
lo_gR_0   =1.0;
u_gR_0    =2.0;
p_gR_0    =0.4;
lo_sR_0   =1.0;
u_sR_0    =2.0;
p_sR_0    =0.4;
phi_sR_0  =0.5;
save test_toro4.mat

%test 5
clear
gama_s = 3.0;
gama_g = 1.4;
p0 = 10;
lo_gL_0   =1.4;
u_gL_0    =0.0;
p_gL_0    =1.0;
lo_sL_0   =1.4;
u_sL_0    =0.0;
p_sL_0    =2.0;
phi_sL_0  =0.6;
lo_gR_0   =1.0;
u_gR_0    =0.0;
p_gR_0    =1.0;
lo_sR_0   =1.0;
u_sR_0    =0.0;
p_sR_0    =3.0;
phi_sR_0  =0.3;
save test_toro5.mat

%test 6
clear
gama_s = 3.0;
gama_g = 1.4;
p0 = 100;
Tend=0.007;
x0=0.8;
lo_gL_0   =1.0;
u_gL_0    =-19.5975;
p_gL_0    =1000;
lo_sL_0   =1.0;
u_sL_0    =-19.5975;
p_sL_0    =1000;
phi_sL_0  =0.7;
lo_gR_0   =1.0;
u_gR_0    =-19.5975;
p_gR_0    =0.01;
lo_sR_0   =1.0;
u_sR_0    =-19.5975;
p_sR_0    =0.01;
phi_sR_0  =0.2;
save test_toro6.mat
