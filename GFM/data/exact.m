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
put_out=linear_GRP_solver_Edir_Q1D([lo_L_0,lo_R_0,0,0,0,0,u_L_0,u_R_0,0,0,0,0,-0,-0,0,0,0,0,p_L_0,p_R_0,0,0,0,0,0,0,0,0,0,0,gama_s,gama_g,1e-9,1e-9]);		
lo_L_M=put_out(15);
lo_R_M=put_out(17);
u_M   =put_out(16);
p_M   =put_out(18);
c_L_0=sqrt(gama_s*p_L_0/lo_L_0);
S_1L  =u_L_0-c_L_0;
S_1M  =u_M  -sqrt(gama_s*p_M  /lo_L_M);
S_3   =put_out(2);
x_min=0;
x_max=1;
x_0=0.2;
N=200;
d_x=(x_max-x_min)/N;
Tend  =0.0411;
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if x(i)<(x_0+S_1L*Tend)
        lo_ex(i)=lo_L_0;
        u_ex(i) =u_L_0;
        p_ex(i) =p_L_0;
    elseif x(i)<(x_0+S_1M*Tend)
        x_t=(x(i)-x_0)/Tend;
        lo_ex(i)=lo_L_0*(2/(gama_s+1)+(gama_s-1)/(gama_s+1)/c_L_0*(u_L_0-x_t))^(2/(gama_s-1));
        u_ex(i) =2/(gama_s+1)*(c_L_0+(gama_s-1)/2*u_L_0+x_t);
        p_ex(i) =p_L_0*(2/(gama_s+1)+(gama_s-1)/(gama_s+1)/c_L_0*(u_L_0-x_t))^(2*gama_s/(gama_s-1));
    elseif x(i)<(x_0+u_M*Tend)
        lo_ex(i)=lo_L_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    elseif x(i)<(x_0+S_3*Tend)
        lo_ex(i)=lo_R_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    else
        lo_ex(i)=lo_R_0;
        u_ex(i) =u_R_0;
        p_ex(i) =p_R_0;
    end
end
save('exact1.mat','lo_ex','u_ex','p_ex');

%% test11
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
put_out=linear_GRP_solver_Edir_Q1D([lo_L_0,lo_R_0,0,0,0,0,u_L_0,u_R_0,0,0,0,0,-0,-0,0,0,0,0,p_L_0,p_R_0,0,0,0,0,0,0,0,0,0,0,gama_s,gama_g,1e-9,1e-9]);		
lo_L_M=put_out(15);
lo_R_M=put_out(17);
u_M   =put_out(16);
p_M   =put_out(18);
c_L_0=sqrt(gama_s*p_L_0/lo_L_0);
S_1L  =u_L_0-c_L_0;
S_1M  =u_M  -sqrt(gama_s*p_M  /lo_L_M);
S_3   =put_out(2);
x_min=0;
x_max=1;
x_0=0.2;
N=200;
d_x=(x_max-x_min)/N;
Tend  =0.0136;
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if x(i)<(x_0+S_1L*Tend)
        lo_ex(i)=lo_L_0;
        u_ex(i) =u_L_0;
        p_ex(i) =p_L_0;
    elseif x(i)<(x_0+S_1M*Tend)
        x_t=(x(i)-x_0)/Tend;
        lo_ex(i)=lo_L_0*(2/(gama_s+1)+(gama_s-1)/(gama_s+1)/c_L_0*(u_L_0-x_t))^(2/(gama_s-1));
        u_ex(i) =2/(gama_s+1)*(c_L_0+(gama_s-1)/2*u_L_0+x_t);
        p_ex(i) =p_L_0*(2/(gama_s+1)+(gama_s-1)/(gama_s+1)/c_L_0*(u_L_0-x_t))^(2*gama_s/(gama_s-1));
    elseif x(i)<(x_0+u_M*Tend)
        lo_ex(i)=lo_L_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    elseif x(i)<(x_0+S_3*Tend)
        lo_ex(i)=lo_R_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    else
        lo_ex(i)=lo_R_0;
        u_ex(i) =u_R_0;
        p_ex(i) =p_R_0;
    end
end
save('exact11.mat','lo_ex','u_ex','p_ex');

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
put_out=linear_GRP_solver_Edir_Q1D([lo_L_0,lo_R_0,0,0,0,0,u_L_0,u_R_0,0,0,0,0,-0,-0,0,0,0,0,p_L_0,p_R_0,0,0,0,0,0,0,0,0,0,0,gama_s,gama_g,1e-9,1e-9]);		
lo_L_M=put_out(15);
lo_R_M=put_out(17);
u_M   =put_out(16);
p_M   =put_out(18);
S_1   =put_out(1);
S_3   =put_out(2);
x_min=0;
x_max=1;
x_0=0.3;
N=200;
d_x=(x_max-x_min)/N;
Tend  =0.033;
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if x(i)<(x_0+S_1*Tend)
        lo_ex(i)=lo_L_0;
        u_ex(i) =u_L_0;
        p_ex(i) =p_L_0;
    elseif x(i)<(x_0+u_M*Tend)
        lo_ex(i)=lo_L_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    elseif x(i)<(x_0+S_3*Tend)
        lo_ex(i)=lo_R_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    else
        lo_ex(i)=lo_R_0;
        u_ex(i) =u_R_0;
        p_ex(i) =p_R_0;
    end
end
save('exact3.mat','lo_ex','u_ex','p_ex');

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
put_out=linear_GRP_solver_Edir_Q1D([lo_L_0,lo_R_0,0,0,0,0,u_L_0,u_R_0,0,0,0,0,-0,-0,0,0,0,0,p_L_0,p_R_0,0,0,0,0,0,0,0,0,0,0,gama_s,gama_g,1e-9,1e-9]);		
lo_L_M=put_out(15);
lo_R_M=put_out(17);
u_M   =put_out(16);
p_M   =put_out(18);
S_1   =put_out(1);
S_3   =put_out(2);
x_min=0;
x_max=1;
x_0=0.2;
N=200;
d_x=(x_max-x_min)/N;
Tend  =0.0534;
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if x(i)<(x_0+S_1*Tend)
        lo_ex(i)=lo_L_0;
        u_ex(i) =u_L_0;
        p_ex(i) =p_L_0;
    elseif x(i)<(x_0+u_M*Tend)
        lo_ex(i)=lo_L_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    elseif x(i)<(x_0+S_3*Tend)
        lo_ex(i)=lo_R_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    else
        lo_ex(i)=lo_R_0;
        u_ex(i) =u_R_0;
        p_ex(i) =p_R_0;
    end
end
save('exact4.mat','lo_ex','u_ex','p_ex');

%% test5
clear
gama_s=1.6667;
mu_s2 =(gama_s-1)/(gama_s+1);
gama_g=1.4;
p_R_0 =1;
lo_R_0=0.1;
u_R_0 =0;
lo_M_0=250;
p_L_0 =1000;
lo_L_0=lo_M_0/(p_R_0+mu_s2*p_L_0)*(p_L_0+mu_s2*p_R_0);
u_L_0 =u_R_0+sqrt((1/lo_L_0-1/lo_M_0)*(p_R_0-p_L_0));
put_out=linear_GRP_solver_Edir_Q1D([lo_L_0,lo_R_0,0,0,0,0,u_L_0,u_R_0,0,0,0,0,-0,-0,0,0,0,0,p_L_0,p_R_0,0,0,0,0,0,0,0,0,0,0,gama_s,gama_g,1e-9,1e-9]);		
lo_L_M=put_out(15);
lo_R_M=put_out(17);
u_M   =put_out(16);
p_M   =put_out(18);
c_L_0=sqrt(gama_s*p_L_0/lo_L_0);
S_1L  =u_L_0-c_L_0;
S_1M  =u_M  -sqrt(gama_s*p_M  /lo_L_M);
S_3   =put_out(2);
x_min=0;
x_max=1;
x_0=0.1;
N=200;
d_x=(x_max-x_min)/N;
Tend  =0.1;
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if x(i)<(x_0+S_1L*Tend)
        lo_ex(i)=lo_L_0;
        u_ex(i) =u_L_0;
        p_ex(i) =p_L_0;
    elseif x(i)<(x_0+S_1M*Tend)
        x_t=(x(i)-x_0)/Tend;
        lo_ex(i)=lo_L_0*(2/(gama_s+1)+(gama_s-1)/(gama_s+1)/c_L_0*(u_L_0-x_t))^(2/(gama_s-1));
        u_ex(i) =2/(gama_s+1)*(c_L_0+(gama_s-1)/2*u_L_0+x_t);
        p_ex(i) =p_L_0*(2/(gama_s+1)+(gama_s-1)/(gama_s+1)/c_L_0*(u_L_0-x_t))^(2*gama_s/(gama_s-1));
    elseif x(i)<(x_0+u_M*Tend)
        lo_ex(i)=lo_L_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    elseif x(i)<(x_0+S_3*Tend)
        lo_ex(i)=lo_R_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    else
        lo_ex(i)=lo_R_0;
        u_ex(i) =u_R_0;
        p_ex(i) =p_R_0;
    end
end
save('exact5.mat','lo_ex','u_ex','p_ex');

%% test55
clear
gama_s=1.6667;
mu_s2 =(gama_s-1)/(gama_s+1);
gama_g=1.4;
p_R_0 =1;
lo_R_0=0.1;
u_R_0 =0;
lo_M_0=2500;
p_L_0 =10000;
lo_L_0=lo_M_0/(p_R_0+mu_s2*p_L_0)*(p_L_0+mu_s2*p_R_0);
u_L_0 =u_R_0+sqrt((1/lo_L_0-1/lo_M_0)*(p_R_0-p_L_0));
put_out=linear_GRP_solver_Edir_Q1D([lo_L_0,lo_R_0,0,0,0,0,u_L_0,u_R_0,0,0,0,0,-0,-0,0,0,0,0,p_L_0,p_R_0,0,0,0,0,0,0,0,0,0,0,gama_s,gama_g,1e-9,1e-9]);		
lo_L_M=put_out(15);
lo_R_M=put_out(17);
u_M   =put_out(16);
p_M   =put_out(18);
c_L_0=sqrt(gama_s*p_L_0/lo_L_0);
S_1L  =u_L_0-c_L_0;
S_1M  =u_M  -sqrt(gama_s*p_M  /lo_L_M);
S_3   =put_out(2);
x_min=0;
x_max=1;
x_0=0.1;
N=200;
d_x=(x_max-x_min)/N;
Tend  =0.1;
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if x(i)<(x_0+S_1L*Tend)
        lo_ex(i)=lo_L_0;
        u_ex(i) =u_L_0;
        p_ex(i) =p_L_0;
    elseif x(i)<(x_0+S_1M*Tend)
        x_t=(x(i)-x_0)/Tend;
        lo_ex(i)=lo_L_0*(2/(gama_s+1)+(gama_s-1)/(gama_s+1)/c_L_0*(u_L_0-x_t))^(2/(gama_s-1));
        u_ex(i) =2/(gama_s+1)*(c_L_0+(gama_s-1)/2*u_L_0+x_t);
        p_ex(i) =p_L_0*(2/(gama_s+1)+(gama_s-1)/(gama_s+1)/c_L_0*(u_L_0-x_t))^(2*gama_s/(gama_s-1));
    elseif x(i)<(x_0+u_M*Tend)
        lo_ex(i)=lo_L_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    elseif x(i)<(x_0+S_3*Tend)
        lo_ex(i)=lo_R_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    else
        lo_ex(i)=lo_R_0;
        u_ex(i) =u_R_0;
        p_ex(i) =p_R_0;
    end
end
save('exact55.mat','lo_ex','u_ex','p_ex');

%% test5555 %same gama
clear
gama_s=1.4;
mu_s2 =(gama_s-1)/(gama_s+1);
gama_g=1.4;
p_R_0 =1;
lo_R_0=0.1;
u_R_0 =0;
lo_M_0=2500;
p_L_0 =10000;
lo_L_0=lo_M_0/(p_R_0+mu_s2*p_L_0)*(p_L_0+mu_s2*p_R_0);
u_L_0 =u_R_0+sqrt((1/lo_L_0-1/lo_M_0)*(p_R_0-p_L_0));
put_out=linear_GRP_solver_Edir_Q1D([lo_L_0,lo_R_0,0,0,0,0,u_L_0,u_R_0,0,0,0,0,-0,-0,0,0,0,0,p_L_0,p_R_0,0,0,0,0,0,0,0,0,0,0,gama_s,gama_g,1e-9,1e-9]);		
lo_L_M=put_out(15);
lo_R_M=put_out(17);
u_M   =put_out(16);
p_M   =put_out(18);
c_L_0=sqrt(gama_s*p_L_0/lo_L_0);
S_1L  =u_L_0-c_L_0;
S_1M  =u_M  -sqrt(gama_s*p_M  /lo_L_M);
S_3   =put_out(2);
x_min=0;
x_max=1;
x_0=0.1;
N=200;
d_x=(x_max-x_min)/N;
Tend  =0.1;
for i=1:N
    x(i)=x_min+(i-0.5)*d_x;
    if x(i)<(x_0+S_1L*Tend)
        lo_ex(i)=lo_L_0;
        u_ex(i) =u_L_0;
        p_ex(i) =p_L_0;
    elseif x(i)<(x_0+S_1M*Tend)
        x_t=(x(i)-x_0)/Tend;
        lo_ex(i)=lo_L_0*(2/(gama_s+1)+(gama_s-1)/(gama_s+1)/c_L_0*(u_L_0-x_t))^(2/(gama_s-1));
        u_ex(i) =2/(gama_s+1)*(c_L_0+(gama_s-1)/2*u_L_0+x_t);
        p_ex(i) =p_L_0*(2/(gama_s+1)+(gama_s-1)/(gama_s+1)/c_L_0*(u_L_0-x_t))^(2*gama_s/(gama_s-1));
    elseif x(i)<(x_0+u_M*Tend)
        lo_ex(i)=lo_L_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    elseif x(i)<(x_0+S_3*Tend)
        lo_ex(i)=lo_R_M;
        u_ex(i) =u_M;
        p_ex(i) =p_M;
    else
        lo_ex(i)=lo_R_0;
        u_ex(i) =u_R_0;
        p_ex(i) =p_R_0;
    end
end
save('exact5555.mat','lo_ex','u_ex','p_ex');
