syms S_gL S_gR S_sL S_sR;
syms lo_gL lo_gR p_gL p_gR u_gL u_gR lo_sL lo_sR p_sL p_sR u_sL u_sR phi_sL phi_sR phi_gL phi_gR;
syms p_s_srL p_s_srR p_g_srL p_g_srR;
syms gama_g;
u_s_srL  = u_sL + (p_s_srL-p_sL)/lo_sL/(S_sL-u_sL);
u_s_srR  = u_sR + (p_s_srR-p_sR)/lo_sR/(S_sR-u_sR);
u_g_srL  = u_gL + (p_g_srL-p_gL)/lo_gL/(S_gL-u_gL);
u_g_srR  = u_gR + (p_g_srR-p_gR)/lo_gR/(S_gR-u_gR);
lo_s_srL = lo_sL^2*(S_sL-u_sL)^2 / (lo_sL*(S_sL-u_sL)^2-p_s_srL+p_sL);
lo_s_srR = lo_sR^2*(S_sR-u_sR)^2 / (lo_sR*(S_sR-u_sR)^2-p_s_srR+p_sR);
lo_g_srL = lo_gL^2*(S_gL-u_gL)^2 / (lo_gL*(S_gL-u_gL)^2-p_g_srL+p_gL);
lo_g_srR = lo_gR^2*(S_gR-u_gR)^2 / (lo_gR*(S_gR-u_gR)^2-p_g_srR+p_gR);
RHL(1) = u_s_srR-u_s_srL;
RHL(2) = phi_gR*(p_g_srR/p_g_srL)^(1/gama_g)*(u_g_srR-u_s_srR)-phi_gL*(u_g_srL-u_s_srL);
RHL(3) = phi_gL*lo_g_srL*(u_g_srL-u_s_srL)*(u_g_srR-u_g_srL)+phi_sR*p_s_srR+phi_gR*p_g_srR-phi_sL*p_s_srL-phi_gL*p_g_srL;
RHL(4) = gama_g/(gama_g-1)*p_g_srR/lo_g_srL*(p_g_srL/p_g_srR)^(1/gama_g)+0.5*(u_g_srR-u_s_srR)^2-gama_g/(gama_g-1)*p_g_srL/lo_g_srL-0.5*(u_g_srL-u_s_srL)^2;
dfun_L=[diff(RHL,'p_s_srL');diff(RHL,'p_s_srR');diff(RHL,'p_g_srL');diff(RHL,'p_g_srR')]
% fun_L  = matlabFunction(RHL,'Vars',[p_s_srL p_s_srR p_g_srL p_g_srR]);
% dfun_L = matlabFunction([diff(RHL,'p_s_srL');diff(RHL,'p_s_srR');diff(RHL,'p_g_srL');diff(RHL,'p_g_srR')],'Vars',[p_s_srL p_s_srR p_g_srL p_g_srR]);
RHR(1) = u_s_srR-u_s_srL;
RHR(2) = phi_gR*(u_g_srR-u_s_srR)-phi_gL*(p_g_srL/p_g_srR)^(1/gama_g)*(u_g_srL-u_s_srL);
RHR(3) = phi_gR*lo_g_srR*(u_g_srR-u_s_srR)*(u_g_srR-u_g_srL)+phi_sR*p_s_srR+phi_gR*p_g_srR-phi_sL*p_s_srL-phi_gL*p_g_srL;
RHR(4) = gama_g/(gama_g-1)*p_g_srR/lo_g_srR+0.5*(u_g_srR-u_s_srR)^2-gama_g/(gama_g-1)*p_g_srL/lo_g_srR*(p_g_srR/p_g_srL)^(1/gama_g)-0.5*(u_g_srL-u_s_srL)^2;
dfun_R=[diff(RHR,'p_s_srL');diff(RHR,'p_s_srR');diff(RHR,'p_g_srL');diff(RHR,'p_g_srR')]^-1
% fun_R  = matlabFunction(RHR,'Vars',[p_s_srL p_s_srR p_g_srL p_g_srR]);
% dfun_R = matlabFunction([diff(RHR,'p_s_srL');diff(RHR,'p_s_srR');diff(RHR,'p_g_srL');diff(RHR,'p_g_srR')],'Vars',[p_s_srL p_s_srR p_g_srL p_g_srR]);
