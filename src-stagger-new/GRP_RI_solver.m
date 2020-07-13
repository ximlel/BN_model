function [F,W_int]=GRP_RI_solver(rho_gL,rho_gR,u_gL,u_gR,p_gL,p_gR,lo_sL,lo_sR,u_sL,u_sR,p_sL,p_sR,phi_sL,phi_sR,D_L,D_R,d_t)
    global ep;
    global gama_s gama_g;
    F=zeros(3,1);
    put_out_g=linear_GRP_solver_Edir_Q1D([rho_gL,rho_gR,0.0,0.0,0,0,u_gL,u_gR,0.0,0.0,0,0,-0,-0,0,0,0,0,p_gL,p_gR,0.0,0.0,0,0,phi_sL,phi_sR,0.0,0.0,0,0,gama_g,gama_g,ep,ep]);		
    rho_g = put_out_g(9);
    u_g   = put_out_g(10);
    p_g   = put_out_g(12);
    put_out_s=linear_GRP_solver_Edir_Q1D([rho_sL,rho_sR,0.0,0.0,0,0,u_sL,u_sR,0.0,0.0,0,0,-0,-0,0,0,0,0,p_sL,p_sR,0.0,0.0,0,0,phi_sL,phi_sR,0.0,0.0,0,0,gama_s,gama_s,ep,ep]);		
    rho_s = put_out_s(9);
    u_s   = put_out_s(10);
    p_s   = put_out_s(12);
    phi_s = put_out_s(14);    p_s   = put_out_g(1
    c_s = sqrt(gamma_s * p_s / rho_s);
    c_g = sqrt(gamma_g * p_g / rho_g);
    Lambda = diag([u_s, u_s-c_s, u_s, u_s+c_s, u_g-c_g, u_g, u_g+c_g]);
    Lambda_v_p = max(Lambda,0.0);
    Lambda_v_m = min(Lambda,0.0);
    GAMMA_g = gama_g-1.0;
    V = u_g-u_s;
    T_g = rho_g^GAMMA_g/GAMMA_g;    
    R = zeros(7,7);
    R(1,1) = 1.0;
    R(2,2) = 1.0/c_s;
    R(2,3) = 1.0;
    R(2,4) = 1.0/c_s;
    R(3,2) =-1.0/rho_s;
    R(3,4) = 1.0/rho_s;	
    R(4,2) = phi_s*c_s + 2.0*phi_g*rho_g*V/rho_s;
    R(4,4) = phi_s*c_s - 2.0*phi_g*rho_g*V/rho_s;
    R(4,5) = u_g-c_g-u_s;
    R(4,6) =-phi_g*rho_g*T_g*GAMMA_g*V*V/c_g/c_g;
    R(4,7) = u_g+c_g-u_s;
    R(5,2) = phi_g*rho_g/rho_s;
    R(5,4) =-phi_g*rho_g/rho_s;
    R(5,5) = 1.0;
    R(5,6) =-phi_g*rho_g*T_g*GAMMA_g*V/c_g/c_g;
    R(5,7) = 1.0;
    R(6,2) = V/rho_s;
    R(6,4) =-V/rho_s;
    R(6,5) =-c_g/phi_g/rho_g;
    R(6,6) = T_g;
    R(6,7) = c_g/phi_g/rho_g;
    R(7,6) = 1.0;	
    
    F(1) = phi*rho_mid*u_mid;
    F(2) = F(1)*u_mid + phi*p_mid;
    F(3) = (gama/(gama-1.0))*p_mid/rho_mid + 0.5*u_mid*u_mid;
    F(3) = F(1)*F(3);
    W_int(1) = put_out(9)  + d_t*put_out(3);	
    W_int(2) = put_out(10) + d_t*put_out(4);			
    W_int(3) = put_out(12) + d_t*put_out(6);   
    W_int(4) = put_out(14) + d_t*put_out(8);
%     phi  = F(1)*phi;
end

