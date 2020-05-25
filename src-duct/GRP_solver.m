function [F,W_int]=GRP_solver(lo_L,lo_R,dlo_L,dlo_R,u_L,u_R,du_L,du_R,p_L,p_R,dp_L,dp_R,phi_L,phi_R,dphi_L,dphi_R,gama,d_t)
    global ep;
    F=zeros(3,1);
    put_out=linear_GRP_solver_Edir_Q1D([lo_L,lo_R,dlo_L,dlo_R,0,0,u_L,u_R,du_L,du_R,0,0,-0,-0,0,0,0,0,p_L,p_R,dp_L,dp_R,0,0,phi_L,phi_R,dphi_L,dphi_R,0,0,gama,gama,ep,ep]);		
    rho_mid = put_out(9)  + 0.5*d_t*put_out(3);	
    u_mid   = put_out(10) + 0.5*d_t*put_out(4);			
    p_mid   = put_out(12) + 0.5*d_t*put_out(6);   
    phi  = put_out(14) + 0.5*d_t*put_out(8);
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

