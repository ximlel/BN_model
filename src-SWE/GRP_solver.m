function [F,W_int]=GRP_solver(h_L,h_R,dh_L,dh_R,u_L,u_R,du_L,du_R,Z,dZ_L,dZ_R,g,d_t)
    global ep;
    F=zeros(2,1);
    put_out=linear_GRP_solver_SWE([h_L,h_R,dh_L,dh_R,u_L,u_R,du_L,du_R,Z,dZ_L,dZ_R,g,ep]);		
    h_mid = put_out(1) + 0.5*d_t*put_out(3);	
    u_mid = put_out(2) + 0.5*d_t*put_out(4);			
    F(1) = h_mid*u_mid;
    F(2) = h_mid*u_mid^2 + 0.5*g*h_mid*h_mid;
    W_int(1) = put_out(1) + d_t*put_out(3);	
    W_int(2) = put_out(2) + d_t*put_out(4);			
end