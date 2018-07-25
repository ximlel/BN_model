function [F,u_MID,W_int_L,W_int_R]=GRP_solver_LAG(lo_L,lo_R,dlo_L,dlo_R,u_L,u_R,du_L,du_R,p_L,p_R,dp_L,dp_R,gamaL,gamaR,d_t)
    global ep;
    F=zeros(3,1);
    put_out=linear_GRP_solver_Edir_star([lo_L,lo_R,dlo_L,dlo_R,0,0,u_L,u_R,du_L,du_R,0,0,-0,-0,0,0,0,0,p_L,p_R,dp_L,dp_R,0,0,-0,-0,0,0,0,0,gamaL,gamaR,ep,ep]);		
    u_MID = put_out(16) + 0.5*d_t*put_out(27);			
    p_MID = put_out(18) + 0.5*d_t*put_out(28);   
    F(1) = 0;
    F(2) = p_MID;
    F(3) = p_MID*u_MID;
    W_int_L(1) = put_out(15) + d_t*put_out(29);	
    W_int_L(2) = put_out(16) + d_t*put_out(27);	
    W_int_L(3) = put_out(18) + d_t*put_out(28);   
    W_int_R(1) = put_out(17) + d_t*put_out(30);	
    W_int_R(2) = put_out(16) + d_t*put_out(27);	
    W_int_R(3) = put_out(18) + d_t*put_out(28);   
end