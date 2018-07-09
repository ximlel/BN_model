function [p_IL,u_IL,lo_IL,dp_IL,du_IL,dlo_IL,p_IR,u_IR,lo_IR,dp_IR,du_IR,dlo_IR,u_mat]=ghost_cal(lo_L,u_L,p_L,dlo_L,du_L,dp_L,gama_L,lo_R,u_R,p_R,dlo_R,du_R,dp_R,gama_R)
    global ep;
    put_out=linear_GRP_solver_Edir_star([lo_L,lo_R,dlo_L,dlo_R,0,0,u_L,u_R,du_L,du_R,0,0,-0,-0,0,0,0,0,p_L,p_R,dp_L,dp_R,0,0,-0,-0,0,0,0,0,gama_L,gama_R,ep,ep]);
    p_IL =put_out(18);
    u_IL =put_out(16);
    p_IR =put_out(18);
    u_IR =put_out(16);
    lo_IR=put_out(15);
    lo_IL=put_out(17);
    dp_IL=put_out(24);
    du_IL=put_out(22);
    dp_IR=put_out(26);
    du_IR=put_out(25);
    dlo_IL=put_out(21);
    dlo_IR=put_out(23);
    u_mat=put_out(27);
end
