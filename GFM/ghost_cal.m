function [p_IL,u_IL,lo_IL,p_IR,u_IR,lo_IR]=ghost_cal(lo_L,u_L,p_L,gama_L,lo_R,u_R,p_R,gama_R)
    global ep;
    put_out=linear_GRP_solver_Edir_Q1D([lo_L,lo_R,0,0,0,0,u_L,u_R,0,0,0,0,-0,-0,0,0,0,0,p_L,p_R,0,0,0,0,-0,-0,0,0,0,0,gama_L,gama_R,ep,ep]);
    p_IL =put_out(18);
    u_IL =put_out(16);
    p_IR =put_out(18);
    u_IR =put_out(16);
    lo_IL=put_out(15);
    lo_IR=put_out(17);
end
