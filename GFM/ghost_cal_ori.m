function [p_IL,u_IL,lo_IL,p_IR,u_IR,lo_IR]=ghost_cal_ori(lo_L,u_L,p_L,gama_L,lo_R,u_R,p_R,gama_R,lo_ML,u_ML,p_ML,lo_MR,u_MR,p_MR)
    global ep;
    put_out=RP_solver([lo_L,lo_R,u_L,u_R,p_L,p_R,gama_L,gama_R,ep]);
    p_IL =put_out(1);
    u_IL =put_out(2);
    p_IR =put_out(1);
    u_IR =put_out(2);
%    lo_IL=put_out(3);
%    lo_IR=put_out(4);
%    lo_IL=(p_IL/p_R)^(1/gama_R)*lo_R;
%    lo_IR=(p_IR/p_L)^(1/gama_L)*lo_L;
    lo_IR=put_out(3);
    lo_IL=put_out(4);
end
