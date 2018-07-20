function [p_M,u_M,lo_M]=value_cal_GRP(lo,u,p,dlo,du,dp,gama,d_x,d_t)
    global ep;
    lo=lo+0.5*d_x*dlo;
    u =u +0.5*d_x*du;
    p =p +0.5*d_x*dp;
    put_out=linear_GRP_solver_Edir_Q1D([lo,lo,dlo,dlo,0,0,u,u,du,du,0,0,-0,-0,0,0,0,0,p,p,dp,dp,0,0,-0,-0,0,0,0,0,gama,gama,ep,ep]);		
    lo_M = put_out(9)  + d_t*put_out(3);	
    u_M  = put_out(10) + d_t*put_out(4);			
    p_M  = put_out(12) + d_t*put_out(6);   
end

