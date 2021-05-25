function [h_mid,u_mid,H_t_mid,F,W_int]=GRP_solver(h_L,h_R,dh_L,dh_R,u_L,u_R,du_L,du_R,Z,dZ_L,dZ_R,d_t)
    global ep;
    global g;
    F=zeros(2,1);
    put_out=linear_GRP_solver_SWE([h_L,h_R,dh_L,dh_R,u_L,u_R,du_L,du_R,Z,dZ_L,dZ_R,g,ep]);
    h_mid = put_out(3) + 0.5*d_t*put_out(1);
    u_mid = put_out(4) + 0.5*d_t*put_out(2);
    H_t_mid = h_mid + 0.5*u_mid^2/g + Z;
    F(1) = h_mid*u_mid;
    F(2) = h_mid*u_mid^2 + 0.5*g*h_mid^2;
    h_int = put_out(3) + d_t*put_out(1);
    u_int = put_out(4) + d_t*put_out(2);
    W_int(1) = h_int+Z;
    W_int(2) = h_int*u_int;
end
