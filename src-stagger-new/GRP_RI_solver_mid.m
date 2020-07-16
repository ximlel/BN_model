function [phi_s_mid, rho_s_mid, u_s_mid, p_g_mid]=GRP_RI_solver_mid(rho_g,u_g,p_g,rho_s,u_s,p_s,phi_s,D,d_t)
    global ep;
    global gama_s gama_g;
    phi_g = 1.0-phi_s;
    c_s = sqrt(gama_s * p_s / rho_s);
    c_g = sqrt(gama_g * p_g / rho_g);
    Lambda = diag([u_s, u_s-c_s, u_s, u_s+c_s, u_g-c_g, u_g, u_g+c_g]);
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
    W_t = -(R*Lambda)/R*D;

    P = phi_g*rho_g*(u_g-u_s)^2+phi_g*p_g+phi_s*p_s;
    Q = phi_g*rho_g*(u_g-u_s);
    H = 0.5*(u_g-u_s)^2+gama_g/(gama_g-1.0)*p_g/rho_g;
    eta_g = p_g/rho_g^gama_g;
    
    phi_s_mid = phi_s + 0.5*d_t*W_t(1);
    rho_s_mid = rho_s + 0.5*d_t*W_t(2);
    u_s_mid = u_s + 0.5*d_t*W_t(3);
    Q_mid = Q + 0.5*d_t*W_t(5);
    H_mid = H + 0.5*d_t*W_t(6);
    eta_g_mid = eta_g + 0.5*d_t*W_t(7);
    phi_g_mid = 1.0-phi_s_mid;
    
    it_max = 500;
    k = 0; errB = 1e50;
    rho_g_mid = rho_g;
    while (k<it_max && errB>ep)
        fun  = H_mid-0.5*(Q_mid/phi_g_mid)^2/rho_g_mid^2-gama_g/(gama_g-1.0)*eta_g_mid*rho_g_mid^(gama_g-1.0);
        if abs(fun) < ep
            break;
        end
        dfun = (Q_mid/phi_g_mid)^2/rho_g_mid^3-gama_g*eta_g_mid*rho_g_mid^(gama_g-2.0);
        if dfun^2 < ep
            dfun = ep/dfun;
            break;
        end
        if abs(dfun) < ep
            fun
            break;
        end
        % 零特征值修正
        [rho_g_mid, errB] = NewtonRapshon(fun,dfun,rho_g_mid,ep);
        rho_g_mid=max(rho_g_mid,ep);
        k=k+1;
        errB = abs(fun);
    end
    if k>=it_max
        errB
    end
    p_g_mid = rho_g_mid^gama_g*eta_g_mid;
end

