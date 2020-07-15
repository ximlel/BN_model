%compute flux
function RI=U2RI_cal(phi_s,lo_g,u_g,p_g,u_s,p_s,lo_s)
global gama_g;
global ep;
    phi_g=1.0-phi_s;
    eta_g=p_g/lo_g^gama_g;
    Q=phi_g*lo_g*(u_g-u_s);
    P=phi_g*lo_g*(u_g-u_s)^2+phi_g*p_g+phi_s*p_s;
    H=0.5*(u_g-u_s)^2+gama_g/(gama_g-1.0)*p_g/lo_g;
    RI(1)=lo_s;
    RI(2)=u_s;
    RI(3)=P;
    RI(4)=Q;
    RI(5)=H;
    RI(6)=eta_g;
end
