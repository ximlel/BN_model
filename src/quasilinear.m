function QUA = quasilinear(lo,u,p,p_g,u_s,phi_s,s_or_g,A_or_C)
%QUASILINEAR -- compute A(U) or C(U)
QUA = zeros(4,4);
global gama_s gama_g p0;
QUA(2,3)=1.0;
if strcmp(s_or_g,'g')==1
    if strcmp(A_or_C,'C')==1
        QUA(1,1)=QUA(1,1)+u_s;
        QUA(3,1)=QUA(3,1)-p_g;
        QUA(4,1)=QUA(4,1)-p_g*u_s;
    end
else
    QUA(3,1)=-gama_s*p0;
    QUA(4,1)=-gama_s*p0*u_s;
    if strcmp(A_or_C,'C')==1
        QUA(1,1)=QUA(1,1)+u_s;
        QUA(3,1)=QUA(3,1)+p_g;
        QUA(4,1)=QUA(4,1)+p_g*u_s;
    end
end
end

