function QUA = quasilinear(lo,u,p,p_g,u_s,s_or_g,A_or_C)
%QUASILINEAR -- compute A(U) or C(U)
QUA = zeros(4,4);
global gama_s gama_g p0;
QUA(2,3)=1.0;
if strcmp(s_or_g,'g')==1 %gas phase
    E=p/(gama_g-1)+0.5*lo*u^2;
    H=(E+p)/lo;
    a=sqrt(gama_g*p/lo);
    QUA(3,2)=(gama_g-1)*H-u^2-a^2;
    QUA(3,3)=(3-gama_g)*u;
    QUA(3,4)=gama_g-1;
    QUA(4,2)=u*(-H+0.5*(gama_g-1)*u^2);
    QUA(4,3)=H-(gama_g-1)*u^2;
    QUA(4,4)=gama_g*u;
    if strcmp(A_or_C,'C')==1
        QUA(1,1)=u_s;
        QUA(3,1)=QUA(3,1)+p_g;
        QUA(4,1)=QUA(4,1)+p_g*u_s;
    elseif strcmp(A_or_C,'A')~=1
        error('NOT A(U) or C(U)!');
    end
elseif strcmp(s_or_g,'s')==1%solid phase
    QUA(3,1)=-gama_s*p0;
    QUA(4,1)=-gama_s*p0*u_s;
    E=(p+gama_s*p0)/(gama_s-1)+0.5*lo*u^2;
    H=(E+p)/lo;
    a=sqrt(gama_s*(p+p0)/lo);
    QUA(3,2)=(gama_s-1)*H-u^2-a^2;
    QUA(3,3)=(3-gama_s)*u;
    QUA(3,4)=gama_s-1;
    QUA(4,2)=u*(-H+0.5*(gama_s-1)*u^2);
    QUA(4,3)=H-(gama_s-1)*u^2;
    QUA(4,4)=gama_s*u;
    if strcmp(A_or_C,'C')==1
        QUA(1,1)=u_s;
        QUA(3,1)=QUA(3,1)-p_g;
        QUA(4,1)=QUA(4,1)-p_g*u_s;
    elseif strcmp(A_or_C,'A')~=1
        error('NOT A(U) or C(U)!');
    end
else
    error('NOT gas or solid in QUASILINEAR!');
end
end

