%t4: Riemann Solver exactly method
function [out_flux out_W_exact]=Riemann_solver_Exact(loL,loR,pL,pR,vel_uL,vel_uR,A,gama,ratio_x_t)
%state constant
both_rarefaction=1;
both_shock=2;
left_shock_right_rarefaction=3;
right_shock_left_rarefaction=4;
%state value
AL=2/(gama+1)/loL;
AR=2/(gama+1)/loR;
BL=(gama-1)/(gama+1)*pL;
BR=(gama-1)/(gama+1)*pR;
%judge wave pattern
p_min=min([pL pR]);
p_max=max([pL pR]);
d_vel_u=vel_uR-vel_uL;
aL=sqrt(gama*pL/loL);
aR=sqrt(gama*pR/loR);
f_min_L=2*aL/(gama-1)*((p_min/pL)^(1/2-1/2/gama)-1);
f_min_R=2*aR/(gama-1)*((p_min/pR)^(1/2-1/2/gama)-1);
f_min=f_min_L+f_min_R+d_vel_u;
f_max_L=(p_max-pL)*sqrt(AL/(p_max+BL));
f_max_R=(p_max-pR)*sqrt(AR/(p_max+BR));
f_max=f_max_L+f_max_R+d_vel_u;
if f_min>0&&f_max>0
    %both rarefaction
    wave_pattern=both_rarefaction;
elseif f_min<0&&f_max<0
    %both shock
    wave_pattern=both_shock;
elseif f_min<=0&&f_max>=0
    %shock and rarefaction
    if p_min==pL
        %left shock,right rarefaction
        wave_pattern=left_shock_right_rarefaction;
    else
        %right shock,left rarefaction
        wave_pattern=right_shock_left_rarefaction;
    end
else
    disp('error:occur while judging wave pattern');
end
%solve pStar by Newton interation
%guess p initial value
global ep
pPV=1/2*(pL+pR)-1/8*d_vel_u*(loL+loR)*(aL+aR);
p0=max([ep pPV]);
p(1)=p0;
k=1;
CHA=1;
while CHA>=ep
    switch wave_pattern
        case both_rarefaction,
            f_L(k)=2*aL/(gama-1)*((p(k)/pL)^(1/2-1/2/gama)-1);
            f_R(k)=2*aR/(gama-1)*((p(k)/pR)^(1/2-1/2/gama)-1);
            diff_f(k)=1/loL/aL*(p(k)/pL)^(-1/2-1/2/gama)+1/loR/aR*(p(k)/pR)^(-1/2-1/2/gama);
        case both_shock,
            f_L(k)=(p(k)-pL)*sqrt(AL/(p(k)+BL));
            f_R(k)=(p(k)-pR)*sqrt(AR/(p(k)+BR));
            diff_f(k)=sqrt(AL/(BL+p(k)))*(1-(p(k)-pL)/2/(BL+p(k)))+sqrt(AR/(BR+p(k)))*(1-(p(k)-pR)/2/(BR+p(k)));
        case left_shock_right_rarefaction,
            f_L(k)=(p(k)-pL)*sqrt(AL/(p(k)+BL));
            f_R(k)=2*aR/(gama-1)*((p(k)/pR)^(1/2-1/2/gama)-1);
            diff_f(k)=sqrt(AL/(BL+p(k)))*(1-(p(k)-pL)/2/(BL+p(k)))+1/loR/aR*(p(k)/pR)^(-1/2-1/2/gama);
        case right_shock_left_rarefaction,
            f_L(k)=2*aL/(gama-1)*((p(k)/pL)^(1/2-1/2/gama)-1);
            f_R(k)=(p(k)-pR)*sqrt(AR/(p(k)+BR));
            diff_f(k)=sqrt(AR/(BR+p(k)))*(1-(p(k)-pR)/2/(BR+p(k)))+1/loL/aL*(p(k)/pL)^(-1/2-1/2/gama);
        otherwise,
            disp('error:occur while computing derivative');
    end
    f(k)=f_L(k)+f_R(k)+d_vel_u;
    p(k+1)=p(k)-f(k)/diff_f(k);
    CHA=2*abs(p(k+1)-p(k))/(p(k+1)+p(k));
    k=k+1;
end
pStar=p(k);
%solve uStar
uStar=1/2*(vel_uL+vel_uR)+1/2*(f_R(k-1)-f_L(k-1));
%solve the complete solution
ratio_pStar_pL=pStar/pL;
ratio_pStar_pR=pStar/pR;
C1=(gama-1)/(gama+1);
C2=(gama-1)/2/gama;
C3=(gama+1)/2/gama;
switch wave_pattern
    case both_shock,
        loStar_L=loL*(ratio_pStar_pL+C1)/(C1*ratio_pStar_pL+1);
        SL=vel_uL-aL*sqrt(C3*ratio_pStar_pL+C2);
        loStar_R=loR*(ratio_pStar_pR+C1)/(C1*ratio_pStar_pR+1);
        SR=vel_uR+aR*sqrt(C3*ratio_pStar_pR+C2);
        if ratio_x_t<=uStar&&ratio_x_t>=SL
            W_0=[loStar_L;uStar;pStar];
        elseif ratio_x_t<SL
            W_0=[loL;vel_uL;pL];
        elseif ratio_x_t<=SR&&ratio_x_t>uStar
            W_0=[loStar_R;uStar;pStar];
        elseif ratio_x_t>SR
            W_0=[loR;vel_uR;pR];
        else
            disp('error:ratio_x_t runs out of real');
        end
    case both_rarefaction,
        loStar_L=loL*ratio_pStar_pL^(1/gama);
        loStar_R=loR*ratio_pStar_pR^(1/gama);
        aStar_L=aL*ratio_pStar_pL^C2;
        aStar_R=aR*ratio_pStar_pR^C2;
        SHL=vel_uL-aL;
        SHR=vel_uR+aR;
        STL=uStar-aStar_L;
        STR=uStar+aStar_R;
        if ratio_x_t<=SHL
            W_0=[loL;vel_uL;pL];
        elseif ratio_x_t>SHL&&ratio_x_t<STL
            W_0=[loL*(1/C3/gama+C1/aL*(vel_uL-ratio_x_t))^(1/C2/gama);1/C3/gama*(aL+C2*gama*vel_uL+ratio_x_t);pL*(1/C3/gama+C1/aL*(vel_uL-ratio_x_t))^(1/C2)];
        elseif ratio_x_t>=STL&&ratio_x_t<=uStar
            W_0=[loStar_L;uStar;pStar];
        elseif ratio_x_t>uStar&&ratio_x_t<=STR
            W_0=[loStar_R;uStar;pStar];
        elseif ratio_x_t>STR&&ratio_x_t<SHR
            W_0=[loR*(1/C3/gama-C1/aR*(vel_uR-ratio_x_t))^(1/C2/gama);1/C3/gama*(-aR+C2*gama*vel_uR+ratio_x_t);pR*(1/C3/gama-C1/aR*(vel_uR-ratio_x_t))^(1/C2)];
        elseif ratio_x_t>=SHR
            W_0=[loR;vel_uR;pR];
        else
            disp('error:ratio_x_t runs out of real');
        end
    case right_shock_left_rarefaction,
        loStar_R=loR*(ratio_pStar_pR+C1)/(C1*ratio_pStar_pR+1);
        SR=vel_uR+aR*sqrt(C3*ratio_pStar_pR+C2);
        loStar_L=loL*ratio_pStar_pL^(1/gama);
        aStar_L=aL*ratio_pStar_pL^C2;
        SHL=vel_uL-aL;
        STL=uStar-aStar_L;
        if ratio_x_t<=SHL
            W_0=[loL;vel_uL;pL];
        elseif ratio_x_t>SHL&&ratio_x_t<STL
            W_0=[loL*(1/C3/gama+C1/aL*(vel_uL-ratio_x_t))^(1/C2/gama);1/C3/gama*(aL+C2*gama*vel_uL+ratio_x_t);pL*(1/C3/gama+C1/aL*(vel_uL-ratio_x_t))^(1/C2)];
        elseif ratio_x_t>=STL&&ratio_x_t<=uStar
            W_0=[loStar_L;uStar;pStar];
        elseif ratio_x_t<=SR&&ratio_x_t>uStar
            W_0=[loStar_R;uStar;pStar];
        elseif ratio_x_t>SR
            W_0=[loR;vel_uR;pR];
        else
            disp('error:ratio_x_t runs out of real');
        end
    case left_shock_right_rarefaction,
        loStar_L=loL*(ratio_pStar_pL+C1)/(C1*ratio_pStar_pL+1);
        SL=vel_uL-aL*sqrt(C3*ratio_pStar_pL+C2);
        loStar_R=loR*ratio_pStar_pR^(1/gama);
        aStar_R=aR*ratio_pStar_pR^C2;
        SHR=vel_uR+aR;
        STR=uStar+aStar_R;
        if ratio_x_t<=uStar&&ratio_x_t>=SL
            W_0=[loStar_L;uStar;pStar];
        elseif ratio_x_t<SL
            W_0=[loL;vel_uL;pL];
        elseif ratio_x_t>uStar&&ratio_x_t<=STR
            W_0=[loStar_R;uStar;pStar];
        elseif ratio_x_t>STR&&ratio_x_t<SHR
            W_0=[loR*(1/C3/gama-C1/aR*(vel_uR-ratio_x_t))^(1/C2/gama);1/C3/gama*(-aR+C2*gama*vel_uR+ratio_x_t);pR*(1/C3/gama-C1/aR*(vel_uR-ratio_x_t))^(1/C2)];
        elseif ratio_x_t>=SHR
            W_0=[loR;vel_uR;pR];
        else
            disp('error:ratio_x_t runs out of real');
        end
    otherwise,
        disp('error:occur while computing complete solutions');
end
out_lo=W_0(1);
out_vel_u=W_0(2);
out_p=W_0(3);
out_E=out_p/(gama-1)+1/2*out_lo*out_vel_u^2;
out_W_exact=[out_lo;out_vel_u;out_p];
out_flux=[out_lo*out_vel_u;out_lo*out_vel_u*out_vel_u+out_p;out_E*out_vel_u+out_p*out_vel_u];
out_flux=A*out_flux;
end
