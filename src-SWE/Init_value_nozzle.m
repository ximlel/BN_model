%% computational options
% N={640*9,640,480,320,240,160,120,80,40,20};
N = {100};
L = 1;
% tau = 0.05*L;
% time = (0:tau:5*L);
time = L*5;

MaxStp = 800000;
boundary = 1;
threshold = 10.0;
thickness = 1.8;
gamma = 1.4;
CFL = 0.1;
eps = 1e-9;
tol = 1e-13;
alp1 = 1.0;
alp2 = 2.0;
adp = 0;
prim = 0;
deri = 0;
lim = 0;
decomp = 1;
throat = L*0.25;

%% inital value parameters setup
mu2 = (gamma+1)/(gamma-1);
M0 = 0.12;
M1 = 3.0;
Ain = mach2section(M0, gamma);
Aex = mach2section(M1, gamma);
syms r
% A0 = Ain*exp(-log(Ain)*sin(pi*r/L).^2);
A0 = Ain*exp(-log(Ain)*sin(2*pi*r/L).^2);
A1 = Aex*exp(-log(Aex)*sin(2*pi*(1-r/L)/3).^2);
DA0 = diff(A0,r);
DA1 = diff(A1,r);

r0 = 1;
p0 = 1;
pin = p0*( 1.0+0.5*(gamma-1.0)*M0*M0 )^(gamma/(1.0-gamma));
rin = r0*( 1.0+0.5*(gamma-1.0)*M0*M0 )^(1.0/(1.0-gamma));
uin = M0*sqrt(gamma*pin/rin);
pex = p0*( 1.0+0.5*(gamma-1.0)*M1*M1 )^(gamma/(1.0-gamma));
rex = r0*( 1.0+0.5*(gamma-1.0)*M1*M1 )^(1.0/(1.0-gamma));
uex = M1*sqrt(gamma*pex/rex);
pex = 0.4;
rex = r0*((pex/p0)^(1/gamma));

%% reference values
nA = N{1};
h = L/nA;
volumeA = zeros(1,nA);
w=0.5*[0.1012285363,0.2223810345,0.3137066459,0.3626837834,0.3626837834,0.3137066459,0.2223810345,0.1012285363];
x = (0:h:L);
xx = zeros(1,8);

rhoA = zeros(1,nA);
uA = rhoA;
pA = rhoA;
err = rhoA;
for j=1:nA
    x0=0.5*(x(j)+x(j+1));
    xx(8) = 0.9602898565*(x(j+1)-x0)+x0;
    xx(7) = 0.7966664774*(x(j+1)-x0)+x0;
    xx(6) = 0.5255324099*(x(j+1)-x0)+x0;
    xx(5) = 0.1834346425*(x(j+1)-x0)+x0;
    xx(4) = 2*x0 - xx(5);
    xx(3) = 2*x0 - xx(6);
    xx(2) = 2*x0 - xx(7);
    xx(1) = 2*x0 - xx(8);
    aa = ones(1,8);
    for itG=1:8
        if xx(itG) > throat
            aa(itG) = double(subs(A1,r,xx(itG)));
        else
            aa(itG) = double(subs(A0,r,xx(itG)));
        end
    end
%     MM = section2mach(aa, gamma, xx > throat);
%     
%     rr = r0*(( 1.0+0.5*(gamma-1.0)*MM.^2).^(1.0/(1.0-gamma)));
%     pp = p0* ((1.0+0.5*(gamma-1.0)*MM.^2).^(gamma/(1.0-gamma)));
%     uu = MM .* sqrt(gamma*pp./rr);
%     
%     rr = rr.*aa;
%     pp = pp.*aa;
%     pp = 0.5*rr.*uu.*uu + pp/(gamma-1.0);
%     uu = uu .* rr;
%     rhoA(j) = dot(rr,w);
%     uA(j) = dot(uu,w);
%     pA(j) = dot(pp,w);
    volumeA(j) = dot(aa,w);
end


% load('~/consv_laws/DATA/nozzle/accurate-5760.mat');
% rhoA = [r0*ones(1,0.25*nA),r0*((pex/p0)^(1/gamma))*ones(1,0.75*nA)];
% uA = zeros(1,nA);
% pA = [p0*ones(1,0.25*nA),pex*ones(1,0.75*nA)];

%% piece-wise constant initial value
for j=1:nA
    if j < 0.25*nA
        rhoA(j) = r0;
        pA(j) = p0;
    else
        rhoA(j) = r0*((pex/p0)^(1/gamma));
        pA(j) = pex;
    end
%     rhoA(ceil(nA*0.25)) = 0.5*(r0+r0*((pex/p0)^(1/gamma)));
%     pA(ceil(nA*0.25)) = 0.5*(p0+pex);
end
if ~mod(nA,4)
    j0 = ceil(nA/4);
    x0 = (j0-1)*h;
    x1 = x0+h;
    rA(j0) = r0*(0.25*L-x0)/h + r0*((pex/p0)^(1/gamma))*(x1-0.25*L)/h;
    pA(j0) = p0*(0.25*L-x0)/h + pex*(x1-0.25*L)/h;
end
uA = zeros(1,nA);
rhoA = rhoA.*volumeA;
pA = 0.5*rhoA.*uA.*uA + pA.*volumeA/(gamma-1);
uA = rhoA.*uA;
volumeA = volumeA*h;
%% inital value
for it_n = 1:length(N)
    n = N{it_n};
    h = L/n;
    section = zeros(1,n+1);
    sec1 = zeros(1,n);
    sec2 = sec1;
    Dsec0 = section;
    Dsec1 = zeros(1,n);
    Dsec2 = Dsec1;
    
    %% sectional areas and Mach numbers
    w=0.5*[0.1012285363,0.2223810345,0.3137066459,0.3626837834,0.3626837834,0.3137066459,0.2223810345,0.1012285363];
    x = (0:h:L);
    xx = zeros(1,8);
    for k=1:n+1
        if k>0.25*n
            section(k) = double(subs(A1,r,x(k)));
            Dsec0(k) = double(subs(DA1,r,x(k)));
        else
            section(k) = double(subs(A0,r,x(k)));
            Dsec0(k) = double(subs(DA0,r,x(k)));
        end
    end
    for j=1:n
        x0=0.5*(x(j)+x(j+1));
%         xx(8) = 0.9602898565*(x(j+1)-x0)+x0;
%         xx(7) = 0.7966664774*(x(j+1)-x0)+x0;
%         xx(6) = 0.5255324099*(x(j+1)-x0)+x0;
%         xx(5) = 0.1834346425*(x(j+1)-x0)+x0;
%         xx(4) = 2*x0 - xx(5);
%         xx(3) = 2*x0 - xx(6);
%         xx(2) = 2*x0 - xx(7);
%         xx(1) = 2*x0 - xx(8);
% 
%         aa = double(subs(A0,r,xx));
%         volume(j) = dot(aa,w);
        xG = x0 - (0.5*h)/(sqrt(5));
        if(xG>throat)
            sec1(j) = double(subs(A1,r,xG));
            Dsec1(j) = double(subs(DA1,r,xG));
        else
            sec1(j) = double(subs(A0,r,xG));
            Dsec1(j) = double(subs(DA0,r,xG));
        end
        xG = x0 + (0.5*h)/(sqrt(5));
        if(xG>throat)
            sec2(j) = double(subs(A1,r,xG));
            Dsec2(j) = double(subs(DA1,r,xG));
        else
            sec2(j) = double(subs(A0,r,xG));
            Dsec2(j) = double(subs(DA0,r,xG));
        end
    end
    MI = section2mach(section, gamma, x > throat);
    xG = 0.5*(x(1:end-1)+x(2:end)) - 0.5*h/sqrt(5);
    MG1 = section2mach(sec1, gamma, xG > throat);
    xG = 0.5*(x(1:end-1)+x(2:end)) + 0.5*h/sqrt(5);
    MG2 = section2mach(sec2, gamma, xG > throat);

    %% initial values
    volume = zeros(1,n);
    rho = volume;
    u = volume;
    p = volume;
    ratio = nA / n;
    for j=1:ratio
        rho = rho + rhoA(j:ratio:end);
        u = u + uA(j:ratio:end);
        p = p + pA(j:ratio:end);
        volume = volume + volumeA(j:ratio:end);
    end
    rho = rho/ratio;
    u = u/ratio;
    p = p/ratio;
    Ma = [MI, MG1, MG2];


    %% output
    root = '~/consv_laws/DATA/';
    name = 'nozzle';
    fid = fopen([root,name,'/','RHO', name, '.txt'],'wt');
    fprintf(fid,'%14.8f\t', rho);
    fclose(fid);
    fid = fopen([root,name,'/','U', name, '.txt'],'wt');
    fprintf(fid,'%14.8f\t', u);
    fclose(fid);
    fid = fopen([root,name,'/','P', name, '.txt'],'wt');
    fprintf(fid,'%14.8f\t', p);
    fclose(fid);
    fid = fopen([root,name,'/','X',name,'.txt'],'wt');
    fprintf(fid,'%14.12f', h);
    fclose(fid);
    fid = fopen([root,name,'/','A0',name,'.txt'],'wt');
    fprintf(fid,'%14.12f\t', section);
    fclose(fid);
    fid = fopen([root,name,'/','A1',name,'.txt'],'wt');
    fprintf(fid,'%14.12f\t', sec1);
    fclose(fid);
    fid = fopen([root,name,'/','A2',name,'.txt'],'wt');
    fprintf(fid,'%14.12f\t', sec2);
    fclose(fid);
    fid = fopen([root,name,'/','DA0',name,'.txt'],'wt');
    fprintf(fid,'%14.12f\t', Dsec0);
    fclose(fid);
    fid = fopen([root,name,'/','DA1',name,'.txt'],'wt');
    fprintf(fid,'%14.12f\t', Dsec1);
    fclose(fid);
    fid = fopen([root,name,'/','DA2',name,'.txt'],'wt');
    fprintf(fid,'%14.12f\t', Dsec2);
    fclose(fid);
    fid = fopen([root,name,'/','VOL',name,'.txt'],'wt');
    fprintf(fid,'%14.12f\t', volume);
    fclose(fid);
    fid = fopen([root,name,'/','TIME',name,'.txt'],'wt');
    fprintf(fid,'%14.12f\n',time);
    fclose(fid);
    fid = fopen([root,name,'/','MACH',name,'.txt'],'wt');
    fprintf(fid,'%14.18f\n',Ma);
    fclose(fid);




    ampl_rho = 0.0;
    ampl_u = 0.0;
    ampl_p = 0.0;


    fid = fopen([root,name,'/','CONF',name,'.txt'],'wt');
    fprintf(fid,'gamma  %f\n', gamma);
    fprintf(fid,'CFL  %f\n', CFL);
    fprintf(fid,'eps  %g\n', eps);
    fprintf(fid,'alp1  %g\n', alp1);
    fprintf(fid,'alp2  %g\n', alp2);
    fprintf(fid,'bet1  120.0\n');
    fprintf(fid,'bet2  20.0\n');
    fprintf(fid,'modifier  0.03\n');
    fprintf(fid,'tol  %g\n', tol);
    fprintf(fid,'threshold  %f\n', threshold);
    fprintf(fid,'thickness  %f\n', thickness);
    fprintf(fid,'scaling  1.0\n');
    fprintf(fid,'MaxStp  %d\n', MaxStp);
    fprintf(fid,'Adp  %d\n', adp);
    fprintf(fid,'Boundary  %d\n', boundary);
    fprintf(fid,'Primative  %d\n', prim);
    fprintf(fid,'Deri  %d\n', deri);
    fprintf(fid,'Limiter  %d\n',lim);
    fprintf(fid,'Decomp   %d\n', decomp);
    fprintf(fid,'DA0      %f\n',double(subs(diff(A0,r),r,0)));
    fprintf(fid,'DA1      %f\n', double(subs(diff(A1,r),r,L)));
    fprintf(fid, 'ampl_rho %f\n', ampl_rho);
    fprintf(fid, 'ampl_u %f\n', ampl_u);
    fprintf(fid, 'ampl_p %f\n', ampl_p);
    fprintf(fid, 'rin %.15f\n', rin);
    fprintf(fid, 'uin %.15f\n', uin);
    fprintf(fid, 'pin %.15f\n', pin);
    fprintf(fid, 'pex %.15f\n', pex);
    fclose(fid);
end