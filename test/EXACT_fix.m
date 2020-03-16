load EXACT.mat
lo_g_E=lo_g;
u_g_E=u_g;
p_g_E=p_g;
eta_g_E=eta_g;
lo_s_E=lo_s;
u_s_E=u_s;
p_s_E=p_s;
eta_s_E=eta_s;
for i=1:2000
   if lo_g(i)>2.67
       lo_g_E(i)=2.933;
   elseif lo_g(i)>2.2267
       lo_g_E(i)=2.4072;
   else
       lo_g_E(i)=2.0462;       
   end
   
   if i>1300 && u_g(i)<0.7955
       u_g_E(i)=0.7114;
   elseif u_g(i)>0.6467
       u_g_E(i)=0.8797;
   else
       u_g_E(i)=0.4136;       
   end
   
   if p_g(i)>2.1981
       p_g_E(i)=2.5;
   elseif p_g(i)>1.7029
       p_g_E(i)=1.8962;
   else
       p_g_E(i)=1.5096;
   end
   
   if eta_g(i)>0.554139177523263 || i<1300
       eta_g_E(i)=0.554246705750993;
   else
       eta_g_E(i)=0.554031649295533;       
   end
   
   
   if i<760
       lo_s_E(i)=0.5476;
   elseif i>=760  && i<=850
       lo_s_E(i)=0.5476+(0.4581-0.5476)*(i-760)/(850-760);
   elseif i>850  && i<1300
       if lo_s(i)<0.8175
           lo_s_E(i)=0.4581;
       else
           lo_s_E(i)=1.1769;       
       end      
   elseif i>=1300 && i<=1340
       lo_s_E(i)=1.04+(1.1769-1.04)*(1340-i)/(1340-1300); 
   else
       lo_s_E(i)=1.04;
   end
   
   if i<760
       u_s_E(i)=0;
   elseif i>=760  && i<=850
       u_s_E(i)=0.1605*(i-760)/(850-760);
   elseif i>850  && i<1300
       u_s_E(i)=0.1605;       
   elseif i>=1300 && i<=1340
       u_s_E(i)=0.1605*(1340-i)/(1340-1300); 
   else
       u_s_E(i)=0;
   end
   
   if i<760
       p_s_E(i)=0.328;
   elseif i>=760  && i<=850
       p_s_E(i)=0.328+(0.2556-0.328)*(i-760)/(850-760);
   elseif i>850  && i<1300
       if p_s(i)<0.8533
           p_s_E(i)=0.2556;
       else
           p_s_E(i)=1.4509;
       end          
   elseif i>=1300 && i<=1340
       p_s_E(i)=1.22+(1.4509-1.22)*(1340-i)/(1340-1300); 
   else
       p_s_E(i)=1.22;
   end

   if eta_s(i)<0.9585
       eta_s_E(i)=0.7621;
   else
       eta_s_E(i)=1.1548;       
   end
end
save EXACT_karni lo_s_E u_s_E p_s_E eta_s_E lo_g_E u_g_E p_g_E eta_g_E