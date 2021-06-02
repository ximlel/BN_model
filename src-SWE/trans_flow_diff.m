function diff = trans_flow_diff(h1,h2)
    diff = 2*0.18^2*(h2-h1)-9.81*h1*h2*(h1+h2);
end

