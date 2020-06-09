function [x_star, err] = NewtonRapshon(fun,dfun,x0,ep)
if norm(fun,inf) <= ep
    d = 0.0;
else
    d = -fun/dfun;
end
x_star = x0 + d;
err = norm(d,inf);
end
