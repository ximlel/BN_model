function C = minmod2(A,B)
%MINMOD FUNCTION
for i=1:size(A,1)
    if A(i)>0 && B(i)>0
        C(i) = min([A(i) B(i)]);
    elseif A(i)<0 && B(i)<0
        C(i) = max([A(i) B(i)]);
    else
        C(i) = 0.0;
    end
end
end

