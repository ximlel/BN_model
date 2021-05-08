function C = minmod(A,C,B)
%MINMOD FUNCTION
for i=1:size(C,1)
    if C(i)>0 && A(i)>0 && B(i)>0
        C(i) = min([A(i) B(i) C(i)]);
    elseif C(i)<0 && A(i)<0 && B(i)<0
        C(i) = max([A(i) B(i) C(i)]);
    else
        C(i) = 0.0;
    end
end
end

