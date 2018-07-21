function C = minmod(Alpha,A,C,B)
%MINMOD FUNCTION
if nargin < 3
    B=sign(A)*1e30;
end
for i=1:size(C,1)
    if C(i)>0 && A(i)>0 && B(i)>0
        C(i) = min([Alpha*A(i) Alpha*B(i) C(i)]);
    elseif C(i)<0 && A(i)<0 && B(i)<0
        C(i) = max([Alpha*A(i) Alpha*B(i) C(i)]);
    else
        C(i) = 0.0;
    end
end
end

