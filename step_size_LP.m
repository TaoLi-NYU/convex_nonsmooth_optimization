function [t] = step_size_LP(x,b)
%This function returns largest t such that x+tb>0
% with two given vectors x and b
if x<0
    disp('xx is not an interio point')
    return
end
% find the b_i<0
minus=find(b<0);
l=length(minus);
% if {b_i<0} is empty
if l==0
    t=inf;
else
    % compute xi/bi for each i in {b_i<0}
    for i=1:l
        t=x(minus)./b(minus);
        t=min(abs(t));
    end
end

end

