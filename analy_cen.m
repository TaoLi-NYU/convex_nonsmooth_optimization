function [f,g,h]=analy_cen(x,A)
%function defined in BV p519 for computing the analiyitc center of linear
%inequalities
[n,m]=size(A);
a=A'*x;
% out of domain
if sum(a>=1)||sum(abs(x)>=1)
    f=inf;
    g=nan*ones(n,1);
    h=nan*ones(n,n);
% compute the function value, gradient and the hessian
else
    %function evaluation
    f=-sum(log(ones(m,1)-a))-sum(log(ones(n,1)-x.*x));
    %gradient
    g=A*(1./(1-A'*x))+1./(1-x)-1./(1+x);
    %hessian
    d=1./(1-A'*x);
    h=A*diag(d.^2)*A'+diag(1./(1+x).^2)+diag(1./(1-x).^2);
end
end