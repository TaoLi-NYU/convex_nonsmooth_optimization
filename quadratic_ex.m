%Quadratic case using grad_descent with backline 
n = 5;
A = hilb(n); % the Hilbert matrix, which is a very ill conditioned matrix
b = ones(n,1);
x0=ones(n,1);
%the true minimizer
p = -A\b;
[fp,]=quad(p,A,b);
[f0,]=quad(x0,A,b);
%f_0-f^*
delta=f0-fp;
% backtracking
alpha = 0.25;
beta = 0.5;
count = 0;
while count<100
    count=count+1;
    t=1;
    [f,g]=quad(x0,A,b);
    while true 
        x1=x0-t*g;
        [f1,]=quad(x1,A,b);
        %armijo condition
        f_arm=f-alpha*t*norm(g)^2;
        if f1>f_arm
            t=t*beta;
        else
            break;
        end
    end
    x0=x1;
end
%reduction factor and its root
ratio=(f1-fp)/delta 
factor=nthroot(ratio, count)
%compute the factor c=1-a*m/M
eigen = eig(A);
m=min(eigen);
M=max(eigen);
c=1-alpha*m/M;
