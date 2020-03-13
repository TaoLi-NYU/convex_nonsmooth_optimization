function [f_all,gnorm_all,h_all] = gradmeth(fun, x0,tol, maxit)
% code for gradient method, including backtracking line search, goes here,
% instead of the following dummy code, which just takes one gradient step
% and does not even check whether fun(x) < fun(x0).
% Note that the input parameter 'fun' is an anonymous function.
%set alpha beta
alpha = 0.25;
beta = 0.5;
count = 0;
e=1;
%pre-allocation
f_all=zeros(maxit,1);
gnorm_all=zeros(maxit,1);
h_all=zeros(length(x0),length(x0),maxit);
%e is above the tolerance and count is below maxit
while count<maxit && e>tol
    count=count+1;
    t=1;
    %get func eval and gradient
    [f,g,h]=fun(x0);
    f_all(count)=f;
    e=norm(g);
    gnorm_all(count)=e;
    h_all(:,:,count)=h;
    %backline tracking descent direction -g
    while true 
        x1=x0-t*g;
        [f1,~,~]=fun(x1);
        %armijo condition
        f_arm=f-alpha*t*e^2;
        if f1>f_arm
            t=t*beta;
        else
            break;
        end
    end
    x0=x1;
end
%we only keep the valid iterations
f_all=f_all(1:count);
gnorm_all=gnorm_all(1:count);
h_all=h_all(:,:,1:count);
p=f_all(end);
fprintf('the minimal value is %f',p)
iter=1:count;
figure(1)
semilogy(iter,f_all-p,'ro');
hold on
semilogy(iter,gnorm_all,'bx');
title('logplot of function values and norms of gradients')
xlabel('#iterations')
legend('func evaluation','norm of the gradient')

end
