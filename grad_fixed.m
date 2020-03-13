function [f_all,g_all] = grad_fixed(fun,x0,t,tol,maxit)
%Gradient descent method with fixed step size
%fun: function to be minimized, x0: starting point
%tol: tolerance, maxit: maximum number of iterations, t: step size
count=0;
e=1;
%pre-allocation
f_all=zeros(maxit,1);
gnorm_all=zeros(maxit,1);
while count<maxit && e>tol
    count=count+1;
    [f,g]=fun(x0);
    x1=x0-t*g;
    f_all(count)=f;
    e=norm(g);
    gnorm_all(count)=e;
    x0=x1;
end
f_all=f_all(1:count);
g_all=gnorm_all(1:count);
end
