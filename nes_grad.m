function [f_all,g_all] = nes_grad(fun,x0,m,M,tol,maxit)
%Nesterov optimal gradient descent
%fun: function to be minimized, x0: initial point
%m,M: lower and upper bound of the eigenvalues
%tol, maxit: tolerance and maximum number of iterations 
count=0;
e=1;
%pre-allocation
f_all=zeros(maxit,1);
g_all=zeros(maxit,1);
%condition number of q in Nesterov
k=M/m;
q=(sqrt(k)-1)/(sqrt(k)+1);
y0=x0;
while count<maxit && e>tol
    count=count+1;
    [f,g]=fun(x0);
    y1=x0-1/M*g;
    x1=(1+q)*y1-q*y0;
    f_all(count)=f;
    e=norm(g);
    g_all(count)=e;
    x0=x1;
    y0=y1;
end
f_all=f_all(1:count);
g_all=gnorm_all(1:count);
end
