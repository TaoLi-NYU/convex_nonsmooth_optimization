function [x, pr_all,dr_all,p_all]=lasso_admm_small(z0,u0,A,b,lambda,rho, tol, maxit)
% using ADMM method to solve lasso problem 1/2|Ax-b|^2+lambda|x|_1
%inputs: x0,z0: primal variable; u0: dual(scaled) variable
%        A,b,lambda: parameters involved in lasso; rho: parameter in the augmented lagrangian
%        tol: tolerance in stopping criterion; maxit: maximum iteration numbers
%     
%outputs: x
% initializaion
tol_abs=tol(1);
tol_rel=tol(2);
count=0;
M=size(A,1);
N=size(A,2);
I=eye(N);
pr_all=zeros(maxit,1);
dr_all=zeros(maxit,1);
p_all=zeros(maxit,1);
while count<maxit
    %ADMM update
    x1=(A'*A+rho*I)\(A'*b+rho*(z0-u0));
    z1=wthresh(x1+u0,'s',lambda/rho);
    u1=u0+x1-z1;
    %primal residual
    pr=norm(x1-z1);
    %dual residual
    dr=rho*norm(z1-z0);
    %collect data
    count=count+1;
    pr_all(count)=pr;
    dr_all(count)=dr;
    p_all(count)=0.5*norm(A*x1-b)^2+lambda*norm(x1,1);
    %stopping criterion
    e_pri=sqrt(N)*tol_abs+tol_rel*max(norm(x1),norm(z1));
    e_dual=sqrt(M)*tol_abs/rho+tol_rel*norm(u1);
    if pr<e_pri && dr<e_dual
        break
    end
    %assign new values 
    z0=z1;
    u0=u1;
end
x=x1;
pr_all=pr_all(1:count);
dr_all=dr_all(1:count);
p_all=p_all(1:count);

% logplot
figure()
%#itertions 
iter=1:length(pr_all);
%optimal values
p=p_all(end);
%logplot
semilogy(iter,p_all-p,'ro');
hold on
semilogy(iter,pr_all,'go');
semilogy(iter,dr_all,'bo');
legend('function evaluation', 'primal residual', 'dual residual')
title(['logplot of quantities involved, lambda= ', num2str(lambda), '  rho= ', num2str(rho)])
xlabel('#iterations')
ylabel('log of different values')
end