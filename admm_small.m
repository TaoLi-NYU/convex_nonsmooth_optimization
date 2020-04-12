% initialization
clear variables 
rng(1)
N=50;
A=sprandn(N,N,3/4);
b=randi([0,10],N,1);
z0=0;
u0=0;
tol=[1e-5,1e-5];
maxit=1e4;
%% ADMM
lambda=10;
rho=100;
[x, pr_all,dr_all,p_all]=lasso_admm_small(z0,u0,A,b,lambda,rho, tol, maxit);

%% CVX
cvx_begin
    variable y(N)
    minimize(0.5*sum_square( A * y - b )+lambda*norm(y,1))
cvx_end
%% Experiment with lambda
lambda=[2,10,50,100,1000,5000];
rho=100;
[X] = admm_para(z0,u0,A,b,lambda,rho, tol, maxit);
for i=1:6
    sum(X(:,i)<1e-5)
end
%% Experiment with rho
clear lambda rho X 
lambda=100;
rho=[1,10,100,1000];
[X] = admm_para(z0,u0,A,b,lambda,rho, tol, maxit);