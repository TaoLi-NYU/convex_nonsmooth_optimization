clear variables
rng(1)
M=100000;
N=10000;
density=2/M;
A=sprandn(M,N,density);
figure(1)
spy(A)
b=randn(M,1);
z0=0;
u0=0;
tol=[1e-12,1e-12];
maxit=1e4;
% Cholesky 
rho=100;
B = A'*A+rho*speye(N); 
L=chol(B,'lower');
figure(2)
spy(L)
%% ADMM
lambda=10;
rho=100;
[x, pr_all,dr_all,p_all]=lasso_admm_large(z0,u0,A,b,lambda,rho, tol, maxit);
%% CVX
cvx_begin
    variable y(N)
    minimize(0.5*sum_square( A * y - b )+lambda*norm(y,1))
cvx_end
%% experiment with lambda
lambda=[2,10,100,1000];
rho=100;
[X] = admm_para_large(z0,u0,A,b,lambda,rho, tol, maxit);
for i=1:4
    sum(X(:,i)==0)
end
%% experiment with  rho
clear lambda rho X 
lambda=100;
rho=[1,10,100,1000, 10000,1e6];
[X] = admm_para_large(z0,u0,A,b,lambda,rho, tol, maxit);
