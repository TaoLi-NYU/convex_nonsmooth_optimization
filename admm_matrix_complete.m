% random
clear 
rank=2;
dim=5;
X=-2+4*rand(dim); %entries are within [-2,2]
[U,S,V]=svd(X);
X= U(:, 1: rank)* S(1: rank, 1: rank)* V(:, 1: rank)';
num_entries=10;
%%
rho=100;
tol=[1e-8,1e-8];
maxit=10000;
disp('ADMM')
[Xcomplete2,norm_diff2,obj2] = admm_complete(X,num_entries,rho, tol, maxit);
disp('CVX')
[Xcomplete1,norm_diff1] = matrix_complete(X,num_entries);
%%
rank=4;
dim=10;
X=-2+4*rand(dim); %entries are within [-2,2]
[U,S,V]=svd(X);
X= U(:, 1: rank)* S(1: rank, 1: rank)* V(:, 1: rank)';
num_entries=30;
%%
disp('ADMM')
[Xcomplete4,norm_diff4,obj4] = admm_complete(X,num_entries,rho, tol, maxit);
disp('CVX')
[Xcomplete3,norm_diff3] = matrix_complete(X,num_entries);
%%
rank=5;
dim=20;
X=-2+4*rand(dim); %entries are within [-2,2]
[U,S,V]=svd(X);
X= U(:, 1: rank)* S(1: rank, 1: rank)* V(:, 1: rank)';
num_entries=50;
%%
disp('ADMM')
[Xcomplete6,norm_diff6,obj6] = admm_complete(X,num_entries,rho, tol, maxit);
disp('CVX')
[Xcomplete5,norm_diff5] = matrix_complete(X,num_entries);