%example 1
clear all
tol = 1e-6;
maxit = 1e4;
n = 5;
H = hilb(n); % the Hilbert matrix, which is a very ill conditioned matrix
b = ones(n,1);
fun1 = @(x)quad(x, H, b);
x0=zeros(n,1);
% m and M for the Hessian matrix
m=min(eig(H));
M=max(eig(H));
% implement the gradient methods
[f1_all1,g1_all1] = grad_fixed(fun1,x0,1/M,tol,maxit);
[f1_all2,g1_all2] = grad_fixed(fun1,x0,2/(M+m),tol,maxit);
[f1_all3,g1_all3] = nes_grad(fun1,x0,m,M,tol,maxit);
%logplot
figure(1)
%#itertions 
iter1=1:length(f1_all1);
iter2=1:length(f1_all1);
iter3=1:length(f1_all3);
%optimal values
p=f1_all3(end);
%function eval logplot
semilogy(iter1,f1_all1-p,'ro');
hold on
semilogy(iter2,f1_all2-p,'bo');
semilogy(iter3,f1_all3-p,'go');
legend('gradient method t=1/M', 'gradient method t=2/(m+M)', 'Optimal gradient')
title('logplot of fucntion values using different methods')
xlabel('#iterations')
ylabel('log of function values')
% gradient norm logplot
figure(2)
semilogy(iter1,g1_all1,'rx');
hold on
semilogy(iter2,g1_all2,'bx');
semilogy(iter3,g1_all3,'gx');
legend('gradient method t=1/M', 'gradient method t=2/(m+M)', 'Optimal gradient')
title('logplot of norm of gradients using different methods')
xlabel('#iterations')
ylabel('log of norm of gradients')
%% example 2
clear all
tol = 1e-6;
maxit = 1e4;
load('Adata.mat');
n=length(A);
fun2=@(x)analy_cen(x,A);
x0=zeros(n,1);
% m and M for the Hessian matrices and we shall use the estimates obtained 
% in last assignment
m=2;
M=282;
% implement gradient methods
[f2_all1,g2_all1] = grad_fixed(fun2,x0,1/M,tol,maxit);
[f2_all2,g2_all2] = grad_fixed(fun2,x0,2/(M+m),tol,maxit);
[f2_all3,g2_all3] = nes_grad(fun2,x0,m,M,tol,maxit);
%logplot
figure(2)
%#itertions 
iter1=1:length(f2_all1);
iter2=1:length(f2_all2);
iter3=1:length(f2_all3);
%optimal values
p=f2_all3(end);
%logplot
figure(3)
semilogy(iter1,f2_all1-p,'ro');
hold on
semilogy(iter2,f2_all2-p,'bo');
semilogy(iter3,f2_all3-p,'go');
legend('gradient method t=1/M', 'gradient method t=2/(m+M)', 'Optimal gradient')
title('logplot of fucntion values using different methods')
xlabel('#iterations')
ylabel('log of function values')
%
figure(4)
semilogy(iter1,g2_all1,'rx');
hold on
semilogy(iter2,g2_all2,'bx');
semilogy(iter3,g2_all3,'gx');
legend('gradient method t=1/M', 'gradient method t=2/(m+M)', 'Optimal gradient')
title('logplot of norm of gradients using different methods')
xlabel('#iterations')
ylabel('log of norm of gradients')
%% example 3
clear all
tol = 1e-6;
maxit = 1e4;
dim=1e4;
m=1;
M=100;
% construct T matrix involved in gradient
T=2*speye(dim);
%second upper diagnoal
T(2:dim+1:end)=-1;
%lower diagnoal
T(dim+1:dim+1:end)=-1;
e1=sparse(zeros(dim,1));
e1(1)=1;
fun3=@(x)Nes_LCB(x,m,M,T,e1);
x0=sparse(zeros(dim,1));
% implement gradient methods
[f3_all1,g3_all1] = grad_fixed(fun3,x0,1/M,tol,maxit);
[f3_all2,g3_all2] = grad_fixed(fun3,x0,2/(M+m),tol,maxit);
[f3_all3,g3_all3] = nes_grad(fun3,x0,m,M,tol,maxit);
%logplot
%#itertions 
iter1=1:length(f3_all1);
iter2=1:length(f3_all2);
iter3=1:length(f3_all3);
%optimal values
p=f3_all3(end);
%logplot
figure(5)
semilogy(iter1,f3_all1-p,'ro');
hold on
semilogy(iter2,f3_all2-p,'bo');
semilogy(iter3,f3_all3-p,'go');
legend('gradient method t=1/M', 'gradient method t=2/(m+M)', 'Optimal gradient')
title('logplot of fucntion values using different methods')
xlabel('#iterations')
ylabel('log of function values')
%
figure(6)
semilogy(iter1,g3_all1,'rx');
hold on
semilogy(iter2,g3_all2,'bx');
semilogy(iter3,g3_all3,'gx');
legend('gradient method t=1/M', 'gradient method t=2/(m+M)', 'Optimal gradient')
title('logplot of norm of gradients using different methods')
xlabel('#iterations')
ylabel('log of norm of gradients')
%%
clear all
tol = 1e-6;
maxit = 1e4;
dim=1e4;
m=1;
M=1e4;
% construct T matrix involved in gradient
T=2*speye(dim);
%second upper diagnoal
T(2:dim+1:end)=-1;
%lower diagnoal
T(dim+1:dim+1:end)=-1;
e1=sparse(zeros(dim,1));
e1(1)=1;
fun4=@(x)Nes_LCB(x,m,M,T,e1);
x0=sparse(zeros(dim,1));
% implement gradient methods
[f4_all1,g4_all1] = grad_fixed(fun4,x0,1/M,tol,maxit);
[f4_all2,g4_all2] = grad_fixed(fun4,x0,2/(M+m),tol,maxit);
[f4_all3,g4_all3] = nes_grad(fun4,x0,m,M,tol,maxit);
%logplot
figure(3)
%#itertions 
iter1=1:length(f4_all1);
iter2=1:length(f4_all2);
iter3=1:length(f4_all3);
%optimal values
p=f4_all3(end);
%logplot
figure(7)
semilogy(iter1,f4_all1-p,'ro');
hold on
semilogy(iter2,f4_all2-p,'bo');
semilogy(iter3,f4_all3-p,'go');
legend('gradient method t=1/M', 'gradient method t=2/(m+M)', 'Optimal gradient')
title('logplot of fucntion values using different methods')
xlabel('#iterations')
ylabel('log of function values')
%
figure(8)
semilogy(iter1,g4_all1,'rx');
hold on
semilogy(iter2,g4_all2,'bx');
semilogy(iter3,g4_all3,'gx');
legend('gradient method t=1/M', 'gradient method t=2/(m+M)', 'Optimal gradient')
title('logplot of norm of gradients using different methods')
xlabel('#iterations')
ylabel('log of norm of gradients')