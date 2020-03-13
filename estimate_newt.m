%estimate m, M, L, eta, gamma
%#iterations
%solve the problem using newton's method
load('Adata.mat')
x0=zeros(100,1);
tol=1e-8;
maxit=1e3;
fun=@(x)analy_cen(x,A);
[f_all,g_all,h_all,x_all]=newtmeth(fun, x0,tol, maxit);
len=size(h_all);
len=len(3);
%initialization
m=inf;
M=-inf;
L=0;
%
for i=1:len
    %for every hessian involved,c compute its largest eig and smallest
    %eigenvalue
    h=h_all(:,:,i);
    eigen=eig(h);
    eigmin=min(eigen);
    eigmax=max(eigen);
    %update m, if s samller one occurs
    if eigmin<m
        m=eigmin;
    end
    %update M, if a larger one occurs
    if eigmax>M
        M=eigmax;
    end
    %for two consecutive hessians, compute the difference quotient
    j=mod(i,len)+1;
    diff_h=h_all(:,:,j)-h;
    diff_x=x_all(:,j)-x_all(:,i);
    diff_quo=norm(diff_h)/norm(diff_x);
    if diff_quo>L
        L=diff_quo;
    end
end
%estimation for L, M, m
M=ceil(M)
L=ceil(L)
m
%estimation for eta gamma
alpha=0.25;
beta=0.5;
eta=3*(1-2*alpha)*m^2/L
gamma=alpha*beta*eta^2*m/M^2
% #iter
est_num=(f_all(1)-f_all(end))/gamma
% #iter for newton phase
e0=2*m^3/L^2;
net_num=log2(log2(e0/tol))