load('Adata.mat')
x0=zeros(100,1);
tol=1e-6;
maxit=1e3;
fun=@(x)analy_cen(x,A);
[f_all,~,h_all]=gradmeth(fun, x0,tol, maxit);
f_diff=f_all-f_all(end);
log_f=log(f_diff);
log_f=log_f(1:end-1);
X=1:length(log_f);
log_c=polyfit(X',log_f,1);
c=exp(log_c(1));
k=0.25/(1-c)
%%
len=size(h_all);
len=len(3);
m=inf;
M=-inf;

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
end

