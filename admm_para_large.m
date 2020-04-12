function [X] = admm_para_large(z0,u0,A,b,lambda,rho, tol, maxit)
%This function provides reults of ADMM using different parameters
%input arguments are the same as lasso_admm_small
% we run through all possible combinations of paramenters
len_rho=length(rho);
len_lam=length(lambda);
N=size(A,2);
X=zeros(N,len_rho*len_lam);
% for each combination of parameters
for i=1:len_rho
    for j=1:len_lam
        r=rho(i);
        l=lambda(j);
        [x, ~,~,~]=lasso_admm_large(z0,u0,A,b,l,r, tol, maxit);
        % record the minimier
        index=(i-1)*len_lam+j;
        X(:,index)=x;
    end
end

end

