function [t] = step_size_SDP(X,B)
%This function returns the largest t such that X+tB >=0
% with positive definite X and summetric matrix B

if find(eig(X)<0)
    disp('X is not positive definite')
    return
end
% Cholesky factorization X=R'R
R=chol(X);
% construct B'=-R'(-1)BR
B_prime=-R'\B;
B_prime=B_prime/R;
b=max(eig(B_prime));
if b>0
    t=1/b;
else 
    t=inf;
end
end
