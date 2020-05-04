function [Norm_all] = matrix_complete_thresh(rank,dim)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% randomly generate a square matirx of dim-dimension
% use SVD 
rng(1)
X=-2+4*rand(dim); %entries are within [-2,2]
[U,S,V]=svd(X);
% alter the rank by SVD 
X= U(:, 1: rank)* S(1: rank, 1: rank)* V(:, 1: rank)';
% tolerance
tol=2e-9;
Norm_all=zeros(1,dim*dim);
for i=1:dim*dim
    %CVX for completion
    [~,norm_diff] = matrix_complete(X,i);
    Norm_all(i)=norm_diff;
    % once the norm of the difference is below tol for ten iterations
    if Norm_all(i)<tol
        if Norm_all(i-9:i)<tol*ones(1,10)
            break
        end
    end
end
% logplot
Norm_all=Norm_all(1:i);
figure()
semilogy(1:i,Norm_all,'o');
title('logplot of norm difference')
xlabel('number of entries')
ylabel('log of norm difference')
end
