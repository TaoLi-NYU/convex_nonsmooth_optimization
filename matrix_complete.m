function [Xcomplete,norm_diff] = matrix_complete(X,num_entries)
%This function performs matrix completion using SDP in CVX 
%Inputs: X: original matrix; num_tries: given known entries
%Output: Xcomplete: the completed matrix
% obtain the dimensions
M=size(X,1);
N=size(X,2);
% if num_entries exceed the total numner
if num_entries>M*N
    disp('error in the input')
    return
end
% randomly generate row and column index pairs
rng(1)% for random experiments comment this 
row_index=randi(M,num_entries,1);
col_index=randi(N,num_entries,1);
% SDP nuclear form
cvx_begin sdp
    variable W1(M,M) symmetric
    variable W2(N,N) symmetric
    variable X1(M,N)
    minimize trace(W1)+trace(W2)
    [W1, X1;X1', W2]>=0;
    for i=1:num_entries
        X1(row_index,col_index(i))==X(row_index,col_index(i));
    end
cvx_end
Xcomplete=X1;
norm_diff=norm(X1-X,'fro');
end

