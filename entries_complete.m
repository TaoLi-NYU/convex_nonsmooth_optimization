function [num,norm_diff,Xcomplete]=entries_complete(X,tol,repeat)
% This function determines the number of entries that are required 
%for approximately recovering  the matirx
%Inputs: X1: given matrix, tol: tolearnce for the norm (Xcomplete-X1)
% repeat: number of repeated experiemnts
%Outputs: num: number of entries are needed; norm_diff: norm of the
%difference; Xcomplete: completed matrix
M=size(X,1);
N=size(X,2);
sum=0;
%for each choice(number of entries), we run 5 experiments and average the
%results
for j=1:repeat
    for i=1:M*N % increase the number until the matrix is completed
        [Xcomplete,norm_diff] = matrix_complete(X,i);
        % matix i completed
        if norm_diff<tol
            num_entries=i;
            break
        end
    end
    sum=sum+num_entries;
end
%average
num=ceil(sum/repeat);