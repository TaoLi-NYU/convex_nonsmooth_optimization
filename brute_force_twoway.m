% brute force in solving two-way partitioning
n=length(W1);
X=ones(n,1);
min=X'*W1*X;
minimizer= X;
for k=1:n
    C=combnk(1:n,k);
    n_rows = size(C,1);
    for i=1:n_rows
        for j=1:k
            X(C(i,j))=-1;
        end
        temp=X'*W1*X;
        if temp==min
            minimizer=[minimizer,X];
        end
        if temp<min
            min=temp;
            minimizer=X;
        end
    end
end


