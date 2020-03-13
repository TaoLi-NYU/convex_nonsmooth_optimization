function [S,S1,cut] = max_cut(W)
%Solving max_cut problem
% W: weight matrix; r: normal vector
n=length(W);
r=ones(n,1)/sqrt(n);
S=[];
S1=[];
cut=0;
% X=V'V where X is from SDP
V=spd_relax_eig(W);
for k=1:n
    vr=V(:,k)'*r;
    % partitioning based on the inner product
    if vr>=0
        S=[S,k];
    else
        S1=[S1,k];
    end
end            
    
for i=1:length(S)
    for j=1:length(S1)
        cut=cut+W(S(i),S1(j));
    end
end

end


