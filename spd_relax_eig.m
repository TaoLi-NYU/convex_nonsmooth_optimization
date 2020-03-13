function [V] = spd_relax_eig(W)
%Solving the SDP relaxation of two-partitioning using X=VV
%input matrix W
n=length(W);
cvx_begin sdp
    variable X(n,n) symmetric
    minimize(trace(W*X))
    X>=0 ;
    diag(X)==ones(n,1) ;
cvx_end
[Q,D]=eig(X);
D_sqrt=diag(sqrt(diag(D)));
V=D_sqrt*Q';
end

