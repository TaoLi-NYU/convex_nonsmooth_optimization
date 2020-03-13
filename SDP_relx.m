function [Z,v,X,tr] = SDP_relx(W)
%Solving the SDP relaxation of two-partitioning
%input matrix W
n=length(W);
cvx_begin sdp
    variable X(n,n) symmetric
    dual variables Z v
    minimize(trace(W*X))
    X>=0 :Z;
    diag(X)==ones(n,1) :v ;
cvx_end
tr=trace(W*X);
end

