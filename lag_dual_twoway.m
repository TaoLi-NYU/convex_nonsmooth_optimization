function [V,X] = lag_dual_twoway(W)
%Solving the Lagrangian dual of two-way partitioning problem
% input: weight matrix W
n = length(W);
cvx_begin sdp
    variable v(n)
    dual variable X
    maximize (-sum(v))
    W+diag(v)>= 0 : X;
cvx_end 
V=v;
end

