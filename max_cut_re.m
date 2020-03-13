function [max_cut_re]=max_cut_re(W,sdp)
% recover the max cut from the optimal of SDP
% from \sum_i\sum_j w_ij x_ij to 1/2 \sum_{i<j}w_{ij} (1-y_i*y_j) 
partial_tr=0.5*(sdp-sum(diag(W)));
partial_sum=0.5*(sum(sum(W))-sum(diag(W)));
max_cut_re=0.5*(partial_sum-partial_tr);

