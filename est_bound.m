function [logc]=est_bound(f_all,range)
% estmiate the slope of the logplot
%input: f_all: all function evaluations, range: number of iterations
%considered in linear fit
%obtain f(x^k)-p^*
f_diff=f_all-f_all(end);
f_diff=f_diff(range);
% obtain the log value
log_f=log(f_diff);
%use least square fit
log_c=polyfit(range',log_f,1);
logc=log_c(1);
end