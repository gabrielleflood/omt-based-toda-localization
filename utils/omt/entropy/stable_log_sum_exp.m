function val_out = stable_log_sum_exp(X,onesN)
%
% For a matrix, this function compute log(sum(exp(.)) of each row and
% returns this as a column vector.
%
% onesN - column vector of all 1's. Length should match row-length of X.
%

if nargin<2
    max_X = max(X);
    val_out = max_X + log(sum(exp(X - max_X)));
else
    max_X = max(X,[],2);
    val_out = max_X+ log(sum(exp(X - max_X*onesN'),2));
end
end