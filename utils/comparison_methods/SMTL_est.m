function [activation_vec_out,ggrid_out] = SMTL_est(ggrid,sensor_pos,TDOA_mat,lambda)


% sensor_pos        -       3-times-M matrix with sensor positions
% ggrid             -       3-times-N matrix with grid coordinates
% TDOA_mat          -       (M-1)-times-K matrix with TDOA measurements for
%                           the K targets. The TDOAs are relative to the
%                           first sensor, i.e., they all have the first
%                           sensor as reference.
%                           Note that each column does not 
%                           have to correspond to a single target (no
%                           association required).
%

N =size(ggrid,2);
M = size(sensor_pos,2);

TDOA_vec = sum(TDOA_mat,2);

TOA_ref_grid = sqrt(sum((ggrid-sensor_pos(:,1)).^2));

TDOA_dictionary = zeros(M-1,N);
for k = 2:M
    TOA_temp = sqrt(sum((ggrid-sensor_pos(:,k)).^2));
    TDOA_dictionary(k-1,:) = TOA_ref_grid-TOA_temp;
end
TDOA_dictionary = -TDOA_dictionary;

%%% Pruning based on plausibility? %%%%
if 0
include_indicator = zeros(N,1);
for k = 1:N
    TDOA_temp = TDOA_dictionary(:,k);
    for ell = 1:size(TDOA_mat,2)
        min_dist = min(min(abs(TDOA_temp - TDOA_mat(:,ell)')));
    end
    dict_threshold = .01;
    if min_dist<dict_threshold
        include_indicator(k) = 1;
    end
end
TDOA_dictionary = TDOA_dictionary(:,include_indicator==1);
N = size(TDOA_dictionary,2);
ggrid_out = ggrid(:,include_indicator==1);
else
    ggrid_out = ggrid;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Proposed in paper: augmentation %%%%
if 1
    L = 5; % 5 in paper
    A = TDOA_dictionary;
    y = TDOA_vec;
    for ell = 2:L
        f = @(x) (x).^ell;
        A = [A;f(TDOA_dictionary)];
        y = [y;sum(f(TDOA_mat),2)];
    end
else
    A = TDOA_dictionary;
    y = TDOA_vec;
end

if 0
    [U,S,V] = svd(A);
    R = pinv(S)*U';

    A = R*A;
    y = R*y;
end


%%%% Pruning? %%%%%
if 0
inner_prods = abs(A'*y);
inner_prods = inner_prods/max(inner_prods);
corr_threshold = 0.00;
inds_above_threshold = find(inner_prods>corr_threshold);
A = A(:,inds_above_threshold);
N = size(A,2);
ggrid_out = ggrid(:,inds_above_threshold);
end
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cvx_begin quiet
variable x(N,1)
minimize(.5*sum_square_abs(A*x-y)+lambda*sum(abs(x)))
%minimize(.5*sum_square_abs(A*x-y))
%subject to
%sum(x)==3
x>=0 % OBS!!!!!!!
%x<=1
cvx_end

activation_vec_out = x;

%alpha = 1.8;
%rho = 10;
%[z, history] = boyd_admm_lasso(A, y, lambda, rho, alpha);


% Sparsity-Aware Multi-Source TDOA Localization