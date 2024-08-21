function [M_out,mu_out] = omt_cluster_solv(nbr_toas,card_grid,nbr_targets,nbr_receivers,cost_mat)

%
% OMT solver for clustering problem with additional trash state
%

% size check
[row_C,cols_C] = size(cost_mat);

if row_C ~=nbr_toas
    error('nbr of rows of cost matrix not equal to nbr ToAs')
end
if cols_C ~=card_grid+1
    error('nbr of columns of cost matrix not equal to grid cardinality +1')
end

S = nbr_targets;
R = nbr_receivers;

one_vec_ToA = ones(nbr_toas,1);
one_vec_grid = ones(card_grid,1);

%%%% Cost function %%%%%
cost_func = @(M) sum(M(:).*cost_mat(:));

%%%%% Constraints %%%%%
margin_toa = @(M) sum(M,2)-one_vec_ToA;
margin_grid = @(M) sum(M(:,1:end-1),1)'-R*(R-1)/2*one_vec_grid;
activation_bound = @(M,my_mu) M(:,1:end-1) - one_vec_ToA*my_mu';
target_constraint = @(my_mu) sum(my_mu)-S;

%%%% CVX %%%%
cvx_begin quiet
variable M(nbr_toas,card_grid+1)
variable my_mu(card_grid,1)
%
minimize(cost_func(M))
%
subject to
%
margin_toa(M) == 0;
margin_grid(M) <= 0;
%
activation_bound(M,my_mu)<=0;
target_constraint(my_mu) == 0;
%
M>=0;
M<=1;
%
my_mu>=0;
my_mu<=1;
%
cvx_end
%%%%%%%%%%%%%%%

%
M_out = M;
mu_out = my_mu;

