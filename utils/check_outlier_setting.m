function [is_outlier_setting,CRLB_source_pos,cond_nbrs] = check_outlier_setting(ind_set,r,s,condition_threshold)

% Threshold for condition number
if nargin < 4 || isempty(condition_threshold)
    condition_threshold = 100;
end

% Threshold CRLB


nbr_pairs = size(ind_set,1);
S = size(s,2);

prop_speed_nominal = 1;
sigma2_noise_nominal = 1;

rec1_rec2_pair_cell = cell(nbr_pairs,1);
for k_tdoa_pair = 1:nbr_pairs
    rec1_rec2_pair_cell{k_tdoa_pair} = ...
        [r(:,ind_set(k_tdoa_pair,1)),r(:,ind_set(k_tdoa_pair,2))];
end

cond_nbrs = zeros(S,1);
CRLB_source_pos = zeros(S,1);
for k_source = 1:S
    [CRLB_source_pos_temp,CRLB_mat,tot_FIM] = ...
        get_CRLB_source_pos(s(:,k_source),rec1_rec2_pair_cell,prop_speed_nominal,sigma2_noise_nominal);
    cond_nbrs(k_source) = cond(tot_FIM);
    CRLB_source_pos(k_source) = CRLB_source_pos_temp;
end

if max(cond_nbrs)>condition_threshold %|| max(CRLB_source_pos)>
    is_outlier_setting = 1;
else
    is_outlier_setting = 0;
end
