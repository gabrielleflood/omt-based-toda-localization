function [weights_out,weights_and_pos] = SMTL_wrap(candidates,r,tdoas_measured,lambda)

ggrid = candidates;
sensor_pos = r;

M = size(sensor_pos,2);

% Handle false/missing TDOAs
max_nbr_targets = 0;
for k = 2:M
    temp_TDOA = tdoas_measured{1,k}.tdoas;
    max_nbr_targets = max(max_nbr_targets,length(temp_TDOA));
end

TDOA_mat = zeros(M-1,max_nbr_targets);

for k = 2:M
    temp_TDOA = tdoas_measured{1,k}.tdoas;
    temp_TDOA = temp_TDOA(:)';
    K_temp = length(temp_TDOA);
    if K_temp<max_nbr_targets
        temp_TDOA = [temp_TDOA,zeros(1,max_nbr_targets-K_temp)];
    end

    TDOA_mat(k-1,:) = temp_TDOA;
end

[activation_vec_out,new_grid] = SMTL_est(ggrid,sensor_pos,TDOA_mat,lambda);

%weights_and_pos = [activation_vec_out';ggrid];

weights_and_pos = [activation_vec_out';new_grid];

weights_out = activation_vec_out;
end


