function [SMTL_mean_distance] = get_SMTL_mean_distance(ggrid,r,s,nbr_targets,tdoas_measured, lambda)

[SMTL_weights_out,SMTL_weights_and_pos] = SMTL_wrap(ggrid,r,tdoas_measured,lambda);
[~,sort_inds] = sort(abs(SMTL_weights_out),'descend');
SMTL_weights_and_pos(:,sort_inds(1:nbr_targets));
SMTL_found_targets_pos = SMTL_weights_and_pos(2:end,sort_inds(1:nbr_targets));
%[est_to_any_source,source_to_any_est,~] = get_dists(SMTL_found_targets_pos,s);
%SMTL_mean_distance = sum(est_to_any_source.^2)/nbr_targets;
SMTL_mean_distance= mean(min(pdist2(s',SMTL_found_targets_pos'),[],2));
