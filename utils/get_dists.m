function [est_to_any_source,source_to_any_est,source_to_grid] = get_dists(est_pos,source_pos,ggrid)
%
% Computes distances betweens sets of estimated target positions and ground
% truth source positions.
%
%

N_ests = size(est_pos,2);
est_to_any_source = zeros(N_ests,1);

N_sources = size(source_pos,2);
source_to_any_est = zeros(N_sources,1);

if nargin>=3
    N_ggrid = size(ggrid,2);
    source_to_grid = zeros(N_sources,1);
else
    N_ggrid = 0;
    source_to_grid = [];
end


for k_est = 1:N_ests
    est_to_any_source(k_est) = min(sqrt(sum((est_pos(:,k_est)-source_pos).^2)));
end
for k_source = 1:N_sources
    source_to_any_est(k_source) = min(sqrt(sum((source_pos(:,k_source)-est_pos).^2)));
end
if N_ggrid>0
    for k_source = 1:N_sources
        source_to_grid(k_source) = min(sqrt(sum((source_pos(:,k_source)-ggrid).^2)));
    end
end
end



