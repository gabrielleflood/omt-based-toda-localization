function [s_refined,mean_distance_refined] = refine_using_association(tdoas_listed, S, M_out, best_pos,  candidates, r,s)


M = M_out(:,best_pos(1:S));
%M_binary = M.*(M>0.7); % conservative, only take measurement if certain.
M_binary = (M>0.7)+0;

%M_binary = create_binary_M(M_out, [best_pos; size(M_out,2)]);
%tdoa_association_vec = zeros(size(M_out,2),1);
%for k_tdoa = 1:size(M_out,2)
%    tdoa_association_vec(k_tdoa) = find(M_binary(k_tdoa,:));
%end
s_refined = zeros(3,S);
% OBS24! Update the for loop s.t. we do not get an error if the number of
% found candidates is smaller than S? No problem, as it is only called if S
% sources are found.
for k_source = 1:S
    %tdoa_indices = find(M_binary(:,best_pos(k_source)));
    tdoa_indices = find(M_binary(:,k_source));
    temp_tdoa_rec_mat = zeros(length(tdoa_indices),3);
    for k_tdoa = 1:length(tdoa_indices)
        temp_struct = tdoas_listed(tdoa_indices(k_tdoa));
        temp_tdoa_rec_mat(k_tdoa,:) = [temp_struct.tdoas,temp_struct.r1,temp_struct.r2];
    end
    % OBS24! Have s_found as an input an change the below line to s_found(:,k_source)?
    s0 = candidates(:,best_pos(k_source));
    s_est_temp = NLS_refine_per_source(temp_tdoa_rec_mat,r,s0);
    s_refined(:,k_source) = s_est_temp;
end
s_eucl_dist_refined = pdist2(s',s_refined');
mean_distance_refined = sum((min(s_eucl_dist_refined,[],1)))/S;
%mean_distance_refined = sum((min(s_eucl_dist_refined,[],1)).^2)/S;

