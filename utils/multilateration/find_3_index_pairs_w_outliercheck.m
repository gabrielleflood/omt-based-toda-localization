function [ind_set,tdoa_pairs,cond_flag,CRLB_source_pos,ind_set_full,sel] = find_3_index_pairs_w_outliercheck(R,r,s,ind_set_ok,tdoa_pairs_ok,ind_set_full,sel_all)
% calls the function find_3_index_pairs to choose tdoa pairs, but also
% checks that the setting is not "too bad" and re-generates new tdoa data
% if the setting is bad.

cond_flag = 1;
cond_iter = 1;
max_cond_iter = 100;
while cond_flag && cond_iter <= max_cond_iter
    % ind_set = find_3_index_pairs(R);
    [ind_set,tdoa_pairs,sel] = find_3_index_pairs(R,sel_all);
    [is_outlier_setting,CRLB_source_pos] = check_outlier_setting(ind_set,r,s);
    if ~is_outlier_setting
        cond_flag = 0;
    end
    cond_iter = cond_iter+1;
end
if cond_flag
    'skumt2'
    ind_set = ind_set_ok;
    tdoa_pairs = tdoa_pairs_ok;
end
[~,CRLB_source_pos] = check_outlier_setting(ind_set,r,s);
CRLB_source_pos = CRLB_source_pos';
ind_set_full = [ind_set_full;ind_set];