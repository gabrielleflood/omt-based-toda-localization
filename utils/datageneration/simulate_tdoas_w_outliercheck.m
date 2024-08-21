function [tdoas_measured_orig,tdoas_true,r,s,nbr_missing, nbr_extra, ind_set_ok,tdoa_pairs_ok] = simulate_tdoas_w_outliercheck(R,S,sigma,r_bounds, s_bounds,P_missing,P_extra,all_tdoas_from_toas)
% calls the function simulate_tdoas to simulate tdoa measurements, but also
% checks that the setting is not "too bad" and re-generates new tdoa data
% if the setting is bad.

is_outlier_setting = 1;
while is_outlier_setting %&& kk == 1
    [tdoas_measured_orig,tdoas_true,r,s,nbr_missing, nbr_extra] = simulate_tdoas(R,S,sigma,r_bounds, s_bounds,P_missing,P_extra,all_tdoas_from_toas);
    max_inner_iter = 100;
    selection_iter = 1;
    while is_outlier_setting && selection_iter<=max_inner_iter
        [ind_set,tdoa_pairs] = find_3_index_pairs(R);
        
        [is_outlier_setting,CRLB_source_pos] = check_outlier_setting(ind_set,r,s);
        
        selection_iter = selection_iter+1;
    end
    
end
ind_set_ok = ind_set;
tdoa_pairs_ok = tdoa_pairs;

if is_outlier_setting
    'skumt'
end
