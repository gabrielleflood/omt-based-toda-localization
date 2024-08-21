function tdoas_listed = tdoa_matrix_to_list(tdoas_measured, tdoa_pairs,all_tdoas)

tdoas_listed = [];

% first, create a collection of tdoas, where we have all tdoas measurements
% in a vector
for i = 1:size(tdoa_pairs,1)
    curr_pair = tdoas_measured{tdoa_pairs(i,1), tdoa_pairs(i,2)};
    start_ind = length(tdoas_listed);
    for j = 1:length(curr_pair.tdoas)
        tdoas_listed(start_ind+j).tdoas = curr_pair.tdoas(j); % tdoa measurements
        tdoas_listed(start_ind+j).r1 = curr_pair.r1; % tdoa receiver 1
        tdoas_listed(start_ind+j).r2 = curr_pair.r2; % tdoa receiver 2
        tdoas_listed(start_ind+j).peak_nbr = j; % which peak in the order the tdoa, to be able to go back to tdoas_measures
        if ~all_tdoas
            tdoas_listed(start_ind+j).s_true = curr_pair.s(j); % the true sender for each tdoa
        else
            tdoas_listed(start_ind+j).s1_true = curr_pair.s1(j); % the true sender 1 for each tdoa
            tdoas_listed(start_ind+j).s2_true = curr_pair.s2(j); % the true sender 2 for each tdoa
        end
    end
end

