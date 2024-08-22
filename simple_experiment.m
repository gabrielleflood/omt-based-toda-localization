% script to localize S sources in an environment with R receivers usign
% TDOA measurements.

%% add paths, define some variables and simulate measurements

disp('Setting up the problem')

% add util path
addpath(genpath('utils'))

% define some variables, according to your specific case (measure or to be
% generated)
R = 12; % nbr of receivers
S = 3; % nbr of senders

% some variables needed for the method
imag_thresh = 0.8; % if the norm of the imaginary part of the multilateration solution is larger than this, the candidate is discarded
best_candidates_per_set = 1; % if 1 then only the best candidates from each multilateration (closest in euclidean measures) are saved. A threshold is set in multilaterate_candidates. Imag_thresh is used either way.
trash_prct = 95; % trash cost is constant, set to be this percentile of all other cost values. Leave emtpy to choose the constant below.
trash_cost = []; % if you want to choose the conastant trash cost, set the constant here and set trash_prct to be empty.
duplicate_threshold = 1e-2; % how close two candidates can be
nbr_candset_iter = 3; % the number of index sets that are used to produce candidates. Recommended to be >1 is there are missing peaks.
all_tdoas_from_toas = 0; % change this to 1 if you want all tdoa peaks given from all combinations of toas
savefile = []; % to save data in a certain file beween runs and also load this data and continue where it ended

cost_matrix_setup_function = @setup_cost_matrix; % the function to create the cost matrix. Can be changed for other cost functions
consider_outlier_settings = 0; % 1 if we should re-generate receiver and sender setups that are bad, 0 otherwise. Recommended to be 0.
use_bestpos_inM = 1; % decides whether only the best columns of M should be used when we make a binary version of M
SMTL_lambda = 0.5; % parameter for the comparison met


%% Things needed for data generation. Skip if you already have data
r_bounds = []; % low and high bounds in x-, y- and z-dir for where receivers can be (size 3x2)
s_bounds = []; % low and high bounds in x-, y- and z-dir for where senders can be (size 3x2)
P_missing = 0; % chance of a tdoa measurement being missed, should be 0 here
P_extra = 0; % chance of an extra tdoa measurement appearing, should be 0 here
sigma = 0.01; % stds of pertubations/noise to add to the measured tdoas

[tdoas_measured,tdoas_true,r,s,nbr_missing, nbr_extra] = simulate_tdoas(R,S,sigma,r_bounds, s_bounds,P_missing,P_extra,all_tdoas_from_toas); % simulate senders and recivers and compute tdoas
exists_gt = 1; % set to 0 if there is not ground truth to compare to

%% Read you data from file or put it in here in some way

% tdoas measured = ...
% r = ...

%  % if you know the ground truth of youre data
% tdoas_true = ...
% s = ...
% exists_gt = 1; % set to 0 if there is not ground truth to compare to


%%

if exists_gt
    % compute CRLB for full setup
    CRLB_mean = get_CRLB_mean(r,s,R,S,sigma);
end 

% find three index pairs for the multilateration. tdoa_pairs
% contains all possible pairs, not only selected.
candidates = [];
candidates_full = [];
candidates_idx = [];
ind_set_full = [];
sel_all = [];

%% Find the candidate set

disp('Creating the candidate set using multilateration solver.')

% create full candidate set by choosing nbr_candset_iter different toda
% triplets
for candset_iter = 1:nbr_candset_iter
    
    % sel contains the indices for the selected pairs
    if consider_outlier_settings
        [ind_set,tdoa_pairs,outlier_indicator,CRLB_temp(candset_iter,:),ind_set_full,sel] = find_3_index_pairs_w_outliercheck(R,r,s,ind_set_ok,tdoa_pairs_ok,ind_set_full,sel_all);
        sel_all = [sel_all sel];
    else
        [ind_set,tdoa_pairs,sel] = find_3_index_pairs(R,sel_all);
        sel_all = [sel_all sel];
    end
    % do multilateration on all combinations of tdoa peaks for the
    % chosen receiver pairs
    [candidates_i, candidates_full_i, candidates_idx_i] = multilaterate_candidates(tdoas_measured, ind_set, r, imag_thresh, best_candidates_per_set);
    candidates = [candidates candidates_i];
    candidates_idx = [candidates_idx candidates_idx_i+size(candidates_full,2)];
    candidates_full = [candidates_full candidates_full_i];
end

if exists_gt
    % Check distance from true source positions to best candidates
    dists_temp = get_dists(s,candidates);
end

% remove 'dupplicate' candidates that are very close
candidates = remove_duplicates(candidates,duplicate_threshold);


%% Run the OMT algorithm

disp('Solving the OMT problem using all TDOA measurements.')

% Re-save the data for easier listing of all tdoa measurements, here all data
tdoas_listed = tdoa_matrix_to_list(tdoas_measured, tdoa_pairs,all_tdoas_from_toas);

% create the cost matrix C
cost_mat = cost_matrix_setup_function(tdoas_listed,candidates,r,trash_prct,trash_cost);

epsilon = 1e-7; % setting for omt solver
eta = 1; % setting for omt solver
for k_ot_trial = 1:5
    [M,cluster_vec,trash_mass] = omt_clustering_w_entropy_reg_penalty_logstable_w_trashstate(cost_mat(:,1:end-1),cost_mat(:,end),epsilon,eta,R);
    M_out = M;
    mu_out = cluster_vec;
    if sum(cluster_vec>.8)==S
        break
    end
end
if sum(trash_mass)>0
    disp('trash')
end

% Omt problem solved. Pick out the top S found senders
[~,best_pos] = sort(mu_out,'descend');
try
    s_found = candidates(:,best_pos(1:S));
catch
    % if not enough senders were found (eg due to large noise)
    s_found = candidates;
end

% Reestimate association 
candidates_estimated = candidates(:,best_pos(1:S));
cost_mat = cost_matrix_setup_function(tdoas_listed,candidates_estimated,r,trash_prct,trash_cost);
[M,trash_mass_2] = omt_clustering_w_entropy_reg_penalty_logstable_reestimate(cost_mat(:,1:end-1),cost_mat(:,end),1e2*epsilon,R);
M_out = 0*M_out;
M_out(:,best_pos(1:S)) = M;

if exists_gt
    % find mapping between gt senders and found senders
    s_eucl_dist = pdist2(s',s_found');
end

%% Evaluate the association and refine solution using local optimization

disp('Running local refinements')

% Refine sender position solution using association
if size(s_found,2) == S
    if exists_gt
        [s_refined,mean_distance_refined] = refine_using_association(tdoas_listed,S,M_out, best_pos,  candidates, r,s);
    else
        [s_refined] = refine_using_association(tdoas_listed,S,M_out, best_pos,  candidates, r);
    end

    % Check distance from true source positions to best candidates
    dists_temp = get_dists(s_refined,candidates);
end

%% Print some results

disp('Done')
disp(' ')

disp('--------RESULTS--------')
disp('The found sender positions (without order) are:')
if size(s_found,2) == S
    disp(s_refined)
else
    disp(s_found)
end

if exists_gt
    disp('While the true sender positions are:')
    disp(s)
    disp('and the mean Euclidean distance from gt to found sources were')
    disp(mean(dists_temp))
    disp('with a CRLB mean value of')
    disp(CRLB_mean)
end


