% script to generate a plot for how the error in euclidean distance from gt
% positions change with increased pertubations in tdoa measurements
% inspired by experiment_error_wrt_tdoapertubation.m and
% experiment_error_wrt_tdoapertubation_CRLBTEST.
% created 2024-08-21

%% add paths, define some variables and simulate measurements
% add util path
clear,clc,close all
addpath(genpath('utils'))

% rng(2*pi) % if a certain rng should be used

% define some variables
R = 12; % nbr of receivers
S = 3; % nbr of senders
r_bounds = []; % low and high bounds in x-, y- and z-dir for where receivers can be (size 3x2)
s_bounds = []; % low and high bounds in x-, y- and z-dir for where senders can be (size 3x2)
P_missing = 0; % chance of a tdoa measurement being missed, should be 0 here
P_extra = 0; % chance of an extra tdoa measurement appearing, should be 0 here
all_tdoas_from_toas = 0; % change this to 1 if you want all tdoa peaks given from all combinations of toas
imag_thresh = 0.8; % if the norm of the imaginary part of the multilateration solution is larger than this, the candidate is discarded
best_candidates_per_set = 1; % if 1 then only the best candidates from each multilateration (closest in euclidean measures) are saved. A threshold is set in multilaterate_candidates. Imag_thresh is used either way. 

trash_prct = 95; % trash cost is constant, set to be this percentile of all other cost values. Leave emtpy to choose the constant below.
trash_cost = []; % if you want to choose the conastant trash cost, set the constant here and set trash_prct to be empty.

sigmas_all = 0:0.002:0.02; % the set of stds of pertubations to add to the measured tdoas
sigmas_all = [0.01 0.02];
nbr_iter = 1; %100; % nbr of noise realization iterations to average over
fail_thresh = 0.8; % when the euclidean distance between any two gt and found source nodes is larger than this, it is considered a fail. for producing plots later.
duplicate_threshold = 1e-2; % how close two candidates can be
nbr_candset_iter = 3; % the number of index sets that are used to produce candidates. Recommended to be >1 is there are missing peaks.
savefile = []; % to save data in a certain file beween runs and also load this data and continue where it ended

cost_matrix_setup_function = @setup_cost_matrix; % the function to create the cost matrix. Can be changed for other cost functions
                                    
consider_outlier_settings = 0; % 1 if we should re-generate receiver and sender setups that are bad, 0 otherwise. Recommended to be 0.
use_bestpos_inM = 1; % decides whether only the best columns of M should be used when we make a binary version of M

SMTL_lambda = 0.5; % parameter for the comparison method SMTL



%% Create some varaibles to save results etc in 

too_few_s_found = []; % keeps track of whether the multilateration gives fewer candidate positions than S

%% Iterate over all sigmas and solve the problem

% if the experiment has already been run and some data has been saved
% already, then continue saving in the same results variable. 
kk_start = 1;
if ~isempty(savefile) && isfile(savefile) 
    load(savefile)
    kk_start = size(results,1);
end

% iterate over number of noise realisations
for kk = kk_start:nbr_iter
  
    % if we want to add different noise on the same estimated sensor
    % positions, first simulate data without noise and then add it here
    sigma = 0;
            
    if consider_outlier_settings
        % Check if there exists minimal sets that dont give outliers
        [tdoas_measured_orig,tdoas_true,r,s,nbr_missing, nbr_extra,ind_set_ok,tdoa_pairs_ok] = simulate_tdoas_w_outliercheck(R,S,sigma,r_bounds, s_bounds,P_missing,P_extra,all_tdoas_from_toas);
    else
       [tdoas_measured_orig,tdoas_true,r,s,nbr_missing, nbr_extra] = simulate_tdoas(R,S,sigma,r_bounds, s_bounds,P_missing,P_extra,all_tdoas_from_toas); % simulate senders and recivers and compute tdoas
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iterate over different noise levels
    for ii = 1:length(sigmas_all)
        disp(['Outer iteration nbr ' num2str(kk) '/' num2str(nbr_iter) ', inner (sigma) nbr ' num2str(ii) '/' num2str(length(sigmas_all))])
        
        % use original tdoa values, to not accumulate pertubations
        tdoas_measured = tdoas_measured_orig;
        % add the noise
        sigma = sigmas_all(ii);
        for i = 1:R
            for j = i+1:R
                tdoas_measured{i,j}.tdoas = tdoas_measured{i,j}.tdoas + sigma*randn(size(tdoas_measured{i,j}.tdoas));
            end
        end
        
        % compute CRLB for full setup for plotting
        CRLB_mean = get_CRLB_mean(r,s,R,S,sigma);
        
        % find three index pairs for the multilateration. tdoa_pairs 
        % contains all possible pairs, not only selected.
        candidates = [];  
        candidates_full = []; 
        candidates_idx = []; 
        ind_set_full = [];
        sel_all = [];
        
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
            [candidates_i, candidates_full_i, candidates_idx_i] = multilaterate_candidates(tdoas_measured, ind_set, r, imag_thresh, best_candidates_per_set,s);
            candidates = [candidates candidates_i]; 
            candidates_idx = [candidates_idx candidates_idx_i+size(candidates_full,2)]; 
            candidates_full = [candidates_full candidates_full_i]; 
        end  
        
        % Check distance from true source positions to best candidates
        dists_temp = get_dists(s,candidates);
        % remove 'dupplicate' candidates that are very close
        candidates = remove_duplicates(candidates,duplicate_threshold);

        % Re-save the data for easier listing of all tdoa measurements, here all data
        tdoas_listed = tdoa_matrix_to_list(tdoas_measured, tdoa_pairs,all_tdoas_from_toas);
        
        %%%%%%%%%%%% Solve OMT problem %%%%%%%%%%%%%%
        disp('Solving the OMT problem using all TDOA measurements.')
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
            too_few_s_found = [too_few_s_found [kk;ii;size(candidates,2) ]];
        end

        %%%%%%% Reestimate association %%%%%%
        candidates_estimated = candidates(:,best_pos(1:S));
        cost_mat = cost_matrix_setup_function(tdoas_listed,candidates_estimated,r,trash_prct,trash_cost);
        [M,trash_mass_2] = omt_clustering_w_entropy_reg_penalty_logstable_reestimate(cost_mat(:,1:end-1),cost_mat(:,end),1e2*epsilon,R);
        M_out = 0*M_out;
        M_out(:,best_pos(1:S)) = M;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % find mapping between gt senders and found senders
        s_eucl_dist = pdist2(s',s_found');
           
        % Evaluate the association
        [association_rate, trash_rate] = evaluate_association(tdoas_listed, S,[M_out trash_mass_2], best_pos(1:S), use_bestpos_inM);
        % Refine sender position solution using association
        if size(s_found,2) == S
            [s_refined,mean_distance_refined] = refine_using_association(tdoas_listed,S,M_out, best_pos,  candidates, r,s);
        end
        
        % Do SMTL calculations. Note that this requires CVX.
        ggrid = candidates; % OBS!!!!!
        [SMTL_mean_distance] = get_SMTL_mean_distance(ggrid,r,s,S,tdoas_measured, SMTL_lambda);
           
        % save some results
        % row (kk) represents iteration, col (ii) represents sigma in sigmas_all(ii)
        results(kk,ii).s_true = s;
        results(kk,ii).mu = mu_out;
        results(kk,ii).mu_above_0_001 = sum(mu_out>0.001);
        results(kk,ii).mu_above_0_2 = sum(mu_out>0.2);
        results(kk,ii).s_found = s_found;
        results(kk,ii).candidates = candidates;
        results(kk,ii).tdoas_measured = tdoas_measured;
        results(kk,ii).ind_set = ind_set;
        results(kk,ii).tdoa_pairs = tdoa_pairs;
        results(kk,ii).mean_distance = sum(min(s_eucl_dist,[],1))/S;
        results(kk,ii).potential_fail = sum(min(s_eucl_dist,[],1) > fail_thresh) >0;
        results(kk,ii).CRLB_mean = CRLB_mean;        
        results(kk,ii).dist_to_candidates = mean(dists_temp);
        results(kk,ii).SMTL_mean_distance = SMTL_mean_distance;
        results(kk,ii).association_rate = association_rate;
        results(kk,ii).trash_rate = trash_rate;
        
        if size(s_found,2) == S
            results(kk,ii).s_refined = s_refined;
            results(kk,ii).mean_distance_refined = mean_distance_refined;
        end

        if consider_outlier_settings
            results(kk,ii).CRLB_over_iterations = mean((min(CRLB_temp)));
            [~,CRLB_source_pos_full] = check_outlier_setting(ind_set_full,r,s);
            results(kk,ii).CRLB_over_iterations_full_set = mean(CRLB_source_pos_full);
            results(kk,ii).outlier_indicator = outlier_indicator;  
        end

        if ~isempty(savefile)
            clear tdoas_measured rec1_rec2_pair_cell counter ...
                prop_speed_nominal sigma2_noise_nominal CRLB_source_pos ...
                CRLB_mat tot_FIM CRLB_mean candidates candidates_full ...
                candidates_idx ind_set tdoa_pairs candidates_i  ...
                candidates_full_i candidates_idx_i tdoas_listed nbr_toas ...
                card_grid nbr_targets nbr_receivers cost_mat M_out ...
                mu_out omt_time s_found s_eucl_dist 
            save(savefile,'results')
        end
    end
end


%% Compute results for plotting

% fail_count = zeros(1,length(sigmas_all));
mean_distance_refined = zeros(1,length(sigmas_all));
mean_distance_SMTL = zeros(1,length(sigmas_all));
mean_association_rate = zeros(1,length(sigmas_all));
% mean_trash_rate = zeros(1,length(sigmas_all));
for ii = 1:length(sigmas_all)
%     [results(:,ii).mean_distance]
    tmp = [results(:,ii).mean_distance];
    mean_distance(1,ii) = mean([results(:,ii).mean_distance]);
    % mean_distance_nofail(1,ii) = mean(tmp(~[results(:,ii).potential_fail]));
    % fail_count(ii) = fail_count(ii) + sum([results(:,ii).potential_fail]);
    mean_distance_refined(1,ii) = mean([results(:,ii).mean_distance_refined]);
    mean_distance_SMTL(1,ii) = mean([results(:,ii).SMTL_mean_distance]);
    mean_CRLB(ii) = mean([results(:,ii).CRLB_mean]);    
    mean_association_rate(ii) = mean([results(:,ii).association_rate]);
    % mean_trash_rate(ii) = mean([results(:,ii).trash_rate]);
end

% if results are to be save, also save setting variables and studd to plot
if ~isempty(savefile)
    save(savefile)
end

%% Plot results

figure(120)
%tcl = tiledlayout(1,2);
% tcl = tiledlayout(2,1);
%nexttile(tcl)

yyaxis left
semilogy(sigmas_all,[mean_distance],'^-','LineWidth',2,'MarkerSize',6)
hold on
semilogy(sigmas_all,mean_distance_refined,'v-','LineWidth',2,'MarkerSize',6)
semilogy(sigmas_all,(mean_CRLB),'--','LineWidth',2,'MarkerSize',6)
semilogy(sigmas_all,mean_distance_SMTL,'o-','LineWidth',2,'MarkerSize',6)

xlabel('$\sigma/\rho \;\;\mathrm{(meters)}$','interpreter','latex')
ylabel('$\bf{E} [ \| \hat{s} - s \|_2  ]$','interpreter','latex')
ylim([min([mean_CRLB mean_distance_refined]) 2*max(mean_distance_SMTL)]);
set(gca,'FontSize',14)
grid on
ax = gca;
ax.GridAlpha = 0.3;

yyaxis right 
plot(sigmas_all,mean_association_rate,'-s','LineWidth',2,'MarkerSize',6)
% plot(sigmas_all,mean_trash_rate,'-x','LineWidth',2.5,'MarkerSize',8)
ylabel('rate')
ylim([0 1])
% xlim([min(sigmas_all) max(sigmas_all)])
xlim([0 max(sigmas_all)])
set(gca,'FontSize',14)
set(gca,'TickLength',[0.05 2.5*0.05])

legend('Proposed','Proposed, refined','Root CRLB','SMTL','Assoc. rate',...
   'location','southoutside','NumColumns',2,'FontSize',14); %,'FontWeight','bold'