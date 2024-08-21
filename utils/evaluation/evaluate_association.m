function [association_rate, trash_rate] = evaluate_association(tdoas_listed, S,M_in, best_pos, use_bestpos_inM)
% function to check how many TDOA that are correctly associated.
% tdoas_listed are the listed tdoas in the same order as in the rows of
% M_in. best_pos are the positions for the best candidate positions. The
% variable use_bestpos_inM tells whether we only should count correct TDOA
% associations if they belong to a selected candidatefrom mu
% association_rate are the rate of tdoas that were associated with the
% correct source (or trash if the measurements are false)
% trash_rate are the rate of false tdoas that are correctly associated with
% trash. If there are no false tdoas, this value is NaN

if ~use_bestpos_inM
    % first, make the M matrix binary. 
    best_pos = []; % use this option if M can be 1 even if this does not correspond to one of the top activations in mu
end

% change M into binary. Add s.t. also trash bin can be activated, not only
% best pos
M = create_binary_M(M_in, [best_pos; size(M_in,2)]);
%M = M.*(M_in>0.7);
M = (M_in>0.7)+0;

for i = 1:S
    % find the tdoa indices that represent tdoas from sender i
    s_ind = ([tdoas_listed.s_true] == i)';
    
    % check how many of these tdoas that are associated with each candidate
    candidate_count = sum(M(s_ind,:));
    
    % find the max in candidate_count, i.e. the candidate that most of the
    % TDOAs for sender i are associated to
    %[~,max_cand] = max(candidate_count(1:end-1));
    [~,max_cand] = max(candidate_count(1:end));
    
    % then compute how many of the candidates that are associated with this
    % candidate and divide by the gt number of TDOA originating from source
    % i
%     correct_association(i) = candidate_count(max_cand)/sum(s_ind);
    correct_association(i) = candidate_count(max_cand);%/sum(s_ind);
    
    % also calculate how many TDOAs that are associated with trahs
    % trash_association_bad(i) = candidate_count(end)/sum(s_ind);
end

% continue counting how many peaks that were correctly classifies as trash 
false_ind = isnan([tdoas_listed.s_true])';
% By this we should have gone through/considered all tdoas, since all come
% from either one of the S true sources, or are false, i.e. have NaN as
% s_true

% check how many of these tdoas that are (correctly) associated with trash
% trash_count = sum(M(false_ind,end));
% trash_count = sum(M(:,end));
% trash_count = sum(M(:,end))/sum(false_ind);


% this is the rate of those tdoas that were false that are associated with
% trash
trash_rate = sum(M(false_ind,end))/sum(false_ind);


% This is the number of tdoas that were false that are (ocorrectly)
% associated with trash
correct_association(end+1) = sum(M(false_ind,end));

% count how many of the  trash classified tdoas that were actually trash


%association_rate = mean(correct_association);

% to get the association rate, divide the sum of all correct associations
% (for the different sources/NaN) with the number of tdoas in total
association_rate = sum(correct_association)/length(tdoas_listed);

% we don't need this anymore
%trash_rate_bad = mean(trash_association_bad);