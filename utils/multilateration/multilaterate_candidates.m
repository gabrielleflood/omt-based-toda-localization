function [candidates, candidates_full, candidates_idx] = multilaterate_candidates(tdoas_measured, ind_set, r, imag_thresh,best_candidates_per_set,s)
% function that takes a RxR matrix of measured tdoas, a seet of 3 index set
% pair (values <R) and a threshold for whether to remove values that have
% to large imaginary part and returns a set of candidate positions for
% senders. The candidates are multilaterated from all tdoa peak
% combinations of the different receiver pairs in ind_set. The elements of
% tdoas_measured should contain the fields r1; r2; s or s1 and s2; N and 
% tdoas for the upper triangular matrix. if imag_thresh is <inf, then
% candidates for which the norm of the imaginary part is larger than
% imag_thresh are discarded. All candidates, with imaginary parts, are
% returned in candidates_full. r are the known receiver positions. the
% output candidates_idx tells which indices in candidates_full that the
% positions in candidates correspond to. 

if nargin<4 || isempty(imag_thresh); imag_thresh = inf; end
if nargin<5 || isempty(best_candidates_per_set); best_candidates_per_set = 0; end

% create matrices to save candidate positions in
candidates = [];
candidates_full = [];
candidates_idx = [];

% create the cartesian product of the tdoa peak index sets
chosen_tdoas(1) = tdoas_measured{ind_set(1,1),ind_set(1,2)}; % the tdoa measurements from the first random microphine pair
chosen_tdoas(2) = tdoas_measured{ind_set(2,1),ind_set(2,2)};
chosen_tdoas(3) = tdoas_measured{ind_set(3,1),ind_set(3,2)};
peak_ind = cartesianProd(1:chosen_tdoas(1).N, 1:chosen_tdoas(2).N, 1:chosen_tdoas(3).N);

% go through all tdoa combinations
count = 0;
for i = 1:size(peak_ind,1)
    curr_ind = peak_ind(i,:);
    udata = [ind_set'; ...
        chosen_tdoas(1).tdoas(curr_ind(1)) chosen_tdoas(2).tdoas(curr_ind(2)) chosen_tdoas(3).tdoas(curr_ind(3))];
    tmp2 = ind_set'; % i, j index
    chosen_r = r(:,tmp2(:));
    tmp4 = [chosen_r(:); udata(3,:)'.^2];
    try
        sols = solver_tdoa_pair_3d(tmp4);
        if best_candidates_per_set 
            % check how well the solutions agree with the three used
            % tdoa measurements
            candidates_full = [candidates_full sols];
            sols = sols(:,vecnorm(imag(sols)) < imag_thresh);
            sols = real(sols);
            % add something that removes doubled positions, that come from
            % conjugates?
            sols = remove_duplicates(sols,[]);
            tmp_tdoa_diff = [];
            for iii = 1:3
                tmp_tdoa=sqrt(sum((sols-repmat(r(:,ind_set(iii,2)),[1 size(sols,2)])).^2,1)) -...
                    sqrt(sum((sols-repmat(r(:,ind_set(iii,1)),[1 size(sols,2)])).^2,1));
                tmp_tdoa_diff(iii,1:size(sols,2)) = ( tmp_tdoa - chosen_tdoas(iii).tdoas(curr_ind(iii)) ).^2;
            end
            best_pos = find(vecnorm(tmp_tdoa_diff)<0.01);
            candidates = [candidates sols(:,best_pos)];
            candidates_idx = [candidates_idx count+best_pos];
            count = count + size(sols,2);
            
        else
            for k = 1:size(sols,2)
                count = count+1;
                % obs! what do we want to do if we get imaginary numbers? Only keep
                % real or check size of imginary part in some way?
                
                if norm(imag(sols(:,k))) < imag_thresh
                    candidates = [candidates real(sols(:,k)) ];
                    candidates_idx = [candidates_idx count];
                end
                candidates_full = [candidates_full (sols(:,k)) ];
            end
        end
    catch
        % something?
        warning('Minimal solver gave an error, continuing...')
    end
end