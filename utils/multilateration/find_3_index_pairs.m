function [ind_set,tdoa_pairs,sel] = find_3_index_pairs(R,sel)
% function that takes an integer R and finds three random index pairs in
% the set 1:Rx1:R. The pairs are chosen s.t. it has to contain at least
% four different indices, giving that the tdoa measurements for the pairs
% are not linearly dependent.
% if sel is sent in as a parameter, the selection is made from all pairs
% except there pairs, otherwise all

if nargin < 2 || isempty(sel) 
    sel = [];
end

% find out which index pairs there are and choose three of these for the
% multilateration
tdoa_pairs = cartesianProd(1:R,1:R);
remove_ind = tdoa_pairs(:,1)>=tdoa_pairs(:,2);
tdoa_pairs = tdoa_pairs(~remove_ind,:);

tdoa_pairs_old = tdoa_pairs;
tdoa_pairs(sel,:) = [];
sel_old = sel;

% choose 3 measurements to be used
sel = randperm(size(tdoa_pairs,1),3);
% Make sure that chosen set is independent and can be solved for
ind_set = tdoa_pairs(sel,:);
ind_set = unique(ind_set(:));
while length(ind_set) < 4
    sel = randperm(size(tdoa_pairs,1),3);
    ind_set = tdoa_pairs(sel,:);
    ind_set = unique(ind_set(:));
end
ind_set = tdoa_pairs(sel,:);

% just fix sel s.t. the correct sels are sent out if some are sent in
sel_old = sort(sel_old);
for i = 1:length(sel_old)
    pos_i = sel>=sel_old(i); 
    sel(pos_i) = sel(pos_i) +1;
end
tdoa_pairs = tdoa_pairs_old;
% OBS! Double check the sel-part before using it
